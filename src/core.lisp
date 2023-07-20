(defpackage :cl-mpm
  (:use :cl
        :cl-mpm/particle
        :cl-mpm/mesh
        :cl-mpm/utils
        :cl-mpm/fastmath
;        :cl-mpm/shape-function
        )
  (:import-from
    :magicl tref .+ .-
    )
  (:import-from
   :cl-mpm/forces det-int-force det-ext-force
   )
  (:export
    #:mpm-sim
    #:make-mpm-sim
    #:make-shape-function-linear
    #:sim-mesh
    #:sim-mps
    #:sim-bcs
    #:sim-dt
    #:sim-damping-factor
    #:sim-mass-filter
    #:post-stress-step
    ))
(declaim (optimize (debug 0) (safety 0) (speed 3)))
;    #:make-shape-function
(in-package :cl-mpm)

(eval-when (:compile-toplevel :load-toplevel :execute)
  #+cl-mpm-pic (print "Compiled with PIC")
  #-cl-mpm-pic (print "Compiled with FLIP")
  #+cl-mpm-fbar (print "Compiled with FBAR")
  #-cl-mpm-fbar (print "Compiled without FBAR")
  #+cl-mpm-special (print "Compiled with special hooks")
  #-cl-mpm-special (print "Compiled without special hooks")
)


(defclass mpm-sim ()
  ((dt
     :accessor sim-dt
     :initarg :dt)
   (mesh
     :accessor sim-mesh
     :initarg :mesh)
   (mps
     :accessor sim-mps
     :initarg :mps)
   (bcs
     :accessor sim-bcs
     :initarg :bcs
     :initform (make-array 0))
   (bcs-force
    :accessor sim-bcs-force
    :initarg :bcs-force
    :initform (make-array 0))
   (damping-factor
     :type double-float
     :accessor sim-damping-factor
     :initarg :damping-factor
     :initform 0d0)
   (mass-scale
    :type double-float
    :accessor sim-mass-scale
    :initarg :mass-scale
    :initform 1d0)
   (allow-mp-split
    :type boolean
    :accessor sim-allow-mp-split
    :initarg :splitting
    :initform nil
    )
   (enable-damage
    :type boolean
    :accessor sim-enable-damage
    :initarg :enable-damage
    :initform nil
    )
   (nonlocal-damage
    :type boolean
    :accessor sim-nonlocal-damage
    :initform t
    )
   (allow-mp-damage-removal
    :type boolean
    :accessor sim-allow-mp-damage-removal
    :initarg :damage-removal
    :initform nil
    )
   (mass-filter
     :type double-float
     :accessor sim-mass-filter
     :initarg :mass-filter
     :initform 1d-15))
  (:documentation "A self contained mpm simulation"))
(defclass mpm-sim-usf (mpm-sim)
  ()
  (:documentation "Explicit simulation with update stress first update"))
(defclass mpm-sim-usl (mpm-sim)
  ()
  (:documentation "Explicit simulation with update stress last update"))

(defun make-mpm-sim (size resolution dt shape-function &key (sim-type 'mpm-sim-usf))
  (make-instance sim-type
                 :dt (coerce dt 'double-float)
                 :mesh (make-mesh size resolution shape-function)
                 :mps '()))
(defun check-mps (mps)
  "Function to check that stresses and positions are sane"
  (loop for mp across mps
        do
        (progn
          (with-accessors ((vel cl-mpm/particle::mp-velocity)) mp
            (loop for i from 0 to 1
                 do (progn
                      (when (equal (abs (magicl:tref vel i 0)) #.sb-ext:double-float-positive-infinity)
                        (print mp)
                        (error "Infinite velocity found"))
                      (when (> (abs (magicl:tref vel i 0)) 1e10)
                        (print mp)
                        (error "High velocity found"))))))))

(defgeneric update-sim (sim)
  (:documentation "Update an mpm simulation by one timestep"))

(defmethod update-sim ((sim mpm-sim))
  (declare (cl-mpm::mpm-sim sim))
  (with-slots ((mesh mesh)
               (mps mps)
               (bcs bcs)
               (bcs-force bcs-force)
               (dt dt)
               (mass-filter mass-filter)
               (split allow-mp-split)
               (enable-damage enable-damage)
               (nonlocal-damage nonlocal-damage)
               (remove-damage allow-mp-damage-removal)
               )
                sim
    (declare (type double-float mass-filter))
                (progn
                    (reset-grid mesh)
                    (p2g mesh mps)
                    (when (> mass-filter 0d0)
                      (filter-grid mesh (sim-mass-filter sim)))
                    (apply-bcs mesh bcs-force dt) ;;Update nodes
                    (update-node-kinematics mesh dt)
                    (p2g-force mesh mps)
                    (update-node-forces mesh (sim-damping-factor sim) dt(sim-mass-scale sim))
                    ;Apply vel bcs
                    (apply-bcs mesh bcs dt)
                    ;;G2p + Particle update
                    (g2p mesh mps dt)
                    (update-stress mesh mps dt)
                    ;; (when enable-damage
                      ;; (cl-mpm/damage::calculate-damage mesh
                      ;;                                  mps
                      ;;                                  dt
                      ;;                                  50d0
                      ;;                                  nonlocal-damage
                      ;;                                  ))
                    ;; (update-particle mps dt)
                    ;;Update stress last
                    (when remove-damage
                      (remove-material-damaged sim))
                    (when split
                      (split-mps sim))
                    (check-mps mps)
                    )))
(defmethod update-sim ((sim mpm-sim-usf))
  (declare (cl-mpm::mpm-sim-usf sim))
  (with-slots ((mesh mesh)
               (mps mps)
               (bcs bcs)
               (bcs-force bcs-force)
               (dt dt)
               (mass-filter mass-filter)
               (split allow-mp-split)
               (enable-damage enable-damage)
               (nonlocal-damage nonlocal-damage)
               (remove-damage allow-mp-damage-removal)
               )
                sim
    (declare (type double-float mass-filter))
                (progn
                    (reset-grid mesh)
                    (p2g mesh mps)
                    (when (> mass-filter 0d0)
                     (filter-grid mesh (sim-mass-filter sim)))
                    (update-node-kinematics mesh dt )
                    (apply-bcs mesh bcs dt)
                    (update-stress mesh mps dt)
                    ;; (when enable-damage
                    ;;  (cl-mpm/damage::calculate-damage mesh
                    ;;                                   mps
                    ;;                                   dt
                    ;;                                   50d0
                    ;;                                   nonlocal-damage
                    ;;                                   ))
                    ;Map forces onto nodes
                    (p2g-force mesh mps)
                    (apply-bcs mesh bcs-force dt)
                    (update-node-forces mesh (sim-damping-factor sim) dt (sim-mass-scale sim))
                    ;Reapply velocity BCs
                    (apply-bcs mesh bcs dt)
                    ;Also updates mps inline
                    (g2p mesh mps dt)

                    (when remove-damage
                     (remove-material-damaged sim))
                    (when split
                     (split-mps sim))
                    (check-mps mps)
                    )))
(defmethod update-sim ((sim mpm-sim-usl))
  (declare (cl-mpm::mpm-sim sim))
  (with-slots ((mesh mesh)
               (mps mps)
               (bcs bcs)
               (bcs-force bcs-force)
               (dt dt)
               (mass-filter mass-filter)
               (split allow-mp-split)
               (enable-damage enable-damage)
               (nonlocal-damage nonlocal-damage)
               (remove-damage allow-mp-damage-removal)
               )
                sim
    (declare (type double-float mass-filter))
                (progn
                    (reset-grid mesh)
                    ;; Map momentum to grid
                    (p2g mesh mps)
                    ;;Reset nodes below our mass-filter
                    (when (> mass-filter 0d0)
                      (filter-grid mesh (sim-mass-filter sim)))
                    ;;Turn momentum into velocity
                    (update-node-kinematics mesh dt)
                    (p2g-force mesh mps)
                    (apply-bcs mesh bcs-force dt)
                    ;;Update our nodes after force mapping
                    (update-node-forces mesh (sim-damping-factor sim) dt(sim-mass-scale sim))
                    ;Apply velocity bcs
                    (apply-bcs mesh bcs dt)
                    ;;Grid to particle mapping
                    (g2p mesh mps dt)

                    ;; ;;2nd round of momentum mapping
                    (reset-grid mesh)
                    (p2g mesh mps)
                    (when (> mass-filter 0d0)
                     (filter-grid mesh (sim-mass-filter sim)))
                    (update-node-kinematics mesh dt)
                    (apply-bcs mesh bcs dt)

                    ;;Update stress last
                    (update-stress mesh mps dt)
                    ;; (when enable-damage
                    ;;   (cl-mpm/damage::calculate-damage mesh
                    ;;                                    mps
                    ;;                                    dt
                    ;;                                    50d0
                    ;;                                    nonlocal-damage
                    ;;                                    ))
                    ;; (cl-mpm/damage::calculate-damage mesh
                    ;;                                  mps
                    ;;                                  dt
                    ;;                                  25d0)

                    (when remove-damage
                      (remove-material-damaged sim))
                    (when split
                      (split-mps sim))
                    (check-mps mps)
                    )))

(defgeneric iterate-over-neighbours-shape (mesh shape-func mp func)
  (:documentation "For a given shape function iterate over relevant nodes and call func with mesh, mp, node, weight and gradients"))

(declaim
(inline iterate-over-neighbours)
         (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle function) (values)) iterate-over-neighbours))
(defun iterate-over-neighbours (mesh mp func)
  "For a given mesh and mp, iterate over all the neighbours with
weight greater than 0, calling func with the mesh, mp, node, svp, and grad"
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (if (> (length (cl-mpm/particle::mp-cached-nodes mp)) 0)
      (iterate-over-neighbours-cached mesh mp func)
      (create-node-cache mesh mp func))
  ;; (iterate-over-neighbours-shape-gimp mesh mp func)
  ;; (iterate-over-neighbours-shape-linear mesh mp func)
  ;; (iterate-over-neighbours-shape mesh (cl-mpm/mesh:mesh-shape-func mesh) mp func)
  (values)
  )

(declaim (inline create-node-cache))
(defun create-node-cache (mesh mp func)
  (declare (function func))
  (with-accessors ((nodes cl-mpm/particle::mp-cached-nodes))
      mp
    (iterate-over-neighbours-shape-gimp
       mesh mp
       (lambda (mesh mp node svp grads fsvp fgrads)
         (vector-push-extend
          (cl-mpm/particle::make-node-cache
           :node node
           :weight svp
           :grads grads
           :weight-fbar fsvp
           :grads-fbar fgrads
           )
          nodes)
         (funcall func mesh mp node svp grads fsvp fgrads)
         ))))

(declaim (inline iterate-over-neighbours-cached))
(defun iterate-over-neighbours-cached (mesh mp func)
  (declare (function func))
  "If a node iteration cache has been generated we loop over the data list"
  (loop for nc across (cl-mpm/particle::mp-cached-nodes mp)
        do
         (funcall func mesh mp
                  (cl-mpm/particle::node-cache-node nc)
                  (cl-mpm/particle::node-cache-weight nc)
                  (cl-mpm/particle::node-cache-grads nc)
                  (cl-mpm/particle::node-cache-weight-fbar nc)
                  (cl-mpm/particle::node-cache-grads-fbar nc)
                  )))

(defmethod iterate-over-neighbours-shape (mesh (shape-func cl-mpm/shape-function:shape-function-linear) mp func)
  (iterate-over-neighbours-shape-linear mesh mp func))

(defmethod iterate-over-neighbours-shape (mesh (shape-func cl-mpm/shape-function:shape-function-bspline) mp func)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/shape-function::shape-function-bspline shape-func))
  (iterate-over-neighbours-shape-bspline mesh mp func))

(declaim (inline iterate-over-neighbours-shape-linear))
(defun iterate-over-neighbours-shape-linear (mesh mp func)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (pos-vec (cl-mpm/particle:mp-position mp))
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
           (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec #'floor)))
      (loop for dx from 0 to 1
            do (loop for dy from 0 to 1
                     do (let* ((id (mapcar #'+ pos-index (list dx dy))))
                          (when (cl-mpm/mesh:in-bounds mesh id)
                            (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                   (node (cl-mpm/mesh:get-node mesh id))
                                   (weights (mapcar (lambda (x) (cl-mpm/shape-function::shape-linear x h)) dist))
                                   (weight (reduce #'* weights))
                                   (grads (mapcar (lambda (d w) (* (cl-mpm/shape-function::shape-linear-dsvp d h) w))
                                                  dist (nreverse weights))))
                              (when (< 0d0 weight)
                                (funcall func mesh mp node weight grads nil nil))))))))))

(declaim (inline iterate-over-neighbours-point-linear)
         (ftype (function (cl-mpm/mesh::mesh magicl:matrix/double-float function) (values))
                iterate-over-neighbours-point-linear)
         )
(defun iterate-over-neighbours-point-linear (mesh position func)
  (declare (cl-mpm/mesh::mesh mesh))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (pos-vec position)
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
           (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec #'floor)))
      (declare (dynamic-extent pos pos-index pos-vec))
      (loop for dx from 0 to 1
            do (loop for dy from 0 to 1
                     do (let* ((id (mapcar #'+ pos-index (list dx dy))))
                          (declare (dynamic-extent id))
                          (when (cl-mpm/mesh:in-bounds mesh id)
                            (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                   (node (cl-mpm/mesh:get-node mesh id))
                                   (weights (mapcar (lambda (x) (cl-mpm/shape-function::shape-linear x h)) dist))
                                   (weight (reduce #'* weights))
                                   (grads (mapcar (lambda (d w) (* (the double-float (cl-mpm/shape-function::shape-linear-dsvp d h))
                                                                   (the double-float w)))
                                                  dist (nreverse weights))))
                              (declare (double-float weight)
                                       (dynamic-extent dist weights))
                              (when (< 0d0 weight)
                                (funcall func mesh node weight grads))))
                          ))))))

(declaim (inline iterate-over-neighbours-shape-gimp)
         (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle function) (values))
                iterate-over-neighbours-shape-gimp))
(defun iterate-over-neighbours-shape-gimp (mesh mp func)
  (declare (type cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp)
           (function func)
           ;; (optimize (speed 0) (safety 3) (debug 3))
           )
  (progn
    (with-accessors ((pos-vec cl-mpm/particle:mp-position)
                     (d0 cl-mpm/particle::mp-domain-size))
        mp
      (let ((h (the double-float (cl-mpm/mesh:mesh-resolution mesh))))
        (flet ((center-diff (x)
                 (declare (double-float x h))
                 (- x (the double-float (* h (the double-float
                                                  (fround (the double-float (/ x h))))))))
               (center-id (x)
                 (round x h))
               )
          (let* ((px (tref pos-vec 0 0))
                 (py (tref pos-vec 1 0))
                 (ix (the fixnum (truncate (center-id px))))
                 (iy (the fixnum (truncate (center-id py))))
                 (cx (center-diff px))
                 (cy (center-diff py))
                 (dox (* 0.5d0 (the double-float (tref d0 0 0))))
                 (doy (* 0.5d0 (the double-float (tref d0 1 0))))
                 (dxf (the fixnum (truncate (ffloor   (- cx dox) h))))
                 (dxc (the fixnum (truncate (fceiling (+ cx dox) h))))
                 (dyf (the fixnum (truncate (ffloor   (- cy doy) h))))
                 (dyc (the fixnum (truncate (fceiling (+ cy doy) h))))
                 )
            ;; (declare (dynamic-extent pos pos-index))
            (declare (type double-float h cx cy dox doy px py)
                     (type integer dxf dxc dyf dyc ix iy)
                     ;; (type list pos pos-index)
                     )
            (loop for dx from dxf
                    to dxc
                  do (loop for dy from dyf
                             to dyc
                           do
                              (let* ((id (list (+ ix dx) (+ iy dy))))
                                (declare (dynamic-extent id))
                                (when (cl-mpm/mesh:in-bounds mesh id)
                                  (let* ((dist (list (- cx (* h dx))
                                                     (- cy (* h dy))))
                                         (domain (list dox doy))
                                         (weights (mapcar (lambda (x l)
                                                            (cl-mpm/shape-function::shape-gimp x l h))
                                                          dist domain))
                                         (weight (* (the double-float (first weights))
                                                    (the double-float (second weights))))

                                         #+cl-mpm-fbar (weights-fbar (mapcar (lambda (x l)
                                                                 (cl-mpm/shape-function::shape-gimp-fbar x l h))
                                                               dist domain))
                                         #+cl-mpm-fbar (weight-fbar (* (the double-float (first weights-fbar))
                                                         (the double-float (second weights-fbar))))
                                         #-cl-mpm-fbar (weight-fbar 0d0)
                                         )
                                    (declare (dynamic-extent domain dist))
                                    (declare ;(type double-float h)
                                     (double-float weight)
                                     (type list domain dist)
                                     (dynamic-extent domain dist))
                                    (when (< 0d0 weight)
                                      (let* ((node (cl-mpm/mesh:get-node mesh id))
                                             (grads (mapcar (lambda (d l w)
                                                              (* (cl-mpm/shape-function::shape-gimp-dsvp d l h)
                                                                 (the double-float w)))
                                                            dist domain (nreverse weights)))
                                             #+cl-mpm-fbar (grads-fbar (mapcar (lambda (d l w)
                                                              (* (cl-mpm/shape-function::shape-gimp-dsvp d l h)
                                                                 (the double-float w)))
                                                            dist domain (nreverse weights-fbar)))
                                             #-cl-mpm-fbar (grads-fbar nil)
                                             )
                                        (funcall func mesh mp node weight grads weight-fbar grads-fbar))))))))))))))
;;This is more consise but half as fast
;; (defun iterate-over-neighbours-shape-gimp (mesh mp func)
;;   (declare (type cl-mpm/mesh::mesh mesh)
;;            (cl-mpm/particle:particle mp)
;;            (function func))
;;   (progn
;;     (with-accessors ((pos-vec cl-mpm/particle:mp-position)
;;                      (d0 cl-mpm/particle::mp-domain-size))
;;         mp
;;       (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
;;              (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
;;              (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec))
;;              )
;;         (declare (dynamic-extent pos pos-index))
;;         (declare (type double-float h)
;;                  (type list pos pos-index))
;;         (loop for dx from -1 to 1
;;               do (loop for dy from -1 to 1
;;                        do (let* ((id (mapcar #'+ pos-index (list dx dy))))
;;                             (declare (dynamic-extent id))
;;                             (when (cl-mpm/mesh:in-bounds mesh id)
;;                               (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
;;                                      (domain (loop for x across (magicl::matrix/double-float-storage d0) collect (* 0.5d0 (the double-float x))))
;;                                      (weights (mapcar (lambda (x l)
;;                                                         (cl-mpm/shape-function::shape-gimp x l h))
;;                                                       dist domain))
;;                                      (weight (reduce #'* weights))
;;                                      )
;;                                 (declare (dynamic-extent domain))
;;                                 (declare (type double-float h)
;;                                          (type list domain))
;;                                 (when (< 0d0 weight)
;;                                   (let* ((node (cl-mpm/mesh:get-node mesh id))
;;                                          (grads (mapcar (lambda (d l w)
;;                                                           (* (cl-mpm/shape-function::shape-gimp-dsvp d l h) w))
;;                                                         dist domain (nreverse weights))))
;;                                     (funcall func mesh mp node weight grads))))))))))))


(defun make-knot-list (mesh pos)
  (list
   ;; (loop for dy from -1 to 1
   ;;       collect)
   (let ((dy 0))
     (loop for dx from -4 to 4
           collect (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos (list dx dy)))))
   ;; (loop for dx from -1 to 1
   ;;       collect)
   (let ((dx 0))
     (loop for dy from -4 to 4
           collect (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos (list dx dy)))))))
(declaim (inline iterate-over-neighbours-shape-bspline))
(defun iterate-over-neighbours-shape-bspline (mesh mp func)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (pos-vec (cl-mpm/particle:mp-position mp))
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
           ;; (pos-index (magicl:scale pos-vec (/ 1 h)))
           (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec #'round))
           (border (not (and (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list 2 2)))
                             (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list -2 -2)))
                             (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list 2 -2)))
                             (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list -2 2)))
                             ))))
      (if border
          (loop for dx from -1 to 1
                do (loop for dy from -1 to 1
                         do (let* ((id (mapcar #'+ pos-index (list dx dy))))
                              (when (cl-mpm/mesh:in-bounds mesh id)
                                ;; Ill defined bspline
                                (when t
                                  (let* (
                                         (knot-list (make-knot-list mesh pos-index))
                                         (local-id-list (list (+ dx 0) (+ dy 0)))
                                         ;;Autogenerated bsplines have different requiremets dist measured from center
                                         (dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh pos-index)))
                                         (node (cl-mpm/mesh:get-node mesh id))
                                         (weights (mapcar (lambda (k x n)
                                                            (cl-mpm/shape-function::nodal-bspline k x n h))
                                                          knot-list dist local-id-list))
                                         (weight (reduce #'* weights))
                                         (grads (mapcar (lambda (k d n w)
                                                          (* (cl-mpm/shape-function::nodal-bspline-dsvp k d n h) w))
                                                        knot-list
                                                        dist
                                                        local-id-list
                                                        (nreverse weights)))
                                         )
                                    (funcall func mesh mp node weight grads nil nil)))))))
          (loop for dx from -1 to 1
                do (loop for dy from -1 to 1
                         do (let* ((id (mapcar #'+ pos-index (list dx dy))))
                              (when (cl-mpm/mesh:in-bounds mesh id)
                                (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                       (node (cl-mpm/mesh:get-node mesh id))
                                       (weights (mapcar (lambda (x) (cl-mpm/shape-function:shape-bspline x h)) dist))
                                       (weight (reduce #'* weights))
                                       (grads (mapcar (lambda (d w)
                                                        (* (cl-mpm/shape-function:shape-bspline-dsvp d h) w))
                                                      dist (nreverse weights))))
                                  (funcall func mesh mp node weight grads nil nil))))))))))

(defgeneric special-p2g (mp node svp dsvp)
  (:documentation "P2G behaviour for specific features")
  (:method (mp node svp dsvp)))

(defmethod special-p2g ((mp cl-mpm/particle::particle-thermal) node svp dsvp)
  (with-accessors ((node-mass  cl-mpm/mesh:node-mass)
                   (node-temp  cl-mpm/mesh:node-temperature)
                   (node-dtemp  cl-mpm/mesh:node-dtemp)
                   (node-volume  cl-mpm/mesh::node-volume)
                   (node-lock  cl-mpm/mesh:node-lock)) node
    (with-accessors ((mp-mass cl-mpm/particle:mp-mass)
                     (mp-temp cl-mpm/particle::mp-temperature)
                     (mp-heat-capaciy cl-mpm/particle::mp-heat-capacity)
                     (mp-thermal-conductivity cl-mpm/particle::mp-thermal-conductivity)
                     (mp-volume cl-mpm/particle:mp-volume)
                     ) mp
            (progn
              (let* ((weighted-temp (* mp-temp mp-mass svp))
                     (weighted-dtemp (* (/ mp-volume
                                           (* mp-mass
                                              mp-heat-capaciy))
                                        mp-thermal-conductivity
                                        mp-temp
                                        (magicl::sum
                                         (magicl:.* dsvp dsvp)))))
                  (sb-thread:with-mutex (node-lock)
                    (setf node-temp
                          (+ node-temp weighted-temp))
                    (setf node-dtemp
                          (+ node-dtemp weighted-dtemp))))))))


(declaim
 (inline p2g-mp)
 (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle) (values))
                p2g-mp))
(defun p2g-mp (mesh mp)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (with-accessors ((mp-vel  cl-mpm/particle:mp-velocity)
                   (mp-mass cl-mpm/particle:mp-mass)
                   (mp-volume cl-mpm/particle:mp-volume)
                   (mp-pmod cl-mpm/particle::mp-p-modulus)
                   (mp-damage cl-mpm/particle::mp-damage)
                   ;; (strain-rate cl-mpm/particle:mp-strain-rate)
                   ) mp
    (declare (type double-float mp-mass mp-volume))
    (iterate-over-neighbours
     mesh mp
     (lambda (mesh mp node svp grads fsvp fgrads)
       (declare
        (cl-mpm/particle:particle mp)
        (cl-mpm/mesh::node node)
        (double-float svp)
        )
       (with-accessors ((node-vel   cl-mpm/mesh:node-velocity)
                        (node-active  cl-mpm/mesh::node-active)
                        (node-mass  cl-mpm/mesh:node-mass)
                        (node-volume  cl-mpm/mesh::node-volume)
                        (node-volume-true  cl-mpm/mesh::node-volume-true)
                        (node-svp-sum  cl-mpm/mesh::node-svp-sum)
                        (node-force cl-mpm/mesh:node-force)
                        (node-p-wave cl-mpm/mesh::node-pwave)
                        (node-damage cl-mpm/mesh::node-damage)
                        (node-lock  cl-mpm/mesh:node-lock)) node
         (declare (type double-float node-mass node-volume mp-volume)
                  (type sb-thread:mutex node-lock))
         (sb-thread:with-mutex (node-lock)
           (setf node-active t)
           (incf node-mass
                 (* mp-mass svp))
           (incf node-volume
                 (* mp-volume svp))
           (incf node-p-wave
                 (* mp-pmod svp))
           (incf node-damage
                 (* mp-damage svp))
           (incf node-svp-sum svp)
           (fast-fmacc node-vel mp-vel (* mp-mass svp))
           )
         ;; (special-p2g mp node svp dsvp)
          )
        ;(p2g-mp-node mp node svp grads)
        )))
  (values))

(declaim (notinline p2g))
(defun p2g (mesh mps)
  (declare (type (array cl-mpm/particle:particle) mps) (cl-mpm/mesh::mesh mesh))
  (lparallel:pdotimes (i (length mps))
    (p2g-mp mesh (aref mps i))))

(declaim (notinline p2g-force-mp)
         (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle) (values)) p2g-force-mp)
         )
(defun p2g-force-mp (mesh mp)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (with-accessors ((mp-vel  cl-mpm/particle:mp-velocity)
                   (mp-mass cl-mpm/particle:mp-mass)
                   ) mp
    (declare (type double-float mp-mass))
    (let ((dsvp (cl-mpm/utils::dsvp-2d-zeros)))
      (iterate-over-neighbours
       mesh mp
       (lambda (mesh mp node svp grads fsvp fgrads)
         (declare
          (cl-mpm/particle:particle mp)
          (cl-mpm/mesh::node node)
          (double-float svp)
          )
         (with-accessors ((node-vel   cl-mpm/mesh:node-velocity)
                          (node-active  cl-mpm/mesh:node-active)
                          (node-mass  cl-mpm/mesh:node-mass)
                          (node-force cl-mpm/mesh:node-force)
                          (node-lock  cl-mpm/mesh:node-lock)) node
           (declare (type double-float node-mass)
                    (type sb-thread:mutex node-lock))
           (when node-active
             (sb-thread:with-mutex (node-lock)
               (det-ext-force mp node svp node-force)
               (cl-mpm/shape-function::assemble-dsvp-2d-prealloc grads dsvp)
               (det-int-force mp dsvp node-force)
               ;; (det-int-force mp (cl-mpm/shape-function::assemble-dsvp-2d grads) node-force)
               ))
           )
         ))))
  (values))

(declaim (notinline p2g-force))
(defun p2g-force (mesh mps)
  (declare (type (array cl-mpm/particle:particle) mps) (cl-mpm/mesh::mesh mesh))
  (lparallel:pdotimes (i (length mps)) 
    (p2g-force-mp mesh (aref mps i))))

(defgeneric special-g2p (mesh mp node svp grads)
  (:documentation "G2P behaviour for specific features")
  (:method (mesh mp node svp grads)))

(defmethod special-g2p (mesh (mp cl-mpm/particle::particle-thermal) node svp grads)
  (declare (ignore mesh))
  "Map grid to particle for one mp-node pair"
  (with-accessors ((node-temp cl-mpm/mesh:node-temperature)) node
    (with-accessors ((temp cl-mpm/particle::mp-temperature)) mp
      (progn
        (setf temp (+ temp (* node-temp svp)))
        ))))

(defmethod special-g2p (mesh (mp cl-mpm/particle::particle-damage) node svp grads)
  (declare (ignore mesh))
  "Map grid to particle for one mp-node pair"
  (with-accessors ((node-damage cl-mpm/mesh::node-damage)) node
    (with-accessors ((temp cl-mpm/particle::mp-damage)) mp
      (progn
        (setf temp (+ temp (* node-temp svp)))
        ))))


(defgeneric reset-mps-g2p (mp)
  (:method (mp)))
(defmethod reset-mps-g2p ((mp cl-mpm/particle::particle-thermal))
  (with-accessors ((temp cl-mpm/particle::mp-temperature)) mp
      (setf temp 0d0)))
(declaim
 (inline g2p-mp)
 (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle double-float) (values))
                g2p-mp))

(defun g2p-mp (mesh mp dt)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp)
           (double-float dt))
  "Map one MP from the grid"
  (with-accessors ((mass mp-mass)
                   (vel mp-velocity)
                   (pos mp-position)
                   (disp cl-mpm/particle::mp-displacement)
                   (acc cl-mpm/particle::mp-acceleration)
                   (temp cl-mpm/particle::mp-boundary)
                   (strain-rate mp-strain-rate)
                   (vorticity cl-mpm/particle:mp-vorticity)
                   (nc cl-mpm/particle::mp-cached-nodes)
                   )
    mp
    (let* ((mapped-vel (make-array 2 :initial-element 0d0 :element-type 'double-float))
          (mapped-acc (make-array 2 :initial-element 0d0 :element-type 'double-float))
           )
      (progn
        ;; (reset-mps-g2p mp)
        (setf temp 0d0)
        )
      ;; Map variables
      (iterate-over-neighbours
       mesh mp
       (lambda (mesh mp node svp grads fsvp fgrads)
         (declare
          (cl-mpm/mesh::node node)
          (cl-mpm/particle:particle mp)
          (type double-float svp))
         (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
                          (node-acc cl-mpm/mesh:node-acceleration)
                          (node-scalar cl-mpm/mesh::node-boundary-scalar)
                          (node-active cl-mpm/mesh:node-active)
                          ) node
           (when node-active
               (cl-mpm/fastmath::simd-fmacc mapped-vel (magicl::matrix/double-float-storage node-vel) svp)
               (cl-mpm/fastmath::simd-fmacc mapped-acc (magicl::matrix/double-float-storage node-acc) svp)
               (incf temp (* svp node-scalar))
               #+cl-mpm-special (special-g2p mesh mp node svp grads)
             )
           )
         ;; (g2p-mp-node mp node svp grads)
         ))
      ;;Update particle
      (progn
        (setf (fill-pointer nc) 0)
        (cl-mpm/fastmath::fast-fmacc-array (magicl::matrix/double-float-storage pos) mapped-vel dt)
        (cl-mpm/fastmath::fast-fmacc-array (magicl::matrix/double-float-storage disp) mapped-vel dt)
        (aops:copy-into (magicl::matrix/double-float-storage acc) mapped-acc)
        ;;FLIP
        #-cl-mpm-pic (cl-mpm/fastmath::simd-fmacc (magicl::matrix/double-float-storage vel)  mapped-acc dt)
        ;;PIC
        #+cl-mpm-pic (aops:copy-into (magicl::matrix/double-float-storage vel) mapped-vel)

        ;;Direct velocity damping
        ;; (magicl:scale! vel (- 1d0 1d-3))
        ))))

(declaim (notinline g2p))
(defun g2p (mesh mps dt)
  (declare (cl-mpm/mesh::mesh mesh) (array mps))
  "Map grid values to all particles"
  (lparallel:pdotimes (i (length mps))
    (g2p-mp mesh (aref mps i) dt)))

(defgeneric special-update-node (mesh dt node damping)
  (:documentation "Update node method")
  (:method (mesh dt node damping)))
(defmethod special-update-node (mesh dt (node cl-mpm/mesh::node-thermal) damping)
  (with-accessors ((mass  node-mass)
                   (temp node-temperature)
                   (dtemp node-dtemp))
      node
    (progn
      (setf temp (+ (/ temp mass) (* dtemp dt)))
      )))

(defun calculate-kinematics (node)
  "Calculate velocity from momentum on a single nod"
  (when (cl-mpm/mesh:node-active node)
      (with-accessors ( (mass  node-mass)
                        (vel   node-velocity))
        node
        (declare (double-float mass))
        (progn
          (magicl:scale! vel (/ 1.0d0 mass))
          ))))

(declaim (notinline calculate-forces)
         (ftype (function (cl-mpm/mesh::node double-float double-float double-float) (vaules)) calculate-forces))
(defun calculate-forces (node damping dt mass-scale)
  (when (cl-mpm/mesh:node-active node)
      (with-accessors ( (mass  node-mass)
                        (vel   node-velocity)
                        (force node-force)
                        (acc   node-acceleration))
        node
        (declare (double-float mass dt damping))
        (progn
          (magicl:scale acc 0d0)
          (cl-mpm/fastmath:fast-fmacc acc force (/ 1d0 (* mass mass-scale)))
          (cl-mpm/fastmath:fast-fmacc acc vel (/ (* damping -1d0) mass-scale))
          ;; (cl-mpm/fastmath:fast-fmacc acc vel (* damping -1d0))
          (cl-mpm/fastmath:fast-fmacc vel acc dt)
          )))
  (values))

(declaim (inline iterate-over-nodes)
         (ftype (function (cl-mpm/mesh::mesh function) (values)) iterate-over-nodes)
         )
(defun iterate-over-nodes (mesh func)
  (declare (type function func))
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
    (declare (type (array cl-mpm/particle:particle) nodes))
    (lparallel:pdotimes (i (array-total-size nodes))
      (let ((node (row-major-aref nodes i)))
        (funcall func node))))
  (values))

(declaim (inline iterate-over-nodes-serial)
         (ftype (function (cl-mpm/mesh::mesh function) (values)) iterate-over-nodes-serial)
         )
(defun iterate-over-nodes-serial (mesh func)
  (declare (type function func))
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
    (declare (type (array cl-mpm/particle:particle) nodes))
    (dotimes (i (array-total-size nodes))
      (let ((node (row-major-aref nodes i)))
        (funcall func node))))
  (values))

(defun update-node (mesh dt node damping)
    (when (cl-mpm/mesh:node-active node)
      (with-accessors ( (mass  node-mass)
                        (vel   node-velocity)
                        (force node-force)
                        (temp node-temperature)
                        (dtemp node-dtemp)
                        (acc   node-acceleration))
        node
        (progn
          (magicl:scale! vel (/ 1.0 mass))
          (setf acc (magicl:scale
                     (magicl:.- force
                                (magicl:scale vel damping))
                                (/ 1.0 mass)))
          (setf vel (magicl:.+ vel (magicl:scale acc dt)))
          (special-update-node mesh dt node damping)
          ))))

(defun update-nodes (mesh dt damping)
  "Explicitly integrate nodes"
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
  (lparallel:pdotimes (i (array-total-size nodes))
    (let ((node (row-major-aref nodes i)))
      (update-node mesh dt node damping)))))

(defun update-node-kinematics (mesh dt)
  (iterate-over-nodes mesh
                      (lambda (node)
                        (calculate-kinematics node))))
(defun update-node-forces (mesh damping dt mass-scale)
  (iterate-over-nodes mesh
                      (lambda (node)
                        (calculate-forces node damping dt mass-scale))))

(defun apply-bcs-old (mesh bcs dt)
  (declare (cl-mpm/mesh::mesh mesh))
  (with-accessors ( (nodes  mesh-nodes)
                    (nD     mesh-nD)
                    (mc     mesh-count)) mesh
    ;each bc is a list (pos value)
    (lparallel:pdotimes (i (length bcs))
      (let ((bc (nth i bcs)))
        (when bc
          (let ((index (cl-mpm/bc:bc-index bc)))
            (when (cl-mpm/mesh:in-bounds mesh index);Potentially throw here
              (cl-mpm/bc:apply-bc bc (cl-mpm/mesh:get-node mesh index) mesh dt))))))))

(defun apply-bcs (mesh bcs dt)
  (declare (cl-mpm/mesh::mesh mesh))
  (with-accessors ( (nodes  mesh-nodes)
                   (nD     mesh-nD)
                   (mc     mesh-count)) mesh
                                        ;each bc is a list (pos value)
    (lparallel:pdotimes (i (array-total-size bcs))
      (let ((bc (aref bcs i)))
        (with-accessors ((node cl-mpm/bc::bc-node)
                         (index cl-mpm/bc::bc-index))
            bc
          (if node
              (cl-mpm/bc:apply-bc bc node mesh dt)
              (progn
                (setf node (cl-mpm/mesh:get-node mesh index))
                (cl-mpm/bc:apply-bc bc node mesh dt)
                )))))))


;; (defun update-particle (mps dt)
;;   (declare (array mps))
;;   (loop for mp across mps
;;     do (with-accessors ((vel cl-mpm/particle:mp-velocity)
;;                         (pos cl-mpm/particle:mp-position)) mp
;;       (progn 
;;       (setf pos (magicl:.+ pos (magicl:scale vel dt)))))))
(declaim (inline update-particle-mp))
(defun update-particle-mp (mp dt)
  (declare (cl-mpm/particle:particle mp) (double-float dt))
  (with-accessors ((vel cl-mpm/particle:mp-velocity)
                   (pos cl-mpm/particle:mp-position)
                   (disp cl-mpm/particle::mp-displacement)
                   (nc cl-mpm/particle::mp-cached-nodes)
                   ) mp
      (progn
        (setf pos (magicl:.+ pos (magicl:scale vel dt)))
        (setf disp (magicl:.+ disp (magicl:scale vel dt)))
        (setf (fill-pointer nc) 0)
        ;; (.+ disp (magicl:scale vel dt) disp)
        )))

;; (defun update-particle (mps dt)
;;   (declare (array mps) (double-float dt))
;;   (lparallel:pdotimes (i (length mps))
;;      (update-particle-mp (aref mps i) dt)))

;Could include this in p2g but idk
(declaim (inline calculate-strain-rate)
         (ftype (function (cl-mpm/mesh::mesh  cl-mpm/particle:particle double-float)) calculate-strain-rate))
(defun calculate-strain-rate (mesh mp dt)
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt))
  (with-accessors ((strain-rate cl-mpm/particle:mp-strain-rate)
                   (vorticity cl-mpm/particle:mp-vorticity)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (jfbar cl-mpm/particle::mp-j-fbar)
                   (velocity-rate cl-mpm/particle::mp-velocity-rate)
                   ) mp
        (progn
          (magicl:scale! strain-rate 0d0)
          (magicl:scale! vorticity 0d0)
          (magicl:scale! stretch-tensor 0d0)
          (let (;(stretch-dsvp (magicl:zeros '(4 2)))
                #+cl-mpm-fbar (stretch-tensor-fbar (magicl:zeros '(2 2) :type 'double-float))
                ;(v-s (magicl:zeros '(2 2)))
                (v-s (matrix-zeros))
                (stretch-dsvp (stretch-dsvp-zeros))
                )
            (declare (magicl:matrix/double-float stretch-dsvp v-s)
                     (dynamic-extent stretch-dsvp v-s))
            (iterate-over-neighbours
             mesh mp
             (lambda (mesh mp node svp grads fsvp fgrads)
               (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
                                (node-active cl-mpm/mesh:node-active))
                   node
                 (declare (double-float))
                 (when node-active
                   (magicl:.+
                    stretch-tensor
                    (cl-mpm/utils::voight-to-stretch-prealloc
                     (magicl:@ (cl-mpm/shape-function::assemble-dstretch-2d-prealloc
                                grads stretch-dsvp) node-vel) v-s)
                    stretch-tensor)
                   #+cl-mpm-fbar (magicl:.+
                                  stretch-tensor-fbar
                                  (cl-mpm/utils::voight-to-stretch-prealloc
                                   (magicl:@ (cl-mpm/shape-function::assemble-dstretch-2d-prealloc
                                              fgrads stretch-dsvp) node-vel) v-s)
                                  stretch-tensor-fbar)
                   ;; (mult (cl-mpm/shape-function::assemble-dsvp-2d grads) node-vel strain-rate)
                   ;; (mult (cl-mpm/shape-function::assemble-vorticity-2d grads) node-vel vorticity)
                   ))))
            #+cl-mpm-fbar (setf jfbar (magicl:det (magicl:.+ (magicl:eye 2) stretch-tensor-fbar)))
            )
          #+cl-mpm-fbar (when (<= jfbar 0d0)
            (error "FBAR volume non-positive"))
            (cl-mpm/fastmath::stretch-to-sym stretch-tensor strain-rate)
            (cl-mpm/fastmath::stretch-to-skew stretch-tensor vorticity)
            (aops:copy-into (magicl::matrix/double-float-storage velocity-rate) (magicl::matrix/double-float-storage strain-rate))
            ;; (setf velocity-rate (magicl:scale strain-rate 1d0))
            (magicl:scale! stretch-tensor dt)
            (magicl:scale! strain-rate dt)
            (magicl:scale! vorticity dt)
            )))


(defun rotation-matrix (degrees)
  (let ((angle (/ (* pi degrees) 180)))
    (magicl:from-list (list (cos angle) (- (sin angle))
                            (sin angle) (cos angle))
                      '(2 2)
                      :type 'double-float
                      )))
(defun update-strain-linear (mesh mp dstrain)
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (volume-0 cl-mpm/particle::mp-volume-0)
                   (strain cl-mpm/particle:mp-strain)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   ) mp
    (let ((df (calculate-df mesh mp)))
                   (progn
                     (setf def (magicl:@ df def))
                     (setf strain (magicl:.+ strain dstrain))
                     ;(setf volume (* volume (det df)))
                     (setf volume (* volume-0 (magicl:det def)))

                     (multiple-value-bind (l v) (magicl:eig (magicl:@ def (magicl:transpose def)))
                       (let ((stretch
                               (magicl:@
                                v
                                (magicl:from-diag (mapcar (lambda (x) (the double-float (sqrt
                                                                                         (the double-float x)))) l) :type 'double-float)
                                (magicl:transpose v))))
                         (declare (type magicl:matrix/double-float stretch))
                         (setf (tref domain 0 0) (* (the double-float (tref domain 0 0))
                                                    (the double-float (tref stretch 0 0))))
                         (setf (tref domain 1 0) (* (the double-float (tref domain 1 0))
                                                    (the double-float (tref stretch 1 1))))
                         ;; (setf (tref domain 0 0) (* (the double-float (tref domain-0 0 0))
                         ;;                            (the double-float (tref stretch 0 0))))
                         ;; (setf (tref domain 1 0) (* (the double-float (tref domain-0 1 0))
                         ;;                            (the double-float (tref stretch 1 1))))
                         ))
                     ))))

(declaim (inline update-strain-kirchoff))
(declaim (ftype (function (cl-mpm/mesh::mesh
                           cl-mpm/particle:particle
                           double-float) (values))
               update-strain-kirchoff))
(defun update-strain-kirchoff (mesh mp dt)
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (volume-0 cl-mpm/particle::mp-volume-0)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   (strain-rate-tensor cl-mpm/particle::mp-strain-rate-tensor)
                   (velocity-rate cl-mpm/particle::mp-velocity-rate)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0)
                   (eng-strain-rate cl-mpm/particle::mp-eng-strain-rate)
                   ) mp
    (declare (type double-float volume)
             (type magicl:matrix/double-float
                   domain
                   ))
    (progn
      (let ((df (calculate-df mesh mp)))
        (progn
          ;; (magicl:mult df def :target def)
          (setf def (magicl:@ df def))
          (let ((initial-strain (magicl::copy-matrix/double-float strain))
                ;(initial-strain (cl-mpm/utils::matrix-zeros))
                (temp-strain-mat-a (cl-mpm/utils::matrix-zeros))
                (temp-strain-mat-b (cl-mpm/utils::matrix-zeros))
                )
            (multiple-value-bind (l v) (magicl:eig (voigt-to-matrix strain))
              (let (;(trail-lgs (cl-mpm/utils::matrix-zeros))
                    (trial-lgs (magicl:@ df
                                         v
                                         (cl-mpm/utils::matrix-from-list
                                          (list
                                           (the double-float (exp (* 2d0 (the double-float (nth 0 l)))))
                                           0d0 0d0
                                           (the double-float (exp (* 2d0 (the double-float (nth 1 l)))))))
                                         (magicl:transpose v)
                                         (magicl:transpose df)))
                    )
                (multiple-value-bind (lf vf) (magicl:eig
                                              trial-lgs)
                  (setf strain (magicl:scale!
                                (matrix-to-voigt
                                 (magicl:@
                                  vf
                                  (cl-mpm/utils::matrix-from-list
                                   (list
                                    (the double-float (log (the double-float (nth 0 lf))))
                                    0d0 0d0
                                    (the double-float (log (the double-float (nth 1 lf))))))
                                  (magicl:transpose vf)))
                                0.5d0))
                  ;; (print strain)
                  ;; (setf (magicl:tref strain 2 0) (* 2d0 (the double-float (magicl:tref strain 2 0))))
                  )))
            (magicl:.- initial-strain strain initial-strain)
            (setf eng-strain-rate initial-strain)
            (magicl:scale! eng-strain-rate (/ 1d0 dt))

            ;; (setf eng-strain-rate (magicl:scale! (magicl:.- strain initial-strain) (/ 1d0 dt)))
            )

          ;;Post multiply to turn to eng strain
          ;; (setf volume (* volume (magicl:det df)))
          (setf volume (* volume-0 (magicl:det def)))
          (when (<= volume 0d0)
            (error "Negative volume"))
          ;;Stretch rate update
          (let ((F (cl-mpm/utils::matrix-zeros)))
            (magicl:mult def def :target F :transb :t)
            (multiple-value-bind (l v) (magicl:eig F)
              (let ((stretch
                      (magicl:@
                       v
                       (cl-mpm/utils::matrix-from-list
                        (list (the double-float (sqrt (the double-float (nth 0 l))))
                              0d0 0d0
                              (the double-float (sqrt (the double-float (nth 1 l))))))
                       (magicl:transpose v)))
                    )
                (declare (type magicl:matrix/double-float stretch))
                (setf (tref domain 0 0) (* (the double-float (tref domain-0 0 0))
                                           (the double-float (tref stretch 0 0))))
                (setf (tref domain 1 0) (* (the double-float (tref domain-0 1 0))
                                           (the double-float (tref stretch 1 1))))
                )))
          ;; (update-domain-corner mesh mp dt)
          ;; (setf (tref domain 0 0) (* (the double-float (tref domain 0 0))
          ;;                             (the double-float (tref df 0 0))))
          ;; (setf (tref domain 1 0) (* (the double-float (tref domain 1 0))
          ;;                             (the double-float (tref df 1 1))))
          )
          )))
  (values))
(defun update-domain-corner (mesh mp dt)
  (with-accessors ((position cl-mpm/particle::mp-position)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0)
                   )
      mp
    (let ((points (list '(-0.5d0 -0.5d0) '(0.5d0 -0.5d0) '(0.5d0 0.5d0) '(-0.5d0 0.5d0)))
          (dim '(0 0 1 1))
          ;(disp (make-array 4 :element-type 'double-float :initial-element 0d0))
          (disp (make-array 4 :initial-element (magicl:zeros '(2 1) :type 'double-float)))
          )
      (loop for point in points
            for d in dim
            for i from 0 below 4
            do
               (let ((corner (magicl:.+ position
                                        (magicl:.*
                                         (magicl:from-list point '(2 1)
                                                           :type 'double-float)
                                         domain
                                         ))))
                 (iterate-over-neighbours-point-linear
                  mesh corner
                  (lambda (mesh node svp grads)
                    (with-accessors ((vel cl-mpm/mesh:node-velocity))
                        node
                      ;; (incf (aref disp i)
                      ;;       (* svp dt (magicl:tref vel d 0)))
                      (setf (aref disp i)
                            (magicl:.+ (aref disp i)
                                       (magicl:scale vel (* dt svp)))))))))
      (incf (magicl:tref domain 0 0) (* 0.5
                                        (-
                                         (+ (magicl:tref (aref disp 1) 0 0)
                                            (magicl:tref (aref disp 2) 0 0))
                                         (+ (magicl:tref (aref disp 0) 0 0)
                                            (magicl:tref (aref disp 3) 0 0))
                                         )))
      (incf (magicl:tref domain 1 0) (* 0.5
                                        (-
                                         (+ (magicl:tref (aref disp 3) 1 0)
                                            (magicl:tref (aref disp 2) 1 0))
                                         (+ (magicl:tref (aref disp 0) 1 0)
                                            (magicl:tref (aref disp 1) 1 0))
                                         )))
      (let* ((jf (magicl:det def))
             (jl (* (magicl:tref domain 0 0) (magicl:tref domain 1 0)))
             (jl0 (* (magicl:tref domain-0 0 0) (magicl:tref domain-0 1 0)))
             (scaling (expt (/ (* jf jl0) jl) 0.5d0)))
        (setf (magicl:tref domain 0 0) (* (magicl:tref domain 0 0) scaling)
              (magicl:tref domain 1 0) (* (magicl:tref domain 1 0) scaling)
              ))
      )

    )
  )



(declaim (inline update-stress-mp)
         (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle double-float) (values)) update-stress-mp))
(defun update-stress-mp (mesh mp dt)
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt))
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (stress-kirchoff cl-mpm/particle::mp-stress-kirchoff)
                     (volume cl-mpm/particle:mp-volume)
                     (strain cl-mpm/particle:mp-strain)
                     (def    cl-mpm/particle:mp-deformation-gradient)
                     (strain-rate cl-mpm/particle:mp-strain-rate)
                        ) mp
      (declare (magicl:matrix/double-float stress stress-kirchoff strain def strain-rate))
      (progn
      ;;   ;;For no FBAR we need to update our strains
          (progn
            (calculate-strain-rate mesh mp dt)

            ;; Turn cauchy stress to kirchoff
            (setf stress stress-kirchoff)

            ;; Update our strains
            (update-strain-kirchoff mesh mp dt)

            ;;Update our kirchoff stress with constitutive model
            (setf stress-kirchoff (cl-mpm/particle:constitutive-model mp strain dt))

            ;; Check volume constraint!
            (when (<= volume 0d0)
              (error "Negative volume"))
            ;; Turn kirchoff stress to cauchy
            (setf stress (magicl:scale stress-kirchoff (/ 1.0d0 (the double-float (magicl:det def)))))
            ))))
(defun calculate-cell-deformation (mesh cell dt)
  (with-accessors ((def cl-mpm/mesh::cell-deformation-gradient)
                   (volume cl-mpm/mesh::cell-volume))
      cell
    (let ((dstrain (magicl:zeros '(3 1))))
      (cl-mpm/mesh::cell-iterate-over-neighbours
       mesh cell
       (lambda (mesh c node svp grads)
         (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
                          (node-active cl-mpm/mesh:node-active)
                          )
             node
           (declare (double-float))
           (when node-active
             (mult (cl-mpm/shape-function::assemble-dsvp-2d grads) node-vel dstrain))
           )))
      (magicl:scale! dstrain dt)
      (let ((df (.+ (magicl:eye 2) (voight-to-matrix dstrain))))
        (progn
          (setf def (magicl:@ df def))
          (setf volume (* volume (det df)))
          )))))
(declaim (ftype (function (cl-mpm/mesh::mesh
                           (array cl-mpm/particle:particle)
                           double-float) (values)) update-stress))
(defun map-jacobian (mesh mp dt)
  (with-accessors ((dstrain cl-mpm/particle::mp-strain-rate)
                   (stretch-rate cl-mpm/particle::mp-stretch-tensor)
                   (volume cl-mpm/particle:mp-volume)
                   (def cl-mpm/particle:mp-deformation-gradient))
      mp
    (let* ((df (.+ (magicl:eye 2) stretch-rate))
           (j-inc (det df))
           (j-n (det def)))
      (iterate-over-neighbours mesh mp
                               (lambda (mesh mp node svp grads fsvp fgrads)
                                 (with-accessors ((node-active cl-mpm/mesh:node-active)
                                                  (node-j-inc cl-mpm/mesh::node-jacobian-inc)
                                                  (node-lock cl-mpm/mesh::node-lock)
                                                  )
                                     node
                                   (when node-active
                                     (sb-thread:with-mutex (node-lock)
                                       (incf node-j-inc (* svp volume j-inc j-n)))
                                     ))
                                 )))))

(declaim (inline calculate-df)
         (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle) magicl:matrix/double-float)
                calculate-df))
(defun calculate-df (mesh mp)
  (with-accessors ((dstrain cl-mpm/particle::mp-strain-rate)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (jfbar cl-mpm/particle::mp-j-fbar)
                   (def cl-mpm/particle:mp-deformation-gradient))
      mp
    (let* ((df (cl-mpm/utils::matrix-from-list '(1d0 0d0 0d0 1d0))))
      (.+ df stretch-tensor df)
      #+cl-mpm-fbar (progn
                      (magicl:scale! df (expt (/ jfbar (magicl:det df)) (/ 1d0 2d0))))
      df)))

(defun update-stress (mesh mps dt)
  (declare ((array cl-mpm/particle:particle) mps) (cl-mpm/mesh::mesh mesh))
  ;;Update stress
  (lparallel:pdotimes (i (length mps))
     (update-stress-mp mesh (aref mps i) dt))
  ;; (lparallel:pmap nil (lambda (mp) (update-stress-mp mesh mp dt)) mps)
  ;; (dotimes (i (length mps))
  ;;   (update-stress-mp mesh (aref mps i) dt))
  (values))

(declaim (notinline reset-grid))
(defun reset-grid (mesh)
  (declare (cl-mpm/mesh::mesh mesh))
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
    (lparallel:pdotimes (i (array-total-size nodes))
      (let ((node (row-major-aref nodes i)))
        (when (cl-mpm/mesh:node-active node)
          (cl-mpm/mesh:reset-node node))))))

(defun filter-grid (mesh mass-thresh)
  (declare (cl-mpm/mesh::mesh mesh) (double-float mass-thresh))
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
  (dotimes (i (array-total-size nodes))
    (let ((node (row-major-aref nodes i)))
      (when (and (cl-mpm/mesh:node-active node)
                 (< (cl-mpm/mesh:node-mass node) mass-thresh))
        (setf (cl-mpm/mesh::node-active node) nil)
        (cl-mpm/mesh:reset-node node)
        )))))

(defun remove-material-damaged (sim)
  "Remove material points that have strain energy density exceeding fracture toughness"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (let ((h (cl-mpm/mesh:mesh-resolution mesh)))
       (setf mps
             (lparallel:premove-if (lambda (mp)
                          (with-accessors ((damage cl-mpm/particle:mp-damage)
                                           (def cl-mpm/particle::mp-deformation-gradient))
                              mp
                            (and (>= damage 1d0)
                                 ;; (or
                                  (split-criteria mp h)
                                 ;;  (< (magicl:det def) 1d-3)
                                 ;;  )
                                 ))) mps)))
      ;(delete-if (lambda (mp)
      ;                        (with-accessors ((damage cl-mpm/particle:mp-damage))
      ;                            mp
      ;                          (and (>= damage 1d0)
      ;                                  ;(split-criteria mp h)
      ;                               ))) mps))
    ))
(defun split-criteria (mp h)
  (with-accessors ((def cl-mpm/particle:mp-deformation-gradient)
                   (lens cl-mpm/particle::mp-domain-size)
                   (lens-0 cl-mpm/particle::mp-domain-size-0)
                   )
      mp
    (let ((l-factor 1.50d0)
          (h-factor (* 0.8d0 h))
          (s-factor 1.5d0))
      (cond
        ;; ((< h-factor (tref lens 0 0)) :x)
        ;; ((< h-factor (tref lens 1 0)) :y)
        ((< (* l-factor (tref lens-0 0 0)) (tref lens 0 0)) :x)
        ((< (* l-factor (tref lens-0 1 0)) (tref lens 1 0)) :y)
        ;; ((< l-factor (/ (tref lens 0 0) (tref lens-0 0 0))) t)
        ;; ((< l-factor (/ (tref lens 1 0) (tref lens-0 1 0))) t)
                                        ;((< 2.0d0 (tref def 1 1)) t)
                                        ;((> 1.5d0 (tref def 0 0)) t)
        (t nil)
        ))))
(defun copy-particle (original &rest initargs &key &allow-other-keys)
  (let* ((class (class-of original))
         (copy (allocate-instance class)))
    (dolist (slot (mapcar #'sb-mop:slot-definition-name (sb-mop:class-slots class)))
      (when (slot-boundp original slot)
        (cond
          ((typep (slot-value original slot) 'magicl:tensor)
           (setf (slot-value copy slot)
                 (magicl:scale (slot-value original slot) 1d0)))
          (t (setf (slot-value copy slot)
                  (slot-value original slot))))))
    (apply #'reinitialize-instance copy initargs)))
(defun split-mp (mp h direction)
  (with-accessors ((def cl-mpm/particle:mp-deformation-gradient)
                   (lens cl-mpm/particle::mp-domain-size)
                   (lens-0 cl-mpm/particle::mp-domain-size-0)
                   (pos cl-mpm/particle:mp-position)
                   (mass cl-mpm/particle:mp-mass)
                   (volume cl-mpm/particle:mp-volume)
                   ;; (volume cl-mpm/particle::mp-volume-0)
                   )
      mp
    (let ((l-factor 1.20d0)
          (h-factor (* 0.8d0 h)))
    (cond
      ((eq direction :x)
       (let ((new-size (magicl:.* lens (magicl:from-list '(0.5d0 1d0) '(2 1))))
             (new-size-0 (magicl:.* lens-0 (magicl:from-list '(0.5d0 1d0) '(2 1))))
             (pos-offset (magicl:.* lens (magicl:from-list '(0.25d0 0d0) '(2 1)))))
         (list
          (copy-particle mp
                         :mass (/ mass 2)
                         :volume (/ volume 2)
                         :size new-size
                         :size-0 new-size-0
                         :position (magicl:.+ pos pos-offset))
          (copy-particle mp
                         :mass (/ mass 2)
                         :volume (/ volume 2)
                         :size new-size
                         :size-0 new-size-0
                         :position (magicl:.- pos pos-offset)))))
      ((eq direction :y)
       (let ((new-size (magicl:.* lens (magicl:from-list '(1d0 0.5d0) '(2 1))))
             (new-size-0 (magicl:.* lens-0 (magicl:from-list '(1d0 0.5d0) '(2 1))))
             (pos-offset (magicl:.* lens (magicl:from-list '(0d0 0.25d0) '(2 1)))))
         (list
          (copy-particle mp
                         :mass (/ mass 2)
                         :volume (/ volume 2)
                         :size new-size
                         :size-0 new-size-0
                         :position (magicl:.+ pos pos-offset))
          (copy-particle mp
                         :mass (/ mass 2)
                         :volume (/ volume 2)
                         :size new-size
                         :size-0 new-size-0
                         :position (magicl:.- pos pos-offset)))))
      ((eq direction :xy)
       (let ((new-size (magicl:.* lens (magicl:from-list '(0.5d0 0.5d0) '(2 1))))
             (new-size-0 (magicl:.* lens-0 (magicl:from-list '(0.5d0 0.5d0) '(2 1))))
             (pos-offset (magicl:.* lens (magicl:from-list '(0.25d0 0.25d0) '(2 1)))))
         (list
          (copy-particle mp
                         :mass (/ mass 4d0)
                         :volume (/ volume 4d0)
                         :size new-size
                         :size-0 new-size-0
                         :position (magicl:.+ pos (magicl:.* pos-offset (magicl:from-list '(1d0 1d0) '(2 1)))))
          (copy-particle mp
                         :mass (/ mass 4d0)
                         :volume (/ volume 4d0)
                         :size new-size
                         :size-0 new-size-0
                         :position (magicl:.+ pos (magicl:.* pos-offset (magicl:from-list '(-1d0 1d0) '(2 1)))))
          (copy-particle mp
                         :mass (/ mass 4d0)
                         :volume (/ volume 4d0)
                         :size new-size
                         :size-0 new-size-0
                         :position (magicl:.+ pos (magicl:.* pos-offset (magicl:from-list '(-1d0 -1d0) '(2 1)))))
          (copy-particle mp
                         :mass (/ mass 4d0)
                         :volume (/ volume 4d0)
                         :size new-size
                         :size-0 new-size-0
                         :position (magicl:.+ pos (magicl:.* pos-offset (magicl:from-list '(1d0 -1d0) '(2 1)))))
          )))
      (t nil)
      ))))
(defun split-mps (sim)
  "Remove material points that have strain energy density exceeding fracture toughness"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (mps-to-split (remove-if-not (lambda (mp) (split-criteria mp h)) mps))
           (split-direction (map 'list (lambda (mp) (split-criteria mp h)) mps-to-split)))
      (setf mps (delete-if (lambda (mp) (split-criteria mp h)) mps))
      (loop for mp across mps-to-split
            for direction in split-direction
            do (loop for new-mp in (split-mp mp h direction)
                     do (vector-push-extend new-mp mps)))
      )))

(defun split-mps-criteria (sim criteria)
  "Remove material points that have strain energy density exceeding fracture toughness"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (mps-to-split (remove-if-not (lambda (mp) (funcall criteria mp h)) mps))
           (split-direction (map 'list (lambda (mp) (funcall criteria mp h)) mps-to-split)))
      (print split-direction)
      (setf mps (delete-if (lambda (mp) (funcall criteria mp h)) mps))
      (loop for mp across mps-to-split
            for direction in split-direction
            do (loop for new-mp in (split-mp mp h direction)
                     do (vector-push-extend new-mp mps)))
      )))

(defun calculate-min-dt (sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mass-scale cl-mpm::sim-mass-scale))
      sim
    (let ((inner-factor most-positive-double-float))
      (iterate-over-nodes-serial
       mesh
       (lambda (node)
         (with-accessors ((node-active  cl-mpm/mesh:node-active)
                          (pmod cl-mpm/mesh::node-pwave)
                          (mass cl-mpm/mesh::node-mass)
                          (svp-sum cl-mpm/mesh::node-svp-sum)
                          (vol cl-mpm/mesh::node-volume)
                          ) node
           (when node-active
             (let ((nf (/ mass (* vol (/ pmod svp-sum)))))
                 (when (< nf inner-factor)
                   ;; (format t "Mass: ~a - Vol: ~a - Pmod: ~a~%" mass vol (/ pmod svp-sum))
                   (setf inner-factor nf)))))))
      (if (< inner-factor most-positive-double-float)
          (* (sqrt mass-scale) (sqrt inner-factor) (cl-mpm/mesh:mesh-resolution mesh))
          0d0))))
#||
(progn
(ql:quickload :cl-mpm/examples/fracture)
(in-package :cl-mpm/examples/fracture))
||#


