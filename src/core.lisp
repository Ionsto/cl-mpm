(defpackage :cl-mpm
  (:use :cl
        :cl-mpm/particle
        :cl-mpm/mesh
        :cl-mpm/utils
        :cl-mpm/fastmath
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
    #:update-sim
    #:make-shape-function-linear
    #:sim-mesh
    #:sim-mps
    #:sim-bcs
    #:sim-dt
    #:sim-damping-factor
    #:sim-mass-filter
    #:sim-mass-scale
    #:sim-allow-mp-split
    #:post-stress-step
    #:iterate-over-nodes
    #:iterate-over-nodes-serial
    #:iterate-over-neighbours
    #:calculate-adaptive-time
    #:iterate-over-mps
    #:iterate-over-neighbours-point-linear
    ))
;; (sb-ext:restrict-compiler-policy 'speed  3 3)
;; (declaim (optimize (debug 3) (safety 3) (speed 3)))
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
  (format t "Compiled with optimize ~A~%" (sb-ext:restrict-compiler-policy))
)

(defun FLIP-status ()
  #+cl-mpm-pic (print "Compiled with PIC")
  #-cl-mpm-pic (print "Compiled with FLIP")
  )


(defclass mpm-sim ()
  ((dt
     :accessor sim-dt
     :initarg :dt)
   (time
    :accessor sim-time
    :initform 0d0
    :initarg :time)
   (stats-mps-removed
    :accessor sim-stats-mps-removed
    :initform 0d0
    :initarg :time)
   (stats-mps-added
    :accessor sim-stats-mps-added
    :initform 0d0
    :initarg :time)
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
   (bcs-force-list
    :accessor sim-bcs-force-list
    :initarg :bcs-force-list
    :initform nil)
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
    :initform nil)
   (enable-damage
    :type boolean
    :accessor sim-enable-damage
    :initarg :enable-damage
    :initform nil)
   (enable-fbar
    :type boolean
    :accessor sim-enable-fbar
    :initarg :enable-fbar
    :initform nil)
   (nonlocal-damage
    :type boolean
    :accessor sim-nonlocal-damage
    :initform t)
   (allow-mp-damage-removal
    :type boolean
    :accessor sim-allow-mp-damage-removal
    :initarg :damage-removal
    :initform nil)
   (mp-damage-removal-instant
    :type boolean
    :accessor sim-mp-damage-removal-instant
    :initform nil)
   (mass-filter
     :type double-float
     :accessor sim-mass-filter
     :initarg :mass-filter
     :initform 1d-15))
  (:documentation "A self contained mpm simulation"))
(defclass mpm-nd-2d ()())
(defclass mpm-nd-3d ()())
(defclass mpm-sf-mpm ()())
(defclass mpm-sf-gimp ()())
(defclass mpm-sim-quasi-static (mpm-sim) ())

(defclass mpm-sim-usf (mpm-sim)
  ()
  (:documentation "Explicit simulation with update stress first update"))
(defclass mpm-sim-usl (mpm-sim)
  ()
  (:documentation "Explicit simulation with update stress last update"))


(defun make-mpm-sim (size resolution dt shape-function &key (sim-type 'mpm-sim-usf))
  "Constructs an mp with critical infomation like mesh and number of dimentions"
  (make-instance sim-type
                 :dt (coerce dt 'double-float)
                 :mesh (make-mesh size resolution shape-function)
                 :mps '()))

;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
(defun check-mps (sim)
  "Function to check that stresses and positions are sane, deleting mps moving very fast"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
  (with-accessors ((h cl-mpm/mesh::mesh-resolution))
      mesh
    ;; loop for mp across mps do
    (cl-mpm::iterate-over-mps
     mps
     (lambda (mp)
       (progn
         (with-accessors ((pos cl-mpm/particle::mp-position)
                          (vel cl-mpm/particle::mp-velocity)
                          (domain cl-mpm/particle::mp-domain-size)
                          ) mp
           (loop for i from 0 to 2
                 do (progn
                      (when (sb-ext:float-nan-p (magicl:tref pos i 0)) 
                        (pprint mp)
                        (error "NaN location found for ~A" mp))
                      (when (equal (abs (magicl:tref vel i 0)) #.sb-ext:double-float-positive-infinity)
                        (pprint mp)
                        (error "Infinite velocity found"))
                      (when (> (abs (magicl:tref vel i 0)) 1e10)
                        (pprint mp)
                        (error "High velocity found"))))))))
    (remove-mps-func
       sim
       (lambda (mp)
         (with-accessors ((damage cl-mpm/particle:mp-damage)
                          (def cl-mpm/particle::mp-deformation-gradient))
             mp
           (or
            (gimp-removal-criteria mp h)
            ;; (cl-mpm/mesh::in-bounds-array mesh (magicl::matrix/double-float-storage
            ;;                                     (cl-mpm/particle::mp-position mp)))
            ))))
    )))
;; (sb-ext:restrict-compiler-policy 'debug  0 0)
;; (sb-ext:restrict-compiler-policy 'speed  3 3)
(defun check-single-mps (sim)
  "Function to check and remove single material points"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
  (with-accessors ((h cl-mpm/mesh::mesh-resolution))
      mesh
    (iterate-over-mps
     mps
     (lambda (mp)
       (setf (cl-mpm/particle::mp-single-particle mp)
             (single-particle-criteria mesh mp))))
      (remove-mps-func
       sim
       #'cl-mpm/particle::mp-single-particle))))


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
               (fbar enable-fbar)
               (bcs-force-list bcs-force-list)
               (time time))
                sim
    (declare (type double-float mass-filter))
                (progn

                  (reset-grid mesh)
                  (p2g mesh mps)
                  (when (> mass-filter 0d0)
                    (filter-grid mesh (sim-mass-filter sim)))
                  (update-node-kinematics mesh dt )
                  (apply-bcs mesh bcs dt)
                  (update-stress mesh mps dt fbar)
                  ;; ;; Map forces onto nodes
                  (p2g-force mesh mps)
                  (loop for bcs-f in bcs-force-list
                        do (apply-bcs mesh bcs-f dt))
                  (update-node-forces sim)
                  ;; Reapply velocity BCs
                  (apply-bcs mesh bcs dt)
                  ;; Also updates mps inline
                  (g2p mesh mps dt)
                  (when split
                    (split-mps sim))
                  (check-mps sim)
                  (incf time dt))))
(defmethod update-sim ((sim mpm-sim-usf))
  "Update stress first algorithm"
  (declare (cl-mpm::mpm-sim-usf sim))
  (with-slots ((mesh mesh)
               (mps mps)
               (bcs bcs)
               (bcs-force bcs-force)
               (bcs-force-list bcs-force-list)
               (dt dt)
               (mass-filter mass-filter)
               (split allow-mp-split)
               (enable-damage enable-damage)
               (nonlocal-damage nonlocal-damage)
               (remove-damage allow-mp-damage-removal)
               (fbar enable-fbar)
               (time time)
               )
                sim
    (declare (double-float mass-filter dt time))
                (progn
                    (reset-grid mesh)
                    (p2g mesh mps)
                    (when (> mass-filter 0d0)
                      (filter-grid mesh (sim-mass-filter sim)))
                    (update-node-kinematics mesh dt )
                    (apply-bcs mesh bcs dt)
                    (update-stress mesh mps dt fbar)
                    ;; Map forces onto nodes
                    (p2g-force mesh mps)
                    (loop for bcs-f in bcs-force-list
                          do (apply-bcs mesh bcs-f dt))
                    (update-node-forces sim)
                    ;; Reapply velocity BCs
                    (apply-bcs mesh bcs dt)
                    ;; Also updates mps inline
                    (g2p mesh mps dt)
                    (when split
                      (split-mps sim))
                    (check-mps sim)
                    (incf time dt))))

(defmethod update-sim ((sim mpm-sim-usl))
  "Update stress last algorithm"
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
               (fbar enable-fbar)
               (bcs-force-list bcs-force-list)
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
                    (loop for bcs-f in bcs-force-list
                          do (apply-bcs mesh bcs-f dt))
                    ;; (apply-bcs mesh bcs-force dt)
                    ;;Update our nodes after force mapping
                    (update-node-forces sim)
                    ;;Apply velocity bcs
                    (apply-bcs mesh bcs dt)
                    ;;Grid to particle mapping
                    (g2p mesh mps dt)

                    ;;2nd round of momentum mapping
                    (reset-grid-velocity mesh)
                    (p2g mesh mps)
                    (when (> mass-filter 0d0)
                     (filter-grid-velocity mesh (sim-mass-filter sim)))
                    (update-node-kinematics mesh dt)
                    (apply-bcs mesh bcs dt)
                    ;;Update stress last
                    (update-stress mesh mps dt fbar)

                    (when remove-damage
                      (remove-material-damaged sim))
                    (when split
                      (split-mps sim))
                    (check-mps sim)
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
  ;; (create-node-cache mesh mp func)
  (if (> (length (cl-mpm/particle::mp-cached-nodes mp)) 0)
      (iterate-over-neighbours-cached mesh mp func)
      (create-node-cache mesh mp func))
  ;;Ideally we should be dispatching over simulation shape function type but its faster to hard code
  ;; (iterate-over-neighbours-shape-gimp mesh mp func)
  ;; (iterate-over-neighbours-shape-linear mesh mp func)
  ;; (iterate-over-neighbours-shape mesh (cl-mpm/mesh:mesh-shape-func mesh) mp func)
  (values)
  )

(declaim (inline create-node-cache))
(defun create-node-cache (mesh mp func)
  "A function iterating over neighbours executing function, while also caching the relevent node, gradients and weights"
  (declare (function func))
  (with-accessors ((nodes cl-mpm/particle::mp-cached-nodes))
      mp
    ;;Simple if statement - we take the hit
    (if (= (the fixnum (cl-mpm/mesh:mesh-nd mesh)) 2)
        (iterate-over-neighbours-shape-gimp-simd;bspline
         ;; iterate-over-neighbours-shape-linear
         mesh mp
         (lambda (mesh mp node svp grads fsvp fgrads)
           (vector-push-extend
            (cl-mpm/particle::make-node-cache
             :node node
             :weight svp
             :grads grads
             :weight-fbar fsvp
             :grads-fbar fgrads)
            nodes)
           (funcall func mesh mp node svp grads fsvp fgrads)))
        (iterate-over-neighbours-shape-gimp-3d
         mesh mp
         (lambda (mesh mp node svp grads fsvp fgrads)
           (vector-push-extend
            (cl-mpm/particle::make-node-cache
             :node node
             :weight svp
             :grads grads
             :weight-fbar fsvp
             :grads-fbar fgrads)
            nodes)
           (funcall func mesh mp node svp grads fsvp fgrads))))))

(declaim (inline iterate-over-neighbours-cached))
(defun iterate-over-neighbours-cached (mesh mp func)
  "If a node iteration cache has been generated we loop over the data list"
  (declare (function func))
  (loop for nc across (the (vector cl-mpm/particle::node-cache *) (cl-mpm/particle::mp-cached-nodes mp))
        do
         (funcall func mesh mp
                  (cl-mpm/particle::node-cache-node nc)
                  (cl-mpm/particle::node-cache-weight nc)
                  (cl-mpm/particle::node-cache-grads nc)
                  (cl-mpm/particle::node-cache-weight-fbar nc)
                  (cl-mpm/particle::node-cache-grads-fbar nc)
                  )))

;;This is one method of dispatching over different types of shape functions
(defmethod iterate-over-neighbours-shape (mesh (shape-func cl-mpm/shape-function:shape-function-linear) mp func)
  (iterate-over-neighbours-shape-linear mesh mp func))

(defmethod iterate-over-neighbours-shape (mesh (shape-func cl-mpm/shape-function:shape-function-bspline) mp func)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/shape-function::shape-function-bspline shape-func))
  (iterate-over-neighbours-shape-bspline mesh mp func))

(declaim (inline iterate-over-neighbours-shape-linear))
(defun iterate-over-neighbours-shape-linear (mesh mp func)
  "The simplest shape function implementation - linear"
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (pos-vec (cl-mpm/particle:mp-position mp))
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
           (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec #'floor)))
      (loop for dx from 0 to 1
            do (loop for dy from 0 to 1
                     do (let* ((id (mapcar #'+ pos-index (list dx dy 0))))
                          (when (cl-mpm/mesh:in-bounds mesh id)
                            (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                   (node (cl-mpm/mesh:get-node mesh id))
                                   (weights (mapcar (lambda (x) (cl-mpm/shape-function::shape-linear x h)) dist))
                                   (weight (reduce #'* weights))
                                   (grads (mapcar (lambda (d w) (* (cl-mpm/shape-function::shape-linear-dsvp d h) w))
                                                  dist (nreverse weights))))
                              (when (< 0d0 weight)
                                (funcall func mesh mp node weight (append grads (list 0d0)) 0d0 (list 0d0 0d0 0d0)))))))))))

(declaim (inline iterate-over-neighbours-shape-linear))
(defun iterate-over-neighbours-shape-linear-3d (mesh mp func)
  "The simplest shape function implementation - linear"
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (pos-vec (cl-mpm/particle:mp-position mp))
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)  (tref pos-vec 2 0)))
           (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec #'floor)))
      (loop for dx from 0 to 1
            do (loop for dy from 0 to 1
                     do
                        (loop for dz from 0 to 1
                              do
                                 (let* ((id (mapcar #'+ pos-index (list dx dy dz))))
                                   (when (cl-mpm/mesh:in-bounds mesh id)
                                     (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                            (node (cl-mpm/mesh:get-node mesh id))
                                            (weights (mapcar (lambda (x) (cl-mpm/shape-function::shape-linear x h)) dist))
                                            (weight (reduce #'* weights))
                                            (grads
                                              (list
                                               (* (cl-mpm/shape-function::shape-linear-dsvp (nth 0 dist) h)
                                                  (nth 1 dist) (nth 2 dist))
                                               (* (cl-mpm/shape-function::shape-linear-dsvp (nth 1 dist) h)
                                                  (nth 0 dist) (nth 2 dist))
                                               (* (cl-mpm/shape-function::shape-linear-dsvp (nth 2 dist) h)
                                                  (nth 1 dist) (nth 0 dist))
                                               )))
                                       (when (< 0d0 weight)
                                         (funcall func mesh mp node weight grads 0d0 nil)))))))))))


(declaim (inline iterate-over-neighbours-point-linear)
         (ftype (function (cl-mpm/mesh::mesh magicl:matrix/double-float function) (values))
                iterate-over-neighbours-point-linear))

(defun iterate-over-neighbours-point-linear (mesh position func)
  "Iterate over neighbours of an arbitrary point - using FEM linear basis"
  (if (= (the fixnum (cl-mpm/mesh:mesh-nd mesh)) 2)
      (iterate-over-neighbours-point-linear-simd mesh position func)
      (iterate-over-neighbours-point-linear-3d mesh position func)))

(declaim (inline iterate-over-neighbours-point-linear-lisp)
         (ftype (function (cl-mpm/mesh::mesh magicl:matrix/double-float function) (values))
                iterate-over-neighbours-point-linear-lisp))

(defun iterate-over-neighbours-point-linear-lisp (mesh position func)
  "Iterating over basis functions in 2D"
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
(defun iterate-over-neighbours-point-linear-3d (mesh position func)
  "Iterating over basis functions in 3D"
  (declare (cl-mpm/mesh::mesh mesh))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (pos-vec position)
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0) (tref pos-vec 2 0)))
           (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec #'floor)))
      (declare (dynamic-extent pos pos-index pos-vec))
      (loop for dx from 0 to 1
            do (loop for dy from 0 to 1
                     do (loop for dz from 0 to 1
                              do (let* ((id (mapcar #'+ pos-index (list dx dy dz))))
                                   (declare (dynamic-extent id))
                                   (when (cl-mpm/mesh:in-bounds mesh id)
                                     (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                            (node (cl-mpm/mesh:get-node mesh id))
                                            (weights (mapcar (lambda (x) (cl-mpm/shape-function::shape-linear x h)) dist))
                                            (weight (reduce #'* weights))
                                            (lin-grads (mapcar (lambda (d)
                                                                 (cl-mpm/shape-function::shape-linear-dsvp d h))
                                                               dist))
                                            (grads (cl-mpm/shape-function::grads-3d weights lin-grads)))
                                       (declare
                                        (double-float weight)
                                        (dynamic-extent dist weights))
                                       (when (< 0d0 weight)
                                         (funcall func mesh node weight grads)))))))))))

(declaim (inline iterate-over-neighbours-point-linear-simd)
         (ftype (function (cl-mpm/mesh::mesh magicl:matrix/double-float function) (values))
                iterate-over-neighbours-point-linear-simd)
         )
(defun iterate-over-neighbours-point-linear-simd (mesh position func)
  "A fast implemenntation of 2D linear basis function iteration"
  (declare (cl-mpm/mesh::mesh mesh)
           (magicl:matrix/double-float position)
           (function func))
  (labels ((simd-abs (vec)
           (sb-simd-avx:f64.2-and vec (sb-simd-avx:f64.2-not -0.0d0)))
         ;; (in-bounds-simd (pos)
         ;;   t)
         (linear-weight-simd (dist h)
           (declare (double-float h))
           ;;Add an abs
           (sb-simd-avx:f64.2- 1d0 (sb-simd-avx:f64.2/ (simd-abs dist) h)))
         (linear-grads-simd (dist h)
           (declare (double-float h))
           (sb-simd-avx:f64.2-if (sb-simd-avx:f64.2> dist 0d0) (/ 1d0 h) (/ -1d0 h))
           dist)
         (in-bounds-simd (mesh dist)
           t)
         )
    (progn
      (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
             (pos-vec (sb-simd-avx:f64.2-aref (magicl::matrix/double-float-storage position) 0))
             (pos-index (sb-simd-avx:f64.2-floor
                         (sb-simd-avx:f64.2/ pos-vec h))))
        (declare (sb-simd-avx:f64.2 pos-vec)
                 (double-float h))
        (loop for dx fixnum from 0 to 1
              do (loop for dy fixnum from 0 to 1
                       do (let* ((id-vec
                                   (sb-simd-avx:f64.2+ pos-index (sb-simd-avx:make-f64.2 dx dy)))
                                 ;; (id (mapcar (lambda (x) (truncate x))
                                 ;;             (multiple-value-list (sb-simd-avx:f64.2-values id-vec))))
                                 (id
                                   (append
                                    (mapcar (lambda (x) (truncate (the double-float x)))
                                            (multiple-value-list (sb-simd-avx:f64.2-values id-vec)))
                                    '(0)))
                                 )
                            (declare (dynamic-extent id))
                            (when (cl-mpm/mesh:in-bounds mesh id)
                                ;(in-bounds-simd mesh id-vec)
                              (let* ((dist (sb-simd-avx:f64.2-
                                            pos-vec
                                            (sb-simd-avx:f64.2* id-vec h)))
                                     (node (cl-mpm/mesh:get-node mesh id))
                                     (weights (linear-weight-simd dist h))
                                     (weight (sb-simd-avx::f64.2-horizontal* weights))
                                     (grads-vec (sb-simd-avx:f64.2*
                                                 (linear-grads-simd dist h)
                                                 (sb-simd-avx:f64.2-shuffle weights weights 1)))
                                     (grads (multiple-value-list (sb-simd-avx:f64.2-values grads-vec)))
                                     )
                                (declare (double-float weight))
                                (when (< 0d0 weight)
                                  (funcall func mesh node weight grads))))
                            )))))))


(declaim ;(inline iterate-over-neighbours-shape-gimp-simd)
         (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle function) (values))
                iterate-over-neighbours-shape-gimp-simd))
(defun iterate-over-neighbours-shape-gimp-simd (mesh mp func)
  "Iterate over a gimp domains neighbours in 2D using SIMD constructs"
  (declare (type cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp)
           (function func)
           (optimize (speed 3) (safety 0) (debug 0))
           )
  (progn
    (with-accessors ((pos-vec cl-mpm/particle:mp-position)
                     (d0 cl-mpm/particle::mp-domain-size))
        mp
      (let ((h (the double-float (cl-mpm/mesh:mesh-resolution mesh))))
        ;; (flet ((center-diff (x)
        ;;          (declare (double-float x h))
        ;;          (- x (the double-float (* h (the double-float
        ;;                                           (fround (the double-float (/ x h))))))))
        ;;        (center-id (x)
        ;;          (round x h))
        ;;        ))
        (let* ((pa (sb-simd-avx:f64.2-aref (magicl::matrix/double-float-storage pos-vec) 0))
               (d0a (sb-simd-avx:f64.2-aref (magicl::matrix/double-float-storage d0) 0))
               (ia (sb-simd-avx:f64.2-round (sb-simd-avx:f64.2*
                                             pa (/ 1d0 h))))
               (ca (sb-simd-avx:f64.2- pa (sb-simd-avx:f64.2* ia h)))
               ;;Lower bound of domain
               (dfa (sb-simd-avx:f64.2-floor
                     (sb-simd-avx:f64.2/
                      (sb-simd-avx:f64.2-
                       ca
                       (sb-simd-avx:f64.2*
                        d0a 0.5d0))
                      h)))
               ;;Upper bound of domain
               (dca (sb-simd-avx:f64.2-ceiling
                     (sb-simd-avx:f64.2/
                      (sb-simd-avx:f64.2+
                       ca
                       (sb-simd-avx:f64.2*
                        d0a 0.5d0)
                       )
                      h)))
               )
          ;; (declare ((simple-array double-float *) pa))
          ;; (declare (dynamic-extent pa))
          (multiple-value-bind (ix iy) (sb-simd-avx:f64.2-values ia)
            (multiple-value-bind (dox doy) (sb-simd-avx:f64.2-values d0a)
              (multiple-value-bind (cx cy) (sb-simd-avx:f64.2-values ca)
                (multiple-value-bind (dxf dyf) (sb-simd-avx:f64.2-values dfa)
                  (multiple-value-bind (dxc dyc) (sb-simd-avx:f64.2-values dca)
                    (declare (type double-float h cx cy dox doy)
                                        ;(type integer dxf dxc dyf dyc ix iy)
                             )
                    (loop for dx fixnum from (the fixnum (truncate dxf))
                            to (the fixnum (truncate dxc))
                          do (loop for dy fixnum from (the fixnum (truncate dyf))
                                     to (the fixnum (truncate dyc))
                                   do
                                      (let* ((id (list (+ (the fixnum (truncate ix)) dx)
                                                       (+ (the fixnum (truncate iy)) dy)
                                                       0
                                                       )))
                                        (declare (dynamic-extent id))
                                        (when (cl-mpm/mesh:in-bounds mesh id)
                                          (let* ((dist (sb-simd-avx:f64.2-
                                                        ca
                                                        (sb-simd-avx:f64.2*
                                                         (sb-simd-avx:make-f64.2 dx dy)
                                                         h))))
                                            (multiple-value-bind (distx disty) (sb-simd-avx:f64.2-values dist)
                                              (declare (type double-float distx disty))
                                              (let* ((weightsx (the double-float (cl-mpm/shape-function::shape-gimp-fast distx (* 0.5d0 dox) h)))
                                                     (weightsy (the double-float (cl-mpm/shape-function::shape-gimp-fast disty (* 0.5d0 doy) h)))
                                                     (weight (* weightsx weightsy))
                                                     ;; #+cl-mpm-fbar
                                                     (weights-fbar-x (the double-float (cl-mpm/shape-function::shape-gimp-fbar distx (* 0.5d0 dox) h)))
                                                     ;; #+cl-mpm-fbar
                                                     (weights-fbar-y (the double-float (cl-mpm/shape-function::shape-gimp-fbar disty (* 0.5d0 doy) h)))
                                                     (weight-fbar (* weights-fbar-x weights-fbar-y))
                                                     )
                                                (when (< 0d0 weight)
                                                  (let* ((node (cl-mpm/mesh:get-node mesh id))
                                                         (gradx (* (cl-mpm/shape-function::shape-gimp-dsvp distx (* 0.5d0 dox) h)
                                                                   (the double-float weightsy)))
                                                         (grady (* (cl-mpm/shape-function::shape-gimp-dsvp disty (* 0.5d0 doy) h)
                                                                   (the double-float weightsx)))
                                                         (gradz 0d0)
                                                         (grads-fbar
                                                           (list (* weights-fbar-y (cl-mpm/shape-function::shape-gimp-dsvp distx (* 0.5d0 dox) h))
                                                                 (* weights-fbar-x (cl-mpm/shape-function::shape-gimp-dsvp disty (* 0.5d0 doy) h))
                                                                 0d0))
                                                         )
                                                    (declare (double-float gradx grady))
                                                    ;; (when (sb-vm:current-float-trap :underflow :overflow :invalid :divide-by-zero)
                                                    ;;   (error "FLOATS FUCKED"))
                                                    (funcall func mesh mp node
                                                             weight (list gradx grady gradz)
                                                             weight-fbar grads-fbar))
                                                  )))
                                            )))))
                    ))))))))))
;;This is more consise but half as fast
(defun iterate-over-neighbours-shape-gimp (mesh mp func)
  "Iterate over a gimp domains neighbours in 3D using simple lisp constructs"
  (declare (type cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp)
           (function func))
  (progn
    (with-accessors ((pos-vec cl-mpm/particle:mp-position)
                     (d0 cl-mpm/particle::mp-domain-size))
        mp
      (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
             (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0) (tref pos-vec 2 0)))
             (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec))
             )
        (declare (dynamic-extent pos pos-index))
        (declare (type double-float h)
                 (type list pos pos-index))
        (loop for dx from -1 to 1
              do (loop for dy from -1 to 1
                       do (loop for dz from -1 to 1
                                do (let* ((id (mapcar #'+ pos-index (list dx dy dz))))
                                     (declare (dynamic-extent id))
                                     (when (cl-mpm/mesh:in-bounds mesh id)
                                       (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                              (domain (loop for x across (magicl::matrix/double-float-storage d0) collect (* 0.5d0 (the double-float x))))
                                              (weights (mapcar (lambda (x l)
                                                                 (cl-mpm/shape-function::shape-gimp x l h))
                                                               dist domain))
                                              (weight (reduce #'* weights))
                                              )
                                         (declare (dynamic-extent domain))
                                         (declare (type double-float h)
                                                  (type list domain))
                                         (when (< 0d0 weight)
                                           (let* ((node (cl-mpm/mesh:get-node mesh id))
                                                  (lin-grads (mapcar (lambda (d l)
                                                                       (cl-mpm/shape-function::shape-gimp-dsvp d l h))
                                                                     dist domain))
                                                  (grads (cl-mpm/shape-function::grads-3d weights lin-grads))
                                                  )
                                             (funcall func mesh mp node weight grads 0d0 nil)))))))))))))

(defun iterate-over-neighbours-shape-gimp-2d (mesh mp func)
  "Iterate over a gimp domains neighbours in 2D unrolled, this version is more performant but worse than SIMD"
  (declare (type cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp)
           (function func)
           (optimize (speed 3) (safety 0) (debug 0))
           )
  (progn
    (with-accessors ((pos-vec cl-mpm/particle:mp-position)
                     (d0 cl-mpm/particle::mp-domain-size))
        mp
      (let ((h (the double-float (cl-mpm/mesh:mesh-resolution mesh))))
        (flet ((center-diff (x)
                 (declare (double-float x h))
                 (- x (the double-float
                           (* h
                              (the double-float
                                   (fround (the double-float (/ x h))))))))
               (center-id (x)
                 (round x h))
               )
          (let* ((pa (magicl::matrix/double-float-storage pos-vec))
                 (da (magicl::matrix/double-float-storage d0))
                 (px (the double-float (aref pa 0)))
                 (py (the double-float (aref pa 1)))
                 (ix (the fixnum (truncate (center-id px))))
                 (iy (the fixnum (truncate (center-id py))))
                 (cx (- px (* h ix)))
                 (cy (- py (* h iy)))
                 (dox (* 0.5d0 (the double-float (aref da 0))))
                 (doy (* 0.5d0 (the double-float (aref da 1))))
                 (dxf (the fixnum (truncate (ffloor   (- cx dox) h))))
                 (dxc (the fixnum (truncate (fceiling (+ cx dox) h))))
                 (dyf (the fixnum (truncate (ffloor   (- cy doy) h))))
                 (dyc (the fixnum (truncate (fceiling (+ cy doy) h))))
                 )
            (declare ((simple-array double-float *) pa))
            ;; (declare (dynamic-extent pa))
            (declare (type double-float h cx cy dox doy px py )
                     (type integer dxf dxc dyf dyc ix iy )
                     )
            (loop for dx from dxf to dxc
                  do (loop for dy from dyf to dyc
                           do
                              (let* ((id (list (the fixnum (+ ix dx))
                                               (the fixnum (+ iy dy))
                                               (the fixnum 0)
                                               )))
                                (declare (dynamic-extent id))
                                (when (cl-mpm/mesh:in-bounds mesh id)
                                  (let* ((distx (- cx (* h dx)))
                                         (disty (- cy (* h dy)))
                                         (weightsx (cl-mpm/shape-function::shape-gimp-fast distx dox h))
                                         (weightsy (cl-mpm/shape-function::shape-gimp-fast disty doy h))
                                         (weight (* weightsx weightsy))
                                         )
                                    (declare ;(type double-float h)
                                     (double-float weight weightsx weightsy distx disty )
                                     )
                                    (when (< 0d0 weight)
                                      (let* ((node (cl-mpm/mesh:get-node mesh id))
                                             (gradx (* (cl-mpm/shape-function::shape-gimp-dsvp distx dox h)
                                                       weightsy 
                                                       ))
                                             (grady (* (cl-mpm/shape-function::shape-gimp-dsvp disty doy h)
                                                       weightsx
                                                       ))
                                             (gradz 0d0))
                                        (declare (double-float gradx grady gradz))
                                        (funcall func mesh mp node weight (list gradx grady gradz) 0d0 nil))))))))))))))

(declaim ;(inline iterate-over-neighbours-shape-gimp)
         (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle function) (values))
                iterate-over-neighbours-shape-gimp-3d))
(defun iterate-over-neighbours-shape-gimp-3d (mesh mp func)
  "Iterate over gimp neighbours in 3D, unrolled for speed"
  (declare (type cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp)
           (function func)
           (optimize (speed 3) (safety 0) (debug 0))
           )
  (progn
    (with-accessors ((pos-vec cl-mpm/particle:mp-position)
                     (d0 cl-mpm/particle::mp-domain-size))
        mp
      (let ((h (the double-float (cl-mpm/mesh:mesh-resolution mesh))))
        (flet ((center-id (x)
                 (round x h))
               )
          (let* ((pa (magicl::matrix/double-float-storage pos-vec))
                 (da (magicl::matrix/double-float-storage d0))
                 (px (the double-float (aref pa 0)))
                 (py (the double-float (aref pa 1)))
                 (pz (the double-float (aref pa 2)))
                 (ix (the fixnum (truncate (center-id px))))
                 (iy (the fixnum (truncate (center-id py))))
                 (iz (the fixnum (truncate (center-id pz))))
                 (cx (- px (* h ix)))
                 (cy (- py (* h iy)))
                 (cz (- pz (* h iz)))
                 (dox (* 0.5d0 (the double-float (aref da 0))))
                 (doy (* 0.5d0 (the double-float (aref da 1))))
                 (doz (* 0.5d0 (the double-float (aref da 2))))
                 (dxf (the fixnum (truncate (ffloor   (- cx dox) h))))
                 (dxc (the fixnum (truncate (fceiling (+ cx dox) h))))
                 (dyf (the fixnum (truncate (ffloor   (- cy doy) h))))
                 (dyc (the fixnum (truncate (fceiling (+ cy doy) h))))
                 (dzf (the fixnum (truncate (ffloor   (- cz doz) h))))
                 (dzc (the fixnum (truncate (fceiling (+ cz doz) h))))
                 )
            (declare ((simple-array double-float *) pa))
            ;; (declare (dynamic-extent pa))
            (declare (type double-float h cx cy dox doy px py pz doz cz)
                     (type integer dxf dxc dyf dyc ix iy dzf dzc iz)
                     )
            (loop for dx from dxf to dxc
                  do (loop for dy from dyf to dyc
                           do
                              (loop for dz from dzf to dzc
                                    do
                              (let* ((id (list (the fixnum (+ ix dx))
                                               (the fixnum (+ iy dy))
                                               (the fixnum (+ iz dz))
                                               )))
                                (declare (dynamic-extent id))
                                (when (cl-mpm/mesh:in-bounds mesh id)
                                  (let* ((distx (- cx (* h dx)))
                                         (disty (- cy (* h dy)))
                                         (distz (- cz (* h dz)))
                                         (weightsx (cl-mpm/shape-function::shape-gimp-fast distx dox h))
                                         (weightsy (cl-mpm/shape-function::shape-gimp-fast disty doy h))
                                         (weightsz (cl-mpm/shape-function::shape-gimp-fast distz doz h))
                                         (weight (* weightsx weightsy weightsz))

                                         (weights-fbar-x (the double-float (cl-mpm/shape-function::shape-gimp-fbar distx (* 0.5d0 dox) h)))
                                         (weights-fbar-y (the double-float (cl-mpm/shape-function::shape-gimp-fbar disty (* 0.5d0 doy) h)))
                                         (weights-fbar-z (the double-float (cl-mpm/shape-function::shape-gimp-fbar distz (* 0.5d0 doz) h)))
                                         (weight-fbar (* weights-fbar-x weights-fbar-y weights-fbar-z))
                                         ;; (weight-fbar 0d0)
                                         )
                                    (declare ;(type double-float h)
                                     (double-float weight weightsx weightsy weightsz distx disty distz)
                                     )
                                    (when (< 0d0 weight)
                                      (let* ((node (cl-mpm/mesh:get-node mesh id))
                                             (gradx (* (cl-mpm/shape-function::shape-gimp-dsvp distx dox h)
                                                       weightsy weightsz
                                                       ))
                                             (grady (* (cl-mpm/shape-function::shape-gimp-dsvp disty doy h)
                                                       weightsx weightsz
                                                       ))
                                             (gradz (* (cl-mpm/shape-function::shape-gimp-dsvp distz doz h)
                                                       weightsx weightsy
                                                       ))
                                             (grads-fbar
                                               (list (* weights-fbar-y weights-fbar-z
                                                        (cl-mpm/shape-function::shape-gimp-dsvp distx (* 0.5d0 dox) h))
                                                     (* weights-fbar-x weights-fbar-z
                                                        (cl-mpm/shape-function::shape-gimp-dsvp disty (* 0.5d0 doy) h))
                                                     (* weights-fbar-x weights-fbar-y
                                                        (cl-mpm/shape-function::shape-gimp-dsvp distz (* 0.5d0 doz) h))))
                                             ;; (grads-fbar (list 0d0 0d0 0d0))
                                             )
                                        (declare (double-float gradx grady gradz))
                                        (funcall func mesh mp node
                                                 weight (list gradx grady gradz)
                                                 weight-fbar grads-fbar)))))))))))))))

(defun iterate-over-neighbours-shape-gimp-simd-3d (mesh mp func)
  "Iterate over gimp neighbours in 3D, unrolled for speed"
  (declare (type cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp)
           (function func)
           (optimize (speed 3) (safety 0) (debug 0))
           )
  (progn
    (with-accessors ((pos-vec cl-mpm/particle:mp-position)
                     (d0 cl-mpm/particle::mp-domain-size))
        mp
      (let ((h (the double-float (cl-mpm/mesh:mesh-resolution mesh))))
        (flet ((center-id (x)
                 (round x h))
               )
          (let* ((pa (magicl::matrix/double-float-storage pos-vec))
                 (da (magicl::matrix/double-float-storage d0))

                 (ia (let ((arr (make-array 3 :initial-element 0d0 :element-type 'double-float)))
                       (loop for i from 0 to 2 do (setf (aref arr i) (fround (aref arr i) h))) arr))
                 (ca (let ((arr (make-array 3 :initial-element 0d0 :element-type 'double-float)))
                       ;; (loop for i from 0 to 2
                       ;;       for pv across pa
                       ;;       for iv across ia
                       ;;       do (setf (aref arr i) (- pv (* h iv))))
                       (setf
                        (sb-simd-avx:f64.2-aref
                         arr 0)
                        (sb-simd-avx:f64.2-
                         (sb-simd-avx:f64.2-aref
                          pa 0)
                         (sb-simd-avx:f64.2*
                          (sb-simd-avx:f64.2-aref
                           ia 0)
                          h))
                        (aref arr 2)
                        (- (aref pa 2) (* h (aref ia 2))))
                       arr))

                 (doa (let ((arr (make-array 3 :initial-element 0d0 :element-type 'double-float)))
                        (setf
                         (sb-simd-avx:f64.2-aref
                          arr 0)
                         (sb-simd-avx:f64.2*
                          (sb-simd-avx:f64.2-aref
                           da 0)
                          0.5d0)
                         (aref arr 2)
                         (* 0.5d0 (aref da 2)))
                       ;; (loop for i from 0 to 2
                       ;;       for dv across da
                       ;;       do (setf (aref arr i) (* 0.5d0 dv)))
                        arr))
                 (dfa (let ((arr (make-array 3 :initial-element 0d0 :element-type 'double-float)))
                       (loop for i from 0 to 2
                             for cv across ca
                             for dov across doa
                             do (setf (aref arr i) (ffloor (- cv dov) h))) arr))
                 (dca (let ((arr (make-array 3 :initial-element 0d0 :element-type 'double-float)))
                        (loop for i from 0 to 2
                              for cv across ca
                              for dov across doa
                              do (setf (aref arr i) (fceiling (+ cv dov) h))) arr))
                 (dfia (let ((arr (make-array 3 :initial-element 0 :element-type 'fixnum)))
                         (loop for v across dfa
                               for i from 0
                               do (setf (aref arr i) (truncate v))) arr))
                 (dcia (let ((arr (make-array 3 :initial-element 0 :element-type 'fixnum)))
                         (loop for v across dca
                               for i from 0
                               do (setf (aref arr i) (truncate v))) arr))

                 (iia (let ((arr (make-array 3 :initial-element 0 :element-type 'fixnum)))
                       (loop for i from 0 to 2 do (setf (aref arr i) (truncate (aref ia i)))) arr))
                 )
            (declare ((simple-array double-float (3)) pa da ia ca doa dfa dca))
            (declare (dynamic-extent pa da ia ca doa dfa dca dfia dcia iia))
            (loop for dx fixnum from (aref dfia 0) to (aref dcia 0)
                  do (loop for dy fixnum from (aref dfia 1) to (aref dcia 1)
                           do
                              (loop for dz fixnum from (aref dfia 2) to (aref dcia 2)
                                    do
                                       (let* ((id (list (the fixnum (+ (the fixnum (aref iia 0)) dx))
                                                        (the fixnum (+ (the fixnum (aref iia 1)) dy))
                                                        (the fixnum (+ (the fixnum (aref iia 2)) dz))
                                                        )))
                                         (declare (dynamic-extent id))
                                         (when (cl-mpm/mesh:in-bounds mesh id)
                                           (let* ((dista
                                                    (let ((arr (make-array 3 :initial-element 0d0 :element-type 'double-float)))
                                                      (declare ((simple-array double-float (3))))
                                                           (setf
                                                            (sb-simd-avx:f64.2-aref
                                                             arr 0)
                                                            (sb-simd-avx:f64.2-
                                                             (sb-simd-avx:f64.2-aref
                                                              ca 0)
                                                             (sb-simd-avx:f64.2*
                                                              (sb-simd-avx:make-f64.2
                                                               dx
                                                               dy)
                                                              h))
                                                            ;; (aref arr 0) (- (aref ca 0) (* h dx))
                                                            ;; (aref arr 1) (- (aref ca 1) (* h dy))
                                                            (aref arr 2) (- (aref ca 2) (* h dz)))
                                                            arr))
                                                  (weightsx (cl-mpm/shape-function::shape-gimp-fast
                                                             (aref dista 0) (aref doa 0) h))
                                                  (weightsy (cl-mpm/shape-function::shape-gimp-fast
                                                             (aref dista 1) (aref doa 1) h))
                                                  (weightsz (cl-mpm/shape-function::shape-gimp-fast
                                                             (aref dista 2) (aref doa 2) h))
                                                  (weight (* weightsx weightsy weightsz))
                                                  )
                                             (declare ;(type double-float h)
                                              (double-float weight weightsx weightsy weightsz)
                                              )
                                             (when (< 0d0 weight)
                                               (let* ((node (cl-mpm/mesh:get-node mesh id))
                                                      (gradx (* (cl-mpm/shape-function::shape-gimp-dsvp (aref dista 0) (aref doa 0) h)
                                                                weightsy weightsz
                                                                ))
                                                      (grady (* (cl-mpm/shape-function::shape-gimp-dsvp (aref dista 1) (aref doa 1) h)
                                                                weightsx weightsz
                                                                ))
                                                      (gradz (* (cl-mpm/shape-function::shape-gimp-dsvp (aref dista 2) (aref doa 2) h)
                                                                weightsx weightsy
                                                                ))
                                                      )
                                                 (declare (double-float gradx grady gradz))
                                                 (funcall func mesh mp node
                                                          weight (list gradx grady gradz)
                                                          0d0 (list 0d0 0d0 0d0))))))))))))))))


(defun make-knot-list (mesh pos)
  "Function for maybe being able to make a bspline knot list"
  (list
   (let ((dy 0))
     (loop for dx from -4 to 4
           collect (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos (list dx dy)))))
   (let ((dx 0))
     (loop for dy from -4 to 4
           collect (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos (list dx dy)))))))
(declaim (inline iterate-over-neighbours-shape-bspline))
(defun iterate-over-neighbours-shape-bspline (mesh mp func)
  "Iterate over 2D bsplines, however boundary adjustment doesn't work"
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (pos-vec (cl-mpm/particle:mp-position mp))
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0) 0))
           ;; (pos-index (magicl:scale pos-vec (/ 1 h)))
           (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec #'round))
           (border (not (and (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list 2 2 0)))
                             (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list -2 -2 0)))
                             ;; (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list 2 -2 0)))
                             ;; (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list -2 2 0)))
                             )))
           )
      (if border
          (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
                 (pos-vec (cl-mpm/particle:mp-position mp))
                 (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
                 (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec #'floor)))
            (loop for dx from 0 to 1
                  do (loop for dy from 0 to 1
                           do (let* ((id (mapcar #'+ pos-index (list dx dy 0))))
                                (when (cl-mpm/mesh:in-bounds mesh id)
                                  (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                         (node (cl-mpm/mesh:get-node mesh id))
                                         (weights (mapcar (lambda (x) (cl-mpm/shape-function::shape-linear x h)) dist))
                                         (weight (reduce #'* weights))
                                         (grads (mapcar (lambda (d w) (* (cl-mpm/shape-function::shape-linear-dsvp d h) w))
                                                        dist (nreverse weights))))
                                    (when (< 0d0 weight)
                                      (funcall func mesh mp node weight (append grads (list 0d0)) 0d0 (list 0d0 0d0 0d0)))))))))

          (loop for dx from -1 to 1
                do (loop for dy from -1 to 1
                         do (let* ((id (mapcar #'+ pos-index (list dx dy 0))))
                              (when (cl-mpm/mesh:in-bounds mesh id)
                                (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                       (node (cl-mpm/mesh:get-node mesh id))
                                       (weights (mapcar (lambda (x) (cl-mpm/shape-function:shape-bspline x h)) dist))
                                       (weight (reduce #'* weights))
                                       (gradx (* (cl-mpm/shape-function::shape-bspline-dsvp (nth 0 dist) h)
                                                 (the double-float (nth 1 weights))))
                                       (grady (* (cl-mpm/shape-function::shape-bspline-dsvp (nth 1 dist) h)
                                                 (the double-float (nth 0 weights))))
                                       (gradz 0d0)
                                       )
                                  (funcall func mesh mp node weight (list gradx grady gradz) 0d0 (list 0d0 0d0 0d0)))))))))))

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
                                         (cl-mpm/fastmath::fast-.* dsvp dsvp)))))
                (sb-thread:with-mutex (node-lock)
                    (setf node-temp
                          (+ node-temp weighted-temp))
                    (setf node-dtemp
                          (+ node-dtemp weighted-dtemp))))))))


(declaim
 (notinline p2g-mp)
 (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle) (values))
                p2g-mp))
(defun p2g-mp (mesh mp)
  "P2G transfer for one MP"
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (with-accessors ((mp-vel  cl-mpm/particle:mp-velocity)
                   (mp-mass cl-mpm/particle:mp-mass)
                   (mp-volume cl-mpm/particle:mp-volume)
                   (mp-pmod cl-mpm/particle::mp-p-modulus)
                   (mp-damage cl-mpm/particle::mp-damage)
                   ) mp
    (declare (type double-float mp-mass mp-volume))
    (let (
          (mp-mass mp-mass)
          (mp-vel mp-vel)
          (mp-volume mp-volume)
          (mp-pmod mp-pmod)
          (mp-damage mp-damage)
          )
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
           (declare (type double-float node-mass node-volume mp-volume mp-pmod mp-damage node-svp-sum svp node-p-wave node-damage)
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
           ;;Ideally we include these generic functions for special mapping operations, however they are slow
           ;; (special-p2g mp node svp dsvp)
           )))))
  (values))

(declaim (notinline p2g))
(defun p2g (mesh mps)
  "Map particle momentum to the grid"
  (declare (type (array cl-mpm/particle:particle) mps) (cl-mpm/mesh::mesh mesh))
  (iterate-over-mps
   mps
   (lambda (mp)
     (p2g-mp mesh mp))))

(declaim (notinline p2g-force-mp)
         (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle) (values)) p2g-force-mp)
         )
(defun p2g-force-mp (mesh mp)
  "Map particle forces to the grid for one mp"
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (with-accessors ((mp-vel  cl-mpm/particle:mp-velocity)
                   (mp-mass cl-mpm/particle:mp-mass)
                   ) mp
    (declare (type double-float mp-mass))
    (let ((dsvp (cl-mpm/utils::dsvp-3d-zeros)))
      (declare (dynamic-extent dsvp))
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
                          (node-int-force cl-mpm/mesh::node-internal-force)
                          (node-ext-force cl-mpm/mesh::node-external-force)
                          (node-lock  cl-mpm/mesh:node-lock)) node
           (declare (double-float node-mass)
                    (boolean node-active)
                    (sb-thread:mutex node-lock)
                    (magicl:matrix/double-float node-vel node-force node-int-force node-ext-force))
           (when node-active
             (sb-thread:with-mutex (node-lock)
               (det-ext-force mp node svp node-ext-force)
               (cl-mpm/shape-function::assemble-dsvp-3d-prealloc grads dsvp)
               (det-int-force mp dsvp node-int-force)

               ;; (det-ext-force mp node svp node-force)
               ;; (cl-mpm/shape-function::assemble-dsvp-3d-prealloc grads dsvp)
               ;; (det-int-force mp dsvp node-force)

               ;; (cl-mpm/fastmath::fast-zero node-force)
               ;; (cl-mpm/fastmath::fast-.+ node-force node-int-force)
               ;; (cl-mpm/fastmath::fast-.+ node-force node-ext-force)

               (cl-mpm/fastmath::fast-.+-vector node-int-force node-ext-force node-force)
               ;; (magicl:.+ node-int-force node-ext-force node-force)
               ))
           )
         ))))
  (values))

(declaim (inline p2g-force))
(defun p2g-force (mesh mps)
  "Map particle forces to the grid"
  (declare (type (array cl-mpm/particle:particle) mps) (cl-mpm/mesh::mesh mesh))
  (iterate-over-mps
   mps
   (lambda (mp)
     (p2g-force-mp mesh mp))))

(defgeneric special-g2p (mesh mp node svp grads)
  (:documentation "G2P behaviour for specific features")
  (:method (mesh mp node svp grads)))

(defmethod special-g2p (mesh (mp cl-mpm/particle::particle-thermal) node svp grads)
  (declare (ignore mesh))
  "Map grid to particle for one mp-node pair"
  (with-accessors ((node-temp cl-mpm/mesh:node-temperature)) node
    (with-accessors ((temp cl-mpm/particle::mp-temperature)) mp
      (progn
        (setf temp (+ temp (* node-temp svp)))))))

(defmethod special-g2p (mesh (mp cl-mpm/particle::particle-damage) node svp grads)
  (declare (ignore mesh))
  "Map grid to particle for one mp-node pair"
  (with-accessors ((node-damage cl-mpm/mesh::node-damage)
                   (node-temp cl-mpm/mesh::node-temp)) node
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
 (notinline g2p-mp)
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
                   (fixed-velocity cl-mpm/particle::mp-fixed-velocity)
                   )
    mp
    (let* (
           (mapped-vel (cl-mpm/utils:vector-zeros))
           )
      ;; (declare (dynamic-extent mapped-vel))
      (progn
        ;;With special operations we need to reset some params for g2p
        ;; (reset-mps-g2p mp)
        (setf temp 0d0)
        )
      (cl-mpm/fastmath::fast-zero acc)
      ;; Map variables
      (iterate-over-neighbours
       mesh mp
       (lambda (mesh mp node svp grads fsvp fgrads)
         (declare
          (ignore mp mesh fsvp fgrads)
          (cl-mpm/mesh::node node)
          (cl-mpm/particle:particle mp)
          ( double-float svp))
         (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
                          (node-acc cl-mpm/mesh:node-acceleration)
                          (node-scalar cl-mpm/mesh::node-boundary-scalar)
                          (node-active cl-mpm/mesh:node-active)
                          ) node
           (declare (double-float node-scalar temp)
                    (boolean node-active))
           (when node-active
             (cl-mpm/fastmath::fast-fmacc mapped-vel node-vel svp)
             (cl-mpm/fastmath::fast-fmacc acc node-acc svp)
             (incf temp (* svp node-scalar))
             ;;With special operations we want to include this operation
             #+cl-mpm-special (special-g2p mesh mp node svp grads)
             )
           )
         ;; (g2p-mp-node mp node svp grads)
         ))
      ;;Update particle
      (progn
        ;;Invalidate shapefunction/gradient cache
        (setf (fill-pointer nc) 0)

        ;;PIC
        #+ :cl-mpm-pic (cl-mpm/utils::voigt-copy-into mapped-vel vel)
        (cl-mpm/fastmath::fast-scale mapped-vel dt)
        (cl-mpm/fastmath::simd-add pos mapped-vel)
        (cl-mpm/fastmath::simd-add disp mapped-vel)
        (unless (cl-mpm/particle::mp-penalty-contact mp)
          (cl-mpm/fastmath:fast-zero (cl-mpm/particle:mp-penalty-frictional-force mp)))
        (setf (cl-mpm/particle::mp-penalty-contact mp) nil)
        ;;FLIP
        #- :cl-mpm-pic (cl-mpm/fastmath:fast-fmacc vel acc dt)
        ))))
(defgeneric pre-particle-update-hook (particle dt)
  )
(defmethod pre-particle-update-hook (particle dt))

(declaim (notinline g2p))
(defun g2p (mesh mps dt)
  (declare (cl-mpm/mesh::mesh mesh) (array mps))
  "Map grid values to all particles"
  (iterate-over-mps
   mps
   (lambda (mp)
     (g2p-mp mesh mp dt))))


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
  "Calculate velocity from momentum on a single node"
  (when (cl-mpm/mesh:node-active node)
      (with-accessors ( (mass  node-mass)
                        (vel   node-velocity))
        node
        (declare (double-float mass))
        (progn
          (cl-mpm/fastmath::fast-scale vel (/ 1.0d0 mass))))))

(declaim (notinline calculate-forces)
         (ftype (function (cl-mpm/mesh::node double-float double-float double-float) (vaules)) calculate-forces))
(defun calculate-forces (node damping dt mass-scale)
  "Update forces and nodal velocities with viscous damping"
  (when (cl-mpm/mesh:node-active node)
      (with-accessors ( (mass  node-mass)
                        (vel   node-velocity)
                        (force node-force)
                        (acc   node-acceleration))
        node
        (declare (double-float mass dt damping mass-scale))
        (progn

          (magicl:scale! acc 0d0)
          ;;Set acc to f/m
          (cl-mpm/fastmath:fast-fmacc acc force (/ 1d0 (* mass mass-scale)))
          (cl-mpm/fastmath:fast-fmacc acc vel (/ (* damping -1d0) mass-scale))
          ;; (cl-mpm/fastmath:fast-fmacc acc vel (/ (* damping -1d0) (* mass mass-scale)))
          (cl-mpm/fastmath:fast-fmacc vel acc dt)
          )))
  (values))
(defun calculate-forces-psudo-viscous (node damping dt mass-scale)
  "Update forces and nodal velocities with viscous damping - except without scaling by mass
This allows for a non-physical but viscous damping scheme that is robust to GIMP domains "
  (when (cl-mpm/mesh:node-active node)
    (with-accessors ( (mass  node-mass)
                     (vel   node-velocity)
                     (force node-force)
                     (acc   node-acceleration))
        node
      (declare (double-float mass dt damping mass-scale))
      (progn
        (magicl:scale! acc 0d0)
        ;;Set acc to f/m
        (cl-mpm/fastmath:fast-fmacc acc force (/ 1d0 (* mass mass-scale)))
        (cl-mpm/fastmath:fast-fmacc acc vel (* damping -1d0))
        (cl-mpm/fastmath:fast-fmacc vel acc dt)
        )))
  (values))

(defun calculate-forces-cundall (node damping dt mass-scale)
  "Apply cundall damping to the system"
  (when (cl-mpm/mesh:node-active node)
    (with-accessors ((mass  node-mass)
                     (vel   node-velocity)
                     (force node-force)
                     (acc   node-acceleration)
                     )
        node
      (declare (double-float mass dt damping mass-scale))
      (progn
        (magicl:scale! acc 0d0)
        ;;Set acc to f/m

        (let* ((vel-sign (cl-mpm/utils:vector-zeros))
               (f-s (magicl::matrix/double-float-storage force))
               (v-s (magicl::matrix/double-float-storage vel))
               (fnorm (cl-mpm/fastmath::mag force))
               (vs-s (magicl::matrix/double-float-storage vel-sign)))
          (loop for i from 0 to 2
                do
                   ;;Possible form of damping
                   (when t;(> (* (aref f-s i) (aref v-s i)) 0)
                     (incf (aref f-s i)
                           (*
                            (signum (aref v-s i))
                            (abs (aref f-s i))
                            damping
                            -1d0))))
          (cl-mpm/fastmath:fast-fmacc acc force (/ 1d0 (* mass mass-scale))))
        (cl-mpm/fastmath:fast-fmacc vel acc dt))))
  (values))
(defun calculate-forces-cundall-conservative (node damping dt mass-scale)
  "Apply cundall damping to the system only when doing negative work"
  (when (cl-mpm/mesh:node-active node)
    (with-accessors ((mass  node-mass)
                     (vel   node-velocity)
                     (force node-force)
                     (acc   node-acceleration)
                     )
        node
      (declare (double-float mass dt damping))
      (progn
        (magicl:scale! acc 0d0)
        ;;Set acc to f/m

        (let* ((vel-sign (cl-mpm/utils:vector-zeros))
               (f-s (magicl::matrix/double-float-storage force))
               (v-s (magicl::matrix/double-float-storage vel))
               (fnorm (cl-mpm/fastmath::mag force))
               (vs-s (magicl::matrix/double-float-storage vel-sign)))
          (loop for i from 0 to 2
                do
                   ;;Possible form of damping
                   (when (> (* (aref f-s i) (aref v-s i)) 0)
                     (incf (aref f-s i)
                           (*
                            (signum (aref v-s i))
                            (abs (aref f-s i))
                            damping
                            -1d0)))
                )
          (cl-mpm/fastmath:fast-fmacc acc force (/ 1d0 (* mass mass-scale))))
        (cl-mpm/fastmath:fast-fmacc vel acc dt))))
  (values))

(declaim (inline iterate-over-nodes)
         (ftype (function (cl-mpm/mesh::mesh function) (values)) iterate-over-nodes)
         )
(defun iterate-over-nodes (mesh func)
  "Helper function for iterating over all nodes in a mesh
   Calls func with only the node"
  (declare (type function func))
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
    (declare (type (array cl-mpm/particle:particle) nodes))
    ;; (dotimes (i (array-total-size nodes))
    ;;   (let ((node (row-major-aref nodes i)))
    ;;     (funcall func node)))
    (lparallel:pdotimes (i (array-total-size nodes))
      (let ((node (row-major-aref nodes i)))
        (when node
          (funcall func node))))
    )
  (values))

(declaim (inline iterate-over-nodes-serial)
         (ftype (function (cl-mpm/mesh::mesh function) (values)) iterate-over-nodes-serial)
         )
(defun iterate-over-nodes-serial (mesh func)
  "Helper function for iterating over all nodes in a mesh - in a serial way
Calls func with only the node"
  (declare (type function func))
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
    (declare (type (array cl-mpm/particle:particle) nodes))
    (dotimes (i (array-total-size nodes))
      (let ((node (row-major-aref nodes i)))
        (when node
          (funcall func node)))))
  (values))

(declaim
 (ftype (function (cl-mpm/mesh::mesh function) (values)) iterate-over-cells))
(defun iterate-over-cells (mesh func)
  "Helper function for iterating over all nodes in a mesh
   Calls func with only the node"
  (declare (type function func))
  (let ((cells (cl-mpm/mesh::mesh-cells mesh)))
    (declare (type (array (or cl-mpm/mesh::cell null)) cells))
    (lparallel:pdotimes (i (array-total-size cells))
      (let ((cell (row-major-aref cells i)))
        (when cell
          (funcall func cell)))))
  (values))

(defun iterate-over-cells-serial (mesh func)
  "Helper function for iterating over all nodes in a mesh
   Calls func with only the node"
  (declare (type function func))
  (let ((cells (cl-mpm/mesh::mesh-cells mesh)))
    (declare (type (array (or cl-mpm/mesh::cell null)) cells))
    (dotimes (i (array-total-size cells))
      (let ((cell (row-major-aref cells i)))
        (when cell
          (funcall func cell)))))
  (values))

(declaim (ftype (function
                 ((vector cl-mpm/particle:particle *)
                  function)
                 (values))
                iterate-over-mps))
(defun iterate-over-mps (mps func)
  "Helper function for iterating over all nodes in a mesh
   Calls func with only the node"
  (declare (type function func)
           (type (array cl-mpm/particle:particle) mps))
  ;; (dotimes (i (length mps))
  ;;   (funcall func (aref mps i)))
  (lparallel:pdotimes (i (length mps))
                      (funcall func (aref mps i)))
  (values))
;; 
;; (defun iterate-over-mps-serial (mps func)
;;   "Helper function for iterating over all nodes in a mesh
;;    Calls func with only the node"
;;   (declare (type function func)
;;            (type (array cl-mpm/particle:particle) mps))
;;   (loop for mp across mps
;;         do (funcall func mp))
;;   (values))


(defun update-node-kinematics (mesh dt)
  (iterate-over-nodes mesh
                      (lambda (node)
                        (calculate-kinematics node))))
(defgeneric update-node-forces (sim)
  (:documentation "Update the acceleration from forces and apply any damping"))
(defmethod update-node-forces ((sim mpm-sim))
  (with-accessors ((damping sim-damping-factor)
                   (mass-scale sim-mass-scale)
                   (mesh sim-mesh)
                   (dt sim-dt))
      sim
    (iterate-over-nodes
     mesh
     (lambda (node)
       (calculate-forces node damping dt mass-scale)))))

(defmethod update-node-forces ((sim mpm-sim-quasi-static))
  (with-accessors ((damping sim-damping-factor)
                   (mass-scale sim-mass-scale)
                   (mesh sim-mesh)
                   (dt sim-dt))
      sim
    ;; (break)
    (iterate-over-nodes
     mesh
     (lambda (node)
       (calculate-forces-cundall node damping dt mass-scale)))))



(defun apply-bcs (mesh bcs dt)
  "Apply all normal bcs onto the mesh"
  (declare (cl-mpm/mesh::mesh mesh))
  (with-accessors ( (nodes  mesh-nodes)
                   (nD     mesh-nD)
                   (mc     mesh-count)) mesh
                                        ;each bc is a list (pos value)
    (lparallel:pdotimes (i (array-total-size bcs))
      (let ((bc (aref bcs i)))
        (when bc
          (with-accessors ((node cl-mpm/bc::bc-node)
                           (index cl-mpm/bc::bc-index))
              bc
            (if (or node
                    (not index))
                (cl-mpm/bc:apply-bc bc node mesh dt)
                (progn
                  (setf node (cl-mpm/mesh:get-node mesh index))
                  (if node
                      (cl-mpm/bc:apply-bc bc node mesh dt)
                      (error "BC attempted to get a nil node ~A ~A" bc index))))))))))


;Could include this in p2g but idk
(declaim (notinline calculate-strain-rate)
         (ftype (function (cl-mpm/mesh::mesh  cl-mpm/particle:particle double-float)) calculate-strain-rate))
(defun calculate-strain-rate (mesh mp dt)
  "Calculate the strain rate, stretch rate and vorticity"
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt)
           (optimize (speed 3) (safety 1))
           )
  (with-accessors ((strain-rate cl-mpm/particle:mp-strain-rate)
                   (vorticity cl-mpm/particle:mp-vorticity)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (stretch-tensor-fbar cl-mpm/particle::mp-stretch-tensor-fbar)
                   (velocity-rate cl-mpm/particle::mp-velocity-rate)
                   ) mp
    (declare (magicl:matrix/double-float strain-rate vorticity stretch-tensor stretch-tensor-fbar velocity-rate))
        (progn
          ;; (cl-mpm/fastmath::fast-zero strain-rate)
          ;; (cl-mpm/fastmath::fast-zero vorticity)
          (cl-mpm/fastmath::fast-zero stretch-tensor)
          (cl-mpm/fastmath::fast-zero stretch-tensor-fbar)
          (let (
                ;; (stretch-dsvp (stretch-dsvp-3d-zeros))
                ;; (temp-mult (cl-mpm/utils::stretch-dsvp-voigt-zeros))
                ;; (temp-add (cl-mpm/utils::matrix-zeros))
                )
            ;; (declare (magicl:matrix/double-float stretch-dsvp temp-mult temp-add)
            ;;          ;; (dynamic-extent stretch-dsvp temp-mult temp-add)
            ;;          )
            (iterate-over-neighbours
             mesh mp
             (lambda (mesh mp node svp grads fsvp fgrads)
               (declare (ignore mp svp fsvp))
               (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
                                (node-active cl-mpm/mesh:node-active))
                   node
                 (declare (magicl:matrix/double-float node-vel)
                          (boolean node-active))
                 (when node-active

                   (cl-mpm/shape-function::@-combi-assemble-dstretch-3d grads node-vel stretch-tensor)
                   (cl-mpm/shape-function::@-combi-assemble-dstretch-3d fgrads node-vel stretch-tensor-fbar)

                   ;; (cl-mpm/shape-function::assemble-dstretch-3d-prealloc grads stretch-dsvp)
                   ;; ;; (ecase (cl-mpm/mesh::mesh-nd mesh)
                   ;; ;;   (2 (cl-mpm/shape-function::assemble-dstretch-2d-prealloc grads stretch-dsvp))
                   ;; ;;   (3 (cl-mpm/shape-function::assemble-dstretch-3d-prealloc grads stretch-dsvp)))

                   ;; (cl-mpm/fastmath::@-stretch-vec stretch-dsvp node-vel temp-mult)
                   ;; (cl-mpm/utils::voight-to-stretch-prealloc temp-mult temp-add)

                   ;; ;; (setf stretch-tensor (magicl:.+ stretch-tensor temp-add))
                   ;; ;; (magicl:.+ stretch-tensor temp-add stretch-tensor)
                   ;; (cl-mpm/fastmath::fast-.+-matrix
                   ;;  stretch-tensor
                   ;;  temp-add
                   ;;  stretch-tensor)

                   ;; ;; (cl-mpm/shape-function::assemble-dstretch-3d-prealloc fgrads stretch-dsvp)
                   ;; ;; (cl-mpm/fastmath::@-stretch-vec stretch-dsvp node-vel temp-mult)
                   ;; ;; (cl-mpm/utils::voight-to-stretch-prealloc temp-mult temp-add)
                   ;; ;; (cl-mpm/fastmath::fast-.+-matrix
                   ;; ;;  stretch-tensor-fbar
                   ;; ;;  temp-add
                   ;; ;;  stretch-tensor-fbar)
                   )))))

            ;; (cl-mpm/utils::stretch-to-sym stretch-tensor strain-rate)
            ;; (cl-mpm/utils::stretch-to-skew stretch-tensor vorticity)
            ;; (aops:copy-into (cl-mpm/utils::fast-storage velocity-rate) (cl-mpm/utils::fast-storage strain-rate))
            ;; (setf velocity-rate (magicl:scale strain-rate 1d0))
            (cl-mpm/fastmath::fast-scale stretch-tensor dt)
            (cl-mpm/fastmath::fast-scale stretch-tensor-fbar dt)
            ;; (cl-mpm/fastmath::fast-scale strain-rate dt)
            ;; (cl-mpm/fastmath::fast-scale vorticity dt)
            )))


(defun rotation-matrix (degrees)
  (let ((angle (/ (* pi degrees) 180)))
    (magicl:from-list (list (cos angle) (- (sin angle))
                            (sin angle) (cos angle))
                      '(2 2)
                      :type 'double-float
                      )))
(defun update-strain-linear (mesh mp dt fbar)
  "Linear small strain update algorithm"
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (volume-0 cl-mpm/particle::mp-volume-0)
                   (strain cl-mpm/particle:mp-strain)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   ) mp
    (let ((df (calculate-df mesh mp fbar))
          (dstrain (magicl:scale strain-rate dt)))
                   (progn
                     (setf def (magicl:@ df def))
                     (setf strain (cl-mpm/fastmath::fast-.+-voigt strain dstrain))
                     (setf volume (* volume-0 (magicl:det def)))
                     ;(update-domain-corner mesh mp dt)
                     (update-domain-midpoint mesh mp dt)))))

(declaim (inline update-strain-kirchoff))
(declaim (ftype (function (cl-mpm/mesh::mesh
                           cl-mpm/particle:particle
                           double-float
                           boolean) (values))
               update-strain-kirchoff))

(defun update-strain-kirchoff (mesh mp dt fbar)
  "Finite strain kirchhoff strain update algorithm"
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
                   domain))
    (progn
      (let ((df (calculate-df mesh mp fbar)))
        (progn
          (setf def (magicl:@ df def))

          ;; (let ((temp (cl-mpm/utils:matrix-zeros)))
          ;;   (magicl:mult df def :target temp)
          ;;   (cl-mpm/utils:matrix-copy-into temp def))

          (let (;(initial-strain (cl-mpm/utils::vector-copy strain))
                )
            ;; (aops:copy-into (magicl::matrix/double-float-storage eng-strain-rate)
            ;;                 (magicl::matrix/double-float-storage strain))
            ;;This can may be in lisp or accelerated
            ;; (when (> (abs (magicl:tref strain 4 0)) 0d0)
            ;;   (error "Pre Nonzero out of plane strain with ~A" (loop for v across (magicl::storage strain) collect v)))
            (cl-mpm/ext:kirchoff-update strain df)
            ;; (when (> (abs (magicl:tref strain 4 0)) 0d0)
            ;;   (error "Post Nonzero out of plane strain with ~A" (loop for v across (magicl::storage strain) collect v)))
            ;;Not sure about this engineering strain calculation
            ;; (magicl:.- eng-strain-rate strain eng-strain-rate)
            ;; (setf eng-strain-rate initial-strain)

            ;; (cl-mpm/fastmath::fast-scale eng-strain-rate (/ 1d0 dt))

            ;; (setf eng-strain-rate (magicl:scale! (magicl:.- strain initial-strain) (/ 1d0 dt)))
            )

          ;;Post multiply to turn to eng strain
          (setf volume (* volume (the double-float (magicl:det df))))
          ;; (setf volume (* volume-0 (magicl:det def)))
          (when (<= volume 0d0)
            (error "Negative volume"))
          ;;Stretch rate update
            (update-domain-corner mesh mp dt)
            ;; (update-domain-midpoint mesh mp dt)
          ;; (update-domain-stretch-rate df domain)

          )
          )))
  (values))
(declaim
 (inline update-strain-kirchoff-damaged)
 (ftype (function (cl-mpm/mesh::mesh
                           cl-mpm/particle:particle
                           double-float
                           boolean) (values))
                update-strain-kirchoff-damaged))
(defun update-strain-kirchoff-damaged (mesh mp dt fbar)
  "Finite strain kirchhoff strain update algorithm"
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   (domain cl-mpm/particle::mp-domain-size)
                   (eng-strain-rate cl-mpm/particle::mp-eng-strain-rate)
                   ) mp
    (declare (type double-float volume)
             (type magicl:matrix/double-float
                   domain))
    (progn
      (let ((df (calculate-df mesh mp fbar)))
        (progn
          (let ((def-n (cl-mpm/utils::matrix-zeros)))
            (magicl:mult df def :target def-n)
            (cl-mpm/utils::matrix-copy-into def-n def))
          ;; (aops:copy-into (magicl::matrix/double-float-storage eng-strain-rate)
          ;;                 (magicl::matrix/double-float-storage strain))
          ;; (when (> (abs (magicl:tref strain 4 0)) 0d0)
          ;;   (error "Pre Nonzero out of plane strain with ~A" (loop for v across (magicl::storage strain) collect v)))
          (cl-mpm/utils:voigt-copy-into strain strain-rate)
          (cl-mpm/ext:kirchoff-update strain df)
          ;; (when (> (abs (magicl:tref strain 4 0)) 0d0)
          ;;   (error "Post Nonzero out of plane strain with ~A" (loop for v across (magicl::storage strain) collect v)))
          (cl-mpm/fastmath:fast-.- strain-rate strain strain-rate)
          ;; (magicl:.- eng-strain-rate strain eng-strain-rate)
          ;; (magicl:scale! eng-strain-rate (/ 1d0 dt))
          ;;Post multiply to turn to eng strain
          (setf volume (* volume (the double-float (magicl:det df))))
          (when (<= volume 0d0)
            (error "Negative volume"))
          ;; (update-domain-stretch-rate-damage stretch-tensor (cl-mpm/particle::mp-damage mp) domain
          ;;                                    (cl-mpm/particle::mp-damage-domain-update-rate mp))
          )))
    )
  (values))

(defun update-domain-deformation-rate (domain df)
  "Update the domain length based on the increment of defomation rate"
  (setf (tref domain 0 0) (* (the double-float (tref domain 0 0))
                              (the double-float (tref df 0 0))))
  (setf (tref domain 1 0) (* (the double-float (tref domain 1 0))
                              (the double-float (tref df 1 1))))
  )
(defun update-domain-stretch-rate (df domain)
  "Update the domain length based on the increment of the stretch rate"
  (let ((F (cl-mpm/utils::matrix-zeros)))
    (magicl:mult df df :target F :transb :t)
    (multiple-value-bind (l v) (cl-mpm/utils::eig F)
      (let* ((stretch
              (magicl:@
               v
               (cl-mpm/utils::matrix-from-list
                (list (the double-float (sqrt (the double-float (nth 0 l)))) 0d0 0d0
                      0d0 (the double-float (sqrt (the double-float (nth 1 l)))) 0d0
                      0d0 0d0 (the double-float (sqrt (the double-float (nth 2 l))))))
               (magicl:transpose v)))
            )
        (declare (type magicl:matrix/double-float stretch))
        (setf (tref domain 0 0) (* (the double-float (tref domain 0 0))
                                   (the double-float (tref stretch 0 0))))
        (setf (tref domain 1 0) (* (the double-float (tref domain 1 0))
                                   (the double-float (tref stretch 1 1))))
        (setf (tref domain 2 0) (* (the double-float (tref domain 2 0))
                                   (the double-float (tref stretch 2 2))))
        ))))


(declaim (ftype (function (magicl:matrix/double-float
                           double-float
                           magicl:matrix/double-float
                           double-float) (values))
                update-domain-stretch-rate-damage))
(defun update-domain-stretch-rate-damage (stretch-rate damage domain
                                          damage-domain-rate)
  "Update the domain length based on the increment of the stretch rate"
  (declare (double-float damage damage-domain-rate))
  (let ((df (cl-mpm/utils::matrix-from-list '(1d0 0d0 0d0
                                              0d0 1d0 0d0
                                              0d0 0d0 1d0)))
        (degredation (- 1d0 (* damage-domain-rate damage)))
        (domain-array (cl-mpm/utils:fast-storage domain))
        )

    (cl-mpm/fastmath:fast-.+ df (magicl:scale stretch-rate degredation) df)
    ;; (cl-mpm/fastmath::fast-.+ df (magicl:scale stretch-rate degredation) df)
    (let ((F (cl-mpm/utils::matrix-zeros)))
      (magicl:mult df df :target F :transb :t)
      (multiple-value-bind (l v) (cl-mpm/utils::eig F)
        (destructuring-bind (l1 l2 l3) l
          (declare (double-float l1 l2 l3))
          (let* ((stretch
                   (magicl:@
                    v
                    (cl-mpm/utils::matrix-from-list
                     (list (the double-float (sqrt l1)) 0d0 0d0
                           0d0 (the double-float (sqrt l2)) 0d0
                           0d0 0d0 (the double-float (sqrt l3))))
                    (magicl:transpose v)))
                 )
            (declare (type magicl:matrix/double-float stretch))
            (setf (aref domain-array 0) (* (the double-float (tref domain 0 0))
                                       (the double-float (tref stretch 0 0))))
            (setf (aref domain-array 1) (* (the double-float (tref domain 1 0))
                                       (the double-float (tref stretch 1 1))))
            (setf (aref domain-array 2) (* (the double-float (tref domain 2 0))
                                       (the double-float (tref stretch 2 2))))
            ))))))

(defun update-domain-stretch (def domain domain-0)
  "Update the domain length based on the total stretch rate"
  (let ((F (cl-mpm/utils::matrix-zeros)))
    (magicl:mult def def :target F :transb :t)
    (multiple-value-bind (l v) (cl-mpm/utils::eig F)
      (let* ((stretch
              (magicl:@
               v
               (cl-mpm/utils::matrix-from-list
                (list (the double-float (sqrt (the double-float (nth 0 l)))) 0d0 0d0
                      0d0 (the double-float (sqrt (the double-float (nth 1 l)))) 0d0
                      0d0 0d0 (the double-float (sqrt (the double-float (nth 2 l))))))
               (magicl:transpose v)))
            )
        (declare (type magicl:matrix/double-float stretch))
        (setf (tref domain 0 0) (* (the double-float (tref domain-0 0 0))
                                   (the double-float (tref stretch 0 0))))
        (setf (tref domain 1 0) (* (the double-float (tref domain-0 1 0))
                                   (the double-float (tref stretch 1 1))))
        (setf (tref domain 2 0) (* (the double-float (tref domain-0 2 0))
                                   (the double-float (tref stretch 2 2))))
        ))))
(defun update-domain-midpoint (mesh mp dt)
  "Use a corner tracking scheme to update domain lengths"
  (with-accessors ((position cl-mpm/particle::mp-position)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0)
                   )
      mp
    (let ((nd (the fixnum (cl-mpm/mesh:mesh-nd mesh))))
      (with-accessors ((mesh-size cl-mpm/mesh::mesh-mesh-size))
          mesh
        (let ((diff (make-array 3 :initial-element 0d0))
              (domain-storage (magicl::matrix/double-float-storage domain)))
          (loop for d from 0 below nd
                do
                   (let ((disp (cl-mpm/utils:vector-zeros)))
                     (loop for direction in (list 1d0 -1d0)
                           do
                              (let ((corner (cl-mpm/utils::vector-copy position)))
                                (incf (magicl:tref corner d 0)
                                      (* direction 0.5d0 (aref domain-storage d)))
                                (iterate-over-neighbours-point-linear
                                 mesh corner
                                 (lambda (mesh node svp grads)
                                   (declare (double-float dt svp))
                                   (with-accessors ((vel cl-mpm/mesh:node-velocity))
                                       node
                                     (cl-mpm/fastmath:fast-fmacc corner vel (* dt svp))
                                     )))
                                (cl-mpm/fastmath:fast-fmacc disp corner direction)
                                ))
                     (setf (aref diff d) (magicl:tref disp d 0))))
          (setf
           (aref domain-storage 0) (aref diff 0)
           (aref domain-storage 1) (aref diff 1)
           (aref domain-storage 2) (aref diff 2))
          ;; (incf (aref domain-storage 0) (aref diff 0))
          ;; (incf (aref domain-storage 1) (aref diff 1))
          ;; (incf (aref domain-storage 2) (aref diff 2))
          (if (= 2 nd)
              (let* ((jf  (magicl:det def))
                     (jl  (* (magicl:tref domain 0 0) (magicl:tref domain 1 0)))
                     (jl0 (* (magicl:tref domain-0 0 0) (magicl:tref domain-0 1 0)))
                     (scaling (expt (/ (* jf jl0) jl) (/ 1d0 2d0))))
                (setf (magicl:tref domain 0 0) (* (magicl:tref domain 0 0) scaling)
                      (magicl:tref domain 1 0) (* (magicl:tref domain 1 0) scaling)))
              (let* ((jf  (magicl:det def))
                     (jl  (* (magicl:tref domain 0 0) (magicl:tref domain 1 0) (magicl:tref domain 2 0)))
                     (jl0 (* (magicl:tref domain-0 0 0) (magicl:tref domain-0 1 0) (magicl:tref domain-0 2 0)))
                     (scaling (expt (/ (* jf jl0) jl) (/ 1d0 3d0))))
                (setf (magicl:tref domain 0 0) (* (magicl:tref domain 0 0) scaling)
                      (magicl:tref domain 1 0) (* (magicl:tref domain 1 0) scaling)
                      (magicl:tref domain 2 0) (* (magicl:tref domain 2 0) scaling)
                      ))))))))

(defun update-domain-corner-2d (mesh mp dt)
  "Use a corner tracking scheme to update domain lengths"
  (with-accessors ((position cl-mpm/particle::mp-position)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0)
                   )
      mp
    (with-accessors ((mesh-size cl-mpm/mesh::mesh-mesh-size))
        mesh
        (let ((diff (make-array 2 :initial-element 0d0))
              (domain-storage (magicl::matrix/double-float-storage domain)))
          (array-operations/utilities:nested-loop (x y) '(2 2)
            (let ((corner (cl-mpm/utils:vector-zeros))
                  (disp (cl-mpm/utils:vector-zeros)))
              (cl-mpm/fastmath::fast-.+-vector
               position
               (magicl:scale!
                (magicl:.*
                 (vector-from-list (mapcar (lambda (x) (- (* 2d0 (coerce x 'double-float)) 1d0)) (list x y 0)))
                 domain
                 ) 0.5d0) corner)
              (loop for i from 0 to 1
                    do (setf (the double-float (magicl:tref corner i 0))
                             (max 0d0 (min
                                       (the double-float (coerce (nth i mesh-size) 'double-float))
                                       (the double-float (magicl:tref corner i 0))))))
              (iterate-over-neighbours-point-linear-simd
               mesh corner
               (lambda (mesh node svp grads)
                 (declare (double-float dt svp))
                 (with-accessors ((vel cl-mpm/mesh:node-velocity))
                     node
                   (cl-mpm/fastmath:fast-fmacc disp vel (* dt svp)))))

              (incf (the double-float (aref diff 0)) (* 1.0d0 (the double-float (magicl:tref disp 0 0)) (- (* 2d0 (coerce x 'double-float)) 1d0)))
              (incf (the double-float (aref diff 1)) (* 1.0d0 (the double-float (magicl:tref disp 1 0)) (- (* 2d0 (coerce y 'double-float)) 1d0)))
              ))
          (incf (the double-float (aref domain-storage 0)) (* 0.5d0 (the double-float (aref diff 0))))
          (incf (the double-float (aref domain-storage 1)) (* 0.5d0 (the double-float (aref diff 1))))
          (let* ((jf  (the double-float (magicl:det def)))
                 (jl  (* (the double-float (magicl:tref domain 0 0))
                         (the double-float (magicl:tref domain 1 0))))
                 (jl0 (* (the double-float (magicl:tref domain-0 0 0))
                         (the double-float (magicl:tref domain-0 1 0))))
                 (scaling (the double-float
                               (expt (the double-float (/ (the double-float (* (the double-float jf) (the double-float jl0))) (the double-float jl)))
                                     (the double-float (/ 1d0 2d0))))))
            (setf (magicl:tref domain 0 0) (* (the double-float (magicl:tref domain 0 0)) scaling)
                  (magicl:tref domain 1 0) (* (the double-float (magicl:tref domain 1 0)) scaling)
                  ))))))

(defun update-domain-corner-3d (mesh mp dt)
  "Use a corner tracking scheme to update domain lengths"
  (with-accessors ((position cl-mpm/particle::mp-position)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0)
                   )
      mp
    (with-accessors ((mesh-size cl-mpm/mesh::mesh-mesh-size))
        mesh
        (let ((diff (make-array 3 :initial-element 0d0))
              (domain-storage (magicl::matrix/double-float-storage domain)))
          (array-operations/utilities:nested-loop (x y z) '(2 2 2)
            (let ((corner (cl-mpm/utils:vector-zeros))
                  (disp (cl-mpm/utils:vector-zeros)))
              (cl-mpm/fastmath::fast-.+-vector position
                                    (magicl:scale!
                                     (magicl:.*
                                      (vector-from-list (mapcar (lambda (x) (- (* 2d0 (coerce x 'double-float)) 1d0)) (list x y z)))
                                      domain
                                      ) 0.5d0) corner)
              (loop for i from 0 to 2
                    do (setf (the double-float (magicl:tref corner i 0))
                             (max 0d0 (min
                                       (the double-float (coerce (nth i mesh-size) 'double-float))
                                       (the double-float (magicl:tref corner i 0))))))
              (iterate-over-neighbours-point-linear
               mesh corner
               (lambda (mesh node svp grads)
                 (declare (double-float dt svp))
                 (with-accessors ((vel cl-mpm/mesh:node-velocity))
                     node
                   (cl-mpm/fastmath:fast-fmacc disp vel (* dt svp)))))

              (incf (the double-float (aref diff 0)) (* (the double-float (magicl:tref disp 0 0)) (- (* 2d0 (coerce x 'double-float)) 1d0)))
              (incf (the double-float (aref diff 1)) (* (the double-float (magicl:tref disp 1 0)) (- (* 2d0 (coerce y 'double-float)) 1d0)))
              (incf (the double-float (aref diff 2)) (* (the double-float (magicl:tref disp 2 0)) (- (* 2d0 (coerce z 'double-float)) 1d0)))
              ))
          (let ((nd (the fixnum (cl-mpm/mesh:mesh-nd mesh))))
            (incf (the double-float (aref domain-storage 0)) (* 0.5d0 (the double-float (aref diff 0))))
            (incf (the double-float (aref domain-storage 1)) (* 0.5d0 (the double-float (aref diff 1))))
            (incf (the double-float (aref domain-storage 2)) (* 0.5d0 (the double-float (aref diff 2))))
            (let* ((jf  (magicl:det def))
                   (jl  (* (magicl:tref domain 0 0) (magicl:tref domain 1 0) (magicl:tref domain 2 0)))
                   (jl0 (* (magicl:tref domain-0 0 0) (magicl:tref domain-0 1 0) (magicl:tref domain-0 2 0)))
                   (scaling (expt (/ (* jf jl0) jl) (/ 1d0 3d0))))
              (setf (magicl:tref domain 0 0) (* (magicl:tref domain 0 0) scaling)
                    (magicl:tref domain 1 0) (* (magicl:tref domain 1 0) scaling)
                    (magicl:tref domain 2 0) (* (magicl:tref domain 2 0) scaling)
                    ))
            )
          ))))
(defun update-domain-corner (mesh mp dt)
  (if (= (the fixnum (cl-mpm/mesh:mesh-nd mesh)) 2)
      (update-domain-corner-2d mesh mp dt)
      (update-domain-corner-3d mesh mp dt)
      )
  )


(defun update-stress-kirchoff (mesh mp dt fbar)
  "Update stress for a single mp"
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt)
           (optimize (speed 3) (safety 0) (debug 0)))
  (with-accessors ((stress cl-mpm/particle:mp-stress)
                   (stress-kirchoff cl-mpm/particle::mp-stress-kirchoff)
                   (volume cl-mpm/particle:mp-volume)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   ) mp
    (declare (magicl:matrix/double-float stress stress-kirchoff strain def strain-rate))
    (progn
      ;;   ;;For no FBAR we need to update our strains
      (progn
        ;; (unless fbar)
        (calculate-strain-rate mesh mp dt)

        ;; (loop for v across (cl-mpm/utils::fast-storage stretch-tensor)
        ;;       do (when (or (sb-ext::float-nan-p v)
        ;;                    (> v 1d5)
        ;;                    (< v -1d5)
        ;;                    )
        ;;            (error "Bad stretch tensor found ~A ~A" mp stretch-tensor)))

        ;; Turn cauchy stress to kirchoff
        ;; (setf stress stress-kirchoff)
        (cl-mpm/utils::voigt-copy-into stress-kirchoff stress)
        ;; (loop for v across (cl-mpm/utils::fast-storage strain)
        ;;       do (when (sb-ext::float-nan-p v)
        ;;            (error "PRE NaN strain found ~A" mp)))
        ;; (pprint stretch-tensor)
        ;; Update our strains
        (update-strain-kirchoff mesh mp dt fbar)
        ;; (loop for v across (cl-mpm/utils::fast-storage strain)
        ;;       do (when (sb-ext::float-nan-p v)
        ;;            (error "POST NaN strain found ~A" mp)))
        ;; Update our kirchoff stress with constitutive model
        ;; (setf stress-kirchoff (cl-mpm/particle:constitutive-model mp strain dt))
        (cl-mpm/utils::voigt-copy-into (cl-mpm/particle:constitutive-model mp strain dt) stress-kirchoff)

        ;; Check volume constraint!
        (when (<= volume 0d0)
          (error "Negative volume"))
        ;; Turn kirchoff stress to cauchy
        (cl-mpm/utils::voigt-copy-into stress-kirchoff stress)
        (cl-mpm/fastmath::fast-scale stress (/ 1.0d0 (the double-float (magicl:det def))))
        ;; (setf stress (magicl:scale stress-kirchoff (/ 1.0d0 (the double-float (magicl:det def)))))
        ))))

(declaim (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle double-float boolean) (values)) update-stress-kirchoff-damaged))
(defun update-stress-kirchoff-damaged (mesh mp dt fbar)
  "Update stress for a single mp"
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
      ;; (unless fbar)
      (calculate-strain-rate mesh mp dt)

      ;; Turn cauchy stress to kirchoff
      (cl-mpm/utils::voigt-copy-into stress-kirchoff stress)
      ;; Update our strains
      (update-strain-kirchoff-damaged mesh mp dt fbar)
      ;; Update our kirchoff stress with constitutive model
      (cl-mpm/utils::voigt-copy-into (cl-mpm/particle:constitutive-model mp strain dt) stress-kirchoff)
      ;; Check volume constraint!
      (when (<= volume 0d0)
        (error "Negative volume"))
      ;; Turn kirchoff stress to cauchy
      (cl-mpm/fastmath::fast-scale stress (/ 1.0d0 (the double-float (magicl:det def))))
      )))

(defun update-stress-linear (mesh mp dt fbar)
  "Update stress for a single mp"
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

        ;; Update our strains
        (update-strain-linear mesh mp dt fbar)

        ;; Update our kirchoff stress with constitutive model
        (setf stress (cl-mpm/particle:constitutive-model mp strain dt))
        ;; Check volume constraint!
        (when (<= volume 0d0)
          (error "Negative volume"))
        ))))

(defgeneric update-stress-mp (mesh mp dt fbar)
  (:documentation "A mp dependent stress update scheme"))

(defmethod update-stress-mp (mesh (mp cl-mpm/particle::particle) dt fbar)
  (update-stress-kirchoff mesh mp dt fbar)
  ;; (update-stress-linear mesh mp dt fbar)
  )

;; (defun calculate-cell-deformation (mesh cell dt)
;;   (with-accessors ((def cl-mpm/mesh::cell-deformation-gradient)
;;                    (volume cl-mpm/mesh::cell-volume))
;;       cell
;;     (let ((dstrain (magicl:zeros '(3 1))))
;;       (cl-mpm/mesh::cell-iterate-over-neighbours
;;        mesh cell
;;        (lambda (mesh c node svp grads)
;;          (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
;;                           (node-active cl-mpm/mesh:node-active)
;;                           )
;;              node
;;            (declare (double-float))
;;            (when node-active
;;              (mult (cl-mpm/shape-function::assemble-dsvp-2d grads) node-vel dstrain))
;;            )))
;;       (magicl:scale! dstrain dt)
;;       (let ((df (cl-mpm/fastmath::fast-.+ (magicl:eye 2) (voight-to-matrix dstrain))))
;;         (progn
;;           (setf def (magicl:@ df def))
;;           (setf volume (* volume (det df)))
;;           )))))
(defun map-jacobian (mesh mp dt)
  (with-accessors ((stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (volume cl-mpm/particle:mp-volume)
                   (def cl-mpm/particle:mp-deformation-gradient))
      mp
    (let* ((df (cl-mpm/utils::matrix-from-list '(1d0 0d0 0d0
                                                 0d0 1d0 0d0
                                                 0d0 0d0 1d0))))
      (cl-mpm/fastmath::fast-.+ df stretch-tensor df)
      (let ((j-inc (magicl:det df))
            (j-n (magicl:det def)))
        (iterate-over-neighbours
         mesh mp
         (lambda (mesh mp node svp grads fsvp fgrads)
           (with-accessors ((node-active cl-mpm/mesh:node-active)
                            (node-j-inc cl-mpm/mesh::node-jacobian-inc)
                            (node-volume cl-mpm/mesh::node-volume)
                            (node-lock cl-mpm/mesh::node-lock)
                            )
               node
             (declare (double-float node-j-inc svp volume j-inc j-n node-volume))
             (when node-active
               (sb-thread:with-mutex (node-lock)
                 (incf node-j-inc (/ (* svp volume j-inc j-n) node-volume)))))))))))

(declaim (inline calculate-df)
         (ftype (function (cl-mpm/mesh::mesh
                           cl-mpm/particle:particle
                           boolean) magicl:matrix/double-float)
                calculate-df))
(defun calculate-df (mesh mp fbar)
  (with-accessors ((dstrain cl-mpm/particle::mp-strain-rate)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (stretch-tensor-fbar cl-mpm/particle::mp-stretch-tensor-fbar)
                   (jfbar cl-mpm/particle::mp-j-fbar)
                   (def cl-mpm/particle:mp-deformation-gradient)
                   (pos cl-mpm/particle:mp-position)
                   )
      mp
    (let* ((df (cl-mpm/utils::matrix-from-list '(1d0 0d0 0d0
                                                 0d0 1d0 0d0
                                                 0d0 0d0 1d0))))
      (cl-mpm/fastmath::fast-.+-matrix df stretch-tensor df)
      ;;Coombs fbar
      (when fbar
        (let* ((df-fbar (cl-mpm/utils::matrix-from-list '(1d0 0d0 0d0
                                                          0d0 1d0 0d0
                                                          0d0 0d0 1d0)))
               (nd (cl-mpm/mesh::mesh-nd mesh)))
          (cl-mpm/fastmath::fast-.+-matrix df-fbar stretch-tensor-fbar df-fbar)
          (setf (cl-mpm/particle::mp-debug-j mp) (magicl:det df)
                (cl-mpm/particle::mp-debug-j-gather mp) (magicl:det df-fbar))
          (magicl:scale! df (expt
                             (the double-float (/ (magicl:det df-fbar)
                                                  (magicl:det df)))
                             (the double-float (/ 1d0 (float nd)))))
          (when (= nd 2)
            (setf (magicl:tref df 2 2) 1d0))
          ))

      ;; (when fbar
      ;;   (let ((j-inc (magicl:det df))
      ;;         (j-n (magicl:det def))
      ;;         (j-n1 0d0)
      ;;         (wsum 0d0)
      ;;         )
      ;;     (declare (double-float j-inc j-n j-n1))
      ;;     (iterate-over-neighbours
      ;;      mesh mp
      ;;      (lambda (mesh mp node svp grads f fb)
      ;;        (declare (ignore mesh mp  grads f fb))
      ;;        (with-accessors ((node-active cl-mpm/mesh:node-active)
      ;;                         (node-j-inc cl-mpm/mesh::node-jacobian-inc))
      ;;            node
      ;;          (declare (double-float svp node-j-inc))
      ;;          (when node-active
      ;;            (incf j-n1 (* svp node-j-inc))
      ;;            (incf wsum svp)
      ;;            ))))
      ;;     (setf j-n1 (/ j-n1 wsum))
      ;;     (when (<= (/ j-n1 (* j-n j-inc)) 0d0)
      ;;       (error "Negative volume"))
      ;;     (let ((nd (cl-mpm/mesh:mesh-nd mesh)))
      ;;       (magicl:scale! df
      ;;                      (expt
      ;;                       (the double-float (/ j-n1 (* j-n j-inc)))
      ;;                       (the double-float (/ 1d0 (float nd)))))
      ;;       (when (= nd 2)
      ;;         (setf (magicl:tref df 2 2) 1d0))
      ;;       )
      ;;     (setf (cl-mpm/particle::mp-debug-j mp) j-inc
      ;;           (cl-mpm/particle::mp-debug-j-gather mp) (magicl:det df))
      ;;     ))
      df)))
(defgeneric post-stress-step (mesh mp dt))
(defmethod post-stress-step (mesh mp dt)
  )
(defmethod post-stress-step (mesh (mp cl-mpm/particle::particle) dt)
  (declare (ignore mesh mp dt)))

(declaim (ftype (function (cl-mpm/mesh::mesh
                           (array cl-mpm/particle:particle)
                           double-float
                           &optional boolean) (values)) update-stress))
(defun update-stress (mesh mps dt &optional (fbar nil))
  "Update all stresses, with optional f-bar"
  (declare ((array cl-mpm/particle:particle) mps) (cl-mpm/mesh::mesh mesh))
  ;;Update stress
  ;; (when fbar
  ;;   (cl-mpm:iterate-over-mps
  ;;    mps
  ;;    (lambda (mp)
  ;;      (calculate-strain-rate mesh mp dt)
  ;;      (map-jacobian mesh mp dt))))

  (iterate-over-mps
   mps
   (lambda (mp)
     (update-stress-mp mesh mp dt fbar)
     (post-stress-step mesh mp dt)
     ))

  ;; (lparallel:pdotimes (i (length mps))
  ;;   (update-stress-mp mesh (aref mps i) dt fbar)
  ;;   (post-stress-step mesh (aref mps i) dt)
  ;;   )
  ;; (dotimes (i (length mps))
  ;;   (update-stress-mp mesh (aref mps i) dt fbar)
  ;;   (post-stress-step mesh (aref mps i) dt)
  ;;   )
  (values))

(declaim (notinline reset-grid))
(defun reset-grid (mesh)
  "Reset all nodes on the grid"
  (declare (cl-mpm/mesh::mesh mesh))
  (iterate-over-nodes
   mesh
   (lambda (node)
     (when (cl-mpm/mesh:node-active node)
       (cl-mpm/mesh:reset-node node)))))
(declaim (notinline reset-grid-velocity))
(defun reset-grid-velocity (mesh)
  "Reset all velocity map on grid for MUSL"
  (declare (cl-mpm/mesh::mesh mesh))
  (iterate-over-nodes
   mesh
   (lambda (node)
     (when (cl-mpm/mesh:node-active node)
       (cl-mpm/mesh::reset-node-velocity node)))))


(defun filter-grid (mesh mass-thresh)
  "Filter out nodes with very small massess"
  (declare (cl-mpm/mesh::mesh mesh) (double-float mass-thresh))
  (iterate-over-nodes
   mesh
   (lambda (node)
     (when (and (cl-mpm/mesh:node-active node)
                (<= (cl-mpm/mesh:node-mass node) mass-thresh))
       (cl-mpm/mesh:reset-node node)))))

(defun filter-grid-velocity (mesh mass-thresh)
  "Filter out nodes with very small massess"
  (declare (cl-mpm/mesh::mesh mesh) (double-float mass-thresh))
  (iterate-over-nodes
   mesh
   (lambda (node)
     (when (and (cl-mpm/mesh:node-active node)
                (<= (cl-mpm/mesh:node-mass node) mass-thresh))
       (cl-mpm/mesh::reset-node-velocity node)))))

(defun filter-grid-volume (mesh volume-ratio)
  (declare (cl-mpm/mesh::mesh mesh) (double-float volume-ratio))
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
    (dotimes (i (array-total-size nodes))
      (let ((node (row-major-aref nodes i)))
        (when (and (cl-mpm/mesh:node-active node)
                   (< (the double-float (cl-mpm/mesh:node-mass node))
                      (* volume-ratio (the double-float (cl-mpm/mesh::node-volume-true node)))))
          (setf (cl-mpm/mesh::node-active node) nil)
          (cl-mpm/mesh:reset-node node)
          )))))

(defgeneric sim-add-mp (sim mp)
  (:documentation "A function to append an mp to a simulation"))
(defmethod sim-add-mp (sim mp)
  (vector-push-extend mp (cl-mpm:sim-mps sim)))

(defgeneric remove-mps-func (sim func)
  (:documentation "A function for removing mps from a sim"))

(defmethod remove-mps-func (sim func)
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (declare ((vector cl-mpm/particle::particle) mps))
    (when (> (length mps) 0)
      (setf mps
            ;;We cant do this in parallel apparently
            (delete-if func mps)
                                        ;(lparallel:premove-if func mps)
            ))
    ;;Sometimes when compacting the array; sbcl will just discard make and unadjustable array in place which is a bit wild
    (when (and (not (adjustable-array-p mps))
               (= (length mps) 0))
      (setf mps (make-array 0 :adjustable t :fill-pointer 0)))
    )
  (values))

(defun remove-material-damaged (sim)
  "Remove material points that have strain energy density exceeding fracture toughness"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)
                   (instant-damage-removal cl-mpm::sim-mp-damage-removal-instant)
                   )
      sim
    (let ((h (cl-mpm/mesh:mesh-resolution mesh)))
      (remove-mps-func
       sim
       (lambda (mp)
         (with-accessors ((damage cl-mpm/particle:mp-damage)
                          (def cl-mpm/particle::mp-deformation-gradient))
             mp
           (and (>= damage 0.9d0)
                (or instant-damage-removal
                    (damage-removal-criteria mp h) )
                )))))
    ))

(defun damage-removal-criteria (mp h)
  "Some criteria for when we should remove fully damaged MPs"
  (with-accessors ((def cl-mpm/particle:mp-deformation-gradient)
                   (lens cl-mpm/particle::mp-domain-size)
                   (strain cl-mpm/particle::mp-strain)
                   (lens-0 cl-mpm/particle::mp-domain-size-0)
                   ;; (ft cl-mpm/particle::mp-initiation-stress)
                   ;; (ductility cl-mpm/particle::mp-ductility)
                   ;; (E cl-mpm/particle::mp-e)
                   )
      mp
    ;; (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix strain)))
    (let* ((l-factor 2.00d0)
           ;; (e0 (/ ft E))
           ;; (ef (/ (* ft (+ ductility 1d0)) (* 2d0 E)))
           ;; (ef (* ef 50d0))
           ;; (e_max (reduce #' max l))
           )
      (cond
        ;; ((< ef e_max) :x)
        ;; ((< ef (tref strain 0 0)) :x)
        ;; ((< ef (tref strain 1 0)) :y)
        ;; ((< ef (tref strain 2 0)) :z)
        ((< (* l-factor (tref lens-0 0 0)) (tref lens 0 0)) :x)
        ((< (* l-factor (tref lens-0 1 0)) (tref lens 1 0)) :y)
        ;; ((and (< (* l-factor (tref lens-0 2 0)) (tref lens 2 0))
        ;;       (> (tref lens-0 2 0) 0d0)
        ;;       ) :z)
        (t nil)
        ))))

(defparameter *max-split-depth* 3)
(defun split-criteria (mp h)
  "Some numerical splitting estimates"
  (with-accessors ((def cl-mpm/particle:mp-deformation-gradient)
                   (lens cl-mpm/particle::mp-domain-size)
                   (lens-0 cl-mpm/particle::mp-domain-size-0)
                   (split-depth cl-mpm/particle::mp-split-depth)
                   )
      mp
    (when (< split-depth *max-split-depth*)
      (let ((l-factor 1.00d0)
            (h-factor (* 0.6d0 h))
            (s-factor 1.5d0))
        (cond
          ((< h-factor (tref lens 0 0)) :x)
          ((< h-factor (tref lens 1 0)) :y)
          ((< h-factor (tref lens 2 0)) :z)
          ;; ((< (* l-factor (tref lens-0 0 0)) (tref lens 0 0)) :x)
          ;; ((< (* l-factor (tref lens-0 1 0)) (tref lens 1 0)) :y)
          ;; ((< (* l-factor (tref lens-0 2 0)) (tref lens 2 0)) :z)
          ;; ((< l-factor (/ (tref lens 0 0) (tref lens-0 0 0))) t)
          ;; ((< l-factor (/ (tref lens 1 0) (tref lens-0 1 0))) t)
                                        ;((< 2.0d0 (tref def 1 1)) t)
                                        ;((> 1.5d0 (tref def 0 0)) t)
          (t nil)
          )))))

(defun single-particle-criteria (mesh mp)
  "Criteria for checking if material point is unconnected to another MP"
  (let ((svp-sum 0d0)
        (alone t))
    (iterate-over-neighbours
      mesh mp
      (lambda (mesh mp node svp grads fsvp fgrads)
        (declare
          (cl-mpm/mesh::node node)
          (cl-mpm/particle:particle mp)
          (type double-float svp))
        (with-accessors ((node-svp cl-mpm/mesh::node-svp-sum)
                         (node-active cl-mpm/mesh:node-active)
                         ) node
          (when node-active
            (incf svp-sum node-svp)
            ;; (when (> node-svp svp)
            ;;   (setf alone nil))
            )
          )))
    (and (< svp-sum 2d0) (not (= svp-sum 0d0)))
    ;; (setf alone t)
    
    ;; alone
    ))
(defun gimp-removal-criteria (mp h)
  "Criteria for removal of gimp mps based on domain length"
  (with-accessors ((lens cl-mpm/particle::mp-domain-size))
      mp
    (declare (double-float h))
    (let ((h-factor (* 1.0d0 h))
          (aspect 0.01d0))
      (cond
        ((< h-factor (the double-float (tref lens 0 0))) :x)
        ((< h-factor (the double-float (tref lens 1 0))) :y)
        ((< h-factor (the double-float (tref lens 2 0))) :z)
        ((> aspect (/ (the double-float (tref lens 0 0))
                      (the double-float (tref lens 1 0))
                      )) :x)
        ((> aspect (/ (the double-float (tref lens 1 0))
                      (the double-float (tref lens 0 0))
                      )) :x)
        (t nil)
        ))))
(defun copy-particle (original &rest initargs &key &allow-other-keys)
  "Function for copying a particle, allowing for re-initialising members"
  (let* ((class (class-of original))
         (copy (allocate-instance class)))
    (dolist (slot (mapcar #'sb-mop:slot-definition-name (sb-mop:class-slots class)))
      (when (slot-boundp original slot)
        (cond
          ((typep (slot-value original slot) 'magicl::abstract-tensor)
           (setf (slot-value copy slot)
                 (magicl:scale (slot-value original slot) 1d0)))
          (t (setf (slot-value copy slot)
                  (slot-value original slot))))))
    (apply #'reinitialize-instance copy initargs)))

(defmacro split-linear (dir direction dimension)
  "Helper macro for single splitting along cartesian directions "
  `((eq ,dir ,direction)
    (let ((new-size (vector-from-list (list 1d0 1d0 1d0)))
          (new-size-0 (vector-from-list (list 1d0 1d0 1d0)))
          (pos-offset (vector-zeros))
          (new-split-depth (+ (cl-mpm/particle::mp-split-depth mp) 1))
          )
      (setf (tref new-size ,dimension 0) 0.5d0)
      (setf (tref new-size-0 ,dimension 0) 0.5d0)
      (setf (tref pos-offset ,dimension 0) 0.25d0)
      (cl-mpm/fastmath::fast-.* lens new-size new-size)
      (cl-mpm/fastmath::fast-.* lens-0 new-size-0 new-size-0)
      (cl-mpm/fastmath::fast-.* lens pos-offset pos-offset)
      (list
       (copy-particle mp
                      :mass (/ mass 2)
                      :volume (/ volume 2)
                      :size (cl-mpm/utils::vector-copy new-size)
                      :size-0 (cl-mpm/utils::vector-copy new-size-0)
                      :position (cl-mpm/fastmath::fast-.+-vector pos pos-offset)
                      :nc (make-array 8 :fill-pointer 0 :element-type 'node-cache)
                      :split-depth new-split-depth
                       )
       (copy-particle mp
                      :mass (/ mass 2)
                      :volume (/ volume 2)
                      :size (cl-mpm/utils::vector-copy new-size)
                      :size-0 (cl-mpm/utils::vector-copy new-size-0)
                      :position (magicl:.- pos pos-offset)
                      :nc (make-array 8 :fill-pointer 0 :element-type 'node-cache)
                      :split-depth new-split-depth
                      )))))
(defmacro split-cases (direction)
  "Another helper macro for splitting mps"
  `(cond
     ,(macroexpand-1 '(split-linear direction :x 0))
     ,(macroexpand-1 '(split-linear direction :y 1))
     ,(macroexpand-1 '(split-linear direction :z 2))
     (t nil)
     ))
(defun split-mp (mp h direction)
  "Function to split an mp across a direction
   Directions should be :x,:y,:z "
  (with-accessors ((def cl-mpm/particle:mp-deformation-gradient)
                   (lens cl-mpm/particle::mp-domain-size)
                   (lens-0 cl-mpm/particle::mp-domain-size-0)
                   (pos cl-mpm/particle:mp-position)
                   (mass cl-mpm/particle:mp-mass)
                   (volume cl-mpm/particle:mp-volume)
                   ;; (volume cl-mpm/particle::mp-volume-0)
                   )
      mp
    (let ((l-factor 1.50d0)
          (h-factor (* 0.8d0 h)))
      (split-cases direction))))
(defun split-mps (sim)
  "Split mps that match the split-criteria"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (mps-to-split (remove-if-not (lambda (mp) (split-criteria mp h)) mps))
           (split-direction (map 'list (lambda (mp) (split-criteria mp h)) mps-to-split)))
      ;; (setf mps (delete-if (lambda (mp) (split-criteria mp h)) mps))
      (remove-mps-func sim (lambda (mp) (split-criteria mp h)))
      (loop for mp across mps-to-split
            for direction in split-direction
            do (loop for new-mp in (split-mp mp h direction)
                     do (sim-add-mp sim new-mp)))
      )))

(defun split-mps-criteria (sim criteria)
  "Split mps that fail an arbritary criteria"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (mps-to-split (remove-if-not (lambda (mp) (funcall criteria mp h)) mps))
           (split-direction (map 'list (lambda (mp) (funcall criteria mp h)) mps-to-split)))
      ;; (setf mps (delete-if (lambda (mp) (funcall criteria mp h)) mps))
      (remove-mps-func sim criteria)
      (loop for mp across mps-to-split
            for direction in split-direction
            do (loop for new-mp in (split-mp mp h direction)
                     do (sim-add-mp sim new-mp)))
      )))

(defgeneric calculate-min-dt (sim)
  (:documentation "A function for calculating an approximate stable timestep"))
(defmethod calculate-min-dt ((sim mpm-sim))
  "Estimate minimum p-wave modulus"
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
                          (vel cl-mpm/mesh::node-velocity)
                          ) node
           (when (and node-active
                      (> vol 0d0)
                      (> pmod 0d0)
                      (> svp-sum 0d0))
             (let ((nf (+ (/ mass (* vol (+ (/ pmod svp-sum) ;; (* svp-sum (cl-mpm/fastmath::mag-squared vel))
                                            )))))) 
                 (when (< nf inner-factor)
                   (setf inner-factor nf)))))))
      (if (< inner-factor most-positive-double-float)
          (* (sqrt mass-scale) (sqrt inner-factor) (cl-mpm/mesh:mesh-resolution mesh))
          (cl-mpm:sim-dt sim)))))
(defun calculate-adaptive-time (sim target-time &key (dt-scale 1d0))
  "Given that the system has previously been mapped to, caluclate an estimated dt and substep for a given target time
The value dt-scale allows for the estimated dt to be scaled up or down
This modifies the dt of the simulation in the process
"
  (let* ((dt-e (* dt-scale (calculate-min-dt sim)))
         (substeps-e (floor target-time dt-e)))
    ;; (format t "CFL dt estimate: ~f~%" dt-e)
    ;; (format t "CFL step count estimate: ~D~%" substeps-e)
    (setf (cl-mpm:sim-dt sim) dt-e)
    ;; (setf substeps substeps-e)
    (values dt-e substeps-e)
    )
  )
(defgeneric add-mps (sim mps-array))
(defmethod add-mps (sim mps-array)
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
      (if (> (length mps) 0)
          (progn
            (loop for mp across mps-array
                  do (vector-push-extend mp mps (length mps-array))))
          (setf (cl-mpm:sim-mps sim) mps-array))))

#||
(progn
(ql:quickload :cl-mpm/examples/fracture)
(in-package :cl-mpm/examples/fracture))
||#




;; (sb-ext:restrict-compiler-policy 'speed  0 0)

