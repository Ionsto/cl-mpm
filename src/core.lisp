(defpackage :cl-mpm
  (:use :cl
        :cl-mpm/particle
        :cl-mpm/mesh
        :cl-mpm/utils
        :cl-mpm/fastmath
;        :cl-mpm/shape-function
        )
  (:import-from
    :magicl tref .+ .-)
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

(defun make-mpm-sim (size resolution dt shape-function &key (sim-type 'mpm-sim))
  (make-instance 'mpm-sim-usf
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
                    (update-node-forces mesh (sim-damping-factor sim) dt)
                    ;Apply vel bcs
                    (apply-bcs mesh bcs dt)
                    ;;G2p + Particle update
                    (g2p mesh mps dt)
                    (update-stress mesh mps dt)
                    (when enable-damage
                      (cl-mpm/damage::calculate-damage mesh
                                                       mps
                                                       dt
                                                       50d0))
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
               (remove-damage allow-mp-damage-removal)
               )
                sim
    (declare (type double-float mass-filter))
                (progn
                    (reset-grid mesh)
                    (p2g mesh mps)
                    (when (> mass-filter 0d0)
                     (filter-grid mesh (sim-mass-filter sim)))
                    (update-node-kinematics mesh dt)
                    (apply-bcs mesh bcs dt)
                    (update-stress mesh mps dt)
                    (when enable-damage
                     (cl-mpm/damage::calculate-damage mesh
                                                      mps
                                                      dt
                                                      50d0))
                    ;Map forces onto nodes
                    (p2g-force mesh mps)
                    (apply-bcs mesh bcs-force dt)
                    (update-node-forces mesh (sim-damping-factor sim) dt)
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
                    (update-node-forces mesh (sim-damping-factor sim) dt)
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
                    (when enable-damage
                      (cl-mpm/damage::calculate-damage mesh
                                                       mps
                                                       dt
                                                       50d0))
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

(defun create-node-cache (mesh mp func)
  (declare (function func))
  (with-accessors ((nodes cl-mpm/particle::mp-cached-nodes))
      mp
      (iterate-over-neighbours-shape-gimp
       mesh mp
       (lambda (mesh mp node svp grads)
         (vector-push-extend
          (cl-mpm/particle::make-node-cache
           :node node
           :weight svp
           :grads grads)
          nodes)
         (funcall func mesh mp node svp grads)
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
                  (cl-mpm/particle::node-cache-grads nc))))

(defmethod iterate-over-neighbours-shape (mesh (shape-func cl-mpm/shape-function:shape-function-linear) mp func)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/shape-function::shape-function-linear shape-func)
           (cl-mpm/particle:particle mp))
  (let* ((nd (cl-mpm/mesh:mesh-nd mesh))
         (h (cl-mpm/mesh:mesh-resolution mesh))
         ;(order (slot-value shape-func 'order))
         (pos-vec (cl-mpm/particle:mp-position mp));Material point position 
         (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
         (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec 'floor)))
    (loop for dx from 0 to 1
          do (loop for dy from 0 to 1
                   do
                   (let* ((id (mapcar #'+ pos-index (list dx dy))))
                     (when (cl-mpm/mesh:in-bounds mesh id)
                       (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                              (node (cl-mpm/mesh:get-node mesh id))
                              (weight (apply (cl-mpm/shape-function::svp shape-func) dist))
                              (grads (apply (cl-mpm/shape-function::dsvp shape-func) dist))) 
                         (funcall func mesh mp node weight (cl-mpm/shape-function::assemble-dsvp nd grads) grads))))))))

(defmethod iterate-over-neighbours-shape (mesh (shape-func cl-mpm/shape-function:shape-function-bspline) mp func) (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/shape-function::shape-function-bspline shape-func)
           (cl-mpm/particle:particle mp))
  (let* ((nd (cl-mpm/mesh:mesh-nd mesh))
         (h (cl-mpm/mesh:mesh-resolution mesh))
         (pos-vec (cl-mpm/particle:mp-position mp));Material point position 
         (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
         (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec #'round))
         (border (not (and (cl-mpm/mesh:in-bounds (mapcar #'+ pos-index (list 1 1)))
                          (cl-mpm/mesh:in-bounds (mapcar #'- pos-index (list 1 1)))
                          ))))
    (if nil ;border
        ;; Do something janky if we are in a border element
        (loop for dx from -1 to 1
              do (loop for dy from -1 to 1
                       do
                          (let* (
                                 (id (mapcar #'+ pos-index (list dx dy)))
                                 (dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                 (weight (apply (cl-mpm/shape-function:svp shape-func) dist))
                                 (grads (apply (cl-mpm/shape-function:dsvp shape-func) dist))
                                 )
                            (when (cl-mpm/mesh:in-bounds mesh id)
                              (funcall func mesh mp (cl-mpm/mesh:get-node mesh id) weight
                                       (cl-mpm/shape-function:assemble-dsvp nd grads)))))))))
;; (declaim (inline dispatch)
;;          (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle function) (values)) iterate-over-neighbours-shape-linear))
(defmacro iterate-over-neighbours-general (mesh mp order body)
  `(progn
    (let* ((h (cl-mpm/mesh:mesh-resolution ,mesh))
           (pos-vec (cl-mpm/particle:mp-position ,mp))
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
           (pos-index (magicl:scale pos-vec (/ 1 h)))
           )
      (loop for dx from (- ,order) to ,order
            do (loop for dy from (- ,order) to ,order
                     do (let* ((id (list (+ (round (/ (tref pos-index 0 0) h)) dx)
                                         (+ (round (/ (tref pos-index 1 0) h)) dy))))
                          (when (cl-mpm/mesh:in-bounds ,mesh id)
                            (funcall ,body
                                     ,mesh
                                     ,mp
                                     id
                                     (mapcar (lambda (p i) (- p (* i h))) pos id)))))))))
(defmacro iterate-over-neighbours-general-pure (mesh mp order args body)
  (let ((h (gensym))
        (pos-vec (gensym))
        (pos (gensym))
        (pos-index (gensym))
        (id (gensym))
        (dx (gensym))
        (dy (gensym))
        (dist (gensym)))
    `(progn
       (let* ((,h (cl-mpm/mesh:mesh-resolution ,mesh))
              (,pos-vec (cl-mpm/particle:mp-position ,mp))
              (,pos (list (tref ,pos-vec 0 0) (tref ,pos-vec 1 0)))
              (,pos-index (magicl:scale ,pos-vec (/ 1 ,h)))
              )
         (loop for ,dx from (- ,order) to ,order
               do (loop for ,dy from (- ,order) to ,order
                        do (let* ((,id (list (+ (round (tref ,pos-index 0 0)) ,dx)
                                             (+ (round (tref ,pos-index 1 0)) ,dy))))
                             (when (cl-mpm/mesh:in-bounds ,mesh ,id)
                               (let ((,dist (mapcar (lambda (p i) (- p (* i ,h))) ,pos ,id)))
                                 (symbol-macrolet ((,(first args) ,mesh)
                                                   (,(second args) ,mp)
                                                   (,(third args) ,id)
                                                   (,(fourth args) ,dist))
                                   ,body))))))))))
(defmacro ion-linear (number body)
  (let ((test (expt number 2)))
    `(+ ,test ,body))
  ;; `(with-gensyms (exp-number-name)
  ;;    (let (exp-number-name ,)))
  )
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
                                (funcall func mesh mp node weight grads))))))))))

(declaim (inline iterate-over-neighbours-shape-gimp)
         (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle function) (values))
                iterate-over-neighbours-shape-gimp))
(defun iterate-over-neighbours-shape-gimp (mesh mp func)
  (declare (type cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp)
           (function func))
  (progn
    (with-accessors ((pos-vec cl-mpm/particle:mp-position)
                     (d0 cl-mpm/particle::mp-domain-size))
        mp
      (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
             (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
             (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec))
             (domain (loop for x across (magicl::storage d0) collect (* 0.5d0 (the double-float x))))
             )
        (declare (type double-float h)
                 (type list pos pos-index domain))
        (loop for dx from -1 to 1
              do (loop for dy from -1 to 1
                       do (let* ((id (mapcar #'+ pos-index (list dx dy))))
                            (when (cl-mpm/mesh:in-bounds mesh id)
                              (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                     (weights (mapcar (lambda (x l)
                                                        (cl-mpm/shape-function::shape-gimp x l h))
                                                      dist domain))
                                     (weight (reduce #'* weights))
                                     )
                                (when (< 0d0 weight)
                                  (let* ((node (cl-mpm/mesh:get-node mesh id))
                                         (grads (mapcar (lambda (d l w)
                                                          (* (cl-mpm/shape-function::shape-gimp-dsvp d l h) w))
                                                        dist domain (nreverse weights))))
                                    (funcall func mesh mp node weight grads))
                                )
                                )))))
        ))))

;; (declaim (inline iterate-over-neighbours-shape-gimp)
;;          (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle function) (values))
;;                 iterate-over-neighbours-shape-gimp))
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
;;              (domain (loop for x across (magicl::storage d0) collect (* 0.5d0 (the double-float x))))
;;              (pminx (floor (- (tref pos-vec 0 0)   (* 0.5d0 (tref d0 0 0)))))
;;              (pminy (floor (- (tref pos-vec 1 0)   (* 0.5d0 (tref d0 1 0)))))
;;              (pmaxx (ceiling (+ (tref pos-vec 0 0) (* 0.5d0 (tref d0 0 0)))))
;;              (pmaxy (ceiling (+ (tref pos-vec 1 0) (* 0.5d0 (tref d0 1 0)))))
;;              )
;;         (loop for dx from pminx to pmaxx
;;               do (loop for dy from pminy to pmaxy
;;                        do (let* ((id (mapcar #'+ pos-index (list dx dy))))
;;                             (when (cl-mpm/mesh:in-bounds mesh id)
;;                               (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
;;                                      (node (cl-mpm/mesh:get-node mesh id))
;;                                      (weights (mapcar (lambda (x l)
;;                                                         (cl-mpm/shape-function::shape-gimp x l h))
;;                                                       dist domain))
;;                                      (weight (reduce #'* weights))
;;                                      (grads (mapcar (lambda (d l w)
;;                                                       (* (cl-mpm/shape-function::shape-gimp-dsvp d l h) w))
;;                                                     dist domain (nreverse weights)))
;;                                      )
;;                                 (when (< 0d0 weight)
;;                                   (funcall func mesh mp node weight grads)
;;                                 )
;;                                 )))))
;;         ))))

(defun make-knot-list (mesh pos)
  (list
    (loop for dy from -1 to 1
          collect
          (loop for dx from -4 to 4
                collect (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos (list dx dy)))))
    (loop for dx from -1 to 1
          collect
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
      (if nil;border
          (let (
                ;; (knots-x
                ;;   (loop for dy from -1 to 1
                ;;         collect
                ;;         (loop for dx from -4 to 4
                ;;                          collect (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list dx dy))))))
                ;; (knots-y
                ;;   (loop for dx from -1 to 1
                ;;         collect
                ;;         (loop for dy from -4 to 4
                ;;               collect (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list dx dy))))))
                )
            (loop for dx from -1 to 1
                  do (loop for dy from -1 to 1
                           do (let* ((id (mapcar #'+ pos-index (list dx dy))))
                                (when (cl-mpm/mesh:in-bounds mesh id)
                                  ;; Ill defined bspline
                                  (when t
                                      ;; (and (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list 0 dy)))
                                      ;;        (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list dx 0)))
                                      ;;        )
                                    (let* (
                                           ;; (knot-list (list (nth (+ dy 1) knots-x)
                                           ;;                  (nth (+ dx 1) knots-y)))
                                           (knot-list (make-knot-list mesh id))
                                           (local-id-list (list (+ dx 0) (+ dy 0)))
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
                                      (funcall func mesh mp node weight grads))))))))
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
                                  (funcall func mesh mp node weight grads))))))))))
(defmacro iterate-over-neighbours-shape-linear-macro (mesh mp func-body)
  `(progn
    (let* ((h (cl-mpm/mesh:mesh-resolution ,mesh))
           (pos-vec (cl-mpm/particle:mp-position ,mp))
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
           (pos-index (magicl:scale pos-vec (/ 1 h)))
           ;; (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec 'floor)))
           )
      (loop for dx from 0 to 1
            do (loop for dy from 0 to 1
                     do (let* (
                               (id (list (+ (floor (tref pos-index 0 0)) dx)
                                         (+ (floor (tref pos-index 1 0)) dy)
                                         ))
                               (dist (mapcar (lambda (p i) (- p (* i h))) pos id))
                               (node (cl-mpm/mesh:get-node ,mesh id))
                               (weights (mapcar #'shape-linear dist))
                               (weight (reduce #'* weights))
                               (dsvp (mapcar (lambda (d w) (* (shape-linear-dsvp d) w)) dist (nreverse weights)))
                               ;; (dsvp (assemble-dsvp-2d (list (* (shape-linear-dsvp (first dist)) wy)
                               ;;                               (* (shape-linear-dsvp (second dist)) wx))
                               ;;        ))
                               )
                          ,func-body
                          )
                     )))))


;(defparameter *mesh* (sim-mesh *sim*))
;(defparameter *sf* (cl-mpm/mesh:mesh-shape-func *mesh*))
;(defparameter *mp* (cl-mpm/particle:make-particle 2 :pos '(0.5 0.1)))
;(iterate-over-neighbours-shape *mesh* *sf* *mp* 
;                               (lambda (mesh mp node w dsvp) (print (cl-mpm/mesh:node-index  node))))

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
                   ;; (strain-rate cl-mpm/particle:mp-strain-rate)
                   ) mp
    (declare (type double-float mp-mass mp-volume))
    (iterate-over-neighbours
     mesh mp
     (lambda (mesh mp node svp grads)
       (declare
        (cl-mpm/particle:particle mp)
        (cl-mpm/mesh::node node)
        (double-float svp)
        )
       (with-accessors ((node-vel   cl-mpm/mesh:node-velocity)
                        (node-active  cl-mpm/mesh::node-active)
                        (node-mass  cl-mpm/mesh:node-mass)
                        (node-volume  cl-mpm/mesh::node-volume)
                        (node-force cl-mpm/mesh:node-force)
                        (node-lock  cl-mpm/mesh:node-lock)) node
         (declare (type double-float node-mass node-volume mp-volume)
                  (type sb-thread:mutex node-lock))
         (sb-thread:with-mutex (node-lock)
           (setf node-active t)
           (incf node-mass
                 (* mp-mass svp))
           (incf node-volume
                 (* mp-volume svp))
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

(declaim (inline p2g-force-mp)
         (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle) (values)) p2g-force-mp)
         )
(defun p2g-force-mp (mesh mp)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (with-accessors ((mp-vel  cl-mpm/particle:mp-velocity)
                   (mp-mass cl-mpm/particle:mp-mass)
                   ) mp
    (declare (type double-float mp-mass))
    (iterate-over-neighbours
     mesh mp
     (lambda (mesh mp node svp grads)
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
             (det-int-force mp (cl-mpm/shape-function::assemble-dsvp-2d grads) node-force)
             ))
          )
        )))
  (values))

(declaim (notinline p2g))
(defun p2g-force (mesh mps)
  (declare (type (array cl-mpm/particle:particle) mps) (cl-mpm/mesh::mesh mesh))
  (lparallel:pdotimes (i (length mps)) 
    (p2g-force-mp mesh (aref mps i))))

(defgeneric special-g2p (mesh mp node svp dsvp grads)
  (:documentation "G2P behaviour for specific features")
  (:method (mesh mp node svp dsvp grads)))

(defmethod special-g2p (mesh (mp cl-mpm/particle::particle-thermal) node svp dsvp grads) 
  (declare (ignore mesh))
  "Map grid to particle for one mp-node pair"
  (with-accessors ((node-temp cl-mpm/mesh:node-temperature)) node
    (with-accessors ((temp cl-mpm/particle::mp-temperature)) mp
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
  "Map one MP onto the grid"
  (with-accessors ((mass mp-mass)
                   (vel mp-velocity)
                   (pos mp-position)
                   (disp cl-mpm/particle::mp-displacement)
                   (acc cl-mpm/particle::mp-acceleration)
                   (temp cl-mpm/particle::mp-temperature)
                   (strain-rate mp-strain-rate)
                   (vorticity cl-mpm/particle:mp-vorticity)
                   (nc cl-mpm/particle::mp-cached-nodes)
                   )
    mp
    (let* ((mapped-vel (make-array 2 :initial-element 0d0 :element-type 'double-float))
          (mapped-acc (make-array 2 :initial-element 0d0 :element-type 'double-float))
          ;; (mapped-vel-mat (magicl::from-storage mapped-vel '(2 1)))
           )
      (progn
        ;; (magicl:scale vel 0d0)
        ;; (reset-mps-g2p mp)
        )
      ;; Map variables
      (iterate-over-neighbours
       mesh mp
       (lambda (mesh mp node svp grads)
         (declare
          (cl-mpm/mesh::node node)
          (cl-mpm/particle:particle mp)
          (type double-float svp))
         (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
                          (node-acc cl-mpm/mesh:node-acceleration)
                          (node-active cl-mpm/mesh:node-active)
                          ) node
           (when node-active
               (cl-mpm/fastmath::simd-fmacc mapped-vel (magicl::storage node-vel) svp)
               (cl-mpm/fastmath::simd-fmacc mapped-acc (magicl::storage node-acc) svp)
             )
           )
         ;; (g2p-mp-node mp node svp grads)
         ))
      ;;Update particle
      (progn
        (setf (fill-pointer nc) 0)
        (cl-mpm/fastmath::fast-fmacc-array (magicl::storage pos) mapped-vel dt)
        (cl-mpm/fastmath::fast-fmacc-array (magicl::storage disp) mapped-vel dt)
        (aops:copy-into (magicl::storage acc) mapped-acc)
        ;;FLIP
        (cl-mpm/fastmath::simd-fmacc (magicl::storage vel)  mapped-acc dt)
        ;;Direct velocity damping
        ;; (magicl:scale! vel (- 1d0 1d-3))
        ;;PIC
        ;; (aops:copy-into (magicl::storage vel) mapped-vel)
        ;; (update-domain mesh mp dt)
        )
        ;; (setf pos (magicl:.+ pos (magicl:scale vel dt)))
        ;; (setf disp (magicl:.+ disp (magicl:scale vel dt))))
      )))

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
          (magicl:scale! vel (/ 1.0d0 mass))))))

(declaim (inline calculate-forces)
         (ftype (function (cl-mpm/mesh::node double-float double-float) (vaules)) calculate-forces))
(defun calculate-forces (node damping dt)
  (when (cl-mpm/mesh:node-active node)
      (with-accessors ( (mass  node-mass)
                        (vel   node-velocity)
                        (force node-force)
                        (acc   node-acceleration))
        node
        (declare (double-float mass dt damping))
        (progn
          (magicl:scale acc 0d0)
          (cl-mpm/fastmath:fast-fmacc acc force (/ 1d0 mass))
          (cl-mpm/fastmath:fast-fmacc acc vel (* damping -1d0))
          ;;Disable for FLIP - not true?
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
(defun update-node-forces (mesh damping dt)
  (iterate-over-nodes mesh
                      (lambda (node)
                        (calculate-forces node damping dt))))

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
                   (velocity-rate cl-mpm/particle::mp-velocity-rate)
                   ) mp
        (progn
          (magicl:scale! strain-rate 0d0)
          (magicl:scale! vorticity 0d0)
;          (setf (cl-mpm/particle::mp-cached-nodes mp) '())
            (iterate-over-neighbours mesh mp
                (lambda (mesh mp node svp grads)
                  (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
                                   (node-active cl-mpm/mesh:node-active)
                                   )
                      node
                    (declare (double-float))
                    (when node-active
                      (mult (cl-mpm/shape-function::assemble-dsvp-2d grads) node-vel strain-rate)
                      (mult (cl-mpm/shape-function::assemble-vorticity-2d grads) node-vel vorticity)))
                  ))
            (setf velocity-rate (magicl:scale strain-rate 1d0))
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
                     (setf volume (* volume (det df)))

                     (multiple-value-bind (l v) (magicl:eig (magicl:@ def (magicl:transpose def)))
                       (let ((stretch
                               (magicl:@
                                v
                                (magicl:from-diag (mapcar (lambda (x) (the double-float (sqrt
                                                                                         (the double-float x)))) l) :type 'double-float)
                                (magicl:transpose v))))
                         (declare (type magicl:matrix/double-float stretch))
                         ;(setf (tref domain 0 0) (* (the double-float (tref domain 0 0))
                         ;                           (the double-float (tref stretch 0 0))))
                         ;(setf (tref domain 1 0) (* (the double-float (tref domain 1 0))
                         ;                           (the double-float (tref stretch 1 1))))
                         (setf (tref domain 0 0) (* (the double-float (tref domain-0 0 0))
                                                    (the double-float (tref stretch 0 0))))
                         (setf (tref domain 1 0) (* (the double-float (tref domain-0 1 0))
                                                    (the double-float (tref stretch 1 1))))
                         ))
                     ))))

(declaim (inline update-strain-kirchoff))
(declaim (ftype (function (cl-mpm/mesh::mesh
                           cl-mpm/particle:particle
                          magicl:matrix/double-float) (values))
               update-strain-kirchoff))
(defun update-strain-kirchoff (mesh mp dstrain)
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   (velocity-rate cl-mpm/particle::mp-velocity-rate)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0)
                   ) mp
    (declare (type double-float volume)
             (type magicl:matrix/double-float
                   domain
                   dstrain
                   ))
    (progn
      (let (;(df (calculate-df mesh mp))
            (df (.+ (magicl:eye 2) (voight-to-matrix dstrain)))
            (prev-strain strain))
        (progn
          (setf def (magicl:@ df def))
          ;;Eng strain to log strain
          ;;(voight-to-matrix )
          ;; (multiple-value-bind (l v) (magicl:eig (voight-to-matrix strain))
          (multiple-value-bind (l v) (magicl:eig (voight-to-matrix
                                                  (magicl:.* strain (magicl:from-list '(1d0 1d0 0.5d0) '(3 1)))))
            ;0.5
            (let ((trial-lgs (magicl:@ df
                                        v
                                        (magicl:from-diag (mapcar (lambda (x)
                                                                    (declare (type double-float x))
                                                                    (the double-float (exp (* x 2d0)))) l)
                                                          :type 'double-float)
                                        (magicl:transpose v)
                                       (magicl:transpose df))))
              (multiple-value-bind (lf vf) (magicl:eig trial-lgs)
                (setf strain (magicl:scale!
                               (matrix-to-voight
                                (magicl:@ vf
                                          (magicl:from-diag (mapcar (lambda (x)
                                                                      (log (the double-float x))) lf)
                                                            :type 'double-float)
                                          (magicl:transpose vf)))
                              0.5d0))
                (setf strain (magicl:.* strain (magicl:from-list '(1d0 1d0 2d0) '(3 1))))
                )))
        ;;Post multiply to turn to eng strain
        ;; (magicl:from-list '(1d0 1d0 2d0) '(3 1))

          ;(setf velocity-rate (magicl:scale strain-rate 1d0))
          ;; (setf strain-rate (magicl:.- strain prev-strain))
          ;(setf velocity-rate (magicl:scale strain-rate (/ 1d0 dt)))
          (setf volume (* volume (det df)))
          (when (<= volume 0d0)
            (error "Negative volume"))
          (multiple-value-bind (l v) (magicl:eig (magicl:@ def (magicl:transpose def)))
            (let ((stretch
                    (magicl:@
                     v
                     (magicl:from-diag (mapcar (lambda (x) (sqrt (the double-float x))) l) :type 'double-float)
                     (magicl:transpose v))))
              (declare (type magicl:matrix/double-float stretch))
              ;; (setf (tref domain 0 0) (* (the double-float (tref domain 0 0))
              ;;                            (the double-float (tref stretch 0 0))))
              ;; (setf (tref domain 1 0) (* (the double-float (tref domain 1 0))
              ;;                            (the double-float (tref stretch 1 1))))

              (setf (tref domain 0 0) (* (the double-float (tref domain-0 0 0))
                                         (the double-float (tref stretch 0 0))))
              (setf (tref domain 1 0) (* (the double-float (tref domain-0 1 0))
                                         (the double-float (tref stretch 1 1))))
              ))
          ))))
  (values))

(defun update-domain (mesh mp dt)
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   (velocity-rate cl-mpm/particle::mp-velocity-rate)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0)
                   ) mp
    (declare (type double-float volume)
             (type magicl:matrix/double-float
                   domain
                   ))
    (progn
      (let ((dstrain (magicl:zeros '(3 1))))
        (iterate-over-neighbours mesh mp
                                 (lambda (mesh mp node svp grads)
                                   (with-accessors ((node-vel cl-mpm/mesh:node-velocity))
                                       node
                                     (mult (cl-mpm/shape-function::assemble-dsvp-2d grads) node-vel dstrain)
                                     )
                                   ))
        (magicl:scale! dstrain dt)
        (let ((df (.+ (magicl:eye 2) (voight-to-matrix dstrain))))
            (multiple-value-bind (l v) (magicl:eig (magicl:@ df (magicl:transpose df)))
              (let ((stretch
                      (magicl:@
                       v
                       (magicl:from-diag (mapcar (lambda (x) (sqrt x)) l) :type 'double-float)
                       (magicl:transpose v))))
                (declare (type magicl:matrix/double-float stretch))
                (setf (tref domain 0 0) (* (the double-float (tref domain 0 0))
                                           (the double-float (tref stretch 0 0))))
                (setf (tref domain 1 0) (* (the double-float (tref domain 1 0))
                                           (the double-float (tref stretch 1 1))))
                ))))))
  (values))


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
                     ;; (damage cl-mpm/particle:mp-damage)
                        ) mp
      (progn
        ;;For normal
        (calculate-strain-rate mesh mp dt)
        (let ((dstrain (cl-mpm/particle:mp-strain-rate mp)))
          (progn
            (progn
              ;;Linear strain update
              ;; (update-strain-linear mesh mp dstrain)
              ;; (setf stress (cl-mpm/particle:constitutive-model mp strain dt))
              ;; (when (<= volume 0d0)
              ;;   (error "Negative volume"))

              ;; Turn cauchy stress to kirchoff
              (setf stress stress-kirchoff)
              ;Update our strains
              (update-strain-kirchoff mesh mp dstrain)
              ;;Update our kirchoff stress with constitutive model
              (setf stress-kirchoff (cl-mpm/particle:constitutive-model mp strain dt))
              ;;Check volume constraint!
              (when (<= volume 0d0)
                (error "Negative volume"))
              ;;Turn kirchoff stress to cauchy
              (setf stress (magicl:scale stress-kirchoff (/ 1.0d0 (magicl:det def))))
              ))))))
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
                   (volume cl-mpm/particle:mp-volume)
                   (def cl-mpm/particle:mp-deformation-gradient))
      mp
    (let* ((df (.+ (magicl:eye 2) (voight-to-matrix dstrain)))
           (j-inc (det df))
           (j-n (det def)))
      (iterate-over-neighbours mesh mp
                               (lambda (mesh mp node svp grads)
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
                   (def cl-mpm/particle:mp-deformation-gradient))
      mp
    (let* ((df (.+ (magicl:eye 2) (voight-to-matrix dstrain)))
           (j-inc (det df))
           (j-n (det def))
           (j-n1 0d0))
      (progn
        ;;If we want to do f bar
        (when nil
          (iterate-over-neighbours mesh mp
                                   (lambda (mesh mp node svp grads)
                                     (with-accessors ((node-active cl-mpm/mesh:node-active)
                                                      (node-j-inc cl-mpm/mesh::node-jacobian-inc)
                                                      (node-volume cl-mpm/mesh::node-volume)
                                                      )
                                         node
                                       (when node-active
                                         (incf j-n1 (/ (* svp node-j-inc) node-volume))
                                         ))))
          (when (<= (/ j-n1 (* j-n j-inc)) 0d0)
            (error "Negative volume"))
          (magicl:scale! df (expt (/ j-n1 (* j-n j-inc)) 1/2)))
        df))))

(defun update-stress (mesh mps dt)
  (declare ((array cl-mpm/particle:particle) mps) (cl-mpm/mesh::mesh mesh))
  ;; (with-accessors ((cells cl-mpm/mesh::mesh-cells))
  ;;     mesh
  ;;   (lparallel:pdotimes (i (array-total-size cells))
  ;;     (calculate-cell-deformation mesh (row-major-aref cells i) dt)
  ;;     )
  ;;   )
  ;;
  ;;Calculate jp calculate weighted value
  ;; (lparallel:pdotimes (i (length mps))
  ;;   (let ((mp (aref mps i)))
  ;;     (calculate-strain-rate mesh mp dt)
  ;;     (map-jacobian mesh mp dt)))

  ;;Update stress
  (lparallel:pdotimes (i (length mps))
     (update-stress-mp mesh (aref mps i) dt))
  (values))

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
                         (with-accessors ((damage cl-mpm/particle:mp-damage))
                             mp
                           (and (>= damage 1d0)
                                (split-criteria mp h)
                                ))) mps)))
      ;; (delete-if (lambda (mp)
      ;;              (with-accessors ((damage cl-mpm/particle:mp-damage))
      ;;                  mp
      ;;                (and (>= damage 1d0)
      ;;                     (split-criteria mp h)
      ;;                     ))) mps))
    ))
(defun split-criteria (mp h)
  (with-accessors ((def cl-mpm/particle:mp-deformation-gradient)
                   (lens cl-mpm/particle::mp-domain-size)
                   (lens-0 cl-mpm/particle::mp-domain-size-0)
                   )
      mp
    (let ((l-factor 1.10d0)
          (h-factor (* 1.5d0 h)))
      (cond
        ;; ((< h-factor (tref lens 0 0)) t)
        ;; ((< h-factor (tref lens 1 0)) t)
        ((< (* l-factor (tref lens-0 0 0)) (tref lens 0 0)) t)
        ((< (* l-factor (tref lens-0 1 0)) (tref lens 1 0)) t)
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
(defun split-mp (mp h)
  (with-accessors ((def cl-mpm/particle:mp-deformation-gradient)
                   (lens cl-mpm/particle::mp-domain-size)
                   (lens-0 cl-mpm/particle::mp-domain-size-0)
                   (pos cl-mpm/particle:mp-position)
                   (mass cl-mpm/particle:mp-mass)
                   (volume cl-mpm/particle:mp-volume))
      mp
    (let ((l-factor 1.25d0)
          (h-factor (* 1.0d0 h)))
    (cond
      ((< (* l-factor (tref lens-0 0 0)) (tref lens 0 0))
       ;; (< h-factor (tref lens 0 0))
       (let ((new-size (magicl:.* lens (magicl:from-list '(0.5d0 1d0) '(2 1))))
             (pos-offset (magicl:.* lens (magicl:from-list '(0.25d0 0d0) '(2 1)))))
         (list
          (copy-particle mp
                         :mass (/ mass 2)
                         :volume (/ volume 2)
                         :size new-size
                         :position (magicl:.+ pos pos-offset))
          (copy-particle mp
                         :mass (/ mass 2)
                         :volume (/ volume 2)
                         :size new-size
                         :position (magicl:.- pos pos-offset)))))
      ((< (* l-factor (tref lens-0 1 0)) (tref lens 1 0))
       (let ((new-size (magicl:.* lens (magicl:from-list '(1 0.5d0) '(2 1))))
             (pos-offset (magicl:.* lens (magicl:from-list '(0 0.25d0) '(2 1)))))
         (list
          (copy-particle mp
                         :mass (/ mass 2)
                         :volume (/ volume 2)
                         :size new-size
                         :position (magicl:.+ pos pos-offset))
          (copy-particle mp
                         :mass (/ mass 2)
                         :volume (/ volume 2)
                         :size new-size
                         :position (magicl:.- pos pos-offset)))))
      ;((< 2.0d0 (tref def 1 1)) t)
                                        ;((> 1.5d0 (tref def 0 0)) t)
      (t nil)
      ))))
(defun split-mps (sim)
  "Remove material points that have strain energy density exceeding fracture toughness"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (mps-to-split (remove-if-not (lambda (mp) (split-criteria mp h)) mps)))
      ;; (format t "~A ~%"h)
      ;(format t "~A~%" mps-to-split)
      ;; (print mps-to-split)
      ;(format t "~A~%" mps-to-split)
      (delete-if (lambda (mp) (split-criteria mp h)) mps)
      (loop for mp across mps-to-split
            do (loop for new-mp in (split-mp mp h)
                     do (vector-push-extend new-mp mps)))
      )
    ;; (delete-if (lambda (mp)
    ;;              (with-accessors ((damage cl-mpm/particle:mp-damage))
    ;;                  mp
    ;;                (>= damage 1d0))) mps))
  ))
#||
(progn
(ql:quickload :cl-mpm/examples/fracture)
(in-package :cl-mpm/examples/fracture))
||#

