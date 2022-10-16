(defpackage :cl-mpm
  (:use :cl 
        :cl-mpm/particle
        :cl-mpm/mesh
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
    ))
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
     :initform '())
   (damping-factor
     :type double-float
     :accessor sim-damping-factor
     :initarg :damping-factor
     :initform 1d0)
   (mass-filter
     :type double-float
     :accessor sim-mass-filter
     :initarg :mass-filter
     :initform 0.01d0))
  (:documentation "A self contained mpm simulation"))

(defun make-mpm-sim (size resolution dt shape-function)
  (make-instance 'mpm-sim
                 :dt (coerce dt 'double-float)
                 :mesh (make-mesh size resolution shape-function)
                 :mps '()))

(defun matrix-to-voight (matrix)
  (let* ( (exx (tref matrix 0 0))
          (eyy (tref matrix 1 1))
          (exy (tref matrix 1 0)))
    (magicl:from-list (list exx eyy
                            exy)
                      '(3 1) :type 'double-float)))
(defun voight-to-matrix (vec)
  (let* ( (exx (tref vec 0 0))
          (eyy (tref vec 1 0))
          (exy (tref vec 2 0)))
    (magicl:from-list (list exx exy 
                            exy eyy)
                      '(2 2) :type 'double-float)))

(defgeneric update-sim (sim)
  (:documentation "Update an mpm simulation by one timestep"))

(defmethod update-sim (sim)
  (with-slots ((mesh mesh)
               (mps mps)
               (bcs bcs)
               (dt dt))
                sim
                (progn
                    (reset-grid mesh)
                    (p2g mesh mps)
                    (when (> (sim-mass-filter sim) 0)
                      (filter-grid mesh (sim-mass-filter sim)))
                    (update-nodes mesh dt (sim-damping-factor sim))
                    (apply-bcs mesh bcs)
                    (g2p mesh mps)
                    (update-particle mps dt)
                    (update-stress mesh mps dt) 
                  )))

(defgeneric iterate-over-neighbours-shape (mesh shape-func mp func)
  (:documentation "For a given shape function iterate over relevant nodes and call func with mesh, mp, node, weight and gradients"))

(defun iterate-over-neighbours (mesh mp func)
  (iterate-over-neighbours-shape mesh (cl-mpm/mesh:mesh-shape-func mesh) mp func))

(defmethod iterate-over-neighbours-shape (mesh shape-func mp func)
  (let* ((nd (cl-mpm/mesh:mesh-nd mesh))
         (h (cl-mpm/mesh:mesh-resolution mesh))
         (order (slot-value shape-func 'order))
         (pos-vec (cl-mpm/particle:mp-position mp));Material point position 
         (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
         (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec 'floor)))
    (loop for dx from 0 to order
          do (loop for dy from 0 to order
                   do
                   (let* (
                          (id (mapcar #'+ pos-index (list dx dy)))
                          (dist (mapcar (lambda (p i) (- p (* i h))) pos id))
                          (node (cl-mpm/mesh:get-node mesh id))
                          (weight (apply (svp shape-func) dist))
                          (grads (apply (dsvp shape-func) dist)) 
                          )
                          (funcall func mesh mp node weight (assemble-dsvp nd grads)))))))

(defmethod iterate-over-neighbours-shape (mesh
                                          (shape-func shape-function-linear)
                                          (mp cl-mpm/particle::particle-damage)
                                          func)
  (let* ((nd (cl-mpm/mesh:mesh-nd mesh))
         (h (cl-mpm/mesh:mesh-resolution mesh))
         (order (slot-value shape-func 'order))
         (pos-vec (cl-mpm/particle:mp-position mp));Material point position 
         (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
         (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec 'floor)))
    (loop for dx from 0 to order
          do (loop for dy from 0 to order
                   do
                   (let* (
                          (id (mapcar #'+ pos-index (list dx dy)))
                          (dist (mapcar (lambda (p i) (- p (* i h))) pos id))
                          (node (cl-mpm/mesh:get-node mesh id))
                          (weight (apply (svp shape-func) dist))
                          (grads (apply (dsvp shape-func) dist))
                          )
                          (funcall func mesh mp node weight (assemble-dsvp nd grads)))))))

(defmethod iterate-over-neighbours-shape (mesh (shape-func shape-function-bspline) mp func)
  (let* ((nd (cl-mpm/mesh:mesh-nd mesh))
         (h (cl-mpm/mesh:mesh-resolution mesh))
         (order (slot-value shape-func 'order))
         (pos-vec (cl-mpm/particle:mp-position mp));Material point position 
         (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
         (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec)))
    (loop for dx from 0 to order
          do (loop for dy from 0 to order
                   do
                   (let* (
                          (id (mapcar #'+ pos-index (list dx dy)))
                          (dist (mapcar (lambda (p i) (- p (* i h))) pos id))
                          (weight (apply (svp shape-func) dist))
                          (grads (apply (dsvp shape-func) dist)) 
                          )
                     (when (cl-mpm/mesh:in-bounds mesh id)
                       (funcall func mesh mp (cl-mpm/mesh:get-node mesh id) weight (assemble-dsvp nd grads))))))))

;(defparameter *mesh* (sim-mesh *sim*))
;(defparameter *sf* (cl-mpm/mesh:mesh-shape-func *mesh*))
;(defparameter *mp* (cl-mpm/particle:make-particle 2 :pos '(0.5 0.1)))
;(iterate-over-neighbours-shape *mesh* *sf* *mp* 
;                               (lambda (mesh mp node w dsvp) (print (cl-mpm/mesh:node-index  node))))

(defun p2g-mp (mesh mp)
    (iterate-over-neighbours mesh mp
      (lambda (mesh mp node svp dsvp) 
      (progn
        (with-accessors (   (node-vel   cl-mpm/mesh:node-velocity)
                            (node-mass  cl-mpm/mesh:node-mass)
                            (node-force cl-mpm/mesh:node-force)
                            (node-lock  cl-mpm/mesh:node-lock)) node
          (with-accessors ( (mp-vel  cl-mpm/particle:mp-velocity)
                            (mp-mass cl-mpm/particle:mp-mass)
                            (strain-rate cl-mpm/particle:mp-strain-rate)) mp 
            (progn
              (setf strain-rate (magicl:scale strain-rate 0))
              (sb-thread:with-mutex (node-lock)
                (setf node-mass 
                      (+ node-mass (* mp-mass svp)))
                (setf node-vel 
                      (magicl:.+ node-vel (magicl:scale mp-vel (* mp-mass svp))))
                (setf node-force 
                      (magicl:.- node-force  (det-int-force mp node dsvp)))
                (setf node-force 
                      (magicl:.+ node-force  (det-ext-force mp node svp)))
                (setf strain-rate (.+ strain-rate (magicl:@ dsvp node-vel))))
              )))))))

(defun p2g (mesh mps)
  (declare (array mps) (cl-mpm/mesh::mesh mesh))
  (lparallel:pdotimes (i (length mps)) 
    (p2g-mp mesh (aref mps i))))

(defun g2p-mp-node (mesh mp node svp dsvp) 
  (declare (cl-mpm/mesh::mesh mesh))
  "Map grid to particle for one mp-node pair"
  (with-accessors ((node-vel cl-mpm/mesh:node-velocity)) node
    (with-accessors ((vel mp-velocity)
                     (mass mp-mass)) mp
      (progn 
        (setf vel (magicl:.+ vel (magicl:scale node-vel svp)))
        ;(setf strain-rate (magicl:.+ strain-rate (magicl:@ dsvp node-vel)))
        ))))

(defun g2p-mp (mesh mp)
  (declare (cl-mpm/mesh::mesh mesh))
  "Map one MP onto the grid"
  (with-accessors ((vel mp-velocity))
    mp
    (progn
      (setf vel (magicl:scale vel 0d0)))
    (iterate-over-neighbours mesh mp #'g2p-mp-node)))

(defun g2p (mesh mps)
  (declare (cl-mpm/mesh::mesh mesh) (array mps))
  "Map grid values to all particles"
  (lparallel:pdotimes (i (length mps)) 
    (g2p-mp mesh (aref mps i))))

(defun update-node (mesh dt node damping)
    (when (> (cl-mpm/mesh:node-mass node) 0)
      (with-accessors ( (mass  node-mass)
                        (vel   node-velocity)
                        (force node-force)
                        (acc   node-acceleration))
        node
        (progn 
          (setf vel (magicl:scale vel (/ 1.0 mass)))
          (setf acc (magicl:scale force (/ 1.0 mass)))
          (setf acc (magicl:.- acc (magicl:scale vel damping)))
          (setf vel (magicl:.+ vel (magicl:scale acc dt)))
          ))))

(defun update-nodes (mesh dt damping)
  "Explicitly integrate nodes"
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
  (lparallel:pdotimes (i (array-total-size nodes))
    (let ((node (row-major-aref nodes i)))
      (update-node mesh dt node damping)))))

(defun apply-bcs (mesh bcs)
  (declare (cl-mpm/mesh::mesh mesh))
  (with-accessors ( (nodes  mesh-nodes)
                    (nD     mesh-nD)
                    (mc     mesh-count)) mesh
    ;each bc is a list (pos value)
    (dolist (bc bcs)
      (let ((index (cl-mpm/bc:bc-index bc)))
        (when (cl-mpm/mesh:in-bounds mesh index);Potentially throw here
          (cl-mpm/bc:apply-bc bc (cl-mpm/mesh:get-node mesh index)))))))


(defun update-particle (mps dt)
  (declare (array mps))
  (loop for mp across mps
    do (with-accessors ((vel cl-mpm/particle:mp-velocity)
                        (pos cl-mpm/particle:mp-position)) mp
      (progn 
      (setf pos (magicl:.+ pos (magicl:scale vel dt)))))))

;Could include this in p2g but idk
(defun calculate-strain-rate (mesh mp dt)
    (let ((dstrain (magicl:zeros '(3 1))))
        (progn
            (iterate-over-neighbours mesh mp 
                (lambda (mesh mp node svp dsvp)
                    (setf dstrain (.+ dstrain (magicl:@ dsvp (cl-mpm/mesh:node-velocity node))))))
            (magicl:scale  dstrain dt))))

(defun update-damage (mp dt)
  (let ((critical-stress 1)
        (damage-increment 0))
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (damage cl-mpm/particle:mp-damage)) mp
        (progn
          (multiple-value-bind (u s v) (magicl:svd (voight-to-matrix stress))
            (loop for i from 0 to 1
                  do (when (< (tref s i i) critical-stress)
                       (incf damage-increment 1)))
            (incf damage damage-increment)
            (setf damage (min 1 (max 0 damage-increment)))
            (loop for i from 0 to 1
                  do (when (< (tref s i i) 0)
                       (setf (tref s i i) (* (tref s i i) (- 1 damage)))))
            (setf stress (matrix-to-voight (magicl:@ u s v)))
            ))
      )))
(defun rotation-matrix (degrees)
  (let ((angle (/ (* pi degrees) 180)))
    (magicl:from-list (list (cos angle) (- (sin angle))
                            (sin angle) (cos angle))
                      '(2 2)
                      :type 'double-float
                      )))
(defun update-stress-mp (mesh mp dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (volume cl-mpm/particle:mp-volume)
                     (strain cl-mpm/particle:mp-strain)
                     (def    cl-mpm/particle:mp-deformation-gradient)
                     (strain-rate cl-mpm/particle:mp-strain-rate)
                     (damage cl-mpm/particle:mp-damage)
                        ) mp
         
           (let ((dstrain (calculate-strain-rate mesh mp dt)
                          ;(magicl:scale strain-rate dt)
                          ))
             (progn
               (setf strain (magicl:.+ strain dstrain))
               (setf stress (cl-mpm/particle:constitutive-model mp strain))
                (let ((df (.+ (magicl:eye 2) (voight-to-matrix dstrain))))
                   (progn
                     (setf def (magicl:@ df def))
                     (setf volume (* volume (magicl:det df)))
                     (setf stress (magicl:scale stress (/ 1.0 (magicl:det def))))
                     (update-damage mp dt)
                     ))))))

(defun update-stress (mesh mps dt)
  (declare (array mps) (cl-mpm/mesh::mesh mesh))
  (lparallel:pdotimes (i (length mps))
     (update-stress-mp mesh (aref mps i) dt)))

(defun reset-grid (mesh)
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
  (dotimes (i (array-total-size nodes))
    (let ((node (row-major-aref nodes i)))
    (when (> (cl-mpm/mesh:node-mass node) 0)
      (cl-mpm/mesh:reset-node node))))))

(defun filter-grid (mesh mass-thresh)
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
  (dotimes (i (array-total-size nodes))
    (let ((node (row-major-aref nodes i)))
      (when (and (> (cl-mpm/mesh:node-mass node) 0) 
                 (< (cl-mpm/mesh:node-mass node) mass-thresh))
        (cl-mpm/mesh:reset-node node))))))

