(defpackage :cl-mpm/penalty
  (:use :cl
        :cl-mpm
        :cl-mpm/particle
        :cl-mpm/mesh
        :cl-mpm/utils
        :cl-mpm/fastmath
        )
  (:import-from
    :magicl tref .+ .-)
  (:export
   #:make-bc-penalty
   #:bc-penalty-friction))
;; (declaim (optimize (debug 0) (safety 0) (speed 3)))
(declaim (optimize (debug 3) (safety 3) (speed 0)))
;    #:make-shape-function
(in-package :cl-mpm/penalty)

(defun ssqrt (a)
  (* (signum a) (sqrt (abs a))))

(declaim
 (ftype (function
         (magicl::matrix/double-float
          double-float
          magicl::matrix/double-float)
         double-float)
        penetration-distance-point))
(defun penetration-distance-point (point datum normal)
  "Get linear penetration distance"
  (let* ((ypos (cl-mpm/fastmath::dot point normal))
         (dist (- datum ypos)))
    (the double-float dist)))

(declaim
 (ftype (function
         (cl-mpm/particle::particle
          double-float
          magicl::matrix/double-float)
         double-float)
        penetration-distance))
(defun penetration-distance (mp datum normal)
  "Get linear penetration distance"
  (let* ((ypos (cl-mpm/fastmath::dot (cl-mpm/particle:mp-position mp) normal))
         (yheight (cl-mpm/fastmath::dot
                   (magicl:scale
                    (cl-mpm/particle::mp-domain-size mp) 0.5d0)
                   (cl-mpm/fastmath::norm
                    (magicl:.* normal normal))))
         (dist (- datum (- ypos yheight))))
    (the double-float dist)))

(defun penetration-point (mp pen datum normal)
  "Get linear penetration distance"
  (let* ((pos (cl-mpm/particle:mp-position mp))
         (domain (cl-mpm/particle::mp-domain-size mp)))
    (magicl:.- pos
               (magicl:.* normal (magicl:scale domain 0.5d0))
               )))

(defun trial-corner (mp normal)
  "Get linear penetration distance"
  (let* ((corner (cl-mpm/fastmath:fast-.+
                  (cl-mpm/particle:mp-position mp)
                  (cl-mpm/fastmath:fast-.*
                   (cl-mpm/fastmath::fast-scale-vector (cl-mpm/particle::mp-domain-size mp) -0.5d0)
                   normal))))
    corner))

(defparameter *debug-mutex* (sb-thread:make-mutex))
(defparameter *debug-force* 0d0)
(defparameter *debug-force-count* 0)
(defparameter *debug-force-mp-count* 0)
;;Only vertical condition
(declaim (notinline apply-force-mps))
(defun apply-force-mps (mesh mps dt normal datum epsilon friction &optional (func-clip (lambda (mp) t))
                                                                    &key (damping 0.1d0))
  "Update force on nodes, with virtual stress field from mps"
  ;;If we lose contact we need to zero out our friction force
  (with-accessors ((nd cl-mpm/mesh::mesh-nd))
      mesh
    (cl-mpm:iterate-over-mps
     mps
     (lambda (mp)
       (let* ((penetration-dist (penetration-distance mp datum normal)))
         (declare (double-float penetration-dist))
         (when (> penetration-dist 0d0)
             (progn
               ;;Contact
               (with-accessors ((volume cl-mpm/particle:mp-volume)
                                (pressure cl-mpm/particle::mp-pressure)
                                (mp-vel cl-mpm/particle::mp-velocity)
                                (mp-mass cl-mpm/particle::mp-mass)
                                (mp-contact cl-mpm/particle::mp-penalty-contact)
                                (mp-friction cl-mpm/particle::mp-penalty-frictional-force)
                                (mp-normal-force cl-mpm/particle::mp-penalty-normal-force)
                                )
                   mp
                 (let* ((pen-point (penetration-point mp penetration-dist datum normal))
                        (normal-force (* (expt penetration-dist 1d0)
                                         epsilon
                                         ;; (expt volume (/ (- nd 1) nd))
                                         (expt volume (/ (- nd 1) nd))
                                         )))
                   (when (funcall func-clip
                                  (trial-corner mp normal)
                                  ;pen-point
                                  )
                     (sb-thread:with-mutex (*debug-mutex*)
                       (incf *debug-force* (* normal-force 1d0))
                       (incf *debug-force-count* 1))
                     (setf mp-contact t)
                     ;;Iterate over neighbour nodes
                     (let* ((force (cl-mpm/utils:vector-zeros))
                            (rel-vel (cl-mpm/fastmath:dot normal mp-vel))
                            (tang-vel (cl-mpm/fastmath:fast-.- mp-vel (magicl:scale normal rel-vel)))
                            (tang-vel-norm-squared (cl-mpm/fastmath::mag-squared tang-vel))
                            ;; (normal-damping (* damping (sqrt (/ epsilon (/ mp-mass volume)))))
                            (normal-damping (* (/ pi 2) damping (sqrt (* epsilon mp-mass))))
                            ;; (normal-damping (* (/ pi 2) damping (sqrt epsilon)))
                            (damping-force (* normal-damping rel-vel))
                            (force-friction mp-friction)
                            (stick-friction (* friction (abs normal-force))))
                       ;; update trial frictional force
                       (when (> tang-vel-norm-squared 0d0)
                         ;; We have sliding behaviour
                         (let* (;(tang-vel-norm (sqrt tang-vel-norm-squared))
                                ;; (tang-normal (cl-mpm/fastmath:norm tang-vel))
                                ;;Trial friction
                                        ;(trial-friction-force (* (/ epsilon 2d0) tang-vel-norm dt))
                                )
                           (cl-mpm/fastmath::fast-fmacc
                            force-friction
                            tang-vel
                            (* -1d0 (/ epsilon 2d0) dt))))
                       (incf mp-normal-force (- normal-force damping-force))
                       (when (> (cl-mpm/fastmath::mag-squared force-friction) 0d0)
                         (if (> (cl-mpm/fastmath::mag force-friction) stick-friction)
                             (progn
                               ;; (cl-mpm/fastmath::fast-scale! force-friction
                               ;;                              (/ (cl-mpm/fastmath::mag force-friction)
                               ;;                                 stick-friction))
                               (setf force-friction
                                     (magicl:scale
                                      (cl-mpm/fastmath:norm force-friction)
                                      stick-friction))
                               (setf (cl-mpm/particle::mp-penalty-friction-stick mp) t)
                               )
                             (progn
                               (setf (cl-mpm/particle::mp-penalty-friction-stick mp) nil)
                               ))
                         (cl-mpm/fastmath::fast-.+ force force-friction force))
                       (setf mp-friction force-friction)
                       (cl-mpm/fastmath::fast-fmacc force
                                                    normal
                                                    (- normal-force damping-force))
                       (cl-mpm::iterate-over-neighbours-point-linear-3d
                        mesh
                        pen-point
                        (lambda (mesh node svp grads)
                          (with-accessors ((node-force cl-mpm/mesh:node-force)
                                           (node-ext-force cl-mpm/mesh::node-external-force)
                                           (node-lock  cl-mpm/mesh:node-lock)
                                           (node-vel  cl-mpm/mesh:node-velocity)
                                           (node-mass  cl-mpm/mesh:node-mass)
                                           (node-boundary cl-mpm/mesh::node-boundary-node)
                                           (node-boundary-scalar cl-mpm/mesh::node-boundary-scalar)
                                           (node-active  cl-mpm/mesh:node-active))
                              node
                            (declare (double-float volume svp))
                            ;;Lock node for multithreading
                            (when node-active
                              ;;Lock node for multithreading
                              (sb-thread:with-mutex (node-lock)
                                ;; (let* ((svp (* svp (expt volume (/ (- nd 1) nd))))))

                                (cl-mpm/fastmath::fast-fmacc node-force
                                                             force
                                                             svp)
                                (cl-mpm/fastmath::fast-fmacc node-ext-force
                                                             force
                                                             svp)
                                ))))))))))))))))


(defun iterate-over-corners-2d (mesh mp func)
  (declare (cl-mpm/particle::particle mp)
           (function func))
  (array-operations/utilities:nested-loop (x y) '(2 2)
    (let ((domain (cl-mpm/particle::mp-domain-size mp))
          (position (cl-mpm/particle:mp-position mp))
          (corner (cl-mpm/utils:vector-zeros)))
      (cl-mpm/fastmath::fast-.+-vector
       position
       (magicl:scale!
        (magicl:.*
         (vector-from-list (mapcar (lambda (x) (- (* 2d0 (coerce x 'double-float)) 1d0)) (list x y 0d0)))
         domain
         ) 0.5d0) corner)
      (funcall func corner))))

(defun iterate-over-corners-3d (mesh mp func)
  (declare (cl-mpm/particle::particle mp)
           (function func))
  (array-operations/utilities:nested-loop (x y z) '(2 2 2)
    (let ((domain (cl-mpm/particle::mp-domain-size mp))
          (position (cl-mpm/particle:mp-position mp))
          (corner (cl-mpm/utils:vector-zeros)))
      (cl-mpm/fastmath::fast-.+-vector
       position
       (magicl:scale!
        (magicl:.*
         (vector-from-list (mapcar (lambda (x) (- (* 2d0 (coerce x 'double-float)) 1d0)) (list x y z)))
         domain
         ) 0.5d0) corner)
      (funcall func corner))))

(defun iterate-over-corners (mesh mp func)
  (if (= (the fixnum (cl-mpm/mesh:mesh-nd mesh)) 2)
      (iterate-over-corners-2d mesh mp func)
      (iterate-over-corners-3d mesh mp func)))

(defun apply-penalty-force-gimp-mps (mesh mps dt normal datum epsilon friction &optional (func-clip (lambda (mp) t))
                                                                    &key (damping 0.1d0))
  "Update force on nodes, with virtual stress field from mps"
  ;;If we lose contact we need to zero out our friction force
  (with-accessors ((nd cl-mpm/mesh::mesh-nd))
      mesh
    (cl-mpm:iterate-over-mps
     mps
     (lambda (mp)
       (iterate-over-corners
        mesh
        mp
        (lambda (corner)
          (let* ((penetration-dist (penetration-distance-point corner datum normal)))
            (declare (double-float penetration-dist))
            (when (> penetration-dist 0d0)
              (progn
                ;;Contact
                (with-accessors ((volume cl-mpm/particle:mp-volume)
                                 (pressure cl-mpm/particle::mp-pressure)
                                 (mp-vel cl-mpm/particle::mp-velocity)
                                 (mp-mass cl-mpm/particle::mp-mass)
                                 (mp-contact cl-mpm/particle::mp-penalty-contact)
                                 (mp-friction cl-mpm/particle::mp-penalty-frictional-force)
                                 (mp-normal-force cl-mpm/particle::mp-penalty-normal-force)
                                 )
                    mp
                  (let* ((pen-point corner)
                         (normal-force (* (expt penetration-dist 1d0)
                                          epsilon
                                          ;; (expt volume (/ (- nd 1) nd))
                                          (expt volume (/ (- nd 1) nd))
                                          )))
                    (when (and (funcall func-clip corner))
                      (sb-thread:with-mutex (*debug-mutex*)
                        (incf *debug-force* (* normal-force 1d0))
                        (incf *debug-force-count* 1))
                      (setf mp-contact t)
                      ;;Iterate over neighbour nodes
                      (let* ((force (cl-mpm/utils:vector-zeros))
                             (rel-vel (cl-mpm/fastmath:dot normal mp-vel))
                             (tang-vel (cl-mpm/fastmath:fast-.- mp-vel (magicl:scale normal rel-vel)))
                             (tang-vel-norm-squared (cl-mpm/fastmath::mag-squared tang-vel))
                             ;; (normal-damping (* damping (sqrt (/ epsilon (/ mp-mass volume)))))
                             (normal-damping (* (/ pi 2) damping (sqrt (* epsilon mp-mass))))
                             ;; (normal-damping (* (/ pi 2) damping (sqrt epsilon)))
                             (damping-force (* normal-damping rel-vel))
                             (force-friction mp-friction)
                             (stick-friction (* friction (abs normal-force))))
                        ;; update trial frictional force
                        (when (> tang-vel-norm-squared 0d0)
                          ;; We have sliding behaviour
                          (let* (;(tang-vel-norm (sqrt tang-vel-norm-squared))
                                 ;; (tang-normal (cl-mpm/fastmath:norm tang-vel))
                                 ;;Trial friction
                                        ;(trial-friction-force (* (/ epsilon 2d0) tang-vel-norm dt))
                                 )
                            (cl-mpm/fastmath::fast-fmacc
                             force-friction
                             tang-vel
                             (* -1d0 (/ epsilon 2d0) dt))))
                        (incf mp-normal-force (- normal-force damping-force))
                        (when (> (cl-mpm/fastmath::mag-squared force-friction) 0d0)
                          (if (> (cl-mpm/fastmath::mag force-friction) stick-friction)
                              (progn
                                ;; (cl-mpm/fastmath::fast-scale! force-friction
                                ;;                              (/ (cl-mpm/fastmath::mag force-friction)
                                ;;                                 stick-friction))
                                (setf force-friction
                                      (magicl:scale
                                       (cl-mpm/fastmath:norm force-friction)
                                       stick-friction))
                                (setf (cl-mpm/particle::mp-penalty-friction-stick mp) t)
                                )
                              (progn
                                (setf (cl-mpm/particle::mp-penalty-friction-stick mp) nil)
                                ))
                          (cl-mpm/fastmath::fast-.+ force force-friction force))
                        (setf mp-friction force-friction)
                        (cl-mpm/fastmath::fast-fmacc force
                                                     normal
                                                     (- normal-force damping-force))
                        (cl-mpm::iterate-over-neighbours-point-linear-3d
                         mesh
                         pen-point
                         (lambda (mesh node svp grads)
                           (with-accessors ((node-force cl-mpm/mesh:node-force)
                                            (node-ext-force cl-mpm/mesh::node-external-force)
                                            (node-lock  cl-mpm/mesh:node-lock)
                                            (node-vel  cl-mpm/mesh:node-velocity)
                                            (node-mass  cl-mpm/mesh:node-mass)
                                            (node-boundary cl-mpm/mesh::node-boundary-node)
                                            (node-boundary-scalar cl-mpm/mesh::node-boundary-scalar)
                                            (node-active  cl-mpm/mesh:node-active))
                               node
                             (declare (double-float volume svp))
                             ;;Lock node for multithreading
                             (when node-active
                               ;;Lock node for multithreading
                               (sb-thread:with-mutex (node-lock)
                                 (cl-mpm/fastmath::fast-fmacc node-force
                                                              force
                                                              svp)
                                 (cl-mpm/fastmath::fast-fmacc node-ext-force
                                                              force
                                                              svp)))))))))))))))))))

(defun collect-contact-points (mesh mps normal datum)
  (loop for mp across mps
        when (> (penetration-distance mp datum normal) 0d0)
          collect
          (let* ((penetration-dist (penetration-distance mp datum normal)))
            (declare (double-float penetration-dist))
            (when (> penetration-dist 0d0)
              (with-accessors ((volume cl-mpm/particle:mp-volume)
                               (pressure cl-mpm/particle::mp-pressure)
                               (mp-vel cl-mpm/particle::mp-velocity)
                               )
             mp
           (penetration-point mp penetration-dist datum normal))))))
(defun collect-contact-points-bc (mesh mps bc)
  (with-accessors ((normal bc-penalty-normal)
                   (datum bc-penalty-datum)
                   )
      bc
    (collect-contact-points mesh mps normal datum)))

(defun apply-penalty (sim normal datum epsilon friction)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm::sim-mps)
                   (dt cl-mpm::sim-dt)
                   )
      sim
    (apply-penalty-force-gimp-mps
     mesh
     mps
     dt
     normal
     datum
     epsilon
     friction)))

(defclass bc-penalty (cl-mpm/bc::bc)
  ((normal
    :accessor bc-penalty-normal
    :initarg :normal
    :initform (cl-mpm/utils:vector-zeros))
   (datum
    :accessor bc-penalty-datum
    :initarg :datum
    :initform 0d0)
   (sim
    :accessor bc-penalty-sim
    :initarg :sim)
   (friction
    :accessor bc-penalty-friction
    :initarg :friction)
   (epsilon
    :accessor bc-penalty-epsilon
    :initarg :epsilon)
   (load
    :accessor bc-penalty-load
    :initform 0d0)
   (load-lock
    :accessor bc-penalty-load-lock
    :initform (sb-thread:make-mutex)))
  (:documentation "A nonconforming neumann bc"))
(defclass bc-penalty-distance (bc-penalty)
  ((center-point
    :accessor bc-penalty-distance-center-point
    :initarg :center-point
    :initform (cl-mpm/utils:vector-zeros))
   (radius
    :accessor bc-penalty-distance-radius
    :initarg :radius
    :initform 0d0)
   ))

(defun make-bc-penalty-distance-point (sim normal point radius epsilon friction damping)
  (let* ((normal (cl-mpm/fastmath::norm normal))
         (datum (- (penetration-distance-point point 0d0 normal))))
    (make-instance 'bc-penalty-distance
                   :index nil
                   :sim sim
                   :datum datum
                   :normal normal
                   :epsilon epsilon
                   :friction friction
                   :center-point point
                   :radius radius)))


(defgeneric penalty-contact-valid (bc point))

(defmethod penalty-contact-valid ((bc bc-penalty) point)
  t)
(defmethod penalty-contact-valid ((bc bc-penalty-distance) point)
  (let* ((normal (bc-penalty-normal bc))
         (diff
          (cl-mpm/fastmath::fast-.-
           point
           (bc-penalty-distance-center-point bc))))
    (<
     (cl-mpm/fastmath::mag
      (cl-mpm/fastmath:fast-.-
       diff
       (cl-mpm/fastmath:fast-.*
        diff
        normal)))
     (bc-penalty-distance-radius bc))))


(defmethod cl-mpm/bc::apply-bc ((bc bc-penalty) node mesh dt)
  (with-accessors ((datum bc-penalty-datum)
                   (rho bc-buoyancy-rho)
                   (normal bc-penalty-normal)
                   (epsilon bc-penalty-epsilon)
                   (friction bc-penalty-friction)
                   (sim bc-penalty-sim))
      bc
    (apply-penalty
     sim
     normal
     datum
     epsilon
     friction)))

(defun make-bc-penalty (sim datum epsilon friction)
  (let ((normal (cl-mpm/fastmath::norm (cl-mpm/utils:vector-zeros))))
    (setf (magicl:tref normal (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)) 0) 1d0)
    (setf normal (cl-mpm/fastmath::norm normal))
    (make-instance 'bc-penalty
                   :index nil
                   :sim sim
                   :datum datum
                   :normal normal
                   :epsilon epsilon
                   :friction friction)))

(defun make-bc-penalty-point-normal (sim normal point epsilon friction)
  (let* ((normal (cl-mpm/fastmath::norm normal))
         ;(datum (cl-mpm/fastmath::dot normal point))
         (datum (- (penetration-distance-point point 0d0 normal)))
         )
    ;; (format t "Normal ~F ~F ~%" (magicl:tref normal 0 0) (magicl:tref normal 1 0))
    (make-instance 'bc-penalty
                   :index nil
                   :sim sim
                   :datum datum
                   :normal normal
                   :epsilon epsilon
                   :friction friction
                   )))

(defun disp-distance (mp datum normal)
  "Get linear penetration distance"
  (let* ((ypos (cl-mpm/fastmath::dot (cl-mpm/particle:mp-position mp) normal))
         (yheight (cl-mpm/fastmath::dot (magicl:scale (cl-mpm/particle::mp-domain-size mp) 0.5d0)
                                        (cl-mpm/fastmath::norm (magicl:.* normal normal))))
         (dist (- datum (- ypos yheight))))
    (the double-float dist)))

(defun disp-penetration-point (mp pen datum normal)
  "Get linear penetration distance"
  (let* ((pos (cl-mpm/particle:mp-position mp))
         (domain (cl-mpm/particle::mp-domain-size mp)))
    (magicl:.- pos
               (magicl:.* normal (magicl:scale domain 0.5d0)))))

(defun apply-displacement-control-mps (mesh mps dt normal datum epsilon friction)
  "Update force on nodes, with virtual stress field from mps"
  (declare (double-float datum dt epsilon friction))
  (cl-mpm:iterate-over-mps
   mps
   (lambda (mp)
     (let* ((penetration-dist (disp-distance mp datum normal)))
       (declare (double-float penetration-dist))
       (when t;(> (abs penetration-dist) 0d0)
         (with-accessors ((volume cl-mpm/particle:mp-volume)
                          (pressure cl-mpm/particle::mp-pressure)
                          (mp-vel cl-mpm/particle::mp-velocity)
                          (mp-mass cl-mpm/particle::mp-mass)
                          )
             mp
           (let* ((pen-point (disp-penetration-point mp penetration-dist datum normal))
                  (h (cl-mpm/mesh:mesh-resolution mesh))
                  (nd (cl-mpm/mesh:mesh-nd mesh))
                  (normal-force (* (signum penetration-dist)
                                   (expt (abs penetration-dist) 1d0)
                                   epsilon
                                   ;; h
                                   (expt volume (/ (- nd 1) nd))
                                   ;(expt h (/ 1 nd))
                                   )))
             (sb-thread:with-mutex (*debug-mutex*)
               (incf *debug-force* (* normal-force 1d0))
               (incf *debug-force-count* 1))
             ;;Iterate over neighbour nodes
             (cl-mpm::iterate-over-neighbours-point-linear
              mesh pen-point
              (lambda (mesh node svp grads)
                (when (cl-mpm/mesh:node-active node)
                  (with-accessors ((node-force cl-mpm/mesh:node-force)
                                   (node-ext-force cl-mpm/mesh::node-external-force)
                                   (node-lock  cl-mpm/mesh:node-lock)
                                   (node-vel  cl-mpm/mesh:node-velocity)
                                   (node-mass  cl-mpm/mesh:node-mass)
                                   (node-boundary cl-mpm/mesh::node-boundary-node)
                                   (node-boundary-scalar cl-mpm/mesh::node-boundary-scalar))
                      node
                    (declare (double-float volume svp))
                    ;;Lock node for multithreading
                    (let* ((force (cl-mpm/utils:vector-zeros))
                           (rel-vel (cl-mpm/fastmath::dot normal mp-vel))
                           ;; (svp (* svp 1))
                           ;; (normal-force (* dt normal-force))
                           (normal-damping 0d0)
                           (damping-force 0d0
                                          ;(* normal-damping rel-vel)
                                          ))
                      (declare (double-float rel-vel normal-damping damping-force))

                      (cl-mpm/fastmath::fast-fmacc force
                                                   normal
                                                   (- normal-force
                                                      damping-force))

                      (sb-thread:with-mutex (node-lock)
                        (cl-mpm/fastmath::fast-fmacc
                         node-force
                         force
                         svp)
                        (cl-mpm/fastmath::fast-fmacc
                         node-ext-force
                         force
                         svp))))))))))))))


(defstruct line-segment
  normal
  datum
  point-start
  point-end)

(defclass bc-penalty-structure (cl-mpm/bc::bc)
  ((sim
    :accessor bc-penalty-sim
    :initarg :sim)
   (friction
    :accessor bc-penalty-friction
    :initarg :friction)
   (epsilon
    :accessor bc-penalty-epsilon
    :initarg :epsilon)
   (damping
    :accessor bc-penalty-damping
    :initarg :damping)
   (load
    :accessor bc-penalty-load
    :initform 0d0)
   (load-lock
    :accessor bc-penalty-load-lock
    :initform (sb-thread:make-mutex))
   (sub-bcs
    :accessor bc-penalty-structure-sub-bcs
    :initarg :sub-bcs
    :initform (list))
   (contact-points
    :accessor bc-penalty-structure-contact-points
    :initform (list)))
  (:documentation "A single multi-surface structure that should resolve contact through a closest-point algorithm"))

(defun make-bc-penalty-structure (sim epsilon friction damping sub-bcs)
  (make-instance 'bc-penalty-structure
                  :index nil
                  :sim sim
                  :epsilon epsilon
                  :friction friction
                  :damping damping
                  :sub-bcs sub-bcs))

(defstruct contact
  point
  datum
  penetration
  normal
  sub-bc
  )

(defmethod cl-mpm/bc::apply-bc ((bc bc-penalty-structure) node mesh dt)
  (with-accessors ((epsilon bc-penalty-epsilon)
                   (friction bc-penalty-friction)
                   (damping bc-penalty-damping)
                   (sub-bcs bc-penalty-structure-sub-bcs)
                   (debug-mutex bc-penalty-load-lock)
                   (debug-force bc-penalty-load)
                   (sim bc-penalty-sim))
      bc
    (setf (bc-penalty-structure-contact-points bc) nil)
    (with-accessors ((mps cl-mpm:sim-mps)
                     (mesh cl-mpm:sim-mesh))
        sim
      (cl-mpm:iterate-over-mps
       mps
       (lambda (mp)
         (cl-mpm/penalty::iterate-over-corners
          mesh
          mp
          (lambda (corner)
            (let ((in-contact nil)
                  (closest-point (make-contact :penetration 0d0)))
              (loop for sub-bc in sub-bcs
                    do (with-accessors ((datum bc-penalty-datum)
                                        (normal bc-penalty-normal))
                           sub-bc
                         (let* ((penetration-dist (penetration-distance-point corner datum normal)))
                           (declare (double-float penetration-dist))
                           (when (and
                                  (> penetration-dist 0d0)
                                  (penalty-contact-valid sub-bc corner))
                             (if in-contact
                                 (when (< (abs penetration-dist) (abs (contact-penetration closest-point)))
                                   (setf closest-point (make-contact
                                                        :point corner
                                                        :datum (bc-penalty-datum sub-bc)
                                                        :penetration penetration-dist
                                                        :normal (bc-penalty-normal sub-bc)
                                                        :sub-bc sub-bc)))
                                 (progn
                                   (setf in-contact t)
                                   (setf closest-point (make-contact
                                                        :point corner
                                                        :datum (bc-penalty-datum sub-bc)
                                                        :penetration penetration-dist
                                                        :normal (bc-penalty-normal sub-bc)
                                                        :sub-bc sub-bc))))))))
              (when in-contact
                (with-accessors ((nd cl-mpm/mesh:mesh-nd))
                    mesh
                  (with-accessors ((datum bc-penalty-datum)
                                   (normal bc-penalty-normal))
                      (contact-sub-bc closest-point)
                    (with-accessors ((volume cl-mpm/particle:mp-volume)
                                     (pressure cl-mpm/particle::mp-pressure)
                                     (mp-vel cl-mpm/particle::mp-velocity)
                                     (mp-mass cl-mpm/particle::mp-mass)
                                     (mp-contact cl-mpm/particle::mp-penalty-contact)
                                     (mp-friction cl-mpm/particle::mp-penalty-frictional-force)
                                     (mp-normal-force cl-mpm/particle::mp-penalty-normal-force)
                                     )
                        mp
                      (push (contact-point closest-point) (bc-penalty-structure-contact-points bc))
                      (let* ((pen-point (contact-point closest-point))
                             (normal-force (* (expt (contact-penetration closest-point) 1d0)
                                              epsilon
                                              (expt volume (/ (- nd 1) nd)))))
                        (sb-thread:with-mutex (debug-mutex)
                          (incf debug-force (* normal-force 1d0)))
                        (let ((sub-bc (contact-sub-bc closest-point)))
                          (sb-thread:with-mutex ((bc-penalty-load-lock sub-bc))
                            (incf (bc-penalty-load sub-bc) normal-force)))
                        (setf mp-contact t)
                        (let* ((force (cl-mpm/utils:vector-zeros))
                               (rel-vel (cl-mpm/fastmath:dot normal mp-vel))
                               (tang-vel (cl-mpm/fastmath:fast-.- mp-vel (magicl:scale normal rel-vel)))
                               (tang-vel-norm-squared (cl-mpm/fastmath::mag-squared tang-vel))
                               (normal-damping (* (/ pi 2) damping (sqrt (* epsilon mp-mass))))
                               (damping-force (* normal-damping rel-vel))
                               (force-friction mp-friction)
                               (stick-friction (* friction (abs normal-force))))
                          ;; update trial frictional force
                          (when (> tang-vel-norm-squared 0d0)
                            (cl-mpm/fastmath::fast-fmacc
                             force-friction
                             tang-vel
                             (* -1d0 (/ epsilon 2d0) dt)))
                          (incf mp-normal-force (- normal-force damping-force))
                          (when (> (cl-mpm/fastmath::mag-squared force-friction) 0d0)
                            (if (> (cl-mpm/fastmath::mag force-friction) stick-friction)
                                (progn
                                  (setf force-friction
                                        (magicl:scale
                                         (cl-mpm/fastmath:norm force-friction)
                                         stick-friction))
                                  (setf (cl-mpm/particle::mp-penalty-friction-stick mp) t)
                                  )
                                (progn
                                  (setf (cl-mpm/particle::mp-penalty-friction-stick mp) nil)
                                  ))
                            (cl-mpm/fastmath::fast-.+ force force-friction force))
                          (setf mp-friction force-friction)
                          (cl-mpm/fastmath::fast-fmacc force
                                                       normal
                                                       (- normal-force damping-force))
                          (cl-mpm::iterate-over-neighbours-point-linear-3d
                           mesh
                           pen-point
                           (lambda (mesh node svp grads)
                             (with-accessors ((node-force cl-mpm/mesh:node-force)
                                              (node-ext-force cl-mpm/mesh::node-external-force)
                                              (node-lock  cl-mpm/mesh:node-lock)
                                              (node-vel  cl-mpm/mesh:node-velocity)
                                              (node-mass  cl-mpm/mesh:node-mass)
                                              (node-boundary cl-mpm/mesh::node-boundary-node)
                                              (node-boundary-scalar cl-mpm/mesh::node-boundary-scalar)
                                              (node-active  cl-mpm/mesh:node-active))
                                 node
                               (declare (double-float volume svp))
                               ;;Lock node for multithreading
                               (when node-active
                                 ;;Lock node for multithreading
                                 (sb-thread:with-mutex (node-lock)
                                   (cl-mpm/fastmath::fast-fmacc node-force
                                                                force
                                                                svp)
                                   (cl-mpm/fastmath::fast-fmacc node-ext-force
                                                                force
                                                                svp))))))))))))))))))))



