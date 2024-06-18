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
(declaim (optimize (debug 0) (safety 0) (speed 3)))
;    #:make-shape-function
(in-package :cl-mpm/penalty)

(defun ssqrt (a)
  (* (signum a) (sqrt (abs a))))

(defun penetration-distance-point (point datum normal)
  "Get linear penetration distance"
  (let* ((ypos (cl-mpm/fastmath::dot point normal))
         (dist (- datum ypos)))
    (the double-float dist)))

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

(defparameter *debug-mutex* (sb-thread:make-mutex))
(defparameter *debug-force* 0d0)
(defparameter *debug-force-count* 0)
(defparameter *debug-force-mp-count* 0)
;;Only vertical condition
(defun apply-force-mps (mesh mps dt normal datum epsilon friction)
  "Update force on nodes, with virtual stress field from mps"
  ;;If we lose contact we need to zero out our friction force
  (with-accessors ((nd cl-mpm/mesh::mesh-nd))
      mesh
    (cl-mpm:iterate-over-mps
     mps
     (lambda (mp)
       (let* ((penetration-dist (penetration-distance mp datum normal)))
         (declare (double-float penetration-dist))
         (if (> penetration-dist 0d0)
             (progn
               ;;Contact
               (with-accessors ((volume cl-mpm/particle:mp-volume)
                                (pressure cl-mpm/particle::mp-pressure)
                                (mp-vel cl-mpm/particle::mp-velocity)
                                (mp-mass cl-mpm/particle::mp-mass)
                                (mp-contact cl-mpm/particle::mp-penalty-contact)
                                (mp-friction cl-mpm/particle::mp-penalty-frictional-force)
                                )
                   mp
                 (let* ((pen-point (penetration-point mp penetration-dist datum normal))
                        (normal-force (* (expt penetration-dist 1d0)
                                         epsilon
                                         ;; (expt volume (/ (- nd 1) nd))
                                         )))
                   (sb-thread:with-mutex (*debug-mutex*)
                     ;; (format t "Pen dist ~E - Normal Force ~E~%" penetration-dist normal-force)
                     (incf *debug-force* (* normal-force 1d0))
                     (incf *debug-force-count* 1))
                   (setf mp-contact t)
                   ;; (format t "Contact point ~A : dist ~F~%" (magicl::storage pen-point) penetration-dist)
                   ;;Iterate over neighbour nodes
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
                            (let* ((force (cl-mpm/utils:vector-zeros))
                                   (rel-vel (cl-mpm/fastmath:dot normal mp-vel))
                                   (tang-vel (cl-mpm/fastmath:fast-.- mp-vel (magicl:scale normal rel-vel)))
                                   (tang-vel-norm-squared (cl-mpm/fastmath::mag-squared tang-vel))
                                   (normal-damping 0d0)
                                   (svp (* svp (expt volume (/ (- nd 1) nd))))
                                   (damping-force (* normal-damping rel-vel))
                                   (force-friction mp-friction)
                                   (stick-friction (* friction (abs normal-force))))

                              ;; update trial frictional force
                              (when (> tang-vel-norm-squared 0d0)
                                ;; We have sliding behaviour
                                (cl-mpm/fastmath::fast-fmacc
                                 force-friction
                                 tang-vel
                                 (* -1d0 (/ epsilon 2d0) dt)))

                              (when (> (cl-mpm/fastmath::mag-squared force-friction) 0d0)
                                (if (> (cl-mpm/fastmath::mag force-friction) stick-friction)
                                    (progn
                                      (setf force-friction
                                            (magicl:scale
                                             (cl-mpm/fastmath:norm force-friction)
                                             stick-friction))
                                      (setf (cl-mpm/particle::mp-penalty-friction-stick mp) t))
                                    (progn
                                      (setf (cl-mpm/particle::mp-penalty-friction-stick mp) nil)))
                                (cl-mpm/fastmath::fast-.+ force force-friction force))
                              (setf mp-friction force-friction)

                              (cl-mpm/fastmath::fast-fmacc force
                                                           normal
                                                           normal-force)

                              (cl-mpm/fastmath::fast-fmacc node-force
                                                           force
                                                           svp)

                              (cl-mpm/fastmath::fast-fmacc node-ext-force
                                                           force
                                                           svp)
                              )))))))))
             (progn
               ;;No contact
               (setf (cl-mpm/particle::mp-penalty-contact mp) nil)
               (cl-mpm/fastmath::fast-zero (cl-mpm/particle::mp-penalty-frictional-force mp)))))))))

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
    (apply-force-mps mesh mps dt
                     normal
                     datum
                     epsilon
                     friction
                     )))

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
