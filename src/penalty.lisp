
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
    ))
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
         (yheight (cl-mpm/fastmath::dot (magicl:scale (cl-mpm/particle::mp-domain-size mp) 0.5d0) (cl-mpm/fastmath::norm (magicl:.* normal normal))))
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
  (lparallel:pdotimes
      (i (length mps))
    (let ((mp (aref mps i)))

      (let* ((penetration-dist (penetration-distance mp datum normal)))
        (declare (double-float penetration-dist))
        (when (> penetration-dist 0d0)
          ;; (break)
          (with-accessors ((volume cl-mpm/particle:mp-volume)
                           (pressure cl-mpm/particle::mp-pressure)
                           (mp-vel cl-mpm/particle::mp-velocity)
                           (mp-mass cl-mpm/particle::mp-mass)
                           )
              mp
            (let* ((pen-point (penetration-point mp penetration-dist datum normal))
                   (normal-force (* (expt penetration-dist 1d0) epsilon)))
              (sb-thread:with-mutex (*debug-mutex*)
                (incf *debug-force* (* normal-force 1d0))
                (incf *debug-force-count* 1))
              ;; (format t "Contact point ~A : dist ~F~%" (magicl::storage pen-point) penetration-dist)
              ;;Iterate over neighbour nodes
              (cl-mpm::iterate-over-neighbours-point-linear-3d
               mesh pen-point
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
                              (rel-vel (magicl::sum (magicl::.* normal mp-vel)))
                              (tang-vel (magicl:.- mp-vel (magicl:scale normal rel-vel)))
                              (normal-damping 0d0)
                              ;; (svp (* svp volume))
                              (damping-force (* normal-damping rel-vel)))

                         ;; (when (> (cl-mpm/fastmath:dot tang-vel tang-vel) 0d0)
                         ;;   (let ((tang-normal (cl-mpm/fastmath:norm tang-vel))
                         ;;         (force-friction (cl-mpm/utils:vector-zeros)))
                         ;;     (cl-mpm/fastmath::fast-fmacc force-friction
                         ;;                                  tang-normal
                         ;;                                  (* -1d0
                         ;;                                     friction
                         ;;                                     (coerce (abs normal-force) 'double-float)
                         ;;                                     ))
                         ;;     (magicl:.+ force force-friction force))
                         ;; )

                         (cl-mpm/fastmath::fast-fmacc force
                                                      normal
                                                      (- normal-force
                                                         damping-force))
                         (cl-mpm/fastmath::fast-fmacc node-force
                                                      force
                                                      svp)
                         (cl-mpm/fastmath::fast-fmacc node-ext-force
                                                      force
                                                      svp)
                         )))
                   ;; (sb-thread:with-mutex (node-lock)
                   ;;   (let* ((force (cl-mpm/utils:vector-zeros))
                   ;;          (rel-vel (magicl::sum (magicl::.* normal mp-vel)))
                   ;;          (tang-vel (magicl:.- mp-vel (magicl:scale normal rel-vel)))
                   ;;          (normal-damping 1d6)
                   ;;          (damping-force (* normal-damping rel-vel)))


                   ;;     (cl-mpm/fastmath::fast-fmacc force
                   ;;                                  normal
                   ;;                                  (- normal-force
                   ;;                                     damping-force))
                   ;;     (when (> (cl-mpm/fastmath:dot tang-vel tang-vel) 0d0)
                   ;;       (let ((tang-normal (cl-mpm/fastmath:norm tang-vel))
                   ;;             (force-friction (cl-mpm/utils:vector-zeros)))
                   ;;         (cl-mpm/fastmath::fast-fmacc force-friction
                   ;;                                      tang-normal
                   ;;                                      (* -1d0
                   ;;                                         friction
                   ;;                                         (coerce (abs normal-force) 'double-float)
                   ;;                                         ))
                   ;;         ;; (break)
                   ;;         (magicl:.+ force force-friction force))
                   ;;     )
                   ;;     (cl-mpm/fastmath::fast-fmacc node-force
                   ;;                                  force
                   ;;                                  svp
                   ;;                                  ;; (* svp volume)
                   ;;                                  )
                   ;;     ))
                   ))))))))))
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

;; (declaim (notinline apply-non-conforming-nuemann))
;; (defun apply-non-conforming-nuemann (sim func-stress func-div)
;;   (with-accessors ((mesh cl-mpm:sim-mesh)
;;                    (mps cl-mpm::sim-mps))
;;       sim
;;     (with-accessors ((h cl-mpm/mesh:mesh-resolution))
;;         ;; mesh
;;       ;; (locate-mps-cells mesh mps)
;;       (populate-nodes-volume mesh)
;;       ;; (populate-nodes-domain mesh)
;;       (apply-force-mps mesh mps
;;                        (lambda (mp) (calculate-val-mp mp func-stress))
;;                        (lambda (mp) (calculate-val-mp mp func-div)))
;;       (apply-force-cells mesh
;;                          func-stress
;;                          func-div
;;                          )
;;       )))
(defun apply-penalty (sim normal datum epsilon friction)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm::sim-mps)
                   (dt cl-mpm::sim-dt)
                   )
      sim
    (apply-force-mps mesh mps dt normal datum epsilon friction)))

(defclass bc-penalty (cl-mpm/bc::bc)
  ((normal
    :accessor bc-penalty-normal
    :initarg :normal
    :initform (magicl:zeros '(3 1))
    )
   (datum
    :accessor bc-penalty-datum
    :initarg :datum
    :initform 0d0
    )
   (sim
    :accessor bc-penalty-sim
    :initarg :sim)
   (friction
    :accessor bc-penalty-friction
    :initarg :friction)
   (epsilon
    :accessor bc-penalty-epsilon
    :initarg :epsilon)

   )
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
     friction)
    ))

(defun make-bc-penalty (sim datum epsilon friction)
  (let ((normal (cl-mpm/fastmath::norm (cl-mpm/utils:vector-zeros))))
    (setf (magicl:tref normal (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)) 0) 1d0)
    (setf normal (cl-mpm/fastmath::norm normal))
    (make-instance 'bc-penalty
                   :index '(0 0 0)
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
                   :index '(0 0 0)
                   :sim sim
                   :datum datum
                   :normal normal
                   :epsilon epsilon
                   :friction friction
                   )))

(defun disp-distance (mp datum normal)
  "Get linear penetration distance"
  (let* ((ypos (cl-mpm/fastmath::dot (cl-mpm/particle:mp-position mp) normal))
         (yheight (cl-mpm/fastmath::dot (magicl:scale (cl-mpm/particle::mp-domain-size mp) 0.5d0) (cl-mpm/fastmath::norm (magicl:.* normal normal))))
         (dist (- datum (- ypos yheight))))
    (the double-float dist)))

(defun disp-penetration-point (mp pen datum normal)
  "Get linear penetration distance"
  (let* ((pos (cl-mpm/particle:mp-position mp))
         (domain (cl-mpm/particle::mp-domain-size mp)))
    (magicl:.- pos
               (magicl:.* normal (magicl:scale domain 0.5d0))
               )))

(defun apply-displacement-control-mps (mesh mps dt normal datum epsilon friction)
  "Update force on nodes, with virtual stress field from mps"
  (declare (double-float datum dt epsilon friction))
  (lparallel:pdotimes
      (i (length mps))
    (let ((mp (aref mps i)))
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
                   (normal-force (* (expt penetration-dist 1d0) epsilon
                                    ;(sqrt volume)
                                    )))
              (sb-thread:with-mutex (*debug-mutex*)
                (incf *debug-force* (* normal-force 1d0))
                (incf *debug-force-count* 1))
              ;;Iterate over neighbour nodes
              (cl-mpm::iterate-over-neighbours-point-linear-3d
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
                            (rel-vel (magicl::sum (magicl::.* normal mp-vel)))
                            (tang-vel (magicl:.- mp-vel (magicl:scale normal rel-vel)))
                            (normal-damping 0d0)
                            (damping-force (* normal-damping rel-vel)))
                       (declare (double-float rel-vel normal-damping damping-force))

                       (cl-mpm/fastmath::fast-fmacc force
                                                    normal
                                                    (- normal-force
                                                       damping-force))

                       (sb-thread:with-mutex (node-lock)
                         (cl-mpm/fastmath::fast-fmacc node-force
                                                      force
                                                      svp)
                         (cl-mpm/fastmath::fast-fmacc node-ext-force
                                                      force
                                                      svp))
                       ))))))))))))
