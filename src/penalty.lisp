
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
         (yheight (cl-mpm/fastmath::dot (magicl:scale (cl-mpm/particle::mp-domain-size mp) 0.5d0) normal))
        (dist (- datum (- ypos yheight))))
    (the double-float dist)))
(defun penetration-point (mp pen datum normal)
  "Get linear penetration distance"
  (let* ((pos (cl-mpm/particle:mp-position mp))
         (domain (cl-mpm/particle::mp-domain-size mp)))
    (magicl:.- pos (magicl.simd::.*-simd normal (magicl:scale domain 0.5d0)))))

;;Only vertical condition
(defun apply-force-mps (mesh mps normal datum epsilon friction)
  "Update force on nodes, with virtual stress field from mps"
  (lparallel:pdotimes
      (i (length mps))
    (let ((mp (aref mps i)))

      (let* ((penetration-dist (penetration-distance mp datum normal)))
        (declare (double-float penetration-dist))
        (when (> penetration-dist 0d0)
          (with-accessors ((volume cl-mpm/particle:mp-volume)
                           (pressure cl-mpm/particle::mp-pressure)
                           (mp-vel cl-mpm/particle::mp-velocity)
                           )
              mp
            (let ((pen-point (penetration-point mp penetration-dist datum normal)))
              ;; (format t "Contact point ~F : dist ~F~%" (magicl::storage pen-point) penetration-dist)
              ;;Iterate over neighbour nodes
              (cl-mpm::iterate-over-neighbours-point-linear-simd
               mesh pen-point
               (lambda (mesh node svp grads)
                 (with-accessors ((node-force cl-mpm/mesh:node-force)
                                  (node-lock  cl-mpm/mesh:node-lock)
                                  (node-vel  cl-mpm/mesh:node-velocity)
                                  (node-boundary cl-mpm/mesh::node-boundary-node)
                                  (node-boundary-scalar cl-mpm/mesh::node-boundary-scalar)
                                  (node-active  cl-mpm/mesh:node-active))
                     node
                   (declare (double-float volume svp))
                   (when node-active
                     ;;Lock node for multithreading
                     (sb-thread:with-mutex (node-lock)
                       ;; (cl-mpm/shape-function::assemble-dsvp-2d-prealloc grads dsvp)
                       (let* ((force (cl-mpm/utils:vector-zeros))
                              (normal-force (* (expt penetration-dist 1) epsilon))
                              ;; (reaction-force (magicl:scale normal normal-force))
                              (rel-vel (magicl::sum (magicl::.* normal mp-vel)))
                              (tang-vel (magicl:.- mp-vel (magicl:scale normal rel-vel)))
                              (normal-damping 0d10)
                              (damping-force (* normal-damping rel-vel))
                              )
                         ;; (cl-mpm/fastmath::fast-add force reaction-force)
                         ;; (when (> rel-vel 1d0))
                         ;; (magicl:scale! node-vel 0d0)
                         ;; (magicl:.-
                         ;;  node-vel
                         ;;  (magicl:scale!
                         ;;   (magicl:.-
                         ;;    node-vel
                         ;;    (magicl:scale
                         ;;     normal
                         ;;     (cl-mpm/fastmath::dot
                         ;;      node-vel
                         ;;      normal)))
                         ;;   (* svp))
                         ;;  node-vel)
                         ;;  (magicl:.-
                         ;;   node-force
                         ;;   (magicl:scale!
                         ;;    (magicl:.-
                         ;;     node-force
                         ;;     (magicl:scale
                         ;;      normal
                         ;;      (cl-mpm/fastmath::dot
                         ;;       node-force
                         ;;       normal)))
                         ;;    (* svp))
                         ;;   node-force)
                         (cl-mpm/fastmath::fast-fmacc force
                                                      normal
                                                      (- normal-force
                                                         damping-force))
                         (cl-mpm/fastmath::fast-fmacc force
                                                      tang-vel
                                                      (* -1d0
                                                         friction
                                                         ))
                         ;; (cl-mpm/fastmath::fast-fmacc node-vel
                         ;;                              tang-vel
                         ;;                              (* -1d0))
                         (cl-mpm/fastmath::fast-fmacc node-force
                                                      force
                                                      svp
                                                      ;; (* svp volume)
                                                      )
                         )))))))))))))
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
                   (mps cl-mpm::sim-mps))
      sim
    (apply-force-mps mesh mps normal datum epsilon friction)))

(defclass bc-penalty (cl-mpm/bc::bc-closure)
  ((normal
    :accessor bc-penalty-normal
    :initarg :normal
    )
   (datum
    :accessor bc-penalty-datum
    :initarg :datum
    ))
  (:documentation "A nonconforming neumann bc"))

(defun make-bc-penalty (sim datum epsilon friction)
  (let ((normal (cl-mpm/fastmath::norm (magicl:from-list '(0d0 1d0) '(2 1)))))
    (make-instance 'cl-mpm/bc::bc-closure
                   :index '(0 0)
                   :func (lambda ()
                           (progn
                             (apply-penalty
                              sim
                              normal
                              datum
                              epsilon
                              friction))))))

(defun make-bc-penalty-point-normal (sim normal point epsilon friction)
  (let* ((normal (cl-mpm/fastmath::norm normal))
         ;(datum (cl-mpm/fastmath::dot normal point))
         (datum (- (penetration-distance-point point 0d0 normal)))
         )
    (format t "Normal ~F ~F ~%" (magicl:tref normal 0 0) (magicl:tref normal 1 0))
    (make-instance 'bc-penalty
                   :index '(0 0)
                   :normal normal
                   :datum datum
                   :func (lambda ()
                           (progn
                             (apply-penalty
                              sim
                              normal
                              datum
                              epsilon
                              friction))))))
