
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


(defun pentration-distance (mp datum)
  (let* ((ypos (magicl:tref (cl-mpm/particle:mp-position mp) 1 0))
         (yheight (* 0.5d0 (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0)))
        (dist (- datum (- ypos yheight))))
    dist
    ))

;;Only vertical condition
(defun apply-force-mps (mesh mps datum epsilon)
  "Update force on nodes, with virtual stress field from mps"
  (lparallel:pdotimes (i (length mps))
    (let ((mp (aref mps i)))

      (let* ((pos (cl-mpm/particle:mp-position mp))
            (ypos (magicl:tref (cl-mpm/particle:mp-position mp) 1 0))
            (yheight (* 1.0d0 0.5d0 (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0)))
            (penetration-dist (- datum (- ypos yheight))));(pentration-distance mp datum)))
        (declare (double-float penetration-dist))
        (when (> penetration-dist 0d0)
          (with-accessors ((volume cl-mpm/particle:mp-volume)
                           (pressure cl-mpm/particle::mp-pressure)
                           )
              mp
            (let ((dsvp (cl-mpm/utils::stretch-dsvp-zeros)))
              ;;Iterate over neighbour nodes
              (cl-mpm::iterate-over-neighbours-point-linear
               mesh (.+ pos (magicl:from-list (list 0d0 (- yheight)) '(2 1)))
               (lambda (mesh node svp grads)
                 (with-accessors ((node-force cl-mpm/mesh:node-force)
                                  (node-lock  cl-mpm/mesh:node-lock)
                                  (node-boundary cl-mpm/mesh::node-boundary-node)
                                  (node-boundary-scalar cl-mpm/mesh::node-boundary-scalar)
                                  (node-active  cl-mpm/mesh:node-active))
                     node
                   (declare (double-float volume svp))
                   (when node-active
                     ;;Lock node for multithreading
                     (sb-thread:with-mutex (node-lock)
                       ;; (cl-mpm/shape-function::assemble-dsvp-2d-prealloc grads dsvp)
                       (let ((force (cl-mpm/utils:vector-zeros)))
                         (setf (tref force 1 0) (+ (* penetration-dist epsilon)
                                                   ;; (- (* pressure volume))
                                                   0d0
                                                   ))
                         (cl-mpm/fastmath::fast-fmacc node-force
                                                      force
                                                      ;; svp
                                                      (* svp volume)
                                                      ))))))))))))))

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
(defun apply-penalty (sim datum epsilon)
    (with-accessors ((mesh cl-mpm:sim-mesh)
                     (mps cl-mpm::sim-mps))
        sim
      (apply-force-mps mesh mps datum epsilon)))

(defclass bc-penalty (bc)
  ()
  (:documentation "A nonconforming neumann bc"))

(defun make-bc-penalty (sim datum epsilon)
  (make-instance 'cl-mpm/bc::bc-closure
                 :index '(0 0)
                 :func (lambda ()
                         (progn
                           (apply-penalty
                            sim
                            datum
                            epsilon)))))
