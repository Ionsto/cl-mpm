(defpackage :cl-mpm/dynamic-relaxation-erosion
  (:use :cl :cl-mpm/dynamic-relaxation)
  (:export))
(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)

(define-condition error-erosion-criteria (non-convergence-error)
  ((max-erosion-inc :initarg :max-erosion-inc :reader max-erosion-inc)))

;; (defparameter *max-erosion-per-step* )

(defmethod convergence-check :after ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr))
  (let ((max-erosion-inc nil))
    (cl-mpm::iterate-over-bcs-force-serial
     sim
     (lambda (bc)
       (when (typep bc 'cl-mpm/erosion::bc-erode)
         (unless max-erosion-inc
           (setf max-erosion-inc (cl-mpm/erosion::bc-erode-max-erosion bc)))
         (setf max-erosion-inc
               (min (cl-mpm/erosion::bc-erode-max-erosion bc)
                    max-erosion-inc)))))
    (when max-erosion-inc
      (cl-mpm:iterate-over-mps
       (cl-mpm:sim-mps sim)
       (lambda (mp)
         (when (typep mp 'cl-mpm/particle::particle-erosion)
           (let ((inc (/ (cl-mpm/particle::mp-eroded-increment mp)
                         (cl-mpm/particle::mp-mass mp))))
             (when (> inc max-erosion-inc)
               (format t "Erosion criteria exceeded~%")
               (error (make-instance 'error-erosion-criteria :max-erosion-inc inc))))))))))
