(defpackage :cl-mpm/dynamic-relaxation-erosion
  (:use :cl :cl-mpm/dynamic-relaxation)
  (:export))
(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)

(define-condition error-erosion-criteria (non-convergence-error)
  ((max-erosion-inc :initarg :max-erosion-inc :reader max-erosion-inc)))

(defmethod convergence-check :after ((sim cl-mpm::mpm-sim))
  (cl-mpm:iterate-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (when (typep mp 'cl-mpm/particle::particle-erosion)
       (let ((inc (/ (- (cl-mpm/particle::mp-eroded-volume mp)
                        (cl-mpm/particle::mp-eroded-volume-n mp))
                     (cl-mpm/particle::mp-mass mp))))
         (when (> inc 0.5d0)
           (format t "Erosion criteria exceeded~%")
           (error (make-instance 'error-erosion-criteria :max-erosion-inc inc))))))))
