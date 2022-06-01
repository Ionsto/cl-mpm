(defpackage :cl-mpm/setup
  (:use :cl)
  (:export
    #:make-column
  ))

(in-package :cl-mpm/setup)

(defun make-column (size)
  (let* ((nD 2)
         (mp-spacing 1d0)
         (sim (cl-mpm:make-mpm-sim (list 1 size) 1 1e-3 
                                                  (cl-mpm::make-shape-function-linear nD))))
    (progn (setf (cl-mpm:sim-mps sim) 
                 (loop for i from 0 to (- size 1) collect 
                       (cl-mpm::make-particle-elastic nD
                                                      1e2
                                                      0
                                                      :pos (list 0.5 (+ 0.5d0 (* mp-spacing i)))
                                                      :volume mp-spacing)))
          (setf (cl-mpm:sim-bcs sim) (cl-mpm/bc:make-outside-bc (cl-mpm:mesh-mesh-size (cl-mpm:sim-mesh sim)))) 
           sim)))
