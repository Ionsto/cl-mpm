(defpackage :cl-mpm/setup
  (:use :cl)
  (:export
    #:make-column
  ))

(in-package :cl-mpm/setup)

(defun make-column (size)
  "Make a 2D column of heigh size, and width 1 - filled with elements"
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

(defun make-block-mps (offset size mps)
  (let*  ((nD 2)
          (spacing (mapcar #'/ size mps)))
    (loop for x from 0 to (- (first mps) 1)
          for y from 0 to (- (second mps) 1)
          collect 
          (let* ((i (+ y (* x (first mps)))))
            (cl-mpm::make-particle-elastic nD
                                           1e2
                                           0
                                           :pos (list (+ (first offset) (* (first spacing) x))
                                                      (+ (second offset) (* (second spacing) i)))
                                           :volume (* (first spacing) (second spacing)))))))
