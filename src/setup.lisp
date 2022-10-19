(defpackage :cl-mpm/setup
  (:use :cl)
  (:export
    #:make-column
  ))

(in-package :cl-mpm/setup)




(defun make-block (res element-count &optional (shape-maker #'cl-mpm::make-shape-function-linear))
  "Make a 2D column of heigh size, and width 1 - filled with elements"
  (let* ((nD 2)
         (size (mapcar (lambda (x) (* x res)) element-count))
         (sim (cl-mpm:make-mpm-sim size res 1e-3 (funcall shape-maker nD res))))
    (progn 
          (setf (cl-mpm:sim-mps sim) #())
          (setf (cl-mpm:sim-bcs sim) (cl-mpm/bc:make-outside-bc (cl-mpm/mesh:mesh-count (cl-mpm:sim-mesh sim)))) 
           sim)))

(defun make-column (height element-count &optional (shape-maker #'cl-mpm::make-shape-function-linear))
  "Make a 2D column of heigh size, and width 1 - filled with elements"
  (let* ((nD 2)
         (mp-spacing (/ height element-count))
         (sim (cl-mpm:make-mpm-sim (list mp-spacing height) mp-spacing 1e-3
                                                  (funcall shape-maker nD mp-spacing))))
    (progn 
          (setf (cl-mpm:sim-mps sim) #())
          (setf (cl-mpm:sim-bcs sim) (cl-mpm/bc:make-outside-bc (cl-mpm/mesh:mesh-count (cl-mpm:sim-mesh sim)))) 
           sim)))

(defun make-block-mps (offset size mps constructor &rest args)
  (let*  ((nD 2)
          (spacing (mapcar #'/ size mps))
          (data (loop for x from 0 to (- (first mps) 1)
                      append
                      (loop
                        for y from 0 to (- (second mps) 1)
                        collect 
                        (let* ((i (+ y (* x (first mps)))))
                          (apply constructor (append '(2) args
                                                     (list
                                                      :pos (list (+ (first offset) (* (first spacing) x))
                                                                 (+ (second offset) (* (second spacing) y)))
                                                      :volume (* (first spacing) (second spacing)))))
                          )))))
    (make-array (length data) :initial-contents data)))


(defun make-column-mps (size mp-spacing constructor &rest args)
  (let* ((mp-spacing-x (first mp-spacing))
        (mp-spacing-y (second mp-spacing))
        (data (loop for i from 0 to (- size 1) collect 
                    (apply constructor (append '(2) args 
                                               (list
                                                 :pos (list (/ mp-spacing-x 2) 
                                                            (+ (/ mp-spacing-y 2) (* mp-spacing-y i)))
                                                 :volume (* mp-spacing-x mp-spacing-y)))))))
    (make-array (length data) :initial-contents data)))

(defun make-column-mps-elastic (element-count spacing E nu)
  (make-column-mps element-count spacing 'cl-mpm::make-particle-elastic E nu))
