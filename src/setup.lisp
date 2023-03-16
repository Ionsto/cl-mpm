(defpackage :cl-mpm/setup
  (:use :cl)
  (:export
    #:make-column
    #:make-block-mps
    #:remove-sdf
    #:ellipse-sdf
    #:circle-sdf
    #:rectangle-sdf
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

(defun make-block-mps-list (offset size mps density constructor &rest args)
  (let*  ((nD 2)
          (spacing (mapcar #'/ size mps))
          (volume (* (first spacing) (second spacing)))
          (data (loop for x from 0 to (- (first mps) 1)
                      append
                      (loop
                        for y from 0 to (- (second mps) 1)
                        collect 
                        (let* ((i (+ y (* x (first mps)))))
                          (apply constructor (append '(2) args
                                                     (list
                                                      :position (list (+ (first offset) (* (first spacing) x))
                                                                 (+ (second offset) (* (second spacing) y)))
                                                      :mass (* density volume)
                                                      :volume volume
                                                      :size (magicl:from-list spacing (list nD 1) :type 'double-float))))
                          )))))
    data))
(defun make-mps-from-list (mp-list)
  (make-array (length mp-list) :initial-contents mp-list :adjustable t :fill-pointer (length mp-list)))
(defun make-block-mps (offset size mps constructor &rest args)
  (let*  ((data (apply #'make-block-mps-list offset size mps constructor args)))
    (make-mps-from-list data)))


(defun make-column-mps (size mp-spacing constructor &rest args)
  (let* ((mp-spacing-x (first mp-spacing))
        (mp-spacing-y (second mp-spacing))
        (data (loop for i from 0 to (- size 1) collect 
                    (apply constructor (append '(2) args 
                                               (list
                                                 :position (list (/ mp-spacing-x 2) 
                                                            (+ (/ mp-spacing-y 2) (* mp-spacing-y i)))
                                                 :volume (* mp-spacing-x mp-spacing-y)
                                                 :size (magicl:from-list (list mp-spacing-x mp-spacing-y) '(2 1))))))))
    (make-array (length data) :initial-contents data)))

(defun make-column-mps-elastic (element-count spacing E nu)
  (make-column-mps element-count spacing 'cl-mpm::make-particle-elastic E nu))

(defun node-to-mp-pos (position res mps)
  (mapcar #'+ (list (* res (+ (/ 1 (* 2 mps))))
                    (* res (+ (/ 1 (* 2 mps)))))
                      position))

(defun remove-sdf (sim sdf)
  (setf (cl-mpm:sim-mps sim)
        (lparallel:premove-if (lambda (mp)
                                (with-accessors ((pos cl-mpm/particle:mp-position)) mp
                                  (>= 0 (funcall sdf pos))
                                  ))
                              (cl-mpm:sim-mps sim))))
(defun rectangle-sdf (position size)
  (lambda (pos)
      (let* ((position (magicl:from-list position '(2 1) :type 'double-float))
             (dist-vec (magicl:.- (magicl:map! #'abs (magicl:.- pos position))
                                  (magicl:from-list size '(2 1) :type 'double-float))))

        (+ (sqrt (magicl::sum
                  (magicl:map! (lambda (x) (* x x))
                               (magicl:map! (lambda (x) (max 0d0 x)) dist-vec))))
           (min (max (magicl:tref dist-vec 0 0)
                     (magicl:tref dist-vec 1 0)
                     ) 0d0)))))
(defun ellipse-sdf (position x-l y-l)
  (let ((aspect (/ x-l y-l)))
    (lambda (pos)
      (let* ((position (magicl:from-list position '(2 1) :type 'double-float))
             (dist-vec (magicl:.* (magicl:.- position pos) (magicl:from-list (list 1 aspect) '(2 1)
                                                                             :type 'double-float)))
             (distance (sqrt (magicl:tref (magicl:@ (magicl:transpose dist-vec)
                                                    dist-vec) 0 0))))
        (- distance x-l)))))
(defun circle-sdf (position radius)
  (ellipse-sdf position radius radius))
