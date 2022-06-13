(defpackage :cl-mpm/bc
  (:use :cl)
  (:export
    #:apply-bc
    #:bc-index
    #:bc-value
    #:make-outside-bc
    )
  )
(in-package :cl-mpm/bc)
(defclass bc ()
  ( (index 
      :accessor bc-index
      :initarg :index))
  (:documentation "A boundary condition that applies some operation at an index"))

(defclass bc-fixed (bc)
  ((value
     :accessor bc-value
     :initarg :value
     :initform '(nil nil)))
  (:documentation "A fixed velocity BC can be 1/0 dof"))

(defun make-bc-fixed (index value)
  (make-instance 'bc-fixed
                 :index index 
                 :value value))

(defgeneric apply-bc (bc node)
  (:documentation "Apply a boundary condition onto a node"))

(defmethod apply-bc (bc node)
  "Natural BC condition is nothing")

(defmethod apply-bc ((bc bc-fixed) node)
  "Fixed velocity BC over some dimensions"
  (with-slots ((value value))
    bc
    (loop for d from 0 to (length value)
            do (when (nth d value)
                 (setf (magicl:tref (cl-mpm/mesh:node-velocity node) d 0) (nth d value))))))

(defun make-outside-bc (mesh-count)
  "Construct fixed bcs over the outside of a mesh"
  (destructuring-bind (xsize ysize) (mapcar (lambda (x) (- x 1)) mesh-count)
    (append 
      (loop for x from 0 to xsize 
            append 
            (list (make-bc-fixed (list x 0)     '(nil 0d0))
                  (make-bc-fixed (list x ysize) '(nil 0d0))))
       (loop for y from 0 to ysize 
            append 
            (list (make-bc-fixed (list 0     y) '(0d0 nil))
                  (make-bc-fixed (list xsize y) '(0d0 nil)))))))
