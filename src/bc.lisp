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
      :initarg :index)))

(defclass bc-fixed (bc)
  ((value
     :accessor bc-value
     :initarg :value
     :initform '(nil nil))))
(defun make-bc-fixed (index value)
  (make-instance 'bc-fixed
                 :index index 
                 :value value))

(defgeneric apply-bc (bc node)
  (:documentation "Apply a boundary condition onto a node"))

(defmethod apply-bc ((bc bc-fixed) node)
  (with-slots ((value value))
    bc
    (loop for d from 0 to (length value)
            do (when (nth d value)
                 (setf (magicl:tref (cl-mpm:node-velocity node) d 0) (nth d value))))))

(defun make-outside-bc (size)
  (destructuring-bind (xsize ysize) size
    (append 
      (loop for x from 0 to xsize 
            append 
            (list (make-bc-fixed (list x 0)     '(nil 0d0))
                  (make-bc-fixed (list x ysize) '(nil 0d0))))
       (loop for y from 0 to ysize 
            append 
            (list (make-bc-fixed (list 0     y) '(0d0 nil))
                  (make-bc-fixed (list xsize y) '(0d0 nil)))))))
