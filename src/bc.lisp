(defpackage :cl-mpm/bc
  (:use :cl)
  (:export
    #:apply-bc
    #:bc-index
    #:bc-value
    #:make-outside-bc
    #:make-outside-bc-nostick
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

(defclass bc-surface (bc)
  ((normal
     :accessor bc-normal
     :initarg :normal
     :initform (magicl:zeros '(2 1))
     :type MAGICL:MATRIX/DOUBLE-FLOAT))
  (:documentation "A BC applied on a normal"))

(defun make-bc-surface (index normal)
  (make-instance 'bc-surface
                 :index index 
                 :normal normal))

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

(defmethod apply-bc ((bc bc-surface) node)
  "Fixed velocity BC over some non-stick surface"
  (with-slots ((normal normal))
    bc
    (with-accessors ((node-vel cl-mpm/mesh:node-velocity)) node
      (let ((rel-vel (magicl:tref (magicl:@ (magicl:transpose normal) node-vel) 0 0)))
        (when (< rel-vel 0)
          (setf node-vel (magicl:.- (magicl:scale normal rel-vel) node-vel)))))))

(defun make-outside-bc (mesh-count)
  "Construct fixed bcs over the outside of a mesh"
  (destructuring-bind (xsize ysize) (mapcar (lambda (x) (- x 1)) mesh-count)
    (append 
      (loop for x from 0 to xsize 
            append 
            (list (make-bc-fixed (list x 0)     '(0d0 0d0))
                  (make-bc-fixed (list x ysize) '(0d0 0d0))))
       (loop for y from 0 to ysize 
            append 
            (list (make-bc-fixed (list 0     y) '(0d0 0d0))
                  (make-bc-fixed (list xsize y) '(0d0 0d0)))))))

(defun make-outside-bc-nostick (mesh-count)
    "Construct nostick bcs over the outside of a mesh"
    (destructuring-bind (xsize ysize) (mapcar (lambda (x) (- x 1)) mesh-count)
        (append 
            (loop for x from 0 to xsize 
                append 
                (list (make-bc-surface (list x 0)     (magicl:from-list '(0d0  1d0) '(2 1)))
                      (make-bc-surface (list x ysize) (magicl:from-list '(0d0 -1d0) '(2 1)))))
            (loop for y from 0 to ysize 
                append 
                (list (make-bc-surface (list 0     y) (magicl:from-list '( 1d0 0d0) '(2 1)))
                      (make-bc-surface (list xsize y) (magicl:from-list '(-1d0 0d0) '(2 1))))))))
