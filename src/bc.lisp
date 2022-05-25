(defpackage :cl-mpm
  (:use :cl))
(in-package :cl-mpm)
(defclass boundary-condition ()
  ( (position 
      :initarg :position)))

(defclass bc-fixed (boundary-condition)
  ((values
     :initarg :value
     :initform '(nil))))
(defun make-bc-fixed (pos value)
  (make-instance 'bc-fixed
                 :position pos
                 :value value))

(defgeneric apply-bc (bc node)
  (:documentation "Apply a boundary condition onto a node"))
(defmethod apply-bc ((bc bc-fixed-vel) node)
  (with-slots ((value value))
    bc
    (with-slots ((vel velocity)) node
      (loop for d from 0 to (length value)
            do (when (nth d value)
                 (setf (magicl:tref d 0) (nth d value)))))))

(defun make-outside-bc (size)
  (let (nD (length size))
    (loop for d from 0 to (- nD 1)
          append (loop for x from 0 to (nth d size)
                   collect (make-fixed-bc )
                   )
          )
    )
  )
