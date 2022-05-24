(defpackage :cl-mpm
  (:use :cl)
  (:export 
    ;#:make-shape-function
    ;#:make-shape-function-linear
    )
  )
(in-package :cl-mpm)

(defmacro shape-linear (x)
  `(- 1d0 (abs ,x)))

(defmacro create-svp (arg form)
  `(lambda (,arg) ,form))

(defmacro create-dsvp (arg form)
  (let ((dx (derive arg form)))
    `(lambda (,arg) ,dx)))


(let* (
       (svp (create-svp x (- 1 (abs x))))
       (dsvp (create-dsvp x (- 1 (abs x))))
       (svp-bi (lambda (x y) (* (funcall svp x) (funcall svp y))))
      ))

(defun svp-1d (svp dsvp)
  svp)
(defun svp-2d (svp dsvp)
  (lambda (x y) (* (funcall svp x) (funcall svp y))))
(defun svp-3d (svp dsvp)
  (lambda (x y z) (* (funcall svp x) (funcall svp y) (funcall svp y))))

(defun dsvp-1d (svp dsvp)
  dsvp)
(defun dsvp-2d (svp dsvp)
  (lambda (x y) (let (
                      (wx (funcall svp x))
                      (wy (funcall svp y))
                      (dx (funcall dsvp x))
                      (dy (funcall dsvp y))
                      )
                  ;(d/dx d/dy)
                  (list (* wy dx) (* wx dy)))))
(defun dsvp-3d (svp dsvp)
  (lambda (x y z) 
    (let (
          (wx (funcall svp x))
          (wy (funcall svp y))
          (wz (funcall svp z))
          (dx (funcall dsvp x))
          (dy (funcall dsvp y))
          (dz (funcall dsvp z))
          )
      ;(d/dx d/dy d/dz)
      (list (* wy wz dx) (* wx wz dy) (* wx wy dz)))))

(defun nd-svp (nD svp dsvp)
  (case  nD
    (1 (svp-1d svp dsvp))
    (2 (svp-2d svp dsvp))
    (3 (svp-3d svp dsvp))))

(defun nd-dsvp (nD svp dsvp)
  (case  nD
    (1 (dsvp-1d svp dsvp))
    (2 (dsvp-2d svp dsvp))
    (3 (dsvp-3d svp dsvp))))

(defparameter *svp* (create-svp x (- 1 (abs x))))
(defparameter *dsvp* (create-dsvp x (- 1 (abs x))))
(funcall (nd-dsvp 3 *svp* *dsvp*) 0 0 0)

(defclass shape-function ()
  ((nD 
     :accessor nD
     :initarg :nD)
   (order 
     :accessor order
     :initarg :order)
   (svp 
     :accessor svp
     :initarg :svp)
   (dsvp 
     :accessor dsvp
     :initarg :dsvp)
   ))

(defclass shape-function-linear (shape-function)
  ((order :initform 2)
   (svp :initform (autodiff:lambda-ad (x y) (* (shape-linear x) (shape-linear y))))))

(defclass shape-function-linear (shape-function)
  ((order :initform 1)))

(defmacro make-shape-function (arg shape-form nD order)
  `(let ((svp (create-svp ,arg ,shape-form))
         (dsvp (create-dsvp ,arg ,shape-form)))
     (make-instance 'shape-function
                    :nD ,nD 
                    :order ,order 
                    :svp (nd-svp ,nD svp dsvp)
                    :dsvp (nd-dsvp ,nD svp dsvp))))
(defun make-shape-function-linear (nD)
  (make-shape-function x (- 1 (abs x)) nD 1))
