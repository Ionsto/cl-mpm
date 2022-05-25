(defpackage :cl-mpm
  (:use :cl)
  )
(in-package :cl-mpm)

(defmacro shape-linear (x)
  `(- 1d0 (abs ,x)))

(defmacro create-svp (arg form)
  `(lambda (,arg) ,form))

(defmacro create-dsvp (arg form)
  (let ((dx (derive arg form)))
    `(lambda (,arg) ,dx)))


(defun svp-1d (svp dsvp)
  svp)
(defun svp-2d (svp dsvp)
  (lambda (x y) (* (funcall svp x) (funcall svp y))))
(defun svp-3d (svp dsvp)
  (lambda (x y z) (* (funcall svp x) (funcall svp y) (funcall svp y))))

(defun dsvp-1d (svp dsvp)
  (lambda (x) (list (funcall dsvp x))))
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

;(defun assemble-dsvp (dsvp)
;  "Assemble d/di to the strain-displacement matrix"
;  (let ((nD 2)
;        (dx (aref dsvp 0))
;        (dy (aref dsvp 1)))
;    (magicl:from-list (list dx 0d0 0d0 dy dy 0d0) '(3 2) :type 'double-float)))
(defun assemble-dsvp-1d (dsvp)
  "Assemble d/di to the strain-displacement matrix"
  (let ((dx (nth 0 dsvp)))
    (magicl:from-list (list dx) '(1 1) :type 'double-float)))

(defun assemble-dsvp-2d (dsvp)
  "Assemble d/di to the strain-displacement matrix"
  (let ((dx (nth 0 dsvp))
        (dy (nth 1 dsvp)))
    (magicl:from-list (list dx 0d0 0d0 dy dy 0d0) '(3 2) :type 'double-float)))

(defun assemble-dsvp-3d (dsvp)
  "Assemble d/di to the strain-displacement matrix"
  (let ((dx (nth 0 dsvp))
        (dy (nth 1 dsvp))
        (dz (nth 2 dsvp)))
    (magicl:from-list (list dx  0d0 0d0;xx
                            0d0 dy  0d0;yy
                            0d0 0d0  dz;zz
                            0d0 dz  0d0;yz
                            dz  0d0 0d0;xz
                            dx  0d0 0d0;xy
                            ) '(6 3) :type 'double-float)))

(defun assemble-dsvp (nD dsvp)
  (case nD
    (1 (assemble-dsvp-1d dsvp))
    (2 (assemble-dsvp-2d dsvp))
    (3 (assemble-dsvp-3d dsvp))))


