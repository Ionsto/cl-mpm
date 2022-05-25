(defpackage :cl-mpm
  (:use :cl))
(in-package :cl-mpm)

(defclass particle ()
  ((mass :initarg :mass 
         :initform 1)
   (nD :initarg :nD)
   (volume :initarg :volume 
           :initform 1)
   (position :initarg :position)
   (velocity 
     :initarg :velocity
     :initform (magicl:zeros '(2 1)))
   (stress 
     :initarg :stress
     :initform (magicl:zeros '(3 1)))
   (strain-rate 
     :initarg :strain-rate
     :initform (magicl:zeros '(3 1)))
   (deformation-matrix 
     :initarg :deformation-gradient
     :initform (magicl:eye 2)))
  (:documentation "A single material point"))

(defun make-particle (nD &key (pos nil) (volume 1))
  (progn
    (print pos)
    (if (eq pos nil)
        (setf pos (magicl:zeros (list nD 1)))
        (setf pos (magicl:from-list pos (list nD 1))))
    (print pos)
    (let ((stress-size (floor (+ nD (/ (- (expt nD 2) nD) 2)))))
      (make-instance 'particle
                     :nD nD
                     :volume volume
                     :velocity (magicl:zeros (list nD 1))
                     :deformation-gradient (magicl:zeros (list nD nD))
                     :stress (magicl:zeros (list stress-size 1))
                     :strain-rate (magicl:zeros (list stress-size 1))
                     :position pos))))

(defgeneric constitutive-model (mp elastic-trial-strain)  

  (:documentation "Compute new stress state given elastic strain"))
