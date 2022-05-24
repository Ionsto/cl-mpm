(defpackage :cl-mpm
  (:use :cl))
(in-package :cl-mpm)

(defclass particle ()
  ((mass :initarg :mass 
         :initform 1)
   (volume :initarg :volume 
           :initform 1)
   (position :initarg :position)
   (velocity :initarg velocity
             :initform (magicl:zeros '(2 1)))
   (stress :initform (magicl:zeros '(3 1)))
   (strain-rate :initform (magicl:zeros '(3 1)))
   (deformation-matrix :initform (magicl:eye 2)))
  (:documentation "A single material point"))

(defun make-particle (&key (x 0) (y 0) (volume 1))
  (make-instance 'particle
                 :volume volume
                 :position (magicl:from-list (list x y) '(2 1))))

(defgeneric constitutive-model (mp elastic-trial-strain)  
  (:documentation "Compute new stress state given elastic strain"))
