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
   (strain 
     :accessor mp-strain
     :initarg :strain
     :initform (magicl:zeros '(3 1)))
   (deformation-matrix 
     :initarg :deformation-gradient
     :initform (magicl:eye 2)))
  (:documentation "A single material point"))

(defclass particle-elastic (particle)
  ((E
     :accessor mp-E
     :initarg :E
     )
   (nu
     :accessor mp-nu
     :initarg :nu)
   )
  )

(defun make-particle (nD &key (pos nil) (volume 1))
  (progn
    (if (eq pos nil)
        (setf pos (magicl:zeros (list nD 1)))
        (setf pos (magicl:from-list pos (list nD 1))))
    (let ((stress-size 3))
      (make-instance 'particle
                     :nD nD
                     :volume volume
                     :velocity (magicl:zeros (list 2 1))
                     :deformation-gradient (magicl:eye 2)
                     :stress (magicl:zeros '(3 1))
                     :position pos))))
(defun make-particle-elastic (nD E nu &key (pos nil) (volume 1))
  (let ((p (make-particle nD :pos pos :volume volume)))
    (progn
      (change-class p 'particle-elastic)
      (setf (mp-E p) E)
      (setf (mp-nu p) nu)
      p)))

(defgeneric constitutive-model (mp elastic-trial-strain)  
 (:documentation "Compute new stress state given elastic strain")
 (:method (mp strain)
  (magicl:scale strain 0)))
