(defpackage :cl-mpm/particle
  (:use :cl)
  (:export
    #:make-particle
    #:make-particle-elastic
    #:mp-mass
    #:mp-nd
    #:mp-volume
    #:mp-position
    #:mp-velocity
    #:mp-stress
    #:mp-strain
    #:mp-strain-rate
    #:mp-gravity
    #:mp-deformation-gradient
    #:constitutive-model
    )
  )
(in-package :cl-mpm/particle)

(defclass particle ()
  ((mass 
     :accessor mp-mass
     :type double-float
     :initarg :mass 
     :initform 1d0)
   (nD 
     :accessor mp-nd
     :type integer
     :initarg :nD)
   (volume 
     :accessor mp-volume
     :type double-float
     :initarg :volume 
     :initform 1d0)
   (position 
     :accessor mp-position
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :initarg :position)
   (velocity 
     :accessor mp-velocity
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :initarg :velocity
     :initform (magicl:zeros '(2 1)))
   (stress 
     :accessor mp-stress
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :initarg :stress
     :initform (magicl:zeros '(3 1)))
   (strain 
     :accessor mp-strain
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :accessor mp-strain
     :initarg :strain
     :initform (magicl:zeros '(3 1)))
   (strain-rate 
     :accessor mp-strain-rate
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :accessor mp-strain-rate
     :initarg :strain-rate
     :initform (magicl:zeros '(3 1)))
   (deformation-gradient
     :accessor mp-deformation-gradient
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :initarg :deformation-gradient
     :initform (magicl:eye 2))
   (gravity
     :type double-float
     :accessor mp-gravity
     :initform -9.8d0))
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
  (:documentation "A linear-elastic material point"))

(defun make-particle (nD &key (pos nil) (volume 1))
  (progn
    (if (eq pos nil)
        (setf pos (magicl:zeros (list nD 1)))
        (setf pos (magicl:from-list (mapcar (lambda (x) (coerce x 'double-float)) pos) (list nD 1))))
    (let ((stress-size 3))
      (make-instance 'particle
                     :nD nD
                     :volume (coerce volume 'double-float)
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

(defmethod constitutive-model ((mp particle-elastic) strain)
    (with-slots ((E E)
                 (nu nu))
                mp
        (cl-mpm/constitutive:linear-elastic strain E nu)))



