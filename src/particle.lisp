(defpackage :cl-mpm/particle
  (:use :cl)
  (:export
    #:make-particle
    #:make-particle-elastic
    #:make-particle-elastic-damage
    #:mp-mass
    #:mp-nd
    #:mp-volume
    #:mp-position
    #:mp-velocity
    #:mp-stress
    #:mp-strain
    #:mp-strain-rate
    #:mp-gravity
    #:mp-body-force
    #:mp-damage
    #:mp-critical-stress
    #:mp-deformation-gradient
    #:constitutive-model
    #:particle
    #:particle-damage
    #:post-stress-step
    )
  )
(in-package :cl-mpm/particle)
(declaim (optimize (debug 0) (safety 0) (speed 3)))

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
   (size
     :accessor mp-domain-size
     :type magicl:matrix/double-float
     :initarg :size
     :initform (magicl:zeros '(2 1)))
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
     :initarg :strain
     :initform (magicl:zeros '(3 1)))
   (strain-plastic
     :accessor mp-strain-plastic
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :initarg :strain-plastic
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
     :initform -9.8d0)
   (body-force
     :accessor mp-body-force
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :initarg :body-force
     :initform (magicl:zeros '(2 1)))
   )
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

(defclass particle-damage (particle)
  (
   (damage
    :accessor mp-damage
    :type DOUBLE-FLOAT
    :initarg :damage
    :initform 0)
   (critical-stress
    :accessor mp-critical-stress
    :type DOUBLE-FLOAT
    :initarg :critical-stress
    :initform 0)
   )
  (:documentation "A material point with a damage tensor"))
(defclass particle-fracture (particle-damage)
  (
   (strain-energy-density
    :accessor mp-strain-energy-density
    :type DOUBLE-FLOAT
    :initform 0)
   (fracture-toughness
    :accessor mp-fracture-toughness
    :type DOUBLE-FLOAT
    :initform 0
    :initarg :fracture-toughness)
   )
  (:documentation "A material point with fracture mechanics"))

(defclass particle-elastic-damage (particle-elastic particle-damage)
  ()
  (:documentation "A mp with damage influanced elastic model"))
(defclass particle-elastic-fracture (particle-elastic particle-fracture)
  ()
  (:documentation "A mp with fracture mechanics"))

(defun make-particle (nD &optional (constructor 'particle) &rest args &key  (position nil) (volume 1) (mass 1) &allow-other-keys)
  (progn
    (if (eq position nil)
        (setf position (magicl:zeros (list nD 1)))
        (setf position (magicl:from-list (mapcar (lambda (x) (coerce x 'double-float)) position) (list nD 1))))
    (let ((stress-size 3))
      (apply #'make-instance constructor
                     :nD nD
                     :volume (coerce volume 'double-float)
                     :mass (coerce mass 'double-float)
                     :position position
                     args))))
;; (defun make-particle (nD &rest args &key (constructor 'particle) (pos nil) (volume 1) (mass 1))
;;   (progn
;;     (if (eq pos nil)
;;         (setf pos (magicl:zeros (list nD 1)))
;;         (setf pos (magicl:from-list (mapcar (lambda (x) (coerce x 'double-float)) pos) (list nD 1))))
;;     (let ((stress-size 3))
;;       (make-instance constructor
;;                      :nD nD
;;                      :volume (coerce volume 'double-float)
;;                      :mass (coerce mass 'double-float)
;;                      :position pos))))
(defun make-particle-elastic (nD E nu &key (pos nil) (volume 1) (mass 1))
    (let ((p (make-particle nD 'particle-elastic :position pos :volume volume :mass mass) ))
        (progn
            (setf (mp-E p) E)
            (setf (mp-nu p) nu)
        p)))

(defun make-particle-elastic-damage (nD E nu &key (pos nil) (volume 1) (mass 1))
    (let ((p (make-particle nD 'particle-elastic-damage :position pos :volume volume :mass mass
                            ) ))
        (progn
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

(defgeneric post-stress-step (mesh mp dt)
  (:documentation "This step gets called after full stress state resolved and allows for other processing"))
(defmethod post-stress-step (mesh mp dt)
  ())


