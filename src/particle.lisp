(defpackage :cl-mpm/particle
  (:use :cl
        :cl-mpm/utils)
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
    #:mp-vorticity
    #:mp-gravity
    #:mp-body-force
    #:mp-damage
    #:mp-temperature
    #:mp-heat-capacity
    #:mp-critical-stress
    #:mp-deformation-gradient
    #:constitutive-model
    #:particle
    #:particle-damage
    #:post-stress-step
    )
  )
(in-package :cl-mpm/particle)
(declaim (optimize (debug 3) (safety 3) (speed 2)))

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
   (vorticity
    :accessor mp-vorticity
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (magicl:zeros '(3 1)))
   (deformation-gradient
     :accessor mp-deformation-gradient
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :initarg :deformation-gradient
     :initform (magicl:eye 2))
   (gravity
     :type double-float
     :accessor mp-gravity
     :initform 0d0;-9.8d0
     :initarg :gravity
     )
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

(defclass particle-fluid (particle)
  ((rest-density
    :accessor mp-rest-density
    :initarg :rest-density
    )
   (stiffness
    :accessor mp-stiffness
    :initarg :stiffness)
   (adiabatic-index
    :accessor mp-adiabatic-index
    :initarg :adiabatic-index)
   (viscosity
    :accessor mp-viscosity
    :initarg :viscosity)
   )
  (:documentation "A fluid material point"))

(defclass particle-viscoelastic (particle)
  ((E
    :accessor mp-E
    :initarg :E
    )
   (nu
    :accessor mp-nu
    :initarg :nu)
   )
  (:documentation "A visco-elastic material point"))

(defclass particle-viscoplastic (particle)
  ((E
    :accessor mp-E
    :initarg :E
    )
   (nu
    :accessor mp-nu
    :initarg :nu)
   (visc-factor
    :accessor mp-visc-factor
    :initarg :visc-factor)
   (visc-power
    :accessor mp-visc-power
    :initarg :visc-power)
   )
  (:documentation "A visco-plastic material point"))

(defclass particle-damage (particle)
  (
   (damage
    :accessor mp-damage
    :type DOUBLE-FLOAT
    :initarg :damage
    :initform 0d0)
   (critical-stress
    :accessor mp-critical-stress
    :type DOUBLE-FLOAT
    :initarg :critical-stress
    :initform 0d0)
   )
  (:documentation "A material point with a damage tensor"))
(defclass particle-thermal (particle)
  (
   (temperature
    :accessor mp-temperature
    :type DOUBLE-FLOAT
    :initarg :temperature
    :initform 0d0)
   (heat-capacity
    :accessor mp-heat-capacity
    :type DOUBLE-FLOAT
    :initarg :heat-capacity
    :initform 0d0)
   )
  (:documentation "A material point with a thermal properties"))
(defclass particle-fracture (particle-damage)
  (
   (strain-energy-density
    :accessor mp-strain-energy-density
    :type DOUBLE-FLOAT
    :initform 0d0)
   (strain-energy-density-local
    :accessor mp-strain-energy-density-local
    :type DOUBLE-FLOAT
    :initform 0d0)
   (fracture-toughness
    :accessor mp-fracture-toughness
    :type DOUBLE-FLOAT
    :initform 0d0
    :initarg :fracture-toughness)
   )
  (:documentation "A material point with fracture mechanics"))

(defclass particle-elastic-damage (particle-elastic particle-damage)
  ()
  (:documentation "A mp with damage influanced elastic model"))


(defclass particle-elastic-fracture (particle-elastic particle-fracture)
  ()
  (:documentation "A mp with fracture mechanics"))

(defclass particle-viscoplastic-damage (particle-viscoplastic particle-fracture)
  ()
  (:documentation "A mp with damage mechanics"))
(defclass particle-thermoviscoplastic-damage (particle-viscoplastic particle-damage particle-thermal)
  ()
  (:documentation "A mp with viscoplastic mechanics with variable thermal fields"))

(defclass particle-viscoelastic-fracture (particle-viscoelastic particle-fracture)
  ()
  (:documentation "A viscoelastic mp with fracture mechanics"))

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

(defgeneric constitutive-model (mp elastic-trial-strain dt)
    (:documentation "Compute new stress state given elastic strain")
    (:method (mp strain dt)
        (magicl:scale strain 0)))

(defmethod constitutive-model ((mp particle-elastic) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
    (with-slots ((E E)
                 (nu nu))
                mp
        (cl-mpm/constitutive:linear-elastic strain E nu)))

(defmethod constitutive-model ((mp particle-fluid) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((viscosity viscosity)
               (mass mass)
               (volume volume)
               (rest-density rest-density)
               (stiffness stiffness)
               (adiabatic-index adiabatic-index)
               )
      mp
    (let* ((density (/ mass volume))
           (pressure (* stiffness (expt (- (/ density rest-density) 1) adiabatic-index))))
      (cl-mpm/constitutive:newtonian-fluid strain pressure viscosity))))

(defmethod constitutive-model ((mp particle-viscoelastic) strain dt)
  "Function for modelling stress intergrated viscoelastic maxwell material"
  (with-slots ((E E)
               (nu nu)
               (strain-rate strain-rate) ;Note strain rate is actually strain increment through dt
               (vorticity vorticity)
               (stress stress))
      mp
    (cl-mpm/constitutive::maxwell strain-rate stress E nu dt vorticity)))

(defmethod constitutive-model ((mp particle-viscoplastic) strain dt)
  "Function for modelling stress intergrated viscoplastic norton-hoff material"
  (with-slots ((E E)
               (nu nu)
               (visc-factor visc-factor)
               (visc-power visc-power)
               (strain-rate strain-rate) ;Note strain rate is actually strain increment through dt
               (strain-plastic strain-plastic)
               (deformation-gradient deformation-gradient)
               (vorticity vorticity)
               (stress stress))
      mp
    ;; (let* ((linear-plastic-strain (cl-mpm/constitutive::norton-hoff-plastic-strain stress visc-factor visc-power dt)))
    ;;   (multiple-value-bind (l v) (magicl:eig (voigt-to-matrix linear-plastic-strain))
    ;;     (let (trial-plastic-strain ))
    ;;     )
    ;;   )
    ;; (setf strain-plastic (cl-mpm/constitutive::norton-hoff-plastic-strain (magicl:scale stress (magicl:det deformation-gradient)) visc-factor visc-power dt))
    ;; (setf strain
    ;;       (magicl:.- strain strain-plastic))
    ;; (cl-mpm/constitutive:linear-elastic strain E nu)
    ;(cl-mpm/constitutive::norton-hoff strain-rate stress E nu visc-factor visc-power dt vorticity)
    (cl-mpm/constitutive::glen-flow strain-rate stress (/ E (* 3d0 (- 1d0 nu nu))) visc-factor visc-power dt vorticity)
    ))
(defmethod constitutive-model ((mp particle-thermoviscoplastic-damage) strain dt)
  "Function for modelling stress intergrated viscoplastic norton-hoff material"
  (with-slots ((E E)
               (nu nu)
               (visc-factor visc-factor)
               (visc-power visc-power)
               (strain-rate strain-rate) ;Note strain rate is actually strain increment through dt
               (strain-plastic strain-plastic)
               (deformation-gradient deformation-gradient)
               (vorticity vorticity)
               (temperature temperature)
               (stress stress))
      mp
    (cl-mpm/constitutive::glen-flow strain-rate stress
                                    (/ E (* 3d0 (- 1d0 nu nu)))
                                    (* (expt (- temperature) 2) visc-factor)
                                    visc-power dt vorticity)
    ))

(defgeneric post-stress-step (mesh mp dt)
  (:documentation "This step gets called after full stress state resolved and allows for other processing"))
(defmethod post-stress-step (mesh mp dt)
  ())


