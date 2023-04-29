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
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(defstruct node-cache
  node
  weight
  grads)

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
   (index
     :accessor mp-index
     :type integer
     :initform 0
     :initarg :index)
   (volume
     :accessor mp-volume
     :type double-float
     :initarg :volume
     :initform 1d0)
   (volume-0
    :accessor mp-volume-0
    :type double-float
    :initarg :volume
    :initform 1d0)
   (size-0
    :accessor mp-domain-size-0
    :type magicl:matrix/double-float
    :initarg :size-0
    :initform (magicl:zeros '(2 1)))
   (size
     :accessor mp-domain-size
     :type magicl:matrix/double-float
     :initarg :size
     :initform (magicl:ones '(2 1)))
   (position
     :accessor mp-position
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :initarg :position)
   (velocity
     :accessor mp-velocity
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :initarg :velocity
     :initform (magicl:zeros '(2 1)))
   (acceleration
    :accessor mp-acceleration
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (magicl:zeros '(2 1)))
   (stress
     :accessor mp-stress
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :initarg :stress
     :initform (magicl:zeros '(3 1)))
   (stress-kirchoff
    :accessor mp-stress-kirchoff
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
   (stretch-tensor
    :accessor mp-stretch-tensor
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (magicl:zeros '(2 2)))
   (strain-rate-tensor
    :accessor mp-strain-rate-tensor
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (magicl:zeros '(2 2)))
   (strain-rate
     :accessor mp-strain-rate
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :accessor mp-strain-rate
     :initarg :strain-rate
     :initform (magicl:zeros '(3 1)))
   (velocity-rate
    :accessor mp-velocity-rate
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :accessor mp-velocity-rate
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
   (pressure
    :type double-float
    :accessor mp-pressure
    :initform 0d0
    )
   (gravity-axis
    :type magicl:matrix/double-float
    :accessor mp-gravity-axis
    :initform (magicl:from-list '(0d0 -1d0) '(2 1)))
   (body-force
     :accessor mp-body-force
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :initarg :body-force
     :initform (magicl:zeros '(2 1)))
   (displacement
    :accessor mp-displacement
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (magicl:zeros '(2 1)))
   (cached-nodes
    :accessor mp-cached-nodes
    :initform (make-array 8 :fill-pointer 0 :element-type 'node-cache)
    )
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
   (elastic-matrix
    :accessor mp-elastic-matrix
    :type magicl:matrix/double-float)
   )
  (:documentation "A linear-elastic material point"))
(defun update-elastic-matrix (particle)
  (with-accessors ((de mp-elastic-matrix)
                   (E  mp-E)
                   (nu mp-nu))
      particle
      (setf de (cl-mpm/constitutive::linear-elastic-matrix E nu))))
(defmethod (setf mp-E) :after (value (p particle-elastic))
  (update-elastic-matrix p))
(defmethod (setf mp-nu) :after (value (p particle-elastic))
  (update-elastic-matrix p))
(defmethod initialize-instance :after ((p particle-elastic) &key)
  (update-elastic-matrix p))

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

(defclass particle-viscoelastic (particle-elastic)
  ((viscosity
    :accessor mp-viscosity
    :initarg :viscosity)
   )
  (:documentation "A visco-elastic material point"))

(defclass particle-glen (particle-elastic)
  ((visc-factor
    :accessor mp-visc-factor
    :initarg :visc-factor)
   (visc-power
    :accessor mp-visc-power
    :initarg :visc-power)
   )
  (:documentation "A glen flow law material point"))

(defclass particle-viscoplastic (particle-elastic)
  ((visc-factor
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
   (local-damage
    :accessor mp-local-damage
    :type DOUBLE-FLOAT
    :initform 0d0)
   (local-damage-increment
    :accessor mp-local-damage-increment
    :type DOUBLE-FLOAT
    :initform 0d0)
   (damage-increment
    :accessor mp-damage-increment
    :type DOUBLE-FLOAT
    :initform 0d0)
   (undamaged-stress
    :accessor mp-undamaged-stress
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (magicl:zeros '(3 1)))
   (initiation-stress
    :accessor mp-initiation-stress
    :type DOUBLE-FLOAT
    :initarg :initiation-stress
    :initform 0d0)
   (critical-stress
    :accessor mp-critical-stress
    :type DOUBLE-FLOAT
    :initarg :critical-stress
    :initform 0d0)
   (damage-rate
    :accessor mp-damage-rate
    :type DOUBLE-FLOAT
    :initarg :damage-rate
    :initform 0d0)
   (critical-damage
    :accessor mp-critical-damage
    :type DOUBLE-FLOAT
    :initarg :critical-damage
    :initform 1d0)
   (damage-ybar
    :accessor mp-damage-ybar
    :type DOUBLE-FLOAT
    :initform 0d0)
   (local-length
    :accessor mp-local-length
    :type DOUBLE-FLOAT
    :initarg :local-length
    :initform 1d0)

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
   (thermal-conductivity
    :accessor mp-thermal-conductivity
    :type DOUBLE-FLOAT
    :initarg :thermal-conductivity
    :initform 1d0)
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

(defclass particle-viscoelastic-damage (particle-viscoelastic particle-damage)
  ()
  (:documentation "A visco-elastic material point with damage"))

(defclass particle-viscoplastic-damage (particle-viscoplastic particle-damage)
  ()
  (:documentation "A mp with damage mechanics"))

(defclass particle-glen-damage (particle-glen particle-damage)
  ()
  (:documentation "A weakly compressible glen flow mp with damage mechanics"))

(defclass particle-thermoelastic-damage (particle-elastic particle-damage particle-thermal)
  ()
  (:documentation "A mp with elastic mechanics with variable thermal fields"))
(defclass particle-thermoviscoplastic-damage (particle-viscoplastic particle-damage particle-thermal)
  ()
  (:documentation "A mp with viscoplastic mechanics with variable thermal fields"))

(defclass particle-thermofluid-damage (particle-fluid particle-damage particle-thermal)
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
      (let ((mp (apply #'make-instance constructor
                      :nD nD
                      :volume (coerce volume 'double-float)
                      :mass (coerce mass 'double-float)
                      :position position
                      args)))
        (progn
          (setf (mp-domain-size-0 mp) (magicl:scale (mp-domain-size-0 mp) 1d0))
          mp)))))
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

(defmethod constitutive-model :before ((mp particle-damage) strain dt)
  "For rate-based equations we want them to see the undamaged stress"
  (with-slots ((stress stress)
               (stress-u undamaged-stress))
      mp
    ;; (setf stress (magicl:scale stress-u 1d0))
    ))


(defmethod constitutive-model ((mp particle-elastic) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (stress stress)
               (strain-rate strain-rate)
               (velocity-rate velocity-rate)
               (vorticity vorticity)
               (def deformation-gradient)
               )
      mp
    ;; Kirchoff stress rate?
    ;; (magicl:.+ stress
    ;; (magicl:.+
    ;;  stress
    ;;  (cl-mpm/constitutive::linear-elastic-mat strain-rate de))
    ;; (magicl:.+
    ;;  (cl-mpm/constitutive:linear-elastic strain-rate E nu)
    ;;  (objectify-stress-kirchoff-truesdale stress vorticity))
    ;;Jaumann rate equation
    ;; (magicl:.+
    ;;  stress
    ;;  (objectify-stress-jaumann
    ;;   (cl-mpm/constitutive:linear-elastic strain-rate E nu)
    ;;   stress
    ;;   vorticity))
    ;; Truesdale rate
    ;; (magicl:.+
    ;;  stress
    ;;  (objectify-stress-kirchoff-truesdale
    ;;   (cl-mpm/constitutive:linear-elastic strain-rate E nu)
    ;;   stress
    ;;   strain-rate))
    (cl-mpm/constitutive::linear-elastic-mat strain de)
    ;; (cl-mpm/constitutive:linear-elastic strain E nu)
    ))

(defmethod constitutive-model ((mp particle-elastic-damage) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (stress stress)
               (stress-undamaged undamaged-stress)
               (strain-rate strain-rate)
               (velocity-rate velocity-rate)
               (vorticity vorticity)
               (def deformation-gradient)
               (damage damage)
               (pressure pressure)
               )
      mp
    ;; Non-objective stress intergration
    (setf stress-undamaged
           ;; (cl-mpm/constitutive::linear-elastic-mat strain de)
          (magicl:.+
           stress-undamaged
           (cl-mpm/constitutive::linear-elastic-mat strain-rate de)))
    (setf stress (magicl:scale stress-undamaged 1d0))

    (when (> damage 0.0d0)
      (multiple-value-bind (l v) (magicl:eig
                                  (voight-to-matrix stress))
        (loop for i from 0 to 1
              do (let* ((sii (nth i l))
                        (esii (- sii
                                 (* pressure 1)))
                        )
                   (when (> esii 0d0)
                     (setf (nth i l)
                           (+
                            (* esii (- 1d0 damage))
                            (* 1 pressure)
                            ))
                     ))
        (setf stress (matrix-to-voight (magicl:@ v
                                                 (magicl:from-diag l :type 'double-float)
                                                 (magicl:transpose v))))
        )))
    stress
    ))

(defclass particle-elastic-inc (particle-elastic)
  ()
  (:documentation "A linear-elastic material point"))
(defmethod constitutive-model ((mp particle-elastic-inc) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (stress stress)
               (strain-rate strain-rate)
               (velocity-rate velocity-rate)
               (vorticity vorticity)
               (def deformation-gradient)
               )
      mp
    ;; Non-objective stress intergration
    (magicl:.+
     stress
     (cl-mpm/constitutive::linear-elastic-mat strain-rate de))
    ))
(defclass particle-elastic-jaumann (particle-elastic)
  ()
  (:documentation "A linear-elastic material point"))
(defmethod constitutive-model ((mp particle-elastic-jaumann) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (stress stress)
               (strain-rate strain-rate)
               (velocity-rate velocity-rate)
               (vorticity vorticity)
               (def deformation-gradient)
               )
      mp
    ;; Kirchoff stress rate?
    ;; (magicl:.+ stress
    ;; (magicl:.+
    ;;  stress
    ;;  (cl-mpm/constitutive::linear-elastic-mat strain-rate de))
    ;; (magicl:.+
    ;;  (cl-mpm/constitutive:linear-elastic strain-rate E nu)
    ;;  (objectify-stress-kirchoff-truesdale stress vorticity))
    ;;Jaumann rate equation
    (magicl:.+
     stress
     (objectify-stress-jaumann
      (cl-mpm/constitutive::linear-elastic-mat strain-rate de)
      stress
      vorticity))
    ;; Truesdale rate
    ;; (magicl:.+
    ;;  stress
    ;;  (objectify-stress-kirchoff-truesdale
    ;;   (cl-mpm/constitutive:linear-elastic strain-rate E nu)
    ;;   stress
    ;;   strain-rate))
    ))
(defclass particle-elastic-truesdale (particle-elastic)
  ()
  (:documentation "A linear-elastic material point"))
(defmethod constitutive-model ((mp particle-elastic-truesdale) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (stress stress)
               (strain-rate strain-rate)
               (velocity-rate velocity-rate)
               (vorticity vorticity)
               (def deformation-gradient)
               )
      mp
    ;; Kirchoff stress rate?
    ;; (magicl:.+ stress
    ;; (magicl:.+
    ;;  stress
    ;;  (cl-mpm/constitutive::linear-elastic-mat strain-rate de))
    ;; (magicl:.+
    ;;  (cl-mpm/constitutive:linear-elastic strain-rate E nu)
    ;;  (objectify-stress-kirchoff-truesdale stress vorticity))
    ;;Jaumann rate equation
    ;; (magicl:.+
    ;;  stress
    ;;  (objectify-stress-jaumann
    ;;   (cl-mpm/constitutive:linear-elastic strain-rate E nu)
    ;;   stress
    ;;   vorticity))
    ;; Truesdale rate
    (magicl:.+
     stress
     (objectify-stress-kirchoff-truesdale
      (cl-mpm/constitutive::linear-elastic-mat strain-rate de)
      stress
      strain-rate))
    ))
(defclass particle-elastic-logspin (particle-elastic)
  ()
  (:documentation "A linear-elastic material point"))
(defmethod constitutive-model ((mp particle-elastic-logspin) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (stress stress)
               (strain-rate strain-rate)
               (D stretch-tensor)
               (velocity-rate velocity-rate)
               (vorticity vorticity)
               (def deformation-gradient)
               )
      mp
    (let (
          ;; (strain-rate (magicl:scale strain-rate (/ 1d0 dt)))
          ;; (vorticity (magicl:scale vorticity (/ 1d0 dt)))
          )
      (magicl:.+
       stress
       (objectify-stress-logspin
        (cl-mpm/constitutive::linear-elastic-mat strain-rate de)
        stress
        def
        vorticity
        D)))
    ))

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
               (viscosity viscosity)
               (strain-rate strain-rate) ;Note strain rate is actually strain increment through dt
               (vorticity vorticity)
               (stress stress))
      mp
    (cl-mpm/constitutive:maxwell strain-rate stress E nu viscosity dt vorticity)))

(defmethod constitutive-model ((mp particle-glen) strain dt)
  "Function for modeling stress intergrated viscoplastic norton-hoff material"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (visc-factor visc-factor)
               (visc-power visc-power)
               (strain-rate strain-rate) ;Note strain rate is actually strain increment through dt
               (velocity-rate velocity-rate) ;Note strain rate is actually strain increment through dt
               (stress stress)
               (strain strain)
               (stretch stretch-tensor)
               )
      mp
    (declare (double-float E visc-factor visc-power))
    (let* ((viscosity (cl-mpm/constitutive::glen-viscosity-strain velocity-rate visc-factor visc-power))
           ;(viscosity (* viscosity 0d0))
           ;; (viscosity 1d6)
           )
      (cl-mpm/constitutive::elasto-glen strain-rate stress E nu de viscosity dt strain)
      )))
(defmethod constitutive-model ((mp particle-viscoplastic) strain dt)
  "Function for modeling stress intergrated viscoplastic norton-hoff material"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (visc-factor visc-factor)
               (visc-power visc-power)
               (strain-rate strain-rate) ;Note strain rate is actually strain increment through dt
               (velocity-rate velocity-rate) ;Note strain rate is actually strain increment through dt
               (strain-plastic strain-plastic)
               (def deformation-gradient)
               (vorticity vorticity)
               (D stretch-tensor)
               (stress stress)
               )
      mp
    (declare (double-float E visc-factor visc-power))
    (let (
          (viscosity (cl-mpm/constitutive::glen-viscosity-strain velocity-rate visc-factor visc-power))
          ;(viscosity (cl-mpm/constitutive::glen-viscosity stress (expt visc-factor (- visc-power)) visc-power))
          )
      ;; stress
      (magicl:.+
       stress
       (objectify-stress-logspin
        (if (> viscosity 0d0)
            (cl-mpm/constitutive::maxwell strain-rate stress E nu de viscosity dt)
            (cl-mpm/constitutive::linear-elastic-mat strain-rate de))
        stress
        def
        vorticity
        D))
      ;; (cl-mpm/constitutive::elasto-glen strain-rate stress E nu de viscosity dt)
      ;; (magicl:.+
      ;;  stress
      ;;   (if (> viscosity 0d0)
      ;;       ;; (cl-mpm/constitutive::maxwell-exp-v strain-rate stress E nu de viscosity dt)
      ;;       (cl-mpm/constitutive::maxwell strain-rate stress E nu de viscosity dt)
      ;;       (cl-mpm/constitutive::linear-elastic-mat strain-rate de)))
      )
    ))

(defmethod constitutive-model ((mp particle-viscoplastic-damage) strain dt)
  "Function for modeling stress intergrated viscoplastic norton-hoff material"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (visc-factor visc-factor)
               (visc-power visc-power)
               (strain-rate strain-rate) ;Note strain rate is actually strain increment through dt
               (velocity-rate velocity-rate) ;Note strain rate is actually strain increment through dt
               (strain-plastic strain-plastic)
               (damage damage)
               (critical-damage critical-damage)
               (def deformation-gradient)
               (D stretch-tensor)
               (vorticity vorticity)
               (stress stress)
               (stress-u undamaged-stress)
               ;; (stress undamaged-stress)
               ;; (stress-damaged stress)
               )
      mp
    (declare (double-float E visc-factor visc-power))
    (let* ((viscosity (cl-mpm/constitutive::glen-viscosity-strain velocity-rate visc-factor visc-power))
           ;; (viscosity (* viscosity (- 1 (* damage (- 1 0.1)))))
           ;; (E-damage (* E (- 1 damage)))
           )
      (setf stress-u
            (magicl:.+
             stress-u
             (objectify-stress-logspin
              (if (> viscosity 0d0)
                  (cl-mpm/constitutive::maxwell-damage strain-rate stress E nu de viscosity dt damage)
                  (cl-mpm/constitutive::linear-elastic-mat strain-rate de))
              stress-u
              def
              vorticity
              D)))
      (setf stress (magicl:scale stress-u 1d0))
      stress
      ;; (cl-mpm/constitutive::elasto-glen strain-rate stress E nu de viscosity dt)
      )))

(defmethod constitutive-model ((mp particle-glen-damage) strain dt)
  "Function for modeling stress intergrated viscoplastic norton-hoff material"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (visc-factor visc-factor)
               (visc-power visc-power)
               (strain-rate strain-rate) ;Note strain rate is actually strain increment through dt
               (velocity-rate velocity-rate) ;Note strain rate is actually strain increment through dt
               (stress stress)
               (strain strain)
               (damage damage)
               )
      mp
    (declare (double-float E visc-factor visc-power))
    (let* ((viscosity (cl-mpm/constitutive::glen-viscosity-strain velocity-rate visc-factor visc-power))
          (viscosity (* viscosity (- 1 (* damage (- 1 0.1)))))
          ;; (E (* E damage))
          )
      (cl-mpm/constitutive::elasto-glen strain-rate stress E nu de viscosity dt strain)
      )
    ))

(defmethod constitutive-model ((mp particle-viscoelastic-damage) strain dt)
  "Function for modelling stress intergrated viscoelastic maxwell material"
  (with-slots ((E E)
               (nu nu)
               (viscosity viscosity)
               (strain-rate strain-rate) ;Note strain rate is actually strain increment through dt
               (vorticity vorticity)
               (def deformation-gradient)
               (stress undamaged-stress)
               ;; (stress stress)
               )
      mp
    ;; (magicl:.+ stress (cl-mpm/constitutive:linear-elastic strain-rate E nu) (objectify-stress-jaumann stress vorticity))
   ;(cl-mpm/constitutive:maxwell strain-rate (magicl:scale stress (magicl:det def)) E nu viscosity dt vorticity)
    (cl-mpm/constitutive:maxwell-exp strain-rate stress E nu viscosity dt vorticity)
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
               (stress-u undamaged-stress))
      mp
    (let ((viscosity (cl-mpm/constitutive::glen-viscosity stress-u visc-factor visc-power)))
            (cl-mpm/constitutive::maxwell strain-rate stress-u E nu viscosity dt vorticity))
    ))

(defmethod constitutive-model ((mp particle-thermofluid-damage) strain dt)
  "Function for modelling stress intergrated viscoplastic norton-hoff material"
  (with-slots ((viscosity viscosity)
               (mass mass)
               (volume volume)
               (rest-density rest-density)
               (stiffness stiffness)
               (temp temperature)
               (adiabatic-index adiabatic-index)
               )
      mp
    (let* ((density (/ mass volume))
           (effective-rest-density rest-density)
           (pressure (* stiffness
                        (expt
                         (- (/ density effective-rest-density)
                            1) adiabatic-index))))
      (cl-mpm/constitutive:newtonian-fluid strain pressure viscosity)))
  )

(defgeneric post-stress-step (mesh mp dt)
  (:documentation "This step gets called after full stress state resolved and allows for other processing"))
(defmethod post-stress-step (mesh mp dt)
  ())

(defun assemble-vorticity-matrix (vorticity)
  (let ((dx (magicl:tref vorticity 0 0))
        (dy (magicl:tref vorticity 1 0))
        (dxdy (magicl:tref vorticity 2 0))
        )
    (magicl:from-list (list dx dxdy (- dxdy) dy) '(2 2))))

(defun objectify-stress (mp)
  (cl-mpm/particle:mp-stress mp))

(defun objectify-stress-jaumann (stress-inc stress vorticity)
  (magicl:.-
   stress-inc
   (matrix-to-voight
    (magicl::.- (magicl:@ (voight-to-matrix stress) (assemble-vorticity-matrix vorticity))
                (magicl:@ (assemble-vorticity-matrix vorticity) (voight-to-matrix stress)))))
             )

(defun objectify-stress-kirchoff-truesdale (stress-inc stress velocity-rate)
  (magicl:.-
   stress-inc
   (matrix-to-voight
    (magicl::.- (magicl:@ (voight-to-matrix stress) (voight-to-matrix velocity-rate))
                (magicl:@ (voight-to-matrix velocity-rate) (voight-to-matrix stress))
                ))))

(declaim (notinline objectify-stress-logspin))
(defun objectify-stress-logspin (stress-inc stress def vorticity D)
    (let ((b (magicl:@ def (magicl:transpose def)))
          (omega (assemble-vorticity-matrix vorticity))
          (D (magicl:scale (magicl:.+ D (magicl:transpose D)) 0.5d0))
          )
        (multiple-value-bind (l v) (magicl:eig b)
          (loop for i from 0 to 1
                do (loop for j from 0 to 1
                         do
                            ;;For all the pairs of eigenvalues
                            (when (not (= i j))
                              (let ((l_i (nth i l))
                                    (l_j (nth j l))
                                    (v_i (magicl:column v i))
                                    (v_j (magicl:column v j))
                                    )
                                ;; (when (< l_i 0d0)
                                ;;   (magicl:scale! v_i -1d0)
                                ;;   (setf l_i (abs l_i))
                                ;;   )
                                ;; (when (< l_j 0d0)
                                ;;   (magicl:scale! v_j -1d0)
                                ;;   (setf l_j (abs l_j))
                                ;;   )
                                ;;When the eigenvalues are distinct
                                (when (and
                                       (> (abs (- l_i l_j)) 1d-6)
                                       ;;When they are nonzero
                                       ;; (> l_i 0d0)
                                       ;; (> l_j 0d0)
                                       )
                                  ;; When we have pairs of unique nonzero eigenvalues
                                  (setf omega (magicl:.+
                                               omega
                                               (magicl:scale
                                                (magicl:@
                                                 (magicl:@
                                                  v_i
                                                  (magicl:transpose v_i))
                                                 D
                                                 (magicl:@
                                                  v_j
                                                  (magicl:transpose v_j)))
                                                (+
                                                 (/ (+ 1 (/ l_i l_j)) (- 1 (/ l_i l_j)))
                                                 (/ 2d0 (log (/ l_i l_j))))
                                                )))
                                  ))))))
      ;; (print omega)
      ;; (break)
      ;; (break)
      (magicl:.-
       stress-inc
       (matrix-to-voight
        (magicl::.- (magicl:@ (voight-to-matrix stress) omega)
                    (magicl:@ omega (voight-to-matrix stress)))))
      ))
;; (defun objectify-stress-logspin (mp)
;;   (with-accessors ((stress cl-mpm/particle:mp-stress)
;;                    (def cl-mpm/particle:mp-deformation-gradient)
;;                    (vorticity cl-mpm/particle:mp-vorticity)
;;                    (strain-rate cl-mpm/particle::mp-strain-rate)
;;                    )
;;       mp
;;     (let ((b (magicl:@ def (magicl:transpose def)))
;;           (omega (assemble-vorticity-matrix vorticity))
;;           (D (voight-to-matrix strain-rate)))
;;         (multiple-value-bind (l v) (magicl:eig b)
;;           (loop for i from 0 to 2
;;                 do (loop for j from 0 to 2
;;                          do
;;                             (when (not (= i j))
;;                               (let ((l_i (nth i l))
;;                                     (l_j (nth j l))
;;                                     )
;;                                 (when (< (abs (- l_i l_j)) 1d-6)
;;                                   ;; When we have unique eigenvalues
;;                                   (setf omega (magicl:.+ omega
;;                                                          (magicl:scale
;;                                                           (magicl:@ (magicl:column v i)
;;                                                                     D
;;                                                                     (magicl:column v j)
;;                                                                     )
;;                                                           (+
;;                                                            (/ (+ l_i l_j) (- l_i l_j))
;;                                                            (/ 2 (- (log l_i) (log l_j))))
;;                                                           ))))))))))))

