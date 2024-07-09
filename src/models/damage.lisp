(defpackage :cl-mpm/models/damage
  (:use
   :cl
   :cl-mpm/utils
   :cl-mpm/particle)
  (:export))
(in-package :cl-mpm/models/damage)

(defclass particle-damage (particle)
  (
   (log-damage
    :accessor mp-log-damage
    :type DOUBLE-FLOAT
    :initform 0d0)
   ;; (damage
   ;;  :accessor mp-damage
   ;;  :type DOUBLE-FLOAT
   ;;  :initarg :damage
   ;;  :initform 0d0)
   (local-damage
    :accessor mp-local-damage
    :type DOUBLE-FLOAT
    :initform 0d0)
   (local-damage-increment
    :accessor mp-local-damage-increment
    :type DOUBLE-FLOAT
    :initarg :damage-y
    :initform 0d0)
   (damage-increment
    :accessor mp-damage-increment
    :type DOUBLE-FLOAT
    :initform 0d0)
   (undamaged-stress
    :accessor mp-undamaged-stress
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils:voigt-zeros))
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
    :initform 0d0
    :initarg :damage-ybar
    )
   (damage-y-local
    :accessor mp-damage-y-local
    :type DOUBLE-FLOAT
    :initform 0d0
    :initarg :damage-y
    )
   (time-averaged-ybar
    :accessor mp-time-averaged-ybar
    :initform 1d0)
   (time-averaged-damage-inc
    :accessor mp-time-averaged-damage-inc
    :initform 1d0)
   (time-averaged-counter
    :accessor mp-time-averaged-counter
    :initform 0d0)
   (local-length
    :accessor mp-local-length
    :type DOUBLE-FLOAT
    :initarg :local-length
    :initform 1d0)
   (local-length-damaged
    :accessor mp-local-length-damaged
    :type DOUBLE-FLOAT
    :initarg :local-length-damaged
    :initform 1d0)
   (true-local-length
    :accessor mp-true-local-length
    :type DOUBLE-FLOAT
    :initarg :local-length-t
    :initform 1d0)
   (damage-position
    :accessor mp-damage-position
    ;:type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform nil
    ;:initform (cl-mpm/utils::vector-zeros)
    )
   (damage-tensor
    :accessor mp-damage-tensor
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils::matrix-zeros)
    )
   (damage-ybar-tensor
    :accessor mp-damage-ybar-tensor
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils::matrix-zeros))
   (damage-model
    :accessor mp-damage-model
    :initarg :damage-model
    )
   (damage-domain-update-rate
    :accessor mp-damage-domain-update-rate
    :initarg :damage-domain-rate
    :initform 0d0))
  (:documentation "A material point with a damage tensor"))

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

(defclass particle-creep-damage (particle-elastic-damage)
  ()
  (:documentation "A mp with isotropic damage - reduced for use in creep test (higher performance)"))
(defclass particle-elastic-fracture (particle-elastic particle-fracture)
  ()
  (:documentation "A mp with fracture mechanics"))

(defclass particle-viscoelastic-damage (particle-viscoelastic particle-damage)
  ()
  (:documentation "A visco-elastic material point with damage"))

(defclass particle-viscoplastic-damage (particle-viscoplastic particle-damage)
  ()
  (:documentation "A mp with damage mechanics"))

;; (defmethod constitutive-model :before ((mp particle-damage) strain dt)
;;   "For rate-based equations we want them to see the undamaged stress"
;;   (with-slots ((stress stress)
;;                (stress-u undamaged-stress))
;;       mp
;;     ;; (setf stress (magicl:scale stress-u 1d0))
;;     ))

(defmethod constitutive-model ((mp particle-elastic-damage) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (stress stress)
               (stress-undamaged undamaged-stress)
               (strain-rate strain-rate)
               (D stretch-tensor)
               (velocity-rate velocity-rate)
               (vorticity vorticity)
               (def deformation-gradient)
               (damage damage)
               (pressure pressure)
               ;; (datum pressure-datum)
               ;; (rho pressure-head)
               (pos position)
               (calc-pressure pressure-func)
               )
      mp
    (declare (double-float pressure damage)
             (function calc-pressure))
    (cl-mpm/constitutive::linear-elastic-mat strain de stress-undamaged)
    (cl-mpm/utils::voigt-copy-into stress-undamaged stress)
    (when (> damage 0d0)
      (cl-mpm/fastmath::fast-scale! stress (- 1d0 (* (- 1d0 1d-9) damage))))
    stress))

(defmethod constitutive-model ((mp particle-creep-damage) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (stress stress)
               (stress-undamaged undamaged-stress)
               (strain-rate strain-rate)
               (def deformation-gradient)
               (damage damage)
               )
      mp
    (declare (double-float damage))
    (setf stress-undamaged (cl-mpm/constitutive::linear-elastic-mat strain de))
    (setf stress (magicl:scale stress-undamaged (- 1d0 damage)))
    stress
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
               (pressure pressure)
               (pos position)
               (eng-strain-rate eng-strain-rate)
               (time-averaged-visc time-averaged-visc)
               (p p-modulus)

               (calc-pressure pressure-func)

               ;; (stress undamaged-stress)
               ;; (stress-damaged stress)
               )
      mp
    (declare (double-float E visc-factor visc-power))
    (let* (;; (eng-strain-rate (cl-mpm/fastmath::fast-.* (magicl:map (lambda (x) (* x (exp x))) strain) velocity-rate
           ;;                                        (cl-mpm/utils:voigt-from-list '(1d0 1d0 1d0 0.5d0 0.5d0 0.5d0))))
           (viscosity (cl-mpm/constitutive::glen-viscosity-strain eng-strain-rate visc-factor visc-power))
           ;; (viscosity 1d-20)
           (visc-u viscosity)
           (viscosity (* viscosity (max 1d-10 (expt (- 1d0 damage) 2))))
           )
      ;; (setf true-visc viscosity)
      ;; (setf p
      ;;       (/ (/ E (* (+ 1 nu) (- 1 nu)))
      ;;          (expt (max 1d-2 (- 1d0 damage)) 2)
      ;;          ;(expt (/ 1d13 viscosity) 1)
      ;;          ))
      ;; (setf stress (magicl:scale stress-undamaged (- 1d0 damage)))

      (incf time-averaged-visc viscosity)
      ;; (setf stress-u
      ;;       (cl-mpm/fastmath::fast-.+
      ;;        stress-u
      ;;        (objectify-stress-logspin
      ;;         (if (> viscosity 0d0)
      ;;             ;; (cl-mpm/constitutive::maxwell-damage strain-rate stress E nu de viscosity dt damage)
      ;;             ;; (cl-mpm/constitutive:maxwell-exp-v strain-rate stress E nu de visc-u dt))
      ;;             (cl-mpm/constitutive::maxwell strain-rate stress E nu de viscosity dt)
      ;;             (cl-mpm/constitutive::linear-elastic-mat strain-rate de))
      ;;         stress-u
      ;;         def
      ;;         vorticity
      ;;         D)))

      ;; (magicl:.+
      ;;  stress-u
      ;;  (cl-mpm/constitutive:maxwell strain-rate stress E nu de viscosity dt)
      ;;  stress-u)
      (setf stress-u
            (cl-mpm/constitutive:maxwell-exp-v
             strain-rate
             stress-u
             E nu de
             viscosity dt))
      ;; (when (> (abs(magicl:tref stress-u 3 0)) 1d-10)
      ;;   (break)
      ;;   )
      ;; (setf (magicl:tref stress-u 3 0) 0d0
      ;;       (magicl:tref stress-u 4 0) 0d0
      ;;       )


      ;; (setf stress-u (cl-mpm/constitutive:maxwell strain-rate stress E nu de viscosity dt))
      ;; (setf stress
      ;;       (cl-mpm/constitutive:maxwell-exp-v strain-rate stress E nu de viscosity dt))
      ;; (setf stress-u (cl-mpm/constitutive::linear-elastic-mat strain de))
      (setf stress (magicl:scale stress-u 1d0))
      (when (> damage 0.0d0)
        (let ((j (magicl:det def)))
          (multiple-value-bind (l v) (cl-mpm/utils::eig
                                      (magicl:scale! (voight-to-matrix stress) (/ 1d0 j)))

            (let* ((tp (funcall calc-pressure pos))
                  (driving-pressure (* tp (expt (min 0.90d0 damage) (magicl:det def))))
                  ;; (driving-pressure 0d0) ;;(* pressure (min 1.00d0 damage)))
                  (degredation (expt (- 1d0 damage) 2d0)))

              (loop for i from 0 to 2
                    do (let* ((sii (nth i l))
                              (esii (- sii driving-pressure)))
                         (when (> esii 0d0)
                           ;;Tensile damage -> unbounded
                           (setf (nth i l) (* esii (max 1d-8 degredation)))
                           (setf (nth i l) (+ (nth i l) driving-pressure))
                           )
                         (when (< esii 0d0)
                           ;;Bounded compressive damage
                           (setf (nth i l) (* esii (max 1d-1 degredation)))
                           (setf (nth i l) (+ (nth i l) driving-pressure))
                           )
                         ;; (setf (nth i l) (* sii (max 0d0 (- 1d0 damage))))
                         ))
              (setf stress (magicl:scale! (matrix-to-voight (magicl:@ v
                                                                      (magicl:from-diag l :type 'double-float)
                                                                      (magicl:transpose v))) j))
              ))))
      ;; (when nil;(> damage 0.0d0)
      ;;   (let ((j 1d0))
      ;;     (multiple-value-bind (l v) (cl-mpm/utils::eig
      ;;                                 (magicl:scale! (voight-to-matrix stress) (/ 1d0 j)))
      ;;       (let ((driving-pressure (* pressure (min 1.00d0 damage)))
      ;;             (degredation (expt (- 1d0 damage) 2d0)))
      ;;         (loop for i from 0 to 2
      ;;               do (let* ((sii (nth i l))
      ;;                         (esii (- sii driving-pressure)))
      ;;                    (when (> esii 0d0)
      ;;                      ;;Tensile damage -> unbounded
      ;;                      (setf (nth i l) (* esii (max 0d-8 degredation)))
      ;;                      (setf (nth i l) (+ (nth i l) driving-pressure))
      ;;                      )
      ;;                    (when (< esii 0d0)
      ;;                      ;;Bounded compressive damage
      ;;                      (setf (nth i l) (* esii (max 1d-2 degredation)))
      ;;                      (setf (nth i l) (+ (nth i l) driving-pressure))
      ;;                      )
      ;;                    ;; (setf (nth i l) (* sii (max 0d0 (- 1d0 damage))))
      ;;                    ))
      ;;         (setf stress (magicl:scale! (matrix-to-voight (magicl:@ v
      ;;                                                                 (magicl:from-diag l :type 'double-float)
      ;;                                                                 (magicl:transpose v))) j))
      ;;         ))))

      stress
      ;; (cl-mpm/constitutive::elasto-glen strain-rate stress E nu de viscosity dt)
      )))


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
    ;; (cl-mpm/fastmath::fast-.+ stress (cl-mpm/constitutive:linear-elastic strain-rate E nu) (objectify-stress-jaumann stress vorticity))
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
