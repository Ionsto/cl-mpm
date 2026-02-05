;; (defpackage :cl-mpm/models/damage
;;   (:use
;;    :cl
;;    :cl-mpm/utils
;;    :cl-mpm/particle)
;;   (:export))
;; (in-package :cl-mpm/models/damage)

(in-package :cl-mpm/particle)



(defclass particle-damage-mpi (particle-damage-fundemental)
    ((volume
      :accessor mp-volume
      :type double-float
      :initarg :volume
      :initform 1d0)
     (position
      :accessor mp-position
      :type MAGICL:MATRIX/DOUBLE-FLOAT
      :initform (cl-mpm/utils:vector-zeros)
      :initarg :position)
     (damage
      :accessor mp-damage
      :type DOUBLE-FLOAT
      :initarg :damage
      :initform 0d0)
     (true-local-length
      :accessor mp-true-local-length
      :type DOUBLE-FLOAT
      :initarg :local-length
      :initform 1d0)
     (local-damage-increment
      :accessor mp-local-damage-increment
      :type DOUBLE-FLOAT
      :initarg :damage-y
      :initform 0d0)
     (damage-y-local
      :accessor mp-damage-y-local
      :type DOUBLE-FLOAT
      :initform 0d0
      :initarg :damage-y)
     (damage-position
      :accessor mp-damage-position
      :initform nil)))



(defclass particle-fracture (particle-damage)
  ((strain-energy-density
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
  ((ductility
    :accessor mp-ductility
    :initarg :ductility
    :initform 1d0)
   (history-stress
    :accessor mp-history-stress
    :initform 0d0)
   (residual-strength
    :initarg :residual-strength
    :accessor mp-residual-strength
    :initform 1d-6)
   (history-stress-n
    :accessor mp-history-stress-n
    :initform 0d0))
  (:documentation "A mp with damage influanced elastic model"))

(defmethod cl-mpm/particle::reset-loadstep-mp ((mp particle-elastic-damage))
  (with-accessors ((k    cl-mpm/particle::mp-history-stress)
                   (k-n    cl-mpm/particle::mp-history-stress-n)
                   (y cl-mpm/particle::mp-damage-y-local)
                   (ybar cl-mpm/particle::mp-damage-ybar)
                   (damage    cl-mpm/particle::mp-damage)
                   (damage-prev    cl-mpm/particle::mp-damage-prev-trial)
                   (damage-inc    cl-mpm/particle::mp-damage-increment)
                   (damage-n    cl-mpm/particle::mp-damage-n))
      mp
    ;; (pprint "Reset load")
    (setf k k-n
          damage damage-n
          damage-prev damage
          damage-inc 0d0
          y 0d0
          ybar 0d0)
    (cl-mpm/damage::compute-damage mp)
    (call-next-method)))

(defmethod new-loadstep-mp ((mp particle-elastic-damage))
  (with-accessors ((k    cl-mpm/particle::mp-history-stress)
                   (k-n    cl-mpm/particle::mp-history-stress-n)
                   (damage    cl-mpm/particle::mp-damage)
                   (damage-n    cl-mpm/particle::mp-damage-n))
      mp
    (setf k-n k
          damage-n damage)
    (call-next-method)))

(defclass particle-elastic-damage-delayed (particle-elastic-damage)
   ((delay-time
     :accessor mp-delay-time
     :initform 1d0
     :initarg :delay-time
     )
    (delay-exponent
     :accessor mp-delay-exponent
     :initform 1d0
     :initarg :delay-exponent
     ))
   (:documentation "A mp with isotropic damage - with delayed damage"))

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
               (def deformation-gradient)
               (damage damage))
      mp
    (declare (double-float damage))
    (setf (cl-mpm/particle::mp-p-modulus mp)
          (*
           (cl-mpm/particle::estimate-log-enhancement mp)
           (cl-mpm/particle::compute-p-modulus mp)))
    (cl-mpm/constitutive::linear-elastic-mat strain de stress-undamaged)
    (cl-mpm/utils::voigt-copy-into stress-undamaged stress)
    stress))

(defun apply-isotropic-degredation (mp)
  (with-accessors ((damage        cl-mpm/particle::mp-damage)
                   (stress        cl-mpm/particle::mp-stress)
                   (enable-damage cl-mpm/particle::mp-enable-damage)
                   (p-mod cl-mpm/particle::mp-p-modulus)
                   (e cl-mpm/particle::mp-e)
                   (nu cl-mpm/particle::mp-nu))
      mp
    (declare (double-float damage))
    (when (and
           enable-damage
           (> damage 0.0d0))
      (cl-mpm/fastmaths:fast-scale! stress (- 1d0 damage))
      (setf (cl-mpm/particle::mp-p-modulus mp)
            (*
             (- 1d0 damage)
             (cl-mpm/particle::estimate-log-enhancement mp)
             (cl-mpm/particle::compute-p-modulus mp)))
      )))

(defun apply-vol-degredation (mp)
  (with-accessors ((damage        cl-mpm/particle::mp-damage)
                   (stress        cl-mpm/particle::mp-stress)
                   (enable-damage cl-mpm/particle::mp-enable-damage)
                   (p-mod cl-mpm/particle::mp-p-modulus)
                   (e cl-mpm/particle::mp-e)
                   (nu cl-mpm/particle::mp-nu)
                   )
      mp
    (declare (double-float damage))
    (when (and
           enable-damage
           (> damage 0.0d0))
      ;; (setf p-mod (cl-mpm/particle::compute-p-modulus mp))
      (let ((p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
            (s (cl-mpm/constitutive::deviatoric-voigt stress)))
        (setf stress
              (cl-mpm/fastmaths:fast-.+
               (cl-mpm/constitutive::voight-eye p)
               (cl-mpm/fastmaths:fast-scale! s (- 1d0 damage))
               stress))
        ;; (let ((K (/ e (* 3 (- 1d0 (* 2 nu)))))
        ;;       (G (* (- 1d0 damage) (/ e (* 2 (+ 1d0 nu))))))
        ;;   (setf p-mod (+ K (* G (/ 4 3)))))
        

        ))))

(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-elastic-damage) dt)
  (apply-isotropic-degredation mp)
  ;; (apply-vol-degredation mp)
  )

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
               ;; (p p-modulus)

               (calc-pressure pressure-func)

               ;; (stress undamaged-stress)
               ;; (stress-damaged stress)
               )
      mp
    (declare (double-float E visc-factor visc-power))
    (let* (;; (eng-strain-rate (cl-mpm/fastmaths::fast-.* (magicl:map (lambda (x) (* x (exp x))) strain) velocity-rate
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
      ;;       (cl-mpm/fastmaths::fast-.+
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
    ;; (cl-mpm/fastmaths::fast-.+ stress (cl-mpm/constitutive:linear-elastic strain-rate E nu) (objectify-stress-jaumann stress vorticity))
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
      (cl-mpm/constitutive:newtonian-fluid strain pressure viscosity))))


(defmethod cl-mpm/damage::update-damage ((mp cl-mpm/particle::particle-elastic-damage-delayed) dt)
  (when (cl-mpm/particle::mp-enable-damage mp)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (damage-n cl-mpm/particle::mp-damage-n)
                     (E cl-mpm/particle::mp-e)
                     (damage-inc cl-mpm/particle::mp-damage-increment)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (ybar-prev cl-mpm/particle::mp-damage-ybar-prev)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (residual-strength cl-mpm/particle::mp-residual-strength)
                     (length cl-mpm/particle::mp-local-length)
                     (k cl-mpm/particle::mp-history-stress)
                     (k-n cl-mpm/particle::mp-history-stress-n)
                     (tau cl-mpm/particle::mp-delay-time)
                     (tau-exp cl-mpm/particle::mp-delay-exponent)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (ductility cl-mpm/particle::mp-ductility))
        mp
      (declare (double-float damage damage-inc k ybar tau dt))
      (when t;(<= damage 1d0)
        ;;Damage increment holds the delocalised driving factor
        (setf (cl-mpm/particle::mp-damage-prev-trial mp) (cl-mpm/particle::mp-damage mp))
        (let ((k0 init-stress))
          (when (or
                 (>= ybar-prev k0)
                 (>= ybar k0))
            (setf k
                  (cl-mpm/damage::huen-integration
                   k-n
                   ybar-prev
                   ybar
                   k0
                   tau
                   tau-exp
                   dt))))
        (let ((new-damage
                (max
                 damage-n
                 (cl-mpm/damage::damage-response-exponential-peerlings-residual k E init-stress ductility (- 1d0 residual-strength)))))
          (declare (double-float new-damage))
          (setf damage new-damage)
          (setf damage-inc (- new-damage damage-n)))
        (when (>= damage 1d0)
          (setf damage-inc 0d0))
        (incf (the double-float (cl-mpm/particle::mp-time-averaged-damage-inc mp)) (* damage-inc dt))
        (incf (the double-float (cl-mpm/particle::mp-time-averaged-ybar mp)) ybar)
        (incf (the double-float (cl-mpm/particle::mp-time-averaged-counter mp)))
        ;;Transform to log damage
        ;; (incf damage damage-inc)
        ;;Transform to linear damage
        (setf damage (max 0d0 (min 1d0 damage)))
        (when (> damage critical-damage)
          (setf damage 1d0)
          (setf damage-inc 0d0)))
      (values))))

(defmethod cl-mpm/damage::compute-damage ((mp cl-mpm/particle::particle-elastic-damage))
  (with-accessors ((damage cl-mpm/particle:mp-damage)
                   (damage-n cl-mpm/particle::mp-damage-n)
                   (E cl-mpm/particle::mp-e)
                   (ybar cl-mpm/particle::mp-damage-ybar)
                   (ybar-prev cl-mpm/particle::mp-damage-ybar-prev)
                   (init-stress cl-mpm/particle::mp-initiation-stress)
                   (k cl-mpm/particle::mp-history-stress)
                   (ductility cl-mpm/particle::mp-ductility))
      mp
    (declare (double-float damage k ybar))
    (setf damage
          (cl-mpm/damage::damage-response-exponential-peerlings-residual
           k
           E init-stress ductility (- 1d0 1d-9)))))

(defmethod cl-mpm/damage::update-damage ((mp cl-mpm/particle::particle-elastic-damage) dt)
  (when (cl-mpm/particle::mp-enable-damage mp)
    (with-accessors ((damage cl-mpm/particle:mp-damage)
                     (damage-n cl-mpm/particle::mp-damage-n)
                     (E cl-mpm/particle::mp-e)
                     (damage-inc cl-mpm/particle::mp-damage-increment)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (ybar-prev cl-mpm/particle::mp-damage-ybar-prev)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (length cl-mpm/particle::mp-local-length)
                     (k cl-mpm/particle::mp-history-stress)
                     (k-n cl-mpm/particle::mp-history-stress-n)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (ductility cl-mpm/particle::mp-ductility))
        mp
      (declare (double-float damage damage-inc k ybar dt))
      (setf (cl-mpm/particle::mp-damage-prev-trial mp) (cl-mpm/particle::mp-damage mp))
      (when t;(<= damage 1d0)
        ;;Damage increment holds the delocalised driving factor
        (setf k
              (max
               k-n
               ybar-prev
               ybar))
        (cl-mpm/damage::compute-damage mp)
        (setf damage-inc (- damage damage-n))
        ;; (when (> damage critical-damage)
        ;;   (setf damage 1d0)
        ;;   (setf damage-inc 0d0))
        )
      (values))))

(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-elastic-damage) dt)
  (with-accessors ((strain cl-mpm/particle::mp-strain)
                   (E cl-mpm/particle::mp-e)
                   (y cl-mpm/particle::mp-damage-y-local)
                   (de cl-mpm/particle::mp-elastic-matrix))
      mp
    (setf y (cl-mpm/damage::tensile-energy-norm strain E de))))

