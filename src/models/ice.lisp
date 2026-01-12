;; (defpackage :cl-mpm/models/ice
;;   (:use :cl
;;    :cl-mpm/utils)
;;   (:export
;;    )
;;   )
;; (in-package :cl-mpm/models/ice)

(in-package :cl-mpm/particle)

(defclass particle-glen (particle-elastic)
  ((visc-factor
    :accessor mp-visc-factor
    :initarg :visc-factor)
   (visc-power
    :accessor mp-visc-power
    :initarg :visc-power)
   (true-visc
    :accessor mp-true-visc
    :initform 0d0)
   )
  (:documentation "A glen flow law material point"))

(defclass particle-m-ice (particle-viscoplastic)
  ((plastic-stress
    :accessor mp-plastic-stress
    :initarg :plastic-stress)
   (visc-plastic
    :accessor mp-visc-plastic
    :initform 0d0)
   (visc-glen
    :accessor mp-visc-glen
    :initform 0d0)
   )
  (:documentation "A visco-plastic material point"))

(defclass particle-visco-elasto-plastic-damage (particle-viscoplastic particle-chalk-delayed)
  ((ductility-mode-2
    :accessor mp-ductility-mode-2
    :initarg :ductility-mode-2
    :initform 1d0)
   (damage-mode-2
    :accessor mp-damage-mode-2
    :initform 0d0)
   (enable-viscosity
    :accessor mp-enable-viscosity
    :initform t
    :initarg :enable-viscosity)
   (visc-water-damping
    :accessor mp-visc-water-damping
    :initarg :water-damping
    :initform 0d0))
  (:documentation "A ice mp with viscoplastic damage and viscoelastic relaxation "))

(defclass particle-ice-brittle (particle-elastic-damage  particle-mc)
  ((trial-elastic-strain
    :accessor mp-trial-strain
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils:voigt-zeros))
   (friction-angle
    :accessor mp-friction-angle
    :initarg :friction-angle
    :initform 30d0)
   (peerlings-damage
    :accessor mp-peerlings-damage
    :initform t
    :initarg :peerlings-damage)
   (shear-residual-ratio
    :accessor mp-shear-residual-ratio
    :initarg :g-res-ratio
    :initform 1d-9
    )
   (k-tensile-residual-ratio
    :accessor mp-k-tensile-residual-ratio
    :initarg :kt-res-ratio
    :initform 1d-9
    )
   (k-compressive-residual-ratio
    :accessor mp-k-compressive-residual-ratio
    :initarg :kc-res-ratio
    :initform 1d-9
    )
   (damage-tension
    :accessor mp-damage-tension
    :initarg :damage-tension
    :initform 0d0)
   (damage-compression
    :accessor mp-damage-compression
    :initarg :damage-compression
    :initform 0d0)
   (damage-shear
    :accessor mp-damage-shear
    :initarg :damage-shear
    :initform 0d0)
   )
  (:documentation "An ice damage model with a quasi-brittle damage model"))

(defclass particle-ice-delayed (particle-ice-brittle)
  (
   (delay-time
    :accessor mp-delay-time
    :initform 1d0
    :initarg :delay-time)
   (delay-exponent
    :accessor mp-delay-exponent
    :initform 1d0
    :initarg :delay-exponent))
  (:documentation "An ice damage model with time dependant damage evolution"))

(defclass particle-glen-damage (particle-glen particle-damage)
  ()
  (:documentation "A weakly compressible glen flow mp with damage mechanics"))



(defmethod constitutive-model ((mp particle-ice-brittle) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-accessors ((de mp-elastic-matrix)
                   (stress mp-stress)
                   (stress-u mp-undamaged-stress)
                   (strain mp-strain)
                   (trial-elastic-strain mp-trial-strain)
                   (damage mp-damage)
                   (damage-t mp-damage-tension)
                   (damage-c mp-damage-compression)
                   (damage-s mp-damage-shear)
                   (coheasion mp-c)
                   (ps-vm mp-strain-plastic-vm)
                   (ps-vm-inc mp-strain-plastic-vm-inc)
                   (ps-vm-1 mp-strain-plastic-vm-1)
                   (plastic-strain mp-strain-plastic)
                   (yield-func mp-yield-func)
                   (enable-plasticity mp-enable-plasticity)
                   (enable-damage mp-enable-damage)
                   (E mp-E)
                   (nu mp-nu)
                   (phi mp-phi)
                   (psi mp-psi)
                   (kc-r mp-k-compressive-residual-ratio)
                   (kt-r mp-k-tensile-residual-ratio)
                   (g-r mp-shear-residual-ratio)
                   (peerlings mp-peerlings-damage)
                   (L                cl-mpm/particle::mp-stretch-tensor)
                   (def mp-deformation-gradient)
                   (pressure mp-pressure)
                   (p-wave cl-mpm/particle::mp-p-modulus-0)
                   (p cl-mpm/particle::mp-pressure)
                   )
      mp
    (declare (magicl:matrix/double-float de stress stress-u strain plastic-strain)
             (double-float coheasion ps-vm-inc ps-vm yield-func E nu phi psi kc-r kt-r g-r damage))
    ;;Train elastic strain - plus trail kirchoff stress
    (setf stress-u (cl-mpm/constitutive::linear-elastic-mat strain de stress-u))
    (cl-mpm/utils:voigt-copy-into strain trial-elastic-strain)
    (setf p-wave (cl-mpm/particle::compute-p-modulus mp))
    (when enable-plasticity
        (let* (;; (epstr (voigt-copy strain))
               (K (/ e (* 3 (- 1d0 (* 2 nu)))))
               (pressure-strain (cl-mpm/utils::voigt-eye (/ pressure (* 3d0 K (/ 1d0 (magicl:det def)))))))
          ;; (cl-mpm/fastmaths:fast-.- strain pressure-strain strain)
          ;; (setf stress-u (cl-mpm/constitutive::linear-elastic-mat strain de stress-u))
          (multiple-value-bind (sig eps-e f inc pmod)
              (cl-mpm/ext::constitutive-mohr-coulomb stress-u
                                                     de
                                                     strain
                                                     E
                                                     nu
                                                     phi
                                                     psi
                                                     coheasion)
              ;; (cl-mpm/constitutive::vm-plastic stress-u
              ;;                                  de
              ;;                                  strain
              ;;                                  coheasion)
              ;; (cl-mpm/constitutive::mc-plastic stress-u
              ;;                                  de
              ;;                                  strain
              ;;                                  E
              ;;                                  nu
              ;;                                  phi
              ;;                                  psi
              ;;                                  coheasion)
            ;; (cl-mpm/fastmaths:fast-.+ eps-e pressure-strain eps-e)
            (setf sig (cl-mpm/constitutive::linear-elastic-mat eps-e de sig))
            (setf
             stress-u sig
             strain eps-e
             yield-func f
             ;; p-wave pmod
             )
            (let ((inc (expt (* 1/3 (max 0d0 (cl-mpm/utils::trace-voigt (cl-mpm/fastmaths:fast-.-
                  trial-elastic-strain strain)))) 1)))
              (setf ps-vm (+ ps-vm-1 inc))
              (setf ps-vm-inc inc))
            ;; (cl-mpm/fastmaths:fast-.+ plastic-strain
            ;;                           (cl-mpm/fastmaths:fast-.- trial-elastic-strain strain)
            ;;                           plastic-strain)
            ;; (cl-mpm/fastmaths:fast-.- trial-elastic-strain strain plastic-strain)
            ;; (setf inc (cl-mpm/utils:trace-voigt plastic-strain))
            ;; (let ()
            ;;   (incf ps-vm inc)
            ;;   (setf ps-vm-inc inc))
            )))
    (cl-mpm/utils:voigt-copy-into stress-u stress)
    stress))



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
               (eng-strain-rate eng-strain-rate)
               )
      mp
    (declare (double-float E visc-factor visc-power))
    (let* (;(viscosity (cl-mpm/constitutive::glen-viscosity-strain velocity-rate visc-factor visc-power))
           (viscosity (cl-mpm/constitutive::glen-viscosity-strain eng-strain-rate visc-factor visc-power)))
      (cl-mpm/constitutive::elasto-glen strain-rate stress E nu de viscosity dt strain)
      )))
(defmethod constitutive-model ((mp particle-m-ice) strain dt)
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
               (temp true-visc)
               (plastic-stress plastic-stress)
               (eng-strain-rate eng-strain-rate)
               (visc-plastic visc-plastic)
               (visc-glen visc-glen)
               (time-averaged-visc time-averaged-visc)
               (p p-modulus)
               )
      mp
    (declare (double-float E visc-factor visc-power))
    (flet ((inv (x) (/ 1d0 (the double-float x))))
      (let* (;(eng-strain-rate (cl-mpm/fastmaths::fast-.* (magicl:map (lambda (x) (* x (exp x))) strain) velocity-rate
                                        ;                            (cl-mpm/utils:stress-from-list '(1d0 1d0 0.5d0))))
             (estrain (cl-mpm/constitutive::effective-strain-rate eng-strain-rate))
             (viscosity-diff (/ 6.5d8 (* 2d0 (expt 10d-3 2))))
             (viscosity-glen (cl-mpm/constitutive::glen-viscosity-strain eng-strain-rate visc-factor visc-power))
             (viscosity-plastic (/ plastic-stress (+ 1d-40 (* 2d0 estrain))))
             (viscosity (+ 1d-40 (inv (+ ;(inv viscosity-diff)
                                         (inv viscosity-glen)
                                         (inv viscosity-plastic)))))
             ;; (viscosity viscosity-glen)
             )
        ;; stress
        ;; (print viscosity-glen)
        ;; ;; (print viscosity-plastic)
        ;; (setf p
        ;;       (* (/ E (* (+ 1 nu) (- 1 nu))) (/ 1d12 viscosity)))
        (setf temp viscosity)
        (setf visc-plastic viscosity-plastic
              visc-glen viscosity-glen)
        (incf time-averaged-visc viscosity)
        ;; (setf stress-u
        ;;       (cl-mpm/constitutive:maxwell-exp-v strain-rate stress E nu de visc-u dt))
        ;; (setf stress
        ;;       (cl-mpm/constitutive:maxwell-exp-v strain-rate stress E nu de viscosity dt))

        (cl-mpm/fastmaths::fast-.+
         stress
         (objectify-stress-logspin
          (if (> viscosity 0d0)
              (cl-mpm/constitutive::maxwell-exp-inc strain-rate stress E nu de viscosity dt)
              (cl-mpm/constitutive::linear-elastic-mat strain-rate de))
          stress
          def
          vorticity
          D
          ))
        ))
    ))

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
          ;(viscosity (* viscosity (- 1 (* damage (- 1 0.1)))))
           (viscosity (* viscosity (max 1d-3 (expt (- 1d0 damage) visc-power))))
          )
      ;; (cl-mpm/constitutive::elasto-glen-damage strain-rate stress E nu de viscosity dt strain damage)
      (cl-mpm/constitutive::elasto-glen strain-rate stress E nu de viscosity dt strain)
      )
    ))

(defmethod constitutive-model ((mp particle-visco-elasto-plastic-damage) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-accessors ((de mp-elastic-matrix)
                   (stress mp-stress)
                   (stress-u mp-undamaged-stress)
                   (strain mp-strain)
                   (strain-rate mp-strain-rate)
                   (damage mp-damage)
                   (pressure mp-pressure)
                   (def mp-deformation-gradient)
                   (pos mp-position)
                   (calc-pressure mp-pressure-func)
                   (coheasion mp-c)
                   (ps-vm mp-strain-plastic-vm)
                   (plastic-strain mp-strain-plastic)
                   (yield-func mp-yield-func)
                   (enable-plasticity mp-enable-plasticity)
                   (enable-viscosity mp-enable-viscosity)
                   (enable-damage mp-enable-damage)
                   (D mp-stretch-tensor)
                   (vorticity mp-vorticity)
                   (initiation-stress mp-initiation-stress)
                   (ductility-I mp-ductility)
                   (ductility-II mp-ductility-mode-2)
                   (E mp-E)
                   (nu mp-nu)
                   (phi mp-phi)
                   (psi mp-psi)
                   (kc-r mp-k-compressive-residual-ratio)
                   (kt-r mp-k-tensile-residual-ratio)
                   (g-r mp-shear-residual-ratio)

                   (eng-strain-rate mp-eng-strain-rate)
                   (visc-factor mp-visc-factor)
                   (visc-power mp-visc-power)
                   (true-visc mp-true-visc)
                   )
      mp
    (declare (function calc-pressure)
             (double-float damage pressure))
    ;;Train elastic strain - plus trail kirchoff stress

    ;; (let* ((p (cl-mpm/utils:trace-voigt strain))
    ;;        (q (cl-mpm/utils:deviatoric-voigt strain))
    ;;        (tau 1d0)
    ;;        (rho (exp (- (/ dt tau)))))
    ;;   (setf strain
    ;;         (magicl:.+
    ;;          (cl-mpm/constitutive::voight-eye (/ p 3d0))
    ;;          (magicl:scale! q rho))))

    ;;Somebody somewhere is destructivly changing stress-u which is very sad
    (setf stress-u (cl-mpm/constitutive::linear-elastic-mat strain de stress-u))
    ;;"Somebody" is me
    (when enable-viscosity
        (progn
          ;;Approximate visc over loadstep by \sigma_{n+1}
          (let ((viscosity (cl-mpm/constitutive::glen-viscosity stress-u visc-factor visc-power)))
            ;(setf true-visc viscosity)
            (setf stress-u (cl-mpm/constitutive:maxwell-exp-v
                            strain-rate
                            stress-u
                            E
                            nu
                            de
                            viscosity
                            dt
                            ))
            (setf strain (magicl:linear-solve de stress-u)))))

    ;; (let ((pressure (* pressure (expt damage 1) (magicl:det def))))
    ;;   (cl-mpm/fastmaths::fast-.+ stress-u
    ;;                             (cl-mpm/utils::voigt-eye (/ pressure 3d0)))
    ;;   (setf strain (magicl:linear-solve de stress-u)))
    (when enable-plasticity
      (multiple-value-bind (sig eps-e f)
          (cl-mpm/constitutive::mc-plastic
           stress-u
           de
           strain
           E
           nu
           phi
           psi
           coheasion)
        (setf stress-u
              sig
              plastic-strain (magicl:.- strain eps-e)
              yield-func f)
        (setf strain eps-e))
      ;; (multiple-value-bind (sig eps-e f)
      ;;     (cl-mpm/ext::constitutive-drucker-prager
      ;;      ;; stress-u
      ;;      strain
      ;;      de
      ;;      E
      ;;      nu
      ;;      phi
      ;;      psi
      ;;      coheasion)
      ;;   (setf stress-u
      ;;         sig
      ;;         plastic-strain (cl-mpm/fastmaths:fast-.+ plastic-strain (magicl:.- strain eps-e) plastic-strain)
      ;;         yield-func f)
      ;;   (setf strain eps-e))
      (incf ps-vm
            (multiple-value-bind (l v)
                (cl-mpm/utils:eig (cl-mpm/utils:voigt-to-matrix (cl-mpm/particle::mp-strain-plastic mp)))
              (destructuring-bind (s1 s2 s3) l
                (sqrt
                 (/ (+ (expt (- s1 s2) 2d0)
                       (expt (- s2 s3) 2d0)
                       (expt (- s3 s1) 2d0)
                       ) 2d0))))))

    (setf stress (magicl:scale stress-u 1d0))

    (when (and enable-damage
               (> damage 0.0d0))
      (let* ((d-p (* pressure damage 0d0 (magicl:det def))))
        (declare (double-float damage))
        (let* ((p (+ (/ (cl-mpm/constitutive::voight-trace stress) 3d0) d-p))
               (s (cl-mpm/constitutive::deviatoric-voigt stress))
               (ex 1)
               )
          ;; (incf p (* -1d0 pressure damage))
          (setf
           p
           (+ d-p
              (if (> p 0d0)
                  (* (expt (- 1d0 (* (- 1d0 kt-r) damage)) ex) p)
                  (* (expt (- 1d0 (* (- 1d0 kc-r) damage)) ex) p))))
          (cl-mpm/fastmaths:fast-.+
           (cl-mpm/constitutive::voight-eye p)
           (magicl:scale! s (expt (- 1d0 (* (- 1d0 g-r) damage)) ex))
           stress))))

    (let ((pressure (* pressure (expt damage 1)
                       (magicl:det def)
                       )))
      (cl-mpm/fastmaths::fast-.+ stress
                                (cl-mpm/utils::voigt-eye pressure)
                                stress))
    stress
    ))

(in-package :cl-mpm/damage)
(defmethod set-mp-damage ((mp cl-mpm/particle::particle-ice-brittle) d)
  (let ((k (cl-mpm/damage::find-k-damage-mp mp d)))
    (setf (cl-mpm/particle::mp-history-stress mp) k
          (cl-mpm/particle::mp-history-stress-n mp) k)
    (compute-damage mp)))

(defun damage-response-linear-ice (stress E Gf length init-stress ductility)
  (declare (double-float stress E Gf length init-stress ductility))
  "Function that controls how damage evolves with principal stresses"
  (let* ((ft init-stress)
         (e0 (/ ft E))
         (ef (* e0 ductility))
         (k (/ stress E)))
    (if (> k e0)
        (min 1d0
             (/ (max 0d0 (- k e0))
                (- ef e0)))
        0d0)))

(defmethod update-damage ((mp cl-mpm/particle::particle-visco-elasto-plastic-damage) dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (E cl-mpm/particle::mp-e)
                     (Gf cl-mpm/particle::mp-Gf)
                     (damage-inc cl-mpm/particle::mp-damage-increment)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (pressure cl-mpm/particle::mp-pressure)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     ;; (length-t cl-mpm/particle::mp-true-local-length)
                     (length cl-mpm/particle::mp-local-length)
                     (k cl-mpm/particle::mp-history-stress)
                     (tau cl-mpm/particle::mp-delay-time)
                     (tau-exp cl-mpm/particle::mp-delay-exponent)
                     (ductility cl-mpm/particle::mp-ductility)
                     (nu cl-mpm/particle::mp-nu)
                     ) mp
      (declare (double-float damage damage-inc critical-damage k ybar tau dt))
        (progn
          ;;Damage increment holds the delocalised driving factor
          (setf damage-inc 0d0)
          ;; (setf k (max k ybar))

          (when (> ybar init-stress)
            (let ((a tau-exp)
                  (k0 init-stress))
              (incf k
                    (the double-float
                         (*
                          dt
                          (/
                           (* k0
                              (expt
                               (/ (the double-float (max 0d0 (- ybar k)))
                                  k0)
                               a))
                           tau))))))

          (let ((new-damage
                  (max
                   damage
                   (damage-response-exponential k E init-stress ductility)
                   ;; (damage-response-linear k E Gf (/ length (the double-float (sqrt 7d0))) init-stress ductility)
                   )))
            (declare (double-float new-damage))
            (setf damage-inc (- new-damage damage))
            )

          (when (>= damage 1d0)
            (setf damage-inc 0d0)
            (setf ybar 0d0))
          (incf (the double-float(cl-mpm/particle::mp-time-averaged-damage-inc mp)) damage-inc)
          (incf (the double-float(cl-mpm/particle::mp-time-averaged-ybar mp)) ybar)
          (incf (the double-float(cl-mpm/particle::mp-time-averaged-counter mp)))
          ;;Transform to log damage
          (incf damage damage-inc)
          ;;Transform to linear damage
          (setf damage (max 0d0 (min 1d0 damage)))
          (when (> damage critical-damage)
            (setf damage 1d0)
            (setf damage-inc 0d0))
          )
  (values)))



(in-package :cl-mpm/damage)
(defmethod damage-model-calculate-y ((mp cl-mpm/particle::particle-ice-brittle) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
                     (strain cl-mpm/particle::mp-strain)
                     (trial-strain cl-mpm/particle::mp-trial-strain)
                     (plastic-strain cl-mpm/particle::mp-strain-plastic)
                     (ps-vm cl-mpm/particle::mp-strain-plastic-vm)
                     (damage cl-mpm/particle:mp-damage)
                     (pressure cl-mpm/particle::mp-pressure)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (angle cl-mpm/particle::mp-friction-angle)
                     (de cl-mpm/particle::mp-elastic-matrix)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     (E cl-mpm/particle::mp-e))
        mp
      (progn
        (setf damage-increment
              ;; (max 0d0
              (+
               (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile
                stress
                (* angle (/ pi 180d0)))))
        (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
        ;; (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
        ))))


(defun deriv-partial (k y k0 tau n)
  (if (> y k0)
      (/
       (* k0
          (expt
           (/ (the double-float (max 0d0 (- y k)))
              k0) n))
       tau)
      0d0))

(defmethod compute-damage ((mp cl-mpm/particle::particle-ice-brittle))
  (with-accessors ((k cl-mpm/particle::mp-history-stress)
                   (k-n cl-mpm/particle::mp-history-stress-n)
                   (E cl-mpm/particle::mp-e)
                   (damage cl-mpm/particle:mp-damage)
                   (init-stress cl-mpm/particle::mp-initiation-stress)
                   (ductility cl-mpm/particle::mp-ductility)
                   (damage-tension cl-mpm/particle::mp-damage-tension)
                   (damage-shear cl-mpm/particle::mp-damage-shear)
                   (damage-compression cl-mpm/particle::mp-damage-compression)
                   (peerlings cl-mpm/particle::mp-peerlings-damage)
                   (kc-r cl-mpm/particle::mp-k-compressive-residual-ratio)
                   (kt-r cl-mpm/particle::mp-k-tensile-residual-ratio)
                   (g-r cl-mpm/particle::mp-shear-residual-ratio))
      mp
    (declare (double-float kt-r kc-r g-r damage k k-n))
    ;; Directly compute the damage from K
    (let ()
      (setf damage (damage-response-exponential k E init-stress ductility))
      (if peerlings
          (setf
           damage-tension (damage-response-exponential-peerlings-residual k E init-stress ductility kt-r)
           damage-shear (damage-response-exponential-peerlings-residual k E init-stress ductility g-r)
           damage-compression (damage-response-exponential-peerlings-residual k E init-stress ductility kc-r))
          (setf
           damage-tension (* kt-r damage)
           damage-compression (* kc-r damage)
           damage-shear (* g-r damage))))))

(defmethod update-damage ((mp cl-mpm/particle::particle-ice-brittle) dt)
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
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (pressure cl-mpm/particle::mp-pressure)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     ;; (length-t cl-mpm/particle::mp-true-local-length)
                     (length cl-mpm/particle::mp-local-length)
                     (k cl-mpm/particle::mp-history-stress)
                     (k-n cl-mpm/particle::mp-history-stress-n)
                     (tau cl-mpm/particle::mp-delay-time)
                     (tau-exp cl-mpm/particle::mp-delay-exponent)
                     (ductility cl-mpm/particle::mp-ductility)
                     (nu cl-mpm/particle::mp-nu)
                     (p cl-mpm/particle::mp-p-modulus)
                     (kc-r cl-mpm/particle::mp-k-compressive-residual-ratio)
                     (kt-r cl-mpm/particle::mp-k-tensile-residual-ratio)
                     (g-r cl-mpm/particle::mp-shear-residual-ratio)
                     (damage-tension cl-mpm/particle::mp-damage-tension)
                     (damage-shear cl-mpm/particle::mp-damage-shear)
                     (damage-compression cl-mpm/particle::mp-damage-compression)
                     (peerlings cl-mpm/particle::mp-peerlings-damage)
                     )
        mp
      (declare (double-float damage damage-inc damage-n critical-damage k ybar tau dt ybar-prev init-stress k-n ybar))
      (when t;(<= damage 1d0)
        ;;Damage increment holds the delocalised driving factor
        ;; (let ((k0 init-stress))
        ;;   (when (or (>= ybar-prev k0)
        ;;             (>= ybar k0))
        ;;     (setf k
        ;;           (max
        ;;            k-n
        ;;            ybar))))
        (setf k
              (max
               k-n
               ybar-prev
               ybar))
        (compute-damage mp)
        (setf damage-inc (- damage damage-n))

        (incf (the double-float (cl-mpm/particle::mp-time-averaged-damage-inc mp)) (* damage-inc dt))
        (incf (the double-float (cl-mpm/particle::mp-time-averaged-ybar mp)) ybar)
        (incf (the double-float (cl-mpm/particle::mp-time-averaged-counter mp)))
        ;;Transform to log damage
        (setf damage (max 0d0 (min 1d0 damage))))
      (values))))

(defmethod update-damage ((mp cl-mpm/particle::particle-ice-delayed) dt)
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
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (pressure cl-mpm/particle::mp-pressure)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     ;; (length-t cl-mpm/particle::mp-true-local-length)
                     (length cl-mpm/particle::mp-local-length)
                     (k cl-mpm/particle::mp-history-stress)
                     (k-n cl-mpm/particle::mp-history-stress-n)
                     (tau cl-mpm/particle::mp-delay-time)
                     (tau-exp cl-mpm/particle::mp-delay-exponent)
                     (ductility cl-mpm/particle::mp-ductility)
                     (nu cl-mpm/particle::mp-nu)
                     (p cl-mpm/particle::mp-p-modulus)
                     (kc-r cl-mpm/particle::mp-k-compressive-residual-ratio)
                     (kt-r cl-mpm/particle::mp-k-tensile-residual-ratio)
                     (g-r cl-mpm/particle::mp-shear-residual-ratio)
                     (damage-tension cl-mpm/particle::mp-damage-tension)
                     (damage-shear cl-mpm/particle::mp-damage-shear)
                     (damage-compression cl-mpm/particle::mp-damage-compression)
                     (peerlings cl-mpm/particle::mp-peerlings-damage)
                     (ps-vm cl-mpm/particle::mp-strain-plastic-vm)
                     )
        mp
      (declare (double-float damage damage-inc damage-n critical-damage k ybar tau dt ybar-prev init-stress k-n ybar))
      (when t;(<= damage 1d0)
        ;;Damage increment holds the delocalised driving factor
        (let ((k0 init-stress)
              (ps-y (sqrt (* E (expt ps-vm 2)))))
          (when (or (>= ybar-prev k0)
                    (>= ybar k0))
            (setf k
                  (max
                   k-n
                   (+
                    ps-y
                    (cl-mpm/damage::huen-integration
                     k-n
                     ybar-prev
                     ybar
                     k0
                     tau
                     tau-exp
                     dt
                     ))))))
        (compute-damage mp)
        (setf damage-inc (- damage damage-n))

        (incf (the double-float (cl-mpm/particle::mp-time-averaged-damage-inc mp)) (* damage-inc dt))
        (incf (the double-float (cl-mpm/particle::mp-time-averaged-ybar mp)) ybar)
        (incf (the double-float (cl-mpm/particle::mp-time-averaged-counter mp)))
        ;;Transform to log damage
        (setf damage (max 0d0 (min 1d0 damage))))
      (values))))

(defun apply-vol-pressure-degredation (mp dt pressure)
  (with-accessors ((damage        cl-mpm/particle::mp-damage)
                   (damage-t      cl-mpm/particle::mp-damage-tension)
                   (damage-c      cl-mpm/particle::mp-damage-compression)
                   (damage-s      cl-mpm/particle::mp-damage-shear)
                   (stress        cl-mpm/particle::mp-stress)
                   (p-mod         cl-mpm/particle::mp-p-modulus-0)
                   (E cl-mpm/particle::mp-e)
                   (nu cl-mpm/particle::mp-nu)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (enable-damage cl-mpm/particle::mp-enable-damage))
      mp
    (declare (double-float damage damage-t damage-c damage-s))
    (when (and
           enable-damage
           (> damage 0.0d0))
      (let* ((exponent 1)
            (p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
            (pind (- p pressure))
            ;; (pind p)
            (s (cl-mpm/constitutive::deviatoric-voigt stress)))
        (declare (double-float damage-t damage-c damage-s))
        (setf p
              (if (> pind 0d0)
                  (* (- 1d0 (expt damage-t exponent)) p)
                  (* (- 1d0 (expt damage-c exponent)) p)))
        ;; (pprint pressure)
        ;; (break)
        (setf stress
              (cl-mpm/fastmaths:fast-.+
               (cl-mpm/constitutive::voight-eye (- p pressure))
               (cl-mpm/fastmaths:fast-scale! s (- 1d0 (expt damage-s exponent)))
               stress))
        (let ((K (/ e (* 3 (- 1d0 (* 2 nu)))))
              (G (/ e (* 2 (+ 1d0 nu)))))
          (setf K
                (if (> pind 0d0)
                    (* (- 1d0 (expt damage-t exponent)) K)
                    (* (- 1d0 (expt damage-c exponent)) K)))
          (setf G (* G (- 1d0 (expt damage-s exponent))))

          (when (> (cl-mpm/particle::mp-yield-func mp) 0d0)
            ;; (setf G 0d0)
            (setf K (* K (cos (cl-mpm/particle::mp-phi mp))))
            (setf G (* G (sin (cl-mpm/particle::mp-phi mp))))
            )
          (setf p-mod (+ K (* 4/3 G)))
          )
        ))))

(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-ice-brittle) dt)
  (with-accessors ((p cl-mpm/particle::mp-pressure)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (stress cl-mpm/particle::mp-stress)
                   (damage cl-mpm/particle::mp-damage)
                   (enable-damage cl-mpm/particle::mp-enable-damage)
                   (p-mod cl-mpm/particle::mp-p-modulus)
                   )
      mp
    ;; (setf p-mod (* (expt (cl-mpm/fastmaths::det def) -2) (cl-mpm/particle::compute-p-modulus mp)))
    (when enable-damage
      ;; (apply-vol-degredation mp dt)
      ;; (apply-vol-pressure-degredation mp dt (* 1d0 (magicl:det def) (/ p 3) damage))
      ;; (apply-vol-pressure-degredation mp dt (* -1d0 (/ 1d0 (magicl:det def)) (/ p 1) damage))
      (apply-vol-pressure-degredation mp dt (* -1d0
                                               ;; (/ 1d0 (magicl:det def))
                                               (/ p 3) damage))
      ;; (setf stress (cl-mpm/constitutive::voight-eye p))
      )
      ;; (cl-mpm/fastmaths:fast-.+ stress
      ;;                           (cl-mpm/constitutive::voight-eye (* ;; (magicl:det def)
      ;;                                                               (/ p 3) damage)) stress)
    (cl-mpm/particle::update-log-p-wave mp)
    )
  
  )

