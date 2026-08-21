(in-package :cl-mpm/particle)
(defclass particle-plastic-damage-frictional (particle-damage-frictional particle-mc)
  ((plastic-damage-evolution
    :accessor mp-plastic-damage-evolution
    :initform nil
    :initarg :plastic-damage-evolution)
   (plastic-strain-tensile
    :accessor mp-plastic-strain-tensile
    :initform nil
    :initarg :plastic-strain-tensile)
   (oversize-scale
    :accessor mp-oversize-scale
    :initarg :oversize
    :initform (- 1d0 1d-3))
   (residual-friction
    :accessor mp-residual-friction
    :initarg :residual-friction
    :initform 0d0)))

(defclass particle-plastic-damage-frictional-delayed (particle-plastic-damage-frictional)
  ((delay-time
    :accessor mp-delay-time
    :initform 1d0
    :initarg :delay-time)
   (delay-exponent
    :accessor mp-delay-exponent
    :initform 1d0
    :initarg :delay-exponent)))


(defmethod initialize-instance :after ((mp particle-plastic-damage-frictional) &key)
  (with-accessors ((ductility cl-mpm/particle::mp-ductility)
                   (angle cl-mpm/particle::mp-friction-angle)
                   (angle-r cl-mpm/particle::mp-residual-friction)
                   (rc mp-k-compressive-residual-ratio)
                   (init-stress mp-initiation-stress)
                   (oversize mp-oversize-scale))
      mp
    (let* ((c (cl-mpm/damage::mohr-coloumb-tensile-to-coheasion init-stress (rad-to-deg angle)))
           (rs (cl-mpm/damage::est-shear-from-angle (rad-to-deg angle) (rad-to-deg angle-r) rc))
           (residual-strength (mp-residual-strength mp))
           (oversize-ratio (cl-mpm/damage::compute-oversize-factor oversize ductility)))
      (setf
       (cl-mpm/particle::mp-phi mp) angle
       (mp-c mp) (* oversize-ratio c)
       (mp-shear-residual-ratio mp) rs)
      ;; (setf (mp-shear-residual-ratio mp) (min (mp-shear-residual-ratio mp) residual-strength)
      ;;       (mp-k-tensile-residual-ratio mp) (min residual-strength (mp-k-tensile-residual-ratio mp)))
      )))

(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-plastic-damage-frictional) dt)
  (cl-mpm/damage::apply-tcs-degredation mp))

(defmethod constitutive-model ((mp particle-plastic-damage-frictional) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-accessors ((E mp-e)
                   (nu mp-nu)
                   (de mp-elastic-matrix)
                   (stress mp-stress)
                   (stress-undamaged mp-undamaged-stress)
                   (def mp-deformation-gradient)
                   (enable-plasticity mp-enable-plasticity)
                   (phi mp-phi)
                   (psi mp-psi)
                   (coheasion mp-c)
                   (model cl-mpm/particle::mp-friction-model)
                   (p-wave cl-mpm/particle::mp-p-modulus-0)
                   (ps-vm mp-strain-plastic-vm)
                   (ps-vm-inc mp-strain-plastic-vm-inc)
                   (ps-vm-1 mp-strain-plastic-vm-1)
                   (yield-func mp-yield-func)
                   (damage mp-damage))
      mp
    (declare (double-float damage))
    (setf p-wave (cl-mpm/particle::compute-p-modulus mp))
    (cl-mpm/constitutive::linear-elastic-mat strain de stress-undamaged)
    (when enable-plasticity
      (let ((trial-elastic-strain (cl-mpm/utils::voigt-copy strain)))
        (multiple-value-bind (sig eps-e f inc pmod)
            (cl-mpm/ext::constitutive-mohr-coulomb
             stress-undamaged
             de
             strain
             E
             nu
             phi
             psi
             coheasion)
            ;; (cl-mpm/ext::constitutive-vm
            ;;  stress-undamaged
            ;;  strain
            ;;  de
            ;;  E
            ;;  nu
            ;;  coheasion
            ;;  ;; (cl-mpm/particle::mp-tangent-stiffness mp)
            ;;  )
          ;; (ecase model
          ;;   (:MC
          ;;    (cl-mpm/ext::constitutive-mohr-coulomb
          ;;     stress-undamaged
          ;;     de
          ;;     strain
          ;;     E
          ;;     nu
          ;;     phi
          ;;     psi
          ;;     coheasion
          ;;     ))
          ;;   (:DP
          ;;    (cl-mpm/ext::constitutive-drucker-prager
          ;;     stress-undamaged
          ;;     strain
          ;;     de
          ;;     E
          ;;     nu
          ;;     phi
          ;;     psi
          ;;     coheasion)))
          ;; (setf sig (cl-mpm/constitutive::linear-elastic-mat eps-e de sig))
          (setf
           stress-undamaged sig
           strain eps-e
           yield-func f
           p-wave (* 1.0d0 pmod)
           )
          (when (cl-mpm/particle::mp-plastic-strain-tensile mp)
            (setf inc (expt (* 1/3 (max 0d0
                                        (- (cl-mpm/utils::trace-voigt trial-elastic-strain)
                                           (cl-mpm/utils::trace-voigt strain)))) 1)))
          (setf ps-vm (+ ps-vm-1 inc))
          (setf ps-vm-inc inc)
          ;; (setf ps-vm (+ ps-vm-1 inc))
          ;; (setf ps-vm-inc inc)
          )))
    (cl-mpm/utils::voigt-copy-into stress-undamaged stress)
    stress))

(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-plastic-damage-frictional) dt)
  (with-accessors ((strain cl-mpm/particle::mp-strain)
                   (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                   (E cl-mpm/particle::mp-e)
                   (y cl-mpm/particle::mp-damage-y-local)
                   (angle cl-mpm/particle::mp-friction-angle)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (model cl-mpm/particle::mp-friction-model)
                   (pd-inc cl-mpm/particle::mp-plastic-damage-evolution)
                   (ps-vm cl-mpm/particle::mp-strain-plastic-vm)
                   (de cl-mpm/particle::mp-elastic-matrix)
                   )
      mp
    (let ((stress
            undamaged-stress
            ;; (cl-mpm/fastmaths:fast-scale-voigt undamaged-stress (/ 1d0 (cl-mpm/fastmaths:det-3x3 def)))
            )
          (ps-y (sqrt (* E (expt ps-vm 2))))
          )
      (setf
       y
       (+
        (if pd-inc ps-y 0d0)
        (ecase model
          (:SE (cl-mpm/damage::tensile-energy-norm strain e de))
          (:DP (cl-mpm/damage::drucker-prager-criterion stress angle))
          (:MC (cl-mpm/damage::criterion-mohr-coloumb-rankine-stress-tensile stress angle))))))))

;;A range of test cases that are used in the thesis



(defclass cl-mpm/particle::particle-fpd-isotropic (cl-mpm/particle::particle-plastic-damage-frictional)
  ())
(defclass cl-mpm/particle::particle-fpd-tcs (cl-mpm/particle::particle-plastic-damage-frictional)
  ())

(defclass cl-mpm/particle::particle-fpd-tcs-strain (cl-mpm/particle::particle-plastic-damage-frictional)
  ())

(defclass cl-mpm/particle::particle-fpd-spectral (cl-mpm/particle::particle-plastic-damage-frictional)
  ())

(defclass cl-mpm/particle::particle-fpd-spectral-strain (cl-mpm/particle::particle-plastic-damage-frictional)
  ())

(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-fpd-tcs) dt)
  (cl-mpm/damage::apply-tcs-degredation mp))

(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-fpd-tcs-strain) dt)
  (cl-mpm/damage::apply-tcs-strain-degredation mp))

(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-fpd-spectral) dt)
  ;; (setf (cl-mpm/particle::mp-damage mp) (cl-mpm/particle::mp-damage-tension mp))
  (cl-mpm/damage::apply-tensile-stress-degredation mp))

(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-fpd-spectral-strain) dt)
  ;; (setf (cl-mpm/particle::mp-damage mp) (cl-mpm/particle::mp-damage-tension mp))
  (cl-mpm/damage::apply-tensile-strain-degredation mp))

(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-fpd-isotropic) dt)
  ;; (setf (cl-mpm/particle::mp-damage mp) (cl-mpm/particle::mp-damage-tension mp))
  (cl-mpm/damage::apply-isotropic-degredation mp))

(defmethod update-damage ((mp cl-mpm/particle::particle-plastic-damage-frictional-delayed) dt)
  (when (cl-mpm/particle::mp-enable-damage mp)
    (with-accessors ((damage cl-mpm/particle:mp-damage)
                     (damage-n cl-mpm/particle::mp-damage-n)
                     (E cl-mpm/particle::mp-e)
                     (Gf cl-mpm/particle::mp-Gf)
                     (damage-inc cl-mpm/particle::mp-damage-increment)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (ybar-prev cl-mpm/particle::mp-damage-ybar-prev)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (length cl-mpm/particle::mp-local-length)
                     (k cl-mpm/particle::mp-history-stress)
                     (k-n cl-mpm/particle::mp-history-stress-n)
                     (tau cl-mpm/particle::mp-delay-time)
                     (tau-exp cl-mpm/particle::mp-delay-exponent)
                     (ductility cl-mpm/particle::mp-ductility))
        mp
      (declare (double-float damage damage-inc k ybar tau dt))
      (when t;(<= damage 1d0)
        ;;Damage increment holds the delocalised driving factor
        (break)
        (setf (cl-mpm/particle::mp-damage-prev-trial mp) (cl-mpm/particle::mp-damage mp))
        (setf damage-inc 0d0)
        (let ((a tau-exp)
              (k0 init-stress))
          (when t;; (or (>= ybar-prev k0)
                ;;     (>= ybar k0))
            (setf k
                  (max
                   k-n
                   (cl-mpm/damage::auto-refine-substepper
                    k-n
                    ybar-prev
                    ybar
                    dt
                    (lambda (k y0 y1 s-dt)
                      (cl-mpm/damage::huen-integration
                       k
                       y0
                       y1
                       k0
                       tau
                       tau-exp
                       s-dt))
                    :tol 1d-1)
                   ))))
        (compute-damage mp)
        (setf damage-inc (- damage damage-n))
        (setf (cl-mpm/particle::mp-time-averaged-damage-inc mp) (/ damage-inc dt))
        (incf (the double-float (cl-mpm/particle::mp-time-averaged-counter mp)))
        (setf damage (max 0d0 (min 1d0 damage))))
      (values))))
