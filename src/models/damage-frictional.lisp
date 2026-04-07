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

(defmethod initialize-instance :after ((mp particle-plastic-damage-frictional) &key)
  (with-accessors ((ductility cl-mpm/particle::mp-ductility)
                   (angle cl-mpm/particle::mp-friction-angle)
                   (angle-r cl-mpm/particle::mp-residual-friction)
                   (rc mp-shear-residual-ratio)
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
  ;; (setf (cl-mpm/particle::mp-damage mp) (cl-mpm/particle::mp-damage-tension mp))
  ;; (cl-mpm/damage::apply-isotropic-degredation mp)
  ;; (pprint "hello")
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
             coheasion
             )
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
