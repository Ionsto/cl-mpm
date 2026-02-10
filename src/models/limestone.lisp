(defpackage :cl-mpm/models/limestone
  (:use
   :cl
   :cl-mpm/utils
   :cl-mpm/particle)
  (:export))
(in-package :cl-mpm/models/limestone)
(in-package :cl-mpm/particle)
(defclass cl-mpm/particle::particle-limestone (cl-mpm/particle::particle-elastic-damage cl-mpm/particle::particle-mc)
  (
   (compression-ratio
    :accessor mp-compression-ratio
    :initarg :compression-ratio
    :initform 1d0)
   ;; (interal-length
   ;;  :accessor mp-internal-length
   ;;  :type DOUBLE-FLOAT
   ;;  :initarg :internal-length
   ;;  :initform 1d0)
   (enable-plasticity
    :accessor mp-enable-plasticity
    :initarg :enable-plasticity
    :initform nil))
  (:documentation "A concrete damage model"))

(defclass cl-mpm/particle::particle-limestone-delayed (particle-limestone)
  ((delay-time
    :accessor mp-delay-time
    :initform 1d0
    :initarg :delay-time
    ))
  (:documentation "A time dependant limestone elastic damage model"))

(defmethod cl-mpm/particle::constitutive-model ((mp particle-limestone) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-accessors
        ((phi mp-phi)
         (psi mp-psi)
         (de mp-elastic-matrix)
         (coheasion mp-c)
         (plastic-strain mp-strain-plastic)
         (yield-func mp-yield-func)
         (enable-plasticity mp-enable-plasticity)
         (ps-vm mp-strain-plastic-vm)
         (damage mp-damage)
         (stress mp-stress)
         (stress-undamaged mp-undamaged-stress)
         )
      mp
    (declare (double-float damage))
    ;; Non-objective stress intergration
    (setf stress-undamaged (cl-mpm/constitutive::linear-elastic-mat strain de stress-undamaged))
    (setf stress (cl-mpm/utils:voigt-copy-into stress-undamaged stress))
    ;; (when (> damage 0.0d0)
    ;;   (cl-mpm/fastmaths::fast-scale! stress (- 1d0 (* (- 1d0 1d-9) damage))))
    stress))
(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-limestone) dt)
  (cl-mpm/damage::apply-isotropic-degredation mp))

(defmethod cl-mpm/damage::update-damage ((mp cl-mpm/particle::particle-limestone) dt)
  (when t;(cl-mpm/particle::mp-enable-damage mp)
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
        (setf damage-inc (- damage damage-n)))
      (values))))

(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-limestone) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((strain cl-mpm/particle::mp-strain)
                     (stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (pressure cl-mpm/particle::mp-pressure)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     (de cl-mpm/particle::mp-elastic-matrix)
                     (E cl-mpm/particle::mp-e)
                     (nu cl-mpm/particle::mp-nu)
                     (k cl-mpm/particle::mp-compression-ratio)
                     ) mp
      (declare (double-float pressure damage E k nu))
        (progn
          (setf damage-increment (* E (cl-mpm/damage::modified-vm-criterion strain nu k)))
          ;; (setf damage-increment (cl-mpm/damage::modified-vm-stress stress nu k))
          )
      (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
      )))

(defmethod cl-mpm/damage::update-damage ((mp cl-mpm/particle::particle-limestone-delayed) dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (strain cl-mpm/particle::mp-strain)
                     (damage cl-mpm/particle:mp-damage)
                     (E cl-mpm/particle::mp-e)
                     (Gf cl-mpm/particle::mp-Gf)
                     (damage-inc cl-mpm/particle::mp-damage-increment)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (pressure cl-mpm/particle::mp-pressure)
                     (volume cl-mpm/particle::mp-volume)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     ;(length cl-mpm/particle::mp-internal-length)
                     (length cl-mpm/particle::mp-local-length)
                     (k cl-mpm/particle::mp-history-stress)
                     (eng-inc cl-mpm/particle::mp-dissipated-energy-inc)
                     (eng-int cl-mpm/particle::mp-dissipated-energy)
                     (p cl-mpm/particle::mp-p-modulus)
                     (nu cl-mpm/particle::mp-nu)
                     (ductility cl-mpm/particle::mp-ductility)
                     (tau cl-mpm/particle::mp-delay-time)
                     ) mp
      (declare (double-float damage damage-inc critical-damage init-stress ductility k E tau ybar damage-inc dt))
        (progn
          ;;Damage increment holds the delocalised driving factor
          (setf (cl-mpm/particle::mp-damage-prev-trial mp) (cl-mpm/particle::mp-damage mp))

          (setf ybar damage-inc)

          (when (> ybar init-stress)
            (incf k (the double-float (* dt (/ (the double-float (max 0d0 (- ybar k))) tau)))))
          (let ((new-damage
                  (max
                   damage
                   (cl-mpm/damage::damage-response-exponential k E  init-stress ductility))))
            (declare (double-float new-damage))
            (setf damage-inc (- new-damage damage)))

          (when (>= damage 1d0)
            (setf damage-inc 0d0)
            (setf ybar 0d0))
          (incf (the double-float (cl-mpm/particle::mp-time-averaged-damage-inc mp)) damage-inc)
          (incf (the double-float (cl-mpm/particle::mp-time-averaged-ybar mp)) ybar)
          (incf (the double-float (cl-mpm/particle::mp-time-averaged-counter mp)))
          ;;Transform to log damage
          (incf damage damage-inc)
          ;;Transform to linear damage
          (setf damage (max 0d0 (min 1d0 damage)))
          (when (> damage critical-damage)
            (setf damage 1d0))
          )
  (values)))
