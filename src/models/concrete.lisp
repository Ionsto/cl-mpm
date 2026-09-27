;; (defpackage :cl-mpm/models/concrete
;;   (:use
;;    :cl
;;    :cl-mpm/utils
;;    :cl-mpm/particle)
;;   (:export))
;; (in-package :cl-mpm/models/concrete)
(in-package :cl-mpm/particle)

(defclass particle-concrete (particle-elastic-damage)
  ((fracture-energy
    :accessor mp-gf
    :initarg :fracture-energy
    :initform 1d0)
   (ductility
    :accessor mp-ductility
    :initarg :ductility
    :initform 0d0)
   (history-stress
    :accessor mp-history-stress
    :initform 0d0)
   (dissipated-energy
    :accessor mp-dissipated-energy
    :initform 0d0)
   (dissipated-energy-inc
    :accessor mp-dissipated-energy-inc
    :initform 0d0)
   )
  (:documentation "A concrete damage model"))
(defmethod constitutive-model ((mp particle-concrete) strain dt)
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
    ;; Non-objective stress intergration
    (setf stress-undamaged (cl-mpm/constitutive::linear-elastic-mat strain de))

    (setf stress (magicl:scale stress-undamaged 1d0))
    (when (> damage 0.0d0)
      (let ((degredation (expt (- 1d0 damage) 1d0)))
        (magicl:scale! stress (max 0d-9 degredation))))
    stress
    ))

(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-concrete) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (pressure cl-mpm/particle::mp-pressure)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     (angle cl-mpm/particle::mp-friction-angle)
                     (c cl-mpm/particle::mp-coheasion)
                     (J cl-mpm/particle::mp-deformation-jacobian-strain)
                     ) mp
      (declare (double-float pressure damage))
        (progn
          (when (< damage 1d0)
            (let ((cauchy-undamaged (magicl:scale stress (/ 1d0 J))))
              (multiple-value-bind (s_1 s_2 s_3) (cl-mpm/utils::principal-stresses-3d cauchy-undamaged)
                (let* (;(s_1 (max 0d0 s_1))
                       )
                  (when (> s_1 0d0)
                    ;; (setf damage-increment s_1)
                    (setf damage-increment (sqrt
                                            (+ (expt s_1 2)
                                               (expt s_2 2)
                                               (expt s_3 2))))
                    )))))
          (when (>= damage 1d0)
            (setf damage-increment 0d0))
          ;;Delocalisation switch
          (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
          (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
          ))))
