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
        (magicl:scale! stress (max 0d-9 degredation)))
      ;; (let* ((j 1d0)
      ;;        (p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
      ;;        (s (cl-mpm/constitutive::deviatoric-voigt stress)))
      ;;   ;; (setf stress (magicl:.+ (cl-mpm/constitutive::voight-eye p)
      ;;   ;;                         (magicl:scale! s (max 1d-3 (expt (- 1d0 damage) 2d0)))
      ;;   ;;                         ))
      ;;   (multiple-value-bind (l v) (cl-mpm/utils::eig
      ;;                               (magicl:scale! (voight-to-matrix stress) (/ 1d0 j)))
      ;;     (let* (;(tp (funcall calc-pressure (magicl:tref pos 1 0) datum rho))
      ;;            (tp (funcall calc-pressure pos))
      ;;            (driving-pressure (* tp 1d0 (expt (min 1.00d0 damage) 1)))
      ;;            (degredation (expt (- 1d0 damage) 2d0)))
      ;;       ;; (setf stress (magicl:scale stress-undamaged (max 1d-3 degredation)))
      ;;       (setf pressure tp)
      ;;       (loop for i from 0 to 2
      ;;             do
      ;;                (let* ((sii (nth i l))
      ;;                         (esii (- sii driving-pressure)))
      ;;                    (when (> esii 0d0)
      ;;                      ;;Tensile damage -> unbounded
      ;;                      (setf (nth i l) (* esii (max 0d-6 degredation)))
      ;;                      (setf (nth i l) (+ (nth i l) driving-pressure))
      ;;                      )
      ;;                    ;; (when (< esii 1d0)
      ;;                    ;;   ;;Bounded compressive damage
      ;;                    ;;   (setf (nth i l) (* esii (max 1d0 degredation)))
      ;;                    ;;   (setf (nth i l) (+ (nth i l) driving-pressure))
      ;;                    ;;   )
      ;;                    ;; (setf (nth i l) (* sii (max 0d0 (- 1d0 damage))))
      ;;                    )
      ;;             )
      ;;         (setf stress (magicl:scale! (matrix-to-voight (magicl:@ v
      ;;                                                                 (magicl:from-diag l :type 'double-float)
      ;;                                                                 (magicl:transpose v))) j))
      ;;         ))
      ;;   )
      )
    stress
    ))
