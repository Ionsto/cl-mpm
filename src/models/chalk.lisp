(defpackage :cl-mpm/models/chalk
  (:use :cl
   :cl-mpm/utils)
  (:export))
(in-package :cl-mpm/models/chalk)
(defclass particle-chalk (particle-elastic-damage)
  (
   (coheasion
    :accessor mp-coheasion
    :initarg :coheasion
    :initform 0d0)
   (friction-angle
    :accessor mp-friction-angle
    :initarg :friction-angle
    :initform 30d0)
   (max-strain
    :initform 0d0
    )
   (ductility
    :accessor mp-ductility
    :initarg :ductility
    :initform 1d0)

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
   )
  (:documentation "A chalk damage model"))

(defclass particle-chalk-creep (particle-chalk)
  ()
  (:documentation "A chalk damage model"))

(defclass particle-chalk-brittle (particle-chalk particle-mc)
  (
   (fracture-energy
    :accessor mp-gf
    :initarg :fracture-energy
    :initform 1d0)
   (history-stress
    :accessor mp-history-stress
    :initform 0d0)
   (delay-time
    :accessor mp-delay-time
    :initform 1d0
    :initarg :delay-time
    )
   (ft
    :accessor mp-ft
    :initform 200d3
    :initarg :ft)
   (fc
    :accessor mp-fc
    :initform 200d3
    :initarg :fc)
   )
  (:documentation "A chalk damage model"))

(defclass particle-chalk-delayed (particle-chalk-brittle)
  (
   (delay-time
    :accessor mp-delay-time
    :initform 1d0
    :initarg :delay-time
    )
   (delay-exponent
    :accessor mp-delay-exponent
    :initform 1d0
    :initarg :delay-exponent
    )
   )
  (:documentation "A chalk damage model"))
(defclass particle-chalk-anisotropic (particle-chalk-delayed)
  ()
  (:documentation "A chalk damage model"))


(defmethod constitutive-model ((mp particle-chalk) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-accessors ((de mp-elastic-matrix)
                   (stress mp-stress)
                   (stress-u mp-undamaged-stress)
                   (strain mp-strain)
                   (damage mp-damage)
                   (pressure mp-pressure)
                   (def mp-deformation-gradient)
                   (pos mp-position)
                   (calc-pressure mp-pressure-func)
                   (rho mp-rho)
                   (ps-vm mp-strain-plastic-vm)
                   (plastic-strain mp-strain-plastic)
                   (yield-func mp-yield-func)
                   (enable-plasticity mp-enable-plasticity)
                   )
      mp
    (declare (function calc-pressure))
    ;;Train elastic strain - plus trail kirchoff stress
    (setf stress-u
          (cl-mpm/constitutive::linear-elastic-mat strain de))
    ;; (call-next-method)

    (if enable-plasticity
        (let* ((rho_0 rho)
               (rho_1 rho_0)
               (rho-d (+ rho_1 (* (- rho_0 rho_1) (- 1d0 damage))))
               ;; (rho-d rho_0)
               )
          (multiple-value-bind (sig eps-e f) (cl-mpm/constitutive::vm-plastic stress-u de strain rho)
            (setf stress
                  sig
                  plastic-strain (magicl:.- strain eps-e)
                  yield-func f
                  )
            (setf strain eps-e)
            )
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
        )
    (when (> damage 0.0d0)
      (let* ((j (magicl:det def))
             (degredation (expt (- 1d0 damage) 2d0))
             )
        ;; (setf stress (magicl:.+ (cl-mpm/constitutive::voight-eye p)
        ;;                         (magicl:scale! s (max 1d-9 degredation))
        ;;                         ))
        ;; (setf stress (magicl:scale stress (max 1d-9 degredation)))
        ;; (multiple-value-bind (l v) (cl-mpm/utils::eig
        ;;                             (magicl:scale! (voight-to-matrix stress) (/ 1d0 j)))
        ;;   (let* ((tp 0d0)
        ;;          ;(tp (funcall calc-pressure pos))
        ;;          (driving-pressure (* tp 1d0 (expt (min 1.00d0 damage) 1)))
        ;;          )
        ;;     (setf pressure tp)
        ;;     (loop for i from 0 to 2
        ;;           do
        ;;              (let* ((sii (nth i l))
        ;;                       (esii (- sii driving-pressure)))
        ;;                  (when (> esii 0d0)
        ;;                    ;;tensile damage -> unbounded
        ;;                    (setf (nth i l) (* sii (max 1d-8 degredation)))
        ;;                    ;; (setf (nth i l) (+ (nth i l) driving-pressure))
        ;;                    ;; (setf (nth i l) (* sii degredation))
        ;;                    )
        ;;                ;; (setf (nth i l) (* sii (max 1d-5 degredation)))
        ;;                  ;; (when (< esii 1d0)
        ;;                  ;;   ;;bounded compressive damage
        ;;                  ;;   (setf (nth i l) (* (nth i l) (max 1d-1 degredation)))
        ;;                  ;;   ;; (setf (nth i l) (+ (nth i l) driving-pressure))
        ;;                  ;;   )
        ;;                ;; (setf (nth i l) 0)
        ;;                  ;; (setf (nth i l) (* sii (max 0d0 (- 1d0 damage))))
        ;;                  )
        ;;           )
        ;;       (setf stress (magicl:scale! (matrix-to-voight (magicl:@ v
        ;;                                                               (magicl:from-diag l :type 'double-float)
        ;;                                                               (magicl:transpose v))) j))
        ;;       ))
        (let ((p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
              (s (cl-mpm/constitutive::deviatoric-voigt stress)))
          (setf stress (magicl:.+ (cl-mpm/constitutive::voight-eye p)
                                  (magicl:scale! s (max 0d-3 degredation))
                                  )))
        ))
    stress
    ))
(defmethod constitutive-model ((mp particle-chalk-brittle) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-accessors ((de mp-elastic-matrix)
                   (stress mp-stress)
                   (stress-u mp-undamaged-stress)
                   (strain mp-strain)
                   (damage mp-damage)
                   (coheasion mp-c)
                   (ps-vm mp-strain-plastic-vm)
                   (ps-vm-inc mp-strain-plastic-vm-inc)
                   (plastic-strain mp-strain-plastic)
                   (yield-func mp-yield-func)
                   (enable-plasticity mp-enable-plasticity)
                   (E mp-E)
                   (nu mp-nu)
                   (phi mp-phi)
                   (psi mp-psi)
                   (kc-r mp-k-compressive-residual-ratio)
                   (kt-r mp-k-tensile-residual-ratio)
                   (g-r mp-shear-residual-ratio))
      mp
    (declare (magicl:matrix/double-float de stress stress-u strain plastic-strain)
             (double-float coheasion ps-vm-inc ps-vm yield-func E nu phi psi kc-r kt-r g-r damage))
    ;;Train elastic strain - plus trail kirchoff stress
    (cl-mpm/constitutive::linear-elastic-mat strain de stress-u)
    ;; (setf stress-u (cl-mpm/constitutive::linear-elastic-mat strain de))

    (if enable-plasticity
        (progn
          (multiple-value-bind (sig eps-e f)
              (cl-mpm/constitutive::mc-plastic stress-u
                                               de
                                               strain
                                               E
                                               nu
                                               phi
                                               psi
                                               coheasion)
            (setf stress sig
                  plastic-strain (magicl:.- strain eps-e)
                  yield-func f)
            (setf strain eps-e))
          (let ((inc (multiple-value-bind (l v)
                         (cl-mpm/utils:eig (cl-mpm/utils:voigt-to-matrix (cl-mpm/particle::mp-strain-plastic mp)))
                       (destructuring-bind (s1 s2 s3) l
                         (sqrt
                          (/ (+ (expt (- s1 s2) 2d0)
                                (expt (- s2 s3) 2d0)
                                (expt (- s3 s1) 2d0)
                                ) 2d0))))))
            (incf ps-vm
                  inc)
            (setf ps-vm-inc inc)))
        (setf stress (cl-mpm/utils:voigt-copy stress-u)))
    (when (> damage 0.0d0)
      (let ((p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
            (s (cl-mpm/constitutive::deviatoric-voigt stress)))
        (setf p
              (if (> p 0d0)
                  (* (expt (- 1d0 (* (- 1d0 kt-r) damage)) 1d0) p)
                  ;; (* (expt (- 1d0 (* (- 1d0 kc-r) damage)) 1d0) p)
                  p
                  ))
        (magicl:.+ (cl-mpm/constitutive::voight-eye p)
                   (magicl:scale! s (expt (- 1d0 (* (- 1d0 g-r) damage)) 1d0))
                   stress)))
    ;; (when (> damage 0.0d0)
    ;;   (let* ()
    ;;     (multiple-value-bind (l v) (cl-mpm/utils::eig (voight-to-matrix stress))
    ;;       (let* ()
    ;;         (loop for i from 0 to 2
    ;;               do
    ;;                  (let* ((sii (nth i l))
    ;;                           (esii sii))
    ;;                      (when (> esii 0d0)
    ;;                        ;;Tensile damage -> unbounded
    ;;                        (setf (nth i l) (* esii (* (- 1d0 1d-15) (- 1d0 damage))))
    ;;                        )
    ;;                      (when (< esii 1d0)
    ;;                        ;;Bounded compressive damage
    ;;                        ;; (setf (nth i l) (* esii (* (- 1d0 1d-9) (- 1d0 damage))))
    ;;                        )
    ;;                      ;; (setf (nth i l) (* sii (max 0d0 (- 1d0 damage))))
    ;;                      )
    ;;               )
    ;;           (setf stress (matrix-to-voight (magicl:@ v
    ;;                                                                   (magicl:from-diag l :type 'double-float)
    ;;                                                                   (magicl:transpose v))))
    ;;           ))
    ;;     ))
    stress))
