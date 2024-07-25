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
   )
  (:documentation "A ice mp with viscoplastic damage and viscoelastic relaxation "))

(defclass particle-glen-damage (particle-glen particle-damage)
  ()
  (:documentation "A weakly compressible glen flow mp with damage mechanics"))

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
      (let* (;(eng-strain-rate (cl-mpm/fastmath::fast-.* (magicl:map (lambda (x) (* x (exp x))) strain) velocity-rate
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

        (cl-mpm/fastmath::fast-.+
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

    (setf stress-u
          (cl-mpm/constitutive::linear-elastic-mat strain de))

    ;; (let ((pressure (* pressure (expt damage 1))))
    ;;   (cl-mpm/fastmath::fast-.+ stress-u
    ;;                             (cl-mpm/utils::voigt-eye (/ pressure 3d0)))
    ;;   (setf strain (magicl:linear-solve de stress-u)))



    (when enable-plasticity
      (multiple-value-bind (sig eps-e f)
          (cl-mpm/constitutive::mc-plastic-terzaghi
           stress-u
           de
           strain
           E
           nu
           phi
           psi
           coheasion
           (* pressure damage 0d0)
           )
        (setf stress-u
              sig
              plastic-strain (magicl:.- strain eps-e)
              yield-func f)
        (setf strain eps-e))
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

    ;; (let ((pressure (* pressure (expt damage 1))))
    ;;   (cl-mpm/fastmath::fast-.+ stress
    ;;                             (cl-mpm/utils::voigt-eye (/ pressure 3d0))
    ;;                             stress))
    ;; ;; (when (> damage 0.0d0)
    ;; ;;   (cl-mpm/fastmath::fast-scale! stress (- 1d0 (* (- 1d0 1d-9) damage))))
    ;; ;; (when (> damage 0.0d0)
    ;; ;;   (multiple-value-bind (l v) (cl-mpm/utils::eig (voight-to-matrix stress))
    ;; ;;     (loop for i from 0 to 2
    ;; ;;           do (let* ((sii (nth i l)))
    ;; ;;                (when (> sii 0d0)
    ;; ;;                  ;;Tensile damage -> unbounded
    ;; ;;                  (setf (nth i l) (* sii (max 0d0 (- 1d0 damage))))
    ;; ;;                  )
    ;; ;;                (when (< sii 0d0)
    ;; ;;                  ;;Bounded compressive damage
    ;; ;;                  (setf (nth i l) (* sii (max 1d0 (- 1d0 (expt damage 2d0))))))))
    ;; ;;     (setf stress (matrix-to-voight (magicl:@ v
    ;; ;;                                              (cl-mpm/utils::matrix-from-diag l)
    ;; ;;                                              (magicl:transpose v))))))
    ;; (setf (cl-mpm/particle::mp-body-force mp)
    ;;       (cl-mpm/utils:vector-from-list (list 0d0 (* -1d0 1d3 -9.8d0) 0d0)))
    ;; (let ((pressure (* pressure damage (magicl:det def))))
    ;;   (cl-mpm/fastmath::fast-.+ stress
    ;;                             (cl-mpm/utils::voigt-eye pressure)
    ;;                             stress))
    (when (> damage 0.0d0)
      (let* ((d-p (* pressure damage 0d0 (magicl:det def))))
        (declare (double-float damage))
        (let* ((p (+ (/ (cl-mpm/constitutive::voight-trace stress) 3d0) d-p))
               (s (cl-mpm/constitutive::deviatoric-voigt stress)))
          ;; (incf p (* -1d0 pressure damage))
          (setf
           p
           (+ d-p
              (if (> p 0d0)
                  (* (- 1d0 (* (- 1d0 kt-r) damage)) p)
                  (* (- 1d0 (* (- 1d0 kc-r) damage)) p)
                  )))
          (cl-mpm/fastmath:fast-.+
           (cl-mpm/constitutive::voight-eye p)
           (magicl:scale! s (- 1d0 (* (- 1d0 g-r) damage)))
           stress)
          )))

    (let ((pressure (* pressure (expt damage 1) (magicl:det def))))
      (cl-mpm/fastmath::fast-.+ stress
                                (cl-mpm/utils::voigt-eye pressure)
                                stress))

    (let ((damping-coeff 1d-5))
      (setf (cl-mpm/particle::mp-body-force mp)
            (cl-mpm/fastmath::fast-scale-vector (cl-mpm/particle::mp-velocity mp)
                                                (*  1d0
                                                    damping-coeff
                                                    (cl-mpm/particle::mp-boundary mp)
                                                    (cl-mpm/particle::mp-mass mp)
                                                    ))))

    ;; (cl-mpm/fastmath::fast-.+ stress
    ;;                           (cl-mpm/utils::voigt-eye (/ (* pressure damage) 3d0)))

    ;; (when (> damage 0.0d0)
    ;;   (let* ((pressure (* pressure damage)))
    ;;     (declare (double-float damage))
    ;;     (let* ((pressure pressure)
    ;;            (p (- (/ (cl-mpm/constitutive::voight-trace stress) 3d0) pressure))
    ;;            (s (cl-mpm/constitutive::deviatoric-voigt stress)))
    ;;       (setf p
    ;;             (+ pressure
    ;;                (if (> p 0d0)
    ;;                    (* (- 1d0 damage) p)
    ;;                    p
    ;;                    ;; (* (- 1d0 (expt damage 2)) p)
    ;;                    )))
    ;;       (cl-mpm/fastmath:fast-.+
    ;;        (cl-mpm/constitutive::voight-eye p)
    ;;        (magicl:scale! s (- 1d0 damage))
    ;;        stress
    ;;        ))))0

    ;; (when (> damage 0.0d0)
    ;;   (let* ((beta-I (cl-mpm/damage::damage-response-exponential-beta E initiation-stress ductility-I))
    ;;          (beta-II (cl-mpm/damage::damage-response-exponential-beta E initiation-stress ductility-II))
    ;;          (damage-II (if (< damage 1d0)
    ;;                         (- 1d0 (exp (* (/ beta-II beta-I) (log (- 1d0 damage)))))
    ;;                         1d0)))
    ;;     (declare (double-float damage))
    ;;     (setf (mp-damage-mode-2 mp) damage-II)
    ;;     (let* ((p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
    ;;            (s (cl-mpm/constitutive::deviatoric-voigt stress)))
    ;;       (setf p
    ;;             (if (> p 0d0)
    ;;                 (* (- 1d0 damage) p)
    ;;                 ;; (* (- 1d0 (* (- 1d0 kc-r) damage)) p)
    ;;                 p
    ;;                 )
    ;;             )
    ;;       (setf stress
    ;;             (magicl:.+
    ;;              (cl-mpm/constitutive::voight-eye p)
    ;;              (magicl:scale! s (- 1d0 damage-II))
    ;;              )))))
    stress
    ))
