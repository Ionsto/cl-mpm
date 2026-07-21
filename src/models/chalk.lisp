(defpackage :cl-mpm/models/chalk
  (:use :cl
   :cl-mpm/utils)
  (:export))
;; (in-package :cl-mpm/models/chalk)
;;We don't actually want to intern anything in chalk i guess
(declaim (optimize (speed 3) (debug 0) (safety 0)))
;; (declaim (optimize (speed 0) (debug 3) (safety 3)))

(in-package :cl-mpm/particle)
(declaim (optimize (speed 3) (debug 0) (safety 0)))
;; (declaim (optimize (speed 0) (debug 3) (safety 3)))
(defclass particle-chalk (particle-plastic-damage-frictional)
  ()
  (:documentation "A chalk damage model"))

(defclass particle-chalk-creep (particle-chalk)
  ()
  (:documentation "A chalk damage model"))

(defclass particle-chalk-brittle (particle-chalk particle-mc)
  ()
  (:documentation "A chalk damage model"))

(defclass particle-chalk-delayed (particle-chalk-brittle)
  ((delay-time
    :accessor mp-delay-time
    :initform 1d0
    :initarg :delay-time)
   (delay-exponent
    :accessor mp-delay-exponent
    :initform 1d0
    :initarg :delay-exponent))
  (:documentation "A chalk damage model"))

(defclass particle-chalk-delayed-grassl (particle-chalk-delayed)
  ()
  (:documentation "A chalk damage model based on grassl"))

(defclass particle-chalk-delayed-linear (particle-chalk-delayed)
  ()
  (:documentation "A chalk damage model based on grassl"))

(defclass particle-chalk-anisotropic (particle-chalk-delayed)
  ((damage-tensor
    :accessor mp-damage-tensor
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils::matrix-zeros))
   (damage-ybar-tensor
    :accessor mp-damage-ybar-tensor
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils::matrix-zeros)))
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
                   (enable-plasticity mp-enable-plasticity))
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
        (let ((p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
              (s (cl-mpm/constitutive::deviatoric-voigt stress)))
          (setf stress (magicl:.+ (cl-mpm/constitutive::voight-eye p)
                                  (magicl:scale! s (max 0d-3 degredation))
                                  )))
        ))
    stress
    ))
;; (defmethod initialize-instance :after ((mp particle-chalk-brittle) &key)
;;   (with-accessors ((ductility cl-mpm/particle::mp-ductility)
;;                    (angle cl-mpm/particle::mp-friction-angle)
;;                    (angle-r cl-mpm/particle::mp-residual-friction)
;;                    (rc mp-shear-residual-ratio)
;;                    (init-stress mp-initiation-stress)
;;                    (oversize mp-oversize-scale))
;;       mp
;;     (let* ((c (cl-mpm/damage::mohr-coloumb-tensile-to-coheasion init-stress (rad-to-deg angle)))
;;            (rs (cl-mpm/damage::est-shear-from-angle (rad-to-deg angle) (rad-to-deg angle-r) rc))
;;            ;; (residual-strength (mp-residual-strength mp))
;;            (oversize-ratio (cl-mpm/damage::compute-oversize-factor oversize ductility)))
;;       (setf
;;        (cl-mpm/particle::mp-phi mp) angle
;;        (mp-c mp) (* oversize-ratio c)
;;        (mp-shear-residual-ratio mp) rs)
;;       ;; (setf (mp-shear-residual-ratio mp) (min (mp-shear-residual-ratio mp) residual-strength)
;;       ;;       (mp-k-tensile-residual-ratio mp) (min residual-strength (mp-k-tensile-residual-ratio mp)))
;;       )))
(defmethod constitutive-model ((mp particle-chalk-brittle) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-accessors ((de mp-elastic-matrix)
                   (stress mp-stress)
                   (stress-u mp-undamaged-stress)
                   (strain mp-strain)
                   (damage mp-damage)
                   (damage-t mp-damage-tension)
                   (damage-c mp-damage-compression)
                   (damage-s mp-damage-shear)
                   (coheasion mp-c)
                   (ps-vm mp-strain-plastic-vm)
                   (ps-vm-1 mp-strain-plastic-vm-1)
                   (ps-vm-inc mp-strain-plastic-vm-inc)
                   (plastic-strain mp-strain-plastic)
                   (yield-func mp-yield-func)
                   (enable-plasticity mp-enable-plasticity)
                   (enable-damage mp-enable-damage)
                   (E mp-E)
                   (nu mp-nu)
                   (phi mp-phi)
                   (psi mp-psi)
                   (p-mod mp-p-modulus-0)
                   )
      mp
    (declare (magicl:matrix/double-float de stress stress-u strain plastic-strain)
             (double-float coheasion ps-vm-inc ps-vm yield-func E nu phi psi kc-r kt-r g-r damage))
    ;;Train elastic strain - plus trail kirchoff stress

    (setf stress-u (cl-mpm/constitutive::linear-elastic-mat strain de stress-u))

    (setf p-mod (cl-mpm/particle::compute-p-modulus mp))

    (when enable-plasticity
        (progn
          (multiple-value-bind (sig eps-e f inc pmod)
              (cl-mpm/ext::constitutive-mohr-coulomb stress-u
                                                     de
                                                     strain
                                                     E
                                                     nu
                                                     phi
                                                     psi
                                                     coheasion)
            (setf
             stress-u sig
             strain eps-e
             yield-func f
             p-mod pmod)
            (let ()
              (setf ps-vm (+ ps-vm-1 inc))
              (setf ps-vm-inc inc)))))
    (cl-mpm/utils:voigt-copy-into stress-u stress)
    stress))

(in-package :cl-mpm/damage)


(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-chalk-brittle) dt)
  (with-accessors ((p cl-mpm/particle::mp-pressure)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (stress cl-mpm/particle::mp-stress)
                   (damage cl-mpm/particle::mp-damage)
                   (enable-damage cl-mpm/particle::mp-enable-damage)
                   )
      mp
    ;; (cl-mpm/damage::apply-isotropic-degredation mp)
    ;; (cl-mpm/damage::apply-vol-tensile-degredation mp)
    (cl-mpm/damage::apply-tcs-degredation mp)
    ))

(defmethod update-damage ((mp cl-mpm/particle::particle-chalk-delayed) dt)
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
        (setf (cl-mpm/particle::mp-damage-prev-trial mp) (cl-mpm/particle::mp-damage mp))
        (setf damage-inc 0d0)
        (let ((a tau-exp)
              (k0 init-stress))
          (when t;; (or (>= ybar-prev k0)
                ;;     (>= ybar k0))
            (setf k
                  (max
                   k-n
                   ;; (cl-mpm/damage::secant-solver
                   ;;  k-n
                   ;;  ybar-prev
                   ;;  ybar
                   ;;  dt
                   ;;  (lambda (kmid ymid)
                   ;;    (cl-mpm/damage::deriv-partial
                   ;;     kmid
                   ;;     ymid
                   ;;     k0
                   ;;     tau
                   ;;     tau-exp)))
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
                       s-dt)))
                   ))))
        (compute-damage mp)
        (setf damage-inc (- damage damage-n))
        (setf (cl-mpm/particle::mp-time-averaged-damage-inc mp) (/ damage-inc dt))
        (incf (the double-float (cl-mpm/particle::mp-time-averaged-counter mp)))
        (setf damage (max 0d0 (min 1d0 damage))))
      (values))))
;; (defmethod update-damage ((mp cl-mpm/particle::particle-chalk-delayed-linear) dt)
;;   (when (cl-mpm/particle::mp-enable-damage mp)
;;     (with-accessors ((stress cl-mpm/particle:mp-stress)
;;                      (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
;;                      (damage cl-mpm/particle:mp-damage)
;;                      (E cl-mpm/particle::mp-e)
;;                      (Gf cl-mpm/particle::mp-Gf)
;;                      (damage-inc cl-mpm/particle::mp-damage-increment)
;;                      (ybar cl-mpm/particle::mp-damage-ybar)
;;                      (init-stress cl-mpm/particle::mp-initiation-stress)
;;                      (damage-rate cl-mpm/particle::mp-damage-rate)
;;                      (critical-damage cl-mpm/particle::mp-critical-damage)
;;                      (pressure cl-mpm/particle::mp-pressure)
;;                      (def cl-mpm/particle::mp-deformation-gradient)
;;                      ;; (length-t cl-mpm/particle::mp-true-local-length)
;;                      (length cl-mpm/particle::mp-local-length)
;;                      (k cl-mpm/particle::mp-history-stress)
;;                      (tau cl-mpm/particle::mp-delay-time)
;;                      (tau-exp cl-mpm/particle::mp-delay-exponent)
;;                      (ductility cl-mpm/particle::mp-ductility)
;;                      (nu cl-mpm/particle::mp-nu)
;;                      (p cl-mpm/particle::mp-p-modulus)
;;                      (kc-r cl-mpm/particle::mp-k-compressive-residual-ratio)
;;                      (kt-r cl-mpm/particle::mp-k-tensile-residual-ratio)
;;                      (g-r cl-mpm/particle::mp-shear-residual-ratio)
;;                      (damage-tension cl-mpm/particle::mp-damage-tension)
;;                      (damage-shear cl-mpm/particle::mp-damage-shear)
;;                      (damage-compression cl-mpm/particle::mp-damage-compression)
;;                      ) mp
;;       (declare (double-float damage damage-inc critical-damage k ybar tau dt
;;                              damage-tension damage-shear damage-compression
;;                              kt-r kc-r g-r
;;                              init-stress ductility))
;;       (when (< damage 1d0)
;;         ;;Damage increment holds the delocalised driving factor
;;         (setf damage-inc 0d0)
;;         (let ((a tau-exp)
;;               (k0 init-stress))
;;           (incf k (the double-float
;;                        (*
;;                         dt
;;                         (/
;;                          (* k0
;;                             (the double-float
;;                                  (expt
;;                                   (/ (the double-float (max 0d0 (- ybar k)))
;;                                      k0) a)))
;;                          tau)))))
;;         (let ((new-damage
;;                 (max
;;                  damage
;;                  ;; (damage-response-exponential k E init-stress ductility)
;;                  (damage-response-linear k E Gf (/ length (the double-float (sqrt 7d0))) init-stress ductility)
;;                  )))
;;           (declare (double-float new-damage))
;;           (setf damage-inc (- new-damage damage)))
;;         (setf
;;          damage-tension (max damage-tension (damage-response-linear-residual k E init-stress ductility kt-r))
;;          damage-shear (max damage-shear (damage-response-linear-residual k E init-stress ductility g-r))
;;          damage-compression (max damage-compression (damage-response-linear-residual k E init-stress ductility kc-r)))

;;         (when (>= damage 1d0)
;;           (setf damage-inc 0d0))
;;         (incf (the double-float(cl-mpm/particle::mp-time-averaged-damage-inc mp)) (* damage-inc dt))
;;         (incf (the double-float(cl-mpm/particle::mp-time-averaged-ybar mp)) ybar)
;;         (incf (the double-float(cl-mpm/particle::mp-time-averaged-counter mp)))
;;         ;;Transform to log damage
;;         (incf damage damage-inc)
;;         ;;Transform to linear damage
;;         (setf damage (max 0d0 (min 1d0 damage)))
;;         (when (> damage critical-damage)
;;           (setf damage 1d0)
;;           (setf damage-inc 0d0))
;;         )
;;       (values))))

(defmethod update-damage ((mp cl-mpm/particle::particle-chalk-anisotropic) dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (damage-tensor cl-mpm/particle::mp-damage-tensor)
                     (ybar-tensor cl-mpm/particle::mp-damage-ybar-tensor)
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
                     (ductility cl-mpm/particle::mp-ductility)
                     (D cl-mpm/particle::mp-stretch-tensor)
                     ) mp
      (declare (double-float damage damage-inc critical-damage))
        (progn
          (magicl:scale! ybar-tensor 0d0)
          (when (< damage 1d0)
            (let ((damage-inc-mat (cl-mpm/utils:matrix-zeros))
                  (cauchy-undamaged (magicl:scale stress (/ 1d0 (magicl:det def))))
                  (anisotropicity 0.0d0)
                  )
              (multiple-value-bind (ls v) (cl-mpm/utils:eig (cl-mpm/utils:voigt-to-matrix cauchy-undamaged))
                (let ((l-y (mapcar (lambda (sii) (if (> sii 0d0)
                                                     damage-inc
                                                     0d0)) ls)))
                  (loop for i from 0 to 2
                        do
                           (let* ((sii (+ (* (- 1d0 anisotropicity) (nth i l-y))
                                          (reduce #'+ (mapcar (lambda (z) (* z anisotropicity (/ 1d0 3d0))) l-y))))
                                  (vii (magicl::column v i))
                                  (vsi (magicl:@ vii (magicl:transpose vii)))
                                  (dii (magicl::trace (magicl:@ damage-tensor vsi)))
                                  )
                             ;; (setf sii ybar)
                             ;; (when (< sii 0d0)
                             ;;   (setf sii 0d0))
                             ;; (when (> sii 0d0)
                             ;;   (setf sii damage-inc))
                             (let* ((new-damage (damage-response-exponential sii E init-stress ductility))
                                    (damage-increment ;(- (max dii new-damage) dii)
                                      (* (/ dt tau) (- 1d0 (exp (- (* 1d0 (abs (- new-damage dii)))))))))
                               (magicl:.+ ybar-tensor
                                          (magicl:scale vsi sii)
                                          ybar-tensor)
                               (magicl:.+ damage-inc-mat
                                          (magicl:scale vsi damage-increment)
                                          damage-inc-mat))))))
              ;; (break)
              (let ((omega (magicl:scale! (magicl:.- D (magicl:transpose D)) 0.5d0))
                    )
                (magicl:.+ damage-tensor
                           (magicl:.+ damage-inc-mat
                                      (magicl:.+ (magicl:@ omega damage-tensor)
                                                 (magicl:@ damage-tensor omega)))
                           damage-tensor))))
          ;;Damage increment holds the delocalised driving factor
          ;; (setf k (max k ybar))
          ;; (setf damage-inc 0d0)

          ;; (incf k (* dt (/ (max 0d0 (- ybar k)) tau)))

          ;; (let ((new-damage (max damage
          ;;                        ;; (test-damage k 1d7 init-stress)
          ;;                        ;; (delayed-damage-y ybar E Gf length init-stress)
          ;;                        (damage-response-exponential k E Gf length init-stress)
          ;;                        )))
          ;;   (setf damage-inc (- new-damage damage))
          ;;   ;; (setf damage-inc (* (/ dt tau) (- 1d0 (exp (- (* 1d0 (abs (- new-damage damage))))))))
          ;;   )
          ;; (setf damage (max damage (brittle-chalk-d k E Gf length init-stress))
          ;;       damage-inc 0d0)
          ;; (when (>= damage 1d0)
          ;;   (setf damage-inc 0d0)
          ;;   (setf ybar 0d0))
          ;; (incf (cl-mpm/particle::mp-time-averaged-damage-inc mp) damage-inc)
          ;; (incf (cl-mpm/particle::mp-time-averaged-ybar mp) ybar)
          ;; (incf (cl-mpm/particle::mp-time-averaged-counter mp))
          ;;Transform to log damage
          ;; (incf damage damage-inc)
          ;;Transform to linear damage
          (multiple-value-bind (ls v) (cl-mpm/utils:eig damage-tensor)
            (loop for i from 0 to 2
                  do
                     (let* ()
                       (setf (nth i ls) (max 0d0 (min 1d0 (nth i ls))))
                       (when (> (nth i ls) critical-damage)
                         (setf (nth i ls) 1d0))
                       (setf damage (max damage (nth i ls)))
                       )
                  )
            (setf damage-tensor (magicl:@ v
                                          (magicl:from-diag ls :type 'double-float)
                                          (magicl:transpose v)))
            )
          )
  (values)
  ))

;; (defmethod update-damage ((mp cl-mpm/particle::particle-chalk) dt)
;;     (with-accessors ((stress cl-mpm/particle:mp-stress)
;;                      (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
;;                      (damage cl-mpm/particle:mp-damage)
;;                      (damage-inc cl-mpm/particle::mp-damage-increment)
;;                      (ybar cl-mpm/particle::mp-damage-ybar)
;;                      (init-stress cl-mpm/particle::mp-initiation-stress)
;;                      (damage-rate cl-mpm/particle::mp-damage-rate)
;;                      (critical-damage cl-mpm/particle::mp-critical-damage)
;;                      (pressure cl-mpm/particle::mp-pressure)
;;                      (def cl-mpm/particle::mp-deformation-gradient)
;;                      ) mp
;;       (declare (double-float damage damage-inc critical-damage))
;;         (progn
;;           ;;Damage increment holds the delocalised driving factor
;;           (when (< damage 1d0)
;;             (setf damage-inc (* dt
;;                                 ;; (/ 1d0 (- 1d0 damage))
;;                                 (damage-rate-profile-chalk damage-inc damage damage-rate init-stress))))
;;           (when (>= damage 1d0)
;;             ;; (setf damage-inc 0d0)
;;             (setf ybar 0d0))
;;           (incf (cl-mpm/particle::mp-time-averaged-damage-inc mp) damage-inc)
;;           (incf (cl-mpm/particle::mp-time-averaged-ybar mp) ybar)
;;           (incf (cl-mpm/particle::mp-time-averaged-counter mp))
;;           ;;Transform to log damage
;;           (incf damage damage-inc)
;;           ;;Transform to linear damage
;;           (setf damage (max 0d0 (min 1d0 damage)))
;;           (when (> damage critical-damage)
;;             (setf damage 1d0)
;;             (setf damage-inc 0d0)))
;;   (values)
;;   ))

(defmethod damage-model-calculate-y ((mp cl-mpm/particle::particle-chalk) dt)
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
                     ) mp
      (declare (double-float pressure damage))
        (progn
          (when (< damage 1d0)
            (let ((cauchy-undamaged (magicl:scale stress (/ 1d0 (magicl:det def)))))
              (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d cauchy-undamaged)
                (let* ((pressure-effective (* 1d0 damage pressure))
                       (j2 (sqrt (cl-mpm/constitutive::voigt-j2
                                  (cl-mpm/utils::deviatoric-voigt cauchy-undamaged))))
                       (p (+ s_1 s_2 s_3))
                       (s_1 (- s_1 pressure-effective))
                       (s_2 (- s_2 pressure-effective))
                       (s_3 (- s_3 pressure-effective))
                       (s_1 (max 0d0 s_1))
                       (s_2 (max 0d0 s_2))
                       (s_3 (max 0d0 s_3))
                       (angle (* angle (/ pi 180d0)))
                       ;; (c 1d4)
                       (A (/ (* 6 c (cos angle))
                             (* (sqrt 3) (- 3 (sin angle)))))
                       (B (/ (* 2d0 (sin angle)) (* (sqrt 3) (- 3 (sin angle)))))
                       ;; (s_1 (+ j2 (- (* B p) A)))
                       ;; (s_1 j2)
                       )
                  (when (> s_1 0d0)
                    ;; (setf damage-increment (* s_1 (expt (- 1d0 damage) -2d0)))
                    (setf damage-increment s_1)
                    )))))
          (when (>= damage 1d0)
            ;; (setf damage-increment 0d0)
            )
          ;;Delocalisation switch
          (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
          ))))

(defmethod damage-model-calculate-y ((mp cl-mpm/particle::particle-chalk-creep) dt)
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
                     ) mp
      (declare (double-float pressure damage))
        (progn
          (when (< damage 1d0)
            (let ((cauchy-undamaged (magicl:scale stress (/ 1d0 (magicl:det def)))))
              (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d cauchy-undamaged)
                (let* ((pressure-effective (* 1d0 damage pressure))
                       (j2 (sqrt (cl-mpm/constitutive::voigt-j2
                                  (cl-mpm/utils::deviatoric-voigt cauchy-undamaged))))
                       (p (+ s_1 s_2 s_3))
                       (s_1 (- s_1 pressure-effective))
                       (s_2 (- s_2 pressure-effective))
                       (s_3 (- s_3 pressure-effective))
                       (s_1 (max 0d0 s_1))
                       (s_2 (max 0d0 s_2))
                       (s_3 (max 0d0 s_3))
                       (angle (* angle (/ pi 180d0)))
                       ;; (c 1d4)
                       (A (/ (* 6 c (cos angle))
                             (* (sqrt 3) (- 3 (sin angle)))))
                       (B (/ (* 2d0 (sin angle)) (* (sqrt 3) (- 3 (sin angle)))))
                       ;; (s_1 (+ j2 (- (* B p) A)))
                       ;; (s_1 j2)
                       )
                  (when (> s_1 0d0)
                    ;; (setf damage-increment (* s_1 (expt (- 1d0 damage) -2d0)))
                    (setf damage-increment s_1)
                    )))))
          (when (>= damage 1d0)
            (setf damage-increment 0d0))
          ;;Delocalisation switch
          (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
          ))))

(defmethod update-damage ((mp cl-mpm/particle::particle-chalk-creep) dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (damage-inc cl-mpm/particle::mp-damage-increment)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (pressure cl-mpm/particle::mp-pressure)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     ) mp
      (declare (double-float damage damage-inc critical-damage))
        (progn
          ;;Damage increment holds the delocalised driving factor
          (when (< damage 1d0)
            (setf damage-inc (* dt
                                ;; (/ 1d0 (- 1d0 damage))
                                (damage-rate-profile-chalk damage-inc damage damage-rate init-stress))))
          (when (>= damage 1d0)
            (setf damage-inc 0d0)
            (setf ybar 0d0))
          (incf (cl-mpm/particle::mp-time-averaged-damage-inc mp) damage-inc)
          (incf (cl-mpm/particle::mp-time-averaged-ybar mp) ybar)
          (incf (cl-mpm/particle::mp-time-averaged-counter mp))
          ;;Transform to log damage
          (incf damage damage-inc)
          ;;Transform to linear damage
          (setf damage (max 0d0 (min 1d0 damage)))
          (when (> damage critical-damage)
            (setf damage 1d0)
            (setf damage-inc 0d0)))
  (values)
  ))

(defmethod damage-model-calculate-y ((mp cl-mpm/particle::particle-chalk-brittle) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
                     (strain cl-mpm/particle::mp-strain)
                     (damage cl-mpm/particle:mp-damage)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (pressure cl-mpm/particle::mp-pressure)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     (angle cl-mpm/particle::mp-friction-angle)
                     (c cl-mpm/particle::mp-coheasion)
                     (nu cl-mpm/particle::mp-nu)
                     (ft cl-mpm/particle::mp-ft)
                     (fc cl-mpm/particle::mp-fc)
                     (E cl-mpm/particle::mp-e)
                     (de cl-mpm/particle::mp-elastic-matrix)
                     (kc-r cl-mpm/particle::mp-k-compressive-residual-ratio)
                     (kt-r cl-mpm/particle::mp-k-tensile-residual-ratio)
                     (g-r cl-mpm/particle::mp-shear-residual-ratio)
                     ) mp
      (declare (double-float pressure damage))
      (progn
        (when (< damage 1d0)
          ;(setf damage-increment (tensile-energy-norm strain E de))
          (setf damage-increment
                (max 0d0
                     (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile
                      stress
                      (* angle (/ pi 180d0)))))
          )
        (when (>= damage 1d0)
          (setf damage-increment 0d0))
        ;;Delocalisation switch
        (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
        (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
        ))))

;; (defun set-damage-directly (mp damage)
;;   (let* ((ft init-stress)
;;          (e0 (/ ft E))
;;          (ef (/ (* ft (+ ductility 1d0)) (* 2d0 E)))
;;          (k (/ stress E))
;;          (beta (/ 1d0 (- ef e0)))
;;          )
;;     (declare (double-float ft e0 ef k beta))
;;     (if (> k e0)
;;         (- 1d0 (* (/ e0 k) (exp (- (* beta (- k e0))))))
;;         0d0))
;;   )

(defmethod update-damage ((mp cl-mpm/particle::particle-chalk-delayed-grassl) dt)
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
                     (p cl-mpm/particle::mp-p-modulus)
                     ) mp
      (declare (double-float damage damage-inc critical-damage k ybar tau dt))
        (progn
          ;;Damage increment holds the delocalised driving factor
          (setf damage-inc 0d0)
          (let ((a tau-exp)
                (k0 init-stress))
            (incf k (the double-float
                         (*
                          dt
                          (/
                           (* k0
                              (expt
                               (/ (the double-float (max 0d0 (- ybar k)))
                                  k0) a))
                             tau)))))
          (let ((new-damage
                  (max
                   damage
                   (plastic-damage-response-exponential k E (/ length (the double-float (sqrt 7d0))) ductility)
                   ;; (damage-response-exponential k E Gf (/ length (the double-float (sqrt 7d0))) init-stress ductility)
                   )))
            (declare (double-float new-damage))
            (setf damage-inc (- new-damage damage))
            )

          (when (>= damage 1d0)
            (setf damage-inc 0d0))
          (incf (the double-float(cl-mpm/particle::mp-time-averaged-damage-inc mp)) (* damage-inc dt))
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

(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-chalk-delayed-grassl) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
                     (strain cl-mpm/particle::mp-strain)
                     (plastic-strain cl-mpm/particle::mp-strain-plastic)
                     (plastic-vm cl-mpm/particle::mp-strain-plastic-vm)
                     (damage cl-mpm/particle:mp-damage)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (pressure cl-mpm/particle::mp-pressure)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     (angle cl-mpm/particle::mp-friction-angle)
                     (c cl-mpm/particle::mp-coheasion)
                     (nu cl-mpm/particle::mp-nu)
                     (ft cl-mpm/particle::mp-ft)
                     (fc cl-mpm/particle::mp-fc)
                     (E cl-mpm/particle::mp-e)
                     (de cl-mpm/particle::mp-elastic-matrix)
                     (kc-r cl-mpm/particle::mp-k-compressive-residual-ratio)
                     (kt-r cl-mpm/particle::mp-k-tensile-residual-ratio)
                     (g-r cl-mpm/particle::mp-shear-residual-ratio))
        mp
      (declare (double-float pressure damage E plastic-vm damage-increment))
      (progn
        (setf damage-increment (* E plastic-vm))
        ;;Delocalisation switch
        (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
        (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
        ))))
