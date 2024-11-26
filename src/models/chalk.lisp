(defpackage :cl-mpm/models/chalk
  (:use :cl
   :cl-mpm/utils)
  (:export))
;; (in-package :cl-mpm/models/chalk)
;;We don't actually want to intern anything in chalk i guess
(declaim (optimize (speed 3) (debug 0) (safety 0)))

(in-package :cl-mpm/particle)
(declaim (optimize (speed 3) (debug 0) (safety 0)))
(defclass particle-chalk (particle-elastic-damage)
  ((coheasion
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
    :initform 0d0))
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
   (peerlings-damage
    :accessor mp-peerlings-damage
    :initform t
    :initarg :peerlings-damage))
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
(defclass particle-chalk-delayed-grassl (particle-chalk-delayed)
  ()
  (:documentation "A chalk damage model based on grassl"))
(defclass particle-chalk-delayed-linear (particle-chalk-delayed)
  ()
  (:documentation "A chalk damage model based on grassl"))
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
                   (damage-t mp-damage-tension)
                   (damage-c mp-damage-compression)
                   (damage-s mp-damage-shear)
                   (coheasion mp-c)
                   (ps-vm mp-strain-plastic-vm)
                   (ps-vm-inc mp-strain-plastic-vm-inc)
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
                   )
      mp
    (declare (magicl:matrix/double-float de stress stress-u strain plastic-strain)
             (double-float coheasion ps-vm-inc ps-vm yield-func E nu phi psi kc-r kt-r g-r damage))
    ;;Train elastic strain - plus trail kirchoff stress
    (setf stress-u (cl-mpm/constitutive::linear-elastic-mat strain de stress-u))
    (when enable-plasticity
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
            (setf stress-u sig
                  plastic-strain (magicl:.- strain eps-e)
                  yield-func f)
            (setf strain eps-e)
            )
          (let ((inc (multiple-value-bind (l v)
                         (cl-mpm/utils:eig (cl-mpm/utils:voigt-to-matrix (cl-mpm/particle::mp-strain-plastic mp)))
                       (destructuring-bind (s1 s2 s3) l
                         (sqrt
                          (/ (+ (expt (- s1 s2) 2d0)
                                (expt (- s2 s3) 2d0)
                                (expt (- s3 s1) 2d0)
                                ) 2d0))))))
            (incf ps-vm inc)
            (setf ps-vm-inc inc))))
    (cl-mpm/utils:voigt-copy-into stress-u stress)
    (when (and
           enable-damage
           (> damage 0.0d0))
      (unless peerlings
        (setf damage-t (* kt-r damage)
              damage-c (* kc-r damage)
              damage-s (* g-r damage)))
      (let ((p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
            (s (cl-mpm/constitutive::deviatoric-voigt stress)))
        (declare (double-float damage-t damage-c damage-s))
        (setf p
              (if (> p 0d0)
                  (* (- 1d0 damage-t) p)
                  (* (- 1d0 damage-c) p)))
        (setf stress
              (cl-mpm/fastmaths:fast-.+
               (cl-mpm/constitutive::voight-eye p)
               (cl-mpm/fastmaths:fast-scale! s (- 1d0 damage-s))
               stress)))
      )
    stress))


(defmethod constitutive-model ((mp particle-chalk-anisotropic) strain dt)
  ;; (declare (optimize (speed 0) (debug 3)))
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
                   (damage-tensor mp-damage-tensor)
                   (E mp-E)
                   (nu mp-nu)
                   (coheasion mp-c)
                   (phi mp-phi)
                   (psi mp-psi)
                   )
      mp
    (declare (function calc-pressure))
    ;;Train elastic strain - plus trail kirchoff stress
    (setf stress-u
          (cl-mpm/constitutive::linear-elastic-mat strain de))
    ;; (call-next-method)

    (if enable-plasticity
        (progn
          (multiple-value-bind (sig eps-e f)
              ;; (cl-mpm/constitutive::vm-plastic stress-u de strain coheasion)
            (cl-mpm/constitutive::mc-plastic stress-u de strain
                                             E
                                             nu
                                             phi
                                             psi
                                             coheasion)
            (setf stress
                  sig
                  plastic-strain (magicl:.- strain eps-e)
                  yield-func f
                  )
            ;; (cl-mpm/utils:voigt-copy-into stress sig)
            ;; (magicl:.- strain eps-e plastic-strain)
            ;; (cl-mpm/utils:voigt-copy-into strain eps-e)
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
    ;; (setf stress (magicl:scale stress-u 1d0))
    (when (> damage 0.0d0)
      (let* ((j (magicl:det def))
             (degredation (expt (- 1d0 damage) 1d0))
             (min-strength 1d-6)
             )

        (let* ((min-damage (- 1d0 0d-6))
               (d1 (* min-damage (magicl:tref damage-tensor 0 0)))
               (d2 (* min-damage (magicl:tref damage-tensor 1 1)))
               (d3 (* min-damage (magicl:tref damage-tensor 2 2)))
               (d11 (- 1d0 d1))
               (d22 (- 1d0 d2))
               (d33 (- 1d0 d3))
               (d12 (sqrt (* d11 d22)))
               (d13 (sqrt (* d11 d33)))
               (d23 (sqrt (* d22 d33)))
               (g12 d12)
               (g13 d13)
               (g23 d23)
               (damage-effect (cl-mpm/utils:tensor-voigt-4th-from-list
                               (list
                                d11 d12 d13 0d0 0d0 0d0
                                d12 d22 d23 0d0 0d0 0d0
                                d13 d23 d33 0d0 0d0 0d0
                                0d0 0d0 0d0 g23 0d0 0d0
                                0d0 0d0 0d0 0d0 g13 0d0
                                0d0 0d0 0d0 0d0 0d0 g12
                                )))
               (new-de (magicl:.* damage-effect de))
               )
          ;; (pprint de)
          ;; (pprint new-de)
          ;; (pprint stress)
          (setf stress (magicl:@ new-de strain))
          ;; (pprint stress)
          ;; (break)
          )

        ;; (let* ((b (magicl:.-
        ;;            (magicl:from-diag (list 1d0 1d0 1d0) '(3 3))
        ;;           damage-tensor))
        ;;        (bsqrt (matrix-sqrt b)))
        ;;   )
        ;; (multiple-value-bind (l v)
        ;;     (cl-mpm/utils:eig damage-tensor)
        ;;   (let* ((d1 (- 1d0 (nth 0 l)))
        ;;          (d2 (- 1d0 (nth 1 l)))
        ;;          (d3 (- 1d0 (nth 2 l)))
        ;;          (m (magicl:from-diag
        ;;              (list d1
        ;;                    d2
        ;;                    d3
        ;;                    (sqrt (* d2 d3))
        ;;                    (sqrt (* d3 d1))
        ;;                    (sqrt (* d1 d2)))))
        ;;          (Q
        ;;            (magicl:transpose!
        ;;             (magicl:block-matrix (list
        ;;                                   (magicl:.* (magicl:column v 0)
        ;;                                              (magicl:column v 0))

        ;;                                   (magicl:.* (magicl:column v 1)
        ;;                                              (magicl:column v 1))
        ;;                                   (magicl:.* (magicl:column v 2)
        ;;                                              (magicl:column v 2))

        ;;                                   (magicl:scale!
        ;;                                    (magicl:.* (magicl:column v 0)
        ;;                                               (magicl:column v 1)) 2d0)
        ;;                                   (magicl:scale!
        ;;                                    (magicl:.* (magicl:column v 1)
        ;;                                               (magicl:column v 2)) 2d0)
        ;;                                   (magicl:scale!
        ;;                                    (magicl:.* (magicl:column v 2)
        ;;                                               (magicl:column v 0)) 2d0)

        ;;                                   (magicl:.* (magicl:column v 0)
        ;;                                              (cl-mpm/constitutive::rotate-vector (magicl:column v 0)))
        ;;                                   (magicl:.* (magicl:column v 1)
        ;;                                              (cl-mpm/constitutive::rotate-vector (magicl:column v 1)))
        ;;                                   (magicl:.* (magicl:column v 2)
        ;;                                              (cl-mpm/constitutive::rotate-vector (magicl:column v 2)))

        ;;                                   (magicl:scale! (magicl:.+
        ;;                                                   (magicl:.* (magicl:column v 0)
        ;;                                                              (cl-mpm/constitutive::rotate-vector (magicl:column v 1)))
        ;;                                                   (magicl:.* (magicl:column v 1)
        ;;                                                              (cl-mpm/constitutive::rotate-vector (magicl:column v 0)))) 1d0)

        ;;                                   (magicl:scale! (magicl:.+
        ;;                                                   (magicl:.* (magicl:column v 1)
        ;;                                                              (cl-mpm/constitutive::rotate-vector (magicl:column v 2)))
        ;;                                                   (magicl:.* (magicl:column v 2)
        ;;                                                              (cl-mpm/constitutive::rotate-vector (magicl:column v 1)))) 1d0)
        ;;                                   (magicl:scale! (magicl:.+
        ;;                                                   (magicl:.* (magicl:column v 2)
        ;;                                                              (cl-mpm/constitutive::rotate-vector (magicl:column v 0)))
        ;;                                                   (magicl:.* (magicl:column v 0)
        ;;                                                              (cl-mpm/constitutive::rotate-vector (magicl:column v 2)))) 1d0)
        ;;                                   ) '(2 6))))
        ;;          )
        ;;     (setf stress (magicl:@
        ;;                   (magicl:transpose Q)
        ;;                   m
        ;;                   Q
        ;;                   stress
        ;;                   )
        ;;           )

        ;;     ))
        
        ;; (magicl:scale! stress degredation)
        ;; (multiple-value-bind (l v) (cl-mpm/utils::eig
        ;;                             (magicl:scale! (voight-to-matrix stress) (/ 1d0 j)))
        ;;   (let* ()
        ;;     (loop for i from 0 to 2
        ;;           do
        ;;              (let* ((sii (nth i l))
        ;;                       (esii sii)
        ;;                     (vii (magicl::column v i))
        ;;                     )
        ;;                  (when t;(> esii 0d0)
        ;;                    ;;tensile damage -> unbounded
        ;;                    (let* (
        ;;                           (vsi (magicl:@ vii (magicl:transpose vii)))
        ;;                           (dii (magicl::trace (magicl:@ damage-tensor vsi)))
        ;;                           (degredation (expt (- 1d0 dii) 2d0)))
        ;;                      (setf (nth i l) (* sii (max 1d-12) degredation)))
        ;;                    ;; (setf (nth i l) (* sii (max 1d-8 (- 1d0 damage))))
        ;;                    )
        ;;                  )
        ;;           )
        ;;       (setf stress (magicl:scale! (matrix-to-voight (magicl:@ v
        ;;                                                               (magicl:from-diag l :type 'double-float)
        ;;                                                               (magicl:transpose v))) j))
        ;;       ))
        ))
    stress
    ))


(in-package :cl-mpm/damage)
(defmethod update-damage ((mp cl-mpm/particle::particle-chalk-brittle) dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (E cl-mpm/particle::mp-e)
                     (Gf cl-mpm/particle::mp-Gf)
                     (log-damage cl-mpm/particle::mp-log-damage)
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
                     ) mp
      (declare (double-float damage damage-inc critical-damage))
        (progn
          ;;Damage increment holds the delocalised driving factor
          (when (< damage 1d0)
            (setf ybar damage-inc)
            (setf k (max k ybar))
            (setf damage-inc 0d0)
            (let ((new-damage (max damage
                                   (damage-response-exponential k E init-stress ductility))))
              (setf damage-inc (- new-damage damage)))
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
              (setf damage-inc 0d0))))
  (values)
  ))
(defun plastic-damage-response-exponential (stress E length ductility)
  (declare (double-float stress E length ductility))
  "Function that controls how damage evolves with principal stresses"
  (let* ((ef (* ductility length E)))
    (- 1d0 (exp (- (/ stress ef))))))

(defun damage-response-linear (stress E Gf length init-stress ductility)
  (declare (double-float stress E Gf length init-stress ductility))
  "Linear softening"
  (let* ((ft init-stress)
         (e0 (/ ft E))
         (ef (* e0 ductility))
         (k (/ stress E)))
    (if (> k e0)
        (min 1d0
             (/ (max 0d0 (- k e0))
                (- ef e0)))
        0d0)))

(defun damage-response-linear-residual (stress E init-stress ductility residual)
  (declare (double-float stress E init-stress ductility))
  "Linear softening with plastic residual behaviour"
  (let* ((ft init-stress)
         (e0 (/ ft E))
         (ef (* e0 ductility))
         (k (/ stress E)))
    (if (> k e0)
        (min
         ;; Linear softening
         (min 1d0
              (/ (max 0d0 (- k e0))
                 (- ef e0)))
         (min 1d0
              (max 0d0
                   (- 1d0 (/ e0 k)))))
        0d0)))




(defmethod update-damage ((mp cl-mpm/particle::particle-chalk-delayed) dt)
  (when (cl-mpm/particle::mp-enable-damage mp)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (E cl-mpm/particle::mp-e)
                     (Gf cl-mpm/particle::mp-Gf)
                     (log-damage cl-mpm/particle::mp-log-damage)
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
                     (kc-r cl-mpm/particle::mp-k-compressive-residual-ratio)
                     (kt-r cl-mpm/particle::mp-k-tensile-residual-ratio)
                     (g-r cl-mpm/particle::mp-shear-residual-ratio)
                     (damage-tension cl-mpm/particle::mp-damage-tension)
                     (damage-shear cl-mpm/particle::mp-damage-shear)
                     (damage-compression cl-mpm/particle::mp-damage-compression)
                     )
        mp
      (declare (double-float damage damage-inc critical-damage k ybar tau dt))
      (when t;(<= damage 1d0)
        ;;Damage increment holds the delocalised driving factor
        (setf ybar damage-inc)
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
                 (damage-response-exponential k E init-stress ductility))))
          (declare (double-float new-damage))
          (setf damage-inc (- new-damage damage)))
        ;; (setf
        ;;  damage-tension (* kt-r damage)
        ;;  damage-compression (* kc-r damage)
        ;;  damage-shear (* g-r damage))
        (setf
         damage-tension (max damage-tension (damage-response-exponential-peerlings-residual k E init-stress ductility kt-r))
         damage-shear (max damage-shear (damage-response-exponential-peerlings-residual k E init-stress ductility g-r))
         damage-compression (max damage-compression (damage-response-exponential-peerlings-residual k E init-stress ductility kc-r)))

        (when (>= damage 1d0)
          (setf damage-inc 0d0))
        (incf (the double-float (cl-mpm/particle::mp-time-averaged-damage-inc mp)) (* damage-inc dt))
        (incf (the double-float (cl-mpm/particle::mp-time-averaged-ybar mp)) ybar)
        (incf (the double-float (cl-mpm/particle::mp-time-averaged-counter mp)))
        ;;Transform to log damage
        (incf damage damage-inc)
        ;;Transform to linear damage
        (setf damage (max 0d0 (min 1d0 damage)))
        (when (> damage critical-damage)
          (setf damage 1d0)
          (setf damage-inc 0d0))
        )
      (values))))
(defmethod update-damage ((mp cl-mpm/particle::particle-chalk-delayed-linear) dt)
  (when (cl-mpm/particle::mp-enable-damage mp)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (E cl-mpm/particle::mp-e)
                     (Gf cl-mpm/particle::mp-Gf)
                     (log-damage cl-mpm/particle::mp-log-damage)
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
                     (kc-r cl-mpm/particle::mp-k-compressive-residual-ratio)
                     (kt-r cl-mpm/particle::mp-k-tensile-residual-ratio)
                     (g-r cl-mpm/particle::mp-shear-residual-ratio)
                     (damage-tension cl-mpm/particle::mp-damage-tension)
                     (damage-shear cl-mpm/particle::mp-damage-shear)
                     (damage-compression cl-mpm/particle::mp-damage-compression)
                     ) mp
      (declare (double-float damage damage-inc critical-damage k ybar tau dt
                             damage-tension damage-shear damage-compression
                             kt-r kc-r g-r
                             init-stress ductility))
      (when (< damage 1d0)
        ;;Damage increment holds the delocalised driving factor
        (setf ybar damage-inc)
        (setf damage-inc 0d0)
        (let ((a tau-exp)
              (k0 init-stress))
          (incf k (the double-float
                       (*
                        dt
                        (/
                         (* k0
                            (the double-float
                                 (expt
                                  (/ (the double-float (max 0d0 (- ybar k)))
                                     k0) a)))
                         tau)))))
        (let ((new-damage
                (max
                 damage
                 ;; (damage-response-exponential k E init-stress ductility)
                 (damage-response-linear k E Gf (/ length (the double-float (sqrt 7d0))) init-stress ductility)
                 )))
          (declare (double-float new-damage))
          (setf damage-inc (- new-damage damage)))
        (setf
         damage-tension (max damage-tension (damage-response-linear-residual k E init-stress ductility kt-r))
         damage-shear (max damage-shear (damage-response-linear-residual k E init-stress ductility g-r))
         damage-compression (max damage-compression (damage-response-linear-residual k E init-stress ductility kc-r)))

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
      (values))))

(defmethod update-damage ((mp cl-mpm/particle::particle-chalk-anisotropic) dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (damage-tensor cl-mpm/particle::mp-damage-tensor)
                     (ybar-tensor cl-mpm/particle::mp-damage-ybar-tensor)
                     (E cl-mpm/particle::mp-e)
                     (Gf cl-mpm/particle::mp-Gf)
                     (log-damage cl-mpm/particle::mp-log-damage)
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
          (setf ybar damage-inc)
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
          ;; (setf ybar damage-inc)
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

(defmethod update-damage ((mp cl-mpm/particle::particle-chalk) dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (log-damage cl-mpm/particle::mp-log-damage)
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
          (setf ybar damage-inc)
          (when (< damage 1d0)
            (setf damage-inc (* dt
                                ;; (/ 1d0 (- 1d0 damage))
                                (damage-rate-profile-chalk damage-inc damage damage-rate init-stress))))
          (when (>= damage 1d0)
            ;; (setf damage-inc 0d0)
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
                     (log-damage cl-mpm/particle::mp-log-damage)
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
          (setf ybar damage-inc)
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
          (setf damage-increment (tensile-energy-norm strain E de))
          ;;            (angle )
          ;(setf damage-increment (criterion-mc strain (* angle (/ pi 180d0)) E nu))
          ;(setf damage-increment (criterion-dp strain (* angle (/ pi 180d0)) E nu))
          ;; (setf damage-increment
          ;;       (max 0d0
          ;;            (drucker-prager-criterion
          ;;             (magicl:scale stress (/ 1d0 (magicl:det def))) (* angle (/ pi 180d0)))))
          ;; (let* ((strain+
          ;;          (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix strain))
          ;;            (loop for i from 0 to 2
          ;;                  do
          ;;                     (setf (nth i l) (max (nth i l) 0d0)))
          ;;            (cl-mpm/utils:matrix-to-voigt (magicl:@ v
          ;;                                                     (magicl:from-diag l :type 'double-float)
          ;;                                                     (magicl:transpose v)))))
          ;;        (strain- (magicl:.- strain strain+))
          ;;        ;(invar (cl-mpm/utils:voigt-from-list (list 1d0 1d0 1d0 2d0 2d0 2d0)))
          ;;        (e+ (sqrt (max 0d0 (* E (cl-mpm/fastmaths::dot strain+ (magicl:@ de strain+))))))
          ;;        (e- (sqrt (max 0d0 (* E (cl-mpm/fastmaths::dot strain- (magicl:@ de strain-))))))
          ;;        (k (/ fc ft)))
          ;;   ;; (format t "Energy real ~A~%" (magicl:@ de strain+))
          ;;   (setf damage-increment
          ;;         ;; (+ e+ (/ e- k))
          ;;         e+
          ;;         ;; (/
          ;;         ;;  (+ (* k e+) e-)
          ;;         ;;  (+ k 1d0)
          ;;         ;;  )
          ;;         ))
          
          ;; (setf damage-increment (smooth-rankine-criterion (magicl:scale stress (/ 1d0 (magicl:det def)))))

          ;; (let ((strain (magicl:.*
          ;;                strain
          ;;                (cl-mpm/utils:voigt-from-list (list 1d0 1d0 1d0 0.5d0 0.5d0 0.5d0)))))
          ;;   (multiple-value-bind (e_1 e_2 e_3) (principal-stresses-3d strain)
          ;;     (let* ((K (/ E (* 3d0 (- 1d0 (* 2d0 nu)))))
          ;;            (mu (* (/ E (* 2d0 (+ 1d0 nu)))))
          ;;            (angle (* angle (/ pi 180d0)))
          ;;            (B (/ (* 2d0 (sin angle))
          ;;                  (* (sqrt 3) (+ 3d0 (sin angle)))))
          ;;            (j2
          ;;              ;; (cl-mpm/constitutive::voigt-j2 (cl-mpm/utils::deviatoric-voigt strain))
          ;;              (* 1/6 (+
          ;;                      (expt (- e_1 e_2) 2)
          ;;                      (expt (- e_2 e_3) 2)
          ;;                      (expt (- e_3 e_1) 2)))
          ;;              )
          ;;            (i1 (+ e_1 e_2 e_3))
          ;;            (idisc (* -6d0 B (sqrt j2)))
          ;;            (jdisc (* 2d0 mu (sqrt j2)))
          ;;            (psi1 (+
          ;;                   (* 0.5d0 K (expt i1 2))
          ;;                   (* 2d0 mu j2)
          ;;                   ))
          ;;            (psi2 (* (/ 1d0 (+ (* 18d0 (expt B 2) K)
          ;;                               (* 2d0 mu)))
          ;;                     (expt
          ;;                      (- (* 2d0 mu (sqrt j2))
          ;;                         (* 3d0 B K i1)) 2)))
          ;;            )
          ;;       (setf damage-increment
          ;;             (sqrt
          ;;              (max 0d0
          ;;                   (* E
          ;;                      (cond
          ;;                        ((< idisc i1) psi1)
          ;;                        ((and (>= idisc i1)
          ;;                              (>= jdisc (* 3d0 B K i1))
          ;;                              ) psi2)
          ;;                        (t 0d0)
          ;;                        ))))
          ;;             ))
          ;;     ))

          ;; (multiple-value-bind (e_1 e_2 e_3) (principal-stresses-3d (magicl:.*
          ;;                                                            strain
          ;;                                                            (cl-mpm/utils:voigt-from-list (list 1d0 1d0 1d0 0.5d0 0.5d0 0.5d0))))
          ;;   (let* ((e1+ (max 0d0 e_1))
          ;;          (e2+ (max 0d0 e_2))
          ;;          (e3+ (max 0d0 e_3))
          ;;          (e1- (min 0d0 e_1))
          ;;          (e2- (min 0d0 e_2))
          ;;          (e3- (min 0d0 e_3))
          ;;          (n 3d0)
          ;;          (mu (/ E (* 2d0 (+ 1d0 nu))))
          ;;          (K (/ E (* 3d0 (- 1d0 (* 2d0 nu)))))
          ;;          (lam (/ (* nu E) (* (+ 1d0 nu) (- 1d0 (* 2d0 nu)))))
          ;;          (ev+ (* 0.5d0 K
          ;;                  (/ 1d0 n)
          ;;                  (expt
          ;;                   (max
          ;;                    0d0
          ;;                    (+ e_1
          ;;                       e_2
          ;;                       e_3)) 2d0)))
          ;;          (ed+
          ;;            (-
          ;;             (* mu (/ (- n 1d0) n) (+ (expt e1+ 2)
          ;;                                      (expt e2+ 2)
          ;;                                      (expt e3+ 2)))
          ;;             (* mu (/ 2 n)
          ;;                (+
          ;;                 (* e1+ e2-)
          ;;                 (* e2+ e1-)
          ;;                 (* e1+ e3-)
          ;;                 (* e3+ e1-)
          ;;                 (* e3+ e2-)
          ;;                 (* e2+ e3-))))
          ;;          )
          ;;          (ed-
          ;;            (* mu (/ (- n 1d0)
          ;;                     n)
          ;;               (+ (expt e1- 2)
          ;;                  (expt e2- 2)
          ;;                  (expt e3- 2))))
          ;;          (angle (* angle (/ pi 180d0)))
          ;;          (f (+
          ;;              (* (/ (+ lam mu)
          ;;                    mu)
          ;;                 (/ (+ e_1 e_3)
          ;;                    (- e_1 e_3))
          ;;                 (sin angle)
          ;;                 )
          ;;              (/ (* -1d0 ft (cos angle))
          ;;                 (* mu (- e_1 e_3)))
          ;;              1d0))
          ;;          )
          ;;     (setf damage-increment (sqrt (max 0d0 (* E
          ;;                                              (+ ev+ ed+
          ;;                                                 ;; (if (> f 0d0)
          ;;                                                 ;;     ed-
          ;;                                                 ;;     0d0)
          ;;                                                 ))))))
          ;;   )

          ;; (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d stress)
          ;;   (multiple-value-bind (e_1 e_2 e_3) (principal-stresses-3d strain)
          ;;     (let* ((pressure-effective (* 1d0 damage pressure))
          ;;            (j2 (cl-mpm/constitutive::voigt-j2
          ;;                 (cl-mpm/utils::deviatoric-voigt cauchy-undamaged)))
          ;;            (p (+ s_1 s_2 s_3))
          ;;            (p (cl-mpm/utils::trace-voigt stress))
          ;;            (p (if (> p 0d0)
          ;;                   (* (- 1d0 kt-r) p)
          ;;                   (* (- 1d0 kc-r) p)))
          ;;            (et (sqrt (max 0d0
          ;;                           (* E
          ;;                              (+
          ;;                               (* (max s_1 0d0) e_1)
          ;;                               (* (max s_2 0d0) e_2)
          ;;                               (* (max s_3 0d0) e_3)
          ;;                               )))))
          ;;            (ec (sqrt (max 0d0
          ;;                           (* E
          ;;                              (+
          ;;                               (* (min s_1 0d0) e_1)
          ;;                               (* (min s_2 0d0) e_2)
          ;;                               (* (min s_3 0d0) e_3)
          ;;                               )))))
          ;;            (k (/ fc ft))
          ;;            ;; (s_1 (max 0d0 s_1))

          ;;            (s_1
          ;;              et
          ;;              ;; (/
          ;;              ;;  (+ (* k et) ec)
          ;;              ;;  (+ k 1d0)
          ;;              ;;  )
          ;;              )

          ;;            ;; (s_1
          ;;            ;;   (sqrt
          ;;            ;;    (* E
          ;;            ;;       (+
          ;;            ;;        (*
          ;;            ;;         (/ 1d0 6d0)
          ;;            ;;         (max 0d0 (cl-mpm/utils::trace-voigt stress))
          ;;            ;;         (cl-mpm/utils::trace-voigt strain))
          ;;            ;;        (* 0.5d0
          ;;            ;;           (cl-mpm/fastmaths:dot (cl-mpm/utils::deviatoric-voigt stress)
          ;;            ;;                                (cl-mpm/utils::deviatoric-voigt strain)))))))
          ;;            ;; (p 0d0)
          ;;            ;; (s_1 (sqrt (* E
          ;;            ;;               (+
          ;;            ;;                (* p (cl-mpm/utils::trace-voigt strain))
          ;;            ;;                (* (- 1d0 g-r)
          ;;            ;;                   (cl-mpm/fastmaths:dot (cl-mpm/utils::deviatoric-voigt stress)
          ;;            ;;                                        (cl-mpm/utils::deviatoric-voigt strain)))))))
          ;;            ;; (s_1
          ;;            ;;   (sqrt (max 0d0
          ;;            ;;              (* E
          ;;            ;;                 (cl-mpm/fastmaths:dot (cl-mpm/utils::deviatoric-voigt stress)
          ;;            ;;                                      (cl-mpm/utils::deviatoric-voigt strain))))))
          ;;            ;; (s_1 (- s_1 pressure-effective))
          ;;            ;; (s_2 (- s_2 pressure-effective))
          ;;            ;; (s_3 (- s_3 pressure-effective))
          ;;            ;; (s_1 (max 0d0 s_1))
          ;;            ;; (s_2 (max 0d0 s_2))
          ;;            ;; (s_3 (max 0d0 s_3))
          ;;            ;; (s_1 (* e
          ;;            ;;         (sqrt
          ;;            ;;          (+ (expt s_1 2)
          ;;            ;;             (expt s_2 2)
          ;;            ;;             (expt s_3 2)))))
          ;;            ;; (k (/ fc ft))
          ;;            ;; (s_1 (* e
          ;;            ;;         (/
          ;;            ;;          (+
          ;;            ;;             (* k
          ;;            ;;                (sqrt
          ;;            ;;                 (+ (expt (max 0d0 s_1) 2)
          ;;            ;;                    (expt (max 0d0 s_2) 2)
          ;;            ;;                    (expt (max 0d0 s_3) 2))))
          ;;            ;;             (sqrt
          ;;            ;;              (+ (expt (max 0d0 (- s_1)) 2)
          ;;            ;;                 (expt (max 0d0 (- s_2)) 2)
          ;;            ;;                 (expt (max 0d0 (- s_3)) 2)))
          ;;            ;;             )
          ;;            ;;          (+ 1d0 k)
          ;;            ;;          ))
          ;;            ;;      )
          ;;            (angle (* angle (/ pi 180d0)))
          ;;            ;; (ft 200d3)
          ;;            ;; (fc 600d3)
          ;;            ;; (fc 500d3)
          ;;            ;; (angle (atan (* 3 (/ (- fc ft) (+ fc ft)))))
          ;;            ;; (fc (*
          ;;            ;;      ft
          ;;            ;;      -3d0
          ;;            ;;      (/ (+ 1d0 (tan angle))
          ;;            ;;         (- 1d0 (tan angle)))))
          ;;            ;; (k (* (/ 2d0 (sqrt 3)) (/ (* ft fc) (+ ft fc))))
          ;;            ;; ;; (c 1d4)
          ;;            ;;another dp

          ;;            ;;smooth rankine
          ;;            ;; (s_1 (sqrt
          ;;            ;;       (+ (expt (max 0d0 s_1) 2)
          ;;            ;;          (expt (max 0d0 s_2) 2)
          ;;            ;;          (expt (max 0d0 s_3) 2))))

          ;;            ;;good drucker-prager
          ;;            ;; (s_1 (* (/ 3d0 (+ 3 (tan angle)))
          ;;            ;;         (+ (sqrt (* 3 j2)) (* 1/3 (tan angle) p))))


          ;;            ;; (k (/ fc ft))
          ;;            ;; (i1 (+ s_1 s_2 s_3))
          ;;            ;; (k-factor (/ (- k 1d0)
          ;;            ;;              (- 1d0 (* 2d0 nu))))
          ;;            ;; (s_1 (* e
          ;;            ;;         (+ (* i1 (/ k-factor (* 2d0 k)))
          ;;            ;;             (* (/ 1d0 (* 2d0 k))
          ;;            ;;                (sqrt (+ (expt (* k-factor i1) 2)
          ;;            ;;                         (* (/ (* 12 k) (expt (- 1d0 nu) 2)) j2)
          ;;            ;;                         ))))))


          ;;            ;; (s_1
          ;;            ;;   (* 0.5d0
          ;;            ;;      (+
          ;;            ;;       (* (/ 1d0 (* k init-stress))
          ;;            ;;          (+
          ;;            ;;           (expt (- s_1 s_2) 2)
          ;;            ;;           (expt (- s_2 s_3) 2)
          ;;            ;;           (expt (- s_3 s_1) 2))
          ;;            ;;          )
          ;;            ;;       (* (/ 2d0 (* k init-stress)) (- (* k init-stress) init-stress)
          ;;            ;;          (+ s_1 s_2 s_3))
          ;;            ;;       )))
          ;;            )
          ;;       ;; (setf damage-increment s_1)
          ;;       (when (> s_1 0d0)
          ;;         ;; (setf damage-increment (* s_1 (expt (- 1d0 damage) -2d0)))
          ;;         (setf damage-increment s_1)
          ;;         )
          ;;       )))
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
                     (log-damage cl-mpm/particle::mp-log-damage)
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
          (setf ybar damage-inc)
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
