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
  ((trial-elastic-strain
    :accessor mp-trial-strain
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils:voigt-zeros))
   (material-damping-factor
    :accessor mp-material-damping-factor
    :initarg :material-damping
    :initform 0d0)
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
  ((damage-tensor
    :accessor mp-damage-tensor
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils::matrix-zeros)
    )
   (damage-ybar-tensor
    :accessor mp-damage-ybar-tensor
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils::matrix-zeros))
   )
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
                   (trial-elastic-strain mp-trial-strain)
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
    (cl-mpm/utils:voigt-copy-into strain trial-elastic-strain)
    (when enable-plasticity
        (progn
          (multiple-value-bind (sig eps-e f inc)
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
             yield-func f)
            ;; (cl-mpm/fastmaths:fast-.+ plastic-strain
            ;;                           (cl-mpm/fastmaths:fast-.- trial-elastic-strain strain)
            ;;                           plastic-strain)
            (let ()
              (setf ps-vm (+ ps-vm-1 inc))
              (setf ps-vm-inc inc)))))
    (cl-mpm/utils:voigt-copy-into stress-u stress)
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
            (;cl-mpm/constitutive::mc-plastic
             cl-mpm/ext::constitutive-mohr-coulomb
             stress-u de strain
                                             E
                                             nu
                                             phi
                                             psi
                                             coheasion)
            (setf stress sig
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
                     (damage-inc cl-mpm/particle::mp-damage-increment)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (ybar-prev cl-mpm/particle::mp-damage-ybar-prev)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (pressure cl-mpm/particle::mp-pressure)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     ;; (length-t cl-mpm/particle::mp-true-local-length)
                     (length cl-mpm/particle::mp-local-length)
                     (k cl-mpm/particle::mp-history-stress)
                     (k-n cl-mpm/particle::mp-history-stress-n)
                     (ductility cl-mpm/particle::mp-ductility)
                     (kc-r cl-mpm/particle::mp-k-compressive-residual-ratio)
                     (kt-r cl-mpm/particle::mp-k-tensile-residual-ratio)
                     (g-r cl-mpm/particle::mp-shear-residual-ratio)
                     (damage-tension cl-mpm/particle::mp-damage-tension)
                     (damage-shear cl-mpm/particle::mp-damage-shear)
                     (damage-compression cl-mpm/particle::mp-damage-compression)
                     (peerlings cl-mpm/particle::mp-peerlings-damage)
                     ) mp
      (declare (double-float damage damage-inc critical-damage
                             damage-tension damage-shear damage-compression
                             kt-r kc-r g-r
                             ))
        (progn
          ;;Damage increment holds the delocalised driving factor
          (when t;(< damage 1d0)
            (setf k (max k-n ybar-prev ybar))
            (let ((new-damage (max damage
                                   (damage-response-exponential k E init-stress ductility))))
              (setf damage-inc (- new-damage damage)))
            (if peerlings
                (setf
                 damage-tension (max damage-tension (damage-response-exponential-peerlings-residual k E init-stress ductility kt-r))
                 damage-shear (max damage-shear (damage-response-exponential-peerlings-residual k E init-stress ductility g-r))
                 damage-compression (max damage-compression (damage-response-exponential-peerlings-residual k E init-stress ductility kc-r)))
                (setf
                 damage-tension (* kt-r damage)
                 damage-compression (* kc-r damage)
                 damage-shear (* g-r damage)))
            ;; (when (>= damage 1d0)
            ;;   (setf damage-inc 0d0)
            ;;   (setf ybar 0d0))
            (incf (cl-mpm/particle::mp-time-averaged-damage-inc mp) damage-inc)
            (incf (cl-mpm/particle::mp-time-averaged-ybar mp) ybar)
            (incf (cl-mpm/particle::mp-time-averaged-counter mp))
            ;;Transform to log damage
            (incf damage damage-inc)
            ;;Transform to linear damage
            (setf damage (max 0d0 (min 1d0 damage)))
            ;; (when (> damage critical-damage)
            ;;   (setf damage 1d0)
            ;;   (setf damage-inc 0d0))
            ))
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



(defun apply-vol-degredation (mp dt)
  (with-accessors ((damage        cl-mpm/particle::mp-damage)
                   (damage-t      cl-mpm/particle::mp-damage-tension)
                   (damage-c      cl-mpm/particle::mp-damage-compression)
                   (damage-s      cl-mpm/particle::mp-damage-shear)
                   (stress        cl-mpm/particle::mp-stress)
                   (p-mod         cl-mpm/particle::mp-p-modulus)
                   (E cl-mpm/particle::mp-e)
                   (nu cl-mpm/particle::mp-nu)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (enable-damage cl-mpm/particle::mp-enable-damage))
      mp
    (declare (double-float damage damage-t damage-c damage-s E nu))
    (when (and
           enable-damage
           (> damage 0.0d0))
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
               stress))
        (let ((K (/ e (* 3 (- 1d0 (* 2 nu)))))
              (G (/ e (* 2 (+ 1d0 nu)))))
          (declare (double-float K G p-mod))
          (setf K
                (if (> p 0d0)
                    (* (- 1d0 damage-t) K)
                    (* (- 1d0 damage-c) K)))
          (setf G (* G (- 1d0 damage-s)))
          ;; (when (> (cl-mpm/particle::mp-yield-func mp) 0d0)
          ;;   (setf G 0d0))
          (when (> (cl-mpm/particle::mp-yield-func mp) 0d0)
            ;; (setf G 0d0)
            (setf K (* K (cos (cl-mpm/particle::mp-phi mp))))
            (setf G (* G (sin (cl-mpm/particle::mp-phi mp))))
            )
          (setf p-mod (* (expt (cl-mpm/fastmaths::det def) -2) (+ K (* 4/3 G)))))))))

(defun apply-tensile-vol-degredation (mp dt)
  (with-accessors ((damage        cl-mpm/particle::mp-damage)
                   (damage-t      cl-mpm/particle::mp-damage-tension)
                   (damage-c      cl-mpm/particle::mp-damage-compression)
                   (damage-s      cl-mpm/particle::mp-damage-shear)
                   (stress        cl-mpm/particle::mp-stress)
                   (enable-damage cl-mpm/particle::mp-enable-damage))
      mp
    (declare (double-float damage damage-t damage-c damage-s))
    (when (and
           enable-damage
           (> damage 0.0d0))
      (let ((p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
            (s (cl-mpm/constitutive::deviatoric-voigt stress)))
        (declare (double-float damage-t damage-c damage-s))
        (setf p
              (if (> p 0d0)
                  (* (- 1d0 damage-t) p)
                  ;; (* (- 1d0 damage-c) p)
                  p
                  ))
        (setf stress
              (cl-mpm/fastmaths:fast-.+
               (cl-mpm/constitutive::voight-eye p)
               (cl-mpm/fastmaths:fast-scale! s (- 1d0 damage-s))
               stress))))))

(defun apply-vol-pressure-degredation (mp dt pressure)
  (with-accessors ((damage        cl-mpm/particle::mp-damage)
                   (damage-t      cl-mpm/particle::mp-damage-tension)
                   (damage-c      cl-mpm/particle::mp-damage-compression)
                   (damage-s      cl-mpm/particle::mp-damage-shear)
                   (stress        cl-mpm/particle::mp-stress)
                   (enable-damage cl-mpm/particle::mp-enable-damage))
      mp
    (declare (double-float damage damage-t damage-c damage-s))
    (when (and
           enable-damage
           (> damage 0.0d0))
      (let ((p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
            (s (cl-mpm/constitutive::deviatoric-voigt stress)))
        (declare (double-float damage-t damage-c damage-s))
        (setf p
              (if (> p 0d0)
                  (* (- 1d0 damage-t) p)
                  (* (- 1d0 damage-c) p)))
        (setf stress
              (cl-mpm/fastmaths:fast-.+
               (cl-mpm/constitutive::voight-eye (- p pressure))
               (cl-mpm/fastmaths:fast-scale! s (- 1d0 damage-s))
               stress))))))

(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-chalk-brittle) dt)
  ;; (apply-tensile-vol-degredation mp dt)
  ;(apply-vol-degredation mp dt)
  (with-accessors ((p cl-mpm/particle::mp-pressure)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (stress cl-mpm/particle::mp-stress)
                   (damage cl-mpm/particle::mp-damage)
                   (enable-damage cl-mpm/particle::mp-enable-damage)
                   (p-mod cl-mpm/particle::mp-p-modulus)
                   )
      mp
    ;; (apply-tensile-vol-degredation mp dt)
    (setf p-mod (* (expt (cl-mpm/fastmaths::det def) -2) (cl-mpm/particle::compute-p-modulus mp)))
    (when enable-damage
      (apply-vol-degredation mp dt)
      ;(apply-vol-pressure-degredation mp dt (* 1d0 (magicl:det def) (/ p 3)))
      ;; (apply-vol-pressure-degredation mp dt (* -1d0 (/ 1d0 (magicl:det def)) (/ p 1)))
      ;; (setf stress (cl-mpm/constitutive::voight-eye p))
      )
      ;; (cl-mpm/fastmaths:fast-.+ stress
      ;;                           (cl-mpm/constitutive::voight-eye (* ;; (magicl:det def)
      ;;                                                               (/ p 3) damage)) stress)
    ))

(defun delay-damage-update (k-0 y-0 y-1))

(defmethod compute-damage ((mp cl-mpm/particle::particle-chalk-delayed))
  (with-accessors ((k cl-mpm/particle::mp-history-stress)
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
    (declare (double-float kt-r kc-r g-r damage))
    (setf
     damage
     (damage-response-exponential k E init-stress ductility))
    (if peerlings
        (setf
         damage-tension (damage-response-exponential-peerlings-residual k E init-stress ductility kt-r)
         damage-shear (damage-response-exponential-peerlings-residual k E init-stress ductility g-r)
         damage-compression (damage-response-exponential-peerlings-residual k E init-stress ductility kc-r))
        (setf
         damage-tension (* kt-r damage)
         damage-compression (* kc-r damage)
         damage-shear (* g-r damage)))))

(defmethod update-damage ((mp cl-mpm/particle::particle-chalk-delayed) dt)
  (when (cl-mpm/particle::mp-enable-damage mp)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (damage-n cl-mpm/particle::mp-damage-n)
                     (E cl-mpm/particle::mp-e)
                     (Gf cl-mpm/particle::mp-Gf)
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
      (declare (double-float damage damage-inc critical-damage k ybar tau dt))
      (when t;(<= damage 1d0)
        ;;Damage increment holds the delocalised driving factor
        (setf damage-inc 0d0)
        (let ((a tau-exp)
              (k0 init-stress))
          (when (or (>= ybar-prev k0)
                    (>= ybar k0))
            (setf
             k
             (cl-mpm/damage::huen-integration
              k-n
              ybar-prev
              ybar
              k0
              tau
              tau-exp
              dt
              ))
            (let ((new-damage
                    (max
                     damage-n
                     (damage-response-exponential k E init-stress ductility))))
              (declare (double-float new-damage))
              (setf damage new-damage)
              (setf damage-inc (- new-damage damage-n))))
          (if peerlings
              (setf
               damage-tension (damage-response-exponential-peerlings-residual k E init-stress ductility kt-r)
               damage-shear (damage-response-exponential-peerlings-residual k E init-stress ductility g-r)
               damage-compression (damage-response-exponential-peerlings-residual k E init-stress ductility kc-r))
              (setf
               damage-tension (* kt-r damage)
               damage-compression (* kc-r damage)
               damage-shear (* g-r damage))))

        ;; (when (>= damage 1d0)
        ;;   (setf damage-inc 0d0))
        (incf (the double-float (cl-mpm/particle::mp-time-averaged-damage-inc mp)) (* damage-inc dt))
        (incf (the double-float (cl-mpm/particle::mp-time-averaged-ybar mp)) ybar)
        (incf (the double-float (cl-mpm/particle::mp-time-averaged-counter mp)))
        ;;Transform to log damage
        ;; (incf damage damage-inc)
        ;;Transform to linear damage
        (setf damage (max 0d0 (min 1d0 damage)))
        ;; (when (> damage critical-damage)
        ;;   (setf damage 1d0)
        ;;   (setf damage-inc 0d0))
        )
      (values))))
(defmethod update-damage ((mp cl-mpm/particle::particle-chalk-delayed-linear) dt)
  (when (cl-mpm/particle::mp-enable-damage mp)
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

(defmethod update-damage ((mp cl-mpm/particle::particle-chalk) dt)
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
