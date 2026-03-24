(in-package :cl-mpm/damage)

(defun apply-isotropic-degredation (mp)
  (with-accessors ((damage        cl-mpm/particle::mp-damage)
                   (undamaged-stress        cl-mpm/particle::mp-undamaged-stress)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (stress        cl-mpm/particle::mp-stress)
                   (enable-damage cl-mpm/particle::mp-enable-damage)
                   (p-mod cl-mpm/particle::mp-p-modulus-0)
                   (e cl-mpm/particle::mp-e)
                   (nu cl-mpm/particle::mp-nu))
      mp
    (declare (double-float damage))
    (when (and
           enable-damage
           (> damage 0.0d0))
      (cl-mpm/utils:voigt-copy-into undamaged-stress stress)
      (cl-mpm/fastmaths:fast-scale! stress (/ (- 1d0 damage) (cl-mpm/fastmaths::det-3x3 def)))
      (setf (cl-mpm/particle::mp-p-modulus mp)
            (*
             (- 1d0 damage)
             (cl-mpm/particle::compute-p-modulus mp))))))

(defun apply-tensile-stress-degredation (mp)
  (with-accessors ((damage cl-mpm/particle::mp-damage)
                   (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                   (stress cl-mpm/particle::mp-stress)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (enable-damage cl-mpm/particle::mp-enable-damage)
                   (p-mod cl-mpm/particle::mp-p-modulus-0)
                   (e cl-mpm/particle::mp-e)
                   (nu cl-mpm/particle::mp-nu)
                   )
      mp
    (let ()
      (when (> damage 0.0d0)
        (cl-mpm/utils::copy-into undamaged-stress stress)
        (cl-mpm/fastmaths:fast-scale! stress (/ 1d0 (cl-mpm/fastmaths:det-3x3 def)))
        (multiple-value-bind (l v) (cl-mpm/utils::eig
                                    (voight-to-matrix
                                     stress))

          (let* ((degredation (expt (- 1d0 damage) 1d0)))
            (loop for i from 0 to 2
                  do (let* ((sii (nth i l)))
                       (when (> sii 0d0)
                         ;;Tensile damage -> unbounded
                         (setf (nth i l) (* sii degredation)))))
            (setf stress
                  (cl-mpm/utils:matrix-to-voight
                   (magicl:@
                    v
                    (cl-mpm/utils::matrix-diag l)
                    (magicl:transpose v))
                   stress))))))))

(defun apply-tensile-strain-degredation (mp)
  (with-accessors ((damage        cl-mpm/particle::mp-damage)
                   (stress        cl-mpm/particle::mp-stress)
                   (strain        cl-mpm/particle::mp-strain)
                   (enable-damage cl-mpm/particle::mp-enable-damage)
                   (p-mod cl-mpm/particle::mp-p-modulus-0)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (e cl-mpm/particle::mp-e)
                   (nu cl-mpm/particle::mp-nu)
                   (de cl-mpm/particle::mp-elastic-matrix))
      mp
    (when (> damage 0.0d0)
      (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils::voigt-to-matrix strain))
        (let* ()
          ;; (pprint l)
          (loop for i from 0 below (length l)
                do (let* ((sii (nth i l)))
                     ;;Tensile damage -> unbounded
                     (when (> sii 0d0)
                       (setf (nth i l) (* (nth i l)
                                          (- 1d0 damage))))))
          ;; (pprint strain)
          (let ((strain+ (cl-mpm/utils:matrix-to-voigt
                          (magicl:@
                           v
                           (magicl:from-diag l :type 'double-float)
                           (magicl:transpose v)))))
            ;; (pprint strain+)
            (setf stress
                  (cl-mpm/constitutive::linear-elastic-mat
                   (cl-mpm/utils:matrix-to-voigt
                    (magicl:@
                     v
                     (cl-mpm/utils::matrix-diag l)
                     (magicl:transpose v)))
                   de
                   stress))
            (cl-mpm/fastmaths:fast-scale! stress (* 1d0 (cl-mpm/fastmaths:det-3x3 def)))))))))


(defun apply-vol-degredation (mp)
  (with-accessors ((damage        cl-mpm/particle::mp-damage)
                   (stress        cl-mpm/particle::mp-stress)
                   (enable-damage cl-mpm/particle::mp-enable-damage)
                   (p-mod cl-mpm/particle::mp-p-modulus-0)
                   (e cl-mpm/particle::mp-e)
                   (nu cl-mpm/particle::mp-nu)
                   )
      mp
    (declare (double-float damage))
    (when (and
           enable-damage
           (> damage 0.0d0))
      (let ((p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
            (s (cl-mpm/constitutive::deviatoric-voigt stress)))
        (setf stress
              (cl-mpm/fastmaths:fast-.+
               (cl-mpm/constitutive::voight-eye p)
               (cl-mpm/fastmaths:fast-scale! s (- 1d0 damage))
               stress))
        (let ((K (/ e (* 3 (- 1d0 (* 2 nu)))))
              (G (* (- 1d0 damage) (/ e (* 2 (+ 1d0 nu))))))
          (setf p-mod (+ K (* G (/ 4 3)))))
        ))))

(defun apply-eigen-degredation (mp)
  
  )

(defun apply-vol-pressure-degredation (mp dt pressure)
  (with-accessors ((damage        cl-mpm/particle::mp-damage)
                   (damage-t      cl-mpm/particle::mp-damage-tension)
                   (damage-c      cl-mpm/particle::mp-damage-compression)
                   (damage-s      cl-mpm/particle::mp-damage-shear)
                   (stress        cl-mpm/particle::mp-stress)
                   (p-mod         cl-mpm/particle::mp-p-modulus-0)
                   (E cl-mpm/particle::mp-e)
                   (nu cl-mpm/particle::mp-nu)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (enable-damage cl-mpm/particle::mp-enable-damage))
      mp
    (declare (double-float damage damage-t damage-c damage-s))
    (when (and
           enable-damage
           (> damage 0.0d0))
      (let* ((exponent 1)
            (p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
            (pind (- p pressure))
            ;; (pind p)
            (s (cl-mpm/constitutive::deviatoric-voigt stress)))
        (declare (double-float damage-t damage-c damage-s))
        (setf p
              (if (> pind 0d0)
                  (* (- 1d0 (expt damage-t exponent)) p)
                  (* (- 1d0 (expt damage-c exponent)) p)))
        ;; (pprint pressure)
        ;; (break)
        (setf stress
              (cl-mpm/fastmaths:fast-.+
               (cl-mpm/constitutive::voight-eye (- p pressure))
               (cl-mpm/fastmaths:fast-scale! s (- 1d0 (expt damage-s exponent)))
               stress))
        (let* ((K (/ e (* 3 (- 1d0 (* 2 nu)))))
               (G (/ e (* 2 (+ 1d0 nu))))
               (P-0 (+ K (* 4/3 G))))
          (setf K
                (if (> pind 0d0)
                    (* (- 1d0 (expt damage-t exponent)) K)
                    (* (- 1d0 (expt damage-c exponent)) K)))
          (setf G (* G (- 1d0 (expt damage-s exponent))))

          (when (> (cl-mpm/particle::mp-yield-func mp) 0d0)
            ;; (setf G 0d0)
            (setf K (* K (cos (cl-mpm/particle::mp-phi mp))))
            (setf G (* G (sin (cl-mpm/particle::mp-phi mp))))
            )
          (setf p-mod (max (* P-0 1d-9) (+ K (* 4/3 G))))
          )
        ))))

(defun apply-tcs-degredation (mp)
  (with-accessors ((damage        cl-mpm/particle::mp-damage)
                   (damage-t      cl-mpm/particle::mp-damage-tension)
                   (damage-c      cl-mpm/particle::mp-damage-compression)
                   (damage-s      cl-mpm/particle::mp-damage-shear)
                   (stress-u      cl-mpm/particle::mp-undamaged-stress)
                   (stress        cl-mpm/particle::mp-stress)
                   (p-mod         cl-mpm/particle::mp-p-modulus-0)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (E cl-mpm/particle::mp-e)
                   (nu cl-mpm/particle::mp-nu)
                   (enable-damage cl-mpm/particle::mp-enable-damage))
      mp
    (declare (double-float damage damage-t damage-c damage-s))
    (when (and
           enable-damage
           (> damage 0.0d0))
      (let* ((exponent 1)
             (stress-uc (cl-mpm/fastmaths:fast-scale-voigt stress-u (/ 1d0 (cl-mpm/fastmaths:det-3x3 def))))
             (p (/ (cl-mpm/constitutive::voight-trace stress-uc) 3d0))
             (pind p)
             (s (cl-mpm/constitutive::deviatoric-voigt stress-uc)))
        (declare (double-float damage-t damage-c damage-s))
        (setf p
              (if (> pind 0d0)
                  (* (- 1d0 (expt damage-t exponent)) p)
                  (* (- 1d0 (expt damage-c exponent)) p)))
        ;; (pprint pressure)
        ;; (break)
        (setf stress
              (cl-mpm/fastmaths:fast-.+
               (cl-mpm/constitutive::voight-eye p)
               (cl-mpm/fastmaths:fast-scale! s (- 1d0 (expt damage-s exponent)))
               stress))
        (let* ((K (/ e (* 3 (- 1d0 (* 2 nu)))))
               (G (/ e (* 2 (+ 1d0 nu))))
               (P-0 (+ K (* 4/3 G))))
          (setf K
                (if (> pind 0d0)
                    (* (- 1d0 (expt damage-t exponent)) K)
                    (* (- 1d0 (expt damage-c exponent)) K)))
          (setf G (* G (- 1d0 (expt damage-s exponent))))
          (setf p-mod (max (* P-0 1d-9) (+ K (* 4/3 G)))))))))

(defun apply-tcs-strain-degredation (mp)
  (with-accessors ((damage        cl-mpm/particle::mp-damage)
                   (damage-t      cl-mpm/particle::mp-damage-tension)
                   (damage-c      cl-mpm/particle::mp-damage-compression)
                   (damage-s      cl-mpm/particle::mp-damage-shear)
                   (stress-u      cl-mpm/particle::mp-undamaged-stress)
                   (stress        cl-mpm/particle::mp-stress)
                   (strain        cl-mpm/particle::mp-strain)
                   (p-mod         cl-mpm/particle::mp-p-modulus-0)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (E cl-mpm/particle::mp-e)
                   (nu cl-mpm/particle::mp-nu)
                   (de cl-mpm/particle::mp-elastic-matrix)
                   (enable-damage cl-mpm/particle::mp-enable-damage))
      mp
    (declare (double-float damage damage-t damage-c damage-s))
    (when (and
           enable-damage
           (> damage 0.0d0))
      (let* ((ep (/ (cl-mpm/constitutive::voight-trace strain) 3d0))
             (pind ep)
             (es (cl-mpm/constitutive::deviatoric-voigt strain)))
        (declare (double-float damage-t damage-c damage-s))
        (setf ep
              (if (> pind 0d0)
                  (* (- 1d0 damage-t) ep)
                  (* (- 1d0 damage-c) ep)))
        (setf stress
              (cl-mpm/constitutive::linear-elastic-mat
               (cl-mpm/fastmaths:fast-.+
                (cl-mpm/constitutive::voight-eye ep)
                (cl-mpm/fastmaths:fast-scale! es (- 1d0 damage-s)))
               de
               stress
               ))
        (let* ((K (/ e (* 3 (- 1d0 (* 2 nu)))))
               (G (/ e (* 2 (+ 1d0 nu))))
               (P-0 (+ K (* 4/3 G))))
          (setf K
                (if (> pind 0d0)
                    (* (- 1d0 damage-t) K)
                    (* (- 1d0 damage-c) K)))
          (setf G (* G (- 1d0 damage-s)))
          (setf p-mod (max (* P-0 1d-9) (+ K (* 4/3 G)))))))))

(defun est-shear-from-angle (angle angle-r rc)
  ;; (assert (< angle (/ pi 2)))
  ;; (assert (< angle-r (/ pi 2)))
  (let* ((angle-plastic (cl-mpm/utils:deg-to-rad angle))
         (angle-plastic-residual (cl-mpm/utils:deg-to-rad angle-r)))
    (let ((rs (- 1d0
                 (* (- 1d0 rc)
                    (/ (tan angle-plastic-residual)
                       (tan angle-plastic))))))
      ;;Special case to avoid issues with fp tan
      (when (= angle-r 0d0)
        (setf rs 1d0))
      rs)))


(defun apply-isotropic-porodamage-degredation (mp pressure)
  (with-accessors ((damage        cl-mpm/particle::mp-damage)
                   (undamaged-stress        cl-mpm/particle::mp-undamaged-stress)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (stress        cl-mpm/particle::mp-stress)
                   (enable-damage cl-mpm/particle::mp-enable-damage)
                   (p-mod cl-mpm/particle::mp-p-modulus-0)
                   (e cl-mpm/particle::mp-e)
                   (nu cl-mpm/particle::mp-nu))
      mp
    (declare (double-float damage))
    (when (and
           enable-damage
           (> damage 0.0d0))
      (cl-mpm/utils:voigt-copy-into undamaged-stress stress)
      (cl-mpm/fastmaths:fast-scale! stress (/ (- 1d0 damage) (cl-mpm/fastmaths::det-3x3 def)))
      (cl-mpm/fastmaths::fast-.-
       stress
       (cl-mpm/constitutive::voight-eye pressure)
       stress)
      (setf (cl-mpm/particle::mp-p-modulus mp)
            (*
             (- 1d0 damage)
             (cl-mpm/particle::compute-p-modulus mp)))))
  )
