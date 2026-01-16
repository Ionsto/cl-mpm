(in-package :cl-mpm/damage)

(defun apply-isotropic-degredation (mp)
  (with-accessors ((damage        cl-mpm/particle::mp-damage)
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
      (cl-mpm/fastmaths:fast-scale! stress (- 1d0 damage))
      (setf (cl-mpm/particle::mp-p-modulus mp)
            (*
             (- 1d0 damage)
             (cl-mpm/particle::compute-p-modulus mp))))))

(defun apply-tensile-stress-degredation (mp)
  (with-accessors ((damage        cl-mpm/particle::mp-damage)
                   (stress        cl-mpm/particle::mp-stress)
                   (enable-damage cl-mpm/particle::mp-enable-damage)
                   (p-mod cl-mpm/particle::mp-p-modulus-0)
                   (e cl-mpm/particle::mp-e)
                   (nu cl-mpm/particle::mp-nu)
                   )
      mp
    (let ()
      (when (> damage 0.0d0)
        (multiple-value-bind (l v) (cl-mpm/utils::eig
                                    (voight-to-matrix
                                     stress))

          (let* ((degredation (expt (- 1d0 damage) 1d0)))
            (loop for i from 0 to 2
                  do (let* ((sii (nth i l)))
                       (when (> sii 0d0)
                         ;;Tensile damage -> unbounded
                         (setf (nth i l) (* sii (max 1d-9 degredation))))))
            (setf stress
                  (cl-mpm/utils:matrix-to-voight
                   (magicl:@
                    v
                    (cl-mpm/utils::matrix-diag l)
                    (magicl:transpose v))))
            ))))))

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
                                          (- 1d0 damage)
                                          ;; (max 1d-6 degredation)
                                          )))
                     ;; (when (< esii 0d0)
                     ;;   (setf (nth i l) (* esii (max 1d-1 degredation)))
                     ;;   )
                     ))
          ;; (pprint strain)
          (let ((strain+ (cl-mpm/utils:matrix-to-voigt
                          (magicl:@
                           v
                           (magicl:from-diag l :type 'double-float)
                           (magicl:transpose v)))))
            ;; (pprint strain+)
            (setf stress
                  ;; (magicl:@
                  ;;  de
                  ;;  strain+
                  ;;  )
                  (cl-mpm/constitutive::linear-elastic-mat
                   (cl-mpm/utils:matrix-to-voigt
                    (magicl:@
                     v
                     (cl-mpm/utils::matrix-diag l)
                     (magicl:transpose v)))
                   de
                   stress)
                  )
            (cl-mpm/fastmaths:fast-scale! stress (* 1d0 (cl-mpm/fastmaths:det-3x3 def))))
          ;; (pprint stress)
          )))))


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
