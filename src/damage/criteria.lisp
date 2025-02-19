(in-package :cl-mpm/damage)

(defun modified-vm-stress (stress k init-stress)
  (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d stress)
    (let ((i1 (+ s_1 s_2 s_3)))
      (* 0.5d0
         (+
          (* (/ 1d0 (* k init-stress))
             (+
              (expt (- s_1 s_2) 2)
              (expt (- s_2 s_3) 2)
              (expt (- s_3 s_1) 2))
             )
          (* 2d0 (- 1 (/ 1d0 k)) i1)))
      )))

(defun criterion-modified-vm (strain k nu)
  (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d strain)
    (let ((i1 (+ s_1 s_2 s_3))
          (j2 (cl-mpm/constitutive::voigt-j2 (cl-mpm/constitutive::deviatoric-voigt strain)))
          (k-factor (/ (- k 1d0)
                       (- 1d0 (* 2d0 nu))))
          )
      (* (/ 1d0 (* 2d0 k))
         (+ (* i1 k-factor)
            (sqrt (+ (expt (* k-factor i1) 2)
                     (* (/ (* 12d0 k) (expt (- 1d0 nu) 2)) j2))))))))

(defun criterion-mc (strain angle E nu)
  (multiple-value-bind (e_1 e_2 e_3) (principal-stresses-3d
                                      ;; strain
                                      (magicl:scale!
                                        ;(cl-mpm/utils::voigt-contra->covar strain)
                                        ;(cl-mpm/utils::voigt-contra->covar strain)
                                       strain
                                       -1d0
                                       )
                                      )
    (let ((a (/ (sin angle) (- 1d0 (* 2d0 nu)))))
      (* E
         (/ 1d0 (cos angle))
         (-
          (* (- 1d0 a) e_1)
          (* (+ 1d0 a) e_2)
          )))))
(defun criterion-enhanced-bi-energy-norm (stress strain E nu k)
  (let* ((strain+
           (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix strain))
             (loop for i from 0 to 2
                   do
                      (setf (nth i l) (max (nth i l) 0d0)))
             (cl-mpm/utils:matrix-to-voigt
              (magicl:@ v
                        (cl-mpm/utils::matrix-from-diag l)
                        (magicl:transpose v)))))

         (i1 (cl-mpm/utils:trace-voigt strain))
         (j2
           (cl-mpm/constitutive::voigt-j2
            (cl-mpm/utils::deviatoric-voigt
             (cl-mpm/utils::voigt-contra->covar strain))))
         (j2-factor (sqrt (* 3d0 j2)))

         (et (+ (/ i1 (* 2d0 (- 1d0 (* 2d0 nu))))
                (/ j2-factor (* 2d0 (+ 1d0 nu)))
                ))
         (ec (+ (/ i1 (* 5d0 (- 1d0 (* 2d0 nu))))
                (/ (* 6d0 j2-factor) (* 5d0 (+ 1d0 nu)))
                ))
         ;; (r
         ;;   (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix stress))
         ;;     (loop for i from 0 to 2
         ;;           do
         ;;              (setf (nth i l) (max (nth i l) 0d0)))
         ;;     (/ )(reduce #'+ (mapcar (lambda (x) (max 0d0 x)) l))
         ;;     )
         ;;   )
         )
    ))
(defun criterion-effective-principal-stress (stress pressure)
  (multiple-value-bind (s1 s2 s3) (principal-stresses-3d stress)
    (- (max s1 s2 s3) pressure))
  ;; (let* ((K (/ E (* 3d0 (- 1d0 (* 2d0 nu)))))
  ;;        (ep (* -1/3 (/ pressure K)))
  ;;        (strain+
  ;;          (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix strain))
  ;;            (loop for i from 0 to 2
  ;;                  do
  ;;                     (setf (nth i l) (max (+ (nth i l) ep) 0d0)))
  ;;            (cl-mpm/utils:matrix-to-voigt
  ;;             (magicl:@
  ;;              v
  ;;              (cl-mpm/utils::matrix-from-diag l)
  ;;              (magicl:transpose v))))))
  ;;   (sqrt (max 0d0 (* E (cl-mpm/fastmaths::dot strain+ (magicl:@ de strain+))))))
  )
;;Chalk brittle

(defun drucker-prager-criterion (stress angle)
  (let ((p (cl-mpm/utils::trace-voigt stress))
        (j2 (cl-mpm/constitutive::voigt-j2
             (cl-mpm/utils::deviatoric-voigt stress))))
    (* (/ 3d0 (+ 3 (tan angle)))
       (+ (sqrt (* 3 j2)) (* 1/3 (tan angle) p)))))


(defun criterion-dp (stress angle)
  (criterion-dp-coheasion stress angle)
  ;; (criterion-dp-tensile stress angle)
  )

(defun criterion-dp-coheasion (stress angle)
  "Return some drucker-prager damage criterion from a stress level and and angle (radians)"
  (let ((p (cl-mpm/utils::trace-voigt stress))
        (j2 (cl-mpm/constitutive::voigt-j2
             (cl-mpm/utils::deviatoric-voigt stress)))
        (A (/ (* 6 (cos angle))
              (* (sqrt 3) (- 3d0 (sin angle)))))
        (B (/ (* 2 (sin angle))
              (* (sqrt 3) (- 3d0 (sin angle))))))
    (* (/ 1d0 A)
       (+ (* B p) (sqrt j2)))))

(defun criterion-dp-middle-circumscribe-coheasion (stress angle)
  "Return some drucker-prager damage criterion from a stress level and and angle (radians)"
  (let ((p (cl-mpm/utils::trace-voigt stress))
        (j2 (cl-mpm/constitutive::voigt-j2
             (cl-mpm/utils::deviatoric-voigt stress)))
        (A (/ (* 6 (cos angle))
              (* (sqrt 3) (+ 3d0 (sin angle)))))
        (B (/ (* 2 (sin angle))
              (* (sqrt 3) (+ 3d0 (sin angle))))))
    (* (/ 1d0 A)
       (+ (* B p) (sqrt j2)))))


(defun criterion-dp-inscribe-coheasion (stress angle)
  "Return some drucker-prager damage criterion from a stress level and and angle (radians)"
  (let* ((p (/ (cl-mpm/utils::trace-voigt stress) 3))
         (j2 (cl-mpm/constitutive::voigt-j2
              (cl-mpm/utils::deviatoric-voigt stress)))
         (factor (sqrt (+ 9 (* 3 (expt (sin angle) 2)))))
         (A (/ (* 3 (cos angle))
               factor))
         (B (/ (* (sin angle))
               factor)))
    (* (/ 1d0 A)
       ;; (/ 1d0 (- (/ 1d0 (sqrt 3)) B))
       (+ (* B p) (sqrt j2)))))


(defun degredation-function (stress damage kc kt ks)
  (let ((p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
        (s (cl-mpm/constitutive::deviatoric-voigt stress))
        (damage-t damage)
        (damage-c (* kc damage))
        (damage-s (* ks damage)))
    (declare (double-float damage-t damage-c damage-s))
    (setf p
          (if (> p 0d0)
              (* (- 1d0 damage-t) p)
              (* (- 1d0 damage-c) p)))
    (cl-mpm/fastmaths:fast-.+
     (cl-mpm/constitutive::voight-eye p)
     (cl-mpm/fastmaths:fast-scale! s (- 1d0 damage-s)))))

(defun criterion-y (stress strain damage
                    E
                    kc
                    kt
                    ks)
  (let ((d_0 damage)
        (ddamage (- (min 1d0 (+ damage 1d-3))
                    damage)))
    (if (> ddamage 0d0)
        (let* ((d_1 (+ damage ddamage))
               (se_0 (* 0.5d0
                        (cl-mpm/fastmaths:dot
                         stress
                         strain)))
               (se_1 (* 0.5d0
                        (cl-mpm/fastmaths:dot
                         (degredation-function stress
                                               d_1
                                               kc
                                               kt
                                               ks)
                         strain))))
          (sqrt (* E (max 0d0 (- (/ (- se_1 se_0) ddamage))))))
        0d0)))

;; (let* ((strain (cl-mpm/utils:voigt-from-list (list
;;                                               1d0 1d0 1d0
;;                                               0d0 0d0 1d0
;;                                               )))
;;        (E 1d0)
;;        (nu 0d0)
;;        (de (cl-mpm/constitutive::linear-elastic-matrix E nu))
;;        (stress (magicl:@ de strain))
;;        (kt-r 1d0)
;;        (kc-r 1d0)
;;        (g-r 1d0))
;;   ;; (pprint (degredation-function stress
;;   ;;                               0.1d0
;;   ;;                               1d0
;;   ;;                               1d0
;;   ;;                               1d0
;;   ;;                               ))
;;   (pprint (* 0.5d0
;;              (cl-mpm/fastmaths:dot
;;               stress
;;               strain)))
;;   ;(pprint (criterion-y stress strain damage E kt-r kc-r g-r))
;;   (pprint (volumetric-deviatoric-norm strain E nu kt-r kc-r g-r))
;;   )


(defun criterion-dp-shear-coheasion (stress angle)
  "Return some drucker-prager damage criterion from a stress level and and angle (radians)"
  (let* ((p (cl-mpm/utils::trace-voigt stress))
         (j2 (cl-mpm/constitutive::voigt-j2
              (cl-mpm/utils::deviatoric-voigt stress)))
         (A (cos angle))
         (B (/ (sin angle) 3d0)))
    (* (/ 1d0 A)
       (+ (* B p) (sqrt j2)))))


(defun dp-tension-from-DP-circumscribed (angle coheasion)
  (declare (double-float angle))
  (let* ((factor (* (sqrt 3) (- 3d0 (sin angle))))
         (A (/ (* 6 coheasion (cos angle))
               factor))
         (B (/ (* 2 (sin angle))
               factor)))
    (declare (double-float A B factor))
    (* (/ A (- (/ 1d0 (sqrt 3)) B)))))

(defun dp-tension-from-DP-inscribed (angle coheasion)
  (declare (double-float angle))
  (let* ((factor (sqrt (+ 9 (* 3 (expt (sin angle) 2)))))
         (A (/ (* 3 coheasion (cos angle)) factor))
         (B (/ (sin angle) factor)))
    (declare (double-float A B factor))
    (* (/ A (- (/ 1d0 (sqrt 3)) B)))))


(defun criterion-dp-tensile (stress angle)
  (declare (double-float angle))
  (let* ((p (cl-mpm/utils::trace-voigt stress))
         (j2 (cl-mpm/constitutive::voigt-j2
              (cl-mpm/utils::deviatoric-voigt stress)))
         (factor (* (sqrt 3) (- 3d0 (sin angle))))
         (A (/ (* 6 (cos angle))
               factor))
         (B (/ (* 2 (sin angle))
               factor)))
    (declare (double-float A B factor j2 p))
    (* (/ 1d0 (- (/ 1d0 (sqrt 3)) B))
       (+ (* B p) (the double-float (sqrt j2))))))

(defun criterion-dp-tensile-circumscribed (stress angle)
  (declare (double-float angle))
  (let* ((p (cl-mpm/utils::trace-voigt stress))
         (j2 (cl-mpm/constitutive::voigt-j2
              (cl-mpm/utils::deviatoric-voigt stress)))
         (factor (* (sqrt 3) (- 3d0 (sin angle))))
         (A (/ (* 6 (cos angle))
               factor))
         (B (/ (* 2 (sin angle))
               factor)))
    (declare (double-float A B factor j2 p))
    (* (/ 1d0 (- (/ 1d0 (sqrt 3)) B))
       (+ (* B p) (the double-float (sqrt j2))))))

(defun criterion-dp-tensile-middle-circumscribed (stress angle)
  (declare (double-float angle))
  (let* ((p (cl-mpm/utils::trace-voigt stress))
         (j2 (cl-mpm/constitutive::voigt-j2
              (cl-mpm/utils::deviatoric-voigt stress)))
         (factor (* (sqrt 3) (+ 3d0 (sin angle))))
         (A (/ (* 6 (cos angle))
               factor))
         (B (/ (* 2 (sin angle))
               factor)))
    (declare (double-float A B factor j2 p))
    (* (/ 1d0 (- (/ 1d0 (sqrt 3)) B))
       (+ (* B p) (the double-float (sqrt j2))))))

(defun criterion-dp-tensile-inscribed (stress angle)
  (declare (double-float angle))
  (let* ((p (cl-mpm/utils::trace-voigt stress))
         (j2 (cl-mpm/constitutive::voigt-j2
              (cl-mpm/utils::deviatoric-voigt stress)))
         (factor (* (sqrt (+ 9 (* 3 (expt (sin angle) 2))))))
         (A (/ (* 3 (cos angle))
               factor))
         (B (/ (sin angle)
               factor)))
    (declare (double-float A B factor j2 p))
    (* (/ 1d0 (- (/ 1d0 (sqrt 3)) B))
       (+ (* B p) (the double-float (sqrt j2))))))


(defun criterion-dp-pressure (stress angle pressure)
  (let ((i1 (+ (cl-mpm/utils::trace-voigt stress) (* 3d0 pressure)))
        (j2 (cl-mpm/constitutive::voigt-j2
             (cl-mpm/utils::deviatoric-voigt stress)))
        (B (/ (* 2 (sin angle))
              (* (sqrt 3) (- 3d0 (sin angle))))))
    (* (+ (* B i1) (sqrt j2)))))

(defun criterion-dp-pressure-tensile (stress angle pressure tension)
  (let ((i1 (+ (cl-mpm/utils::trace-voigt stress) (* 3d0 pressure)))
        (j2 (cl-mpm/constitutive::voigt-j2
             (cl-mpm/utils::deviatoric-voigt stress)))
        (B (/ (* 2 (sin angle))
              (* (sqrt 3) (- 3d0 (sin angle)))))
        )
    (if nil;(> i1 (* tension 3))
        (/ i1 3)
        (*  (+ (* B i1) (sqrt j2)))
        )))



;; (defun criterion-dp-strain (stress ft fc)
;;   (let ((p (cl-mpm/utils::trace-voigt stress))
;;         (j2 (cl-mpm/constitutive::voigt-j2
;;              (cl-mpm/utils::deviatoric-voigt
;;               stress
;;               )))
;;         )
;;     (*
;;      (/ 1d0 (* 2 fc))
;;      (- (* (+ fc ft) (sqrt (* 3d0 j2)))
;;         (* (- ft fc) p)))))

(defun criterion-dp-strain (stress k)
  (let ((p (cl-mpm/utils::trace-voigt stress))
        (j2 (cl-mpm/constitutive::voigt-j2
             (cl-mpm/utils::deviatoric-voigt
              stress
              ))))
    (* (/ (sqrt 3) (* 2 (+ k 1)))
       (-
        (sqrt j2)
        (*
         (/ (- k 1)
            (+ k 1))
         p)))
    ))

(defun criterion-j2 (stress)
  (let ((j2 (cl-mpm/constitutive::voigt-j2
             (cl-mpm/utils::deviatoric-voigt stress))))
     (sqrt j2)))

(defun modified-vm-criterion (stress nu k)
  (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d stress)
    (declare (double-float nu k s_1 s_2 s_3))
    (let* ((j2 (cl-mpm/constitutive::voigt-j2
                (cl-mpm/utils::deviatoric-voigt stress)))
           (i1 (+ s_1 s_2 s_3))
           (k-factor (/ (- k 1d0)
                        (- 1d0 (* 2d0 nu))))
           (s_1 (+ (* i1 (/ k-factor (* 2d0 k)))
                   (* (/ 1d0 (* 2d0 k))
                      (sqrt (+ (expt (* k-factor i1) 2)
                               (* (/ (* 12 k) (expt (- 1d0 nu) 2))j2)
                               ))))))
      s_1
      )))
(defun smooth-rankine-criterion (stress)
  (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d stress)
     (sqrt
      (+ (expt (max 0d0 s_1) 2)
         (expt (max 0d0 s_2) 2)
         (expt (max 0d0 s_3) 2)))
     ))
(defun criterion-max-principal-stress (stress)
  (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d stress)
    (declare (double-float s_1 s_2 s_3))
    (max 0d0 s_1 s_2 s_3)))

(defun criterion-max-principal-strain (strain E)
  (multiple-value-bind (s_1 s_2 s_3) (principal-strains strain)
    (declare (double-float s_1 s_2 s_3))
    (* E (max 0d0 s_1 s_2 s_3))))
;; (defun modified-vm-strain (strain nu k E)
;;   (multiple-value-bind (e1 e2 e3)
;;       (principal-stresses-3d
;;        (magicl:.*
;;         strain
;;         (cl-mpm/utils:voigt-from-list
;;          (list 1d0 1d0 1d0
;;                0.5d0
;;                0.5d0
;;                0.5d0))))
;;     (let* ((i1 (+ e1 e2 e3))
;;            (j2
;;              ;; (cl-mpm/constitutive::voigt-j2 (cl-mpm/utils::deviatoric-voigt strain))
;;              (* 1/6 (+
;;                      (expt (- e1 e2) 2)
;;                      (expt (- e2 e3) 2)
;;                      (expt (- e3 e1) 2)))
;;              )
;;            (k-factor (/ (- k 1d0)
;;                         (- 1d0 (* 2d0 nu))))
;;            (s_1 (* E
;;                    (+ (* i1 (/ k-factor (* 2d0 k)))
;;                       (* (/ 1d0 (* 2d0 k))
;;                          (sqrt (+ (expt (* k-factor i1) 2)
;;                                   (* (/ (* 12 k) (expt (- 1d0 nu) 2)) j2)
;;                                   )))))))
;;       s_)))

(defun principal-strains (strains)
  (multiple-value-bind (l v) (cl-mpm/utils::eig (voigt-to-matrix strains))
    (declare (ignore v))
    (setf l (sort l #'>))
    (values (nth 0 l) (nth 1 l) (nth 2 l))))

(defun criterion-mohr-coloumb (strain angle E nu)
  (multiple-value-bind (e1 e2 e3) (principal-strains strain)
    (let ((a (/ (sin angle)
                (- 1d0 (* 2 nu))))
          (G (/ E (* 2 (- 1d0 nu))))
          )
      (* G
         (- (* (- 1d0 a) e3)
            (* (+ 1d0 a) e1))
         (/ 1d0 (cos angle)))
      ;; (/
      ;;  (* E (- (* (- 1d0 a) e1)
      ;;          (* (+ 1d0 a) e3)
      ;;          )
      ;;     (/ 1d0 (cos angle)))
      ;;  (* (- 1d0 a) (/ 1d0 (cos angle))))
      )))

(defun criterion-mohr-coloumb-3d (stress angle)
  (multiple-value-bind (s1 s2 s3) (principal-stresses-3d stress)
    (declare (double-float angle s1 s2 s3))
    (let ((tang (tan angle))
          (cang (cos angle)))
      (* 0.5d0
         (max
          (- (* (+ s3 s1) tang) (/ (- s3 s1) cang))
          ;; (- (- (* (+ s2 s3) tang) (/ (abs (- s2 s3)) cang)))
          ;; (- (- (* (+ s1 s2) tang) (/ (abs (- s1 s2)) cang)))
          )))))

(defun criterion-mohr-coloumb-stress (stress angle)
  (multiple-value-bind (s1 s2 s3) (principal-stresses-3d stress)
    (declare (double-float angle s1 s2 s3))
    (let ()
      (* 0.5d0
         (-
          (* (+ s3 s1) (tan angle))
          (/ (- s3 s1) (cos angle)))))))

(defun criterion-mohr-coloumb-stress-tensile (stress angle)
  (declare (double-float angle))
  (multiple-value-bind (s1 s2 s3) (principal-stresses-3d stress)
    (declare (double-float angle s1 s2 s3))
    (let ((k (/ (+ 1d0 (sin angle))
                (- 1d0 (sin angle)))))
      (max 0d0
           (/ (- (* k s1) s3)
              k)))))

(defun criterion-mohr-coloumb-stress-will (stress angle)
  (multiple-value-bind (s1 s2 s3) (principal-stresses-3d stress)
    (declare (double-float angle s1 s2 s3))
    (let ((k (/ (+ 1d0 (sin angle))
                (- 1d0 (sin angle))
                )))
      (* 0.5d0
         (/ (- (* k s1) s3)
            (sqrt k))))))

(defun mohr-coloumb-coheasion-to-tensile (coheasion angle)
  (let (
        ;; (a (/ (sin angle)
        ;;       (- 1d0 (* 2 nu))))
        ;; (G (/ E (* 2 (- 1d0 nu))))
        (k (/ (+ 1d0 (sin angle))
              (- 1d0 (sin angle))
              )))
    ;; (/ (* 2 coheasion) (sqrt k))
    (the double-float
         (/
          (* 2 coheasion (cos angle))
          (+ 1d0 (sin angle))))
    ;; (* E
    ;;    0.5d0
    ;;    (/
    ;;     (/ coheasion G)
    ;;     (* (- 1d0 a) (/ 1d0 (cos angle)))))
    ))


;; (let* ((stress (cl-mpm/utils:voigt-from-list (list 1d0 0d0 0d0
;;                                                    0d0 0d0 1d0)))
;;        (angle 42d0))
;;   (print (cl-mpm/damage::criterion-mohr-coloumb-stress stress (* angle (/ pi 180))))
;;   (print (cl-mpm/damage::criterion-mohr-coloumb-3d stress (* angle (/ pi 180))))
;;   ;; (print (cl-mpm/damage::criterion-mohr-coloumb-stress-will stress (* angle (/ pi 180))))
;;   ;; (print (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile stress (* angle (/ pi 180))))
;;   )

;;   (let ((mc (cl-mpm/damage::criterion-mohr-coloumb-stress stress (* angle (/ pi 180)))))
;;     (format t "MC: ~%~E~%" mc)
;;     (format t "Insc: ~E~%"  (/ (cl-mpm/damage::criterion-dp-inscribe-coheasion stress (* angle (/ pi 180))) mc))
;;     (format t "Shear: ~E~%"  (/ (cl-mpm/damage::criterion-dp-shear-coheasion stress (* angle (/ pi 180))) mc))
;;     (format t "Mid: ~E~%"  (/ (cl-mpm/damage::criterion-dp-middle-circumscribe-coheasion stress (* angle (/ pi 180))) mc))
;;     (format t "Circ: ~E~%"  (/ (cl-mpm/damage::criterion-dp-coheasion stress (* angle (/ pi 180))) mc)))
;;   )

(defun volumetric-deviatoric-norm (strain E nu
                                   res-t res-c res-s)
  (let* ((etr (cl-mpm/utils:trace-voigt strain))
         (e-tension (max 0d0 etr))
         (e-compression (- (max -0d0 (- etr))))
         (ed (cl-mpm/constitutive::deviatoric-voigt strain))
         (K (/ E (* 3d0 (- 1d0 (* 2d0 nu)))))
         (G (/ E (* 2 (+ 1d0 nu))))
         (energy
           (* 0.5d0
              (+
               (* res-t K (expt e-tension 2))
               (* res-c K (expt e-compression 2))
               (*
                G
                (cl-mpm/fastmaths::dot ed ed))))))
    (sqrt (max 0d0 (* E energy)))
    ))

(defun tensile-energy-norm (strain E de)
  (let* ((strain+
           (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix strain))
             (loop for i from 0 to 2
                   do
                      (setf (nth i l) (max (nth i l) 0d0)))
             (cl-mpm/utils:matrix-to-voigt
              (magicl:@ v
                        (cl-mpm/utils::matrix-from-diag l)
                        (magicl:transpose v))))))
    (sqrt (max 0d0 (* E (cl-mpm/fastmaths::dot strain+ (magicl:@ de strain+)))))))

(defun tensile-energy-norm-pressure (strain E nu de pressure)
  (let* ((K (/ E (* 3d0 (- 1d0 (* 2d0 nu)))))
         (ep (* -1/3 (/ pressure K)))
         (strain+
           (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix strain))
             (loop for i from 0 to 2
                   do
                      (setf (nth i l) (max (+ (nth i l) ep) 0d0)))
             (cl-mpm/utils:matrix-to-voigt
              (magicl:@
               v
               (cl-mpm/utils::matrix-from-diag l)
               (magicl:transpose v))))))
    (sqrt (max 0d0 (* E (cl-mpm/fastmaths::dot strain+ (magicl:@ de strain+)))))))

(defun criterion-mohr-coloumb-rankine-stress-tensile (stress angle)
  (declare (double-float angle))
  (multiple-value-bind (s1 s2 s3) (principal-stresses-3d stress)
    (declare (double-float angle s1 s2 s3))
    (max
     0d0
     s1
     s2
     s3
     (let ((k (/ (+ 1d0 (sin angle))
                 (- 1d0 (sin angle)))))
       (/ (- (* k s1) s3)
          k))
     )))


