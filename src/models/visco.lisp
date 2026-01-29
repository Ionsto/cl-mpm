(defpackage :cl-mpm/models/visco
  (:use :cl :cl-mpm/utils :cl-mpm/fastmaths)
  (:export))
(declaim (optimize (debug 3) (safety 3) (speed 0)))
(in-package :cl-mpm/models/visco)
(defun finite-strain-linear-viscous (stress strain de e nu dt viscosity)
  (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix strain))
      (let* ((dev (cl-mpm/fastmaths:fast-.-
                   (cl-mpm/utils:matrix-from-list (list
                                                   1d0 0d0 0d0
                                                   0d0 1d0 0d0
                                                   0d0 0d0 1d0))
                   (cl-mpm/fastmaths:fast-scale!
                    (cl-mpm/utils:matrix-from-list
                     (list 1d0 1d0 1d0
                           1d0 1d0 1d0
                           1d0 1d0 1d0))
                    (/ 1d0 3d0))))
             (K (/ e (* 3 (- 1d0 (* 2 nu)))))
             (G (/ E (* 2 (+ 1d0 nu))))
             (identity (cl-mpm/utils:matrix-eye 1d0))
             (d-neq
               (cl-mpm/fastmaths::fast-scale!
                (cl-mpm/utils:matrix-from-list (list
                                                (- 1d0 nu) nu nu
                                                nu (- 1d0 nu) nu
                                                nu nu (- 1d0 nu)))
                (/ E (* (+ 1d0 nu) (- 1d0 (* 2d0 nu)))))
               ;; (cl-mpm/fastmaths:fast-.+
               ;;  (fast-scale identity (/ K 3))
               ;;  (fast-scale
               ;;   dev
               ;;   (* 2 G)))
               ;; (fast-scale
               ;;  dev
               ;;  (* 2 G))
               )
             (epsTr (cl-mpm/utils:vector-from-list l))
             (en (cl-mpm/utils::deep-copy epsTr))
             (f-tol 1d-5)
             (beta (magicl:@ d-neq en))
             (C (magicl:inv d-neq))
             (a
               (magicl:@
                C
                (magicl:inv
                 (cl-mpm/fastmaths:fast-.+
                  C
                  (magicl:scale
                   dev
                   (/ dt viscosity))))))
             (f f-tol))
        (pprint d-neq)
        (pprint de)
        (loop for i from 0 to 100
              while (>= f f-tol)
              do
                 (progn
                   (setf beta (magicl:@ d-neq en))
                   (let* ((r (magicl:.-
                              epsTr
                              (magicl:.+ en
                                         (magicl:scale
                                          (magicl:@ dev beta)
                                          (/ dt (* 1d0 viscosity)))))))
                     (setf f (magicl:norm r))
                     (when (>= f f-tol)
                       (setf
                        en
                        (magicl:.+
                         en
                         (magicl:@ a r)
                         ;; (magicl:linear-solve d-neq
                         ;;                      (magicl:@ a r))
                         ))))))
        (when (> f f-tol)
          (pprint beta)
          (pprint en)
          (pprint epsTr)
          (pprint viscosity)
          (pprint (/ dt viscosity))
          (error "Didn't Converged!~%"))
        (let ((out-strain (cl-mpm/utils:matrix-to-voigt
                           (magicl:@ v
                                     (cl-mpm/utils::matrix-diag (list (varef en 0) (varef en 1) (varef en 2)))
                                     (magicl:transpose v)))))
          (cl-mpm/utils:voigt-copy-into out-strain strain)
          (cl-mpm/constitutive::linear-elastic-mat strain de stress)
          out-strain
          ;; (pprint de)
          ;; (magicl:@ de out-strain)
          ))))


(defun voight-trace-2d (s)
  (+ (varef s 0)
     (varef s 1)))
(defun voigt-eye-2d (v)
  (cl-mpm/utils:voigt-from-list
   (list
    v
    v
    0d0
    0d0
    0d0
    0d0)))

(defun dev-exp-v (stress strain-n strain e nu de viscosity dt)
  "A stress increment form of a viscoelastic maxwell material"
  (cl-mpm::voigt-copy-into
   (dev-exp-v-nostrain stress strain-n strain e nu de viscosity dt)
   stress)
  (cl-mpm::voigt-copy-into
   (magicl:linear-solve de stress)
   strain)
  stress)
(defun dev-exp-v-nostrain (stress strain-n strain e nu de viscosity dt)
  "A stress increment form of a viscoelastic maxwell material"
  (let* ((stress-n (cl-mpm/constitutive:linear-elastic-mat strain-n de))
         (stress-inc (cl-mpm/constitutive:linear-elastic-mat (cl-mpm/fastmaths:fast-.- strain strain-n) de))
         (pressure (/ (trace-voigt stress-n) 3d0))
         (pressure-inc (/ (trace-voigt stress-inc) 3d0))
         (pmat (voigt-eye pressure))
         (dev (cl-mpm/fastmaths:fast-.- stress-n pmat))
         (dev-inc (cl-mpm/fastmaths:fast-.- stress-inc (voigt-eye pressure-inc)))
         (rho (/ viscosity e))
         (dy (/ dt rho))
         (exp-rho (exp (- dy)))
         (lam (/ (- 1 exp-rho) dy)))
    (cl-mpm/fastmaths:fast-.+
     (voigt-eye (+ pressure pressure-inc))
     (cl-mpm/fastmaths:fast-.+
      (cl-mpm/fastmaths:fast-scale dev exp-rho)
      (cl-mpm/fastmaths:fast-scale dev-inc lam))
     ;; (cl-mpm/fastmaths:fast-scale
     ;;  (cl-mpm/fastmaths:fast-.- stress pmat)
     ;;  exp-rho)
     )))

;; (defun dev-exp-v (stress strain-n strain e nu de viscosity dt)
;;   "A stress increment form of a viscoelastic maxwell material"
;;   (let* (

;;          (pressure (/ (voight-trace-2d strain) 2d0))
;;          (pmat (voigt-eye-2d pressure))
;;          ;; (rho (/ (* 2d0 (- 1d0 nu) viscosity) e))
;;          ;; (exp-rho (exp (- (/ dt rho))))
;;          (exp-rho (exp (- (/ dt viscosity)))))
;;     (cl-mpm::voigt-copy-into
;;           (cl-mpm/fastmaths:fast-.+
;;            pmat
;;            (cl-mpm/fastmaths:fast-scale!
;;             (cl-mpm/fastmaths:fast-.- strain pmat)
;;             exp-rho))
;;           strain)
;;     (cl-mpm/constitutive::linear-elastic-mat strain de stress)))

(defun test (dt)
  (let* ((E 1d0)
         (nu 0d0)
         ;; (dt 4d0)
         (viscosity 1d0)
         (strain-n (cl-mpm/utils:voigt-from-list (list 0d0 0d0 0d0 0d0 0d0 0d0)))
         (strain-0 (cl-mpm/utils:voigt-from-list (list 0d0 0d0 0d0 1d0 0d0 0d0)))
         (de (cl-mpm/constitutive:linear-elastic-matrix E nu)))
    (let* ((strain (cl-mpm/utils::deep-copy strain-0))
           (stress (magicl:@ de strain)))
                                        ;(pprint (finite-strain-linear-viscous stress strain de e nu dt viscosity))
      ;; (pprint (dev-exp-v stress strain-n strain e nu de viscosity dt))
      ;; (pprint stress)
      (pprint (cl-mpm/ext::constitutive-viscoelastic stress strain de e nu dt viscosity))
      )

    (let* ((strain (cl-mpm/utils::deep-copy strain-0))
           (stress (magicl:@ de strain)))
      ;; (pprint stress)
      ;; (pprint (cl-mpm/ext::constitutive-viscoelastic stress de strain e nu dt viscosity))
      (pprint (finite-strain-linear-viscous stress strain de e nu dt viscosity))
      )
    )

  )

;; (defun matrix-to-column-major (matrix)
;;   (if nil;(= (magicl:layout matrix) :column-major)
;;       matrix
;;       (let* ((s (cl-mpm/utils:fast-storage matrix))
;;              (mt (cl-mpm/utils::empty-copy matrix))
;;              (rows (magicl::matrix/double-float-nrows matrix))
;;              (cols (magicl::matrix/double-float-ncols matrix))
;;              (st (cl-mpm/utils:fast-storage mt)))
;;         (declare (fixnum rows cols))
;;         (loop for i fixnum from 0 below rows
;;               do (loop for j fixnum from 0 below cols
;;                        do
;;                           (progn
;;                             (setf (aref st (+ i (* j rows)))
;;                                   (aref s (+ j (* i rows)))))))
;;         mt)))

;; (declaim (optimize (debug 3) (safety 3) (speed 0)))
;; (defun testt ()
;;   (let* ((de (cl-mpm/constitutive:linear-elastic-matrix 1d0 0d0))
;;          (C (magicl:inv de)))
;;     (magicl:.+
;;      C
;;      (magicl:scale
;;       (magicl:.-
;;        (magicl:eye 6)
;;        (magicl:scale
;;         (magicl:from-list
;;          (list 1d0 1d0 1d0 0d0 0d0 0d0
;;                1d0 1d0 1d0 0d0 0d0 0d0
;;                1d0 1d0 1d0 0d0 0d0 0d0
;;                0d0 0d0 0d0 0d0 0d0 0d0
;;                0d0 0d0 0d0 0d0 0d0 0d0
;;                0d0 0d0 0d0 0d0 0d0 0d0
;;                )
;;          '(6 6)) (/ 1d0 3d0)))
;;       1d0))))

(defun finite-glen-flow (stress strain de e nu dt viscosity)
  (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix strain))
      (let* ((dev (cl-mpm/fastmaths:fast-.-
                   (cl-mpm/utils:matrix-from-list (list
                                                   1d0 0d0 0d0
                                                   0d0 1d0 0d0
                                                   0d0 0d0 1d0))
                   (cl-mpm/fastmaths:fast-scale!
                    (cl-mpm/utils:matrix-from-list
                     (list 1d0 1d0 1d0
                           1d0 1d0 1d0
                           1d0 1d0 1d0))
                    (/ 1d0 3d0))))
             (De3
               (cl-mpm/fastmaths::fast-scale!
                (cl-mpm/utils:matrix-from-list (list
                                                (- 1d0 nu) nu nu
                                                nu (- 1d0 nu) nu
                                                nu nu (- 1d0 nu)))
                (/ E (* (+ 1d0 nu) (- 1d0 (* 2d0 nu))))))
             (power 3d0)
             (epsTr (cl-mpm/utils:vector-from-list l))
             (en (cl-mpm/utils::deep-copy epsTr))
             (f-tol 1d-5)
             (beta (magicl:@ De3 en))
             (C (magicl:inv De3))
             (a
               (magicl:inv
                (cl-mpm/fastmaths:fast-.+
                 C
                 (magicl:scale
                  dev
                  (*
                   (expt (/ 1d0 (sqrt (* 0.5d0 (max (cl-mpm/fastmaths::mag-squared (magicl:@ dev beta))
                                                    1d-20
                                                    )))) power)
                   (/ dt viscosity))))))
             (f f-tol))
        (loop for i from 0 to 1000
              while (>= f f-tol)
              do
                 (progn
                   (setf beta (magicl:@ De3 en))
                   (let* ((r (magicl:.-
                              (magicl:.+ en
                                         (magicl:scale
                                          (magicl:@ dev beta)
                                          (*
                                           (expt (/ 1d0 (sqrt (* 0.5d0 (max (cl-mpm/fastmaths::mag-squared (magicl:@ dev beta))
                                                                            1d-20
                                                                            )))) power)
                                           (/ dt (* 2d0 viscosity)))))
                              epsTr)))
                     ;; (pprint (sqrt (* 0.5d0 (expt (+ (varef beta 0) (varef beta 1) (varef beta 2)) 2))))
                     (setf f (magicl:norm r))
                     (when (>= f f-tol)
                       (setf
                        en
                        (magicl:.-
                         en
                         (magicl:@
                          C
                          (magicl:@ a r)
                          ;; (magicl:linear-solve a r)
                          )))))))
        (when (> f f-tol)
          (error "Didn't Converged!~%"))
        (let ((out-strain (cl-mpm/utils:matrix-to-voigt
                           (magicl:@ v
                                     (cl-mpm/utils::matrix-diag (list (varef en 0) (varef en 1) (varef en 2)))
                                     (magicl:transpose v)))))
          (cl-mpm/utils:voigt-copy-into out-strain strain)
          (cl-mpm/constitutive::linear-elastic-mat strain de stress)
          out-strain))))

(defun finite-strain-full-viscous (stress strain k u k-neq u-neq dt viscosity)
  (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix strain))
      (let* ((dev (cl-mpm/fastmaths:fast-.-
                   (cl-mpm/utils:matrix-from-list (list
                                                   1d0 0d0 0d0
                                                   0d0 1d0 0d0
                                                   0d0 0d0 1d0))
                   (cl-mpm/fastmaths:fast-scale!
                    (cl-mpm/utils:matrix-from-list
                     (list 1d0 1d0 1d0
                           1d0 1d0 1d0
                           1d0 1d0 1d0))
                    (/ 1d0 3d0))))
             ;; (G (/ E (* 2 (+ 1d0 nu))))
             (identity (cl-mpm/utils:matrix-eye 1d0))
             (d-neq
               (cl-mpm/fastmaths:fast-.+
                (cl-mpm/fastmaths:fast-scale identity k-neq)
                (cl-mpm/fastmaths:fast-scale dev u-neq)))
             (d
               (cl-mpm/fastmaths:fast-.+
                (cl-mpm/fastmaths:fast-scale identity k)
                (cl-mpm/fastmaths:fast-scale dev u)))

             (epsTr (cl-mpm/utils:vector-from-list l))
             (en (cl-mpm/utils::deep-copy epsTr))
             (f-tol 1d-5)
             (beta (magicl:@ d-neq en))
             (C (magicl:inv d-neq))
             (a
               (magicl:@
                C
                (magicl:inv
                 (cl-mpm/fastmaths:fast-.+
                  C
                  (magicl:scale
                   dev
                   (/ dt viscosity))))))
             (f f-tol))
        (loop for i from 0 to 100
              while (>= f f-tol)
              do
                 (progn
                   (setf beta (magicl:@ d-neq en))
                   (let* ((r (magicl:.-
                              (magicl:.+ en
                                         (magicl:scale
                                          (magicl:@ dev beta)
                                          (/ dt (* 2d0 viscosity))))
                              epsTr)))
                     (setf f (magicl:norm r))
                     (when (>= f f-tol)
                       (setf
                        en
                        (magicl:.-
                         en
                         (magicl:@ a r)))))))
        (when (> f f-tol)
          (error "Didn't Converged!~%"))
        (let ((out-strain (cl-mpm/utils:matrix-to-voigt
                           (magicl:@ v
                                     (cl-mpm/utils::matrix-diag (list (varef en 0) (varef en 1) (varef en 2)))
                                     (magicl:transpose v)))))
          (cl-mpm/utils:voigt-copy-into out-strain strain)
          (cl-mpm/constitutive::linear-elastic-mat strain d stress)
          out-strain))))


(in-package :cl-mpm/particle)
(defclass particle-finite-viscoelastic (particle-elastic)
  ((viscosity
    :accessor mp-viscosity
    :initarg :viscosity
    :initform 1d0
    )
   (enable-viscosity
    :accessor mp-enable-viscosity
    :initarg :enable-viscosity
    :initform t)))

(defclass particle-finite-viscoelastic-ice (particle-finite-viscoelastic)
  ((visc-factor
    :accessor mp-visc-factor
    :initarg :visc-factor)
   (visc-power
    :accessor mp-visc-power
    :initarg :visc-power)))
(defmethod constitutive-model ((mp particle-finite-viscoelastic) strain dt)
  (with-accessors ((de mp-elastic-matrix)
                   (e mp-e)
                   (nu mp-nu)
                   (p-wave-0 mp-p-modulus-0)
                   (stress mp-stress)
                   (strain mp-strain)
                   (strain-n mp-strain-n)
                   (viscosity mp-viscosity)
                   (enable-viscosity mp-enable-viscosity))
      mp
    ;;Train elastic strain - plus trail kirchoff stress
    (if enable-viscosity
        ;; (cl-mpm/models/visco::finite-strain-linear-viscous stress strain de e nu dt viscosity)
        (cl-mpm/ext::constitutive-viscoelastic stress strain de e nu dt viscosity)
        (cl-mpm/constitutive:linear-elastic-mat strain de stress))
    (when (and enable-viscosity
               (not (= dt 0d0)))
      (let* ((K (/ e (* 3 (- 1d0 (* 2 nu)))))
             (G (/ e (* 2 (+ 1d0 nu))))
             (rho (/ viscosity G))
             (dy (/ dt rho))
             (exp-rho (exp (- dy)))
             (lam (/ (- 1 exp-rho) dy))
             ;; (lam 1d0)
             )
        (declare (double-float K G p-wave-0))
        (setf p-wave-0 (* (+ K (* 4/3 G lam)))))
      ;; (cl-mpm/particle::update-p-modulus mp)
      )
    stress))

;; (let* ((effective-stress 1d4)
;;        (visc-factor 111d6)
;;        (visc-power 3d0)
;;        (visc-factor (expt visc-factor (- visc-power))))
;;   (declare (double-float effective-stress visc-factor))
;;   (pprint
;;    (/ 1d0
;;       (+ 1d-38
;;          (* 2d0 visc-factor (expt effective-stress (/ (- visc-power 1) 2))))))
;;   )

(defun d-glen-visco (stress visc-factor visc-power)
  (let* ((effective-stress (* 1/2 (cl-mpm/fastmaths::voigt-j2 (cl-mpm/utils:deviatoric-voigt stress))))
         (visc-factor (expt visc-factor (- visc-power))))
    (declare (double-float effective-stress visc-factor))
    (/ 1d0
       (+ 1d-38
          (* 2d0 visc-factor (expt effective-stress (/ (- visc-power 1) 2)))))
    ;; (/
    ;;    (expt effective-stress (/ (- 1 visc-power) 2))
    ;;    (* 2d0 visc-factor ))
    ))

(defun glen-visco (stress visc-factor visc-power)
  (let* ((effective-stress (* ;1/2
                              (cl-mpm/fastmaths::voigt-j2 (cl-mpm/utils:deviatoric-voigt stress))))
         (visc-factor (expt visc-factor (- visc-power))))
    (declare (double-float effective-stress visc-factor))
    (/ 1d0
       (+ 1d-38
          (* 2d0 visc-factor (expt effective-stress (/ (- visc-power 1) 2)))))
    ;; (/
    ;;    (expt effective-stress (/ (- 1 visc-power) 2))
    ;;    (* 2d0 visc-factor ))
    ))

;; (defun visco-predictor (mp dt)
;;   (with-accessors ((de mp-elastic-matrix)
;;                    (e mp-e)
;;                    (nu mp-nu)
;;                    (stress mp-stress)
;;                    (def mp-deformation-gradient)
;;                    (def-n mp-deformation-gradient-0)
;;                    (strain mp-strain)
;;                    (strain-n mp-strain-n)
;;                    (true-visc mp-viscosity)
;;                    (visc-factor mp-visc-factor)
;;                    (visc-power mp-visc-power)
;;                    (enable-viscosity mp-enable-viscosity)
;;                    )
;;       mp
;;     (cl-mpm/constitutive:linear-elastic-mat strain de stress)
;;     (let* ((visc-n (glen-visco
;;                     (cl-mpm/fastmaths:fast-scale!
;;                      (cl-mpm/constitutive:linear-elastic-mat strain-n de)
;;                      (/ 1d0 (cl-mpm/fastmaths:det def-n)))
;;                     visc-factor
;;                     visc-power))
;;            (visc-n1 (glen-visco
;;                      (cl-mpm/fastmaths:fast-scale!
;;                       stress
;;                       (/ 1d0 (cl-mpm/fastmaths:det def)))
;;                      visc-factor
;;                      visc-power))
;;            (visc (expt (+ (expt visc-n -1) (expt visc-n1 -1)) -1)))
;;       (if enable-viscosity
;;           (cl-mpm/models/visco::dev-exp-v stress strain e nu de visc dt)
;;           (cl-mpm/constitutive:linear-elastic-mat strain de stress)))))
(defun visco-iterator (mp dt)
  (with-accessors ((de mp-elastic-matrix)
                   (e mp-e)
                   (nu mp-nu)
                   (stress mp-stress)
                   (def mp-deformation-gradient)
                   (def-n mp-deformation-gradient-0)
                   (strain mp-strain)
                   (strain-n mp-strain-n)
                   (true-visc mp-viscosity)
                   (visc-factor mp-visc-factor)
                   (visc-power mp-visc-power)
                   (enable-viscosity mp-enable-viscosity)
                   (p-mod mp-p-modulus)
                   (p-wave-0 mp-p-modulus-0)
                   )
      mp
    (let* ((iJ-n1 (/ 1d0 (cl-mpm/fastmaths:det def)))
           (stress-tr (cl-mpm/constitutive:linear-elastic-mat strain de))
           (stress-n1 (cl-mpm/utils:voigt-copy stress-tr))
           (visc-prev
             (glen-visco
              (cl-mpm/fastmaths:fast-scale
               stress-tr
               iJ-n1)
              visc-factor
              visc-power))
           (visc visc-prev)
           (visc-0 visc)
           (tol 1d-6)
           (err tol)
           )
      (loop for i from 0 to 1000
            while (>= err tol)
            do
               (progn
                 (setf stress-n1 (cl-mpm/models/visco::dev-exp-v-nostrain stress-tr strain-n strain e nu de visc dt))
                 (let ((new-visc
                         (glen-visco
                          (cl-mpm/fastmaths:fast-scale
                           stress-n1
                           iJ-n1)
                          visc-factor
                          visc-power)))
                   (setf err (/ (abs (- new-visc visc)) visc))
                   ;; (format t "Visc iter ~D ~E ~E - err ~E~%" i new-visc visc err)
                   (setf
                    visc-prev visc)
                   (setf
                    visc new-visc))))
      ;; (format t "~E ~E ~%" visc-0 visc)
      visc)))

(defun visco-iterator-midpoint (mp dt)
  (with-accessors ((de mp-elastic-matrix)
                   (e mp-e)
                   (nu mp-nu)
                   (stress mp-stress)
                   (def mp-deformation-gradient)
                   (def-n mp-deformation-gradient-0)
                   (strain mp-strain)
                   (strain-n mp-strain-n)
                   (true-visc mp-viscosity)
                   (visc-factor mp-visc-factor)
                   (visc-power mp-visc-power)
                   (enable-viscosity mp-enable-viscosity)
                   (p-mod mp-p-modulus)
                   (p-wave-0 mp-p-modulus-0)
                   )
      mp
    (let* ((iJ-n1 (/ 1d0 (cl-mpm/fastmaths:det def)))
           (iJ-n (/ 1d0 (cl-mpm/fastmaths:det def-n)))
           (stress-tr (cl-mpm/constitutive:linear-elastic-mat strain de))
           (stress-n (cl-mpm/constitutive:linear-elastic-mat strain-n de))
           ;; (stress-n1 (cl-mpm/utils:voigt-copy stress-tr))
           (visc-n (glen-visco
                    (cl-mpm/fastmaths:fast-scale
                     stress-n
                     iJ-n)
                    visc-factor
                    visc-power))
           (visc-prev
             visc-n)
           (visc visc-prev)
           (tol 1d-9)
           (err tol)
           )
      (flet ((visc-compute (visc-est)
               (* 0.5d0
                  (+ visc-n
                     (glen-visco
                      (cl-mpm/fastmaths:fast-scale
                       (cl-mpm/models/visco::dev-exp-v-nostrain stress-tr strain-n strain e nu de visc-est dt)
                       iJ-n1)
                      visc-factor
                      visc-power)))))
        (loop for i from 0 to 10
              while (>= err tol)
              do
                 (progn
                   ;; (setf stress-n1 (cl-mpm/models/visco::dev-exp-v-nostrain stress-tr strain-n strain e nu de visc dt))
                   ;; (pprint stress-n1)
                   (let ((new-visc (visc-compute visc)))
                     (setf err (/ (abs (- new-visc visc)) visc))
                     ;; (format t "Visc iter ~D ~E ~E - err ~E~%" i new-visc visc err)
                     (setf
                      visc-prev visc)
                     (setf
                      visc new-visc)))))
      visc)))

;; (defun cl-mpm/particle::visco-nr-midpoint (mp dt)
;;   (with-accessors ((de mp-elastic-matrix)
;;                    (e mp-e)
;;                    (nu mp-nu)
;;                    (stress mp-stress)
;;                    (def mp-deformation-gradient)
;;                    (def-n mp-deformation-gradient-0)
;;                    (strain mp-strain)
;;                    (strain-n mp-strain-n)
;;                    (true-visc mp-viscosity)
;;                    (visc-factor mp-visc-factor)
;;                    (visc-power mp-visc-power)
;;                    (enable-viscosity mp-enable-viscosity)
;;                    (p-mod mp-p-modulus)
;;                    (p-wave-0 mp-p-modulus-0)
;;                    )
;;       mp
;;     (let* ((iJ-n1 (/ 1d0 (cl-mpm/fastmaths:det def)))
;;            (iJ-n (/ 1d0 (cl-mpm/fastmaths:det def-n)))
;;            (stress-tr (cl-mpm/constitutive:linear-elastic-mat strain de))
;;            (stress-n (cl-mpm/constitutive:linear-elastic-mat strain-n de))
;;            (stress-n1 (cl-mpm/utils:voigt-copy stress-tr))
;;            (visc-n (glen-visco
;;                     (cl-mpm/fastmaths:fast-scale
;;                      stress-n
;;                      iJ-n)
;;                     visc-factor
;;                     visc-power))
;;            (visc-prev
;;              visc-n)
;;            (visc visc-prev)
;;            (tol 1d-9)
;;            (err tol)
;;            )
;;       (
;;flet ((visc-compute (visc-est)
;;                (* 0.5d0
;;                   (+ visc-n
;;                      (glen-visco
;;                       (cl-mpm/fastmaths:fast-scale
;;                        (cl-mpm/models/visco::dev-exp-v-nostrain stress-tr strain-n strain e nu de visc-est dt)
;;                        iJ-n1)
;;                       visc-factor
;;                       visc-power)))))
;;         ;; (setf visc (visc-compute visc-n))
;;         (pprint visc)
;;         (loop for i from 0 to 10000
;;               ;; while (>= err tol)
;;               do
;;                  (progn
;;                    ;; (setf stress-n1 (cl-mpm/models/visco::dev-exp-v-nostrain stress-tr strain-n strain e nu de visc dt))
;;                    (let* ((dx (* -0.01d0 visc))
;;                           (tangent 100000)
;;                           ;; (tangent (/ (-
;;                           ;;              (visc-compute (+ visc dx))
;;                           ;;              (visc-compute visc))
;;                           ;;             dx))
;;                           )
;;                      (format t "Residual ~E - tangent ~E~%" (- (visc-compute visc) visc) tangent)
;;                      ;; (format t "Tangent points ~E ~E~%" (visc-compute visc) (visc-compute (+ visc dx)))
;;                      (let ((new-visc (+ visc (/ (- (visc-compute visc) visc) tangent))))
;;                        (setf err (/ (abs (- new-visc visc)) visc))
;;                        (format t "Visc iter ~D ~E ~E - err ~E~%" i new-visc visc err)
;;                        (setf
;;                         visc-prev visc)
;;                        (setf
;;                         visc new-visc))))))
;;       visc)))


;; (defparameter *mp* (aref (cl-mpm:sim-mps *sim*) 100))

(defmethod cl-mpm/particle::constitutive-model ((mp particle-finite-viscoelastic-ice) strain dt)
  (with-accessors ((de mp-elastic-matrix)
                   (e mp-e)
                   (nu mp-nu)
                   (stress mp-stress)
                   (def mp-deformation-gradient)
                   (def-n mp-deformation-gradient-0)
                   (strain mp-strain)
                   (strain-n mp-strain-n)
                   (true-visc mp-viscosity)
                   (visc-factor mp-visc-factor)
                   (visc-power mp-visc-power)
                   (enable-viscosity mp-enable-viscosity)
                   (p-mod mp-p-modulus)
                   (p-wave-0 mp-p-modulus-0)
                   )
      mp
    (cl-mpm/constitutive:linear-elastic-mat strain de stress)
    (let* ((visc-n (glen-visco
                    (cl-mpm/fastmaths:fast-scale
                     (cl-mpm/constitutive:linear-elastic-mat strain-n de)
                     (/ 1d0 (cl-mpm/fastmaths:det def-n)))
                    visc-factor
                    visc-power))
           (visc-n1 (glen-visco
                     (cl-mpm/fastmaths:fast-scale
                      stress
                      (/ 1d0 (cl-mpm/fastmaths:det def)))
                     visc-factor
                     visc-power))
           ;; (visc (expt (+ (expt visc-n -1) (expt visc-n1 -1)) -1))
           ;; (visc visc-n)
           (visc (* 0.5d0 (+ visc-n visc-n1)))
           ;; (visc (visco-iterator mp dt))
           ;; (visc 1d12)
           )
      ;; (pprint visc-n)
      ;; (pprint visc-n1)
      ;; (break)
      ;; (break)
      (when enable-viscosity
        (setf visc (visco-iterator-midpoint mp dt))
        ;; (setf visc (visco-iterator mp dt))
        )
      ;; (pprint visc)
      (setf true-visc visc)
      (if (and enable-viscosity
               (not (= dt 0d0)))
          (progn
            (cl-mpm/models/visco::dev-exp-v stress strain-n strain e nu de visc dt))
          (cl-mpm/constitutive:linear-elastic-mat strain de stress))
      (when (and enable-viscosity
               (not (= dt 0d0)))
          (let* ((K (/ e (* 3 (- 1d0 (* 2 nu)))))
                 (G (/ e (* 2 (+ 1d0 nu))))
                 (rho (/ visc e))
                 (dy (/ dt rho))
                 (exp-rho (exp (- dy)))
                 ;; (lam (/ (- 1 exp-rho) dy))
                 (lam 1d0)
                )
            (declare (double-float K G p-mod))
            (setf p-wave-0 (* (+ K (* 4/3 G lam))))
            )
          ;; (cl-mpm/particle::update-p-modulus mp)
          )
      )
    ;; (cl-mpm/particle::update-log-p-wave mp)
    stress))

;; (let ((test (cl-mpm/utils:voigt-from-list (list 0d0 2d0 3d0 2d0 0d0 3d0))))
;;   (pprint (cl-mpm/utils:deviatoric-voigt test))
;;   (pprint (magicl:@ (magicl:.- (magicl:eye 6)
;;                                (magicl:scale
;;                                 (magicl:from-list
;;                                  (list 1d0 1d0 1d0 0d0 0d0 0d0
;;                                        1d0 1d0 1d0 0d0 0d0 0d0
;;                                        1d0 1d0 1d0 0d0 0d0 0d0
;;                                        0d0 0d0 0d0 0d0 0d0 0d0
;;                                        0d0 0d0 0d0 0d0 0d0 0d0
;;                                        0d0 0d0 0d0 0d0 0d0 0d0
;;                                        )
;;                                  '(6 6)) (/ 1d0 3d0))
;;                                ) test))
;;   )
