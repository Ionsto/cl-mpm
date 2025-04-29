
(defpackage :cl-mpm/models/visco
  (:use :cl :cl-mpm/utils :cl-mpm/fastmaths)
  (:export
   )
  )
(declaim (optimize (debug 3) (safety 3) (speed 0)))
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
                   (stress mp-stress)
                   (strain mp-strain)
                   (strain-n mp-strain-n)
                   (viscosity mp-viscosity)
                   (enable-viscosity mp-enable-viscosity)
                   )
      mp
    ;;Train elastic strain - plus trail kirchoff stress
    (if enable-viscosity
        (cl-mpm/models/visco::finite-strain-linear-viscous stress strain de e nu dt viscosity)
        ;; (cl-mpm/ext::constitutive-viscoelastic stress strain de e nu dt viscosity)
        (cl-mpm/constitutive:linear-elastic-mat strain de stress))
    stress))

(defun glen-visco (stress visc-factor visc-power)
  (let* ((dev (cl-mpm/utils:deviatoric-voigt stress))
         (effective-stress (sqrt (* 1/2 (cl-mpm/fastmaths::voigt-j2 dev))))
         (visc-factor (expt visc-factor (- visc-power)))
         )
    (/ 1d0
       (+ 1d-20
          (* visc-factor 2 (expt effective-stress (- visc-power 1)))))))

(defmethod constitutive-model ((mp particle-finite-viscoelastic-ice) strain dt)
  (with-accessors ((de mp-elastic-matrix)
                   (e mp-e)
                   (nu mp-nu)
                   (stress mp-stress)
                   (def mp-deformation-gradient)
                   (strain mp-strain)
                   (strain-n mp-strain-n)
                   (true-visc mp-viscosity)
                   (visc-factor mp-visc-factor)
                   (visc-power mp-visc-power)
                   (enable-viscosity mp-enable-viscosity)
                   )
      mp
    (let ((visc (glen-visco (cl-mpm/fastmaths:fast-scale-voigt
                             stress
                             (/ 1d0 (magicl:det def)))
                            visc-factor
                            visc-power)))
      (setf true-visc visc)
      (if enable-viscosity
          (cl-mpm/models/visco::finite-strain-linear-viscous stress strain de e nu dt visc)
          ;; (cl-mpm/ext::constitutive-viscoelastic stress strain de e nu dt visc)
          (cl-mpm/constitutive:linear-elastic-mat strain de stress))
      )
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
             (G (/ E (* 2 (+ 1d0 nu))))
             (identity (cl-mpm/utils:matrix-eye 1d0))
             (d-neq
               (cl-mpm/fastmaths::fast-scale!
                (cl-mpm/utils:matrix-from-list (list
                                                (- 1d0 nu) nu nu
                                                nu (- 1d0 nu) nu
                                                nu nu (- 1d0 nu)))
                (/ E (* (+ 1d0 nu) (- 1d0 (* 2d0 nu))))))
             ;; (d-neq (cl-mpm/fastmaths::fast-scale
             ;;         dev
             ;;         (* 2 G)))
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
          (cl-mpm/constitutive::linear-elastic-mat strain de stress)
          out-strain))))

(defun test ()
  (let* ((E 1d0)
         (nu 0d0)
         (dt 1d0)
         (viscosity 1d0)
         (strain-0 (cl-mpm/utils:voigt-from-list (list 1d0 0d0 0d0 0d0 0d0 0d0)))
         (de (cl-mpm/constitutive:linear-elastic-matrix E nu)))
    (let* ((strain (cl-mpm/utils::deep-copy strain-0))
           (stress (magicl:@ de strain)))
      (pprint (finite-strain-linear-viscous stress strain de e nu dt viscosity))
      (pprint stress))

    (let* ((strain (cl-mpm/utils::deep-copy strain-0))
           (stress (magicl:@ de strain)))
      (pprint (cl-mpm/ext::constitutive-viscoelastic stress de strain e nu dt viscosity))
      (pprint stress))
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
          (cl-mpm/constitutive::linear-elastic-mat strain de stress)
          out-strain))))
