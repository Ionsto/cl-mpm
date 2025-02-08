(defpackage :cl-mpm/ext
  (:use :cl
        :cffi
   )
  (:import-from
   :magicl tref .+ .-)
  (:export
   #:kirchoff-update
   #:kirchoff-expt-step
   #:kirchoff-expt-step-lisp
   ))
(in-package :cl-mpm/ext)
(declaim (optimize (debug 0) (safety 0) (speed 3)))
;; (declaim (optimize (debug 3) (safety 3) (speed 0)))


(pushnew (asdf:system-relative-pathname 'cl-mpm "./libs/")  *foreign-library-directories* :test #'equal)
(declaim (ftype (function (magicl:matrix/double-float magicl:matrix/double-float)
                          (values))
                kirchoff-expt-step)
         (notinline kirchoff-expt-step)
         )
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)
(defun kirchoff-expt-step-lisp (strain df)
  ;; (loop for s across (cl-mpm/utils::fast-storage strain)
  ;;       do (when (or
  ;;                 (sb-ext:float-nan-p s)
  ;;                 (< (the double-float s) -1d2)
  ;;                 (> (the double-float s) 1d2))
  ;;            (error "Bad strain! ~A" strain)))
  ;; (loop for s across (cl-mpm/utils::fast-storage df)
  ;;       do (when (or
  ;;                 (sb-ext:float-nan-p s)
  ;;                 (< (the double-float s) -1d2)
  ;;                 (> (the double-float s) 1d2))
  ;;            (error "Bad df! ~A" df)))
  (multiple-value-bind (l v) (cl-mpm/utils::eig
                              (cl-mpm/utils:voigt-to-matrix strain))
    ;; (loop for eigenvalue in l
    ;;       do (when (or
    ;;                 (sb-ext:float-nan-p eigenvalue)
    ;;                 ;; (< (the double-float eigenvalue) -1d2)
    ;;                 (> (the double-float eigenvalue) 1d2))
    ;;            (error "Bad eigenvalue! ~A - ~A" l v)))
    (let ((trial-lgs (magicl:@ df
                               v
                               (cl-mpm/utils::matrix-from-list
                                (list
                                 (the double-float (exp (* 2d0 (the double-float (nth 0 l))))) 0d0 0d0
                                 0d0 (the double-float (exp (* 2d0 (the double-float (nth 1 l))))) 0d0
                                 0d0 0d0 (the double-float (exp (* 2d0 (the double-float (nth 2 l)))))
                                 ))
                               (magicl:transpose v)
                               (magicl:transpose df))))

      ;;Enforce symmetry
      (multiple-value-bind (lf vf)
          (cl-mpm/utils::eig trial-lgs)
        (loop for eigenvalue in lf
              do (when (<= (the double-float eigenvalue) 0d0)
                   (error "Negative eigenvalue! ~A - ~A" lf vf)))
        ;; (loop for i from 0 to 2
        ;;       do (setf (nth i lf) (max 1d-20 (nth i lf))))
        (destructuring-bind (lf0 lf1 lf2) lf
          (declare (double-float lf0 lf1 lf2)
                   )
          (setf lf0 (max 1d-20 lf0)
                lf1 (max 1d-20 lf1)
                lf2 (max 1d-20 lf2))
          (cl-mpm/utils:voigt-copy-into
                          (magicl:scale!
                           ;;Note that this is taking care of the shear scaling factor
                           (cl-mpm/utils:matrix-to-voigt
                            (magicl:@
                             vf
                             (cl-mpm/utils::matrix-from-list
                              (list
                               (the double-float (log (the double-float lf0))) 0d0 0d0
                               0d0 (the double-float (log (the double-float lf1))) 0d0
                               0d0 0d0 (the double-float (log (the double-float lf2))))
                              )
                             (magicl:transpose vf)))
                           0.5d0)
                          strain
                          ))))))

;; (sb-ext:restrict-compiler-policy 'speed  3 3)
;; (sb-ext:restrict-compiler-policy 'debug  0 0)
;; (sb-ext:restrict-compiler-policy 'safety 0 0)

(handler-case
    (progn
      (define-foreign-library cl-mpm-cpp
        (:unix "kirchoff.so")
        (t (:default "kirchoff")))

      (use-foreign-library cl-mpm-cpp)


      (defcfun "Kirchoff_Strain_Update" :bool
        (strain-ptr :pointer)
        (df-ptr :pointer))


      (defcfun "CppDruckerPrager" :bool
        (strain-ptr :pointer)
        (E :double)
        (nu :double)
        (phi :double)
        (psi :double)
        (c :double))

      (defcfun "CppMohrCoulomb" :bool
        (strain-ptr :pointer)
        (f-ptr :pointer)
        (ps-ptr :pointer)
        (E :double)
        (nu :double)
        (phi :double)
        (psi :double)
        (c :double))

      (defcfun "MatrixSqrt" :bool
        (input-ptr :pointer)
        (output-ptr :pointer))

      (defcfun "test" :void
        (flags :pointer))
      (format t "~&Using accelerated kirchoff update~%")
      (defun kirchoff-expt-step (strain-matrix df-matrix)
        (magicl.cffi-types:with-array-pointers ((sp (magicl::matrix/double-float-storage strain-matrix))
                                                (dfp (magicl::matrix/double-float-storage df-matrix)))
          (unless (kirchoff-strain-update sp dfp)
            (error "Kirchoff strain update failed - eigensolver error ~A ~A" strain-matrix df-matrix)))
        (values))

      (defun constitutive-drucker-prager (strain de E nu phi psi c)
        (declare (double-float E nu phi psi c))
        (let ((str (cl-mpm/utils::voigt-copy strain)))
          (magicl.cffi-types:with-array-pointers ((sp (cl-mpm/utils:fast-storage str)))
            (unless (CppDruckerPrager sp E nu phi psi c)
              (error "Drucker-Prager failed")))
          (values (magicl:@ de str) str 0d0 t)))
      (defun constitutive-mohr-coulomb (stress de strain E nu phi psi c)
        "Mohr-coulomb, in-place update strain, return a new stress, yield function and ps inc"
        (declare (double-float E nu phi psi c))
        (if t;(> (cl-mpm/constitutive::fast-mc stress phi c) 0d0)
            (let ((str
                    strain
                    ;; (cl-mpm/utils::voigt-copy strain)
                    )
                  )
              (static-vectors:with-static-vector (f-arr 1 :element-type 'double-float)
                (static-vectors:with-static-vector (ps-arr 1 :element-type 'double-float)
                  (magicl.cffi-types:with-array-pointers ((sp (cl-mpm/utils:fast-storage str))
                                                          (f-arr-p f-arr)
                                                          (ps-arr-p ps-arr))
                    (if (CppMohrCoulomb sp f-arr-p ps-arr-p E nu phi psi c)
                        (values (cl-mpm/fastmaths::fast-@-tensor-voigt de str stress)
                                str
                                (aref f-arr 0)
                                (aref ps-arr 0))
                        (values stress strain (aref f-arr 0) 0d0))
                    )))
              )
            (values stress strain 0d0 0d0)
            ))
      (defun matrix-sqrt (mat)
        (let ((output (cl-mpm/utils:matrix-zeros)))
          (magicl.cffi-types:with-array-pointers ((sp (cl-mpm/utils:fast-storage mat))
                                                  (sp-out (cl-mpm/utils:fast-storage output)))
            (MatrixSqrt sp sp-out))
          output))
      )

    (cffi::load-foreign-library-error (c)
      (progn
        (format t "~&No accelerated kirchoff update~%")
        ;;Use a lisp fallback
        (defun kirchoff-expt-step (strain df)
          (kirchoff-expt-step-lisp strain df))
        (defun constitutive-drucker-prager (stress de strain E nu phi psi c)
          (cl-mpm/constitutive::mc-plastic stress
                                           de
                                           strain
                                           E
                                           nu
                                           phi
                                           psi
                                           c))
        (defun matrix-sqrt (mat)
          (multiple-value-bind (l v) (cl-mpm/utils::eig mat)
            (magicl:@
             v
             (cl-mpm/utils::matrix-from-list
              (list (the double-float (sqrt (the double-float (nth 0 l)))) 0d0 0d0
                    0d0 (the double-float (sqrt (the double-float (nth 1 l)))) 0d0
                    0d0 0d0 (the double-float (sqrt (the double-float (nth 2 l))))))
             (magicl:transpose v))
            )))))

(declaim (ftype (function (magicl:matrix/double-float magicl:matrix/double-float) (values)) kirchoff-update))
(defun kirchoff-update (strain df)
  (kirchoff-expt-step strain df)
  ;; (kirchoff-expt-step-lisp strain df)
  ;; (kirchoff-svd-lisp strain df)
  )


;; (defun test-drucker-prager ()
;;   (let* ((E 1d0)
;;          (nu 0.1d0)
;;          (phi 0.1d0)
;;          (psi 0.0d0)
;;          (c 1d0)
;;          (de (cl-mpm/constitutive::linear-elastic-matrix E nu))
;;          (strain (cl-mpm/constitutive::swizzle-coombs->voigt (cl-mpm/utils:voigt-from-list (list -1d0 -1d0 -1d0 4d0 0d0 0d0 ))))
;;          )
;;     ;; (time
;;     ;;  (lparallel:pdotimes (i 1000000)
;;     ;;    (progn)))
;;     (constitutive-drucker-prager strain de E nu phi psi c)
;;     (multiple-value-bind (sig eps f) (constitutive-drucker-prager strain de E nu phi psi c)
;;       ;; (pprint sig)
;;       (pprint (cl-mpm/constitutive::swizzle-voigt->coombs eps)))))
 
;; (defun test-mc ()
;;   (let* ((E 1d0)
;;          (nu 0.2d0)
;;          (phi 0.1d0)
;;          (psi 0.0d0)
;;          (c 1d0)
;;          (de (cl-mpm/constitutive::linear-elastic-matrix E nu))
;;          (strain (cl-mpm/utils:voigt-from-list (list 1d0 2d0 3d0 4d0 5d0 6d0)))
;;          (stress (magicl:@ de strain))
;;          )
;;     ;; (pprint stress)
;;     ;; (pprint de)
;;     (time
;;      (lparallel:pdotimes (i 1000000)
;;        (progn
;;          ;; (cl-mpm/utils:voigt-copy-into strain-0 strain)
;;          (cl-mpm/constitutive::mc-plastic stress de strain de E nu phi psi c))))
;;     (multiple-value-bind (sig eps f) (cl-mpm/constitutive::mc-plastic stress de strain E nu phi psi c)
;;       (pprint f)
;;       (pprint sig)
;;       (pprint eps))))


