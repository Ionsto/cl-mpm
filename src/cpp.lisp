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
  (declare (optimize (debug 3) (speed 0)))
  (loop for s across (cl-mpm/utils::fast-storage strain)
        do (when (or
                  (sb-ext:float-nan-p s)
                  (< (the double-float s) -1d2)
                  (> (the double-float s) 1d2))
             (error "Bad strain! ~A" strain)))
  (loop for s across (cl-mpm/utils::fast-storage df)
        do (when (or
                  (sb-ext:float-nan-p s)
                  (< (the double-float s) -1d2)
                  (> (the double-float s) 1d2))
             (error "Bad df! ~A" df)))
  (multiple-value-bind (l v) (cl-mpm/utils::eig
                              (cl-mpm/utils:voigt-to-matrix strain))
    (loop for eigenvalue in l
          do (when (or
                    (sb-ext:float-nan-p eigenvalue)
                    ;; (< (the double-float eigenvalue) -1d2)
                    (> (the double-float eigenvalue) 1d2))
               (error "Bad eigenvalue! ~A - ~A" l v)))
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

      ;; (magicl:scale! (magicl:.+ trial-lgs (magicl:transpose trial-lgs)) 0.5d0)
      ;;Enforce symmetry
      ;; (setf (magicl:tref trial-lgs 0 1) (magicl:tref trial-lgs 1 0)
      ;;       (magicl:tref trial-lgs 0 2) (magicl:tref trial-lgs 2 0)
      ;;       (magicl:tref trial-lgs 1 2) (magicl:tref trial-lgs 2 1))
      (multiple-value-bind (lf vf)
          (cl-mpm/utils::eig trial-lgs)
        ;; (loop for i from 0 to 2
        ;;       do (setf (nth i lf) (max (nth i lf) 0d0)))
        (loop for eigenvalue in lf
              do (when (<= (the double-float eigenvalue) 0d0)
                   (error "Negative eigenvalue! ~A - ~A" lf vf)))
        (loop for i from 0 to 2
              do (setf (nth i lf) (max 0d0 (nth i lf))))
        (aops:copy-into (magicl::matrix/double-float-storage strain)
                        (magicl::matrix/double-float-storage
                         (magicl:scale!
                          ;;Note that this is taking care of the shear scaling factor
                          (cl-mpm/utils:matrix-to-voigt
                           (magicl:@
                            vf
                            (cl-mpm/utils::matrix-from-list
                             (list
                              (the double-float (log (the double-float (nth 0 lf)))) 0d0 0d0
                              0d0 (the double-float (log (the double-float (nth 1 lf)))) 0d0
                              0d0 0d0 (the double-float (log (the double-float (nth 2 lf))))
                              )
                             )
                            (magicl:transpose vf)))
                          0.5d0)))
        )))
  )

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
      (defcfun "test" :void
        (flags :pointer))
      (format t "~&Using accelerated kirchoff update~%")
      (defun kirchoff-expt-step (strain-matrix df-matrix)
        (magicl.cffi-types:with-array-pointers ((sp (magicl::matrix/double-float-storage strain-matrix))
                                                (dfp (magicl::matrix/double-float-storage df-matrix)))
          (unless (kirchoff-strain-update sp dfp)
            (error "Kirchoff strain update failed - eigensolver error ~A ~A" strain-matrix df-matrix)))
        (values)))

    (cffi::load-foreign-library-error (c)
      (progn
        (format t "~&No accelerated kirchoff update~%")
        ;;Use a lisp fallback
        (defun kirchoff-expt-step (strain df)
          (kirchoff-expt-step-lisp strain df))
        )))
(declaim (ftype (function (magicl:matrix/double-float magicl:matrix/double-float) (values)) kirchoff-update))
(defun kirchoff-update (strain df)
  (kirchoff-expt-step strain df)
  ;; (kirchoff-svd-lisp strain df)
  )

