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


(pushnew (asdf:system-relative-pathname 'cl-mpm "./libs/")  *foreign-library-directories* :test #'equal)
(define-foreign-library cl-mpm-cpp
    (:unix "kirchoff.so")
  (t (:default "kirchoff")))

(use-foreign-library cl-mpm-cpp)


(defcfun "Kirchoff_Strain_Update" :bool
  (strain-ptr :pointer)
  (df-ptr :pointer))
(defcfun "test" :void
  (flags :pointer))


(defun kirchoff-update (strain df)
  (kirchoff-expt-step strain df))
(declaim (ftype (function (magicl:matrix/double-float magicl:matrix/double-float)
                          (values))
                kirchoff-expt-step)
         (inline kirchoff-expt-step)
         )
(defun kirchoff-expt-step (strain-matrix df-matrix)
  (magicl.cffi-types:with-array-pointers ((sp (magicl::matrix/double-float-storage strain-matrix))
                                          (dfp (magicl::matrix/double-float-storage df-matrix)))
    (kirchoff-strain-update sp dfp)
    )
  (values)
  )
(defun kirchoff-expt-step-lisp (strain df)
  (multiple-value-bind (l v) (cl-mpm/utils::eig
                                  (cl-mpm/utils:voigt-to-matrix strain))
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
          (multiple-value-bind (lf vf)
              (cl-mpm/utils::eig (magicl:scale! (magicl:.+ trial-lgs (magicl:transpose trial-lgs)) 0.5d0))
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
                              0.5d0)))))))










