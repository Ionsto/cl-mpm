(defpackage :cl-mpm/ext
  (:use :cl
        :cffi
   )
  (:import-from
   :magicl tref .+ .-)
  (:export
   #:kirchoff-expt-step
   ))
(in-package :cl-mpm/ext)


(define-foreign-library cl-mpm-cpp
  (:unix (:or "./libs/kirchoff.so"))
  (t (:default "kirchoff")))

(use-foreign-library cl-mpm-cpp)


(defcfun "Kirchoff_Strain_Update" :bool
  (strain-ptr :pointer)
  (df-ptr :pointer))
(defcfun "test" :void
  (flags :pointer))
(let ((strain (make-array 6 :element-type 'double-float :initial-element 0d0))
      (df (make-array 9 :element-type 'double-float :initial-contents (list 1d0 0d0 0d0
                                                                            0d0 1d0 0d0
                                                                            0d0 0d0 3d0))))
  (format t "strain: ~A ~%" strain)
  (magicl.cffi-types:with-array-pointers ((sp strain)
                                          (dfp df))
    (kirchoff-strain-update sp dfp)
    )
  (format t "strain: ~A ~%" strain)
  )
(let ((strain (magicl:zeros '(6 1)))
      (df (cl-mpm/utils:matrix-from-list (list 2d0 0d0 0d0
                                               0d0 1d0 0d0
                                               0d0 0d0 1d0))))


  (kirchoff-expt-step strain df)
  (format t "strain:~A~%" strain)
  )


(defun kirchoff-expt-step (strain-matrix df-matrix)
  (let ((sm-s (magicl::storage strain-matrix))
        (dfm-s (magicl::storage df-matrix)))
    ;; (declare ((simple-array double-float *) sm-s dfm-s))
    (magicl.cffi-types:with-array-pointers ((sp (magicl::matrix/double-float-storage strain-matrix))
                                            (dfp (magicl::matrix/double-float-storage df-matrix)))
      ;; (format t "~A ~%" (type-of ))
      ;; (format t "~A ~%" (type-of ))
      ;; (test sp)
      ;; (format t "~A~%" strain-matrix)
      (kirchoff-strain-update sp dfp)
      ))
  (values)
  )
