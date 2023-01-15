(defpackage :cl-mpm/fastmath
  (:use :cl)
  (:import-from
    :magicl tref .+ .-)
  (:export
   #:fast-add
   #:fast-fmacc
   #:mult
   #:det
    ))

(declaim (optimize (debug 0) (safety 0) (speed 3)))
(in-package :cl-mpm/fastmath)

(push :sb-simd *features*)
;; #+:sb-simd
(progn
  (require 'sb-simd)
  (declaim
   (inline simd-accumulate)
   (ftype (function ((simple-array double-float) (simple-array double-float)) (values)) simd-accumulate))
  (defun simd-accumulate (a b)
    (declare (type sb-simd:f64vec a b))
    (setf (sb-simd-avx:f64.2-aref a 0)
          (sb-simd-avx:f64.2+
           (sb-simd-avx:f64.2-aref a 0)
           (sb-simd-avx:f64.2-aref b 0)
           ))
    (values))

  (declaim
   (inline simd-fmacc)
   (ftype (function ((simple-array double-float) (simple-array double-float) double-float) (values)) simd-fmacc))
  (defun simd-fmacc (a b scale)
    (declare (type sb-simd:f64vec a b)
             (type double-float scale))
    (setf (sb-simd-avx:f64.2-aref a 0)
          (sb-simd-avx:f64.2+
           (sb-simd-avx:f64.2-aref a 0)
           (sb-simd-avx:f64.2*
            (sb-simd-avx:f64.2-aref b 0)
            (sb-simd-avx:f64.2 scale))
           ))
    (values))
  (declaim
   (inline simd-add)
   (ftype (function (magicl:matrix/double-float magicl:matrix/double-float) (values)) simd-add))
  (defun simd-add (a b)
    (simd-accumulate (magicl::storage a) (magicl::storage b))
    (values)))

(declaim
 (inline mult)
 (ftype (function (magicl:matrix/double-float
                   magicl:matrix/double-float
                   magicl:matrix/double-float) (values)
                  ) mult))
(defun mult (a b res)
  (declare (type magicl:matrix/double-float a b res))
  (let ((a-s (magicl::storage a))
        (b-s (magicl::storage b))
        (res-s (magicl::storage res))
        )
    (declare (type (simple-array double-float) a-s b-s res-s))
    (loop for i from 0 to 2
          do (loop for j from 0 to 1
                   do (incf (aref res-s i) (* (aref a-s (+ j (* 2 i)))
                                              (aref b-s j)))))))

(declaim (inline det)
         (ftype (function (magicl:matrix/double-float) (values double-float)) det)
         )
(defun det (x)
  (let ((x-a (magicl::storage x)))
    (declare (type (simple-array double-float) x-a))
    (values (- (* (aref x-a 0) (aref x-a 3))
               (* (aref x-a 1) (aref x-a 2))))))


(declaim
 (inline fast-fmacc)
 (ftype (function (magicl:matrix/double-float
                   magicl:matrix/double-float
                   double-float) (values)) fast-fmacc))
(defun fast-fmacc (a b d)
  #+:sb-simd (simd-fmacc (magicl::storage a)
                          (magicl::storage b)
                          d)
  #-:sb-simd (magicl:.+ a (magicl:scale b d) a))

(declaim
   (inline fast-add)
   (ftype (function (magicl:matrix/double-float magicl:matrix/double-float) (values)) fast-add))
(defun fast-add (a b)
  #+:sb-simd (simd-add a b)
  #-:sb-simd (magicl:.+ a b a)
  )
