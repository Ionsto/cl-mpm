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
(eval-when
    (:compile-toplevel)
  (push :sb-simd *features*)
  #+:sb-simd (require 'sb-simd)
  )

#+:sb-simd
(progn
  ;(require 'sb-simd)
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
    (simd-accumulate (magicl::matrix/double-float-storage a) (magicl::matrix/double-float-storage b))
    (values)))

(declaim
 (inline mult)
 (ftype (function (magicl:matrix/double-float
                   magicl:matrix/double-float
                   magicl:matrix/double-float) (values)
                  ) mult))
(defun mult (a b res)
  (declare (type magicl:matrix/double-float a b res))
  (let ((a-s (magicl::matrix/double-float-storage a))
        (b-s (magicl::matrix/double-float-storage b))
        (res-s (magicl::matrix/double-float-storage res))
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
  (let ((x-a (magicl::matrix/double-float-storage x)))
    (declare (type (simple-array double-float) x-a))
    (values (- (* (aref x-a 0) (aref x-a 3))
               (* (aref x-a 1) (aref x-a 2))))))

(declaim
 (inline fast-fmacc-array)
 (ftype (function ((simple-array double-float)
                   (simple-array double-float)
                   double-float) (values)) fast-fmacc))
(defun fast-fmacc-array (a b d)
  #+:sb-simd (simd-fmacc a
                         b
                         d)
  #-:sb-simd (setf a (aops:reduce #'+ a (aops:each (lambda (x) (* x d)) b))))

(declaim
 (inline fast-fmacc)
 (ftype (function (magicl:matrix/double-float
                   magicl:matrix/double-float
                   double-float) (values)) fast-fmacc))
(defun fast-fmacc (a b d)
  #+:sb-simd (simd-fmacc (magicl::matrix/double-float-storage a)
                          (magicl::matrix/double-float-storage b)
                          d)
  #-:sb-simd (magicl:.+ a (magicl:scale b d) a))

(declaim
   (inline fast-add)
   (ftype (function (magicl:matrix/double-float magicl:matrix/double-float) (values)) fast-add))
(defun fast-add (a b)
  #+:sb-simd (simd-add a b)
  #-:sb-simd (magicl:.+ a b a)
  )


(declaim (inline voigt-tensor-reduce-lisp)
         (ftype (function (magicl:matrix/double-float) (values double-float)) voigt-tensor-reduce-lisp))
(let ((second-invar (magicl:from-array (make-array 3 :initial-contents '(1d0 1d0 1/4)) '(3 1) :type 'double-float :layout :column-major)))
  (defun voigt-tensor-reduce-lisp (a)
     (values (magicl::sum (magicl:.* a a second-invar)))))

(declaim (inline voigt-tensor-reduce-simd)
         (ftype (function (magicl:matrix/double-float) (values double-float)) voigt-tensor-reduce-simd))
(defun voigt-tensor-reduce-simd (a)
  "Calculate the product A_{ij}A_{ij}"
  (let ((arr (magicl::matrix/double-float-storage a)))
    (declare ((simple-array double-float) arr))
    (values (+
            (* (aref arr 0) (aref arr 0))
            (* (aref arr 1) (aref arr 1))
            (* (aref arr 2) (aref arr 2) 1/4)))))

(defun voigt-tensor-reduce (a)
  #+:sb-simd (voigt-tensor-reduce-simd a)
  #-:sb-simd (voigt-tensor-reduce-lisp a)
  )

(declaim (inline stretch-to-sym)
         (ftype (function (magicl:matrix/double-float &optional magicl:matrix/double-float) (values)) stretch-to-sym))
(defun stretch-to-sym (stretch &optional (result nil))
  (unless result
    (setf result (cl-mpm/utils:voigt-zeros)))
  (progn
    (declaim (magicl:matrix/double-float result))

    (setf (magicl:tref result  0 0)
          (magicl:tref stretch 0 0))

    (setf (magicl:tref result  1 0)
          (magicl:tref stretch 1 1))

    (setf (magicl:tref result  2 0)
          (+ (the double-float (magicl:tref stretch 0 1))
             (the double-float (magicl:tref stretch 1 0)))))
  (values))

(defun stretch-to-skew (stretch &optional (result nil))
  (unless result
    (setf result (cl-mpm/utils:voigt-zeros)))
  (setf (magicl:tref result  0 0)
        0)
  (setf (magicl:tref result  1 0)
        0)
  (setf (magicl:tref result  2 0)
        (- (the double-float (magicl:tref stretch 0 1))
           (the double-float (magicl:tref stretch 1 0)))))

(defun mult-transpose-accumulate (a b scale res)
  "(incf res (scale! (@ (tranpose a) b) scale))"
  (declare (type magicl:matrix/double-float a b res)
           (type double-float scale))
  (declare (optimize (safety 3)))
  (let (
        ;; (a-s (magicl::matrix/double-float-storage a))
        ;; (b-s (magicl::matrix/double-float-storage b))
        ;; (res-s (magicl::matrix/double-float-storage res))
        )
    (loop for i from 0 to 1
          do (loop for j from 0 to 2
                   do (incf (the double-float (magicl:tref res i 0)) (* (the double-float (magicl:tref a j i))
                                                                        (the double-float (magicl:tref b j 0)) scale))))))
