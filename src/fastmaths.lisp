(defpackage :cl-mpm/fastmaths
  (:use :cl)
  (:import-from
   :magicl tref .+ .-
   )
  (:import-from
   :cl-mpm/utils varef)
  (:export
   #:fast-fmacc
   #:mult
   #:det
   #:det-3x3
   #:dot
   #:norm
   #:fast-.+
   #:fast-.-
   #:fast-.*
   #:fast-scale-vector
   #:fast-scale-voigt
   #:fast-scale!
   #:fast-scale
   #:fast-zero
   #:fast-sum
   #:mag
   ))

(declaim (optimize (debug 0) (safety 0) (speed 3)))
;; (declaim (optimize (debug 3) (safety 3) (speed 0)))
(in-package :cl-mpm/fastmaths)

(require 'sb-simd)
(pushnew :sb-simd *features*)

(eval-when
    (:compile-toplevel)
  (pushnew :sb-simd *features*)
  (require 'sb-simd)
  ;; #+:sb-simd (require 'sb-simd)
  )
#+:sb-simd (format t "~%Built with sb-simd~%")

#+:sb-simd
(progn
  (declaim
   (inline simd-accumulate)
   (ftype (function ((simple-array double-float)
                     (simple-array double-float)) (values)) simd-accumulate))
  (defun simd-accumulate (a b)
    ;; (declare (type (simple-array double-float 3) a b))
    ;; (declare (type sb-simd:f64vec a b))

    (setf (sb-simd-avx:f64.2-aref a 0)
          (sb-simd-avx:f64.2+
           (sb-simd-avx:f64.2-aref a 0)
           (sb-simd-avx:f64.2-aref b 0)
           ))
    (incf (aref a 2) (aref b 2))

    ;; (loop for i from 0 to 2
    ;;       do (incf (aref a i) (aref b i)))
    (values))

  (declaim
   (inline simd-fmacc)
   (ftype (function ((simple-array double-float (3)) (simple-array double-float (3)) double-float) (values)) simd-fmacc))
  (defun simd-fmacc (result source scale)
    (declare
     ((simple-array double-float (3)) result source)
                                        ;(type sb-simd:f64vec result source)
     (type double-float scale))
    (setf (sb-simd-avx:f64.2-aref result 0)
          (sb-simd-avx:f64.2+
           (sb-simd-avx:f64.2-aref result 0)
           (sb-simd-avx:f64.2*
            (sb-simd-avx:f64.2-aref source 0)
            scale
            ;; (sb-simd-avx:f64.2 scale)
            )
           ))
    (incf (aref result 2) (* scale (aref source 2)))

    ;; (loop for i from 0 to 2
    ;;       do (incf (aref result i) (* scale (aref source i))))
    (values))
  (declaim
   (inline simd-add)
   (ftype (function (magicl:matrix/double-float magicl:matrix/double-float) (values)) simd-add))
  (defun simd-add (a b)
    (simd-accumulate (magicl::matrix/double-float-storage a)
                     (magicl::matrix/double-float-storage b))
    (values))

  (defun test-simd-accumulate ()
    (let* ((a (make-array 3 :element-type 'double-float  :initial-contents (list 1d0 2d0 3d0)))
           (b (make-array 3 :element-type 'double-float  :initial-contents (list 3d0 5d0 9d0))))
      (loop for i from 0 to 2
            do (incf (aref a i) (aref b i)))
      (pprint a)
      )
    (let* ((a (make-array 3 :element-type 'double-float  :initial-contents (list 1d0 2d0 3d0)))
           (b (make-array 3 :element-type 'double-float  :initial-contents (list 3d0 5d0 9d0))))
      (simd-accumulate a b)
      (pprint a)
      ))

  (defun test-simd-fmacc ()
    (let ((scale 2d0))
      (let* ((a (make-array 3 :element-type 'double-float  :initial-contents (list 1d0 2d0 3d0)))
             (b (make-array 3 :element-type 'double-float  :initial-contents (list 3d0 5d0 9d0))))
        (loop for i from 0 to 2
              do (incf (aref a i) (* scale(aref b i))))
        (pprint a)
        )
      (let* ((a (make-array 3 :element-type 'double-float  :initial-contents (list 1d0 2d0 3d0)))
             (b (make-array 3 :element-type 'double-float  :initial-contents (list 3d0 5d0 9d0))))
        (simd-fmacc a b scale)
        (pprint a)
        ))))

(declaim
 (inline simd-diff-norm)
 (ftype (function ((simple-array double-float (3))
                   (simple-array double-float (3))) double-float) simd-diff-norm))
(defun simd-diff-norm (a b)
  (declare (type (simple-array double-float (3)) a b))
  (declare (optimize (speed 3) (safety 0)))
  (+
   (the double-float (expt (- (aref a 0) (aref b 0)) 2))
   (the double-float (expt (- (aref a 1) (aref b 1)) 2))
   (the double-float (expt (- (aref a 2) (aref b 2)) 2))))

(defun diff-norm (a b)
  (simd-diff-norm (cl-mpm/utils:fast-storage a)
                  (cl-mpm/utils:fast-storage b)))

(defun diff-mag (a b)
  (sqrt
   (simd-diff-norm (cl-mpm/utils:fast-storage a)
                   (cl-mpm/utils:fast-storage b))))

(declaim
 (inline lisp-fmacc)
 (ftype (function ((simple-array double-float)
                   (simple-array double-float)
                   double-float)
                  (values)) lisp-fmacc))
(defun lisp-fmacc (a b d)
  (declare ((simple-array double-float (*)) a b)
           (double-float d)
           )
  (loop for i fixnum from 0 below (length a)
        do (setf (aref a i) (+ (aref a i) (* (aref b i) d))))
  (values))

(declaim
 (inline fast-fmacc-array)
 (ftype (function ((simple-array double-float)
                   (simple-array double-float)
                   double-float) (values)) fast-fmacc-array))
(defun fast-fmacc-array (a b d)
  #+:sb-simd (simd-fmacc a
                         b
                         d)
  #-:sb-simd (lisp-fmacc a b d))

(declaim
 (inline fast-fmacc)
 (ftype (function (magicl:matrix/double-float
                   magicl:matrix/double-float
                   double-float) (values)) fast-fmacc))
(defun fast-fmacc (a b d)
  #+:sb-simd (simd-fmacc (magicl::matrix/double-float-storage a)
                          (magicl::matrix/double-float-storage b)
                          d)
  ;#-:sb-simd (fast-.+ a (magicl:scale b d) a)
  #-:sb-simd (lisp-fmacc (magicl::matrix/double-float-storage a)
                         (magicl::matrix/double-float-storage b)
                         d))

#+:sb-simd
(progn
  (declaim
   (inline simd-any+)
   (ftype (function ((simple-array double-float)
                     (simple-array double-float)
                     (simple-array double-float))
                    (values)) simd-any+))
  (defun simd-any+ (a b target)
    ;; (declare (optimize (speed 0) (debug 0) (safety 3)))
    (declare ((simple-array double-float (*)) a b target))
    (let ((offset 0))
      (declare (type sb-simd:f64vec a b target)
               (fixnum offset))
      (multiple-value-bind (iter remain) (floor (length a) 2)
        (declare (fixnum iter remain))
        (dotimes (i iter)
          (setf (sb-simd-avx:f64.2-aref target offset)
                (sb-simd-avx:f64.2+
                 (sb-simd-avx:f64.2-aref a offset)
                 (sb-simd-avx:f64.2-aref b offset)))
          (incf offset 2))
        (unless (eq remain 0)
          (dotimes (i remain)
            (setf (aref target offset)
                  (+ (aref a offset) (aref b offset)))
            (incf offset 1)))))
    target))

(declaim
 (inline lisp-any+)
 (ftype (function ((simple-array double-float)
                   (simple-array double-float)
                   (simple-array double-float))
                  (values)) lisp-any+))
(defun lisp-any+ (a b target)
  (declare ((simple-array double-float (*)) a b target))
  (loop for i fixnum from 0 below (length a)
        do (setf (aref target i) (+ (aref a i) (aref b i))))
  target)



#+:sb-simd
(progn
  (declaim
   (inline simd-any-)
   (ftype (function ((simple-array double-float)
                     (simple-array double-float)
                     (simple-array double-float))
                    (values)) simd- any-))
  (defun simd-any- (a b target)
    (declare ((simple-array double-float (*)) a b target))
    (let ((offset 0))
      (declare (type sb-simd:f64vec a b target)
               (fixnum offset))
      (multiple-value-bind (iter remain) (floor (length a) 2)
        (declare (fixnum iter remain))
        (dotimes (i iter)
          (setf (sb-simd-avx:f64.2-aref target offset)
                (sb-simd-avx:f64.2-
                 (sb-simd-avx:f64.2-aref a offset)
                 (sb-simd-avx:f64.2-aref b offset)))
          (incf offset 2))
        (unless (eq remain 0)
          (dotimes (i remain)
            (setf (aref target offset)
                  (- (aref a offset) (aref b offset)))
            (incf offset 1)))))
    target))

(declaim
 (inline lisp-any-)
 (ftype (function ((simple-array double-float)
                   (simple-array double-float)
                   (simple-array double-float))
                  (values)) lisp-any-))
(defun lisp-any- (a b target)
  (declare ((simple-array double-float (*)) a b target))
  (loop for i fixnum from 0 below (length a)
        do (setf (aref target i) (- (aref a i) (aref b i))))
  target)



#+:sb-simd
(progn
  (defun simd-any* (a b target)
                                        ;(declare (optimize (speed 0) (debug 0) (safety 3)))
    (declare ((simple-array double-float (*)) a b target))
    (let ((offset 0))
      (declare (type sb-simd:f64vec a b target)
               (fixnum offset))
      (multiple-value-bind (iter remain) (floor (length a) 2)
        (declare (fixnum iter remain))
        (dotimes (i iter)
          (setf (sb-simd-avx:f64.2-aref target offset)
                (sb-simd-avx:f64.2*
                 (sb-simd-avx:f64.2-aref a offset)
                 (sb-simd-avx:f64.2-aref b offset)))
          (incf offset 2))
        (unless (eq remain 0)
          (dotimes (i remain)
            (setf (aref target offset)
                  (* (aref a offset) (aref b offset)))
            (incf offset 1))
          )
        ))

    ;; (loop for i from 0 below (length a)
    ;;       do (setf (aref target i) (* (aref a i) (aref b i))))
    target))


(declaim
 (inline lisp-any*)
 (ftype (function ((simple-array double-float)
                   (simple-array double-float)
                   (simple-array double-float))
                  (values)) lisp-any*))
(defun lisp-any* (a b target)
  (declare ((simple-array double-float (*)) a b target))
  (loop for i fixnum from 0 below (length a)
        do (setf (aref target i) (* (aref a i) (aref b i))))
  target)

#+:sb-simd
(progn
  (defun test-simd-any+ ()
    (let ((a (cl-mpm/utils::vector-from-list (list 1d0 2d0 3d0)))
          (b (cl-mpm/utils::vector-from-list (list 3d0 5d0 9d0))))
      (let ((res (cl-mpm/utils::vector-zeros)))
        (magicl:.+ a b res)
        (pprint res))
      (let ((res (cl-mpm/utils::vector-zeros)))
        (simd-any+ (magicl::storage a)
                   (magicl::storage b)
                   (magicl::storage res))
        (pprint res)))
    )
  (defun test-simd-any* ()
    (let ((a (cl-mpm/utils::vector-from-list (list 1d0 2d0 3d0)))
          (b (cl-mpm/utils::vector-from-list (list 3d0 5d0 9d0))))
      (let ((res (cl-mpm/utils::vector-zeros)))
        (magicl:.* a b res)
        (pprint res))
      (let ((res (cl-mpm/utils::vector-zeros)))
        (simd-any* (magicl::storage a)
                   (magicl::storage b)
                   (magicl::storage res))
        (pprint res)))
    )

  (declaim
   (inline simd-any+-4)
   (ftype (function ((simple-array double-float)
                     (simple-array double-float)
                     (simple-array double-float))
                    (values)) simd-any+-4))
  (defun simd-any+-4 (a b target)
    (declare (optimize (speed 0) (debug 0) (safety 3)))
    (declare ((simple-array double-float (*)) a b target))
    (let ((offset 0))
      (declare (type sb-simd:f64vec a b target)
               (fixnum offset))
      (multiple-value-bind (iter remain) (floor (length a) 4)
        (declare (fixnum iter remain))
        (dotimes (i iter)
          (setf (sb-simd-avx:f64.4-aref target offset)
                (sb-simd-avx:f64.4+
                 (sb-simd-avx:f64.4-aref a offset)
                 (sb-simd-avx:f64.4-aref b offset)))
          (incf offset 4))
        (unless (eq remain 0)
          (dotimes (i remain)
            (setf (aref target offset)
                  (+ (aref a offset) (aref b offset)))
            (incf offset 1)))))
    target)

  (defun test-simd-any+-4 ()
    (let ((a (cl-mpm/utils::voigt-from-list (list 1d0 2d0 3d0 9d0 2d0 10d0)))
          (b (cl-mpm/utils::voigt-from-list (list 3d0 5d0 9d0 -1d0 3d0 1d0))))
      (let ((res (cl-mpm/utils::voigt-zeros)))
        (magicl:.+ a b res)
        (pprint res))
      (let ((res (cl-mpm/utils::voigt-zeros)))
        (simd-any+ (magicl::storage a)
                   (magicl::storage b)
                   (magicl::storage res))
        (pprint res)))
    ))

(defun @-m-v (matrix vector result-vector)
  "Multiply a 3x3 matrix with a 3x1 vector to calculate a 3x1 vector in place"
  (declare (magicl:matrix/double-float matrix vector result-vector)
           ;; (optimize (speed 0) (safety 3) (debug 0))
           )
  (let ((a (magicl::matrix/double-float-storage matrix))
        (b (magicl::matrix/double-float-storage vector))
        (c (magicl::matrix/double-float-storage result-vector)))
    (declare ((simple-array double-float (9)) a)
             ((simple-array double-float (3)) b c))
    (policy-cond:with-expectations (> speed safety)
        ((assertion (eq (magicl::matrix/double-float-layout matrix) :column-major))
         (assertion (= 3 (magicl::matrix/double-float-nrows matrix)))
         (assertion (= 3 (magicl::matrix/double-float-ncols matrix)))
         (assertion (= 3 (length b)))
         (assertion (= 3 (length c))))
      (flet ((tref (m x y)
               (aref m (+ x (* 3 y)))))
        (loop for i fixnum from 0 below 3
              do
                 (setf
                  (aref c i)
                  (+
                   (* (aref b 0) (tref a i 0))
                   (* (aref b 1) (tref a i 1))
                   (* (aref b 2) (tref a i 2))))))))
  result-vector)


#+:sb-simd
(progn
  (declaim
   (inline @-stretch-vec-simd)
   (ftype (function (magicl:matrix/double-float magicl:matrix/double-float magicl:matrix/double-float)
                    magicl:matrix/double-float) @-stretch-vec-simd))
  (defun @-stretch-vec-simd (matrix vector result-vector)
    "Multiply a 9x3 matrix with a 3x1 vector to calculate a 9x1 vector in place - SIMD implementation"
    (declare (magicl:matrix/double-float matrix vector result-vector)
             ;; (optimize (speed 0) (safety 3) (debug 0))
             )
    (let ((a (magicl::matrix/double-float-storage matrix))
          (b (magicl::matrix/double-float-storage vector))
          (c (magicl::matrix/double-float-storage result-vector))
          )
      (declare ((simple-array double-float (27)) a)
               ((simple-array double-float (3)) b)
               ((simple-array double-float (9)) c)
               )
      ;; (flet ((tref (m x y)
      ;;          (aref m (+ (* 9 x)  y))))
      ;;   (loop for i fixnum from 0 below 9
      ;;         do
      ;;            (setf (aref c i) 0d0)
      ;;            (loop for j fixnum from 0 below 3
      ;;                  do (incf (aref c i) (the double-float (* (aref b j) (tref a j i)))))
      ;;         ))
      (macrolet ((simd-component (i)
                   (declare (fixnum i))
                   `(setf
                     (sb-simd-avx:f64.4-aref c ,(the fixnum (* 4 i)))
                     (sb-simd-avx:f64.4+
                      (sb-simd-avx:f64.4*
                       (sb-simd-avx:f64.4-aref a ,(the fixnum (+ (* i 4) (* 9 0))))
                       (aref b 0))
                      (sb-simd-avx:f64.4*
                       (sb-simd-avx:f64.4-aref a ,(the fixnum (+ (* i 4) (* 9 1))))
                       (aref b 1))
                      (sb-simd-avx:f64.4*
                       (sb-simd-avx:f64.4-aref a ,(the fixnum (+ (* i 4) (* 9 2))))
                       (aref b 2))))
                   ))
        (simd-component 0)
        (simd-component 1)
        (setf
         (aref c 8)
         (+
          (* (aref a (+ 8 0)) (aref b 0))
          (* (aref a (+ 8 9)) (aref b 1))
          (* (aref a (+ 8 18)) (aref b 2))))
        )
      )
    result-vector))
(progn
  (declaim
   (inline @-stretch-vec-lisp)
   (ftype (function (magicl:matrix/double-float magicl:matrix/double-float magicl:matrix/double-float)
                    magicl:matrix/double-float) @-stretch-vec-lisp))
  (defun @-stretch-vec-lisp (matrix vector result-vector)
    "Multiply a 9x3 matrix with a 3x1 vector to calculate a 9x1 vector in place - LISP implementation"
    (declare (magicl:matrix/double-float matrix vector result-vector))
    (let ((a (magicl::matrix/double-float-storage matrix))
          (b (magicl::matrix/double-float-storage vector))
          (c (magicl::matrix/double-float-storage result-vector))
          )
      (declare ((simple-array double-float (27)) a)
               ((simple-array double-float (3)) b)
               ((simple-array double-float (9)) c)
               )
      (dotimes (i 9)
        (setf
         (aref c i)
         (+
          (* (aref a (+ i 0))  (aref b 0))
          (* (aref a (+ i 9))  (aref b 1))
          (* (aref a (+ i 18)) (aref b 2))))))
    result-vector))

(declaim
 (inline @-stretch-vec)
 (ftype (function (magicl:matrix/double-float magicl:matrix/double-float magicl:matrix/double-float)
                          magicl:matrix/double-float) @-stretch-vec))
(defun @-stretch-vec (matrix vector result-vector)
  (declare (magicl:matrix/double-float matrix vector result-vector))
  #+:sb-simd (@-stretch-vec-simd matrix vector result-vector)
  #-:sb-simd (@-stretch-vec-lisp matrix vector result-vector)
  result-vector)

;; (defun test-@-stretch-vec ()
;;   (let ((stretch (cl-mpm/shape-function::assemble-dstretch-3d (list 1d0 2d0 3d0)))
;;         (vel (cl-mpm/utils::vector-from-list (list 0.1d0 9d0 2d0)))
;;         (res-t (cl-mpm/utils::stretch-dsvp-voigt-zeros))
;;         (res (cl-mpm/utils::stretch-dsvp-voigt-zeros))
;;         )
;;     (magicl:mult stretch vel :target res-t)
;;     (pprint res-t)
;;     (setf res (cl-mpm/utils::stretch-dsvp-voigt-zeros))
;;     (cl-mpm/fastmaths::@-stretch-vec stretch vel res)
;;     (pprint res)
;;     (format t "~%Pass?: ~A~%"
;;      (every #'identity (loop for a across (magicl::storage (magicl:.- res res-t))
;;                              collect (< (abs a) 1d-15))))))

(declaim
 (inline @-dsvp-vec)
 (ftype (function (magicl:matrix/double-float magicl:matrix/double-float double-float magicl:matrix/double-float)
                  magicl:matrix/double-float) @-dsvp-vec))
(defun @-dsvp-vec (matrix vector scale result-vector)
  "Multiply a 3x9 matrix with a 3x1 vector to calculate a 3x1 vector in place"
  (declare (magicl:matrix/double-float matrix vector result-vector)
           (double-float scale)
           ;; (optimize (speed 0) (safety 3) (debug 0))
           )
  (let ((a (magicl::matrix/double-float-storage matrix))
        (b (magicl::matrix/double-float-storage vector))
        (c (magicl::matrix/double-float-storage result-vector))
        )
    (declare ((simple-array double-float (18)) a)
             ((simple-array double-float (3)) c)
             ((simple-array double-float (6)) b)
             )
    (flet ((tref (m x y)
             (aref m (+ (* 6 x) y))))
      (loop for i fixnum from 0 below 3
            do
               (setf (aref c i) 0d0)
               (loop for j fixnum from 0 below 6
                     do (incf (aref c i) (the double-float (* (aref b j) (tref a i j) scale))))
            )))
  result-vector)


;; (declaim
;;  (inline fast-.+)
;;  (ftype (function (magicl:matrix/double-float magicl:matrix/double-float) magicl:matrix/double-float) fast-.+))
;; (defun fast-.+ (a b)
;;   (declare (optimize (speed 0) (space 0)))
;;   (let ((res (cl-mpm/utils::empty-copy a)))
;;     (declare (magicl:matrix/double-float a b res))
;;     (simd-any+ (magicl::matrix/double-float-storage a)
;;                (magicl::matrix/double-float-storage b)
;;                (magicl::matrix/double-float-storage res))
;;     res))

(declaim
 (inline fast-.+)
 (ftype (function (magicl:matrix/double-float magicl:matrix/double-float &optional magicl:matrix/double-float) magicl:matrix/double-float) fast-.+))
(defun fast-.+ (a b &optional res)
  (let ((res (if res
                 res
                 (cl-mpm/utils::empty-copy a))))
    (declare (magicl:matrix/double-float a b res))

    #+:sb-simd
    (simd-any+ (magicl::matrix/double-float-storage a)
               (magicl::matrix/double-float-storage b)
               (magicl::matrix/double-float-storage res))
    #-:sb-simd
    (lisp-any+ (magicl::matrix/double-float-storage a)
               (magicl::matrix/double-float-storage b)
               (magicl::matrix/double-float-storage res))
    res))

(declaim
 (inline fast-.-)
 (ftype (function (magicl:matrix/double-float magicl:matrix/double-float &optional magicl:matrix/double-float) magicl:matrix/double-float) fast-.-))
(defun fast-.- (a b &optional res)
  (let ((res (if res
                 res
                 (cl-mpm/utils::empty-copy a))))
    (declare (magicl:matrix/double-float a b res))
    #+:sb-simd
    (simd-any- (magicl::matrix/double-float-storage a)
               (magicl::matrix/double-float-storage b)
               (magicl::matrix/double-float-storage res))
    #-:sb-simd
    (lisp-any- (magicl::matrix/double-float-storage a)
               (magicl::matrix/double-float-storage b)
               (magicl::matrix/double-float-storage res))
    res))

(macrolet ((def-fast-.--type (name maker length underlying-add)
             (declare (fixnum length))
             `(progn
                (declaim
                 (inline ,name)
                 (ftype (function (magicl:matrix/double-float magicl:matrix/double-float &optional magicl:matrix/double-float) magicl:matrix/double-float) ,name))
                (defun ,name (a b &optional res)
                  (let ((res (if res
                                 res
                                 (cl-mpm/utils::empty-copy a)
                                 ;; (,maker)
                                 )))
                    (declare (magicl:matrix/double-float a b res))
                    (,underlying-add (the (simple-array double-float (,length)) (magicl::matrix/double-float-storage a))
                                     (the (simple-array double-float (,length)) (magicl::matrix/double-float-storage b))
                                     (the (simple-array double-float (,length)) (magicl::matrix/double-float-storage res)))
                    res)
                  ))))
  #+:sb-simd
  (progn
    (def-fast-.--type fast-.--vector cl-mpm/utils:vector-zeros 3 simd-any-)
    (def-fast-.--type fast-.--voigt cl-mpm/utils:voigt-zeros 6 simd-any-)
    (def-fast-.--type fast-.--matrix cl-mpm/utils:matrix-zeros 9 simd-any-))
  #-:sb-simd
  (progn
    (def-fast-.--type fast-.--vector cl-mpm/utils:vector-zeros 3 lisp-any-)
    (def-fast-.--type fast-.--voigt cl-mpm/utils:voigt-zeros   6 lisp-any-)
    (def-fast-.--type fast-.--matrix cl-mpm/utils:matrix-zeros 9 lisp-any-))
  )

(defun fast-.* (a b &optional res)
  (let ((res (if res
                 res
                 (cl-mpm/utils::empty-copy a))))
    (declare (magicl:matrix/double-float a b res))
    #+:sb-simd
    (simd-any* (magicl::matrix/double-float-storage a)
               (magicl::matrix/double-float-storage b)
               (magicl::matrix/double-float-storage res))
    #-:sb-simd
    (lisp-any* (magicl::matrix/double-float-storage a)
               (magicl::matrix/double-float-storage b)
               (magicl::matrix/double-float-storage res))
    res))

;; Here we define some specialised fast adds for different length vectors -> this should almost certianly use the type system to resolve it but who knows
(macrolet ((def-fast-.+-type (name maker length underlying-add)
             (declare (fixnum length))
             `(progn
                (declaim
                 (inline ,name)
                 (ftype (function (magicl:matrix/double-float magicl:matrix/double-float &optional magicl:matrix/double-float) magicl:matrix/double-float) ,name))
                (defun ,name (a b &optional res)
                  (let ((res (if res
                                 res
                                 (cl-mpm/utils::empty-copy a)
                                 ;; (,maker)
                                 )))
                    (declare (magicl:matrix/double-float a b res))
                    (,underlying-add (the (simple-array double-float (,length)) (magicl::matrix/double-float-storage a))
                                     (the (simple-array double-float (,length)) (magicl::matrix/double-float-storage b))
                                     (the (simple-array double-float (,length)) (magicl::matrix/double-float-storage res)))
                    res)
                  ))))

  #+:sb-simd
  (progn
    (def-fast-.+-type fast-.+-vector cl-mpm/utils::vector-zeros 3 simd-any+)
    (def-fast-.+-type fast-.+-voigt cl-mpm/utils::voigt-zeros 6 simd-any+)
    (def-fast-.+-type fast-.+-matrix cl-mpm/utils::matrix-zeros 9 simd-any+-4)
    (def-fast-.+-type fast-.+-stretch cl-mpm/utils::stretch-dsvp-3d-zeros 27 simd-any+-4))
  #-:sb-simd
  (progn
    (def-fast-.+-type fast-.+-vector cl-mpm/utils::vector-zeros           3 lisp-any+)
    (def-fast-.+-type fast-.+-voigt cl-mpm/utils::voigt-zeros             6 lisp-any+)
    (def-fast-.+-type fast-.+-matrix cl-mpm/utils::matrix-zeros           9 lisp-any+)
    (def-fast-.+-type fast-.+-stretch cl-mpm/utils::stretch-dsvp-3d-zeros 27 lisp-any+))
  )

(defun test-.+-vector ()
  (let ((a (cl-mpm/utils::vector-from-list (list 1d0 2d0 3d0)))
        (b (cl-mpm/utils::vector-from-list (list 3d0 5d0 9d0))))
    (let ((res (cl-mpm/utils::vector-zeros)))
      (magicl:.+ a b res)
      (pprint res))
    (let ((res (cl-mpm/utils::vector-zeros)))
      (cl-mpm/fastmaths::fast-.+ a b res)
      (pprint res))
    (let ((res (cl-mpm/utils::vector-zeros)))
      (cl-mpm/fastmaths::fast-.+-vector a b res)
      (pprint res))))
(defun test-.+-voigt ()
  (let ((a (cl-mpm/utils::voigt-from-list (list 1d0 2d0 3d0 2d0 5d0 7d0)))
        (b (cl-mpm/utils::voigt-from-list (list 3d0 5d0 9d0 3d0 4d0 9d0))))
    (let ((res (cl-mpm/utils::voigt-zeros)))
      (magicl:.+ a b res)
      (pprint res))
    (let ((res (cl-mpm/utils::voigt-zeros)))
      (cl-mpm/fastmaths::fast-.+ a b res)
      (pprint res))
    (let ((res (cl-mpm/utils::voigt-zeros)))
      (cl-mpm/fastmaths::fast-.+-voigt a b res)
      (pprint res))))

(defun test-.+-matrix ()
  (let ((a (cl-mpm/utils::matrix-from-list (loop for i from 0 below 9 collect (random 1d0))))
        (b (cl-mpm/utils::matrix-from-list (loop for i from 0 below 9 collect (random 1d0)))))
    (let ((res (cl-mpm/utils::matrix-zeros)))
      (magicl:.+ a b res)
      (pprint res))
    (let ((res (cl-mpm/utils::matrix-zeros)))
      (cl-mpm/fastmaths::fast-.+ a b res)
      (pprint res))
    (let ((res (cl-mpm/utils::matrix-zeros)))
      (cl-mpm/fastmaths::fast-.+-matrix a b res)
      (pprint res))))



;; (declaim
;;  (inline fast-.+-voigt)
;;  (ftype (function (magicl:matrix/double-float magicl:matrix/double-float &optional magicl:matrix/double-float) magicl:matrix/double-float) fast-.+-voigt))
;; (defun fast-.+-voigt (a b &optional res)
;;   (declare (optimize (speed 0) (space 0)))
;;   (let ((res (if res
;;                  res
;;                  (cl-mpm/utils::empty-copy a))))
;;     (declare (magicl:matrix/double-float a b res))
;;     (simd-any+-4 (the (simple-array double-float (6)) (magicl::matrix/double-float-storage a))
;;                  (the (simple-array double-float (6)) (magicl::matrix/double-float-storage b))
;;                  (the (simple-array double-float (6)) (magicl::matrix/double-float-storage res)))
;;     res))




(declaim
 (inline fast-scale!)
 (ftype (function (magicl:matrix/double-float double-float) magicl:matrix/double-float) fast-scale!))
(defun fast-scale! (m scale)
  (let ((m-s (magicl::matrix/double-float-storage m)))
    (loop for i fixnum from 0 below (length m-s)
          do (setf (aref m-s i) (* (aref m-s i) scale))))
  m)

(declaim (ftype (function (magicl:matrix/double-float double-float) magicl:matrix/double-float) fast-scale))
(defun fast-scale (mi scale)
  (declare (double-float scale))
  (let* ((m (cl-mpm/utils::deep-copy mi))
         (m-s (magicl::matrix/double-float-storage m)))
    (loop for i fixnum from 0 below (length m-s)
          do (setf (aref m-s i) (* (aref m-s i) scale)))
    m))

(declaim
 (inline fast-scale-voigt)
 (ftype (function (magicl:matrix/double-float double-float) magicl:matrix/double-float) fast-scale-voigt))
(defun fast-scale-voigt (m scale)
  (fast-scale! (cl-mpm/utils:voigt-copy m) scale))

(defun fast-scale-vector (m scale &optional result)
  (let ((result (if result
                    (cl-mpm/utils:vector-copy-into m result)
                    (cl-mpm/utils:vector-copy m))))
    (fast-scale! result scale)
    result))

(declaim
 (inline fast-zero)
 (ftype (function (magicl:matrix/double-float) magicl:matrix/double-float) fast-zero))
(defun fast-zero (m)
  (let ((m-s (magicl::matrix/double-float-storage m)))
    (loop for i fixnum from 0 below (length m-s)
          do (setf (aref m-s i) 0d0)))
  m)
(declaim
 (inline fast-zero-vector)
 (ftype (function (magicl:matrix/double-float) magicl:matrix/double-float) fast-zero-vector))
(defun fast-zero-vector (m)
  #+:sb-simd
  (let ((m-s (magicl::matrix/double-float-storage m)))
    (declare ((simple-array double-float (3)) m-s))
    (setf (sb-simd-avx:f64.2-aref m-s 0) 0d0)
    (setf (aref m-s 2) 0d0))
  #-:sb-simd (fast-zero m)
  m)
(defun fast-zero-voigt (m)
  #+:sb-simd
  (let ((m-s (magicl::matrix/double-float-storage m)))
    (declare ((simple-array double-float (6)) m-s))
    (setf (sb-simd-avx:f64.4-aref m-s 0) 0d0)
    (setf (sb-simd-avx:f64.2-aref m-s 4) 0d0)
    )
  #-:sb-simd (fast-zero m)
  m)
(defun fast-zero-matrix (m)
  #+:sb-simd
  (let ((m-s (magicl::matrix/double-float-storage m)))
    (declare ((simple-array double-float (9)) m-s))
    (setf (sb-simd-avx2:f64.4-aref m-s 0) 0d0)
    (setf (sb-simd-avx2:f64.4-aref m-s 4) 0d0)
    (setf (aref m-s 8) 0d0)
    )
  #-:sb-simd (fast-zero m)
  m)
(defun test-fast-zero ()
  (let ((iters 1000000000)
        (vec (cl-mpm/utils:matrix-zeros)))
    (time
     (dotimes (i iters)
       (fast-zero vec)))
    (time
     (dotimes (i iters)
       (fast-zero-matrix vec))))
  )


(declaim (inline voigt-tensor-reduce-lisp)
         (ftype (function (magicl:matrix/double-float) (values double-float)) voigt-tensor-reduce-lisp))
(let ((second-invar (magicl:from-array (make-array 6 :initial-contents '(1d0 1d0 1d0 0.5d0 0.5d0 0.5d0)) '(6 1) :type 'double-float :layout :column-major)))
  (defun voigt-tensor-reduce-lisp (a)
     (values (magicl::sum (fast-.* a a second-invar)))))

(declaim (inline voigt-tensor-reduce-simd)
         (ftype (function (magicl:matrix/double-float) double-float) voigt-tensor-reduce-simd))
(defun voigt-tensor-reduce-simd (a)
  "Calculate the product A_{ij}A_{ij}"
  (let ((arr (magicl::matrix/double-float-storage a)))
    (declare ((simple-array double-float) arr))
    (+
     (* (aref arr 0) (aref arr 0))
     (* (aref arr 1) (aref arr 1))
     (* (aref arr 2) (aref arr 2))
     (* (aref arr 3) (aref arr 3) 0.5d0)
     (* (aref arr 4) (aref arr 4) 0.5d0)
     (* (aref arr 5) (aref arr 5) 0.5d0)
     )))

(defun voigt-tensor-reduce (a)
  (voigt-tensor-reduce-simd a)
  ;; #+:sb-simd (voigt-tensor-reduce-simd a)
  ;; #-:sb-simd (voigt-tensor-reduce-lisp a)
  )

(declaim
 (ftype
  (function
   (magicl::matrix/double-float
    magicl::matrix/double-float
    )
   double-float)
  dot))
(defun dot (a b)
  (the double-float (fast-sum (fast-.* a b)))
  )
(declaim
 (ftype
  (function
   (magicl::matrix/double-float
    magicl::matrix/double-float
    )
   double-float)
  dot))
(defun dot-vector (a b)
  (let ((a-s (cl-mpm/utils:fast-storage a))
        (b-s (cl-mpm/utils:fast-storage b))
        )
    (declare ((simple-array double-float (3)) a-s b-s))
    (+
     (the double-float (expt (- (aref a-s 0) (aref b-s 0)) 2))
     (the double-float (expt (- (aref a-s 1) (aref b-s 1)) 2))
     (the double-float (expt (- (aref a-s 2) (aref b-s 2)) 2)))))

(declaim
 (ftype
  (function
   (magicl::matrix/double-float)
   double-float)
  mag-squared))
(defun mag-squared (a)
  "Calculate the magnitude - 2-norm"
  (dot a a))

(declaim
 (ftype
  (function
   (magicl::matrix/double-float)
   double-float)
  mag))
(defun mag (a)
  "Calculate the magnitude - 2-norm"
  (sqrt (mag-squared a)))

(declaim
 (ftype
  (function
   (magicl::matrix/double-float)
   magicl::matrix/double-float)
  norm))
(defun norm (a)
  "Normalise a vector"
  (let ((m (mag-squared a)))
    (declare (double-float m))
    (if (> m 0d0)
      (cl-mpm/fastmaths:fast-scale a (/ 1d0 (sqrt m)))
      a)))

(defun cross-product (a b)
  "Calculate the cross product of column vectors a and b, returning a new vector"
  (let ((as (magicl::matrix/double-float-storage a))
        (bs (magicl::matrix/double-float-storage b))
        )
    (cl-mpm/utils:vector-from-list (list
                                    (- (* (aref as 1) (aref bs 2)) (* (aref as 2) (aref bs 1)))
                                    (- (* (aref as 2) (aref bs 0)) (* (aref as 0) (aref bs 2)))
                                    (- (* (aref as 0) (aref bs 1)) (* (aref as 1) (aref bs 0)))))))


(declaim
 (inline fast-sum)
 (ftype (function (magicl::matrix/double-float)
                  double-float) fast-sum))
(defun fast-sum (matrix)
  (let ((s (cl-mpm/utils:fast-storage matrix)))
    (declare ((simple-array double-float (*)) s))
    (loop for v double-float across s
          sum v)))
(declaim
 (inline fast-sum-vector)
 (ftype (function (magicl::matrix/double-float)
                  double-float) fast-sum-vector))
(defun fast-sum-vector (matrix)
  (let ((s (cl-mpm/utils:fast-storage matrix)))
    (declare ((simple-array double-float (3)) s))
    (loop for v double-float across s
          sum v)))
(declaim
 (inline fast-sum-voigt)
 (ftype (function (magicl::matrix/double-float)
                  double-float) fast-sum-voigt))
(defun fast-sum-voigt (matrix)
  (declare (optimize (speed 3)))
  (let ((s (cl-mpm/utils:fast-storage matrix)))
    (declare ((simple-array double-float (6)) s))
    (loop for v double-float across s
          sum v)))

(declaim (inline det)
         (ftype (function (magicl::matrix/double-float) (double-float)) det))
(defun det (m)
  (magicl:det m))


(declaim (inline det-3x3)
         (ftype (function (magicl::matrix/double-float) (double-float)) det-3x3))
(defun det-3x3 (m-mat)
  ;; (magicl:det m-mat)
  (let ((a (cl-mpm/utils:fast-storage m-mat)))
    (declare ((simple-array double-float (9)) a))
    (macrolet ((tref (m x y)
                 `(the double-float (aref ,m ,(the fixnum (+ (the fixnum (* 3 (the fixnum y))) (the fixnum x)))))))
      (+
       (* (tref a 0 0) (- (* (tref a 1 1) (tref a 2 2)) (* (tref a 2 1) (tref a 1 2))))
       (* -1d0 (tref a 1 0) (- (* (tref a 0 1) (tref a 2 2)) (* (tref a 2 1) (tref a 0 2))))
       (* (tref a 2 0) (- (* (tref a 0 1) (tref a 1 2)) (* (tref a 1 1) (tref a 0 2)))))))
  )


(defun test-sum ()
  (let ((iters 100000)
        (data (cl-mpm/utils:voigt-from-list (list 1d0 2d0 3d0 4d0 5d0 6d0))))
    (format t "Magicl: ~F~%" (magicl::sum data))
    (format t "fastmaths: ~F~%" (fast-sum data))
    (time
     (dotimes (i iters)
       (magicl::sum data)))
    (time
     (dotimes (i iters)
       (fast-sum data)))
    (time
     (dotimes (i iters)
       (fast-sum-voigt data)))))


(defun vector-principal-j2 (s)
  "Calculate j2 invarient from deviatoric stress"
  (let ((storage (magicl::matrix/double-float-storage s)))
    (* 0.5d0
       (+
        (expt (- (aref storage 0) (aref storage 1)) 2)
        (expt (- (aref storage 1) (aref storage 2)) 2)
        (expt (- (aref storage 2) (aref storage 0)) 2)))))

(defun voigt-j2 (s)
  "Calculate j2 invarient from deviatoric stress"
  (let ((storage (magicl::matrix/double-float-storage s)))
    (/ (+ (the double-float (dot s s))
          (the double-float (expt (aref storage 3) 2))
          (the double-float (expt (aref storage 4) 2))
          (the double-float (expt (aref storage 5) 2))
          ) 2d0)))

(defun voigt-von-mises (stress)
  (sqrt (* 3d0 (voigt-j2 (cl-mpm/utils::deviatoric-voigt stress)))))


(defun vector-principal-von-mises (stress)
  (the double-float
       (sqrt
        (the double-float
             (* 1d0
                (vector-principal-j2
                 (cl-mpm/utils::deviatoric-vector stress)
                 ))))))

(defun eig-values-voigt (voigt)
  (let ((p1 (+ (expt (varef voigt 3) 2)
               (expt (varef voigt 4) 2)
               (expt (varef voigt 5) 2))))
    (if (= p1 0d0)
        (progn
          ;;Diag
          (values (varef voigt 0)
                  (varef voigt 1)
                  (varef voigt 2)))
        (progn
          (let* ((q (/ (cl-mpm/utils:trace-voigt voigt) 3))
                 (p2 (+
                      (expt (- (varef voigt 0) q) 2)
                      (expt (- (varef voigt 1) q) 2)
                      (expt (- (varef voigt 2) q) 2)
                      (* 2 p1)))
                 (p (sqrt (/ p2 6)))
                 (B (cl-mpm/fastmaths:fast-scale!
                     (cl-mpm/fastmaths:fast-.- voigt (cl-mpm/utils::voigt-eye q))
                                                  (/ 1d0 p)))
                 (r (/ (magicl:det (cl-mpm/utils:voigt-to-matrix B)) 2d0)))
            (let ((phi
                    (cond
                      ((<= r -1d0)
                       (/ pi 3))
                      ((>= r 1d0)
                       0d0)
                      (t
                       (/ (acos r) 3)))))
              (let* ((s1 (+ q (* 2 p (cos phi))))
                     (s3 (+ q (* 2 p (cos (+ phi (* pi 2/3))))))
                     (s2 (- (* 3 q) s1 s3)))
                (values
                 s1 s2 s3))))))))

(defun test-eig ()
  (let* ((s (cl-mpm/utils:voigt-from-list (loop repeat 6 collect (- (random 100d0) 50d0))))
        ;(s (cl-mpm/utils:voigt-from-list (list 1d0 1d0 0d0 1d0 1d0 0d0)))
         (iters 100000)
        )
    (time (dotimes (i iters)
            (cl-mpm/utils::principal-stresses-3d s)))
    (time (dotimes (i iters)
            (eig-values-voigt s)))
    (format t "~A~%" (multiple-value-list (cl-mpm/utils::principal-stresses-3d s)))
    (format t "~A~%" (multiple-value-list (eig-values-voigt s)))
    ;; (pprint (eig-values-voigt s))
    ))

(defun @-matrix-matrix-lisp (matrix-a matrix-b result-matrix)
  "Multiply a 3x3 matrix with a 3x3 matrix to calculate a 3x3 vector in place"
  (declare (magicl:matrix/double-float matrix-a matrix-b result-matrix)
           ;; (optimize (speed 0) (safety 3) (debug 3))
           )
  (let ((a (cl-mpm/utils:fast-storage matrix-a))
        (b (cl-mpm/utils:fast-storage matrix-b))
        (c (cl-mpm/utils:fast-storage result-matrix)))
    (declare ((simple-array double-float (9)) a)
             ((simple-array double-float (9)) b)
             ((simple-array double-float (9)) c))
    (macrolet ((tref (m x y)
                 `(aref ,m (+ (* 3 ,y) ,x))))
      (loop for x fixnum from 0 below 3
            do (loop for y fixnum from 0 below 3
                     do (loop for i fixnum from 0 below 3
                              do (incf
                                  (tref c x y)
                                  (the double-float
                                       (*
                                        (tref a x i)
                                        (tref b i y)))
                                  ))))))
  (values))

(defun fast-@-matrix-matrix (mat vec &optional res)
  (let ((res (if res
                 (fast-zero-matrix res)
                 (cl-mpm/utils::matrix-zeros))))
    (@-matrix-matrix-lisp mat vec res)
    res))


(defun @-arb-vector-lisp (matrix vector result-vector)
  "Multiply a 3x9 matrix with a 3x1 vector to calculate a 3x1 vector in place"
  (declare (magicl:matrix/double-float matrix vector result-vector)
           (optimize (speed 3) (safety 0) (debug 0)))
  (let ((a (magicl::matrix/double-float-storage matrix))
        (b (magicl::matrix/double-float-storage vector))
        (c (magicl::matrix/double-float-storage result-vector)))
    (declare ((simple-array double-float (*)) a)
             ((simple-array double-float (*)) b c ))
    (let ((rows (magicl::nrows a))
          (cols (magicl::ncols a)))
      (loop for i fixnum from 0 below rows
            do
               (loop for j fixnum from 0 below cols
                     do (incf (aref c i) (the double-float (* (aref b j) (cl-mpm/utils:mtref matrix i j)))))
            )))
  (values))

(defun @-matrix-vector-lisp (matrix vector scale result-vector)
  "Multiply a 3x9 matrix with a 3x1 vector to calculate a 3x1 vector in place"
  (declare (magicl:matrix/double-float matrix vector result-vector)
           (double-float scale)
           (optimize (speed 3) (safety 0) (debug 0)))
  (let ((a (magicl::matrix/double-float-storage matrix))
        (b (magicl::matrix/double-float-storage vector))
        (c (magicl::matrix/double-float-storage result-vector))
        )
    (declare ((simple-array double-float (9)) a)
             ((simple-array double-float (3)) c)
             ((simple-array double-float (3)) b)
             )
    (flet ((tref (m x y)
             (aref m (+ (* 3 y) x))))
      (loop for i fixnum from 0 below 3
            do
               (loop for j fixnum from 0 below 3
                     do (incf (aref c i) (the double-float (* (aref b j) (tref a i j) scale))))
            )))
  (values))

#+:sb-simd
(defun @-matrix-vector-simd (matrix vector scale result-vector)
  "Multiply a 3x9 matrix with a 3x1 vector to calculate a 3x1 vector in place"
  (declare (magicl:matrix/double-float matrix vector result-vector)
           (double-float scale)
           (optimize (speed 3) (safety 0) (debug 0)))
  (let ((a (magicl::matrix/double-float-storage matrix))
        (b (magicl::matrix/double-float-storage vector))
        (c (magicl::matrix/double-float-storage result-vector))
        )
    (declare ((simple-array double-float (9)) a)
             ((simple-array double-float (3)) c)
             ((simple-array double-float (3)) b)
             )
    (flet ((tref (m x y)
             (aref m (+ (* 3 y) x))))
      (loop for j fixnum from 0 below 3
            do
               (progn
                 (setf
                  (sb-simd-avx:f64.2-aref c 0)
                  (sb-simd-avx:f64.2+
                   (sb-simd-avx:f64.2-aref c 0)
                   (sb-simd-avx:f64.2*
                    (sb-simd-avx:f64.2-aref a (* j 3))
                    (aref b j))))
                 (incf
                  (aref c 2)
                  (the double-float (* (aref b j) (tref a 2 j) scale))
                  )))))
  (values))

(defun fast-@-matrix-vector (mat vec &optional res)
  (let ((res (if res
                 (fast-zero-vector res)
                 (cl-mpm/utils::vector-zeros))
             ))
    #-:sb-simd (@-matrix-vector-lisp mat vec 1d0 res)
    #+:sb-simd (@-matrix-vector-simd mat vec 1d0 res)
    res))
(defun fast-@-matrix-scaled-vector (mat vec scale &optional res)
  (let ((res (if res
                 (fast-zero-vector res)
                 (cl-mpm/utils::vector-zeros))))
    #-:sb-simd (@-matrix-vector-lisp mat vec scale res)
    #+:sb-simd (@-matrix-vector-simd mat vec scale res)
    res))
(defun test-fast-@-matrix-vector ()
  (let ((B (cl-mpm/utils:vector-from-list (list 1d0 2d0 3d0)))
        (A (cl-mpm/utils::matrix-from-list (list 1d0 2d0 3d0
                                                 4d0 5d0 6d0
                                                 7d0 8d0 9d0)))
        )
    (pprint (magicl:@ A B))
    (pprint (fast-@-matrix-vector A B))))

(defun @-tensor-voigt-lisp (matrix-a matrix-b result-matrix)
  "Multiply a 3x3 matrix with a 3x3 matrix to calculate a 3x3 vector in place"
  (declare (magicl:matrix/double-float matrix-a matrix-b result-matrix)
           ;; (optimize (speed 0) (safety 3) (debug 3))
           )
  (let ((a (cl-mpm/utils:fast-storage matrix-a))
        (b (cl-mpm/utils:fast-storage matrix-b))
        (c (cl-mpm/utils:fast-storage result-matrix)))
    (declare ((simple-array double-float (36)) a)
             ((simple-array double-float (6)) b)
             ((simple-array double-float (6)) c))
    (macrolet ((tref (m x y)
                 `(aref ,m (+ (* 6 ,y) ,x))))
      (loop for x fixnum from 0 below 6
            do (loop for i fixnum from 0 below 6
                               do (incf
                                   (aref c x)
                                   (the double-float
                                        (*
                                         (tref a x i)
                                         (aref b i)))
                                   )))))
  (values))
(declaim (ftype (function (magicl:matrix/double-float
                           magicl:matrix/double-float
                           &optional magicl:matrix/double-float)
                          magicl:matrix/double-float
                          )
                fast-@-tensor-voigt))
(defun fast-@-tensor-voigt (mat vec &optional res)
  (let ((res (if res
                 (fast-zero-voigt res)
                 (cl-mpm/utils::voigt-zeros))))
    (@-tensor-voigt-lisp mat vec res)
    res))

(defun fast-@-arb-vector (mat vec &optional res)
  (let ((res (if res
                 (fast-zero res)
                 (cl-mpm/utils::deep-copy vec))))
    (@-arb-vector-lisp mat vec res)
    res))

(defun @-arb-arb-lisp (a b result)
  "Multiply a 3x9 matrix with a 3x1 vector to calculate a 3x1 vector in place"
  (declare (magicl:matrix/double-float a b result)
           (optimize (speed 3) (safety 0) (debug 0)))
  (let ((rows (magicl::ncols b))
        (cols (magicl::nrows a))
        (inter (magicl::nrows a))
        )
    (loop for i fixnum from 0 below rows
          do (loop for j fixnum from 0 below cols
                   do (loop for k fixnum from 0 below inter
                            do (incf
                                (cl-mpm/utils::mtref result i j)
                                (the double-float (* (cl-mpm/utils:mtref a i k)
                                                     (cl-mpm/utils:mtref b k j))))))
          ))
  (values))

(defun fast-@-arb-arb (a b &optional res)
  (let ((res (if res
                 (fast-zero res)
                 (cl-mpm/utils::arb-matrix
                  (magicl::nrows a)
                  (magicl::ncols b)))))
    (@-arb-arb-lisp a b res)
    res))

(defun %lambert-w-0 (w z)
  (- w (/ ( - (* w (exp w)) z)
          (+ (exp w) (* w (exp w))))))

(defun %lambert-w-log-0 (w z)
  (log (/ (- z) (- w))))

(defun lambert-w-0-exp (z)
  (let ((w 1d0))
    (loop for i from 0 to 1000
          while (> w 0d0)
          do (setf w (%lambert-w-0 w z)))
    w))

(defun lambert-w-0-log (z)
  (let ((w 1d0)
        (tol 1d-1))
    (loop for i from 0 to 10000
          while (and
                 (> w 0d0)
                 (> (abs (- z (* w (exp w)))) tol))
          do (setf w (%lambert-w-log-0 w z)))
    (when (> (abs (- z (* w (exp w)))) tol)
      (error "Lambert W failed to solve"))
    w))

(defun lambert-w-0 (z)
  (lambert-w-0-halley z))


(defun lambert-w-0-halley (z)
  (let ((w 1d0)
        (tol 1d-1))
    (loop for i from 0 to 1000
          while (and
                 (> w 0d0)
                 (> (abs (- z (* w (exp w)))) tol))
          do (let* ((W1 (* w (exp w)))
                   (W2 (* (1+ w) (exp w)))
                   )
               ;(setf w (%lambert-w-log-0 w z))
               (decf w (/ (- W1 z) (- W2 (/ (* (+ w 2) (- W1 z)) (+ 2 (* 2 w))))))
               ))
    (when (> (abs (- z (* w (exp w)))) tol)
      (error "Lambert W failed to solve"))
    w))

(defun voigt-det (voigt)
  (let ((a (varef voigt 0))
        (d (varef voigt 1))
        (f (varef voigt 2))
        (b (* 0.5d0 (varef voigt 5)))
        (c (* 0.5d0 (varef voigt 4)))
        (e (* 0.5d0 (varef voigt 3)))
        )
    (+
     (* a (- (* d f) (* e e)))
     (* b (- (* c e) (* b f)))
     (* c (- (* b e) (* d c))))))

(defun voight-det (voigt)
  (let ((a (varef voigt 0))
        (b (varef voigt 5))
        (c (varef voigt 4))
        (d (varef voigt 1))
        (e (varef voigt 3))
        (f (varef voigt 2))
        )
    (+
     (* a (- (* d f) (* e e)))
     (* b (- (* c e) (* b f)))
     (* c (- (* b e) (* d c))))))

(declaim (inline fast-inv-3x3))
(defun fast-inv-3x3 (mat &optional res)
  (let* ((mat-s (cl-mpm/utils:fast-storage mat))
         (res (if res res (cl-mpm/utils:matrix-zeros)))
         (res-s (cl-mpm/utils:fast-storage res)))
    (declare ((simple-array double-float (9)) mat-s res-s))
    (macrolet ((mref (x y) `(the double-float (aref mat-s (+ (* ,y 3) ,x)))))
      (let* ((det
               (+
                (* (mref 0 0) (- (* (mref 1 1) (mref 2 2)) (* (mref 2 1) (mref 1 2))))
                (* -1d0 (mref 1 0) (- (* (mref 0 1) (mref 2 2)) (* (mref 2 1) (mref 0 2))))
                (* (mref 2 0) (- (* (mref 0 1) (mref 1 2)) (* (mref 1 1) (mref 0 2)))))))
        (declare (double-float det))
        (when (= det 0d0)
          (error "Zero determinate"))
        (let ((inv-det (/ 1d0 det)))
          (declare (double-float inv-det))
          (setf
           (aref res-s 0) (* inv-det      (- (* (mref 1 1) (mref 2 2)) (* (mref 2 1) (mref 1 2))))
           (aref res-s 1) (* -1d0 inv-det (- (* (mref 1 0) (mref 2 2)) (* (mref 1 2) (mref 2 0))))
           (aref res-s 2) (* inv-det      (- (* (mref 1 0) (mref 2 1)) (* (mref 2 0) (mref 1 1))))
           (aref res-s 3) (* -1d0 inv-det (- (* (mref 0 1) (mref 2 2)) (* (mref 0 2) (mref 2 1))))
           (aref res-s 4) (* inv-det      (- (* (mref 0 0) (mref 2 2)) (* (mref 0 2) (mref 2 0))))
           (aref res-s 5) (* -1d0 inv-det (- (* (mref 0 0) (mref 2 1)) (* (mref 2 0) (mref 0 1))))
           (aref res-s 6) (* inv-det      (- (* (mref 0 1) (mref 1 2)) (* (mref 0 2) (mref 1 1))))
           (aref res-s 7) (* -1d0 inv-det (- (* (mref 0 0) (mref 1 2)) (* (mref 1 0) (mref 0 2))))
           (aref res-s 8) (* inv-det      (- (* (mref 0 0) (mref 1 1)) (* (mref 1 0) (mref 0 1)))))
          res)))))

(defun linear-solve-3x3-voigt (mat voigt &optional result)
  (let ((result (if result result (cl-mpm/utils:voigt-zeros))))
    (let ((mat-s (cl-mpm/utils:fast-storage mat))
          (voigt-s (cl-mpm/utils:fast-storage voigt))
          (res-s (cl-mpm/utils:fast-storage result)))
      (macrolet ((mref (x y) `(aref mat-s (+ (* ,y 3) ,x))))
        (let* ((det
                 (+
                  (* (mref 0 0) (- (* (mref 1 1) (mref 2 2)) (* (mref 2 1) (mref 1 2))))
                  (* -1d0 (mref 1 0) (- (* (mref 0 1) (mref 2 2)) (* (mref 2 1) (mref 0 2))))
                  (* (mref 2 0) (- (* (mref 0 1) (mref 1 2)) (* (mref 1 1) (mref 0 2)))))))
          (declare (double-float det))
          (when (= det 0d0)
            (error "Zero determinate"))
          (let ((inv-det (/ 1d0 det)))
            (declare (double-float inv-det))
            (let ((a00 (* inv-det      (- (* (mref 1 1) (mref 2 2)) (* (mref 2 1) (mref 1 2)))))
                  (a01 (* -1d0 inv-det (- (* (mref 1 0) (mref 2 2)) (* (mref 1 2) (mref 2 0)))))
                  (a02 (* inv-det      (- (* (mref 1 0) (mref 2 1)) (* (mref 2 0) (mref 1 1)))))
                  (a10 (* -1d0 inv-det (- (* (mref 0 1) (mref 2 2)) (* (mref 0 2) (mref 2 1)))))
                  (a11 (* inv-det      (- (* (mref 0 0) (mref 2 2)) (* (mref 0 2) (mref 2 0)))))
                  (a12 (* -1d0 inv-det (- (* (mref 0 0) (mref 2 1)) (* (mref 2 0) (mref 0 1)))))
                  (a20 (* inv-det      (- (* (mref 0 1) (mref 1 2)) (* (mref 0 2) (mref 1 1)))))
                  (a21 (* -1d0 inv-det (- (* (mref 0 0) (mref 1 2)) (* (mref 1 0) (mref 0 2)))))
                  (a22 (* inv-det      (- (* (mref 0 0) (mref 1 1)) (* (mref 1 0) (mref 0 1))))))
              (declare (double-float a00 a10 a20 a01 a11 a21 a02 a12 a22))
              ;; (pprint
              ;;  (cl-mpm/utils:matrix-from-list (list a00 a10 a20 a01 a11 a21 a02 a12 a22)))
              (setf
               (aref res-s 0)
               (+
                (* (aref voigt-s 0) a00)
                (* (aref voigt-s 5) a10)
                (* (aref voigt-s 4) a20))
               (aref res-s 1)
               (+
                (* (aref voigt-s 5) a01)
                (* (aref voigt-s 1) a11)
                (* (aref voigt-s 3) a21))
               (aref res-s 2)
               (+
                (* (aref voigt-s 4) a02)
                (* (aref voigt-s 3) a12)
                (* (aref voigt-s 2) a22))

               (aref res-s 3)
               (+
                (* (aref voigt-s 4) a01)
                (* (aref voigt-s 3) a11)
                (* (aref voigt-s 2) a21))

               ;; (aref res-s 4)
               ;; (+
               ;;  (* (aref voigt-s 4) a00)
               ;;  (* (aref voigt-s 3) a10)
               ;;  (* (aref voigt-s 2) a20))

               ;; (aref res-s 5)
               ;; (+
               ;;  (* (aref voigt-s 5) a00)
               ;;  (* (aref voigt-s 1) a10)
               ;;  (* (aref voigt-s 3) a20))
               )
              )))))
    result))



(declaim
 (inline lisp-any/)
 (ftype (function ((simple-array double-float)
                   (simple-array double-float)
                   (simple-array double-float))
                  (values)) lisp-any/))
(defun lisp-any/ (a b target)
  (declare ((simple-array double-float (*)) a b target))
  (loop for i fixnum from 0 below (length a)
        do (setf (aref target i) (/ (aref a i) (aref b i))))
  target)

(defun fast-./ (a b &optional res)
  (let ((res (if res
                 res
                 (cl-mpm/utils::empty-copy a))))
    (declare (magicl:matrix/double-float a b res))
    (lisp-any/ (magicl::matrix/double-float-storage a)
               (magicl::matrix/double-float-storage b)
               (magicl::matrix/double-float-storage res))
    res))


(defun matrix-reset-identity (matrix)
  (fast-zero matrix)
  (setf (varef matrix 0) 1d0
        (varef matrix 4) 1d0
        (varef matrix 8) 1d0)
  matrix)


(defun element-wise-map (vec func)
  (declare (function func))
  (let ((s (cl-mpm/utils:fast-storage vec)))
    (loop for i from 0 below (length s)
          do (setf (aref s i) (funcall func (aref s i)))))
  vec)


(defun det-2x2 (mat)
  (- (* (cl-mpm/utils::mtref-2x2 mat 0 0)
        (cl-mpm/utils::mtref-2x2 mat 1 1))
     (* (cl-mpm/utils::mtref-2x2 mat 0 1)
        (cl-mpm/utils::mtref-2x2 mat 1 0))))

(defun eigenvalues-2x2 (mat)
  (let* ((tr (+ (cl-mpm/utils:varef mat 0) (cl-mpm/utils:varef mat 3)))
         (det (det-2x2 mat))
         (gap (* 0.5d0 (sqrt (- (* tr tr) (* 4 det)))))
         (mid (/ tr 2))
         (l1 (+ mid gap))
         (l2 (- mid gap)))
    (values l1 l2)))


(defstruct eigen-work-object
  work a-tensor w)


(defparameter *test-pool* nil)
(defparameter *work-pool* (cl-mpm/utils::make-object-pool
                           :constructor (lambda ()
                                          (make-eigen-work-object
                                           :work (make-array 102 :element-type 'double-float)
                                           :a-tensor (cl-mpm/utils:matrix-zeros)
                                           :w (make-array 3 :element-type 'double-float)))))
(let ((work-pool *work-pool*))
  ;; (setf (cl-mpm/utils::object-pool-lock work-pool) (sb-thread:make-mutex)
  ;;       (cl-mpm/utils::object-pool-pool work-pool) (make-array 0)
  ;;       )
  ;; (pprint (cl-mpm/utils::object-pool-grab work-pool))
  ;; (setf *test-pool* work-pool)
  ;; (let ((work (make-array 102 :element-type 'double-float))
  ;;       (a-tensor (cl-mpm/utils:matrix-zeros))
  ;;       (w (make-array 3 :element-type 'double-float))))
  (defun magicl-eigen-values-3x3 (m)
    (policy-cond:with-expectations (> speed safety)
        ((assertion (magicl:square-matrix-p m)))
      (let ((rows 3)
            ;; (a-tensor (cl-mpm/utils:matrix-copy m))
            )
        (let* ((wo (cl-mpm/utils::object-pool-grab work-pool))
               (a-tensor (eigen-work-object-a-tensor wo))
               (w (eigen-work-object-w wo))
               (work (eigen-work-object-work wo)))
          (cl-mpm/utils:matrix-copy-into m a-tensor)
          ;; (when (eql :row-major (magicl::matrix/double-float-layout m)) (magicl:transpose! a-tensor))
          (let ((jobz "N")
                (uplo "U")
                (n rows)
                (a (cl-mpm/utils:fast-storage a-tensor))
                (lda rows)
                ;; (w (make-array rows :element-type 'double-float))
                ;; (work (make-array 1 :element-type 'double-float))
                (lwork -1)
                (info 0))
            ;; run it once as a workspace query
            ;; (magicl.lapack-cffi::%dsyev jobz uplo n a lda w work lwork info)
            ;; (setf lwork (truncate (realpart (row-major-aref work 0))))
            ;; (pprint lwork)
            (setf lwork 102)
            ;; (setf work (make-array (max 1 lwork) :element-type 'double-float))
            ;; run it again with optimal workspace size
            (magicl.lapack-cffi::%dsyev jobz uplo n a lda w work lwork info)
            ;; (coerce w 'list)
            (list (aref w 0) (aref w 1) (aref w 2))
            ;; (values (coerce w 'list) a-tensor)
            ))))))

;;Algorithm for closed-form eigenvalues from "A Robust Eigensolver for 3 × 3 Symmetric Matrices" - David Eberly
(defun eigenvalues-3x3 (voigt)
  (let* ((s (cl-mpm/utils::fast-storage voigt))
         (a00 (aref s 0))
         (a11 (aref s 1))
         (a22 (aref s 2))
         (a01 (* 0.5d0 (aref s 5)))
         (a02 (* 0.5d0 (aref s 4)))
         (a12 (* 0.5d0 (aref s 3)))
         )
    (declare (double-float a00 a01 a02 a11 a12 a22))
    (let ((maxAbsElement (max (abs a00) (abs a01) (abs a02) (abs a12) (abs a22) (abs a11))))
      (if (> maxAbsElement 0d0)
          (let* ((invMaxAbsElement (/ 1d0 maxAbsElement))
                 (a00 (* a00 invMaxAbsElement))
                 (a01 (* a01 invMaxAbsElement))
                 (a02 (* a02 invMaxAbsElement))
                 (a11 (* a11 invMaxAbsElement))
                 (a12 (* a12 invMaxAbsElement))
                 (a22 (* a22 invMaxAbsElement))
                 (norm (+ (* a01 a01) (* a02 a02) (* a12 a12))))
            (declare (double-float a00 a01 a02 a11 a12 a22 norm))
            (if (> norm 0d0)
                (let* ((q (/ (+ a00 a11 a22) 3d0))
                       (b00 (- a00 q))
                       (b11 (- a11 q))
                       (b22 (- a22 q))
                       (p (the double-float
                               (sqrt
                                (/
                                 (the double-float (+ (* b00 b00) (* b11 b11) (* b22 b22) (* norm 2d0))) 6d0))))
                       (c00 (- (* b11 b22) (* a12 a12)))
                       (c01 (- (* a01 b22) (* a12 a02)))
                       (c02 (- (* a01 a12) (* b11 a02)))
                       (det (/ (+ (* b00 c00) (* -1d0 a01 c01) (* a02 c02)) (* p p p)))
                       (halfDet (min 1d0 (max -1d0 (/ det 2d0))))
                       ;; (halfDet (/ det 2d0))
                       (angle (/ (the double-float (acos halfDet)) 3d0))
                       ;; The number of digits in twoThirdsPi is chosen so that, whether float or double,
                       ;; the floating-point number is the closest to theoretical 2*pi/3.
                       (twoThirdsPi 2.09439510239319549d0)
                       (beta2 (* 2d0 (cos angle)))
                       (beta0 (* 2d0 (cos (+ angle twoThirdsPi))))
                       (beta1 (- (+ beta0 beta2))))
                  ;; The eigenvalues of A are ordered as alpha0 <= alpha1 <= alpha2.
                  ;; eval[0] = q + p * beta0;
                  ;; eval[1] = q + p * beta1;
                  ;; eval[2] = q + p * beta2;
                  ;; The preconditioning scaled the matrix A, which scales the
                  ;; eigenvalues.
                  (values (the double-float (* (+ q (* p beta0)) maxAbsElement))
                          (the double-float (* (+ q (* p beta1)) maxAbsElement))
                          (the double-float (* (+ q (* p beta2)) maxAbsElement))))
                ;;Diagonal matrix
                (values (the double-float (* maxAbsElement a00))
                        (the double-float (* maxAbsElement a11))
                        (the double-float (* maxAbsElement a22))))
            )
          (values 0d0 0d0 0d0)))))
;; (let ((m (cl-mpm/utils:voigt-from-list (list 1d0 2d0 3d0 4d0 2d0 1d0)))
;;       (iters 1000))
;;   (pprint (multiple-value-list (eigenvalues-3x3 m)))
;;   (pprint (magicl-eigen-values-3x3 (cl-mpm/utils::voigt-to-matrix m)))
;;   ;; (time (dotimes (i iters) (eigenvalues-3x3 m)))
;;   ;; (time (dotimes (i iters) (magicl-eigen-values-3x3 (cl-mpm/utils::voigt-to-matrix m))))
;;   )


;; (let ((m (cl-mpm/utils::matrix-from-list (loop repeat 9 collect (random 1d0))))
;;       (iters 16))
;;   (time
;;    (lparallel:pdotimes (i iters)
;;      (cl-mpm/fastmaths::magicl-eigen-values-3x3 m)))
;;   (time
;;    (lparallel:pdotimes (i iters)
;;      (cl-mpm/utils::eig m)))
;;   )

(defun fast-@-sparse-mat-dense-vec (mat vec &optional res)
  (let ((res (if res
                 (cl-mpm/fastmaths::fast-zero res)
                 (cl-mpm/utils::arb-matrix (cl-mpm/utils::sparse-matrix-nrows mat) 1))))
    (declare (magicl:matrix/double-float vec res)
             (cl-mpm/utils::sparse-matrix mat))

    (let ((rowindex (cl-mpm/utils::sparse-matrix-rowindex mat))
          (cols (cl-mpm/utils::sparse-matrix-cols mat))
          (values (cl-mpm/utils::sparse-matrix-values mat))
          (res-s (cl-mpm/utils::fast-storage res))
          (vec-s (cl-mpm/utils::fast-storage vec)))
      (dotimes (r (cl-mpm/utils::sparse-matrix-nrows mat))
        (let ((col-0 (aref rowindex r))
              (col-1 (aref rowindex (1+ r))))
          (loop for c from col-0 below col-1 do
            (incf (aref res-s r)
                  (the double-float
                       (*
                        (aref values c)
                        (aref vec-s (aref cols c)))))))))
    res))

(defun fast-@-sparse-mat-dense-vec-multithread (mat vec &optional res)
  (let ((res (if res
                 (cl-mpm/fastmaths::fast-zero res)
                 (cl-mpm/utils::arb-matrix (cl-mpm/utils::sparse-matrix-nrows mat) 1))))
    (declare (magicl:matrix/double-float vec res)
             (cl-mpm/utils::sparse-matrix mat))

    (let ((rowindex (cl-mpm/utils::sparse-matrix-rowindex mat))
          (cols (cl-mpm/utils::sparse-matrix-cols mat))
          (values (cl-mpm/utils::sparse-matrix-values mat))
          (res-s (cl-mpm/utils::fast-storage res))
          (vec-s (cl-mpm/utils::fast-storage vec)))
      (lparallel:pdotimes (r (cl-mpm/utils::sparse-matrix-nrows mat))
        (let ((col-0 (aref rowindex r))
              (col-1 (aref rowindex (1+ r))))
          (loop for c from col-0 below col-1 do
            (incf (aref res-s r)
                  (the double-float
                       (*
                        (aref values c)
                        (aref vec-s (aref cols c)))))))))
    res))

(defun fast-@-sparse-mat-dense-vec-masked (mat vec bcs-r bcs-c &optional res)
  (let ((res (if res
                 res
                 (cl-mpm/utils::arb-matrix (cl-mpm/utils::sparse-matrix-nrows mat) 1))))
    (declare (magicl:matrix/double-float vec res)
             (cl-mpm/utils::sparse-matrix mat))
    (let ()
      (let ()
        ;;When we are solving a fully fixed system - i.e. out of plane dimensions
        (let* ((rowindex (cl-mpm/utils::sparse-matrix-rowindex mat))
               (cols (cl-mpm/utils::sparse-matrix-cols mat))
               (values (cl-mpm/utils::sparse-matrix-values mat))
               (res-s (cl-mpm/utils::fast-storage res))
               (vec-s (cl-mpm/utils::fast-storage vec)))
          (declare ((simple-array double-float *) values)
                   ((simple-array fixnum *) rowindex cols)
                   ((simple-array double-float *) res-s vec-s))
          (dotimes (r (length res-s))
            (when (> (varef bcs-r r) 0d0)
              (let* ((col-0 (aref rowindex r))
                     (col-1 (aref rowindex (1+ r))))
                (loop for c from col-0 below col-1
                        do
                           (when (> (varef bcs-c (aref cols c)) 0d0)
                             (incf (aref res-s r)
                                   (the double-float
                                        (*
                                         (aref values c)
                                         (aref vec-s (aref cols c)))))))))))))
    (fast-.* res bcs-r res)
    res))
(defun fast-@-sparse-mat-dense-vec-masked-multithread (mat vec bcs-r bcs-c &optional res)
  (let ((res (if res
                 res
                 (cl-mpm/utils::arb-matrix (cl-mpm/utils::sparse-matrix-nrows mat) 1))))
    (declare (magicl:matrix/double-float vec res)
             (cl-mpm/utils::sparse-matrix mat))
    (let ()
      (let ()
        ;;When we are solving a fully fixed system - i.e. out of plane dimensions
        (let* ((rowindex (cl-mpm/utils::sparse-matrix-rowindex mat))
               (cols (cl-mpm/utils::sparse-matrix-cols mat))
               (values (cl-mpm/utils::sparse-matrix-values mat))
               (res-s (cl-mpm/utils::fast-storage res))
               (vec-s (cl-mpm/utils::fast-storage vec)))
          (declare ((simple-array double-float *) values)
                   ((simple-array fixnum *) rowindex cols)
                   ((simple-array double-float *) res-s vec-s))
          (lparallel::pdotimes (r (length res-s))
            (when (> (varef bcs-r r) 0d0)
              (let* ((col-0 (aref rowindex r))
                     (col-1 (aref rowindex (1+ r))))
                (loop for c from col-0 below col-1
                        do
                           (when (> (varef bcs-c (aref cols c)) 0d0)
                             (incf (aref res-s r)
                                   (the double-float
                                        (*
                                         (aref values c)
                                         (aref vec-s (aref cols c)))))))))))))
    (fast-.* res bcs-r res)
    res))
