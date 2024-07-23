(defpackage :cl-mpm/fastmath
  (:use :cl)
  (:import-from
    :magicl tref .+ .-)
  (:export
   #:fast-add
   #:fast-fmacc
   #:mult
   #:det
   #:dot
   #:norm
   :fast-.+
   :fast-.-
   :fast-.*
   :fast-scale-vector
   :fast-scale-voigt
   :fast-scale!
   :fast-zero
   :fast-sum
   :mag
    ))

(declaim (optimize (debug 0) (safety 0) (speed 3)))
;; (declaim (optimize (debug 3) (safety 3) (speed 0)))
(in-package :cl-mpm/fastmath)

(pushnew :sb-simd *features*)
(eval-when
    (:compile-toplevel)
  (pushnew :sb-simd *features*)
  #+:sb-simd (require 'sb-simd)
  )

#+:sb-simd
(progn
  ;(require 'sb-simd)
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
    (values)))

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
      )))

(declaim
 (inline simd-diff-norm)
 (ftype (function ((simple-array double-float (3))
                   (simple-array double-float (3))) (values)) simd-diff-norm))
(defun simd-diff-norm (a b)
  (declare (type (simple-array double-float (3)) a b))
  (+
   (the double-float (expt (- (aref a 0) (aref b 0)) 2))
   (the double-float (expt (- (aref a 1) (aref b 1)) 2))
   (the double-float (expt (- (aref a 2) (aref b 2)) 2))
   )
  ;; (let ((temp
  ;;         (sb-simd-avx:f64.2-
  ;;          (sb-simd-avx:f64.2-aref a 0)
  ;;          (sb-simd-avx:f64.2-aref b 0))))
  ;;   (+
  ;;    (sb-simd-avx:f64.2-horizontal+ (sb-simd-avx:f64.2*
  ;;                                    temp
  ;;                                    temp))
  ;;    (the double-float (expt (- (aref a 2) (aref b 2)) 2))))
  )

(declaim
 (inline fast-fmacc-array)
 (ftype (function ((simple-array double-float)
                   (simple-array double-float)
                   double-float) (values)) fast-fmacc-array))
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
  #-:sb-simd (fast-.+ a (magicl:scale b d) a)
  ;; (fast-.+ a (magicl:scale b d) a)
  )

(declaim
   (inline fast-add)
   (ftype (function (magicl:matrix/double-float magicl:matrix/double-float) (values)) fast-add))
(defun fast-add (a b)
  #+:sb-simd (simd-add a b)
  #-:sb-simd (magicl:.+ a b a)
  )

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
      ;; (loop for offset fixnum from 0 by 2
      ;;       repeat iter
      ;;       do
      ;;          (setf (sb-simd-avx:f64.2-aref target offset)
      ;;                (sb-simd-avx:f64.2+
      ;;                 (sb-simd-avx:f64.2-aref a offset)
      ;;                 (sb-simd-avx:f64.2-aref b offset))))
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
          (incf offset 1))
        )
      ))

  ;; (loop for i from 0 below (length a)
  ;;       do (setf (aref target i) (+ (aref a i) (aref b i))))
  target)
(declaim
 (inline simd-any-)
 (ftype (function ((simple-array double-float)
                           (simple-array double-float)
                           (simple-array double-float))
                          (values)) simd-any-))
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
  target)


(defun simd-any* (a b target)
  ;(declare (optimize (speed 0) (debug 0) (safety 3)))
  (declare ((simple-array double-float (*)) a b target))
  (let ((offset 0))
    (declare (type sb-simd:f64vec a b target)
             (fixnum offset))
    (multiple-value-bind (iter remain) (floor (length a) 2)
      (declare (fixnum iter remain))
      ;; (loop for offset fixnum from 0 by 2
      ;;       repeat iter
      ;;       do
      ;;          (setf (sb-simd-avx:f64.2-aref target offset)
      ;;                (sb-simd-avx:f64.2+
      ;;                 (sb-simd-avx:f64.2-aref a offset)
      ;;                 (sb-simd-avx:f64.2-aref b offset))))
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
  target)

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
          (incf offset 1))
        )
      ))

  ;; (loop for i from 0 below (length a)
  ;;       do (setf (aref target i) (+ (aref a i) (aref b i))))
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
  )

(defun @-m-v (matrix vector result-vector)
  "Multiply a 3x3 matrix with a 3x1 vector to calculate a 3x1 vector in place"
  (declare (magicl:matrix/double-float matrix vector result-vector)
           ;; (optimize (speed 0) (safety 3) (debug 0))
           )
  (let ((a (magicl::matrix/double-float-storage matrix))
        (b (magicl::matrix/double-float-storage vector))
        (c (magicl::matrix/double-float-storage result-vector))
        )
    (declare ((simple-array double-float (9)) a)
             ((simple-array double-float (3)) b c))
    (flet ((tref (m x y)
             (aref m (+ x (* 3 y)))))
      (loop for i fixnum from 0 below 3
            do
               (setf
                (aref c i)
                (+
                 (* (aref b 0) (tref a i 0))
                 (* (aref b 1) (tref a i 1))
                 (* (aref b 2) (tref a i 2))
                 )
                )
            )))
  result-vector)

(declaim
 (inline @-stretch-vec)
 (ftype (function (magicl:matrix/double-float magicl:matrix/double-float magicl:matrix/double-float)
                          magicl:matrix/double-float) @-stretch-vec))
;; (defun @-stretch-vec (matrix vector result-vector)
;;   "Multiply a 3x9 matrix with a 3x1 vector to calculate a 3x1 vector in place"
;;   (declare (magicl:matrix/double-float matrix vector result-vector)
;;            (optimize (speed 0) (safety 3) (debug 0)))
;;   (let ((a (magicl::matrix/double-float-storage matrix))
;;         (b (magicl::matrix/double-float-storage vector))
;;         (c (magicl::matrix/double-float-storage result-vector))
;;         )
;;     (declare ((simple-array double-float (27)) a)
;;              ((simple-array double-float (9)) c)
;;              ((simple-array double-float (3)) b)
;;              )
;;     (flet ((tref (m x y)
;;              (aref m (+ (* 9 x)  y))))
;;       (loop for i fixnum from 0 below 9
;;             do
;;                (setf (aref c i) 0d0)
;;                (loop for j fixnum from 0 below 3
;;                      do (incf (aref c i) (the double-float (* (aref b j) (tref a j i)))))
;;             )))
;;   result-vector)
(defun @-stretch-vec (matrix vector result-vector)
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
;;     (cl-mpm/fastmath::@-stretch-vec stretch vel res)
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
;;   (let ((res (cl-mpm/utils::deep-copy a)))
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
                 (cl-mpm/utils::deep-copy a))))
    (declare (magicl:matrix/double-float a b res))
    (simd-any+ (magicl::matrix/double-float-storage a)
               (magicl::matrix/double-float-storage b)
               (magicl::matrix/double-float-storage res))
    res))

(declaim
 (inline fast-.-)
 (ftype (function (magicl:matrix/double-float magicl:matrix/double-float &optional magicl:matrix/double-float) magicl:matrix/double-float) fast-.-))
(defun fast-.- (a b &optional res)
  (let ((res (if res
                 res
                 (cl-mpm/utils::deep-copy a))))
    (declare (magicl:matrix/double-float a b res))
    (simd-any- (magicl::matrix/double-float-storage a)
               (magicl::matrix/double-float-storage b)
               (magicl::matrix/double-float-storage res))
    res))

(defun fast-.* (a b &optional res)
  (let ((res (if res
                 res
                 (cl-mpm/utils::deep-copy a))))
    (declare (magicl:matrix/double-float a b res))
    (simd-any* (magicl::matrix/double-float-storage a)
               (magicl::matrix/double-float-storage b)
               (magicl::matrix/double-float-storage res))
    res))

;; Here we define some specialised fast adds for different length vectors -> this should almost certianly use the type system to resolve it but who knows
(macrolet ((def-fast-.+-type (name length underlying-add)
             (declare (fixnum length))
             `(progn
                (declaim
                 (inline ,name)
                 (ftype (function (magicl:matrix/double-float magicl:matrix/double-float &optional magicl:matrix/double-float) magicl:matrix/double-float) ,name))
                (defun ,name (a b &optional res)
                  (let ((res (if res
                                 res
                                 (cl-mpm/utils::deep-copy a))))
                    (declare (magicl:matrix/double-float a b res))
                    (,underlying-add (the (simple-array double-float (,length)) (magicl::matrix/double-float-storage a))
                                     (the (simple-array double-float (,length)) (magicl::matrix/double-float-storage b))
                                     (the (simple-array double-float (,length)) (magicl::matrix/double-float-storage res)))
                    res)
                  ))))
  (def-fast-.+-type fast-.+-vector 3 simd-any+)
  (def-fast-.+-type fast-.+-voigt 6 simd-any+)
  (def-fast-.+-type fast-.+-matrix 9 simd-any+-4)
  (def-fast-.+-type fast-.+-stretch 27 simd-any+-4)
  )
(defun test-.+-vector ()
  (let ((a (cl-mpm/utils::vector-from-list (list 1d0 2d0 3d0)))
        (b (cl-mpm/utils::vector-from-list (list 3d0 5d0 9d0))))
    (let ((res (cl-mpm/utils::vector-zeros)))
      (magicl:.+ a b res)
      (pprint res))
    (let ((res (cl-mpm/utils::vector-zeros)))
      (cl-mpm/fastmath::fast-.+ a b res)
      (pprint res))
    (let ((res (cl-mpm/utils::vector-zeros)))
      (cl-mpm/fastmath::fast-.+-vector a b res)
      (pprint res))))
(defun test-.+-voigt ()
  (let ((a (cl-mpm/utils::voigt-from-list (list 1d0 2d0 3d0 2d0 5d0 7d0)))
        (b (cl-mpm/utils::voigt-from-list (list 3d0 5d0 9d0 3d0 4d0 9d0))))
    (let ((res (cl-mpm/utils::voigt-zeros)))
      (magicl:.+ a b res)
      (pprint res))
    (let ((res (cl-mpm/utils::voigt-zeros)))
      (cl-mpm/fastmath::fast-.+ a b res)
      (pprint res))
    (let ((res (cl-mpm/utils::voigt-zeros)))
      (cl-mpm/fastmath::fast-.+-voigt a b res)
      (pprint res))))

(defun test-.+-matrix ()
  (let ((a (cl-mpm/utils::matrix-from-list (loop for i from 0 below 9 collect (random 1d0))))
        (b (cl-mpm/utils::matrix-from-list (loop for i from 0 below 9 collect (random 1d0)))))
    (let ((res (cl-mpm/utils::matrix-zeros)))
      (magicl:.+ a b res)
      (pprint res))
    (let ((res (cl-mpm/utils::matrix-zeros)))
      (cl-mpm/fastmath::fast-.+ a b res)
      (pprint res))
    (let ((res (cl-mpm/utils::matrix-zeros)))
      (cl-mpm/fastmath::fast-.+-matrix a b res)
      (pprint res))))


;; (declaim
;;  (inline fast-.+-voigt)
;;  (ftype (function (magicl:matrix/double-float magicl:matrix/double-float &optional magicl:matrix/double-float) magicl:matrix/double-float) fast-.+-voigt))
;; (defun fast-.+-voigt (a b &optional res)
;;   (declare (optimize (speed 0) (space 0)))
;;   (let ((res (if res
;;                  res
;;                  (cl-mpm/utils::deep-copy a))))
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

(declaim
 (inline fast-scale-voigt)
 (ftype (function (magicl:matrix/double-float double-float) magicl:matrix/double-float) fast-scale-voigt))
(defun fast-scale-voigt (m scale)
  (fast-scale! (cl-mpm/utils:voigt-copy m) scale))

(defun fast-scale-vector (m scale)
  (fast-scale! (cl-mpm/utils:vector-copy m) scale))

(declaim
 (inline fast-zero)
 (ftype (function (magicl:matrix/double-float) magicl:matrix/double-float) fast-zero))
(defun fast-zero (m)
  (let ((m-s (magicl::matrix/double-float-storage m)))
    (loop for i fixnum from 0 below (length m-s)
          do (setf (aref m-s i) 0d0))) m)
(declaim
 (inline fast-zero-vector)
 (ftype (function (magicl:matrix/double-float) magicl:matrix/double-float) fast-zero-vector))
(defun fast-zero-vector (m)
  (let ((m-s (magicl::matrix/double-float-storage m)))
    (declare ((simple-array double-float (3)) m-s))
    ;; (setf (sb-simd-avx:f64.2-aref m-s 0) (sb-simd-avx:make-f64.2 0d0 0d0))
    (setf (sb-simd-avx:f64.2-aref m-s 0) 0d0)
    (setf (aref m-s 2) 0d0)
    ;; (loop for i fixnum from 0 below (length m-s)
    ;;       do (setf (aref m-s i) 0d0))
    ) m)
(defun fast-zero-voigt (m)
  (let ((m-s (magicl::matrix/double-float-storage m)))
    (declare ((simple-array double-float (6)) m-s))
    (setf (sb-simd-avx:f64.4-aref m-s 0) 0d0)
    (setf (sb-simd-avx:f64.2-aref m-s 4) 0d0)
    ) m)
(defun fast-zero-matrix (m)
  (let ((m-s (magicl::matrix/double-float-storage m)))
    (declare ((simple-array double-float (9)) m-s))
    (setf (sb-simd-avx2:f64.4-aref m-s 0) 0d0)
    (setf (sb-simd-avx2:f64.4-aref m-s 4) 0d0)
    (setf (aref m-s 8) 0d0)
    ) m)
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


;; (declaim
;;  (inline fast-.+)
;;  (ftype (function (magicl:matrix/double-float magicl:matrix/double-float) (magicl:matrix/double-float)) fast-.+))
;; (defun fast-.+ (a b)
;;   (let ((res ))
;;     )
;;   )
;; (declaim
;;  (inline fast-.+)
;;  (ftype (function (magicl:matrix/double-float magicl:matrix/double-float magicl:matrix/double-float) (magicl:matrix/double-float)) fast-.+))
;; (defun fast-.+ (a b res)
;;   )


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
  #+:sb-simd (voigt-tensor-reduce-simd a)
  #-:sb-simd (voigt-tensor-reduce-lisp a)
  )

;; (declaim (inline flat-tensor-reduce-simd)
;;          (ftype (function (magicl:matrix/double-float) double-float) flat-tensor-reduce-simd))
;; (defun flat-tensor-reduce-simd (a)
;;   "Calculate the product A_{ij}A_{ij} but without voigt notation"
;;   (let ((arr (magicl::matrix/double-float-storage a)))
;;     (declare ((simple-array double-float) arr))
;;     (+
;;      (* (aref arr 0) (aref arr 0))
;;      (* (aref arr 1) (aref arr 1))
;;      (* (aref arr 2) (aref arr 2) 0.5d0))))

;; (declaim (inline stretch-to-sym)
;;          (ftype (function (magicl:matrix/double-float magicl:matrix/double-float) (values)) stretch-to-sym))
;; (defun stretch-to-sym (stretch result)
;;   ;; (unless result
;;   ;;   (setf result (cl-mpm/utils:voigt-zeros)))
;;   (progn
;;     (declaim (magicl:matrix/double-float result))

;;     (let ((res (matrix-to)(magicl:.+ stretch (magicl:transpose stretch)))))
;;     ;; (loop for i from 0 below 3
;;     ;;       do
;;     ;;          (setf (magicl:tref result  i 0)
;;     ;;                (magicl:tref stretch i 0)))
;;     (setf (magicl:tref result  0 0)
;;           (magicl:tref stretch 0 0))

;;     (setf (magicl:tref result  1 0)
;;           (magicl:tref stretch 1 1))
;;     ;;Since off diagonal components get halved, then voigt doubles them this is net 1d0
;;     (setf (magicl:tref result 5 0)
;;           (* 1d0 (+ (the double-float (magicl:tref stretch 0 1))
;;                     (the double-float (magicl:tref stretch 1 0)))))
;;     )
;;   (values))


;; (defun mult-transpose-accumulate (a b scale res)
;;   "(incf res (scale! (@ (tranpose a) b) scale))"
;;   (declare (type magicl:matrix/double-float a b res)
;;            (type double-float scale))
;;   (declare (optimize (safety 3)))
;;   (let (
;;         ;; (a-s (magicl::matrix/double-float-storage a))
;;         ;; (b-s (magicl::matrix/double-float-storage b))
;;         ;; (res-s (magicl::matrix/double-float-storage res))
;;         )
;;     (loop for i from 0 to 1
;;           do (loop for j from 0 to 2
;;                    do (incf (the double-float (magicl:tref res i 0)) (* (the double-float (magicl:tref a j i))
;;                                                                         (the double-float (magicl:tref b j 0)) scale))))))
(declaim
 (ftype
  (function
   (magicl::matrix/double-float
    magicl::matrix/double-float
    )
   double-float)
  dot))
(defun dot (a b)
  ;; (let ((as (cl-mpm/utils:fast-storage a))
  ;;       (bs (cl-mpm/utils:fast-storage b)))
  ;;   (declare ((simple-array double-float (*)) as bs))
  ;;   (the double-float
  ;;        (loop for va across as
  ;;              for vb across bs
  ;;              sum (the double-float (* va vb)))))
  (the double-float (fast-sum (fast-.* a b)))
  )

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
      (magicl::scale a (/ 1d0 (sqrt m)))
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


(defun test-sum ()
  (let ((iters 100000)
        (data (cl-mpm/utils:voigt-from-list (list 1d0 2d0 3d0 4d0 5d0 6d0))))
    (format t "Magicl: ~F~%" (magicl::sum data))
    (format t "fastmath: ~F~%" (fast-sum data))
    (time
     (dotimes (i iters)
       (magicl::sum data)))
    (time
     (dotimes (i iters)
       (fast-sum data)))
    (time
     (dotimes (i iters)
       (fast-sum-voigt data)))))
