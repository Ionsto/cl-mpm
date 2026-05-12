(defpackage :cl-mpm/utils
  (:use :cl)
  (:export
   #:voigt-zeros
   #:matrix-zeros
   #:tensor-voigt-4th-zeros
   #:stretch-dsvp-zeros
   #:stretch-dsvp-3d-zeros
   #:voigt-from-list
   #:matrix-from-list
   #:tensor-voigt-4th-from-list
   #:matrix-to-voight
   #:voight-to-matrix
   #:voight-to-stretch
   #:voight-to-stretch-3d
   #:matrix-to-voight-strain
   #:voight-to-matrix-strain
   #:voigt-to-voight-strain
   #:voigt-strain-to-voight
   #:matrix-to-mandel
   #:mandel-to-matrix
   #:voigt-to-mandel
   #:mandel-to-voigt
   #:stress-from-list
   #:deviatoric-voigt
   #:trace-voigt
   #:matrix-to-voigt
   #:voigt-to-matrix
   #:vector-zeros
   #:vector-from-list
   #:eig
   #:@-mat-vec
   #:@-tensor-voigt
   #:get-stress
   #:get-vector
   #:matrix-copy
   #:vector-copy
   #:voigt-copy
   #:matrix-copy-into
   #:vector-copy-into
   #:voigt-copy-into
   #:voigt-eye
   #:matrix-eye
   #:fast-storage
   #:varef
   #:mtref
   #:deg-to-rad
   #:rad-to-deg
   #:set-workers
   ))
(in-package :cl-mpm/utils)
(declaim (optimize (debug 0) (safety 0) (speed 3)))
;; (declaim (optimize (debug 3) (safety 3) (speed 0)))

(declaim (inline eig)
         (ftype (function (magicl:matrix/double-float)
                          (values list magicl:matrix/double-float)) eig))
(defun eig (mat)
  "Real eigen-decomposition"
  (magicl:self-adjoint-eig mat)
  ;; (magicl:hermitian-eig mat)
  ;; (multiple-value-bind (l v) (magicl:eig mat)
  ;;   (values l (magicl:.realpart v)))
  )


(defmacro check-nan-matrix (mat &rest body)
  `(loop for vsdadasdsd across (cl-mpm/utils::fast-storage ,mat)
         do (when (or
                   (sb-ext::float-nan-p vsdadasdsd)
                   (= (abs vsdadasdsd) #.sb-ext:double-float-positive-infinity)
                   (> (abs vsdadasdsd) 1d10))
             ,@body)))

(declaim (inline fast-storage)
         (ftype (function (magicl:matrix/double-float)
                          (simple-array double-float (*))) fast-storage))
(defun fast-storage (m)
  (magicl::matrix/double-float-storage m))

(declaim (inline varef)
         (ftype (function (magicl:matrix/double-float fixnum)
                          double-float) varef))
(defun varef (m i)
  (declare (magicl:matrix/double-float m)
           (fixnum i))
  (aref (the (simple-array double-float (*)) (fast-storage m)) i))

(declaim (inline mtref)
         (ftype (function (magicl:matrix/double-float fixnum fixnum)
                          double-float) mtref))
(defun mtref (m row col)
  (declare (magicl:matrix/double-float m)
           (fixnum row col))
  (policy-cond:with-expectations (> speed safety)
      ((assertion (eq (magicl::matrix/double-float-layout m) :column-major))
       (assertion (< row (magicl::matrix/double-float-nrows m)))
       (assertion (< col (magicl::matrix/double-float-ncols m))))
    (let ((numrows (magicl::matrix-nrows m)))
      (declare (type fixnum numrows))
      (aref (fast-storage m) (+ row (the fixnum (* col numrows)))))))

(declaim (inline mtref-3x3)
         (ftype (function (magicl:matrix/double-float fixnum fixnum)
                          double-float) mtref-3x3))
(defun mtref-3x3 (m row col)
  (declare (magicl:matrix/double-float m)
           (fixnum row col))
  (policy-cond:with-expectations (> speed safety)
      ((assertion (eq (magicl::matrix/double-float-layout m) :column-major))
       (assertion (= 3 (magicl::matrix/double-float-nrows m)))
       (assertion (= 3 (magicl::matrix/double-float-ncols m))))
    (aref (fast-storage m) (the fixnum (+ row (the fixnum (* 3 col)))))))
(declaim (inline mtref-2x2)
         (ftype (function (magicl:matrix/double-float fixnum fixnum)
                          (double-float)) mtref-2x2))
(defun mtref-2x2 (m row col)
  (declare (magicl:matrix/double-float m)
           (fixnum row col))
  (policy-cond:with-expectations (> speed safety)
      ((assertion (eq (magicl::matrix/double-float-layout m) :column-major))
       (assertion (= 2 (magicl::matrix/double-float-nrows m)))
       (assertion (= 2 (magicl::matrix/double-float-ncols m))))
    (aref (fast-storage m) (+ row (the fixnum (* 2 col))))))


(defun (setf varef) (new-value m i)
  (setf (aref (fast-storage m) i) new-value))

(defun (setf mtref) (value m row col)
  ;; (setf (magicl:tref m row col) value)
  (declare (magicl:matrix/double-float m)
           (fixnum row col))
  (assert (eq (magicl::matrix/double-float-layout m) :column-major))
  (let ((numrows (magicl::matrix-nrows m)))
    (declare (type fixnum numrows))
    (setf (aref (fast-storage m) (+ row (the fixnum (* col numrows)))) value))
  )

(declaim (inline vector-tref)
         (ftype (function (magicl::matrix/double-float fixnum) double-float)
                vector-tref))
(defun vector-tref (v index)
  (the double-float (aref (fast-storage v) index)))

;; (dotimes (i 100))
;; (let ((mat (cl-mpm/utils:voigt-to-matrix (cl-mpm/utils:voigt-from-list (loop repeat 6 collect (random 1d0))))))
;;   (multiple-value-bind (l v) (magicl:hermitian-eig mat)
;;     (pprint l)
;;     (pprint v)
;;     ;; (sort l #'>)
;;     ;; (multiple-value-bind (lu vu) (eig mat)
;;     ;;   (sort lu #'>)
;;     ;;   (when (not (every (lambda (x) (< (abs x) 1d-10)) (mapcar #'- l lu)))
;;     ;;     (format t "~A - ~A ~%" l lu))
;;     ;;   (pprint (every (lambda (x) (< (abs x) 1d-10)) (mapcar #'- l lu))))
;;     ))


(defun @-mat-vec (mat vec)
  (let ((result (vector-zeros)))
    (magicl:mult mat vec :target result)
    result))

(defun @-tensor-voigt (mat vec)
  (let ((result (voigt-zeros)))
    (magicl:mult mat vec :target result)
    result))

(let ((vector-components '(:x :y :z)))
  (defun get-vector (vector component)
    "Get stress component from voigt notation stress vector"
    (policy-cond:with-expectations (> speed safety)
        ((assertion (position component vector-components))
         (assertion (= (length (fast-storage vector)) 3))))
    (aref (magicl::matrix/double-float-storage vector) (position component vector-components))))

(let ((stress-components '(:xx :yy :zz
                           :yz :xz :xy)))
  (defun get-stress (stress component)
    "Get stress component from voigt notation stress vector"
    (policy-cond:with-expectations (> speed safety)
        ((assertion (position component stress-components))
         (assertion (= (length (fast-storage stress)) 6))))
    (aref (magicl::matrix/double-float-storage stress) (position component stress-components))))


(declaim (inline deep-copy)
         (ftype (function (magicl:matrix/double-float) magicl:matrix/double-float) deep-copy))
(defun deep-copy (m)
  "Deep copy a tensor"
  (let* ((rows (magicl::matrix/double-float-nrows m))
         (cols (magicl::matrix/double-float-ncols m))
         (total (* rows cols)))
    (declare (fixnum rows cols total))
    (magicl::make-matrix/double-float rows cols total :column-major (make-array total :element-type 'double-float
                                                                                      :initial-contents (fast-storage m)))))

(declaim (inline empty-copy)
         (ftype (function (magicl:matrix/double-float) magicl:matrix/double-float) empty-copy))
(defun empty-copy (m)
  "Deep copy a tensor"
  (let* ((rows (magicl::matrix/double-float-nrows m))
         (cols (magicl::matrix/double-float-ncols m))
         (total (* rows cols)))
    (declare (fixnum rows cols total))
    (magicl::make-matrix/double-float rows cols total :column-major
                                      (make-array total :element-type 'double-float
                                                        :initial-element 0d0))))

(declaim (inline vector-zeros)
         (ftype (function ()
                          magicl:matrix/double-float) vector-zeros))
(defun vector-zeros ()
  "Make a vector of 3x1 zeros"
  (magicl::make-matrix/double-float 3 1 3 :column-major (make-array 3 :element-type 'double-float)))

(declaim (inline voigt-zeros)
         (ftype (function ()
                          magicl:matrix/double-float) voigt-zeros))
(defun voigt-zeros ()
  "Make a stress-like vector 6x1 of zeros"
  (magicl::make-matrix/double-float 6 1 6 :column-major (make-array 6 :element-type 'double-float)))

(declaim (inline matrix-zeros)
         (ftype (function ()
                          magicl:matrix/double-float) matrix-zeros))
(defun matrix-zeros ()
  "Make a stress matrix 3x3 of zeros"
  (magicl::make-matrix/double-float 3 3 9 :column-major (make-array 9 :element-type 'double-float)))

(declaim (inline tensor-voigt-4th-zeros)
         (ftype (function ()
                          magicl:matrix/double-float) tensor-voigt-4th-zeros))
(defun tensor-voigt-4th-zeros ()
  "Make a stress matrix 3x3 of zeros"
  (magicl::make-matrix/double-float 6 6 36 :column-major (make-array 36 :element-type 'double-float)))

(declaim (inline matrix-copy)
         (ftype (function (magicl:matrix/double-float)
                          magicl:matrix/double-float) matrix-copy))
(defun matrix-copy (mat)
  (matrix-copy-into mat (matrix-zeros)))

(defun vector-copy (vec)
  (vector-copy-into vec (vector-zeros)))

(declaim (ftype (function (magicl::matrix/double-float
                           magicl::matrix/double-float) magicl::matrix/double-float) vector-copy-into))
(defun vector-copy-into (source target)
  "Copy a vector from source into target"
  (let ((s-s (cl-mpm/utils::fast-storage source))
        (t-s (cl-mpm/utils::fast-storage target)))
    (declare ((simple-array double-float (3)) s-s)
             ((simple-array double-float (3)) t-s))
    (loop for i from 0 below 3
          do (setf (aref t-s i) (aref s-s i))))
  target)

(declaim (ftype (function (magicl::matrix/double-float) magicl::matrix/double-float) voigt-copy))
(defun voigt-copy (vec)
  (voigt-copy-into vec (voigt-zeros)))

(declaim (ftype (function (magicl::matrix/double-float
                           magicl::matrix/double-float
                           ) magicl::matrix/double-float) voigt-copy-into))
(defun voigt-copy-into (source target)
  "Copy a vector from source into target"
  (let ((s-s (cl-mpm/utils::fast-storage source))
        (t-s (cl-mpm/utils::fast-storage target)))
    (declare ((simple-array double-float (6)) s-s)
             ((simple-array double-float (6)) t-s))
    (loop for i from 0 below 6
          do (setf (aref t-s i) (aref s-s i))))
   target)

(defun matrix-copy-into (source target)
  "Copy a vector from source into target"
  (let ((s-s (cl-mpm/utils::fast-storage source))
        (t-s (cl-mpm/utils::fast-storage target)))
    (declare ((simple-array double-float (9)) s-s)
             ((simple-array double-float (9)) t-s))
    (loop for i from 0 below 9
          do (setf (aref t-s i) (aref s-s i))))
   target)

(defun copy-into (source target)
  "Copy a matrix from source into target"
  (let ((s-s (cl-mpm/utils::fast-storage source))
        (t-s (cl-mpm/utils::fast-storage target)))
    (declare ((simple-array double-float (*)) s-s)
             ((simple-array double-float (*)) t-s))
    (loop for i from 0 below (length t-s)
          do (setf (aref t-s i) (aref s-s i))))
   target)

(defun voigt-contra->covar (vec)
  (let* ((v (voigt-copy vec))
        (vs (magicl::matrix/double-float-storage v)))
    (setf
     (aref vs 3) (* 0.5d0 (aref vs 3))
     (aref vs 4) (* 0.5d0 (aref vs 4))
     (aref vs 5) (* 0.5d0 (aref vs 5)))
    v
    ))
(defun voigt-covar->contra (vec)
  (let* ((v (voigt-copy vec))
         (vs (magicl::matrix/double-float-storage v)))
    (setf
     (aref vs 3) (* 2d0 (aref vs 3))
     (aref vs 4) (* 2d0 (aref vs 4))
     (aref vs 5) (* 2d0 (aref vs 5)))
    v))

(declaim (inline stretch-dsvp-zeros)
         (ftype (function ()
                          magicl:matrix/double-float) stretch-dsvp-zeros))
(defun stretch-dsvp-zeros ()
  (magicl::make-matrix/double-float 4 2 8 :column-major (make-array 8 :element-type 'double-float)))

(declaim (inline stretch-dsvp-3d-zeros)
         (ftype (function () magicl:matrix/double-float) stretch-dsvp-3d-zeros))
(defun stretch-dsvp-3d-zeros ()
  (magicl::make-matrix/double-float 9 3 27 :column-major (make-array 27 :element-type 'double-float)))

(declaim (inline stretch-dsvp-voigt-zeros)
         (ftype (function ()
                          magicl:matrix/double-float) stretch-dsvp-voigt-zeros))
(defun stretch-dsvp-voigt-zeros ()
  (magicl::make-matrix/double-float 9 1 9 :column-major (make-array 9 :element-type 'double-float)))

(declaim (inline dsvp-2d-zeros)
         (ftype (function ()
                          magicl:matrix/double-float) dsvp-2d-zeros))
(defun dsvp-2d-zeros ()
  (magicl::make-matrix/double-float 3 2 6 :column-major (make-array 6 :element-type 'double-float)))

(defun dsvp-3d-zeros ()
  (magicl::make-matrix/double-float 6 3 18 :column-major (make-array 18 :element-type 'double-float)))

(declaim (inline voigt-from-list)
         (ftype (function (list)
                          magicl:matrix/double-float) voigt-from-list))
(defun voigt-from-list (elements)
  "Voigt notation exx eyy ezz eyz ezx exy"
  (magicl::make-matrix/double-float 6 1 6 :column-major
                                    (make-array 6 :element-type 'double-float :initial-contents elements)))

(declaim (inline vector-from-list)
         (ftype (function (list)
                          magicl:matrix/double-float) vector-from-list))
(defun vector-from-list (elements)
  (magicl::make-matrix/double-float 3 1 3 :column-major
                                    (make-array 3 :element-type 'double-float :initial-contents elements)))


(declaim (inline matrix-from-list)
         (ftype (function (list)
                          magicl:matrix/double-float) matrix-from-list))
(defun matrix-from-list (elements)
  (magicl::make-matrix/double-float 3 3 9 :column-major
                                    (make-array 9 :element-type 'double-float :initial-contents elements)))

(declaim (inline tensor-voigt-4th-zeros)
         (ftype (function (list)
                          magicl:matrix/double-float) tensor-voigt-4th-from-list))
(defun tensor-voigt-4th-from-list (elements)
  "Make a stress matrix 3x3 of zeros"
  (magicl::make-matrix/double-float 6 6 36 :column-major (make-array 36 :element-type 'double-float
                                                                        :initial-contents elements
                                                                        )))
(defun tensor-voigt-4th-from-diag (elements)
  "Take a 3 element list and form a diagonal matrix"
  (destructuring-bind (i1 i2 i3 i4 i5 i6) elements
    (tensor-voigt-4th-from-list (list i1 0d0 0d0 0d0 0d0 0d0
                                      0d0 i2 0d0 0d0 0d0 0d0
                                      0d0 0d0 i3 0d0 0d0 0d0
                                      0d0 0d0 0d0 i4 0d0 0d0
                                      0d0 0d0 0d0 0d0 i5 0d0
                                      0d0 0d0 0d0 0d0 0d0 i6))))

(declaim (inline arb-matrix-from-list)
         (ftype (function (list fixnum fixnum)
                          magicl:matrix/double-float) arb-matrix-from-list))
(defun arb-matrix-from-list (elements x y)
  (declare (fixnum x y))
  (magicl::make-matrix/double-float
   x
   y
   (* x y)
   :column-major
   (make-array (* x y) :element-type 'double-float :initial-contents elements)))

(declaim (inline arb-matrix)
         (ftype (function (fixnum fixnum)
                          magicl:matrix/double-float) arb-matrix))
(defun arb-matrix (x y)
  (magicl::make-matrix/double-float x y (* x y) :column-major
                                    (make-array (* x y) :element-type 'double-float :initial-element 0d0)))

(declaim (inline matrix-from-diag)
         (ftype (function (list)
                          magicl:matrix/double-float) matrix-from-diag))
(defun matrix-from-diag (elements)
  "Take a 3 element list and form a diagonal matrix"
  (destructuring-bind (i1 i2 i3) elements
    (matrix-from-list (list i1 0d0 0d0
                            0d0 i2 0d0
                            0d0 0d0 i3))))


(declaim (inline matrix-to-voight)
         (ftype (function (magicl:matrix/double-float
                           &optional (or null magicl:matrix/double-float))
                          magicl:matrix/double-float
                          ) matrix-to-voight))
(defun matrix-to-voight (matrix &optional (result nil))
  "Stress matrix to voigt"
  (let ((result (if result result (voigt-zeros))))
        (let* ((exx (mtref matrix 0 0))
               (eyy (mtref matrix 1 1))
               (ezz (mtref matrix 2 2))
               (eyz (mtref matrix 2 1))
               (exy (mtref matrix 1 0))
               (ezx (mtref matrix 2 0)))
          (setf
           (varef result 0) exx
           (varef result 1) eyy
           (varef result 2) ezz
           (varef result 3) eyz
           (varef result 4) ezx
           (varef result 5) exy
           )
          result
          ;; (voigt-from-list (list exx eyy ezz eyz ezx exy))
          )))

(declaim (inline voight-to-matrix)
         (ftype (function (magicl:matrix/double-float)
                          magicl:matrix/double-float) voight-to-matrix))

(defun voight-to-matrix (vec)
  "Stress format voight to matrix"
  (let* ((exx (mtref vec 0 0))
         (eyy (mtref vec 1 0))
         (ezz (mtref vec 2 0))
         (eyz (mtref vec 3 0))
         (ezx (mtref vec 4 0))
         (exy (mtref vec 5 0)))
    (matrix-from-list (list exx exy ezx
                            exy eyy eyz
                            ezx eyz ezz))))

(declaim (inline voigt-to-matrix)
         (ftype (function (magicl:matrix/double-float)
                          magicl:matrix/double-float) voigt-to-matrix))
(defun voigt-to-matrix (vec)
  (let* ( (exx (mtref vec 0 0))
          (eyy (mtref vec 1 0))
          (ezz (mtref vec 2 0))
          (eyz (* 0.5d0 (the double-float (mtref vec 3 0))))
          (ezx (* 0.5d0 (the double-float (mtref vec 4 0))))
          (exy (* 0.5d0 (the double-float (mtref vec 5 0))))
          )
    (matrix-from-list (list exx exy ezx
                            exy eyy eyz
                            ezx eyz ezz))))

(declaim (inline matrix-to-voigt)
         (ftype (function (magicl:matrix/double-float)
                          magicl:matrix/double-float) matrix-to-voigt))
(defun matrix-to-voigt (matrix)
  (let* ( (exx (mtref matrix 0 0))
          (eyy (mtref matrix 1 1))
          (ezz (mtref matrix 2 2))
          (exy (mtref matrix 1 0))
          (exz (mtref matrix 2 0))
          (ezy (mtref matrix 2 1))
          )
    (voigt-from-list (list exx eyy ezz
                           (* 2d0 ezy)
                           (* 2d0 exz)
                           (* 2d0 exy)))))
(defun matrix-to-voigt-inplace (matrix vec)
  (let* ( (exx (mtref matrix 0 0))
          (eyy (mtref matrix 1 1))
          (ezz (mtref matrix 2 2))
          (exy (mtref matrix 1 0))
          (exz (mtref matrix 2 0))
          (ezy (mtref matrix 2 1))
          (s (magicl::matrix/double-float-storage vec))
          )
    (setf (aref s 0) exx)
    (setf (aref s 1) eyy)
    (setf (aref s 2) ezz)
    (setf (aref s 3) (* 2 ezy))
    (setf (aref s 4) (* 2 exz))
    (setf (aref s 5) (* 2 exy))))

(defun stretch-to-voight (matrix)
  (let* ( (exx (mtref matrix 0 0))
          (eyy (mtref matrix 1 1))
          (exy (mtref matrix 1 0))
          (eyx (mtref matrix 1 0))
          )
    (magicl:from-list (list exx eyy
                            exy eyx)
                      '(4 1) :type 'double-float)))

(defun voight-to-stretch (vec)
  (voight-to-stretch-prealloc vec (matrix-zeros)))

(declaim
 (inline voight-to-stretch-prealloc)
 (ftype (function (magicl:matrix/double-float magicl:matrix/double-float)
                  magicl:matrix/double-float) voight-to-stretch-prealloc))
(defun voight-to-stretch-prealloc (vec result)
  "Take a voigt matrix of stretches to "
  (let ((vecs (magicl::matrix/double-float-storage vec)))
    (let* ((exx (aref vecs 0))
           (eyy (aref vecs 1))
           (ezz (aref vecs 2))
           (exy (aref vecs 3))
           (eyx (aref vecs 4))
           (exz (aref vecs 5))
           (ezx (aref vecs 6))
           (eyz (aref vecs 7))
           (ezy (aref vecs 8))
           (s result)
           )

      (declare (double-float exx eyy ezz exy eyx ezy ezy exz ezx))
      (setf
       (mtref s 0 0) exx
       (mtref s 0 1) exy
       (mtref s 1 0) eyx
       (mtref s 1 1) eyy

       (mtref s 0 2) exz
       (mtref s 2 0) ezx

       (mtref s 2 1) ezy
       (mtref s 1 2) eyz

       (mtref s 2 2) ezz
       )))
  result)

(defun matrix-to-voight-strain (matrix)
  (let* ( (exx (mtref matrix 0 0))
          (eyy (mtref matrix 1 1))
          (exy (mtref matrix 1 0)))
    (magicl:from-list (list exx exy exy eyy)
                      '(4 1) :type 'double-float)))

(defun voight-to-matrix-strain (vec)
  (let* ( (exx (mtref vec 0 0))
          (eyy (mtref vec 1 0))
          (exy (* 0.5d0 (mtref vec 2 0))))
    (magicl:from-list (list exx exy
                            exy eyy)
                      '(2 2) :type 'double-float)))

(defun voigt-to-voight-strain (matrix)
  (let* ( (exx (mtref matrix 0 0))
          (eyy (mtref matrix 1 0))
          (exy (mtref matrix 2 0)))
    (magicl:from-list (list exx exy exy eyy)
                      '(4 1) :type 'double-float)))

(defun voigt-strain-to-voight (matrix)
  (let* ( (e0 (mtref matrix 0 0))
          (e1 (mtref matrix 1 0))
          (e3 (mtref matrix 3 0)))
    (magicl:from-list (list (- e0 e3) (- e1 e3) e3)
                      '(3 1) :type 'double-float)))

(defun voigt-to-mandel (voigt)
  (let* ((scale (sqrt 2)))
    (voigt-from-list
     (list 
      (varef voigt 0)
      (varef voigt 1)
      (varef voigt 2)

      (* scale (varef voigt 3))
      (* scale (varef voigt 4))
      (* scale (varef voigt 5))))))
(defun mandel-to-voigt (voigt)
  (let* ( (exx (mtref voigt 0 0))
          (eyy (mtref voigt 1 0))
          (exy (* (/ 1 (sqrt 2)) (mtref voigt 2 0))))
    (magicl:from-list (list exx eyy exy)
                      '(3 1) :type 'double-float)))

(defun matrix-to-mandel (matrix)
  (let* ( (exx (mtref matrix 0 0))
          (eyy (mtref matrix 1 1))
          (exy (* (sqrt 2) (mtref matrix 1 0))))
    (magicl:from-list (list exx eyy exy)
                      '(3 1) :type 'double-float)))
(defun mandel-to-matrix (vec)
  (let* ( (exx (mtref vec 0 0))
          (eyy (mtref vec 1 0))
          (exy (/ (mtref vec 2 0) (sqrt 2))))
    (magicl:from-list (list exx exy
                            exy eyy)
                      '(2 2) :type 'double-float)))

(defun stress-from-list (list)
  (voigt-from-list list))

(declaim (inline trace-voigt)
         (ftype (function (magicl:matrix/double-float) double-float) trace-voigt))
(defun trace-voigt (a)
  "Calculate the trace of A stored as a voigt tensor"
  (let ((arr (cl-mpm/utils:fast-storage a)))
    (declare ((simple-array double-float (6)) arr))
    (+ (aref arr 0) (aref arr 1) (aref arr 2))))

(defun trace-vector (a)
  "Calculate the trace of A as a vector of principal values"
  (let ((arr (cl-mpm/utils:fast-storage a)))
    (declare ((simple-array double-float (3)) arr))
    (+ (aref arr 0) (aref arr 1) (aref arr 2))))

(declaim (inline deviatoric-voigt)
         (ftype (function (magicl:matrix/double-float)
                          magicl:matrix/double-float) deviatoric-voigt))
(defun deviatoric-voigt (a)
  "Calculate the product A_{ij}A_{ij}"
  (let* ((tr (/ (trace-voigt a) 3d0))
         (arr (magicl::matrix/double-float-storage a)))
    (declare ((simple-array double-float (6)) arr)
             (double-float tr))
    (voigt-from-list (list (- (aref arr 0) tr)
                           (- (aref arr 1) tr)
                           (- (aref arr 2) tr)
                           (aref arr 3)
                           (aref arr 4)
                           (aref arr 5)))))
(declaim (inline deviatoric-vector)
         (ftype (function (magicl:matrix/double-float)
                          magicl:matrix/double-float) deviatoric-vector))
(defun deviatoric-vector (a)
  "Calculate the product A_{ij}A_{ij}"
  (let* ((tr (/ (trace-vector a) 3d0))
         (arr (magicl::matrix/double-float-storage a)))
    (declare ((simple-array double-float (3)) arr)
             (double-float tr))
    (vector-from-list (list (- (aref arr 0) tr)
                            (- (aref arr 1) tr)
                            (- (aref arr 2) tr)))))

(defun plane-strain-transform (stress)
  (magicl:from-list (list (mtref stress 0 0)
                          (mtref stress 1 0)
                          (mtref stress 5 0))
                    '(3 1)
                    :type 'double-float))


(declaim (inline stretch-to-sym)
         (ftype (function (magicl:matrix/double-float &optional magicl:matrix/double-float) (values)) stretch-to-sym))
(defun stretch-to-sym (stretch &optional result)
  (let ((result (if result result (voigt-zeros))))
    (declare (magicl:matrix/double-float result))
    (loop for i from 0 below 3
          do
             (setf (mtref result  i 0)
                   (mtref stretch i i)))
    (setf (mtref result 3 0)
          (* 2d0 0.5d0 (+ (the double-float (mtref stretch 2 1))
                      (the double-float (mtref stretch 1 2)))))
    (setf (mtref result 4 0)
          (* 2d0 0.5d0 (+ (the double-float (mtref stretch 0 2))
                      (the double-float (mtref stretch 2 0)))))
    (setf (mtref result 5 0)
          (* 2d0 0.5d0 (+ (the double-float (mtref stretch 0 1))
                      (the double-float (mtref stretch 1 0)))))
    result))


(declaim (inline stretch-to-skew)
         (ftype (function (magicl:matrix/double-float magicl:matrix/double-float) (values)) stretch-to-skew))
(defun stretch-to-skew (stretch result)
  ;; (magicl:.- stretch (magicl:transpose stretch) result)
  ;; (loop for i from 0 below 3
  ;;       do
  ;;          (setf (mtref result  i 0) 0d0))
  ;; (unless result
  ;;   (setf result (cl-mpm/utils:voigt-zeros)))
  (setf (mtref result  0 0)
        0)
  (setf (mtref result  1 0)
        0)
  (setf (mtref result  2 0)
        0)
  ;; Since off diagonal components get halved, then voigt doubles them this is net 1d0
  (setf (mtref result 3 0)
        (- (the double-float (mtref stretch 2 1))
                  (the double-float (mtref stretch 1 2))))
  (setf (mtref result 4 0)
        (- (the double-float (mtref stretch 0 2))
                  (the double-float (mtref stretch 2 0))))
  (setf (mtref result 5 0)
        (- (the double-float (mtref stretch 0 1))
                  (the double-float (mtref stretch 1 0))))

  (values))


(declaim
 (ftype (function (double-float) magicl::matrix/double-float) matrix-eye))
(defun matrix-eye (value)
  (matrix-from-list (list value 0d0 0d0
                          0d0 value 0d0
                          0d0 0d0 value)))
(declaim
 (ftype (function (list) magicl::matrix/double-float) matrix-diag))
(defun matrix-diag (values)
  (destructuring-bind (v1 v2 v3) values
    (declare (double-float v1 v2 v3))
    (matrix-from-list (list v1 0d0 0d0
                            0d0 v2 0d0
                            0d0 0d0 v3))))

(declaim
 (ftype (function (double-float) magicl::matrix/double-float) voigt-eye))
(defun voigt-eye (value)
  (voigt-from-list (list value value value 0d0 0d0 0d0)))

(defun rotation-matrix (degrees)
  (declare (double-float degrees))
  (let ((angle (/ (* pi degrees) 180)))
    (matrix-from-list (list (cos angle) (- (sin angle)) 0d0
                            (sin angle) (cos angle) 0d0
                            0d0 0d0 1d0))))

(defmacro time-form (it form)
  `(progn
     (declaim (optimize speed))
     (let* ((iterations ,it)
            (start (get-internal-real-time)))
       (dotimes (i ,it)
         ,form)
       (let* ((end (get-internal-real-time))
              (units internal-time-units-per-second)
              (dt (/ (- end start) (* iterations units)))
              )
         (format t "Total time: ~f ~%" (/ (- end start) units)) (format t "Time per iteration: ~f~%" (/ (- end start) (* iterations units)))
         (format t "Throughput: ~f~%" (/ 1 dt))
         dt))))

(defun principal-stresses-3d (stress)
  (multiple-value-bind (l v) (eig (voight-to-matrix stress))
    (declare (ignore v))
    (setf l (sort l #'>))
    (values (nth 0 l) (nth 1 l) (nth 2 l))))

(defun principal-strains-3d (strain)
  (multiple-value-bind (l v) (eig (voigt-to-matrix strain))
    (declare (ignore v))
    (setf l (sort l #'>))
    (values (nth 0 l) (nth 1 l) (nth 2 l))))

;; (defmacro time-form (it form)
;;   `(progn
;;      (declaim (optimize speed))
;;      (let* ((iterations ,it)
;;             (start (get-internal-real-time)))
;;        (time
;;         (dotimes (i ,it)
;;           ,form))
;;        (let* ((end (get-internal-real-time))
;;               (units internal-time-units-per-second)
;;               (dt (/ (- end start) (* iterations units)))
;;               )
;;          (format t "Total time: ~f ~%" (/ (- end start) units)) (format t "Time per iteration: ~f~%" (/ (- end start) (* iterations units)))
;;          (format t "Throughput: ~f~%" (/ 1 dt))
;;          (format t "Time per MP: ~E~%" (/ dt (length (cl-mpm:sim-mps *sim*))))
;;          dt))))

(defmacro track-time ((var) &body body)
  (alexandria:with-gensyms (timespan result)
    `(lret (,result)
           (let ((,timespan
                   (measure-time (_ (setq ,result (progn ,@body))) 1 t)))
             (when ,var
               (incf ,var ,timespan))))))

(defun d2r (degrees)
  (* degrees (/ pi 180)))
(defun r2d (radians)
  (* radians (/ 180 pi)))

(defun deg-to-rad (degrees)
  (* degrees (/ pi 180)))
(defun rad-to-deg (radians)
  (* radians (/ 180 pi)))


(defun slice-matrix-2d (mat)
  (let* ((nd 2)
         (m (arb-matrix nd nd)))
    (loop for x from 0 below nd
          do
             (loop for y from 0 below nd
                   do (setf (mtref m x y) (mtref mat x y))))
    m))
(defun pad-matrix-2d (mat)
  (let* ((nd 3)
         (m (arb-matrix nd nd)))
    (loop for x from 0 below 2
          do
             (loop for y from 0 below 2
                   do (setf (mtref m x y) (mtref mat x y))))
    m))
(defun pad-matrix-1d (mat)
  (let* ((nd 3)
         (m (arb-matrix nd nd)))
    (loop for x from 0 below 1
          do
             (loop for y from 0 below 1
                   do (setf (mtref m x y) (mtref mat x y))))
    m))

(defun slice-matrix-nd (mat nd)
  (ecase nd
    (1 (arb-matrix-from-list (list (mtref mat 0 0)) 1 1))
    (2 (slice-matrix-2d mat))
    (3 mat)))

(defun pad-matrix-nd (mat nd)
  (ecase nd
    (1 (pad-matrix-1d mat))
    (2 (pad-matrix-2d mat))
    (3 mat)))

(defun matrix-column (mat i)
  (magicl:vector->column-matrix (magicl:column mat i)))


(declaim (ftype (function (double-float double-float) double-float) calculate-bulk-modulus))
(defun calculate-bulk-modulus (E nu)
  (/ E (* 3 (- 1d0 (* 2 nu)))))
(declaim (ftype (function (double-float double-float) double-float) calculate-shear-modulus))
(defun calculate-shear-modulus (E nu)
  (/ E (* 2 (+ 1d0 nu))))
(declaim (ftype (function (double-float double-float) double-float) calculate-p-wave-modulus))
(defun calculate-p-wave-modulus (E nu)
  (/ (* (- 1d0 nu) E) (* (+ 1d0 nu) (- 1d0 (* 2d0 nu)))))


(defstruct object-pool
  (lock (sb-thread:make-mutex) :type sb-thread:mutex)
  (pool (make-array 0 :element-type t) :type (array t (*)))
  (constructor #'identity :type function))

(defun rebuild-object-pool (pool)
  ;;Trigger a re-build
  (sb-thread:with-mutex ((object-pool-lock pool))
    (let ((ws (get-worker-size)))
      (when (> ws (length (object-pool-pool pool)))
        (let* ((prev-pool (object-pool-pool pool))
               (new-pool (make-array ws :element-type t)))
          (loop for i from 0 below (length new-pool)
                do (setf (aref new-pool i)
                         (if (< i (length prev-pool))
                             (aref prev-pool i)
                             (funcall (object-pool-constructor pool)))))

          (setf (object-pool-pool pool) new-pool)
          )))))

(defun object-pool-grab (pool)
  (let ((thread-index (get-worker-index)))
    (progn
      (when (>= (lparallel:kernel-worker-count) (length (object-pool-pool pool)))
        (rebuild-object-pool pool))
      (aref (object-pool-pool pool) thread-index))))

(defun get-worker-size ()
  (let ((ls (lparallel:kernel-worker-count)))
    (if ls
        (the fixnum (1+ (the fixnum ls)))
        1)))

(defun object-pool-ensure-size (pool)
  (let ((ws (get-worker-size)))
    (when (> ws (length (object-pool-pool pool)))
      ;;Trigger a re-build
      (sb-thread:with-mutex ((object-pool-lock pool))
        (when (> ws (length (object-pool-pool pool)))
          (let* ((prev-pool (object-pool-pool pool))
                 (wc (the fixnum ws))
                 (new-pool (make-array (1+ wc) :element-type t)))
            (loop for i from 0 below (length new-pool)
                  do (setf (aref new-pool i)
                           (if (< i (length prev-pool))
                               (aref prev-pool i)
                               (funcall (object-pool-constructor pool)))))
            (setf (object-pool-pool pool) new-pool)))))))

(defun get-worker-index ()
  (let ((ls (or (lparallel:kernel-worker-index)
                cl-mpm/utils::*worker-index*)))
    (if ls
        (1+ ls)
        1)))

(defun object-pool-grab-unsafe (pool)
  (let ((thread-index (get-worker-index)))
    (aref (object-pool-pool pool) thread-index)))

(defstruct sparse-matrix
  "Sparse self-adjoint matrix"
  values
  cols
  rowindex
  nrows
  ncols)

(defstruct sparse-matrix-ccs
  "Sparse matrix compressed column storage"
  values
  rows
  colindex
  nrows
  ncols)


(defun sparse-matrix-aref (mat row col)
  (let* ((rowindex (cl-mpm/utils::sparse-matrix-rowindex mat))
         (cols (cl-mpm/utils::sparse-matrix-cols mat))
         (values (cl-mpm/utils::sparse-matrix-values mat))
         (result 0d0)
         )
    (let* ((r row)
           (col-0 (aref rowindex r))
           (col-1 (aref rowindex (1+ r))))
      (loop for c from col-0 below col-1
            while (<= (aref cols c) col)
            do
               (when (= col (aref cols c))
                 (setf result
                       (the double-float
                            (aref values c))))))
    result))


(declaim (notinline sort-and-sum))
(defun sort-and-sum (values rows cols)
  (let ((len (length values)))
    (let* ((permutation (make-array len :element-type 'fixnum)))
      (loop for i from 0 below (length values)
            do (setf (aref permutation i) i))
      ;; (format t "Compacting ~D fds~%" len)
      ;; (if (> len 10000)
      ;;   (lparallel:psort permutation
      ;;                    (lambda (i j)
      ;;                      (or (< (aref rows i)
      ;;                             (aref rows j))
      ;;                          (and
      ;;                           (= (aref rows i)
      ;;                              (aref rows j))
      ;;                           (< (aref cols i)
      ;;                              (aref cols j)))))
      ;;                    :granularity (floor (* len 0.1)))
      ;;   (lparallel:psort permutation
      ;;                    (lambda (i j)
      ;;                      (or (< (aref rows i)
      ;;                             (aref rows j))
      ;;                          (and
      ;;                           (= (aref rows i)
      ;;                              (aref rows j))
      ;;                           (< (aref cols i)
      ;;                              (aref cols j)))))))
      (setf permutation
            (sort permutation
                  (lambda (i j)
                    (or (< (aref rows i)
                           (aref rows j))
                        (and
                         (= (aref rows i)
                            (aref rows j))
                         (< (aref cols i)
                            (aref cols j)))))))
      ;; (sort
      ;;  permutation
      ;;  (lambda (i j)
      ;;    (or (< (aref rows i)
      ;;           (aref rows j))
      ;;        (and
      ;;         (= (aref rows i)
      ;;            (aref rows j))
      ;;         (< (aref cols i)
      ;;            (aref cols j))))))
      (let* ((sorted-values (make-array len :element-type 'double-float :fill-pointer len))
             (sorted-cols   (make-array len :element-type 'fixnum :fill-pointer len))
             (sorted-rows   (make-array len :element-type 'fixnum :fill-pointer len)))
        (loop for i from 0 below len
              for ip across permutation
              do
                 (setf
                  (aref sorted-values i) (aref values ip)
                  (aref sorted-cols   i) (aref cols   ip)
                  (aref sorted-rows   i) (aref rows   ip)))
        (let ((rolling-index -1)
              (row -1)
              (col -1)
              (sum 0d0))
          (loop for i from 0 below len
                for v across sorted-values
                for c across sorted-cols
                for r across sorted-rows
                do
                   (if (and (= row r) (= col c))
                       (incf sum v)
                       (progn
                         (unless (= rolling-index -1)
                           (setf (aref sorted-values rolling-index) sum
                                 (aref sorted-cols   rolling-index) col
                                 (aref sorted-rows   rolling-index) row))
                         (setf sum v)
                         (setf row r
                               col c)
                         (incf rolling-index))))
          (setf (aref sorted-values rolling-index) sum
                (aref sorted-cols   rolling-index) col
                (aref sorted-rows   rolling-index) row)
          (incf rolling-index)
          (setf (fill-pointer sorted-values) rolling-index
                (fill-pointer sorted-cols)   rolling-index
                (fill-pointer sorted-rows)   rolling-index))
        (values
         (subseq sorted-values 0)
         (subseq sorted-rows 0)
         (subseq sorted-cols 0))))))

(defun build-sparse-matrix (values rows cols nrows ncols)
  "Take a triplet of vectors of (value row col) and make a sparse matrix in compressed row storage"
  (multiple-value-bind (v cr cc) (cl-mpm/utils::sort-and-sum values rows cols)
    (let ((rolling-index 0)
          (row 0)
          (rowindex (make-array (1+ nrows) :element-type 'fixnum :initial-element 0)))
      (loop for i from 0
            for r across cr
            do
               (progn
                 (unless (= row r)
                   (loop for j from (1+ row) below r
                         when (> j 0)
                           do (setf (aref rowindex j) rolling-index))
                   (setf (aref rowindex r) rolling-index)
                   (setf row r))
                 (incf rolling-index)))
      (setf (aref rowindex (1+ row)) rolling-index)
      (make-sparse-matrix :values v
                          :cols cc
                          :rowindex rowindex
                          :nrows nrows
                          :ncols ncols))))
(defun mat-to-sparse (mat)
  (let ((v (make-array 0 :adjustable t :fill-pointer 0 :element-type 'double-float))
        (r (make-array 0 :adjustable t :fill-pointer 0 :element-type 'fixnum))
        (c (make-array 0 :adjustable t :fill-pointer 0 :element-type 'fixnum)))
    (loop for x below (magicl:nrows mat) do
      (loop for y below (magicl:ncols mat)
            when (> (abs (cl-mpm/utils::mtref mat x y)) 0d0)
              do (progn
                   (vector-push-extend (cl-mpm/utils::mtref mat x y) v)
                   (vector-push-extend x r)
                   (vector-push-extend y c))))
    (build-sparse-matrix v r c (magicl:nrows mat) (magicl:ncols mat))))
(defun sparse-to-mat (mat)
  (let ((rowindex (cl-mpm/utils::sparse-matrix-rowindex mat))
        (cols (cl-mpm/utils::sparse-matrix-cols mat))
        (values (cl-mpm/utils::sparse-matrix-values mat))
        (res (cl-mpm/utils::arb-matrix (cl-mpm/utils::sparse-matrix-nrows mat) (cl-mpm/utils::sparse-matrix-ncols mat))))
    (dotimes (r (cl-mpm/utils::sparse-matrix-nrows mat))
      (let ((col-0 (aref rowindex r))
            (col-1 (aref rowindex (1+ r))))
        (loop for c from col-0 below col-1 do
          (setf (mtref res r (aref cols c))
                (the double-float
                     (aref values c))))))
    res))

(defun sparse-to-coordinates (mat)
  (let ((coords (list)))
    (let ((rowindex (cl-mpm/utils::sparse-matrix-rowindex mat))
          (cols (cl-mpm/utils::sparse-matrix-cols mat))
          (values (cl-mpm/utils::sparse-matrix-values mat)))
      (dotimes (r (cl-mpm/utils::sparse-matrix-nrows mat))
        (let ((col-0 (aref rowindex r))
              (col-1 (aref rowindex (1+ r))))
          (loop for c from col-0 below col-1 do
            (push
             (list
              (the double-float
                   (aref values c))
              r
              (aref cols c)
              )
             coords))))
      ;; (sort coords
      ;;       (lambda (i j)
      ;;         (or (< (nth 2 i)
      ;;                (nth 2 j))
      ;;             (and
      ;;              (= (nth 2 i)
      ;;                 (nth 2 j))
      ;;              (< (nth 1 i)
      ;;                 (nth 1 j))))))
      coords
      ;; (sort coords
      ;;       (lambda (i j)
      ;;         (or (< (nth 1 i)
      ;;                (nth 1 j))
      ;;             (and
      ;;              (= (nth 1 i)
      ;;                 (nth 1 j))
      ;;              (< (nth 2 i)
      ;;                 (nth 2 j))))))
      )))

(defun build-sparse-matrix-ccs (values rows cols nrows ncols)
  "Take a triplet of vectors of (value row col) and make a sparse matrix in compressed row storage"
  (multiple-value-bind (v cc cr) (cl-mpm/utils::sort-and-sum values cols rows)
    (let ((rolling-index 0)
          (col -1)
          (colindex (make-array (1+ ncols) :element-type 'fixnum :initial-element 0)))
      (loop for i from 0
            for r across cc
            do
               (progn
                 (unless (= col r)
                   (setf (aref colindex r) rolling-index)
                   (setf col r))
                 (incf rolling-index)))
      (setf (aref colindex (1+ col)) rolling-index)
      (make-sparse-matrix-ccs :values v
                          :rows cr
                          :colindex colindex
                          :nrows nrows
                          :ncols ncols))))



(declaim (inline make-gradients))
(defstruct (gradients (:constructor make-gradients (dx dy dz)))
  (dx 0d0 :type double-float)
  (dy 0d0 :type double-float)
  (dz 0d0 :type double-float))

(declaim (inline gradients-to-list))
(defun gradients-to-list (grads)
  (list (gradients-dx grads)
        (gradients-dy grads)
        (gradients-dz grads)))

(declaim (inline gradients-to-vector))
(defun gradients-to-vector (grads)
  (cl-mpm/utils:vector-from-list (list (gradients-dx grads)
                                       (gradients-dy grads)
                                       (gradients-dz grads))))

;; (defstruct mutex-vector
;;   (length 0)
;;   (vector (make-vector 2 :fill-pointer 0 :adjustable nil :element-type t))
;;   (mutex (sb-thread:make-mutex)))

;; (defun mutex-vector-reset (vec))
;; (defun mutex-vector-push-extend (value vec)
;;   ;; (when (< (+ (mutex-vector-length vec) 1)))
;;   ;; ()
;;   )





(defconstant +thread-parts-scale+ 1)
(defun get-parts ()
  (the fixnum (* (the fixnum +thread-parts-scale+) *worker-count*)))

(defparameter *workers* nil)
(defparameter *worker-count* 0)
(declaim (fixnum *worker-count*))
(defparameter *workers-array* nil)
(defparameter *work-queue* (sb-concurrency:make-queue))
(defparameter *workers-run* (sb-thread:make-semaphore))
(defparameter *workers-finish* (sb-thread:make-semaphore))
(defparameter *workers-kill* nil)
(defparameter *workers-func* (lambda (i len)))
(declaim (function *workers-func*))
(defparameter *workers-chunk* 1)
(defparameter *workers-chunk-count* 0)
(defparameter *workers-array-length* 0)
(defparameter *workers-nesting* nil)
(defparameter *workers-pool-age* 0)
(defparameter *worker-index* nil)
(declaim (fixnum *workers-chunk* *workers-chunk-count* *workers-array-length* *workers-pool-age*))
(defparameter *workers-counter* (make-array 1 :element-type '(unsigned-byte 64)))
(declaim ((simple-array (unsigned-byte 64) 1) *workers-counter*))
(progn
  ;; (kill-workers)
  (defun make-workers (&key (worker-count nil))
    (declare (function *workers-func*)
             (fixnum *workers-chunk-count*))
    (when *workers*
      (kill-workers))
    (let ((thread-count (if worker-count worker-count (lparallel:kernel-worker-count))))
      (incf *workers-pool-age*)
      (setf *workers-nesting* nil)
      (setf *worker-count* thread-count)
      (setf *workers-run* (sb-thread:make-semaphore))
      (setf *workers-finish* (sb-thread:make-semaphore))
      (let ((current-pool-age *workers-pool-age*))
        (setf *workers*
              (loop for i fixnum from 0 below thread-count
                    collect
                    (sb-thread:make-thread
                     (lambda (thread-number current-age)
                       (declare (fixnum current-age))
                       (let (;; (lparallel.kernel::*worker*
                             ;;   (lparallel.kernel::make-worker-instance
                             ;;    :thread nil
                             ;;    :tasks (lparallel.kernel::make-spin-queue)
                             ;;    :index thread-number))
                             (*worker-index* thread-number)
                             )
                         (loop
                           while (not *workers-kill*)
                           do
                              (progn
                                (sb-thread:wait-on-semaphore *workers-run*)
                                (if (= *workers-pool-age* current-age)
                                  (progn
                                    (unless *workers-kill*
                                      ;; (unwind-protect)
                                      (let ((iter (sb-ext:atomic-incf (aref *workers-counter* 0))))
                                        (when (< iter *workers-chunk-count*)
                                          (funcall *workers-func* iter)))
                                      ;; (print "cleanup unexpected thread termination")
                                      ;; (setf *workers-nesting* nil)
                                      ;; (setf *workers-kill* t)
                                      ))
                                  (sb-thread:signal-semaphore *workers-run*))
                                (sb-thread:signal-semaphore *workers-finish*)
                                )))
                       (values))
                     :arguments (list i current-pool-age))))))
    ))
(defun kill-workers ()
  (when *workers*
    (setf *workers-nesting* nil)
    (setf *workers-kill* t)
    (sb-thread:signal-semaphore *workers-run* *worker-count*)
    (ignore-errors
     (loop for worker in *workers* do (sb-thread:join-thread worker)))
    (setf *workers* nil)
    (setf *workers-kill* nil)
    ))
(defun set-workers (threads)
  (kill-workers)
  (lparallel:end-kernel)
  (setf lparallel::*kernel* (lparallel::make-kernel threads))
  (make-workers :worker-count threads)
  )
(declaim (ftype (function (fixnum function &key (:parts (or null fixnum))) (values))
                 omp))
(defun omp (total-length func &key (parts nil))
  (declare (fixnum total-length)
           (function func))
  ;; (make-workers)
  (unless *workers*
    (error "no workers"))
  (if *workers-nesting*
      (progn
        (loop for j fixnum from 0 below total-length
              do (funcall func j)))
      (progn
        (setf *workers-nesting* t)
        (setf *workers-array-length* total-length)
        (let* ((length (if parts parts *worker-count*))
               (block-size (ceiling total-length length)))
          (declare (fixnum length block-size))
          (setf *workers-func*
                (lambda (chunk-number)
                  (declare (fixnum chunk-number))
                  (let* ((start (* chunk-number block-size))
                         (end (min (the fixnum (* (1+ chunk-number) block-size)) total-length)))
                    (declare (fixnum block-size start end))
                    (loop for j fixnum from start below end
                          do (funcall func j))
                    (values))))
          (setf *workers-chunk-count* length)
          (sb-ext:atomic-update (aref *workers-counter* 0) (lambda (a) (declare (ignore a)) 0))
          (sb-thread:signal-semaphore *workers-run* length)
          (sb-thread:wait-on-semaphore *workers-finish* :n length))
        (setf *workers-nesting* nil)))
  (values))

(defmacro bpdotimes ((i length) &body func)
  ;; `(lparallel:pdotimes (,i ,length)
  ;;    ,@func
  ;;    )
  (let ((job-sym (gensym)))
    `(progn
       (let* ((total-size ,length)
              (,job-sym (lambda (,i)
                          (declare (fixnum ,i))
                          ,@func
                          )))
         (declare (function ,job-sym))
         (cl-mpm/utils::omp
          total-size
          ,job-sym
          :parts (get-parts)
          ))))
  )
