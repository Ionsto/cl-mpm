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
   #:matrix-copy
   #:vector-copy
   #:voigt-copy
   #:matrix-copy-into
   #:vector-copy-into
   #:voigt-copy-into
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
  ;; (multiple-value-bind (l v) 
  ;;   (values l v))
  ;;
  )

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


(let ((stress-components '(:xx :yy :zz
                           :yz :xz :xy)))
  (defun get-stress (stress component)
    "Get stress component from voigt notation stress vector"
    (aref (magicl::matrix/double-float-storage stress) (position component stress-components))))


(declaim (inline deep-copy)
         (ftype (function ()
                          magicl:matrix/double-float) vector-zeros))
(defun deep-copy (m)
  "Make a vector of 3x1 zeros"
  (let* ((rows (magicl::matrix/double-float-nrows m))
         (cols (magicl::matrix/double-float-ncols m))
         (total (* rows cols)))
    (declare (fixnum rows cols total))
    (magicl::make-matrix/double-float rows cols total :column-major (make-array total :element-type 'double-float))))

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
  (let ((v (matrix-zeros)))
    (aops:copy-into (magicl::matrix/double-float-storage v)
                    (magicl::matrix/double-float-storage mat))
    v))

(defun vector-copy (vec)
  (let ((v (vector-zeros)))
    (aops:copy-into (magicl::matrix/double-float-storage v)
                    (magicl::matrix/double-float-storage vec))
    v)
  )

(defun voigt-copy (vec)
  (let ((v (voigt-zeros)))
    (aops:copy-into (magicl::matrix/double-float-storage v)
                    (magicl::matrix/double-float-storage vec))
    v))

(defun voigt-copy-into (source target)
  (aops:copy-into (magicl::matrix/double-float-storage target)
                  (magicl::matrix/double-float-storage source)
                  ) target)

(defun matrix-copy-into (source target)
  (aops:copy-into (magicl::matrix/double-float-storage target)
                  (magicl::matrix/double-float-storage source)
                  ) target)

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


(declaim (inline matrix-to-voight)
         (ftype (function (magicl:matrix/double-float)
                          magicl:matrix/double-float) matrix-to-voight))
(defun matrix-to-voight (matrix)
  "Stress matrix to voigt"
  (let* ( (exx (magicl:tref matrix 0 0))
          (eyy (magicl:tref matrix 1 1))
          (ezz (magicl:tref matrix 2 2))
          (eyz (magicl:tref matrix 2 1))
          (exy (magicl:tref matrix 1 0))
          (ezx (magicl:tref matrix 2 0))
          )
    (voigt-from-list (list exx eyy ezz eyz ezx exy))))

(declaim (inline voight-to-matrix)
         (ftype (function (magicl:matrix/double-float)
                          magicl:matrix/double-float) voight-to-matrix))

(defun voight-to-matrix (vec)
  "Stress format voight to matrix"
  (let* ((exx (magicl:tref vec 0 0))
         (eyy (magicl:tref vec 1 0))
         (ezz (magicl:tref vec 2 0))
         (eyz (magicl:tref vec 3 0))
         (ezx (magicl:tref vec 4 0))
         (exy (magicl:tref vec 5 0)))
    (matrix-from-list (list exx exy ezx
                            exy eyy eyz
                            ezx eyz ezz))))

(declaim (inline voigt-to-matrix)
         (ftype (function (magicl:matrix/double-float)
                          magicl:matrix/double-float) voigt-to-matrix))
(defun voigt-to-matrix (vec)
  (let* ( (exx (magicl:tref vec 0 0))
          (eyy (magicl:tref vec 1 0))
          (ezz (magicl:tref vec 2 0))
          (eyz (* 0.5d0 (the double-float (magicl:tref vec 3 0))))
          (ezx (* 0.5d0 (the double-float (magicl:tref vec 4 0))))
          (exy (* 0.5d0 (the double-float (magicl:tref vec 5 0))))
          )
    (matrix-from-list (list exx exy ezx
                            exy eyy eyz
                            ezx eyz ezz))))

(declaim (inline matrix-to-voigt)
         (ftype (function (magicl:matrix/double-float)
                          magicl:matrix/double-float) matrix-to-voigt))
(defun matrix-to-voigt (matrix)
  (let* ( (exx (magicl:tref matrix 0 0))
          (eyy (magicl:tref matrix 1 1))
          (ezz (magicl:tref matrix 2 2))
          (exy (magicl:tref matrix 1 0))
          (exz (magicl:tref matrix 2 0))
          (ezy (magicl:tref matrix 2 1))
          )
    (voigt-from-list (list exx eyy ezz
                           (* 2d0 ezy)
                           (* 2d0 exz)
                           (* 2d0 exy)))))
(defun matrix-to-voigt-inplace (matrix vec)
  (let* ( (exx (magicl:tref matrix 0 0))
          (eyy (magicl:tref matrix 1 1))
          (ezz (magicl:tref matrix 2 2))
          (exy (magicl:tref matrix 1 0))
          (exz (magicl:tref matrix 2 0))
          (ezy (magicl:tref matrix 2 1))
          (s (magicl::matrix/double-float-storage vec))
          )
    (setf (aref s 0) exx)
    (setf (aref s 1) eyy)
    (setf (aref s 2) ezz)
    (setf (aref s 3) (* 2 ezy))
    (setf (aref s 4) (* 2 exz))
    (setf (aref s 5) (* 2 exy))))

(defun stretch-to-voight (matrix)
  (let* ( (exx (magicl:tref matrix 0 0))
          (eyy (magicl:tref matrix 1 1))
          (exy (magicl:tref matrix 1 0))
          (eyx (magicl:tref matrix 1 0))
          )
    (magicl:from-list (list exx eyy
                            exy eyx)
                      '(4 1) :type 'double-float)))

(defun voight-to-stretch (vec)
  (voight-to-stretch-prealloc vec (matrix-zeros)))

(declaim
 (notinline voight-to-stretch-prealloc)
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

      (declare (double-float exx eyy exy eyx ezy ez))
      (setf
       (magicl:tref s 0 0) exx
       (magicl:tref s 0 1) exy
       (magicl:tref s 1 0) eyx
       (magicl:tref s 1 1) eyy

       (magicl:tref s 0 2) exz
       (magicl:tref s 2 0) ezx

       (magicl:tref s 2 1) ezy
       (magicl:tref s 1 2) eyz

       (magicl:tref s 2 2) ezz
       )))
  result)

(defun matrix-to-voight-strain (matrix)
  (let* ( (exx (magicl:tref matrix 0 0))
          (eyy (magicl:tref matrix 1 1))
          (exy (magicl:tref matrix 1 0)))
    (magicl:from-list (list exx exy exy eyy)
                      '(4 1) :type 'double-float)))

(defun voight-to-matrix-strain (vec)
  (let* ( (exx (magicl:tref vec 0 0))
          (eyy (magicl:tref vec 1 0))
          (exy (* 0.5d0 (magicl:tref vec 2 0))))
    (magicl:from-list (list exx exy
                            exy eyy)
                      '(2 2) :type 'double-float)))

(defun voigt-to-voight-strain (matrix)
  (let* ( (exx (magicl:tref matrix 0 0))
          (eyy (magicl:tref matrix 1 0))
          (exy (magicl:tref matrix 2 0)))
    (magicl:from-list (list exx exy exy eyy)
                      '(4 1) :type 'double-float)))

(defun voigt-strain-to-voight (matrix)
  (let* ( (e0 (magicl:tref matrix 0 0))
          (e1 (magicl:tref matrix 1 0))
          (e3 (magicl:tref matrix 3 0)))
    (magicl:from-list (list (- e0 e3) (- e1 e3) e3)
                      '(3 1) :type 'double-float)))

(defun voigt-to-mandel (voigt)
  (let* ( (exx (magicl:tref voigt 0 0))
          (eyy (magicl:tref voigt 1 0))
          (exy (* (/ (sqrt 2) 1) (magicl:tref voigt 2 0))))
    (magicl:from-list (list exx eyy exy)
                      '(3 1) :type 'double-float)))
(defun mandel-to-voigt (voigt)
  (let* ( (exx (magicl:tref voigt 0 0))
          (eyy (magicl:tref voigt 1 0))
          (exy (* (/ 1 (sqrt 2)) (magicl:tref voigt 2 0))))
    (magicl:from-list (list exx eyy exy)
                      '(3 1) :type 'double-float)))

(defun matrix-to-mandel (matrix)
  (let* ( (exx (magicl:tref matrix 0 0))
          (eyy (magicl:tref matrix 1 1))
          (exy (* (sqrt 2) (magicl:tref matrix 1 0))))
    (magicl:from-list (list exx eyy exy)
                      '(3 1) :type 'double-float)))
(defun mandel-to-matrix (vec)
  (let* ( (exx (magicl:tref vec 0 0))
          (eyy (magicl:tref vec 1 0))
          (exy (/ (magicl:tref vec 2 0) (sqrt 2))))
    (magicl:from-list (list exx exy
                            exy eyy)
                      '(2 2) :type 'double-float)))

(defun stress-from-list (list)
  (voigt-from-list list))

(declaim (inline trace-voigt)
         (ftype (function (magicl:matrix/double-float) double-float) trace-voigt))
(defun trace-voigt (a)
  "Calculate the product A_{ij}A_{ij}"
  (let ((arr (magicl::matrix/double-float-storage a)))
    (declare ((simple-array double-float) arr))
    (+ (aref arr 0) (aref arr 1) (aref arr 2) )))

(declaim (inline deviatoric-voigt)
         (ftype (function (magicl:matrix/double-float)
                          magicl:matrix/double-float) deviatoric-voigt))
(defun deviatoric-voigt (a)
  "Calculate the product A_{ij}A_{ij}"
  (let* ((tr (/ (trace-voigt a) 3d0))
         (arr (magicl::matrix/double-float-storage a)))
    (declare ((simple-array double-float *) arr)
             (double-float tr))
    (voigt-from-list (list (- (aref arr 0) tr)
                           (- (aref arr 1) tr)
                           (- (aref arr 2) tr)
                           (aref arr 3)
                           (aref arr 4)
                           (aref arr 5)))))

(defun plane-strain-transform (stress)
  (magicl:from-list (list (magicl:tref stress 0 0)
                          (magicl:tref stress 1 0)
                          (magicl:tref stress 5 0))
                    '(3 1)
                    :type 'double-float))


(declaim (inline stretch-to-sym)
         (ftype (function (magicl:matrix/double-float magicl:matrix/double-float) (values)) stretch-to-sym))
(defun stretch-to-sym (stretch result)
  (progn
    (declaim (magicl:matrix/double-float result))
    ;; (let ((res (matrix-to-voight (magicl:scale! (magicl:.+ stretch (magicl:transpose stretch)) 0.5d0))))
    ;;   (aops:copy-into (magicl::matrix/double-float-storage result)
    ;;                   (magicl::matrix/double-float-storage res))
    ;;   result)
    (loop for i from 0 below 3
          do
             (setf (magicl:tref result  i 0)
                   (magicl:tref stretch i i)))
    (setf (magicl:tref result 3 0)
          (* 1d0 (+ (the double-float (magicl:tref stretch 2 1))
                    (the double-float (magicl:tref stretch 1 2)))))
    (setf (magicl:tref result 4 0)
          (* 1d0 (+ (the double-float (magicl:tref stretch 0 2))
                    (the double-float (magicl:tref stretch 2 0)))))
    (setf (magicl:tref result 5 0)
          (* 1d0 (+ (the double-float (magicl:tref stretch 0 1))
                    (the double-float (magicl:tref stretch 1 0)))))
    ;; (setf (magicl:tref result  i 0) (magicl:tref result  i 0))
    ;; (setf (magicl:tref result  0 0)
    ;;       (magicl:tref stretch 0 0))

    ;; (setf (magicl:tref result  1 0)
    ;;       (magicl:tref stretch 1 1))
    ;; ;;Since off diagonal components get halved, then voigt doubles them this is net 1d0
    ;; (setf (magicl:tref result 5 0)
    ;;       (* 1d0 (+ (the double-float (magicl:tref stretch 0 1))
    ;;                 (the double-float (magicl:tref stretch 1 0)))))
    )
  (values))


(declaim (inline stretch-to-skew)
         (ftype (function (magicl:matrix/double-float magicl:matrix/double-float) (values)) stretch-to-skew))
(defun stretch-to-skew (stretch result)
  ;; (magicl:.- stretch (magicl:transpose stretch) result)
  ;; (loop for i from 0 below 3
  ;;       do
  ;;          (setf (magicl:tref result  i 0) 0d0))
  ;; (unless result
  ;;   (setf result (cl-mpm/utils:voigt-zeros)))
  (setf (magicl:tref result  0 0)
        0)
  (setf (magicl:tref result  1 0)
        0)
  (setf (magicl:tref result  2 0)
        0)
  ;; Since off diagonal components get halved, then voigt doubles them this is net 1d0
  (setf (magicl:tref result 3 0)
        (- (the double-float (magicl:tref stretch 2 1))
                  (the double-float (magicl:tref stretch 1 2))))
  (setf (magicl:tref result 4 0)
        (- (the double-float (magicl:tref stretch 0 2))
                  (the double-float (magicl:tref stretch 2 0))))
  (setf (magicl:tref result 5 0)
        (- (the double-float (magicl:tref stretch 0 1))
                  (the double-float (magicl:tref stretch 1 0))))

  (values))

