(defpackage :cl-mpm/utils
  (:use :cl)
  (:export
   #:voigt-zeros
   #:matrix-zeros
   #:stretch-dsvp-zeros
   #:stretch-dsvp-3d-zeros
   #:voigt-from-list
   #:matrix-from-list
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
   #:get-stress
   ))
(in-package :cl-mpm/utils)
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(declaim (inline eig)
         (ftype (function (magicl:matrix/double-float)
                          (values list magicl:matrix/double-float)) eig))
(defun eig (mat)
  (magicl:self-adjoint-eig mat)
  ;; (multiple-value-bind (l v) (magicl:eig mat)
  ;;   (values l (magicl:.realpart v)))
  )
(let ((stress-components '(:xx :yy :zz :yz :xz :xy)))
  (defun get-stress (stress component)
    "Get stress component from voigt notation stress vector"
    (aref (magicl::matrix/double-float-storage stress) (position component stress-components))))

(declaim (inline vector-zeros)
         (ftype (function ()
                          magicl:matrix/double-float) vector-zeros))
(defun vector-zeros ()
  (magicl::make-matrix/double-float 3 1 3 :column-major (make-array 3 :element-type 'double-float)))

(declaim (inline voigt-zeros)
         (ftype (function ()
                          magicl:matrix/double-float) voigt-zeros))
(defun voigt-zeros ()
  (magicl::make-matrix/double-float 6 1 6 :column-major (make-array 6 :element-type 'double-float)))

(declaim (inline matrix-zeros)
         (ftype (function ()
                          magicl:matrix/double-float) matrix-zeros))
(defun matrix-zeros ()
  (magicl::make-matrix/double-float 3 3 9 :column-major (make-array 9 :element-type 'double-float)))

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

(declaim (inline stretch-dsvp-zeros)
         (ftype (function ()
                          magicl:matrix/double-float) stretch-dsvp-zeros))
(defun stretch-dsvp-zeros ()
  (magicl::make-matrix/double-float 4 2 8 :column-major (make-array 8 :element-type 'double-float)))

(defun stretch-dsvp-3d-zeros ()
  (magicl::make-matrix/double-float 9 3 27 :column-major (make-array 27 :element-type 'double-float)))


(declaim (inline dsvp-2d-zeros)
         (ftype (function ()
                          magicl:matrix/double-float) dsvp-2d-zeros))
(defun dsvp-2d-zeros ()
  (magicl::make-matrix/double-float 3 2 6 :column-major (make-array 6 :element-type 'double-float)))

(defun dsvp-3d-zeros ()
  (magicl::make-matrix/double-float 6 2 12 :column-major (make-array 12 :element-type 'double-float)))

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
    (matrix-from-list (list exx exy eyz
                            exy eyy ezx
                            eyz ezx ezz))))

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

(defun stretch-to-voight (matrix)
  (let* ( (exx (magicl:tref matrix 0 0))
          (eyy (magicl:tref matrix 1 1))
          (exy (magicl:tref matrix 1 0))
          (eyx (magicl:tref matrix 1 0))
          )
    (magicl:from-list (list exx eyy
                            exy eyx)
                      '(4 1) :type 'double-float)))

(defun voight-to-stretch-3d (vec)
  (let* ((exx (magicl:tref vec 0 0))
         (eyy (magicl:tref vec 1 0))
         (ezz (magicl:tref vec 2 0))
         (exy (magicl:tref vec 3 0))
         (eyx (magicl:tref vec 4 0))
         (exz (magicl:tref vec 5 0))
         (ezx (magicl:tref vec 6 0))
         (ezy (magicl:tref vec 7 0))
         (eyz (magicl:tref vec 8 0))
         )
    (magicl:from-list (list exx exy exz
                            eyx eyy eyz
                            ezx ezy ezz)
                      '(3 3) :type 'double-float)))

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
           (ezy (aref vecs 7))
           (eyz (aref vecs 8))
           (s result)
           )
      (declare (double-float exx eyy exy eyx))
      (setf (magicl:tref s 0 0) exx
            (magicl:tref s 0 1) exy
            (magicl:tref s 1 0) eyx
            (magicl:tref s 1 1) eyy
            (magicl:tref s 0 2) exz
            (magicl:tref s 2 0) ezx
            (magicl:tref s 1 2) eyz
            (magicl:tref s 2 1) ezy
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
