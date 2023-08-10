(defpackage :cl-mpm/utils
  (:use :cl)
  (:export
   #:voigt-zeros
   #:matrix-zeros
   #:stretch-dsvp-zeros
   #:voigt-from-list
   #:matrix-from-list
   #:matrix-to-voight
   #:voight-to-matrix
   #:voight-to-stretch
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
   ))
(in-package :cl-mpm/utils)
(declaim (optimize (debug 3) (safety 0) (speed 3)))

(declaim (inline eig)
         (ftype (function (magicl:matrix/double-float)
                          (values list magicl:matrix/double-float)) eig))
(defun eig (mat)
  (magicl:self-adjoint-eig mat)
  ;; (multiple-value-bind (l v) (magicl:eig mat)
  ;;   (values l (magicl:.realpart v)))
  )

(declaim (inline vector-zeros)
         (ftype (function ()
                          magicl:matrix/double-float) vector-zeros))
(defun vector-zeros ()
  (magicl::make-matrix/double-float 2 1 2 :column-major (make-array 2 :element-type 'double-float)))

(declaim (inline voigt-zeros)
         (ftype (function ()
                          magicl:matrix/double-float) voigt-zeros))
(defun voigt-zeros ()
  (magicl::make-matrix/double-float 3 1 3 :column-major (make-array 3 :element-type 'double-float)))

(declaim (inline matrix-zeros)
         (ftype (function ()
                          magicl:matrix/double-float) matrix-zeros))
(defun matrix-zeros ()
  (magicl::make-matrix/double-float 2 2 4 :column-major (make-array 4 :element-type 'double-float)))

(declaim (inline matrix-copy)
         (ftype (function (magicl:matrix/double-float)
                          magicl:matrix/double-float) matrix-copy))
(defun matrix-copy (mat)
  ;; (let ((m (magicl::make-matrix/double-float 2 2 4 :column-major (make-array 4 :element-type 'double-float))))
  ;;   (magicl::copy-matrix/double-float m mat)
  ;;   m)
  (let ((v (matrix-zeros)))
    (aops:copy-into (magicl::matrix/double-float-storage v)
                    (magicl::matrix/double-float-storage mat))
    v)
  ;; (magicl::copy-matrix/double-float mat)
  )
(defun vector-copy (vec)
  (let ((v (vector-zeros)))
    (aops:copy-into (magicl::matrix/double-float-storage v)
                    (magicl::matrix/double-float-storage vec))
    v))

(declaim (inline stretch-dsvp-zeros)
         (ftype (function ()
                          magicl:matrix/double-float) stretch-dsvp-zeros))
(defun stretch-dsvp-zeros ()
  (magicl::make-matrix/double-float 4 2 8 :column-major (make-array 8 :element-type 'double-float)))

(declaim (inline dsvp-2d-zeros)
         (ftype (function ()
                          magicl:matrix/double-float) dsvp-2d-zeros))
(defun dsvp-2d-zeros ()
  (magicl::make-matrix/double-float 3 2 6 :column-major (make-array 6 :element-type 'double-float)))

(declaim (inline voigt-from-list)
         (ftype (function (list)
                          magicl:matrix/double-float) voigt-from-list))
(defun voigt-from-list (elements)
  (magicl::make-matrix/double-float 3 1 3 :column-major
                                    (make-array 3 :element-type 'double-float :initial-contents elements)))

(declaim (inline matrix-from-list)
         (ftype (function (list)
                          magicl:matrix/double-float) matrix-from-list))
(defun matrix-from-list (elements)
  (magicl::make-matrix/double-float 2 2 4 :column-major
                                    (make-array 4 :element-type 'double-float :initial-contents elements)))


(declaim (inline matrix-to-voight)
         (ftype (function (magicl:matrix/double-float)
                          magicl:matrix/double-float) matrix-to-voight))
(defun matrix-to-voight (matrix)
  "Stress matrix to voigt"
  (let* ( (exx (magicl:tref matrix 0 0))
          (eyy (magicl:tref matrix 1 1))
          (exy (magicl:tref matrix 1 0)))
    (voigt-from-list (list exx eyy exy))))
(declaim (inline voight-to-matrix)
         (ftype (function (magicl:matrix/double-float)
                          magicl:matrix/double-float) voight-to-matrix))
(defun voight-to-matrix (vec)
  "Stress format voight to matrix"
  (let* ( (exx (magicl:tref vec 0 0))
          (eyy (magicl:tref vec 1 0))
          (exy (magicl:tref vec 2 0)))
    (matrix-from-list (list exx exy exy eyy))
    ))

(declaim (inline voigt-to-matrix)
         (ftype (function (magicl:matrix/double-float)
                          magicl:matrix/double-float) voigt-to-matrix))
(defun voigt-to-matrix (vec)
  (let* ( (exx (magicl:tref vec 0 0))
          (eyy (magicl:tref vec 1 0))
          (exy (* 0.5d0 (the double-float (magicl:tref vec 2 0)))))
    (matrix-from-list (list exx exy exy eyy))
    ))

(declaim (inline matrix-to-voigt)
         (ftype (function (magicl:matrix/double-float)
                          magicl:matrix/double-float) matrix-to-voigt))
(defun matrix-to-voigt (matrix)
  (let* ( (exx (magicl:tref matrix 0 0))
          (eyy (magicl:tref matrix 1 1))
          (exy (magicl:tref matrix 1 0)))
    (voigt-from-list (list exx eyy (* 2d0 (the double-float exy))))))

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
  (let* ( (exx (magicl:tref vec 0 0))
          (eyy (magicl:tref vec 1 0))
          (exy (magicl:tref vec 2 0))
          (eyx (magicl:tref vec 3 0))
          )
    (magicl:from-list (list exx exy
                            eyx eyy)
                      '(2 2) :type 'double-float)))

(declaim
 (inline voight-to-stretch-prealloc)
 (ftype (function (magicl:matrix/double-float magicl:matrix/double-float)
                  magicl:matrix/double-float) voight-to-stretch-prealloc))
(defun voight-to-stretch-prealloc (vec result)

  (let* ((exx (magicl:tref vec 0 0))
         (eyy (magicl:tref vec 1 0))
         (exy (magicl:tref vec 2 0))
         (eyx (magicl:tref vec 3 0))
         ;(s (magicl::matrix/double-float-storage result))
         (s result)
         )
    (declare (double-float exx eyy exy eyx))
    (setf (magicl:tref s 0 0) exx
          (magicl:tref s 0 1) exy
          (magicl:tref s 1 0) eyx
          (magicl:tref s 1 1) eyy
          ))
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
  (magicl:from-list list '(3 1) :type 'double-float))

(declaim (inline trace-voigt)
         (ftype (function (magicl:matrix/double-float) double-float) trace-voigt))
(defun trace-voigt (a)
  "Calculate the product A_{ij}A_{ij}"
  (let ((arr (magicl::matrix/double-float-storage a)))
    (declare ((simple-array double-float) arr))
    (+ (aref arr 0) (aref arr 1))))

(declaim (inline deviatoric-voigt)
         (ftype (function (magicl:matrix/double-float)
                          magicl:matrix/double-float) deviatoric-voigt))
(defun deviatoric-voigt (a)
  "Calculate the product A_{ij}A_{ij}"
  (let* ((tr (/ (trace-voigt a) 2d0))
         (arr (magicl::matrix/double-float-storage a)))
    (declare ((simple-array double-float) arr))
    (stress-from-list (list (- (aref arr 0) tr)
                            (- (aref arr 1) tr)
                            (aref arr 2)))))
