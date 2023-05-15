(defpackage :cl-mpm/utils
  (:use :cl)
  (:export
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
   ))
(in-package :cl-mpm/utils)
(declaim (optimize (debug 0) (safety 0) (speed 3)))
(defun matrix-to-voight (matrix)
  (let* ( (exx (magicl:tref matrix 0 0))
          (eyy (magicl:tref matrix 1 1))
          (exy (magicl:tref matrix 1 0)))
    (magicl:from-list (list exx eyy
                            exy)
                      '(3 1) :type 'double-float)))
(defun voight-to-matrix (vec)
  (let* ( (exx (magicl:tref vec 0 0))
          (eyy (magicl:tref vec 1 0))
          (exy (magicl:tref vec 2 0)))
    (magicl:from-list (list exx exy
                            exy eyy)
                      '(2 2) :type 'double-float)))

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
         ;(s (magicl::storage result))
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
  (let ((arr (magicl::storage a)))
    (declare ((simple-array double-float) arr))
    (+ (aref arr 0)
               (aref arr 1))))

(declaim (inline deviatoric-voigt)
         (ftype (function (magicl:matrix/double-float)
                          magicl:matrix/double-float) deviatoric-voigt))
(defun deviatoric-voigt (a)
  "Calculate the product A_{ij}A_{ij}"
  (let* ((tr (/ (trace-voigt a) 2d0))
         (arr (magicl::storage a)))
    (declare ((simple-array double-float) arr))
    (stress-from-list (list (- (aref arr 0) tr)
                            (- (aref arr 1) tr)
                            (aref arr 2)))))
