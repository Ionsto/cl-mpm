(defpackage :cl-mpm/utils
  (:use :cl)
  (:export
   #:matrix-to-voight
   #:voight-to-matrix
   #:voight-to-stretch
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
