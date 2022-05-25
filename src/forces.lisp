(defpackage :cl-mpm
  (:use :cl))
(in-package :cl-mpm)

(defun det-int-force (mp node dsvp)
  "Calculate internal force contribution from mp at node"
  (with-slots ( (stress stress)
                (volume volume)) mp

    (magicl:@  (magicl:transpose dsvp) (magicl:scale stress volume))
    ))

(defun det-ext-force (mp node svp)
  "Calculate external force contribution from mp at node"
  (with-slots ( (mass mass)) mp
    (magicl:scale (magicl:from-list '(0d0 -9.8d0) '(2 1)) (* mass svp))))
