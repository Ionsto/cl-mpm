(defpackage :cl-mpm
  (:use :cl))
(in-package :cl-mpm)

(defun det-int-force (mp node dsvp)
  "Calculate internal force contribution from mp at node"
  (with-accessors ( (stress cl-mpm/particle:mp-stress)
                    (volume cl-mpm/particle:mp-volume)) mp

    (magicl:@  (magicl:transpose dsvp) (magicl:scale stress volume))
    ))

(defun det-ext-force (mp node svp)
  "Calculate external force contribution from mp at node"
  (with-accessors ( (mass cl-mpm/particle:mp-mass)) mp
    (magicl:scale (magicl:from-list (list 0d0 (cl-mpm/particle:mp-gravity mp)) '(2 1)) (* mass svp))))
