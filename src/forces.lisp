(defpackage :cl-mpm
  (:use :cl))
(in-package :cl-mpm)
(declaim (optimize (debug 3) (safety 3) (speed 3)))

(defun det-int-force (mp node dsvp)
  "Calculate internal force contribution from mp at node"
  (with-accessors ( (stress cl-mpm/particle:mp-stress)
                    (volume cl-mpm/particle:mp-volume)) mp

    (magicl:@  (magicl:transpose dsvp) (magicl:scale stress volume))
    ))

(defun det-ext-force (mp node svp)
  "Calculate external force contribution from mp at node"
  (with-accessors ( (mass cl-mpm/particle:mp-mass)) mp
    (magicl:scale
     (magicl:.+
      (magicl:from-list (list 0d0 (* mass (cl-mpm/particle:mp-gravity mp))) '(2 1))
      (magicl:scale (cl-mpm/particle:mp-body-force mp) mass)
      ) svp)))
