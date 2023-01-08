(defpackage :cl-mpm
  (:use :cl))
(in-package :cl-mpm)
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(declaim
 (inline det-int-force)
 (ftype (function (cl-mpm/particle::particle cl-mpm/mesh::node magicl:matrix/double-float) magicl:matrix/double-float)
        det-int-force))
(defun det-int-force (mp node dsvp)
  "Calculate internal force contribution from mp at node"
  (with-accessors ( (stress cl-mpm/particle:mp-stress)
                    (volume cl-mpm/particle:mp-volume)) mp
    (magicl:@ (magicl:transpose dsvp) (magicl:scale stress volume))
    ))

(declaim
 (inline det-ext-force)
 (ftype (function (cl-mpm/particle::particle cl-mpm/mesh::node double-float) magicl:matrix/double-float)
        det-ext-force))
(defun det-ext-force (mp node svp)
  "Calculate external force contribution from mp at node"
  (with-accessors ( (mass cl-mpm/particle:mp-mass)) mp
    (magicl:scale
     (magicl:.+
      ;; (magicl:from-array (make-array 2 :initial-contents (list 0d0 (* mass (cl-mpm/particle:mp-gravity mp)))) '(2 1) :type 'double-float)
      (magicl:from-list (list 0d0 (* mass (cl-mpm/particle:mp-gravity mp))) '(2 1) :type 'double-float)
      (magicl:scale (cl-mpm/particle:mp-body-force mp) mass)
      ) svp)))
