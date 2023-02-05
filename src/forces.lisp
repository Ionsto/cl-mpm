(defpackage :cl-mpm
  (:use :cl))
(in-package :cl-mpm)
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(declaim
 (inline mult-force)
 (ftype (function (magicl:matrix/double-float
                   magicl:matrix/double-float
                   double-float
                   magicl:matrix/double-float) (values)
                  ) mult-force))
(defun mult-force (a b scale res)
  (declare (type magicl:matrix/double-float a b res)
           (type double-float scale))
  (let ((a-s (magicl::storage a))
        (b-s (magicl::storage b))
        (res-s (magicl::storage res))
        )
    (declare (type (simple-array double-float) a-s b-s res-s))
    (loop for i from 0 to 1
          do (loop for j from 0 to 2
                   do (decf (aref res-s i) (* (aref a-s (+ i (* 2 j)))
                                              (aref b-s j) scale))))
  ))
(declaim
 (inline det-int-force)
 (ftype (function (cl-mpm/particle::particle magicl:matrix/double-float
                                             &optional magicl:matrix/double-float) magicl:matrix/double-float)
        det-int-force))
(defun det-int-force (mp dsvp &optional f-out)
  "Calculate internal force contribution from mp at node"
  (let* ((f-out (if f-out f-out (magicl:zeros '(2 1)))))
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (volume cl-mpm/particle:mp-volume)) mp
      (declare (type double-float volume))
      ;; (mult-force dsvp stress volume f-out))
      ;; (magicl:.- f-out (magicl:scale! (magicl:@ (magicl:transpose dsvp) stress) volume) f-out))
      (magicl:.- f-out (magicl:scale! (magicl:@ (magicl:transpose dsvp) stress) volume) f-out))
    f-out))

(declaim
 (inline det-ext-force)
 (ftype (function (cl-mpm/particle::particle cl-mpm/mesh::node double-float &optional magicl:matrix/double-float) magicl:matrix/double-float)
        det-ext-force))
(defun det-ext-force (mp node svp &optional f-out)
  "Calculate external force contribution from mp at node"
  (with-accessors ((mass cl-mpm/particle:mp-mass)
                   (gravity cl-mpm/particle:mp-gravity)
                   (volume cl-mpm/particle:mp-volume)
                   (body-force cl-mpm/particle:mp-body-force)
                   (gravity-axis cl-mpm/particle::mp-gravity-axis)
                   ) mp
    (declare (type double-float mass gravity)
             (type magicl:matrix/double-float body-force))
    (let* ((f-out (if f-out f-out (magicl:zeros '(2 1))))
           (f-s (magicl::storage f-out))
           (b-s (magicl::storage body-force))
           ;; (g-s (magicl::storage (magicl:scale gravity-axis gravity)))
           )
      (declare (type (simple-array double-float *) f-s b-s))
      ;;Manually unrolled

      (incf (aref f-s 0)
            (* (aref b-s 0) volume svp))
      (incf (aref f-s 1)
            (* (+ (* mass gravity) (* volume (aref b-s 1))) svp)
            )

          ;; (magicl:scale!
          ;;   ;; (magicl:.+ (magicl:from-array (make-array 2 :initial-contents (list 0d0 (* mass gravity)))
          ;;   ;;                               '(2 1) :type 'double-float :layout :column-major)
          ;;   (magicl:.+ (magicl:from-list (list 0d0 (* mass gravity)) '(2 1) :type 'double-float)
          ;;              (magicl:scale body-force mass))
          ;;   svp))
      f-out)))
