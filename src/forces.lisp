(defpackage :cl-mpm/forces
  (:use :cl)
  (:export
   #:det-int-force
   #:det-ext-force
   ))
(in-package :cl-mpm/forces)
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
  (declare (optimize (safety 3)))
  (let (
        ;; (a-s (magicl::matrix/double-float-storage a))
        ;; (b-s (magicl::matrix/double-float-storage b))
        ;; (res-s (magicl::matrix/double-float-storage res))
        )
    ;; (declare (type (simple-array double-float) a-s b-s res-s))
    (loop for i from 0 to 1
          do (loop for j from 0 to 3
                   do (decf (the double-float (magicl:tref res i 0)) (* (the double-float (magicl:tref a j i))
                                                     (the double-float (magicl:tref b j 0)) scale))))))

(declaim
 (inline mult-force-plane-strain)
 (ftype (function (magicl:matrix/double-float
                   magicl:matrix/double-float
                   double-float
                   magicl:matrix/double-float) (values)
                  ) mult-force-plane-strain))
(defun mult-force-plane-strain (a b scale res)
  (declare (type magicl:matrix/double-float a b res)
           (type double-float scale))
  (declare (optimize (speed 3) (safety 0)))
  (let ((rs (magicl::matrix/double-float-storage res))
        (bs (magicl::matrix/double-float-storage b))
        )
    ;; (declare (type (simple-array double-float) a-s b-s res-s))
    (loop for i from 0 to 1
          do
             ;; (decf (the double-float (aref rs i)) (* (the double-float (magicl:tref a 0 i))
             ;;                                                   (the double-float (magicl:tref b 0 0)) scale))
             ;; (decf (the double-float (aref rs i)) (* (the double-float (magicl:tref a 1 i))
             ;;                                                   (the double-float (magicl:tref b 1 0)) scale))
             ;; (decf (the double-float (aref rs i)) (* (the double-float (magicl:tref a 2 i))
             ;;                                                   (the double-float (magicl:tref b 2 0)) scale))

             ;; (loop for j in '(0 1 5)
             ;;       do (decf (the double-float (magicl:tref res i 0)) (* (the double-float (magicl:tref a j i))
             ;;                                                            (the double-float (magicl:tref b j 0)) scale)))
             (decf (the double-float (aref rs i)) (* (the double-float (magicl:tref a 0 i))
                                                     (the double-float (aref bs 0)) scale))
             (decf (the double-float (aref rs i)) (* (the double-float (magicl:tref a 1 i))
                                                     (the double-float (aref bs 1)) scale))
             (decf (the double-float (aref rs i)) (* (the double-float (magicl:tref a 5 i))
                                                     (the double-float (aref bs 5)) scale))
             ;; (loop for j in '(0 1 5)
             ;;       do (decf (the double-float (aref rs i)) (* (the double-float (magicl:tref a j i))
             ;;                                                  (the double-float (aref bs j)) scale)))

          )))

(defun plane-strain-transform (stress)
  (magicl:from-list (list (magicl:tref stress 0 0)
                          (magicl:tref stress 1 0)
                          (magicl:tref stress 5 0))
                    '(3 1)
                    :type 'double-float))

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
      ;; (print dsvp)
      ;(mult-force dsvp stress volume f-out)
      ;; (mult-force dsvp (plane-strain-transform stress) volume f-out)
      (mult-force-plane-strain dsvp stress volume f-out)
      ;; (magicl:.- f-out (magicl:scale! (magicl:@ (magicl:transpose dsvp) (plane-strain-transform stress)) volume) f-out)
      ;; (magicl:.- f-out (magicl:scale! (magicl:@ (magicl:transpose dsvp) stress) volume) f-out))
      )
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
    (declare (type double-float mass gravity volume)
             (type magicl:matrix/double-float body-force))
    (let* ((f-out (if f-out f-out (magicl:zeros '(2 1))))
           (f-s (magicl::matrix/double-float-storage f-out))
           (b-s (magicl::matrix/double-float-storage body-force))
           (g-s (magicl::matrix/double-float-storage gravity-axis))
           )
      (declare (type (simple-array double-float *) f-s b-s g-s))
      ;;Manually unrolled

      (incf (aref f-s 0)
            (* (+ (* mass gravity (aref g-s 0)) (* volume (aref b-s 0))) svp))
      (incf (aref f-s 1)
            (* (+ (* mass gravity (aref g-s 1)) (* volume (aref b-s 1))) svp)
            )

          ;; (magicl:scale!
          ;;   ;; (magicl.simd::.+-simd (magicl:from-array (make-array 2 :initial-contents (list 0d0 (* mass gravity)))
          ;;   ;;                               '(2 1) :type 'double-float :layout :column-major)
          ;;   (magicl.simd::.+-simd (magicl:from-list (list 0d0 (* mass gravity)) '(2 1) :type 'double-float)
          ;;              (magicl:scale body-force mass))
          ;;   svp))
      f-out)))
