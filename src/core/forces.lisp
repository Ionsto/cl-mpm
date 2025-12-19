(defpackage :cl-mpm/forces
  (:use :cl
        :cl-mpm/utils)
  (:export
   #:det-int-force
   #:det-ext-force
   #:det-ext-force-2d
   #:det-int-force-unrolled
   #:det-int-force-unrolled-2d
   ))
(in-package :cl-mpm/forces)
(declaim (optimize (debug 0) (safety 0) (speed 3)))
;; (declaim #.cl-mpm/settings:*optimise-setting*)
;; (declaim (optimize (debug 3) (safety 3) (speed 0)))

(declaim
 (inline @-dsvp-vec)
 (ftype (function (magicl:matrix/double-float magicl:matrix/double-float double-float magicl:matrix/double-float)
                  (values)) @-dsvp-vec))
(defun @-dsvp-vec (matrix vector scale result-vector)
  "Multiply a 3x9 matrix with a 3x1 vector to calculate a 3x1 vector in place"
  (declare (magicl:matrix/double-float matrix vector result-vector)
           (double-float scale)
           (optimize (speed 3) (safety 0) (debug 0)))
  (let ((a (magicl::matrix/double-float-storage matrix))
        (b (magicl::matrix/double-float-storage vector))
        (c (magicl::matrix/double-float-storage result-vector))
        )
    (declare ((simple-array double-float (18)) a)
             ((simple-array double-float (3)) c)
             ((simple-array double-float (6)) b)
             )
    (flet ((tref (m x y)
             (aref m (+ (* 6 x) y))))
      (loop for i fixnum from 0 below 3
            do
               ;; (setf (aref c i) 0d0)
               (loop for j fixnum from 0 below 6
                     do (decf (aref c i) (the double-float (* (aref b j) (tref a i j) scale))))
            )))
  (values))

(defun @-dsvp-vec-simd (matrix vector scale result-vector)
  "Multiply a 3x6 matrix with a 6x1 voigt-vector to calculate a 3x1 vector in place"
  (declare (magicl:matrix/double-float matrix vector result-vector)
           (double-float scale)
           (optimize (speed 3) (safety 0) (debug 0)))
  (let ((a (magicl::matrix/double-float-storage matrix))
        (b (magicl::matrix/double-float-storage vector))
        (c (magicl::matrix/double-float-storage result-vector)))
    (declare ((simple-array double-float (18)) a)
             ((simple-array double-float (3)) c)
             ((simple-array double-float (6)) b)
             )
    (declare (type sb-simd:f64vec a b c))
    (macrolet
        ((result-component (i)
           (declare (fixnum i))
           `(decf
             (aref c ,i)
             (* scale
                (+
                 (sb-simd-avx:f64.4-horizontal+
                  (sb-simd-avx:f64.4*
                   (sb-simd-avx:f64.4-aref a ,(the fixnum (* i 6)))
                   (sb-simd-avx:f64.4-aref b 0)))
                 (sb-simd-avx:f64.2-horizontal+
                  (sb-simd-avx:f64.2*
                   (sb-simd-avx:f64.2-aref a ,(the fixnum (+ (* i 6) 4)))
                   (sb-simd-avx:f64.2-aref b 4))))))))
      (result-component 0)
      (result-component 1)
      (result-component 2)))
  (values))


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
    (loop for i from 0 below 3
          do (loop for j from 0 below 6
                   do (decf (the double-float (magicl:tref res i 0))
                            (* (the double-float (magicl:tref a j i))
                               (the double-float (magicl:tref b j 0)) scale))))))

(defun test-@-dsvp-vec-simd ()
  (let ((volume 0.1d0)
        (dsvp (cl-mpm/shape-function::assemble-dsvp-3d (list 6d0 2d0 1.5d0)))
        (stress (cl-mpm/utils:voigt-from-list (list 9d0 3d0 5d0 3d0 5d0 8d0))))
                                        ;(mult-force dsvp stress volume f-out)
    ;; (cl-mpm/fastmaths::@-dsvp-vec dsvp stress volume f-out)
    ;; (@-dsvp-vec-simd dsvp stress volume f-out)

    (let ((res (cl-mpm/utils:vector-zeros)))
      (magicl:.- res (magicl:scale! (magicl:@ (magicl:transpose dsvp) stress) volume) res)
      (pprint res))
    (let ((res (cl-mpm/utils:vector-zeros)))
      (@-dsvp-vec-simd dsvp stress volume res)
      ;; (magicl:.- res (magicl:scale! (magicl:@ (magicl:transpose dsvp) stress) volume) res)
      (pprint res))
    (let ((iters 1000000))
      (let ((res (cl-mpm/utils:vector-zeros)))
        (time 
         (dotimes (i iters)
           (magicl:.- res (magicl:scale! (magicl:@ (magicl:transpose dsvp) stress) volume) res)))
        )
      (let ((res (cl-mpm/utils:vector-zeros)))
        (time 
         (dotimes (i iters)
           (@-dsvp-vec-simd dsvp stress volume res)))
        ))
    )
  )

(declaim
 (inline mult-force-plane-strain)
 (ftype (function (magicl:matrix/double-float
                   magicl:matrix/double-float
                   double-float
                   magicl:matrix/double-float) (values)
                  ) mult-force-plane-strain))



(declaim
 (inline det-int-force)
 (ftype (function (cl-mpm/particle::particle magicl:matrix/double-float
                                             &optional magicl:matrix/double-float) magicl:matrix/double-float)
        det-int-force))
(defun det-int-force (mp dsvp volume &optional f-out)
  "Calculate internal force contribution from mp at node"
  (let* ((f-out (if f-out f-out (cl-mpm/utils:vector-zeros))))
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     ;; (volume cl-mpm/particle:mp-volume)
                     (vel cl-mpm/particle:mp-velocity)
                     ) mp
      (declare (type double-float volume))
      (@-dsvp-vec-simd dsvp stress volume f-out))
    f-out))

(defun det-stress-force-unrolled (stress grads volume &optional f-out)
  "Calculate internal force contribution from mp at node"
  (let* ((f-out (if f-out f-out (cl-mpm/utils:vector-zeros))))
    (declare (type double-float volume))
     (destructuring-bind (dx dy dz) grads
       (declare (double-float dx dy dz))
       (let ((f-storage (cl-mpm/utils:fast-storage f-out))
             (stress-storage (cl-mpm/utils:fast-storage stress)))
         (incf (aref f-storage 0)
               (* -1d0 volume
                  (+
                   (* dx (aref stress-storage 0))
                   (* dz (aref stress-storage 4))
                   (* dy (aref stress-storage 5)))))
         (incf (aref f-storage 1)
               (* -1d0 volume
                  (+
                   (* dy (aref stress-storage 1))
                   (* dz (aref stress-storage 3))
                   (* dx (aref stress-storage 5)))))
         (incf (aref f-storage 2)
               (* -1d0 volume
                  (+
                   (* dz (aref stress-storage 2))
                   (* dy (aref stress-storage 3))
                   (* dx (aref stress-storage 4)))))))
    f-out))

(declaim
 (inline det-int-force-unrolled)
 (ftype (function (cl-mpm/particle::particle list double-float &optional magicl:matrix/double-float) magicl:matrix/double-float)
        det-int-force-unrolled))
(defun det-int-force-unrolled (mp grads volume &optional f-out)
  "Calculate internal force contribution from mp at node"
  (let* ((f-out (if f-out f-out (cl-mpm/utils:vector-zeros))))
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     ;; (volume-ac cl-mpm/particle:mp-volume)
                     )
        mp
      (let ()
        (declare (type double-float volume))
        (destructuring-bind (dx dy dz) grads
          (declare (double-float dx dy dz))
          (let ((f-storage (cl-mpm/utils:fast-storage f-out))
                (stress-storage (cl-mpm/utils:fast-storage stress)))
            (incf (aref f-storage 0)
                  (* -1d0 volume
                     (+
                      (* dx (aref stress-storage 0))
                      (* dz (aref stress-storage 4))
                      (* dy (aref stress-storage 5)))))
            (incf (aref f-storage 1)
                  (* -1d0 volume
                     (+
                      (* dy (aref stress-storage 1))
                      (* dz (aref stress-storage 3))
                      (* dx (aref stress-storage 5)))))
            (incf (aref f-storage 2)
                  (* -1d0 volume
                     (+
                      (* dz (aref stress-storage 2))
                      (* dy (aref stress-storage 3))
                      (* dx (aref stress-storage 4)))))))))
    f-out))
(declaim
 (inline det-int-force-unrolled-2d)
 (ftype (function (cl-mpm/particle::particle list double-float &optional magicl:matrix/double-float) magicl:matrix/double-float)
        det-int-force-unrolled-2d))
(defun det-int-force-unrolled-2d (mp grads volume &optional f-out)
  "Calculate internal force contribution from mp at node"
  (let* ((f-out (if f-out f-out (cl-mpm/utils:vector-zeros))))
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     ;; (volume-ac cl-mpm/particle:mp-volume)
                     ) mp
      (let (;; (volume volume-ac)
            )
        (declare (type double-float volume))
        (destructuring-bind (dx dy dz) grads
          (declare (double-float dx dy dz))
          (let ((f-storage (cl-mpm/utils:fast-storage f-out))
                (stress-storage (cl-mpm/utils:fast-storage stress)))
            (incf (aref f-storage 0)
                  (* -1d0 volume
                     (+
                      (* dx (aref stress-storage 0))
                      (* dz (aref stress-storage 4))
                      (* dy (aref stress-storage 5)))))
            (incf (aref f-storage 1)
                  (* -1d0 volume
                     (+
                      (* dy (aref stress-storage 1))
                      (* dz (aref stress-storage 3))
                      (* dx (aref stress-storage 5)))))))))
    f-out))


(declaim
 (inline det-ext-force)
 (ftype (function (cl-mpm/particle::particle cl-mpm/mesh::node double-float double-float double-float &optional magicl:matrix/double-float) magicl:matrix/double-float)
        det-ext-force))
(defun det-ext-force (mp node svp gravity volume &optional f-out)
  "Calculate external force contribution from mp at node"
  (with-accessors ((mass cl-mpm/particle:mp-mass)
                   ;; (gravity cl-mpm/particle:mp-gravity)
                   ;; (volume cl-mpm/particle:mp-volume)
                   (body-force cl-mpm/particle:mp-body-force)
                   (gravity-axis cl-mpm/particle::mp-gravity-axis)
                   ) mp
    (declare (type double-float svp mass gravity volume)
             (type magicl:matrix/double-float body-force))
    (let* ((f-out (if f-out f-out (cl-mpm/utils::vector-zeros)))
           (f-s (fast-storage f-out))
           (b-s (fast-storage body-force))
           (g-s (fast-storage gravity-axis))
           )
      (declare (type (simple-array double-float (3)) f-s b-s g-s))
      ;;Manually unrolled
      (setf
       (sb-simd-avx:f64.2-aref f-s 0)
       (sb-simd-avx:f64.2+
        (sb-simd-avx:f64.2-aref f-s 0)
        (sb-simd-avx:f64.2*
         (sb-simd-avx:f64.2+
          (sb-simd-avx:f64.2*
           (sb-simd-avx:f64.2-aref g-s 0)
           (* mass gravity))
          (sb-simd-avx:f64.2*
           (sb-simd-avx:f64.2-aref b-s 0)
           volume))
         svp)))
      (incf (aref f-s 2)
            (*
             (+
              (* mass gravity (aref g-s 2))
              (* volume (aref b-s 2)))
             svp))
      f-out)))

(declaim
 (inline det-ext-force-2d)
 (ftype (function (cl-mpm/particle::particle cl-mpm/mesh::node double-float double-float double-float &optional magicl:matrix/double-float) magicl:matrix/double-float)
        det-ext-force-2d))
(defun det-ext-force-2d (mp node svp gravity volume &optional f-out)
  (declare (double-float svp gravity))
  "Calculate external force contribution from mp at node"
  (with-accessors ((mass cl-mpm/particle:mp-mass)
                   ;; (gravity cl-mpm/particle:mp-gravity)
                   ;; (volume cl-mpm/particle:mp-volume)
                   (body-force cl-mpm/particle:mp-body-force)
                   (gravity-axis cl-mpm/particle::mp-gravity-axis)
                   ) mp
    (declare (type double-float svp mass gravity volume)
             (type magicl:matrix/double-float body-force))
    (let* ((f-out (if f-out f-out (cl-mpm/utils::vector-zeros)))
           (f-s (fast-storage f-out))
           (b-s (fast-storage body-force))
           (g-s (fast-storage gravity-axis))
           )
      (declare (type (simple-array double-float (3)) f-s b-s g-s))
      ;;Manually unrolled
      (setf
       (sb-simd-avx:f64.2-aref f-s 0)
       (sb-simd-avx:f64.2+
        (sb-simd-avx:f64.2-aref f-s 0)
        (sb-simd-avx:f64.2*
         (sb-simd-avx:f64.2+
          (sb-simd-avx:f64.2*
           (sb-simd-avx:f64.2-aref g-s 0)
           (* mass gravity))
          (sb-simd-avx:f64.2*
           (sb-simd-avx:f64.2-aref b-s 0)
           volume))
         svp)))
      f-out)))
