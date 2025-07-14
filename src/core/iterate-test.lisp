
(defpackage :cl-mpm/iter-test
  (:use #:coalton #:coalton-prelude
        #:coalton-library/math))

(cl:in-package #:cl-mpm/iter-test)

(named-readtables:in-readtable coalton:coalton)

(cl:declaim (cl:optimize (cl:speed 3) (cl:safety 0) (cl:debug 0)))

(coalton-toplevel
  ;; (declare f (Single-Float -> Single-Float))
  ;; (define (f x)
  ;;   (- (* x x) (* 2.0f0 x)))
  (declare shape-gimp-grad (Double-Float -> Double-Float -> Double-Float))
  (define (shape-gimp-grad x h)
    (if (> x 0d0)
        (/ -1d0 h)
        (/ 1d0 h)))
  (declare shape-gimp (Double-Float -> Double-Float -> Double-Float -> Double-Float))
  (define (shape-gimp x l h)
    (cond
      ((and (< (negate (+ h l)) x) (<= x (- l h)))
       (/ (^ (+ (+ x h) l) 2) (* 4d0 (* h l)))
       )
      ((and (< (- l h) x) (<= x (negate l)))
       (+ 1d0 (/ x h)))
      ((and (< (negate l) x) (<= x l))
       (- 1d0 (/ (+ (* x x) (* l l)) (* 2d0 (* h l)))))
      ((and (< l x) (<= x (- h l)))
       (- 1d0 (/ x h)))
      ((and (< (- h l) x) (<= x (+ h l)))
       (/ (pow (- (+ h l) x) 2d0) (* 4d0 (* h l))))
      (True
       0d0)))

  ;; (define-struct Vector
  ;;   (array Vector)
  ;;   )
  (define-struct mp
    (pos (Vector Double-Float))
    (domain-size (Vector Double-Float)))

  (define-struct mesh
    (h Double-Float))

  (define (iterate-over-neighbours-shape-gimp-test mesh mp func)
    "Iterate over a gimp domains neighbours in 2D unrolled, this version is more performant but worse than SIMD"
    (progn
      ;; (with-accessors ((pos-vec cl-mpm/particle:mp-position)
      ;;                  (d0 cl-mpm/particle::mp-domain-size))
      ;;     mp)
      (let ((pos-vec (.pos mp))
            (d0 (.domain-size mp))
            )
        (let ((h (.h mesh)))
          (let* ((pa (cl-mpm/utils:fast-storage pos-vec))
                 (da (cl-mpm/utils:fast-storage d0))
                 (px (aref pa 0))
                 (py (aref pa 1))
                 (ix (the fixnum (truncate (round px h))))
                 (iy (the fixnum (truncate (round py h))))
                 (cx (- px (* h ix)))
                 (cy (- py (* h iy)))
                 (dox (* 0.5d0 (the double-float (aref da 0))))
                 (doy (* 0.5d0 (the double-float (aref da 1))))
                 ;; (dxf (the fixnum (truncate (ffloor   (- cx dox) h))))
                 ;; (dxc (the fixnum (truncate (fceiling (+ cx dox) h))))
                 ;; (dyf (the fixnum (truncate (ffloor   (- cy doy) h))))
                 ;; (dyc (the fixnum (truncate (fceiling (+ cy doy) h))))
                 )
            ;; (loop for dx fixnum from -2 to 2;dxf to dxc
            ;;       do (loop for dy fixnum from -2 to 2;dyf to dyc
            ;;                do
            ;;                   (let* (;; (id (list (+ ix dx)
            ;;                          ;;           (+ iy dy)
            ;;                          ;;           0))
            ;;                          )
            ;;                     ;; (declare (dynamic-extent id))
            ;;                     (when t;(cl-mpm/mesh:in-bounds mesh id)
            ;;                       (let* ((distx (- cx (* h dx)))
            ;;                              (disty (- cy (* h dy)))
            ;;                              (weightsx (cl-mpm/shape-function::shape-gimp distx dox h))
            ;;                              (weightsy (cl-mpm/shape-function::shape-gimp disty doy h))
            ;;                              (weights-fbar-x (the double-float (cl-mpm/shape-function::shape-gimp-fbar distx dox h)))
            ;;                              (weights-fbar-y (the double-float (cl-mpm/shape-function::shape-gimp-fbar disty doy h)))
            ;;                              (weight-fbar (* weights-fbar-x weights-fbar-y))
            ;;                              (weight (* weightsx weightsy)))
            ;;                         (declare ;(type double-float h)
            ;;                          (double-float weight weightsx weightsy distx disty))
            ;;                         (when (or (> weight 0d0) (> weight-fbar 0d0))
            ;;                           (let* ((node
            ;;                                    nil
            ;;                                    ;; (cl-mpm/mesh:get-node mesh id)
            ;;                                        )
            ;;                                  (gradx
            ;;                                    (* (cl-mpm/shape-function::shape-gimp-dsvp distx dox h)
            ;;                                            weightsy))
            ;;                                  (grady
            ;;                                    (* (cl-mpm/shape-function::shape-gimp-dsvp disty doy h)
            ;;                                            weightsx))
            ;;                                  (grads-fbar
            ;;                                    nil
            ;;                                    ;; (list (* weights-fbar-y (cl-mpm/shape-function::shape-gimp-dsvp distx dox h))
            ;;                                    ;;       (* weights-fbar-x (cl-mpm/shape-function::shape-gimp-dsvp disty doy h))
            ;;                                    ;;       0d0)
            ;;                                    )
            ;;                                  )
            ;;                             (declare (double-float gradx grady))
            ;;                             (funcall func mesh mp node
            ;;                                      weight nil
            ;;                                      weight-fbar
            ;;                                      grads-fbar))))))))
            )))))

  )
