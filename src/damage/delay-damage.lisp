(in-package :cl-mpm/damage)
(declaim #.cl-mpm/settings:*optimise-setting*)



(defun deriv-partial (k y k0 tau n)
  (declare (double-float k y k0 tau n))
  (/
   (* k0
      (the double-float
           (expt
            (/
             (the double-float (max 0d0 (- y (max k k0))))
             ;; (the double-float (max 0d0 (- y k)))
             (max
              k0
              k)) n)))
   tau))

(defun huen-integration (k y-0 y-1 k0 tau n dt)
  (declare (double-float k y-0 y-1 k0 tau n dt))
  (let* ((dk-0 (deriv-partial k y-0 k0 tau n))
         (dk-1 (deriv-partial (+ k (* dt dk-0)) y-1 k0 tau n)))
    (declare (double-float dk-0 dk-1))
    (the double-float (+ k (* (/ dt 2) (+ dk-0 dk-1))))))

(defun forwards-integration (k y-0 y-1 k0 tau n dt)
  (declare (double-float k y-0 y-1 k0 tau n dt))
  (+ k (* dt (deriv-partial k y-0 k0 tau n))))


(defun integrate-substep (k y-0 y-1 dt iters function)
  (declare (double-float k y-0 y-1 dt)
           (function function)
           ;; ((integer 1 1000000) iters)
           )
  (let ((kprev k)
        (yprev y-0)
        (ycurrent y-0)
        (yinc (/ (- y-1 y-0) iters)))
    (loop for i from 0 below iters
          do
             (incf ycurrent yinc)
             (setf kprev
                   (funcall function
                            kprev
                            yprev
                            ycurrent
                            (/ dt iters)))
             (setf yprev ycurrent))
    kprev))

(defun analytic-trim (k y-0 y-1 k0 dt function)
  (if (and (< y-0 k0)
           (> y-1 y-0))
      (let* ((ratio (/ (- y-1 k0)
                       (- y-1 y-0))))
        (funcall function k k0 y-1 (* dt ratio)))
      (funcall function k y-0 y-1 dt)))

(defun auto-refine-substepper (k y-0 y-1 dt function &key (tol 1d-3))
  (declare (double-float k y-0 y-1 dt)
           (function function))
  (let* ((r 0)
         (kn0 (integrate-substep k y-0 y-1 dt (expt 2 r) function))
         (kn1 (integrate-substep k y-0 y-1 dt (expt 2 (1+ r)) function))
         (err tol))
    (incf r 1)
    (when (> (max kn0 kn1) 0d0)
      (loop for i from 0 to 100
            while
            (and
             ;; (> (max kn0 kn1) 0d0)
             (> (max kn1) 0d0)
             (>= err tol))
            do
               (progn
                 ;; (format t "~D ~E ~E - ~E~%" i kn0 kn1 err)
                 (incf r 1)
                 (setf kn0 kn1)
                 (setf
                  ;; (integrate-substep k y-0 y-1 dt (expt 2 r) function)
                  kn1 (integrate-substep k y-0 y-1 dt (expt 2 (1+ r)) function))
                 ;; (setf err (/ (abs (- kn0 kn1)) (max kn0 kn1)))
                 (setf err (/ (abs (- kn0 kn1)) (max kn1)))
                 ))
      (when (> err tol)
        (format t "Damage integration failed to hit bounds ~E ~E~%" err tol)))
    kn1))

(defun secant-solver (k0 y0 y1 dt func)
  (let ((ymid (* 0.5 (+ y0 y1))))
    (labels ((kmid (k)
               (/ (+ k0 k) 2)
               )
             (fkn (k)
               (+
                k0
                (max 0d0 (* dt (funcall func (kmid k) ymid)))))
             (f (k)
               (-
                (fkn k)
                k)))
      (let* ((kn k0)
             (kn1 (+ (fkn k0) 1d-9))
             (rn (f kn))
             (rn1 (f kn1)))
        ;; (format t "kn ~E - kn1 ~E~%" kn kn1)
        ;; (format t "rn ~E - rn1 ~E~%" rn rn1)
        (loop for i from 0 to 100
              while (and (> (abs (- rn rn1)) 1d-9)
                         (> (abs (- kn kn1)) 0d0)
                         (> (abs rn) 0d0)
                         (> (abs rn1) 0d0)
                         )
              do
                 ;; (format t "iter ~D - error ~E~%" i (abs (- rn rn1)))
                 ;; (format t "kn ~E - kn1 ~E~%" kn kn1)
                 ;; (format t "rn ~E - kn1 ~E~%" rn rn1)
                 (when (> (abs (- kn kn1)) 0d0)
                   (let ((inc (*
                               rn
                               (/ (- kn kn1)
                                  (- rn rn1)
                                  ))))
                     (setf kn1 kn
                           rn1 rn)
                     (setf
                      kn (- kn inc))
                     (setf
                      rn (f kn))
                     ;; (setf k (- kn inc))
                     )))
        (fkn kn)))))



(defun integration-test ()
  (let ((k 0d0)
        (k0 5d0)
        (y0 4.2d0)
        (y1 8d0)
        (dt 0.1d0)
        (tau 1d0)
        (tau-exp 6d0)
        (substeps 2))
    (format t "~%")
    (format t "auto-refine ~E~%"
            (cl-mpm/damage::auto-refine-substepper
             k
             y0
             y1
             dt
             (lambda (k y0 y1 s-dt)
               (cl-mpm/damage::huen-integration k
                                                y0
                                                y1
                                                k0
                                                tau
                                                tau-exp
                                                s-dt))
             :tol 0.5d-4))
    (time (dotimes (i 100000)
            (cl-mpm/damage::analytic-trim
             k
             y0
             y1
             k0
             dt
             (lambda (k y0 y1 dt)
               (cl-mpm/damage::auto-refine-substepper
                k
                y0
                y1
                dt
                (lambda (k y0 y1 s-dt)
                  (cl-mpm/damage::huen-integration k
                                                   y0
                                                   y1
                                                   k0
                                                   tau
                                                   tau-exp
                                                   s-dt))
                :tol 0.5d-4)))))
    (format t "trim auto-refine ~E~%"
            (cl-mpm/damage::analytic-trim
             k
             y0
             y1
             k0
             dt
             (lambda (k y0 y1 dt)
               (cl-mpm/damage::auto-refine-substepper
                k
                y0
                y1
                dt
                (lambda (k y0 y1 s-dt)
                  (cl-mpm/damage::huen-integration k
                                                   y0
                                                   y1
                                                   k0
                                                   tau
                                                   tau-exp
                                                   s-dt))
                :tol 0.5d-4))))
    (time (dotimes (i 100000)
            (cl-mpm/damage::analytic-trim
             k
             y0
             y1
             k0
             dt
             (lambda (k y0 y1 dt)
               (cl-mpm/damage::secant-solver
                k
                y0
                y1
                dt
                (lambda (kmid ymid)
                  (cl-mpm/damage::deriv-partial
                   kmid
                   ymid
                   k0
                   tau
                   tau-exp)))))))
    (format t "substeps ~D ~E~%" substeps
            (cl-mpm/damage::integrate-substep
             k
             y0
             y1
             dt
             substeps
             (lambda (k y0 y1 s-dt)
               (cl-mpm/damage::huen-integration k
                                                y0
                                                y1
                                                k0
                                                tau
                                                tau-exp
                                                s-dt))))
    (format t "analytic substeps ~D ~E~%"
            substeps
            (cl-mpm/damage::analytic-trim
             k
             y0
             y1
             k0
             dt
             (lambda (k y0 y1 dt)
               (cl-mpm/damage::integrate-substep
                k
                y0
                y1
                dt
                substeps
                (lambda (k y0 y1 s-dt)
                  (cl-mpm/damage::huen-integration k
                                                   y0
                                                   y1
                                                   k0
                                                   tau
                                                   tau-exp
                                                   s-dt))))))

    ))
