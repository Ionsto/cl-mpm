(in-package :cl-mpm/damage)

(defun deriv-partial (k y k0 tau n)
  (if (> y k0)
      (/
       (* k0
          (expt
           (/ (the double-float (max 0d0 (- y k)))
              k0) n))
       tau)
      0d0))

(defun deriv-grad-partial (k y k0 tau n)
  (if (and ;(> y k0) ;(> y k0)
           )
      (/
       (* n
          (if (= 0d0 (max 0d0 (- y k)))
              1d0
              (expt
               (/ (the double-float (max 0d0 (- y k)))
                  k0) (- n 1))))
       tau)
      0d0))

(defun backwards-integration (k y-0 y-1 k0 tau n dt)
  (let* ((r k)
         (f-tol 1d-5)
         (f f-tol)
         ;; (grad (deriv-grad-partial k y-0 y-1 k0 tau n dt))
         (grad 0d0)
         (deriv 0d0)
         )
    (loop for i from 0 to 100
          while (>= (abs f) f-tol)
          do
             (progn
               (setf deriv
                     (* 1d0 dt (deriv-grad-partial r y-1 k0 tau n dt)))
               (setf grad
                     (* -1d0 deriv))
               (setf f (- r (+ k grad)))
               (when (not (= grad 0d0))
                 (incf r
                       (/ f deriv)
                       ))
               ;; (setf r (max r k))
               ;; (setf f (- r (+ k grad)))
               ;; (setf f (abs f))

               (format t "Step ~D - R ~E - f ~E - grad ~E~%" i r f grad)
               )
          )
    (when (> f f-tol)
      (error "non-convergence"))
    r))

(defun huen-integration (k y-0 y-1 k0 tau n dt)
  (let* ((dk-0 (deriv-partial k y-0 k0 tau n))
         (dk-1 (deriv-partial (+ k (* dt dk-0)) y-1 k0 tau n)))
    (+ k (* (/ dt 2) (+ dk-0 dk-1)))))

(defun forwards-integration (k y-0 y-1 k0 tau n dt)
  (+ k (* dt (deriv-partial k y-0 k0 tau n))))


(defun test (int-scheme)
  (let ((data-k (list 0d0))
        (data-y (list 0d0))
        (data-step (list))
        (data-t (list))
        (dt 0.1d0))
    (loop for time from 0d0 to 10d0 by dt
          for step from 0
          do
             (progn
               (let* ((y-forcing (+ time (sin time)))
                      (kn (funcall int-scheme
                                   (first data-k)
                                   (first data-y)
                                   y-forcing
                                   dt
                                   ))
                      )
                 (push y-forcing data-y)
                 (push kn data-k)
                 (push time data-t))
               )
          )
    (vgplot:plot data-t data-k "K"
                 data-t data-y "Y"
                 )
    )
  )

(defun test-forwards ()
  (let ((t0 1d0)
        (tau 0.1d0)
        (n 1d0))
    ;; (test (lambda (k y0 y1 dt) (forwards-integration k y0 y1 t0 tau n dt)))
    (test (lambda (k y0 y1 dt) (huen-integration k y0 y1     t0 tau n dt)))
    )

 )
