(in-package :cl-mpm/damage)

(defun deriv-partial (k y k0 tau n)
  (if (and (>= y k0) (>= y k))
      (/
       (* k0
          (expt
           (/ (the double-float (max 0d0 (- y k)))
              k0) n))
       tau)
      0d0))

(defun deriv-grad-partial (k y k0 tau n)
  (if (and (>= y k0) (>= y k))
      (* 1d0
         (/
          (* n
             (if nil;(= 0d0 (max 0d0 (- y k)))
                 1d0
                 (expt
                  (/ (the double-float (max 1d-15 (- y k)))
                     k0) (- n 1))))
          tau))
      1d-15))

(defun backwards-integration (k y-0 y-1 k0 tau n dt)
  (let* ((r y-1)
         (dr 0d0)
         (f-tol 1d-10)
         (f f-tol)
         ;; (grad (deriv-grad-partial k y-0 y-1 k0 tau n dt))
         (grad 0d0)
         (flow 0d0)
         (alpha 1d0)
         (r-alpha k)
         (y-alpha (+ y-0))
         )
    (if (>= y-1 k0)
        (loop for i from 0 to 10
              while (>= (abs f) f-tol)
              do
                 (progn
                   ;; (setf f (+ (- r)
                   ;;            k
                   ;;            (* dt
                   ;;               (deriv-partial r y-1 k0 tau n))
                   ;;            ))
                   ;; (setf
                   ;;  flow
                   ;;  (deriv-partial r y-1 k0 tau n)
                   ;;  grad
                   ;;  (deriv-grad-partial r y-1 k0 tau n)
                   ;;  )
                   ;; (format t "Step ~D - R ~E - f ~E - grad ~E - flow ~E~%" i r f grad flow)
                   (setf r
                         (+
                          r
                          (/
                           (-
                            (+ k
                               (* dt
                                  (deriv-partial r y-1 k0 tau n)))
                            r
                            )
                           (+ 1d0 (* dt (deriv-grad-partial r y-1 k0 tau n)))
                           ))
                         )
                   ;; (setf y-alpha (deriv-partial dr y-1 k0 tau n))
                   ;; (setf r (max r k))
                   (setf f (+ (- r)
                              k
                              (* dt
                                 (deriv-partial r y-1 k0 tau n))
                              ))

                   ;; (format t "Step ~D - R ~E - f ~E - grad ~E~%" i r f dr)
                   )
              )
        (progn
          (setf f 0d0
                r k
                ))
        )
    ;; (when )
    (when (> (abs f) f-tol)
      (error "non-convergence"))
    r))
(defun midpoint-integration (k y-0 y-1 k0 tau n dt)
  (let* ((r y-1)
         (dr 0d0)
         (f-tol 1d-10)
         (f f-tol)
         ;; (grad (deriv-grad-partial k y-0 y-1 k0 tau n dt))
         (grad 0d0)
         (flow 0d0)
         (alpha 0.5d0)
         (r-1 k)
         (y (+ (* (- 1d0 alpha) y-0)
               (* alpha y-1)))
         )
    (if (or (>= y-1 k0) (>= y-0 k0))
        (loop for i from 0 to 10
              while (>= (abs f) f-tol)
              do
                 (progn
                   ;; (setf
                   ;;  r
                   ;;  (+ (* (- 1d0 alpha) k)
                   ;;     (* alpha r-1))
                   ;;  )
                   ;; (setf f (+ (- r)
                   ;;            k
                   ;;            (* dt
                   ;;               (deriv-partial r y k0 tau n))
                   ;;            ))
                   ;; (setf
                   ;;  flow
                   ;;  (deriv-partial r y k0 tau n)
                   ;;  grad
                   ;;  (deriv-grad-partial r y k0 tau n)
                   ;;  )
                   ;; (format t "Step ~D - R ~E - f ~E - grad ~E - flow ~E~%" i r f grad flow)
                   (setf r
                         (+
                          r
                          (/
                           (-
                            (+ k
                               (* dt
                                  alpha
                                  (deriv-partial r y k0 tau n)))
                            r
                            )
                           (+ 1d0 (* dt alpha (deriv-grad-partial r y k0 tau n)))
                           )))
                   (setf f (+ (- r)
                              k
                              (* dt
                                 alpha
                                 (deriv-partial r y k0 tau n))
                              ))
                   )
              )
        (progn
          (setf f 0d0
                r k
                )))
    (setf r
          (+
           k
           (* dt
              (deriv-partial r y k0 tau n))))
    ;; (when )
    (when (> (abs f) f-tol)
      (error "non-convergence"))
    r))

(defun huen-integration (k y-0 y-1 k0 tau n dt)
  (let* ((dk-0 (deriv-partial k y-0 k0 tau n))
         (dk-1 (deriv-partial (+ k (* dt dk-0)) y-1 k0 tau n)))
    (+ k (* (/ dt 2) (+ dk-0 dk-1)))))

(defun forwards-integration (k y-0 y-1 k0 tau n dt)
  (+ k (* dt (deriv-partial k y-0 k0 tau n))))


(defun test (int-scheme data-time data-y)
  (let ((data-k (list 0d0))
        (data-step (list))
        (y-n 0d0)
        (dt (- (second data-time) (first data-time)))
        )
    (loop for time in data-time
          for y-n1 in data-y
          for step from 0
          do
             (progn
               (let* ((kn (funcall int-scheme
                                   (first data-k)
                                   y-n
                                   y-n1
                                   dt
                                   ))
                      )
                 (setf y-n y-n1)
                 (push kn data-k))
               )
          )
    (reverse data-k)))

;; (defun test-forwards ()
;;   (let* ((t0 1d0)
;;          (tau 1d0)
;;          (dt 1d0)
;;          (n 2d0)
;;          (data-coarse-time (loop for time from 0d0 to 10d0 by 0.1d0 collect time))
;;          (data-coarse-y (loop for time in data-coarse-time collect (+ time (sin time))))
;;          (data-time (loop for time from 0d0 to 10d0 by dt collect time))
;;          (data-y (loop for time in data-time collect (+ time (sin time)))))
;;     (let ((test-f (test
;;                    (lambda (k y0 y1 dt) (forwards-integration k y0 y1 t0 tau n dt))
;;                    data-time
;;                    data-y))
;;           (test-huen (test
;;                    (lambda (k y0 y1 dt) (huen-integration k y0 y1 t0 tau n dt))
;;                    data-time
;;                    data-y))
;;           (test-back (test
;;                       (lambda (k y0 y1 dt) (backwards-integration k y0 y1 t0 tau n dt))
;;                       data-time
;;                       data-y))
;;           (test-mid (test
;;                       (lambda (k y0 y1 dt) (midpoint-integration k y0 y1 t0 tau n dt))
;;                       data-time
;;                       data-y))
;;           (test-coarse-euler (test
;;                               (lambda (k y0 y1 dt) (forwards-integration k y0 y1 t0 tau n dt))
;;                              data-coarse-time
;;                              data-coarse-y))
;;           (test-coarse-huen (test
;;                              (lambda (k y0 y1 dt) (huen-integration k y0 y1 t0 tau n dt))
;;                              data-coarse-time
;;                              data-coarse-y))
;;           (test-coarse-back (test
;;                              (lambda (k y0 y1 dt) (backwards-integration k y0 y1 t0 tau n dt))
;;                              data-coarse-time
;;                              data-coarse-y))
;;           )
;;       (vgplot:plot
;;        ;; data-time data-y "Y"
;;        data-coarse-time data-coarse-y "Y coarse"
;;        data-time test-f "K-euler"
;;        data-time test-huen "K-huen"
;;        data-time test-back "K-back"
;;        data-time test-mid "K-mid"
;;        data-coarse-time test-coarse-euler "K-euler-coarse"
;;        ;; data-coarse-time test-coarse-huen "K-huen-coarse"
;;        ;; data-coarse-time test-coarse-back "K-back-coarse"
;;        ))
    
;;     ;; (test (lambda (k y0 y1 dt) (huen-integration k y0 y1     t0 tau n dt)))
;;     ;; (test (lambda (k y0 y1 dt) (backwards-integration k y0 y1     t0 tau n dt)))
;;     )

;;  )
