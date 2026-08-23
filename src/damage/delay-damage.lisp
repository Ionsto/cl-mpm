(in-package :cl-mpm/damage)

;; (defun deriv-partial (k y k0 tau n)
;;   (if (and (>= y k0) (>= y k))
;;       (/
;;        (* k0
;;           (expt
;;            (/ (the double-float (max 0d0 (- y k)))
;;               k0) n))
;;        tau)
;;       0d0))

;; (defun deriv-grad-partial (k y k0 tau n)
;;   (if (and (>= y k0) (>= y k))
;;       (* 1d0
;;          (/
;;           (* n
;;              (if nil;(= 0d0 (max 0d0 (- y k)))
;;                  1d0
;;                  (expt
;;                   (/ (the double-float (max 1d-15 (- y k)))
;;                      k0) (- n 1))))
;;           tau))
;;       1d-15))

;; (defun backwards-integration (k y-0 y-1 k0 tau n dt)
;;   (let* ((r y-1)
;;          (dr 0d0)
;;          (f-tol 1d-10)
;;          (f f-tol)
;;          ;; (grad (deriv-grad-partial k y-0 y-1 k0 tau n dt))
;;          (grad 0d0)
;;          (flow 0d0)
;;          (alpha 1d0)
;;          (r-alpha k)
;;          (y-alpha (+ y-0))
;;          )
;;     (if (>= y-1 k0)
;;         (loop for i from 0 to 10
;;               while (>= (abs f) f-tol)
;;               do
;;                  (progn
;;                    ;; (setf f (+ (- r)
;;                    ;;            k
;;                    ;;            (* dt
;;                    ;;               (deriv-partial r y-1 k0 tau n))
;;                    ;;            ))
;;                    ;; (setf
;;                    ;;  flow
;;                    ;;  (deriv-partial r y-1 k0 tau n)
;;                    ;;  grad
;;                    ;;  (deriv-grad-partial r y-1 k0 tau n)
;;                    ;;  )
;;                    ;; (format t "Step ~D - R ~E - f ~E - grad ~E - flow ~E~%" i r f grad flow)
;;                    (setf r
;;                          (+
;;                           r
;;                           (/
;;                            (-
;;                             (+ k
;;                                (* dt
;;                                   (deriv-partial r y-1 k0 tau n)))
;;                             r
;;                             )
;;                            (+ 1d0 (* dt (deriv-grad-partial r y-1 k0 tau n)))
;;                            ))
;;                          )
;;                    ;; (setf y-alpha (deriv-partial dr y-1 k0 tau n))
;;                    ;; (setf r (max r k))
;;                    (setf f (+ (- r)
;;                               k
;;                               (* dt
;;                                  (deriv-partial r y-1 k0 tau n))
;;                               ))

;;                    ;; (format t "Step ~D - R ~E - f ~E - grad ~E~%" i r f dr)
;;                    )
;;               )
;;         (progn
;;           (setf f 0d0
;;                 r k
;;                 ))
;;         )
;;     ;; (when )
;;     (when (> (abs f) f-tol)
;;       (error "non-convergence"))
;;     r))
;; (defun midpoint-integration (k y-0 y-1 k0 tau n dt)
;;   (let* ((r y-1)
;;          (dr 0d0)
;;          (f-tol 1d-10)
;;          (f f-tol)
;;          ;; (grad (deriv-grad-partial k y-0 y-1 k0 tau n dt))
;;          (grad 0d0)
;;          (flow 0d0)
;;          (alpha 0.5d0)
;;          (r-1 k)
;;          (y (+ (* (- 1d0 alpha) y-0)
;;                (* alpha y-1)))
;;          )
;;     (if (or (>= y-1 k0) (>= y-0 k0))
;;         (loop for i from 0 to 10
;;               while (>= (abs f) f-tol)
;;               do
;;                  (progn
;;                    ;; (setf
;;                    ;;  r
;;                    ;;  (+ (* (- 1d0 alpha) k)
;;                    ;;     (* alpha r-1))
;;                    ;;  )
;;                    ;; (setf f (+ (- r)
;;                    ;;            k
;;                    ;;            (* dt
;;                    ;;               (deriv-partial r y k0 tau n))
;;                    ;;            ))
;;                    ;; (setf
;;                    ;;  flow
;;                    ;;  (deriv-partial r y k0 tau n)
;;                    ;;  grad
;;                    ;;  (deriv-grad-partial r y k0 tau n)
;;                    ;;  )
;;                    ;; (format t "Step ~D - R ~E - f ~E - grad ~E - flow ~E~%" i r f grad flow)
;;                    (setf r
;;                          (+
;;                           r
;;                           (/
;;                            (-
;;                             (+ k
;;                                (* dt
;;                                   alpha
;;                                   (deriv-partial r y k0 tau n)))
;;                             r
;;                             )
;;                            (+ 1d0 (* dt alpha (deriv-grad-partial r y k0 tau n)))
;;                            )))
;;                    (setf f (+ (- r)
;;                               k
;;                               (* dt
;;                                  alpha
;;                                  (deriv-partial r y k0 tau n))
;;                               ))
;;                    )
;;               )
;;         (progn
;;           (setf f 0d0
;;                 r k
;;                 )))
;;     (setf r
;;           (+
;;            k
;;            (* dt
;;               (deriv-partial r y k0 tau n))))
;;     ;; (when )
;;     (when (> (abs f) f-tol)
;;       (error "non-convergence"))
;;     r))

;; (defun huen-integration (k y-0 y-1 k0 tau n dt)
;;   (let* ((dk-0 (deriv-partial k y-0 k0 tau n))
;;          (dk-1 (deriv-partial (+ k (* dt dk-0)) y-1 k0 tau n)))
;;     (+ k (* (/ dt 2) (+ dk-0 dk-1)))))


(defun deriv-partial (k y k0 tau n)
  (declare (double-float k y k0 tau n))
  (if t;(>= y k0)
      (/
       (* k0
          (the double-float
               (expt
                (/
                 (the double-float (max 0d0 (- y (max k k0))))
                 ;; (the double-float (max 0d0 (- y k)))
                 (max
                  k0
                  k
                  )) n)))
       tau)
      0d0)
  ;; (if t;(> y k0)
  ;;     (/
  ;;      (* k0
  ;;         (expt
  ;;          (/ (the double-float (max 0d0 (- y k)))
  ;;             k0) n))
  ;;      tau)
  ;;     0d0)
  )

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
           ((integer 1 1000000) iters))
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

(defun auto-refine-substepper (k y-0 y-1 dt function &key (tol 1d-3))
  (declare (double-float k y-0 y-1 dt)
           (function function)
           )
  (let* ((r 0)
         (kn0 (integrate-substep k y-0 y-1 dt (expt 2 r) function))
         (kn1 (integrate-substep k y-0 y-1 dt (expt 2 (1+ r)) function))
         (err tol))
    (when (> (max kn0 kn1) 0d0)
      (loop for i from 0 to 100
            while
            (and
             (> (max kn0 kn1) 0d0)
             (>= err tol))
            do
               (progn
                 ;; (format t "~D ~E ~E - ~E~%" i kn0 kn1 err)
                 (incf r 1)
                 (setf
                  kn0 kn1
                  ;; (integrate-substep k y-0 y-1 dt (expt 2 r) function)
                  kn1 (integrate-substep k y-0 y-1 dt (expt 2 (1+ r)) function))
                 (setf err (/ (abs (- kn0 kn1)) (max kn0 kn1)))))
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


;; (let ((k 5d0)
;;       (y0 0d0)
;;       (y1 100d0)
;;       (n 2d0)
;;       (tau 1d0)
;;       (dt 100d0)
;;       (k0 1d0)
;;       )
;;   ;; (pprint (forwards-integration k y0 y1 k0 tau n dt))
;;   (pprint (huen-integration k y0 y1 k0 tau n dt))
;;   (pprint (cl-mpm/damage::integrate-substep
;;            k
;;            y0
;;            y1
;;            dt
;;            1000
;;            (lambda (k y0 y1 s-dt)
;;              (cl-mpm/damage::huen-integration k
;;                                               y0
;;                                               y1
;;                                               k0
;;                                               tau
;;                                               n
;;                                               s-dt))
;;            ))
;;   (pprint (cl-mpm/damage::auto-refine-substepper
;;            k
;;            y0
;;            y1
;;            dt
;;            (lambda (k y0 y1 s-dt)
;;              (cl-mpm/damage::huen-integration k
;;                                               y0
;;                                               y1
;;                                               k0
;;                                               tau
;;                                               n
;;                                               s-dt))
;;            ))
;;   ;; (time 
;;   ;;  (pprint (cl-mpm/damage::integrate-substep
;;   ;;           k
;;   ;;           y0
;;   ;;           y1
;;   ;;           dt
;;   ;;           100
;;   ;;           (lambda (k y0 y1 s-dt)
;;   ;;             (cl-mpm/damage::huen-integration k
;;   ;;                                              y0
;;   ;;                                              y1
;;   ;;                                              k0
;;   ;;                                              tau
;;   ;;                                              n
;;   ;;                                              s-dt))
;;   ;;           )))

;;   ;; (pprint (integrate-substep k y0 y1 k0 tau n dt 10000

;;   ;;                            ))

;;   ;; (pprint (secant-solver k y0 y1 dt
;;   ;;                        (lambda (k y)
;;   ;;                          (deriv-partial k y k0 tau n))))
;;   )
