(in-package :cl-mpm/damage)

(declaim (ftype (function (double-float double-float double-float double-float )
                          double-float) damage-response-exponential))
(defun damage-response-exponential (stress E init-stress ductility)
  (declare (double-float stress E init-stress ductility))
  "Function that controls how damage evolves with principal stresses"
  (let* ((ft init-stress)
         (e0 (/ ft E))
         (ef (/ (* ft (+ ductility 1d0)) (* 2d0 E)))
         (k (/ stress E))
         (beta (/ 1d0 (- ef e0)))
         )
    (declare (double-float ft e0 ef k beta))
    (if (> k e0)
        (- 1d0 (* (/ e0 k) (exp (- (* beta (- k e0))))))
        0d0)))

(defun find-k-damage-mp (mp damage)

  (with-accessors
        ((E cl-mpm/particle::mp-e)
         (init-stress cl-mpm/particle::mp-initiation-stress)
         (ductility cl-mpm/particle::mp-ductility))
      mp
    (find-k-damage E init-stress ductility damage))
  )
(defun find-k-damage (E init-stress ductility damage)
  (let* ((e0 init-stress)
         (d-est 0d0)
         (k-est e0)
         (k-max nil)
         (k-prev e0)
         )
    (loop for i from 0 to 10000
          while (> (abs (- damage d-est)) 1d-8)
          do
             (progn
               (setf d-est (damage-response-exponential k-est E init-stress ductility))
               ;; (format t "~F ~F~%" d-est k-est)
               (when (> damage d-est)
                 (setf k-prev k-est)
                 (setf k-est
                       (if k-max
                           (/ (+ k-prev k-max) 2d0)
                           (* k-est 2))))
               (when (< damage d-est)
                 (setf k-max k-est)
                 (setf k-est (/ (+ k-prev k-max) 2d0)))
               ))
    k-est))

(declaim (ftype (function (double-float double-float double-float double-float double-float)
                          double-float) damage-response-exponential-peerlings-residual))
(defun damage-response-exponential-peerlings-residual (stress E init-stress ductility residual)
  (declare (double-float stress E init-stress ductility residual))
  "Function that controls how damage evolves with principal stresses"
  (let* ((ft init-stress)
         (e0 (/ ft E))
         (ef (/ (* ft (+ ductility 1d0)) (* 2d0 E)))
         (k (/ stress E))
         (beta (/ 1d0 (- ef e0))))
    (declare (double-float ft e0 ef k beta))
    (if (> k e0)
        (- 1d0 (* (/ e0 k) (+ (- 1d0 residual) (* residual (exp (- (* beta (- k e0))))))))
        0d0)))
(declaim (ftype (function (double-float double-float double-float double-float double-float)
                          double-float) damage-response-exponential-residual))
(defun damage-response-exponential-residual (stress E init-stress ductility residual)
  (declare (double-float stress E init-stress ductility residual))
  "Function that controls how damage evolves with principal stresses"
  (let* ((ft init-stress)
         (e0 (/ ft E))
         (ef (/ (* ft (+ ductility 1d0)) (* 2d0 E)))
         (k (/ stress E))
         (beta (/ 1d0 (- ef e0))))
    (declare (double-float ft e0 ef k beta))
    (if (> k e0)
        (- 1d0 (+ (- 1d0 residual) (* residual (exp (- (* beta (- k e0)))))))
        0d0)))

(defun damage-response-exponential-beta (E init-stress ductility)
  (declare (double-float E init-stress ductility))
  "Function that controls how damage evolves with principal stresses"
  (let* ((ft init-stress)
         (e0 (/ ft E))
         (ef (/ (* ft (+ ductility 1d0)) (* 2d0 E)))
         (beta (/ 1d0 (- ef e0))))
    beta))


(defun estimate-ductility-jirsek2004 (GF R ft E &optional (k 1d0))
  "Simple estimate of ductility from key elastic parameters"
  (declare (double-float GF R ft E k))
  (let* ((e0 (/ ft E))
         (ef (+ (/ GF (* k R E e0)) (/ e0 2d0))))
    (- (* 2d0 (/ ef e0)) 1d0)))

(defun gf-from-ductility (ductility R ft E &optional (k 1d0))
  (declare (double-float ductility R ft E k))
  "Sanity check for estimating fracture energy from ductility and elastic parameters"
  (let ((e0 (/ ft E)))
    (/ (* ductility R E (expt e0 2))
       2d0)))

(defun compute-oversize-factor-residual (damage-final E init-stress ductility residual)
  "Compute the estimated oversize factor (i.e. w * init-stress) required to have plasticity kick in at (d=damage-final)"
  (let* ((ft init-stress)
         (e0 (/ ft E))
         (ef (/ (* ft (+ ductility 1d0)) (* 2d0 E)))
         (beta (/ 1d0 (- ef e0)))
         (y (- 1d0 damage-final))
         )
    (* (/ 1d0 e0)
       (/ (+
           (* y (- ef e0) (cl-mpm/fastmaths::lambert-w-0 (/ (* -1d0 e0 residual (exp (/ (* -1d0 e0 (+ residual y -1d0))
                                                                                          (* y (- e0 ef)))))
                                                              (* y (- e0 ef)))))
             (* e0 -1d0 residual)
             e0
             )
          y))))


(defun compute-oversize-factor (damage-final ductility)
  "Compute the estimated oversize factor (i.e. w * init-stress) required to have plasticity kick in at (d=damage-final)"
  (declare (optimize (speed 0) (debug 3)))
  (when (< ductility 1d0)
    (error "Ductility ~A cannot be less than 1" ductility))
  (let* ((eta ductility)
         (value (/ (* 2d0 (exp (/ 2d0 (- eta 1d0))))
                   (* (- 1d0 damage-final) (- eta 1d0)))))
    (* 0.5d0 (- eta 1d0) (cl-mpm/fastmaths::lambert-w-0 value))))

(defun brittle-concrete-linear-d (stress E Gf length init-stress)
  (let* ((hsl (/ (expt init-stress 2) (* 2 E Gf)))
         (hs (/ (* hsl length) (- 1 (* hsl length)))))
    (if (> stress init-stress)
        (* (+ 1d0 hs) (- 1d0 (/ init-stress stress)))
        0d0
        )))

(defun brittle-concrete-d (stress E Gf length init-stress)
  (declare (double-float stress E Gf length init-stress))
  "Function that controls how damage evolves with principal stresses"
  (let* (
         (ft init-stress)
         (e0 (/ ft E))
         (ef (+ (/ e0 2) (/ Gf (* length ft))))
         (k (/ stress E))
         )

    (when (> length (/ (* 2 Gf E) (expt ft 2)))
      (error "Length scale is too long"))
    (if (> k e0)
        (- 1d0 (* (/ e0 k) (exp (- (/ (- k e0) (- ef e0)))))) 0d0)))

(defun delayed-damage-y (stress E Gf length init-stress)
  (declare (double-float stress E Gf length init-stress))
  "Function that controls how damage evolves with principal stresses"
  (let* ((ft init-stress)
         (e0 (/ ft E))
         (ef (+ (/ e0 2) (/ Gf (* length ft))))
         (k (/ stress E)))
    (if (> k e0)
        (/ (- stress init-stress) (- (* ef E) init-stress))
        0d0
        )))

(defun brittle-chalk-d (stress E Gf length init-stress)
  (declare (double-float stress E Gf length init-stress))
  "Function that controls how damage evolves with principal stresses"
  ;; (if (> stress init-stress)
  ;;     (* (expt (max 0d0 (- stress init-stress)) 0.5d0) rate) 0d0)
  ;; (when (> length (/ (* 2 Gf E) (expt init-stress 2)))
  ;;   (error "Length scale is too long"))
  (let* ((hs (/ (expt init-stress 2) (* 2d0 E Gf)))
         (hsl (/ (* hs length) (- 1d0 (* hs length)))))
    (if (> stress init-stress)
        (min 1d0 (max 0d0
                      (- 1d0 (* (/ init-stress stress) (exp (* -2d0 hs (/ (- stress init-stress) stress)))))
                      ))
        0d0)
    ))

(defun damage-rate-profile-chalk (stress damage rate init-stress)
  (declare (double-float stress damage rate init-stress))
  "Function that controls how damage evolves with principal stresses"
  (if (> stress init-stress)
      (* (expt (max 0d0 (- stress init-stress)) 0.5d0) rate)
      0d0))


(defun deriv-partial (k y k0 tau n)
  (if (> y k0)
      (/
       (* k0
          (expt
           (/ (the double-float (max 0d0 (- y k)))
              k0) n))
       tau)
      0d0))

(defun huen-integration (k y-0 y-1 k0 tau n dt)
  (let* ((dk-0 (deriv-partial k y-0 k0 tau n))
         (dk-1 (deriv-partial (+ k (* dt dk-0)) y-0 k0 tau n)))
    (+ k (* (/ dt 2) (+ dk-0 dk-1)))))

(defun forwards-integration (k y-0 y-1 k0 tau n dt)
  (+ k (* dt (deriv-partial k y-0 k0 tau n))))
