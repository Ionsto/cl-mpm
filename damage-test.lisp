(ql:quickload "lisp-stat")
(ql:quickload "vgplot")

(defun stress (s0 tau time)
  s0)

(defun damage-rate (alpha beta sc d s)
  (* alpha
     (expt (max 0 (- (/ (max s 0) (- 1 d)) sc)) beta)))

(defun sim (tfinal dt s0 sc tau alpha beta)
  (let ((tn 0)
        (dn 0)
        (t-data (list))
        (s-data (list))
        (d-data (list))
        (critical-damage 0.56)
        )
    (loop for i from 0 to (floor tfinal dt)
          when (< dn critical-damage)
            do
               (progn
                 (push tn t-data)
                 (push dn d-data)
                 (push (/ (stress s0 tau tn)
                          (- 1 dn)) s-data)
                 (incf tn dt)
                 (incf
                  dn
                  (* dt
                     (damage-rate alpha beta sc
                                  dn
                                  (stress s0 tau tn))))
                 ))
    (format t "Time failure ~a at ~f~%" tn s0)
    (format t "Time failure ~a h at ~f MPa~%" (/ tn (* 60.0 60)) (/ s0 1d6))
    (format t "Cycles to failure ~a at ~f~%" (/ tn tau) s0)
    (values tn t-data d-data s-data)))

(defun last-array (a)
  (aref a (- (length a) 1)))

(defun plot-sim (s0 &optional (alpha 2d-24))
  (multiple-value-bind (tn t-data d-data s-data)
      (sim (* 800 60 60) 1 s0 0.33e6 8 alpha 3.25d0)
    (vgplot:close-all-plots)
    (vgplot:figure)
    (vgplot:title "Time vs damage")
    (vgplot:plot t-data s-data)
    (let ((df (lisp-stat:read-csv
	             (uiop:read-file-string #P"creep.csv"))))
      (vgplot:plot
       (mapcar (lambda (x) (/ x (* 60 60))) t-data) d-data "mpm"
       (lisp-stat:column df 'time) (lisp-stat:column df 'damage) "data"))
    (vgplot:figure)
    (vgplot:title "Normalised time vs damage")
    ;(vgplot:plot t-data s-data)
    (let ((df (lisp-stat:read-csv
	             (uiop:read-file-string #P"creep.csv"))))
      (vgplot:plot
       (mapcar (lambda (x) (/ x (first t-data))) t-data) d-data "mpm"
       (aops:each (lambda (x) (/ x (last-array (lisp-stat:column df 'time))))
                  (lisp-stat:column df 'time)) (lisp-stat:column df 'damage) "data"))
    ))

;; (defun find-alpha ()
;;   (let* (
;;          (s0 275e3)
;;          (tau 8.0)
;;          (n 200)
;;          (above 0)
;;          (below 10)
;;          (last-above t)
;;          (last-n 0)
;;          (alpha above)
;;          )
;;     (loop for i from 0 to 1000
;;           when (> (abs (- last-n n)) 0.1)
;;             do
;;                (progn
;;                  (if last-above
;;                      (setf above alpha)
;;                      (setf below alpha))
;;                  (setf alpha (/ (+ below above) 2))
;;                  (let ((tnf (sim 10000 0.01 s0 250e3 tau alpha)))
;;                    (setf last-n (/ tnf tau))
;;                    (setf last-above (> last-n n))
;;                    (format t "Iteration ~D: alpha ~f above?: ~a n: ~a ~%" i alpha last-above (/ tnf tau))
;;                    )))
;;     alpha))
(defun find-alpha-time ()
  (let* (
         (s0 0.93e6)
         (beta 3.25d0)
         (tau 8.0)
         (n (* 115.22 60 60))
         (above 1e-28)
         (below 1e-20)
         (last-above t)
         (last-n 0)
         (alpha above)
         )
    (loop for i from 0 to 1000
          when (> (abs (- last-n n)) 0.01)
            do
               (progn
                 (if last-above
                     (setf above alpha)
                     (setf below alpha))
                 (setf alpha (/ (+ below above) 2))
                 (let ((tnf
                         (sim (* 200 60 60) 1 s0 0.33e6 8 alpha beta)
                         ;; (sim 10000 0.01 s0 250e3 tau alpha)
                         ))
                   ;; (setf last-n (/ tnf tau))
                   (setf last-n tnf)
                   (setf last-above (> last-n n))
                   (format t "Iteration ~D: alpha ~f above?: ~a n: ~a ~%" i alpha last-above (/ tnf tau))
                   )))
    (plot-sim s0 alpha)
    alpha))

(defparameter *alpha* (find-alpha-time))

(defun make-sn-curve ()
  (let* ((s (loop for s from 260e3 to 500e3 by 10e3
                  collect s))
         (n
           (loop for s0 in s
                 collect
                 (multiple-value-bind (tn a b c)
                     (sim 100000 0.01 s0 250e3 8 *alpha*)
                   (/ tn 8)))))
    (vgplot:figure)
    (let ((df (lisp-stat:read-csv
	             (uiop:read-file-string #P"fatigue_curve.csv"))))
      (vgplot:xlabel "N")
      (vgplot:ylabel "Stress")
      (vgplot:semilogx n (mapcar (lambda (x) (/ x 1d3)) s) "mpm"
                       (lisp-stat:column df 'n) (lisp-stat:column df 's) "data"))
    ))
(defun plot-fatigue-curve ()
  (let ((df (lisp-stat:read-csv
	          (uiop:read-file-string #P"fatigue_curve.csv"))))
    (vgplot:figure)
    (vgplot:semilogx (lisp-stat:column df 'n) (lisp-stat:column df 's) ))
  )

(defun find-beta-time ()
  (let* ((s0 0.93e6)
         (n (* 115.22 60 60))
         (alpha 1d-05)
         (betas (loop for b from 2.5 to 3.5 by 0.25d0 collect b))
         (df (lisp-stat:read-csv (uiop:read-file-string #P"creep.csv")))
         )
    (loop for b in betas
            do
               (progn
                 (format t "Beta value: ~f~%" b)
                 (multiple-value-bind (tn t-data d-data s-data)
                     ;; '(1 1 1 1)
                     (sim (* 1000 60 60) 10 s0 0.33e6 8
                          (/ alpha (expt 1000 (* 2 b))) b)
                   (vgplot:figure)
                   (vgplot:plot
                    (mapcar (lambda (x) (/ x (first t-data))) t-data) d-data "mpm"
                    (aops:each (lambda (x) (/ x (last-array (lisp-stat:column df 'time))))
                               (lisp-stat:column df 'time)) (lisp-stat:column df 'damage) "data")            (vgplot:title b))
                 (sleep 1)
                 ))
    ;; (plot-sim s0 alpha)
    alpha))
