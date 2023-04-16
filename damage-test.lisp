
(defun stress (s0 tau time)
  (* s0 (sin (* 2 pi (/ time tau)))))
(defun damage-rate (alpha sc d s)
  (* alpha (max 0 (- (/ s (- 1 d)) sc))))
(defun sim (tfinal dt s0 sc tau alpha)
  (let ((tn 0)
        (dn 0)
        (t-data (list))
        (s-data (list))
        (d-data (list))
        (critical-damage 0.5)
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
                     (damage-rate alpha sc
                                  dn
                                  (stress s0 tau tn))))
                 ))
    (format t "Cycles to failure ~a~%" (/ tn tau))
    (values tn t-data d-data s-data)))
(defun plot-sim (s0)
  (multiple-value-bind (tn t-data d-data s-data)
      (sim 10000 0.01 s0 250e3 8 3d-8)
    (vgplot:figure)
    (vgplot:title "Time vs stress")
    (vgplot:plot t-data s-data)
    (vgplot:figure)
    (vgplot:title "Time vs damage")
    (vgplot:plot t-data d-data)
    ))

(defun find-alpha ()
  (let* (
         (s0 275e3)
         (tau 8.0)
         (n 100)
         (above 0)
         (below 10)
         (last-above t)
         (last-n 0)
         (alpha above)
         )
    (loop for i from 0 to 1000
          when (> (abs (- last-n n)) 1)
            do
               (progn
                 (if last-above
                     (setf above alpha)
                     (setf below alpha))
                 (setf alpha (/ (+ below above) 2))
                 (let ((tnf (sim 10000 0.1 s0 250e3 tau alpha)))
                   (setf last-n (/ tnf tau))
                   (setf last-above (> last-n n))
                   (format t "Iteration ~D: alpha ~f above?: ~a n: ~a ~%" i alpha last-above (/ tnf tau))
                   )))))
