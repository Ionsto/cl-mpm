(defpackage :cl-mpm/examples/notch
  (:use :cl))
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)
(in-package :cl-mpm/examples/notch)
;(declaim (optimize (debug 0) (safety 0) (speed 3)))


(defun find-max-cfl (sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps)
                   (dt cl-mpm:sim-dt))
      sim
    (let ((max-v (loop for mp across mps
                       maximize (with-accessors ((vel cl-mpm/particle:mp-velocity)) mp
                                  (magicl::sum (magicl:map #'abs vel))))))
      (* dt (/ max-v (cl-mpm/mesh:mesh-resolution mesh))))))

(defun pescribe-velocity (sim load-mps vel)
  (let ((mps (cl-mpm:sim-mps sim)))
    (loop for mp in load-mps
          do
             (progn
               (loop for v in vel
                     for i from 0
                     do
                     (when v
                       (setf (magicl:tref (cl-mpm/particle:mp-velocity mp) i 0) v)
                       )
                     )))))
(defun increase-load (sim load-mps amount)
  (loop for mp in load-mps
        do (with-accessors ((pos cl-mpm/particle:mp-position)
                            (force cl-mpm/particle:mp-body-force)) mp
             ;(incf (magicl:tref force 0 0) amount)
             (magicl:.+ force amount force)
             )))

(defun length-from-def (sim mp dim)
  (let* ((mp-scale 2)
         (h-initial (magicl:tref (cl-mpm/particle::mp-domain-size mp) dim 0))
         ;(h-initial  (/ (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) mp-scale))
         )
    h-initial
    ;(* (magicl:tref (cl-mpm::mp-deformation-gradient mp) dim dim) h-initial)
    ))
(defun max-stress (mp)
  (multiple-value-bind (l v) (magicl:eig (cl-mpm::voight-to-matrix (cl-mpm/particle:mp-stress mp)))
    (apply #'max l)
    (magicl:tref (cl-mpm/particle:mp-stress mp) 2 0)))
(defun plot (sim &optional (plot :point))
  (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  (multiple-value-bind (x y c stress-y lx ly e density temp vx)
    (loop for mp across (cl-mpm:sim-mps sim)
          collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
          collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
          collect (magicl:tref (cl-mpm::mp-velocity mp) 0 0) into vx
          collect (length-from-def sim mp 0) into lx
          collect (length-from-def sim mp 1) into ly
          collect (if (slot-exists-p mp 'cl-mpm/particle::damage) (cl-mpm/particle:mp-damage mp) 0) into c
          collect (if (slot-exists-p mp 'cl-mpm/particle::temperature) (cl-mpm/particle:mp-temperature mp) 0) into temp
          collect (if (slot-exists-p mp 'cl-mpm/particle::strain-energy-density) (cl-mpm/particle::mp-strain-energy-density mp) 0) into e
          collect (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)) into density
          ;; collect (cl-mpm/particle:mp-volume mp) into density
          collect (max-stress mp) into stress-y
          finally (return (values x y c stress-y lx ly e density temp vx)))
    (cond
      ((eq plot :point)
       ;; (vgplot:format-plot t "set cbrange [0:1]")
       (vgplot:plot x y ";;with points pt 7")
       )
      ((eq plot :damage)
       (vgplot:format-plot t "set cbrange [0:1]")
       ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 0.01 (apply #'max c)))
       (vgplot:plot x y c ";;with points pt 7 lc palette")
       )
      ((eq plot :velocity)
       (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min vx) (+ 0.01 (apply #'max vx)))
       (vgplot:plot x y vx ";;with points pt 7 lc palette")
       )
      ((eq plot :temperature)
       (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min temp) (+ 0.01 (apply #'max temp)))
       ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 0.01 (apply #'max c)))
       (vgplot:plot x y temp ";;with points pt 7 lc palette")
       )
      ((eq plot :energy)
       (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min e) (+ 0.01 (apply #'max e)))
       (vgplot:plot x y e ";;with points pt 7 lc palette")
       )
      (
       (eq plot :stress)
       (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 0.01 (apply #'max stress-y)))
       (vgplot:plot x y stress-y ";;with points pt 7 lc palette"))
      ((eq plot :deformed)
       ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 0.01 (apply #'max stress-y)))
       (vgplot:plot x y lx ly ";;with ellipses"))
      ((eq plot :density)
       (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min density) (+ 0.01 (apply #'max density)))
       (vgplot:plot x y e ";;with points pt 7 lc palette")
       ))
    )
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:axis (list 0 ms-x
                       0 ms-y))
    (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
    (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
      (vgplot:format-plot t "set ytics ~f" h)
      (vgplot:format-plot t "set xtics ~f" h))
  (vgplot:replot))

(defun remove-sdf (sim sdf)
      (setf (cl-mpm:sim-mps sim)
            (remove-if (lambda (mp)
                         (with-accessors ((pos cl-mpm/particle:mp-position)) mp
                           (>= 0 (funcall sdf pos))
                           ))
                       (cl-mpm:sim-mps sim))))

(defun damage-sdf (sim sdf)
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
      (loop for mp across mps
            do (with-accessors ((pos cl-mpm/particle:mp-position)
                                 (damage cl-mpm/particle:mp-damage)) mp
                  (when (>= 0 (funcall sdf pos))
                    (setf damage 1d0))))))

(defun rectangle-sdf (position size)
  (lambda (pos)
      (let* ((position (magicl:from-list position '(2 1) :type 'double-float))
             (dist-vec (magicl:.- (magicl:map! #'abs (magicl:.- pos position))
                                  (magicl:from-list size '(2 1) :type 'double-float))))

        (+ (sqrt (magicl::sum
                  (magicl:map! (lambda (x) (* x x))
                               (magicl:map! (lambda (x) (max 0d0 x)) dist-vec))))
           (min (max (magicl:tref dist-vec 0 0)
                     (magicl:tref dist-vec 1 0)
                     ) 0d0)))))
(defun ellipse-sdf (position x-l y-l)
  (let ((aspect (/ x-l y-l)))
    (lambda (pos)
      (let* ((position (magicl:from-list position '(2 1) :type 'double-float))
             (dist-vec (magicl:.* (magicl:.- position pos) (magicl:from-list (list 1 aspect) '(2 1)
                                                                             :type 'double-float)))
             (distance (sqrt (magicl:tref (magicl:@ (magicl:transpose dist-vec)
                                                    dist-vec) 0 0))))
        (- distance x-l)))))

(defun make-mps-bread (density-cheese
                       density-bread
                       cheese-size
                       bread-size
                       e-scale
                       mp-scale
                       h)
  (let* ((mass-bread (* density-bread (expt (/ h mp-scale) 2)))
         (mass (* density-cheese (expt (/ h mp-scale) 2)))
         (cheese-pos (list 0 (second bread-size)))
         (bread-1-pos (list 0 0))
         (bread-2-pos (list 0 (+ (second bread-size)
                                 (second cheese-size))))
         )
    (cl-mpm/setup::make-mps-from-list
     (append
      (cl-mpm/setup::make-block-mps-list
       (cl-mpm/setup::node-to-mp-pos cheese-pos (list h h) mp-scale)
       cheese-size
       (mapcar (lambda (e) (* e e-scale mp-scale)) cheese-size)
       'cl-mpm::make-particle
       'cl-mpm/particle::particle-thermoviscoplastic-damage
       :E 1e7
       :nu 0.325d0
       :visc-factor 1d2
       :visc-power 3d0
       :temperature 0d0
       :heat-capacity 1d0
       :thermal-conductivity 1e0
       :mass mass
       :critical-stress 1d20
       :gravity -9.8d0
       :index 0
       )
      (cl-mpm/setup::make-block-mps-list
       (cl-mpm/setup::node-to-mp-pos bread-1-pos
                                     (list h h)
                                     mp-scale)
       bread-size
       (mapcar (lambda (e) (* e e-scale mp-scale))
               bread-size)
       'cl-mpm::make-particle
       'cl-mpm/particle::particle-thermoelastic-damage
       :E 1e7
       :nu 0.325
       :temperature 0d0
       :heat-capacity 1d0
       :thermal-conductivity 1e0
       :mass mass-bread
       :critical-stress 1d20
       :gravity -9.8d0
       :index 1
       )
      (cl-mpm/setup::make-block-mps-list
       (cl-mpm/setup::node-to-mp-pos bread-2-pos
                                     (list h h)
                                     mp-scale)
       bread-size
       (mapcar (lambda (e) (* e e-scale mp-scale))
               bread-size)
       'cl-mpm::make-particle
       'cl-mpm/particle::particle-thermoelastic-damage
       :E 1e7
       :nu 0.325
       :temperature 0d0
       :heat-capacity 1d0
       :thermal-conductivity 1e0
       :mass mass-bread
       :critical-stress 1d20
       :gravity -9.8d0
       :index 1
       )
      ))))

(defun setup-test-column (size block-size block-offset &optional (e-scale 1d0) (mp-scale 1d0)
                          &rest mp-args)
  (let* ((sim (cl-mpm/setup::make-block (/ 1 e-scale)
                                        (mapcar (lambda (s) (* s e-scale)) size)
                                        #'cl-mpm/shape-function:make-shape-function-bspline)) 
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         ;(e-scale 1)
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 900)
         ;; (mass (/ (* 900 h-x h-y) (expt mp-scale 2)))
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size)))
    (progn
      (let ((block-position
              (mapcar #'+ (list (* h-x (- (+ (/ 1 (* 2 mp-scale))) 0))
                                (* h-y (+ (/ 1d0 (* 2d0 mp-scale)))))
                      block-offset)))
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-block-mps
               block-position
               block-size
               (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
               density
               'cl-mpm::make-particle
               ;'cl-mpm/particle::particle-viscoplastic-damage
               ;; 'cl-mpm/particle::particle-viscoplastic
                'cl-mpm/particle::particle-elastic-damage
                 :E 1d9
                 :nu 0.3250d0
                 ;; :visc-factor 1d6
                 ;; :visc-power 3d0
                 :critical-stress 1d10
                 :gravity -9.8d0
                 ;; :gravity-axis (magicl:from-list '(0.5d0 0.5d0) '(2 1))
                 :index 0
               )))
      (setf (cl-mpm:sim-damping-factor sim) 0.5d0)
      (setf (cl-mpm:sim-mass-filter sim) 1d-15)
      (setf (cl-mpm::sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm:sim-dt sim) 1d-2)
      ;; (lambda (i) (cl-mpm/bc:make-bc-friction i
             ;; (magicl:from-list '(0d0 1d0) '(2 1)) 0.25d0))
      ;; (setf (cl-mpm::sim-bcs-force sim)
      ;;        (cl-mpm/bc::make-outside-bc-var
      ;;         (cl-mpm:sim-mesh sim)
      ;;         nil
      ;;         nil
      ;;         nil
      ;;        (lambda (i)
      ;;          ;; (cl-mpm/bc:make-bc-friction i (magicl:from-list '(0d0 1d0) '(2 1)) 0.0d0)
      ;;          (if (< (* h (nth 0 i)) 250d0)
      ;;                (cl-mpm/bc:make-bc-friction i (magicl:from-list '(0d0 1d0) '(2 1)) 0.6d0)
      ;;                (cl-mpm/bc:make-bc-friction i (magicl:from-list '(0d0 1d0) '(2 1)) 0.01d0)
      ;;                )
      ;;                ;; (cl-mpm/bc:make-bc-friction i (magicl:from-list '(0d0 1d0) '(2 1)) 0.01d0))
      ;;          )
      ;;         ))
            ;; (append
            ;;  (cl-mpm/bc::make-domain-bcs
            ;;   (cl-mpm:sim-mesh sim)
            ;;   (lambda (i) (cl-mpm/bc::make-bc-ambient-temp i
            ;;                                                0d0
            ;;                                                0d0))))
            ;; )
      (setf (cl-mpm:sim-bcs sim)
            (append
             (cl-mpm/bc::make-outside-bc-var
              (cl-mpm:sim-mesh sim)
              (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil)))
              (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil)))
              (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0)))
              (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0)))
             ;; (lambda (i) (cl-mpm/bc:make-bc-friction i
             ;; (magicl:from-list '(0d0 1d0) '(2 1)) 0.25d0))
              )
             ;; (cl-mpm/bc::make-domain-bcs
             ;;  (cl-mpm:sim-mesh sim)
             ;;  (lambda (i) (cl-mpm/bc::make-bc-ambient-temp i
             ;;                                               0d0
             ;;                                               0d0)))
             ;; (cl-mpm/bc::make-outside-bc-var
             ;;  (cl-mpm:sim-mesh sim)
             ;;  (lambda (i) (cl-mpm/bc:make-bc-fixed-temp i nil))
             ;;  (lambda (i) (cl-mpm/bc:make-bc-fixed-temp i nil))
             ;;  (lambda (i) (cl-mpm/bc:make-bc-fixed-temp i nil))
             ;;  (lambda (i) (cl-mpm/bc:make-bc-fixed-temp i 0d0)))
             )
            )
      ;; (setf (cl-mpm::sim-bcs-force sim)
      ;;       (list (cl-mpm/bc::make-bc-closure '(0 0)
      ;;                                         (lambda ()
      ;;                                           (cl-mpm/buoyancy::apply-bouyancy (cl-mpm:sim-mesh sim)
      ;;                                                                ;; (let ((ocean-x 0)
      ;;       (ocean-y 200))
      ;;   (setf (cl-mpm::sim-bcs-force sim)
      ;;         (loop for x from (floor ocean-x h) to (floor (first size) h)
      ;;               append (loop for y from 0 to (floor ocean-y h)
      ;;                            collect (cl-mpm/bc::make-bc-buoyancy
      ;;                                     (list x y)
      ;;                                     (magicl:from-list (list 0d0 (* 9.8d0 (* 1.0d0 1000))) '(2 1))))))
      ;;   )
                  ;; (cl-mpm:sim-mps sim))
      ;;                                           ))))

      (let ((ocean-x 0)
            (ocean-y 200))
        (setf (cl-mpm::sim-bcs-force sim)
              (loop for x from (floor ocean-x h) to (floor (first size) h)
                    append (loop for y from 0 to (floor ocean-y h)
                                 collect (cl-mpm/bc::make-bc-buoyancy
                                          (list x y)
                                          (magicl:from-list (list 0d0 (* 9.8d0 (* 1.0d0 1000))) '(2 1))))))
        )
      sim)))

;Setup
(defun setup ()
  (defparameter *sim* (setup-test-column '(600 400) '(500 100) '(000 100) (/ 1 25) 2))
  ;; (defparameter *sim* (setup-test-column '(1 1) '(1 1) '(0 0) 1 1))
  (remove-sdf *sim* (rectangle-sdf (list 500 200) '(100 50)))
  ;; (remove-sdf *sim* (ellipse-sdf (list 250 100) 20 40))
  ;; (remove-sdf *sim* (rectangle-sdf '(250 100) '(25 25)))
  ;; (remove-sdf *sim* (ellipse-sdf (list 1.5 3) 0.25 0.5))
  (print (cl-mpm:sim-dt *sim*))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  (defparameter *x* 0d0)
  (defparameter *x-pos* '())
  (defparameter *sim-step* 0)
  (defparameter *run-sim* nil)
  (defparameter *run-sim* t)
  (defparameter *load-mps*
    (let* ((mps (cl-mpm:sim-mps *sim*))
           (least-pos
              (apply #'max (loop for mp across mps
                                 collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))))
      (loop for mp across mps when (>= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (- least-pos 0.001))
              collect mp)))
  ;; (increase-load *sim* *load-mps* 1)
  ;; (increase-load *sim* *load-mps* 100)
  )

;; (defun sample-point-gimp (sim point get-value size)
;;   (let ((sample-mp (cl-mpm/particle:make-particle 2
;;                                                   :pos point
;;                                                   :size size)))
;;         (funcall get-value (cl-mpm:sim-mesh sim) sample-mp)))
;; (defun sample-point-volume (sim point size)
;;     (sample-point-gimp sim point size
;;         (scalar-average cl-mpm/mesh::node-volume)))
;; (defun get-mp-)
(defun add-snow-layer (sim);;; Printing methods

(defmethod print-object ((obj node) stream)
  (print-unreadable-object (obj stream :type t)
    (with-accessors ((index node-index))
        obj
      (format stream "index: ~a" (list (magicl:tref index 0 0)
                                       (magicl:tref index 1 0)
                                       )))))

  (let ((snow-height 5))
    (with-accessors ((mps cl-mpm:sim-mps))
        sim
      (loop for mp across mps
            collect (sample-point-volume
                     (.+ (cl-mpm/particle:mp-position)))))))
(defun add-particles (sim position size density)
  (vector-push-extend
                      (cl-mpm::make-particle
                       2
                       'cl-mpm/particle::particle-thermoviscoplastic-damage
                       :position (magicl:from-list position '(2 1))
                       :volume (reduce #'* size)
                       :size (magicl:from-list size '(2 1))
                       :E 1e7
                       :nu 0.325
                       :visc-factor 1d-10
                       :visc-power 3
                       :temperature 0d0
                       :heat-capacity 1d0
                       :thermal-conductivity 1e0
                       :mass (* (reduce #'* size) density)
                       :critical-stress 1d20
                       :gravity -9.8d0
                       :index 0
                       )

                      (cl-mpm:sim-mps sim)
                      )
  )

(defparameter *run-sim* nil)

(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)
  (defparameter *run-sim* t)
    ;; (vgplot:close-all-plots)
  ;; (sleep 1)
  (vgplot:figure)
  ;; (vgplot:plot '(10 10))
  (sleep 1)
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:axis (list 0 ms-x
                       0 ms-y))
    (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
    (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
      (vgplot:format-plot t "set ytics ~f" h)
      (vgplot:format-plot t "set xtics ~f" h))
  (with-open-file (stream (merge-pathnames "output/terminus_position.csv") :direction :output :if-exists :supersede)
    (format stream "Time (s),Terminus position~%")
    (loop for tim in (reverse *time*)
          for x in (reverse *x-pos*)
          do (format stream "~f, ~f ~%" tim x)))
  (defparameter *notch-position* 0.1d0)
  (time (loop for steps from 0 to 100
                while *run-sim*
                do
                (progn
                  (format t "Step ~d ~%" steps)
                  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)

                  (push *t* *time*)
                  (setf *x*
                        (loop for mp across (cl-mpm:sim-mps *sim*)
                              maximize (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0)))
                  (push
                   *x*
                   *x-pos*)
                  (let ((max-cfl 0))
                    (time (dotimes (i 1000)
                            ;; (increase-load *sim* *load-mps* (magicl:from-list (list (* (cl-mpm:sim-dt *sim*)
                                                                                       ;; 5d0) 0d0) '(2 1)))
                            ;; (pescribe-velocity *sim* *load-mps* '(1d0 nil))
                           (cl-mpm::update-sim *sim*)
                            ;; (remove-sdf *sim* (ellipse-sdf (list 200 225) *notch-position* 50))
                            ;; (incf *notch-position* (* (cl-mpm:sim-dt *sim*) 5d0))
                            ;; (cl-mpm/damage::calculate-damage (cl-mpm:sim-mesh *sim*)
                            ;;                                  (cl-mpm:sim-mps *sim*)
                            ;;                                  (cl-mpm:sim-dt *sim*)
                            ;;                                  50d0)
                           (setf *t* (+ *t* (cl-mpm::sim-dt *sim*))))))
                  (incf *sim-step*)
                  (plot *sim*)
                  (swank.live:update-swank)
                  (sleep .01)
                  (with-open-file (stream (merge-pathnames "output/terminus_position.csv") :direction :output :if-exists :append)
                          (format stream "~f, ~f ~%" *t* *x*)
                    )

                  )))
    (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*))
                                          *sim*)
    (vgplot:figure)
    (vgplot:title "Terminus over time")
    (vgplot:plot *time* *x-pos*)

  (with-open-file (stream (merge-pathnames "output/terminus_position.csv") :direction :output :if-exists :supersede)
    (format stream "Time (s),Terminus position~%")
    (loop for tim in (reverse *time*)
          for x in (reverse *x-pos*)
          do (format stream "~f, ~f ~%" tim x)))
  )

(setf lparallel:*kernel* (lparallel:make-kernel 8 :name "custom-kernel"))

(defun profile ()
  (sb-profile:unprofile)
  (sb-profile:reset)
  (sb-profile:profile cl-mpm::update-sim
                      cl-mpm::reset-grid
                      cl-mpm::p2g
                      cl-mpm::filter-grid
                      cl-mpm::update-nodes
                      cl-mpm::apply-bcs
                      cl-mpm::g2p
                      cl-mpm::update-particle
                      cl-mpm::update-stress
                      cl-mpm::iterate-over-neighbours-shape
                      cl-mpm::iterate-over-neighbours-shape-linear
                      cl-mpm::p2g-mp
                      cl-mpm::g2p-mp
                      cl-mpm::p2g-mp-node
                      cl-mpm::g2p-mp-node
                      ;; cl-mpm::update-strain-kirchoff
                      cl-mpm/damage::calculate-damage
                      cl-mpm/damage::apply-damage
                      cl-mpm/damage::delocalise-damage
                      cl-mpm/damage::create-delocalisation-list
                      ;; cl-mpm/eigenerosion:update-fracture
                      ;; cl-mpm/eigenerosion::remove-material-damaged
                      ;; cl-mpm/eigenerosion::find-neighbours
                      )
  (loop repeat 100
        do (progn
             (cl-mpm::update-sim *sim*)
             ;; (cl-mpm/damage::calculate-damage (cl-mpm:sim-mesh *sim*)
             ;;                                  (cl-mpm:sim-mps *sim*)
             ;;                                  (cl-mpm:sim-dt *sim*)
             ;;                                  25d0)
             ;; (cl-mpm/eigenerosion:update-fracture *sim*)
             ))
  (sb-profile:report))
(defun simple-time ()
  (setf lparallel:*kernel* (lparallel:make-kernel 8 :name "custom-kernel"))
  (let ((mp (cl-mpm/particle:make-particle 2 'cl-mpm/particle:particle
                                           :size (magicl:from-list '(1d0 1d0) '(2 1)))))
    (format t "Old")
      ;; (time
      ;;  (let ((m (magicl:zeros '(2 2))))
      ;;  (lparallel:pdotimes (i 100000)
      ;;    (magicl::eighth m)))
      ;; )
  (time
   ;; (lparallel:pdotimes (i 1000000)
   (dotimes (i 100)
     (with-accessors ((mesh cl-mpm:sim-mesh)
                      (mps cl-mpm:sim-mps)
                      (dt cl-mpm:sim-dt))
         *sim*
       (cl-mpm::update-sim *sim*)
       ;; (cl-mpm::update-stress-mp mesh (aref mps 0) dt)
         )
     ))
    (format t "new")
    ))
;(let ((stress (magicl:zeros '(3 1)))
;      (strain (magicl:zeros '(3 1)))
;      (de (cl-mpm/constitutive::linear-elastic-matrix 1d0 0d0))
;      )
;  (format t "mat")
;  (time
;   (dotimes (i 100000)
;       (cl-mpm/constitutive::maxwell-exp strain stress 1d0 0d0 1d0 1d0)))
;  (format t "voight")
;  (time
;   (dotimes (i 100000)
;     (cl-mpm/constitutive::maxwell-exp-v strain stress 1d0 0d0 1d0 1d0)))
;  (format t "simd")
;  (time
;   (dotimes (i 100000)
;     (cl-mpm/constitutive::maxwell-exp-v-simd strain stress 1d0 0d0 de 1d0  1d0)))
;  )
