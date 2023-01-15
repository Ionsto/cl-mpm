(defpackage :cl-mpm/examples/beam
  (:use :cl))
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)

(in-package :cl-mpm/examples/beam)
(declaim (optimize (debug 2) (safety 2) (speed 2)))


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
               (setf (cl-mpm/particle:mp-velocity mp) vel)))))
(defun increase-load (sim load-mps amount)
  (loop for mp in load-mps
        do (with-accessors ((pos cl-mpm/particle:mp-position)
                            (force cl-mpm/particle:mp-body-force)) mp
             (incf (magicl:tref force 0 0) amount))))

(defun length-from-def (sim mp dim)
  (let* ((mp-scale 2)
         (h-initial (magicl:tref (cl-mpm/particle::mp-domain-size mp) dim 0))
         ;(h-initial  (/ (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) mp-scale))
         )
    h-initial))
(defun max-stress (mp)
  (multiple-value-bind (l v) (magicl:eig (cl-mpm::voight-to-matrix (cl-mpm/particle:mp-stress mp)))
    ;; (apply #'max l)
    (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0)
    ))
(defun plot (sim &optional (plot :deformed))
  (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  (multiple-value-bind (x y c stress-y lx ly e density)
    (loop for mp across (cl-mpm:sim-mps sim)
          collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
          collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
          collect (length-from-def sim mp 0) into lx
          collect (length-from-def sim mp 1) into ly
          collect (if (slot-exists-p mp 'cl-mpm/particle::damage) (cl-mpm/particle:mp-damage mp) 0) into c
          collect (if (slot-exists-p mp 'cl-mpm/particle::strain-energy-density) (cl-mpm/particle::mp-strain-energy-density mp) 0) into e
          collect (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)) into density
          ;; collect (cl-mpm/particle:mp-volume mp) into density
          collect (max-stress mp) into stress-y
          finally (return (values x y c stress-y lx ly e density)))
    (cond
      ((eq plot :damage)
       (vgplot:format-plot t "set cbrange [0:1]")
       ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 0.01 (apply #'max c)))
       (vgplot:plot x y c ";;with points pt 7 lc palette")
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
(defun ellipse-sdf (position x-l y-l)
  (let ((aspect (/ x-l y-l)))
    (lambda (pos)
      (let* ((position (magicl:from-list position '(2 1) :type 'double-float))
             (dist-vec (magicl:.* (magicl:.- position pos) (magicl:from-list (list 1 aspect) '(2 1)
                                                                             :type 'double-float)))
             (distance (sqrt (magicl:tref (magicl:@ (magicl:transpose dist-vec)
                                                    dist-vec) 0 0))))
        (- distance x-l)))))

(defun setup-test-column (size block-size block-offset &optional (e-scale 1d0) (mp-scale 1d0))
  (let* ((sim (cl-mpm/setup::make-block (/ 1 e-scale)
                                        (mapcar (lambda (s) (* s e-scale)) size)
                                        #'cl-mpm/shape-function:make-shape-function-linear)) 
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         ;(e-scale 1)
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 100)
         ;(mass (/ (* 100 h-x h-y) (expt mp-scale 2)))
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
               ;; 'cl-mpm/particle::particle-viscoelastic-fracture

               'cl-mpm/particle::particle-elastic-damage
               :E 1d8
               :nu 0.325d0

               ;; :E 1e6 :nu 1d8

               ;; 'cl-mpm/particle::particle-viscoplastic
               ;; :E 1e7
               ;; :nu 0.325
               ;; :visc-factor 1e6
               ;; :visc-power 3

               ;; Fluid
               ;; 'cl-mpm/particle::particle-fluid
               ;; :stiffness 1e5
               ;; :adiabatic-index 6
               ;; :rest-density 900d0
               ;; :viscosity 1.02e-3

               ;; :E 1e6 :nu 0.33
               :critical-stress 1d10
               ;; :fracture-toughness 5d0
               :gravity -9.8d0
               )))
      (let ((prev-mps (cl-mpm:sim-mps sim)))
        (setf (cl-mpm:sim-mps sim) (make-array 1000 :fill-pointer 0 :adjustable t))
        (loop for mp across prev-mps
              do (vector-push-extend mp (cl-mpm:sim-mps sim))))
      (setf (cl-mpm:sim-damping-factor sim) 0.5d0)
      (setf (cl-mpm:sim-mass-filter sim) 1d-10)
      (setf (cl-mpm:sim-dt sim) 1d-2)

      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-var (cl-mpm:sim-mesh sim)
                                            (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 0)))
                                            (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil)))
                                            (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0)))
                                            (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0)))
                                            ;; (lambda (i) (cl-mpm/bc:make-bc-friction i
                                            ;;                                         (magicl:from-list '(0d0 1d0)
                                            ;;                                                           '(2 1))
                                            ;;                                         0.25d0))
                                            ))

      ;; (let ((step-x 200)
      ;;       (step-y 400)
      ;;       )

      ;;   (loop for x from 0 to (round step-x h)
      ;;         do (loop for y from (- (round step-y h) 1) to (round step-y h)
      ;;                  do (push
      ;;                      (cl-mpm/bc:make-bc-fixed (list x y) '(nil 0))
      ;;                           (cl-mpm:sim-bcs sim))))
      ;;   )
      sim)))

;Setup
(defun setup ()
  (declare (optimize (speed 0)))
  (defparameter *sim* (setup-test-column '(700 600) '(500 100) '(000 400) (/ 1 50) 2))
  ;; (defparameter *sim* (setup-test-column '(1 1) '(1 1) '(0 0) 1 1))
  ;;(remove-sdf *sim* (ellipse-sdf (list 400 100) 10 40))
  ;; (remove-sdf *sim* (ellipse-sdf (list 1.5 3) 0.25 0.5))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  (defparameter *sim-step* 0)
  (defparameter *run-sim* nil)
  (defparameter *run-sim* t)
  (defparameter *load-mps*
    (let* ((mps (cl-mpm:sim-mps *sim*))
           (least-pos
              (apply #'max (loop for mp across mps
                                 collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))))
      (loop for id from 0 to (- (length mps) 1) when (>= (magicl:tref (cl-mpm/particle:mp-position (aref mps id)) 0 0) (- least-pos 0.001))
              collect (aref mps id))))
  (defparameter *deflection-mps*
    (let* ((mps (cl-mpm:sim-mps *sim*))
           (least-pos
             (apply #'max (loop for mp across mps
                                collect (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)))))
      (loop for id from 0 to (- (length mps) 1) when (>= (magicl:tref (cl-mpm/particle:mp-position (aref mps id)) 1 0) (- least-pos 0.001))
            collect (aref mps id))))
  ;; (increase-load *sim* *load-mps* 1)
  ;; (increase-load *sim* *load-mps* 100)
  )

(defun remove-material-damaged (sim)
  "Remove material points that have strain energy density exceeding fracture toughness"
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (delete-if (lambda (mp)
                 (with-accessors ((damage cl-mpm/particle:mp-damage))
                     mp
                   (>= damage 1d0))) mps)))
(defun mag-sq (a)
  (magicl::sum (magicl:.* a a)))

(defun calculate-energy-gravity (sim)
  (loop for mp across (cl-mpm:sim-mps sim)
        sum (* 9.8 (magicl:tref (cl-mpm/particle::mp-displacement mp) 1 0) (cl-mpm/particle:mp-mass mp))))

(defun calculate-energy-kinetic (sim)
  (loop for mp across (cl-mpm:sim-mps sim)
        sum (* 0.5 (cl-mpm/particle:mp-mass mp) (mag-sq (cl-mpm/particle:mp-velocity mp)))))

(defun calculate-energy-strain (sim)
  (loop for mp across (cl-mpm:sim-mps sim)
        sum
        (with-accessors ((volume cl-mpm/particle:mp-volume)
                         (stress cl-mpm/particle:mp-stress)
                         (strain cl-mpm/particle:mp-strain)
                         (def cl-mpm/particle::mp-deformation-gradient)
                         )
            mp
          (* 0.5 volume (magicl::sum (magicl:.* strain (magicl:scale stress (magicl:det def)))))
          ;; (* 0.5 volume (magicl::sum (magicl:.* strain (magicl:scale stress 1d0))))
          )))

(defparameter *run-sim* nil)

(defun run ()
  (declare (optimize (speed 1) (debug 2)))
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)
  (defparameter *run-sim* t)
    (vgplot:close-all-plots)
    (vgplot:figure)
  (format t "MP count ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (defparameter *energy-gpe* '())
  (defparameter *energy-ke* '())
  (defparameter *energy-se* '())
  (sleep 1)
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:axis (list 0 ms-x
                       0 ms-y))
    (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
    ;; (vgplot:axis (list 0 (nth 0 (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*))) 
    ;;                    0 (nth 1 (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))))
    (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
      ;; (vgplot:format-plot t "set ytics ~f" h)
      (vgplot:format-plot t "set xtics ~f" h)
      (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*))
                              *sim*)
      (incf *sim-step*)
      (time (loop for steps from 0 to 100
                            while *run-sim*
                            do
                               (progn
                                 (format t "Step ~d ~%" steps)
                  ;;; Adaptive timestepping
                                 ;; (let* ((target-dt 1e-1)
                                 ;;        (courent (find-max-cfl *sim*))
                                 ;;        (c-value 0.001d0)
                                 ;;        (cfl-dt (if (> courent 0d0) (/ c-value (/ courent (cl-mpm:sim-dt *sim*))) nil))
                                 ;;        (new-dt (/ c-value (if (> courent 0d0)
                                 ;;                               (/ courent (cl-mpm:sim-dt *sim*))
                                 ;;                               (/ c-value 1e-6))))
                                 ;;        (sub-steps (max (min (floor (/ target-dt new-dt)) 100) 1)))
                                 ;;   (format t "C: ~f - steps: ~D - %dt: ~f~%" courent sub-steps new-dt)
                                 ;;   (format t "Cfl derived dt:~f~%" cfl-dt)
                                 ;;   (setf (cl-mpm:sim-dt *sim*) new-dt))
                                 ;; (break)
                                 (let ((max-cfl 0))
                                   (time (loop for i from 0 to 100
                                               while *run-sim*
                                               do
                                                  (progn

                                                    ;; (pescribe-velocity *sim* *load-mps* (magicl:from-list '(0.5d0 0d0) '(2 1)))
                                                    (cl-mpm::update-sim *sim*)
                                                    ;; (cl-mpm/damage::calculate-damage (cl-mpm:sim-mesh *sim*)
                                                    ;;                               (cl-mpm:sim-mps *sim*)
                                                    ;;                               (cl-mpm:sim-dt *sim*)
                                                    ;;                               25d0
                                                    ;;                               )
                                                    (setf *t* (+ *t* (cl-mpm::sim-dt *sim*))))
                                           ))
                                   ;; (format t "Max cfl: ~f~%" max-cfl)
                                   )
                                 (push *t* *time*)
                                 (push (calculate-energy-strain *sim*) *energy-se*)
                                 (push (calculate-energy-kinetic *sim*) *energy-ke*)
                                 (push (calculate-energy-gravity *sim*) *energy-gpe*)
                                 (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*))
                                                         *sim*)
                                 (incf *sim-step*)
                                 (plot *sim*)
                                 ;; (vgplot:print-plot (asdf:system-relative-pathname "cl-mpm" (format nil "output/frame_~5,'0d.png" steps)))
                                 (swank.live:update-swank)
                                 (sleep .01)

                                 )))
      (plot-deflection (cl-mpm/particle::mp-e (first *deflection-mps*))
                       100 100 500)
      (with-open-file (stream (merge-pathnames "output/energy.csv") :direction :output :if-exists :supersede)
        (format stream "Time (s),SE,KE,GPE~%")
        (loop for tim in (reverse *time*)
              for se in (reverse *energy-se*)
              for ke in (reverse *energy-ke*)
              for gpe in (reverse *energy-gpe*)
              do (format stream "~f, ~f, ~f, ~f ~%" tim se ke gpe)))
    ;; (vgplot:figure)
    ;; (vgplot:title "Velocity over time")
    ;; (vgplot:plot *time* *velocity*)
      (plot-energy)
    ))
(defun plot-energy ()
  (vgplot:figure)
  (vgplot:title "Velocity over time")
  (vgplot:plot *time* *energy-gpe* "GPE"
               *time* *energy-ke* "KE"
               *time* *energy-SE* "SE"
               *time* (mapcar #'+ *energy-gpe* *energy-ke* *energy-se*) "total energy"
               ))

(defmacro time-form (form it)
  `(progn
     (declaim (optimize speed))
     (let* ((iterations ,it)
            (start (get-internal-real-time)))
       (dotimes (i ,it)
         ,form)
       (let* ((end (get-internal-real-time))
              (units internal-time-units-per-second)
              (dt (/ (- end start) (* iterations units)))
              )
         (format t "Total time: ~f ~%" (/ (- end start) units)) (format t "Time per iteration: ~f~%" (/ (- end start) (* iterations units)))
         (format t "Throughput: ~f~%" (/ 1 dt))
         dt))))

(defun plot-deflection (E h density l)
  (let* (
        (x (loop for mp in *deflection-mps*
                 collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
        (ymin (magicl:tref (cl-mpm/particle:mp-position (first *deflection-mps*)) 1 0))
        (y (loop for mp in *deflection-mps*
                 collect (- (magicl:tref (cl-mpm/particle:mp-position mp) 1 0) ymin)))
        (f (/ (* -12 density 9.8) (* 24 E (* h h))))
        (d-udl (mapcar (lambda (k) (* f (expt k 2) (+ (expt k 2) (* 6 l l) (* -4 l k)))) x))
        )
    (vgplot:figure)
    (vgplot:plot x y "Experimental"
                 x d-udl "Analyticl"
                 )
    )
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
                      cl-mpm::calculate-strain-rate
                      cl-mpm::update-strain-kirchoff
                      cl-mpm::iterate-over-neighbours-shape
                      cl-mpm::iterate-over-neighbours-shape-linear
                      cl-mpm::p2g-mp
                      cl-mpm::g2p-mp
                      cl-mpm::p2g-mp-node
                      cl-mpm::g2p-mp-node
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
