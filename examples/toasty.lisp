(defpackage :cl-mpm/examples/slump
  (:use :cl))
(sb-ext:restrict-compiler-policy 'speed 3 3)
(sb-ext:restrict-compiler-policy 'debug 0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(in-package :cl-mpm/examples/slump)
(declaim (optimize (debug 0) (safety 0) (speed 3)))


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
    (* (magicl:tref (cl-mpm::mp-deformation-gradient mp) dim dim) h-initial)))
(defun max-stress (mp)
  (multiple-value-bind (l v) (magicl:hermitian-eig (cl-mpm::voight-to-matrix (cl-mpm/particle:mp-stress mp)))
    (apply #'max l)
    (magicl:tref (cl-mpm/particle:mp-stress mp) 2 0)))
(defun plot (sim &optional (plot :temperature))
  (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  (multiple-value-bind (x y c stress-y lx ly e density temp)
    (loop for mp across (cl-mpm:sim-mps sim)
          collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
          collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
          collect (length-from-def sim mp 0) into lx
          collect (length-from-def sim mp 1) into ly
          collect (if (slot-exists-p mp 'cl-mpm/particle::damage) (cl-mpm/particle:mp-damage mp) 0) into c
          collect (if (slot-exists-p mp 'cl-mpm/particle::temperature) (cl-mpm/particle:mp-temperature mp) 0) into temp
          collect (if (slot-exists-p mp 'cl-mpm/particle::strain-energy-density) (cl-mpm/particle::mp-strain-energy-density mp) 0) into e
          collect (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)) into density
          ;; collect (cl-mpm/particle:mp-volume mp) into density
          collect (max-stress mp) into stress-y
          finally (return (values x y c stress-y lx ly e density temp)))
    (cond
      ((eq plot :damage)
       (vgplot:format-plot t "set cbrange [0:1]")
       ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 0.01 (apply #'max c)))
       (vgplot:plot x y c ";;with points pt 7 lc palette")
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
       :nu 0.325
       :visc-factor 1d5
       :visc-power 3
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
(defun setup-test-column (size block-size block-offset &optional (e-scale 1d0) (mp-scale 1d0))
  (let* ((sim (cl-mpm/setup::make-block (/ 1 e-scale)
                                        (mapcar (lambda (s) (* s e-scale)) size)
                                        #'cl-mpm/shape-function:make-shape-function-bspline)) 
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         ;(e-scale 1)
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (mass (/ (* 9 h-x h-y) (expt mp-scale 2)))
         (mass-bread (/ (* 100 h-x h-y) (expt mp-scale 2)))
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size)))
    (progn
      (let ((block-position
              (mapcar #'+ (list (* h-x (- (+ (/ 1 (* 2 mp-scale))) 0))
                                (* h-y (+ (/ 1d0 (* 2d0 mp-scale)))))
                      block-offset)))
        ;; (setf (cl-mpm:sim-mps sim)
        ;;       (cl-mpm/setup:make-block-mps
        ;;        block-position
        ;;        block-size
        ;;        (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
        ;;        'cl-mpm::make-particle
        ;;        'cl-mpm/particle::particle-thermoviscoplastic-damage
        ;;          :E 1e7
        ;;          :nu 0.325
        ;;          :visc-factor 1d0
        ;;          :visc-power 3
        ;;          :temperature 0d0
        ;;          :heat-capacity 1d0
        ;;          :thermal-conductivity 1e0
        ;;          :mass mass
        ;;          :critical-stress 1d20
        ;;          :gravity -9.8d0
        ;;          :index 0
        ;;        )
        ;;       )
        (setf (cl-mpm:sim-mps sim)
              (make-mps-bread 9d0
                              90d0
                              '(100 100)
                              '(300 100)
                              e-scale
                              mp-scale
                              h)
              )
        )
      (setf (cl-mpm:sim-damping-factor sim) 1d-9)
      (setf (cl-mpm:sim-mass-filter sim) 1d-5)
      (setf (cl-mpm:sim-dt sim) 1d-3)
             ;; (lambda (i) (cl-mpm/bc:make-bc-friction i
             ;; (magicl:from-list '(0d0 1d0) '(2 1)) 0.25d0))
      (setf (cl-mpm::sim-bcs-force sim)
            (append
             (cl-mpm/bc::make-domain-bcs
              (cl-mpm:sim-mesh sim)
              (lambda (i) (cl-mpm/bc::make-bc-ambient-temp i
                                                           0d0
                                                           0d0))))
            )
      (setf (cl-mpm:sim-bcs sim)
            (append
             (cl-mpm/bc::make-outside-bc-var
              (cl-mpm:sim-mesh sim)
              (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil)))
              (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil)))
              (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0)))
              (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0)))
             ;; (lambda (i) (cl-mpm/bc:make-bc-friction i
             ;; (magicl:from-list '(0d0 1d0) '(2 1)) 0.2d0))
              )
             ;; (cl-mpm/bc::make-domain-bcs
             ;;  (cl-mpm:sim-mesh sim)
             ;;  (lambda (i) (cl-mpm/bc::make-bc-ambient-temp i
             ;;                                               0d0
             ;;                                               0d0)))
             (cl-mpm/bc::make-outside-bc-var
              (cl-mpm:sim-mesh sim)
              (lambda (i) (cl-mpm/bc:make-bc-fixed-temp i nil))
              (lambda (i) (cl-mpm/bc:make-bc-fixed-temp i nil))
              (lambda (i) (cl-mpm/bc:make-bc-fixed-temp i nil))
              (lambda (i) (cl-mpm/bc:make-bc-fixed-temp i 10d0)))
             )
            )
      sim)))

;Setup
(defun setup ()
  (defparameter *sim* (setup-test-column '(600 400) '(100 100) '(0 0) (/ 1 20) 2))
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
  ;; (increase-load *sim* *load-mps* 1)
  ;; (increase-load *sim* *load-mps* 100)
  )


(defparameter *run-sim* nil)

(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)
  (defparameter *run-sim* t)
    (vgplot:close-all-plots)
    (vgplot:figure)
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
    (time (loop for steps from 0 to 100
                while *run-sim*
                do
                (progn
                  (format t "Step ~d ~%" steps)
                  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*))
                                          *sim*)
                  (let ((max-cfl 0))
                    (time (dotimes (i 1000)
                           (cl-mpm::update-sim *sim*)
                           (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))
                           ))
                    )
                  (incf *sim-step*)
                  (plot *sim*)
                  (swank.live:update-swank)
                  (sleep .01)

                  )))
    (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*))
                                          *sim*)
    ;; (vgplot:figure)
    ;; (vgplot:title "Velocity over time")
    ;; (vgplot:plot *time* *velocity*)
    )

(setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))
