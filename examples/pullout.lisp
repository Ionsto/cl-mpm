(defpackage :cl-mpm/examples/pullout
  (:use :cl))
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)
(in-package :cl-mpm/examples/pullout)
;; (pushnew :cl-mpm-pic *features*)
(delete :cl-mpm-pic *features*)
;; (asdf:compile-system :cl-mpm :force T)

(defun max-v-sum (mp)
  (with-accessors ((vel cl-mpm/particle:mp-velocity))
      mp
    (magicl::sum (magicl:map #'abs vel))))

(defun find-max-cfl (sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps)
                   (dt cl-mpm:sim-dt))
      sim
    (let ((max-v (lparallel:preduce #'max
                                    (lparallel:pmapcar #'max-v-sum
                                                    mps))))
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
                            (vol cl-mpm/particle::mp-volume)
                            (force cl-mpm/particle:mp-body-force)) mp
             ;(incf (magicl:tref force 0 0) amount)
             (magicl:.+ force (magicl:scale amount (/ 1d0 vol)) force)
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
    (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0)))
(defun plot (sim &optional (plot :damage))
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
      ((eq plot :damage)
       ;; (vgplot:format-plot t "set cbrange [0:1]")
       ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 1e-6 (apply #'max c)))
       (vgplot:format-plot t "set cbrange [~f:~f]" 0d0 (+ 1e-6 (apply #'max c)))
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

(defun damage-sdf (sim sdf &optional (new-damage 1d0))
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
      (loop for mp across mps
            do (with-accessors ((pos cl-mpm/particle:mp-position)
                                 (damage cl-mpm/particle:mp-damage)) mp
                  (when (>= 0 (funcall sdf pos))
                    (setf damage new-damage))))))

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
               'cl-mpm/particle::particle-elastic-damage
               ;; 'cl-mpm/particle::particle-viscoplastic
                ;; 'cl-mpm/particle::particle-viscoplastic-damage
               :E 1d8
               :nu 0.3250d0
               ;; ;; 'cl-mpm/particle::particle-elastic-damage
               ;; :E 1d9
               ;; :nu 0.3250d0

               ;; :visc-factor 11d6
               ;; :visc-power 3d0

               :initiation-stress 0.2d6
               ;; :damage-rate 1d-9
               ;; :damage-rate 1d-8
               :damage-rate 1d-9
               :critical-damage 1.0d0
               :local-length 1000d0
               :gravity 0d0;-9.8d0

                 ;; :gravity-axis (magicl:from-list '(0.5d0 0.5d0) '(2 1))
                 :index 0
               )))
      (setf (cl-mpm:sim-damping-factor sim) 1.5d0)
      (setf (cl-mpm:sim-mass-filter sim) 1d-15)
      (setf (cl-mpm::sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) t)
      (setf (cl-mpm:sim-dt sim) 1d-1)
      ;; (setf (cl-mpm::sim-bcs-force sim)
      ;;       (cl-mpm/bc:make-bcs-from-list
      ;;        (list
      ;;         (cl-mpm/bc::make-bc-closure '(0 0)
      ;;                                     (lambda ()
      ;;                                       (cl-mpm/buoyancy::apply-bouyancy sim 200d0))))))
      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-var
             (cl-mpm:sim-mesh sim)
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0)))
             ;; (lambda (i) (cl-mpm/bc:make-bc-friction i
             ;; (magicl:from-list '(0d0 1d0) '(2 1)) 0.25d0))
             )
            )
      (defparameter *pressure-inc-rate* 0d4)
      (defparameter *shear-rate* 0.1d0)
      ;; (setf (cl-mpm:sim-bcs sim)
      ;;       (cl-mpm/bc:make-bcs-from-list
      ;;        (append
      ;;         (map 'list #'identity (cl-mpm:sim-bcs sim))
      ;;         (list
      ;;          (cl-mpm/bc::make-bc-closure
      ;;           '(0 0)
      ;;           (lambda ()
      ;;             (apply-pullout
      ;;              sim
      ;;              *terminus-mps*
      ;;              *shear-rate*)
      ;;             )
      ;;           )))))
      (defparameter *load-bc*
        (cl-mpm/buoyancy::make-bc-pressure
         sim
         4d5
         0d0
         ))
      (setf (cl-mpm::sim-bcs-force sim)
            (cl-mpm/bc:make-bcs-from-list
             (list *load-bc*)))

      sim)))

;Setup
(defun setup ()
  (defparameter *run-sim* nil)
  (let ((mesh-size 50)
        (mps-per-cell 2))
    ;;Setup notched pullout
    ;; (defparameter *sim* (setup-test-column '(1000 300) '(500 125) '(000 0) (/ 1 mesh-size) mps-per-cell))
    ;; (remove-sdf *sim* (rectangle-sdf '(250 125) '(10 10)))
    ;;Setup 1d pullout
    (defparameter *sim* (setup-test-column (list 1000 mesh-size)
                                           (list 500 mesh-size) '(000 0) (/ 1 mesh-size) mps-per-cell))
    ;; (damage-sdf *sim* (rectangle-sdf '(250 0) (list
    ;;                                            (* 1 mesh-size)
    ;;                                            mesh-size)) 0.1d0)
    )
  ;; (remove-sdf *sim* (ellipse-sdf (list 1.5 3) 0.25 0.5))
  (print (cl-mpm:sim-dt *sim*))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  (defparameter *x* 0d0)
  (defparameter *x-pos* '())

  (defparameter *max-stress* '())
  (defparameter *max-damage* '())
  (defparameter *max-x* '())

  (defparameter *cfl-max* '())
  (defparameter *sim-step* 0)
  ;(defparameter *run-sim* t)
  (defparameter *load-mps*
    (let* ((mps (cl-mpm:sim-mps *sim*))
           (least-pos
              (apply #'max (loop for mp across mps
                                 collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))))
      (loop for mp across mps when (>= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (- least-pos 0.001))
              collect mp)))
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (let ((x-min (loop for mp across mps
                       minimize (magicl:tref
                                 (cl-mpm/particle:mp-position mp)
                                 0 0))))
      (defparameter *far-field-mps*
        (loop for mp across mps
              when (= x-min (magicl:tref
                         (cl-mpm/particle:mp-position mp)
                         0 0))
                collect mp)))
    (let ((y-min (loop for mp across mps
                       minimize (magicl:tref
                                 (cl-mpm/particle:mp-position mp)
                                 1 0))))
      (defparameter *bottom-mps*
        (loop for mp across mps
              when (= y-min (magicl:tref
                             (cl-mpm/particle:mp-position mp)
                             1 0))
                collect mp)))
    (let ((x-max (loop for mp across mps
                       maximize (magicl:tref
                                 (cl-mpm/particle:mp-position mp)
                                 0 0))))
      (defparameter *terminus-mps*
        (loop for mp across mps
              when (= x-max (magicl:tref
                             (cl-mpm/particle:mp-position mp)
                             0 0))
                collect mp)))
    ;; (increase-load *sim* *terminus-mps*
    ;;                (magicl:from-list (list (* 1d5) 0d0) '(2 1)))
    )
  ;; (increase-load *sim* *load-mps* 1)
  ;; (increase-load *sim* *load-mps* 100)
  )


(defparameter *run-sim* nil)
(defun calculate-dt (courent c-value target-step)
  (let* (
         (cfl-dt (if (> courent 0d0) (/ c-value (/ courent (cl-mpm:sim-dt *sim*))) nil))
         (new-dt (/ c-value (if (> courent 0d0)
                                (/ courent (cl-mpm:sim-dt *sim*))
                                (/ c-value 1e-6))))
         (max-steps 1000)
         (sub-steps (max (min (floor (/ target-step new-dt)) max-steps) 1)))
    (when (> (floor (/ target-step new-dt)) max-steps)
        (format t "CFL requires more steps than max-steps~%"))
    (format t "C: ~f - steps: ~D - %dt: ~f~%" courent sub-steps new-dt)
    (format t "Cfl derived dt:~f~%" cfl-dt)
                                        ;(setf (cl-mpm:sim-dt *sim*) new-dt)
    (values new-dt sub-steps)))


(defun apply-pullout (sim load-mps shear-rate)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (loop for mp in load-mps
          do
             (cl-mpm::iterate-over-neighbours
              mesh
              mp
              (lambda (mesh mp node svp grad)
                (with-accessors (
                                 ;; (vel cl-mpm/particle:mp-velocity)
                                 )
                    mp
                  (with-accessors ((pos cl-mpm/mesh::node-position)
                                   (vel cl-mpm/mesh::node-velocity)
                                   (active cl-mpm/mesh::node-active)
                                   (acc cl-mpm/mesh::node-acceleration)
                                   )
                      node
                    (when active
                      (setf (magicl:tref vel 0 0)
                            shear-rate
                            (magicl:tref vel 1 0)
                            0d0
                            )
                      ;; (setf (magicl:tref acc 0 0) 0d0
                      ;;       (magicl:tref acc 1 0) 0d0
                      ;;       )
                      )
                    ))))
             ;; (progn
             ;;   (setf (cl-mpm/particle:mp-velocity mp) vel))
          )))

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
  (with-open-file (stream (merge-pathnames "output/terminus_position.csv") :direction :output :if-exists :supersede)
    (format stream "Time (s),Terminus position~%")
    (loop for tim in (reverse *time*)
          for x in (reverse *x-pos*)
          do (format stream "~f, ~f ~%" tim x)))
  ;; (dotimes (i 1000)
  ;;   (cl-mpm::update-sim *sim*))

  (let* ((target-time 10d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt)))
    (format t "Substeps ~D~%" substeps)
    (time (loop for steps from 0 to 500
                while *run-sim*
                do
                   (progn
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     ;(cl-mpm/output:save-csv (merge-pathnames (format nil "output/sim_csv_~5,'0d.vtk" *sim-step*)) *sim*)

                     (push *t* *time*)
                     ;; (let ((cfl (find-max-cfl *sim*)))
                     ;;   (push cfl *cfl-max*))
                     (with-accessors ((mps cl-mpm:sim-mps))
                         *sim*
                       (setf *x*
                             (loop for mp across mps
                                   maximize (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0)))
                       (push
                        *x*
                        *x-pos*)

                       (push
                        (loop for mp across mps
                              maximize (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0))
                        *max-x*)
                       (push
                        (/ (loop for mp across mps
                                 sum (magicl:tref (cl-mpm/particle::mp-stress mp) 0 0))
                           (length mps))
                        *max-stress*)
                       (push
                        (/ (loop for mp across mps
                                 sum (cl-mpm/particle::mp-damage mp))
                           (length mps))
                        *max-damage*)
                       ;; (let ((m-d (loop for mp across mps maximize (cl-mpm/particle::mp-damage mp))))
                       ;;   (when (>= m-d 1d0)
                       ;;     (setf *run-sim* nil)))
                       )

                     (let ((cfl 0))
                       (time (dotimes (i substeps)
                               (incf (first (cl-mpm/buoyancy::bc-pressure-pressures *load-bc*))
                                     (* (cl-mpm:sim-dt *sim*) *pressure-inc-rate*))
                               (cl-mpm::update-sim *sim*)
                               (setf cfl (max cfl (find-max-cfl *sim*)))
                               (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))
                       ;; (setf cfl (find-max-cfl *sim*))
                       (format t "CFL: ~f~%" cfl)
                       (push cfl *cfl-max*)
                         (multiple-value-bind (dt-e substeps-e) (calculate-dt cfl 1d-3 target-time)
                           (format t "CFL dt estimate: ~f~%" dt-e)
                           (format t "CFL step count estimate: ~D~%" substeps-e)
                           (push cfl *cfl-max*)
                           ;; (setf (cl-mpm:sim-dt *sim*) dt-e)
                           ;; (setf substeps substeps-e)
                           )
                         )
                     (with-accessors ((mps cl-mpm:sim-mps))
                         *sim*
                         (let ((m-d (loop for mp across mps maximize
                                                            (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0))))
                           (when (>= m-d 100)
                             (setf *run-sim* nil))))
                     ;; (setf (cl-mpm:sim-damping-factor *sim*) (max 0.1d0 (/ (cl-mpm:sim-damping-factor *sim*) 1.1d0)))
                     ;; (setf (cl-mpm::sim-enable-damage *sim*) t)
                     (incf *sim-step*)
                     (plot *sim*)
                     (swank.live:update-swank)
                     (sleep .01)
                     (with-open-file (stream (merge-pathnames "output/terminus_position.csv") :direction :output :if-exists :append)
                       (format stream "~f, ~f ~%" *t* *x*)
                       )

                     ))))
  (plot-stress-damage-time)
  )
(defun plot-damage ()
    (let ((x (loop for mp in *bottom-mps*
                   collect (magicl:tref
                            (cl-mpm/particle:mp-position mp)
                            0 0)))
          (d (loop for mp in *bottom-mps*
                   collect (cl-mpm/particle:mp-damage mp)))
          )
      (vgplot:figure)
      (vgplot:plot
       x d "damage"
       )
    ))
(defun plot-stress-damage-time ()
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (let ((s-y 0.2d6)
          (s-x (first *max-x*))
          (len 500d0))
      (vgplot:figure)
      (vgplot:plot
       *time* (mapcar (lambda (x) (/ x s-x)) *max-x*) "Far field extension"
       *time* (mapcar (lambda (x) (/ x s-y)) *max-stress*) "max stress"
       *time* *max-damage* "max damage"
       )
      (vgplot:figure)
      (vgplot:plot
       (mapcar (lambda (x) (/ x len)) *max-x*) (mapcar (lambda (x) (/ x s-y)) *max-stress*) "max stress"
       (mapcar (lambda (x) (/ x len)) *max-x*) *max-damage* "max damage"
       )
      )
    ))
(defun plot-s-xx (mps H)
  (let* (
         (rho-ice 900)
         (E (cl-mpm/particle::mp-e (first mps)))
         (nu (cl-mpm/particle::mp-nu (first mps)))
         (g (cl-mpm/particle::mp-gravity (first mps)))
         (s-xx (loop for mp in mps
                    collect
                    (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0)))
         (y (loop for mp in mps
                    collect
                    (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)))
         (s-an (mapcar
                (lambda (y)
                  (* (/ nu (- 1 nu)) rho-ice g -1d0 (- y (/ H 2d0))))
                y)))
    (vgplot:figure)
    (vgplot:title "S_{xx} over height")
    (vgplot:xlabel "Height (m)")
    (vgplot:ylabel "Longitudinal stress (MPa)")
    (vgplot:plot y (mapcar (lambda (x) (* 1d-6 x)) s-xx) "s_{xx} mpm"
                 y (mapcar (lambda (x) (* 1d-6 x)) s-an) "s_{xx} analytic")))
(defun plot-cfl ()
  (vgplot:figure)
  (vgplot:title "CFL over time")
  (vgplot:plot *time* *cfl-max*))

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

(defmacro time-form (it form)
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
(defun simple-time ()
  (setf lparallel:*kernel* (lparallel:make-kernel 8 :name "custom-kernel"))
  (let ((mp (cl-mpm/particle:make-particle 2 'cl-mpm/particle:particle
                                           :size (magicl:from-list '(1d0 1d0) '(2 1))))
        (mesh (cl-mpm::sim-mesh *sim*)))
  (time-form 1000
     (with-accessors ((mesh cl-mpm:sim-mesh)
                      (mps cl-mpm:sim-mps)
                      (dt cl-mpm:sim-dt))
         *sim*
       (cl-mpm::update-sim *sim*)
         ))))
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

;; (defun dot (x)
;;   (sqrt (magicl::sum (magicl:.* x x))))
;; (defun mp-sdf (mp x)
;;   (with-accessors ((pos cl-mpm/particle:mp-position)
;;                    (size cl-mpm/particle::mp-domain-size))
;;       mp
;;     (let ((r (max (magicl:tref size 0 0) (magicl:tref size 1 0))))
;;       (- (dot (magicl:.- pos x)) r))))

;; (defun draw-state (file)
;;   (declare (optimize (safety 3) (debug 3)))
;;   (let* ((resolution 10)
;;          (size (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
;;          (width (round (first size) resolution))
;;          (height (round (second size) resolution))
;;          (png (make-instance 'zpng:png
;;                               :color-type :truecolor
;;                               :width width
;;                               :height height))
;;          (image (zpng:data-array png))
;;          (max 255))
;;     (lparallel:pdotimes (y height)
;;       (dotimes (x width)
;;         ;; (print (type-of image))
;;         (setf (aref image y x 0) 255
;;               )
;;         (let ((d 1e10)
;;               (ipos (magicl:scale (magicl:from-list (list x y) '(2 1) :type 'double-float) resolution)))
;;           (loop for mp across (cl-mpm:sim-mps *sim*)
;;                 do
;;                    (setf d (min d (mp-sdf mp ipos)))
;;                 )
;;           (when (< d 0d0)
;;             (setf (aref image y x 0) 255
;;                   (aref image y x 1) 0
;;                   (aref image y x 2) 0)
;;             )))
;;       (zpng:write-png png file)
;;       )))
