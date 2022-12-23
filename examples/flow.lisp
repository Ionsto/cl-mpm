(defpackage :cl-mpm/examples/flow
  (:use :cl))
(sb-ext:restrict-compiler-policy 'speed 3 3)
(sb-ext:restrict-compiler-policy 'debug 0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(in-package :cl-mpm/examples/flow)
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
    h-initial)
    ;(* (magicl:tref (cl-mpm::mp-deformation-gradient mp) dim dim) h-initial))
  )
(defun max-stress (mp)
  (multiple-value-bind (l v) (magicl:eig (cl-mpm::voight-to-matrix (cl-mpm/particle:mp-stress mp)))
    (apply #'max l)
    (magicl:tref (cl-mpm/particle:mp-stress mp) 2 0)))
(defun plot (sim &optional (plot :deformed))
  (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  (when (> (length (cl-mpm:sim-mps sim)) 0)
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
      (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x))))
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

(defun create-new-mp (pos h-x h-y mass)
  (cl-mpm::make-particle 2 
                         'cl-mpm/particle::particle-viscoplastic-damage
                         :E 1e7
                         :nu 0.45
                         ;; :visc-factor 1e-20
                         :visc-factor 1e7
                         :visc-power 3
                         ;; 'cl-mpm/particle::particle-elastic-damage
                         ;; :E 1e9
                         ;; :nu 0.325
                         ;; :visc-factor 1e-24
                         ;; :visc-power 3
                         :mass mass
                         :critical-stress 1d7
                         :gravity -9.8d0
                         :position pos
                         :volume (* h-x h-y)
                         :size (magicl:from-list (list h-x h-y) '(2 1) :type 'double-float)))
(defun create-cell-mps (sim pos h-x h-y mass height mps-per-cell)
  (let* (
         (dy (/ h-y mps-per-cell))
         (dy-inc (/ h-y (+ 0 mps-per-cell)))
         (dy-off (/ h-y (+ 1 mps-per-cell))))
    (loop for i from 0 to (- (* height mps-per-cell) 1)
          do
             (vector-push-extend
              (create-new-mp (mapcar #'+ pos (list 0 (+ dy-off (* i dy-inc)))) dy dy mass)
              (cl-mpm:sim-mps sim)))
        ;; (setf (cl-mpm:sim-mps sim)
        ;;       (concatenate '(vector cl-mpm/particle:particle)
        ;;                    (cl-mpm:sim-mps sim)
        ;;                    (loop for i from 0 to (- (* height mps-per-cell) 1)
        ;;                          collect
        ;;                          (create-new-mp (mapcar #'+ pos (list 0 (+ dy-off (* i dy-inc)))) dy dy mass)))
        ;;       )
    ))
(defun setup-test-column (size block-size block-offset &optional (e-scale 1d0) (mp-scale 1d0))
  (let* ((sim (cl-mpm/setup::make-block (/ 1 e-scale)
                                        (mapcar (lambda (s) (* s e-scale)) size)
                                        #'cl-mpm/shape-function:make-shape-function-linear
                                        ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         ;(e-scale 1)
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (mass (/ (* 900 h-x h-y) (expt mp-scale 2)))
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size)))
    (progn
      (setf (cl-mpm:sim-mps sim) (make-array 1000 :fill-pointer 0))
      (when t
        (let ((block-position
                (mapcar #'+ (list (* h-x (- (+ (/ 1 (* 2 mp-scale))) 0))
                                  (* h-y (+ (/ 1d0 (* 2d0 mp-scale)))))
                        block-offset)))
          (setf (cl-mpm:sim-mps sim)
                (cl-mpm/setup::make-block-mps
                 block-position
                 block-size
                 (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                 'cl-mpm::make-particle
                 ;; 'cl-mpm/particle::particle-viscoelastic-fracture

                 ;; 'cl-mpm/particle::particle-elastic
                 ;; :E 9.5d9
                 ;; :nu 0.499

                 ;; :E 1e6 :nu 1d8

                                        ;'cl-mpm/particle::particle-viscoplastic-damage
                 'cl-mpm/particle::particle-viscoplastic-damage
                 :E 1e7
                 :nu 0.45
                 ;; :visc-factor 1e-20
                 :visc-factor 1e7
                 :visc-power 3

                 ;; Fluid
                 ;; 'cl-mpm/particle::particle-fluid
                 ;; :stiffness 1e5
                 ;; :adiabatic-index 6
                 ;; :rest-density 900d0
                 ;; :viscosity 1.02e-3

                 ;; :E 1e6 :nu 0.33
                 :mass mass
                 :critical-stress 1d7
                 ;; :fracture-toughness 5d0
                 :gravity -9.8d0
                 )))
        (let ((prev-mps (cl-mpm:sim-mps sim)))
          (setf (cl-mpm:sim-mps sim) (make-array 1000 :fill-pointer 0 :adjustable t))
          (loop for mp across prev-mps
                do (vector-push-extend mp (cl-mpm:sim-mps sim)))))
      (setf (cl-mpm:sim-damping-factor sim) 1d-4)
      (setf (cl-mpm:sim-mass-filter sim) 1d-8)
      (setf (cl-mpm:sim-dt sim) 1d-2)

      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-var (cl-mpm:sim-mesh sim)
                                            (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil nil)))
                                            (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil)))
                                            (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0)))
                                            (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0)))
                                            ;; (lambda (i) (cl-mpm/bc:make-bc-friction i
                                            ;;                                         (magicl:from-list '(0d0 1d0)
                                            ;;                                                           '(2 1))
                                            ;;                                         0.25d0))
                                            ))
      (let ((step-x 200)
            (step-y 200)
            (flow-rate 5))
        ;; (loop for x from 0 to (- (round step-x h) 1)
        ;;       do (loop for y from 0 to (- (round step-y h) 1)
        ;;                do (push (list x y) (cl-mpm/mesh::mesh-boundary-shapes (cl-mpm:sim-mesh sim)))))

        (loop for x from 0 to (round step-x h)
              do (loop for y from 0 to (round step-y h)
                       do (push
                           ;(cl-mpm/bc:make-bc-fixed (list x y) '(nil 0))
                           (cl-mpm/bc:make-bc-friction (list x y)
                                                                   (magicl:from-list '(0d0 1d0)
                                                                                     '(2 1))
                                                                   0.25d0)
                                (cl-mpm:sim-bcs sim))))
        (when t
          (loop for x from 1 to 1
                do (loop for y from (round 200 h) to (round 200 h)
                         do (push (cl-mpm/bc::make-bc-inflow (list x y)
                                                             (list flow-rate 0)
                                        ;(* h-x h-y (/ 0.125 2))
                                                             (* h-x h-y (* 1/16 (expt mp-scale 1)))
                                                             (lambda (pos)
                                                               (create-cell-mps sim pos h-x h-y mass
                                                                                (round (- 200 100) h)
                                                                                mp-scale)))
                                  (cl-mpm:sim-bcs sim))))
          (loop for x from 0 to 1
                do (loop for y from (round 200 h) to (round 300 h)
                         do (push (cl-mpm/bc::make-bc-fixed (list x y)
                                                            (list flow-rate nil))
                                  (cl-mpm:sim-bcs sim)))))
        ;; (loop for y from 0 to (round step-y h)
        ;;       do (push (cl-mpm/bc:make-bc-surface (list (round step-x h) y)
        ;;                                           (magicl:from-list '(1d0 0d0) '(2 1)))
        ;;                (cl-mpm:sim-bcs sim)))
        )
      (let ((ocean-x 200)
            (ocean-y 200))
        (setf (cl-mpm::sim-bcs-force sim)
              (loop for x from (floor ocean-x h) to (floor (first size) h)
                    append (loop for y from 0 to (floor ocean-y h)
                                 collect (cl-mpm/bc::make-bc-buoyancy
                                          (list x y)
                                          (magicl:from-list (list 1d1 (* 9.8d0 (* 1d0 1000))) '(2 1)))))))
      sim)))

;Setup
(defun setup ()
  (defparameter *sim* (setup-test-column '(2000 400) '(200 100) '(000 200) (/ 1 25) 2))
  ;; (defparameter *sim* (setup-test-column '(1 1) '(1 1) '(0 0) 1 1))
  ;;(remove-sdf *sim* (ellipse-sdf (list 400 100) 10 40))
  ;; (remove-sdf *sim* (ellipse-sdf (list 1.5 3) 0.25 0.5))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  (defparameter *sim-step* 0)
  (defparameter *run-sim* nil)
  (defparameter *run-sim* t)
  ;; (defparameter *load-mps*
  ;;   (let* ((mps (cl-mpm:sim-mps *sim*))
  ;;          (least-pos
  ;;             (apply #'max (loop for mp across mps
  ;;                                collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))))
  ;;     (loop for id from 0 to (- (length mps) 1) when (>= (magicl:tref (cl-mpm/particle:mp-position (aref mps id)) 0 0) (- least-pos 0.001))
  ;;             collect (aref mps id))))
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
    ;; (vgplot:axis (list 0 (nth 0 (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*))) 
    ;;                    0 (nth 1 (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))))
    (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
      (vgplot:format-plot t "set ytics ~f" h)
      (vgplot:format-plot t "set xtics ~f" h))
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
                    (time (dotimes (i 100)
                           ;; (pescribe-velocity *sim* *load-mps* (magicl:from-list '(0.5d0 0d0) '(2 1)))
                           (cl-mpm::update-sim *sim*)
                            (with-accessors ((mps cl-mpm:sim-mps))
                                  *sim*
                                (delete-if (lambda (mp)
                                             (with-accessors ((pos cl-mpm/particle:mp-position))
                                                 mp
                                               (>= (magicl:tref pos 0 0) 1800d0))) mps))
                           ;; (setf max-cfl (max max-cfl (find-max-cfl *sim*)))
                           ;; (cl-mpm/eigenerosion:update-fracture *sim*)
                           (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))
                           ))
                    ;; (format t "Max cfl: ~f~%" max-cfl)
                    )
                  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*))
                                          *sim*)
                  (incf *sim-step*)
                  (plot *sim*)
                  ;; (vgplot:print-plot (asdf:system-relative-pathname "cl-mpm" (format nil "output/frame_~5,'0d.png" steps)))
                  (swank.live:update-swank)
                  (sleep .01)

                  )))
    ;; (vgplot:figure)
    ;; (vgplot:title "Velocity over time")
    ;; (vgplot:plot *time* *velocity*)
    )

(setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))
;; (require :sb-sprof)
;; (setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))
;; (setf sb-sprof:*max-samples* 10000)
;; (defparameter *sim* (setup-test-column '(100 100) '(100 100) '(0 0) (/ 1 20) 2))
;; (progn
;;   (format t "MP count:~D~%" (length (cl-mpm:sim-mps *sim*)))
;;   ;; (sb-sprof:reset)
;;   ;; (sb-sprof:start-profiling :mode :time)
;;   ;; (time (lparallel:pdotimes (i 10000000)
;;   ;;                           (magicl:from-list '(1d0 1d0) '(2 1))))
;;   (time (dotimes (i 1000)
;;           (cl-mpm::update-sim *sim*)))
;;   ;; (sb-sprof:stop-profiling)
;;   )
;; (sb-sprof:report :type :flat)
;; (time (dotimes (i 1)
;;         (cl-mpm::calculate-strain-rate-test *test-mesh* *test-mp* 1)))

;; (defun profile-funcs ()
;;   (sb-profile:unprofile)
;;   (sb-profile:reset)
;;   (format t "MP count:~D~%" (length (cl-mpm:sim-mps *sim*)))
;;   (sb-profile:profile cl-mpm::update-sim
;;                       cl-mpm::p2g
;;                       cl-mpm::g2p
;;                       cl-mpm::calculate-strain-rate
;;                       cl-mpm::update-stress
;;                       ;; cl-mpm::reset-grid
;;                       ;; cl-mpm::p2g
;;                       ;; cl-mpm::filter-grid
;;                       ;; cl-mpm::update-nodes
;;                       ;; cl-mpm::apply-bcs
;;                       ;; cl-mpm::g2p
;;                       ;; cl-mpm::update-particle
;;                       ;; cl-mpm::update-stress-mp
;;                       ;; cl-mpm::update-stress
;;                       ;; cl-mpm::calculate-strain-rate
;;                       ;; cl-mpm::iterate-over-neighbours-shape
;;                       ;; cl-mpm::iterate-over-neighbours-shape-linear
;;                       ;; cl-mpm::p2g-mp
;;                       ;; cl-mpm::g2p-mp
;;                       ;; cl-mpm::p2g-mp-node
;;                       ;; cl-mpm::g2p-mp-node
;;                       )
;;   (time (loop repeat 100
;;              do (progn
;;                   (cl-mpm::update-sim *sim*))))
;;   (sb-profile:report))
