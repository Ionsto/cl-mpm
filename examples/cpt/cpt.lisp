(defpackage :cl-mpm/examples/cpt
  (:use :cl
   :cl-mpm/example
   :cl-mpm/utils)
  )
(in-package :cl-mpm/examples/cpt)
(declaim (notinline plot-domain))
(defun plot-domain ()
  (when *sim*
    ;; (cl-mpm::g2p (cl-mpm:sim-mesh *sim*)
    ;;              (cl-mpm:sim-mps *sim*)
    ;;              (cl-mpm:sim-dt *sim*)
    ;;              0d0
    ;;              :QUASI-STATIC)
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :trial t
     :colour-func (lambda (mp) (cl-mpm/particle::mp-index mp)))))

(defun setup (&key (refine 1) (mps 3))
  (let* ((L 10d0)
         (d 1d0)
         (domain-width 10d0)
         (h (/ 1d0 refine))
         (height 10d0)
         (width 5d0)
         (domain-height (+ height 4d0))
         (density 16d3)
         (E 22d9)
         (domain-size (list domain-width domain-height))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list height height)))

    (setf *sim* (cl-mpm/setup::make-simple-sim
                 h
                 element-count
                 :sim-type 'cl-mpm/aggregate::mpm-sim-agg-usf
                 ;; :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
                 :args-list (list :enable-aggregate t
                                  :enable-fbar nil
                                  :enable-split nil
                                  :vel-algo :FLIP
                                  )))
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0 0)
      block-size
      (mapcar (lambda (e) (* (/ e h) mps)) block-size)
      density
      ;; 'cl-mpm/particle::particle-elastic
      ;; :E E
      ;; :nu 0.3d0
      ;; :rho 20d3
      ;; 'cl-mpm/particle::particle-vm
      ;; :E E
      ;; :nu 0.3d0
      ;; :rho 20d3
      'cl-mpm/particle::particle-mc
      :E E
      :nu 0.3d0
      :phi (cl-mpm/utils:deg-to-rad 32.8d0)
      :phi (cl-mpm/utils:deg-to-rad 2.8d0)
      :c 0.3d3
      ;; 'cl-mpm/particle::particle-elastic
      ;; :E E
      ;; :nu 0.2d0
      ;; :index 0
      ;; :gravity-axis (cl-mpm/utils:vector-zeros)
      ))
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-5)
    ;; (setf (cl-mpm::sim-ghost-factor *sim*) (* E 1d-4))
    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 nil nil)
     :bottom '(0 0 nil)
     :right '(0 0 nil)
     )

    (defparameter *penalty*
      (make-cpt
       *sim*
       E
       :epsilon-scale 1d1
       )
      ;; (cl-mpm/penalty::make-bc-penalty-distance-point
      ;;  *sim*
      ;;  (cl-mpm/utils:vector-from-list '(0d0 -1d0 0d0))
      ;;  (cl-mpm/utils:vector-from-list (list (* 0.5d0 domain-width)
      ;;                                       height
      ;;                                       0d0))
      ;;  (/ domain-width 2)
      ;;  (* E 1d1)
      ;;  0.5d0
      ;;  0d0)
      )
    (defparameter *current-inc* 0d0)
    (cl-mpm:add-bcs-force-list
     *sim*
     *penalty*)
    ;; (setf (cl-mpm:sim-dt *sim*)
    ;;       (* 0.5d0 (cl-mpm/setup:estimate-elastic-dt *sim*)))
    (setf *run-sim* t))
  (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
  )
(defun make-cpt (sim E &key (epsilon-scale 1d2))
  (let ((epsilon (* E epsilon-scale))
        (width 1d0)
        (height 10d0)
        (cone-height 1d0)
        (offset -1d-5)
        (offset-y 10d0)
        (damping 0d0)
        (friction 0d0)
        )
    (let (
          (A
            (cl-mpm/penalty::make-bc-penalty-line-segment
             sim
             (cl-mpm/utils:vector-from-list (list (+ offset width) (+ offset-y cone-height) 0d0))
             (cl-mpm/utils:vector-from-list (list offset offset-y 0d0))
             epsilon
             friction
             damping
             ))
          (B
            (cl-mpm/penalty::make-bc-penalty-line-segment
             sim
             (cl-mpm/utils:vector-from-list (list (+ offset width) (+ offset-y height) 0d0))
             (cl-mpm/utils:vector-from-list (list (+ offset width) (+ offset-y cone-height) 0d0))
             epsilon
             friction
             damping
             ))
          (C
            (cl-mpm/penalty::make-bc-penalty-line-segment
             sim
             (cl-mpm/utils:vector-from-list (list offset offset-y 0d0))
             (cl-mpm/utils:vector-from-list (list (- offset width) (+ offset-y cone-height) 0d0))
             epsilon
             friction
             damping
             ))
          (D
            (cl-mpm/penalty::make-bc-penalty-line-segment
             sim
             (cl-mpm/utils:vector-from-list (list (- offset width) (+ offset-y cone-height) 0d0))
             (cl-mpm/utils:vector-from-list (list (- offset width) (+ offset-y height) 0d0))
             epsilon
             friction
             damping
             ))
          (FLAT
            (cl-mpm/penalty::make-bc-penalty-line-segment
             sim
             (cl-mpm/utils:vector-from-list (list (+ offset width) (+ offset-y cone-height) 0d0))
             (cl-mpm/utils:vector-from-list (list (- offset width) (+ offset-y cone-height) 0d0))
             epsilon
             friction
             damping
             ))
          )
      (cl-mpm/penalty::make-bc-penalty-structure
       sim
       epsilon
       friction
       damping
       (list
        A
        B
        ;; C
        ;; D
        ))
      )))

(defun run ()
  (defparameter *data-disp* (list))
  (defparameter *data-load* (list))
  (let* ((lstps 50)
         (total-disp -5d0)
         (delta (/ total-disp lstps))
         (current-disp 0d0)
         )
    (cl-mpm/dynamic-relaxation::run-load-control
     *sim*
     :output-dir (format nil "./output/")
     :plotter (lambda (sim) (plot-domain))
     :loading-function
     (lambda (i)
       (setf current-disp (* i total-disp))
       (cl-mpm/penalty::bc-increment-center
        *penalty*
        (cl-mpm/utils:vector-from-list (list 0d0 delta 0d0))))
     :load-steps lstps
     :kinetic-damping nil
     :damping 1d0
     :substeps 50
     :criteria 1d-3
     :save-vtk-dr t
     :save-vtk-loadstep t
     :dt-scale 0.05d0
     :post-conv-step
     (lambda (sim)
       (push current-disp *data-disp*)
       (push (cl-mpm/penalty::resolve-load-direction *penalty* (cl-mpm/utils:vector-from-list (list 0d0 -1d0 0d0)))
             *data-load*)
       )
     ))
  (vgplot:figure)
  (vgplot:plot (reverse *data-disp*) (reverse *data-load*)))

(defun run ()
  (defparameter *data-disp* (list))
  (defparameter *data-load* (list))
  (let* ((lstps 50)
         (total-disp -5d0)
         (delta (/ total-disp lstps))
         (current-disp 0d0)
         )
    )
  (vgplot:figure)
  (vgplot:plot (reverse *data-disp*) (reverse *data-load*)))

(defun test ()
  (setup :mps 4 :refine 1)
  ;; (run)
  (run-real-time)
  )


(defun run-real-time (&key (output-dir "./output/")
              (ms 1d0)
              (dt-scale 0.5d0)
              )
  (uiop:ensure-all-directories-exist (list output-dir))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (cl-mpm/output:save-vtk-mesh (merge-pathnames output-dir "mesh.vtk")
                          *sim*)
  (cl-mpm/dynamic-relaxation::elastic-static-solution
   *sim*
   :crit 1d-3
   )
  (cl-mpm/dynamic-relaxation::set-mp-plastic-damage *sim* :enable-plastic t)
  (let* ((target-time 0.05d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-min (cl-mpm:sim-dt *sim*)))
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 0.1d0
             (cl-mpm/setup:estimate-critical-damping *sim*)))
    (setf (cl-mpm:sim-mass-scale *sim*) ms)
    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup:estimate-elastic-dt *sim*))
    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (setf substeps (round target-time (cl-mpm:sim-dt *sim*)))
    (cl-mpm/dynamic-relaxation::save-timestep-preamble output-dir)
    (cl-mpm/output::save-simulation-parameters (merge-pathnames "settings.json" output-dir)
                                               *sim*
                                               (list :dt target-time))
    (format t "Substeps ~D~%" substeps)
    (setf (cl-mpm::sim-stats-work *sim*) 0d0)
    (let* ((total-disp -5d0)
           (total-time 10d0)
           (delta (/ total-disp total-time)))
      (time (loop for steps from 0 to (round total-time target-time)
                  while *run-sim*
                  do
                     (let ((energy 0d0)
                           (work 0d0)
                           (oobf 0d0))
                       (cl-mpm/dynamic-relaxation::save-timestep *sim* output-dir steps :DYNAMIC)
                       (format t "Step ~d ~%" steps)
                       (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" steps)) *sim*)
                       (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" steps)) *sim*)
                       (cl-mpm/output::save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d.vtk" steps)) *sim*)
                       (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim_p_~5,'0d.vtk" steps)) *sim* )
                       ;; (cl-mpm/output::save-vtk-mesh-nodes (merge-pathnames output-dir (format nil "sim_nodes_p_~5,'0d.vtk" *sim-step*)) (cl-mpm::sim-mesh-p *sim*))
                       (setf dt-min (cl-mpm::calculate-min-dt *sim*))
                       (time
                        (dotimes (i substeps)
                          (cl-mpm/penalty::bc-increment-center
                           *penalty*
                           (cl-mpm/utils:vector-from-list (list 0d0 (* (cl-mpm:sim-dt *sim*) delta) 0d0)))
                          (cl-mpm::update-sim *sim*)
                          ;; (setf dt-min (min (cl-mpm::calculate-min-dt *sim*) dt-min))
                          ;; (setf (cl-mpm:sim-dt *sim*) (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
                          (incf oobf (cl-mpm::sim-stats-oobf *sim*))
                          (incf energy (cl-mpm::sim-stats-energy *sim*))
                          (incf work (cl-mpm::sim-stats-power *sim*))
                          ))
                       (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                         (format t "CFL dt estimate: ~f~%" dt-e)
                         (format t "CFL step count estimate: ~D~%" substeps-e)
                         (setf substeps substeps-e))
                       (format t "Substeps ~D~%" substeps)
                       (format t "~E ~E~%" energy oobf)

                       (plot-domain)
                       (swank.live:update-swank)
                       ))))))
