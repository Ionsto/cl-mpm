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
     :colour-func (lambda (mp) (cl-mpm/particle::mp-index mp))
     )
    (let* ((disp *disp*)
           (x0 0d0)
           (y0 (+ disp *cpt-offset-y*))
           (x1 *cpt-width*)
           (y1 (+ disp *cpt-offset-y* *cpt-cone-height*))
           (x2 *cpt-width*)
           (y2 (+ disp *cpt-offset-y* *cpt-height*))
           (x3 0d0)
           (y3 (+ disp *cpt-offset-y* *cpt-height*))
           )
      (vgplot:format-plot t "set object 1 polygon from ~f,~f to ~f,~f to ~f,~f to ~f,~f fc rgb 'black' fs transparent solid 1 noborder behind"
                          x0 y0
                          x1 y1
                          x2 y2
                          x3 y3))))

(defun setup (&key (refine 1) (mps 3))
  (let* ((L 10d0)
         (d 1d0)
         (domain-width 10d0)
         (h (/ 1d0 refine))
         (height 10d0)
         (domain-height (+ height 4d0))
         (density 16d3)
         (E 22d9)
         (domain-size (list domain-width domain-height))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list domain-width height)))

    (setf *sim* (cl-mpm/setup::make-simple-sim
                 h
                 element-count
                 ;; :sim-type 'cl-mpm/aggregate::mpm-sim-agg-usf
                 ;; :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic
                 :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
                 :args-list (list :enable-aggregate t
                                  :enable-fbar t
                                  :enable-split t
                                  :split-factor (* 1.2d0 (sqrt 2) (/ 1d0 mps))
                                  :vel-algo :QUASI-STATIC
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
      ;; 'cl-mpm/particle::particle-vm
      ;; :E E
      ;; :nu 0.3d0
      ;; :rho 0.2d3
      'cl-mpm/particle::particle-mc
      :E E
      :nu 0.3d0
      :phi (cl-mpm/utils:deg-to-rad 32.8d0)
      :psi 0d0;(cl-mpm/utils:deg-to-rad 2.8d0)
      :c 0.3d3
      ;; 'cl-mpm/particle::particle-elastic
      ;; :E E
      ;; :nu 0.2d0
      ;; :index 0
      ;; :gravity-axis (cl-mpm/utils:vector-zeros)
      ))
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-9)
    ;; (setf (cl-mpm::sim-ghost-factor *sim*) (* E 1d-3))
    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 nil nil)
     :bottom '(0 0 nil)
     :right '(0 nil nil)
     )

    (defparameter *penalty*
      (make-cpt
       *sim*
       E
       :epsilon-scale 1d1
       :offset-y height
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
(defun make-cpt (sim E &key (epsilon-scale 1d2)
                         (offset-y 20d0)
                         )
  (let ((epsilon (* E epsilon-scale))
        (width 1d0)
        (height 10d0)
        (cone-height 1d0)
        (offset -1d-5)
        (damping 0d0)
        (friction 0d0))
    (defparameter *cpt-width* width)
    (defparameter *cpt-height* height)
    (defparameter *cpt-cone-height* cone-height)
    (defparameter *cpt-offset-y* offset-y)
    (defparameter *cpt-offset* offset)
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
         (current-disp 0d0))
    (defparameter *disp* 0d0)
    (let ((step 0))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
      (cl-mpm/dynamic-relaxation::run-adaptive-load-control
       *sim*
       :output-dir (format nil "./output/")
       :plotter (lambda (sim) (plot-domain))
       :loading-function
       (lambda (i)
         (setf (cl-mpm:sim-gravity *sim*) -9.8d0)
         (setf current-disp (* i total-disp))
         (setf *disp* current-disp)
         (cl-mpm/penalty::bc-set-displacement
          *penalty*
          (cl-mpm/utils:vector-from-list (list 0d0 current-disp 0d0))))
       :post-conv-step (lambda (sim)
                         (plot-domain)
                         (vgplot:title (format nil "Step ~D" step))
                         (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step)) :terminal "png size 1920,1080")
                         (incf step))
       :load-steps lstps
       :kinetic-damping nil
       :damping 1d0
       :substeps 50
       :criteria 1d-3
       :save-vtk-dr t
       :save-vtk-loadstep t
       :dt-scale 0.5d0
       :post-conv-step
       (lambda (sim)
         (push current-disp *data-disp*)
         (push (cl-mpm/penalty::resolve-load-direction *penalty* (cl-mpm/utils:vector-from-list (list 0d0 -1d0 0d0)))
               *data-load*)))))
  (vgplot:figure)
  (vgplot:plot (reverse *data-disp*) (reverse *data-load*)))

;; (defun run ()
;;   (defparameter *data-disp* (list))
;;   (defparameter *data-load* (list))
;;   (let* ((lstps 50)
;;          (total-disp -5d0)
;;          (delta (/ total-disp lstps))
;;          (current-disp 0d0)
;;          )
;;     )
;;   (vgplot:figure)
;;   (vgplot:plot (reverse *data-disp*) (reverse *data-load*)))

(defun test ()
  (setup :mps 3 :refine 2)
  (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static)
  (run)
  ;; (run-real-time)
  )

(defun test-implicit ()
  (setup :mps 2 :refine 4)
  (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic)
  (run-real-time :dt-scale 100d0 :target-time 0.1d0)
  )


(defun run-real-time (&key (output-dir "./output/")
                        (ms 1d0)
                        (dt-scale 0.1d0)
                        (target-time 0.1d0)
                        )
  (uiop:ensure-all-directories-exist (list output-dir))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (cl-mpm/output:save-vtk-mesh (merge-pathnames output-dir "mesh.vtk")
                          *sim*)
  (cl-mpm/dynamic-relaxation::save-timestep-preamble output-dir)
  (cl-mpm/dynamic-relaxation::save-conv-preamble output-dir)
  (cl-mpm/dynamic-relaxation::elastic-static-solution
   *sim*
   :crit 1d-3
   )
  (cl-mpm/dynamic-relaxation::set-mp-plastic-damage *sim* :enable-plastic t)
  (let* (
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-min (cl-mpm:sim-dt *sim*)))
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 0.1d0
             (cl-mpm/setup:estimate-critical-damping *sim*)))
    (setf (cl-mpm:sim-mass-scale *sim*) ms)
    ;; (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup:estimate-elastic-dt *sim*))
    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))

    (setf substeps (round target-time (cl-mpm:sim-dt *sim*)))
    ;; (setf (cl-mpm:sim-dt *sim*) target-time)
    ;; (setf substeps 1)

    (cl-mpm/dynamic-relaxation::save-timestep-preamble output-dir)
    (cl-mpm/output::save-simulation-parameters (merge-pathnames "settings.json" output-dir)
                                               *sim*
                                               (list :dt target-time))
    (format t "Substeps ~D~%" substeps)
    (setf (cl-mpm::sim-stats-work *sim*) 0d0)
    (let* ((total-disp -10d0)
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
                       ;; (setf dt-min (cl-mpm::calculate-min-dt *sim*))
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
                       ;; (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       ;;   (format t "CFL dt estimate: ~f~%" dt-e)
                       ;;   (format t "CFL step count estimate: ~D~%" substeps-e)
                       ;;   (setf substeps substeps-e))
                       (format t "Substeps ~D~%" substeps)
                       (format t "~E ~E~%" energy oobf)

                       (plot-domain)
                       (swank.live:update-swank)
                       ))))))
