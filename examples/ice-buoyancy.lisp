(defpackage :cl-mpm/examples/ice-buoyancy
  (:use :cl
   :cl-mpm/example))
(in-package :cl-mpm/examples/ice-buoyancy)

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-mc) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  ;; (cl-mpm::co-domain-corner-2d mesh mp dt)
  (cl-mpm::update-domain-polar-2d mesh mp dt)
  (cl-mpm::scale-domain-size mesh mp))

(defun plot-domain ()
  (when *sim*
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     ;:colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle:mp-stress mp) :xx))
     :colour-func #'cl-mpm/particle::mp-strain-plastic-vm
     )))
(defun setup (&key (refine 1) (mps 2))
  (let* ((density 916.7d0)
         (mesh-resolution (/ 10d0 refine))
         (offset (* mesh-resolution 0d0))
         (datum (+ 10d0 offset))
         (domain-size (list 300d0 300d0))
         (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
         (block-size (list 200d0 200d0)))
    (setf *sim* (cl-mpm/setup::make-simple-sim mesh-resolution element-count :sim-type 'cl-mpm::mpm-sim-usf))

    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0 0)
      block-size
      (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
      density
      'cl-mpm/particle::particle-mc
      :E 1d9
      :nu 0.24d0
      :gravity -9.8d0
      :psi 0d0
      :phi (* 40d0 (/ pi 180))
      :phi-r (* 20d0 (/ pi 180))
      :psi-r 0d0
      :c 100d3
      :c-r 1d3
      :softening 10d0
      :enable-plasticity nil
      ))
    (setf
     (cl-mpm:sim-bcs *sim*)
     (cl-mpm/bc::make-outside-bc-varfix
      (cl-mpm:sim-mesh *sim*)
      '(0 nil 0)
      '(0 nil 0)
      '(nil 0 0)
      '(nil 0 0)
      '(nil nil 0)
      '(nil nil 0)))
    (setf (cl-mpm:sim-mass-scale *sim*) 1d4)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 0.1d0
             (sqrt 1d4)
             (cl-mpm/setup:estimate-critical-damping *sim*)))
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-2)
    (setf (cl-mpm::sim-enable-fbar *sim*) t)
    (setf (cl-mpm::sim-allow-mp-split *sim*) t)



    (setf (cl-mpm:sim-dt *sim*)
          (* 0.5d0 (cl-mpm/setup:estimate-elastic-dt *sim*)))
    ;; (setf (cl-mpm::sim-enable-damage *sim*) t)
    (setf *run-sim* t)
    (cl-mpm:add-bcs-force-list
     *sim*
     (cl-mpm/buoyancy::make-bc-buoyancy-clip
      *sim*
      100d0
      1000d0
      (lambda (pos) t)))
    ))



(defun run (&key (output-dir "./output/"))
  (uiop:ensure-all-directories-exist (list (uiop:merge-pathnames* output-dir)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (cl-mpm/dynamic-relaxation:converge-quasi-static
   *sim*
   :oobf-crit 1d-1
   :energy-crit 1d-1
   )
  (cl-mpm::iterate-over-mps
   (cl-mpm:sim-mps *sim*)
   (lambda (mp)
     (setf (cl-mpm/particle::mp-enable-plasticity mp) t)))

  (setf (cl-mpm:sim-mass-scale *sim*) 1d4)
  (setf (cl-mpm:sim-damping-factor *sim*)
        (* 0.01d0
           ;; (sqrt 1d4)
           (cl-mpm/setup:estimate-critical-damping *sim*)))

  (setf (cl-mpm:sim-dt *sim*)
        (* 0.8d0 (cl-mpm/setup:estimate-elastic-dt *sim*)))
  (let* ((dt-scale 1d0)
         (target-time 1d0)
         (substeps (ceiling target-time (cl-mpm:sim-dt *sim*))))
    (format t "Substeps ~D~%" substeps)
    (loop for step from 0 below 100
          while *run-sim*
          do
             (format t "Step ~D~%" step)
             (time
              (dotimes (i substeps)
                (cl-mpm:update-sim *sim*)
                (cl-mpm::remove-mps-func
                 *sim*
                 (lambda (mp)
                   (> (magicl:det (cl-mpm/particle::mp-deformation-gradient mp)) 1.5d0)))
                ;; (cl-mpm::split-mps-eigenvalue *sim*)
                ))
             (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir (format nil "sim_~5,'0d.vtk" step)) *sim* )
             (cl-mpm/output:save-vtk-nodes (uiop:merge-pathnames* output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) *sim* )
             (plot-domain)
             (swank.live:update-swank))))

(defun test-buoy ()

  (let ((step 0)
        (output-dir "./output/"))
    (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir (format nil "sim_~5,'0d.vtk" step)) *sim* )
    (cl-mpm/output:save-vtk-nodes (uiop:merge-pathnames* output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) *sim* )
    (cl-mpm:update-sim *sim*)
    (incf step)
    (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir (format nil "sim_~5,'0d.vtk" step)) *sim* )
    (cl-mpm/output:save-vtk-nodes (uiop:merge-pathnames* output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) *sim* )))
