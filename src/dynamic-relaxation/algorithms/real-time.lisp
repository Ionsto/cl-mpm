(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)

(defun run-time (sim
                 &key (output-dir "./output/")
                   (post-conv-step (lambda (sim)))
                   (post-iter-step (lambda (sim)))
                   (plotter (lambda (sim)))
                   (dt-scale 0.5d0)
                   (dt 1d0)
                   (total-time 1d0)
                   (substeps 50)
                   (enable-damage t)
                   (enable-plastic t)
                   (mass-scale 1d0)
                   ;; (max-adaptive-steps 5)
                   ;; (min-adaptive-steps -1)
                   (save-vtk-dr t)
                   (save-vtk-loadstep t)
                   (damping 1d0)
                   (initial-quasi-static t)
                   (initial-lstps 1)
                   (elastic-solver 'mpm-sim-dr-ul)
                   (adaptive-mass-enabled nil)
                   (adaptive-mass-scale 1d0)
                   (adaptive-crit 1d-3)
                   (adaptive-hist 1d0)
                   (conv-criteria 1d-3))
  (let ((result t))
    (uiop:ensure-all-directories-exist (list output-dir))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
    (cl-mpm/output::save-simulation-parameters (merge-pathnames output-dir "settings.json")
                                               sim
                                               (list :dt dt
                                                     :criteria-energy conv-criteria
                                                     :criteria-oobf conv-criteria
                                                     :criteria-hist 1d0
                                                     ))
    (save-timestep-preamble sim output-dir)
    (save-conv-preamble sim output-dir)
    (defparameter *total-iter* 0)
    (when initial-quasi-static
      (cl-mpm:iterate-over-mps
       (cl-mpm:sim-mps sim)
       (lambda (mp)
         (when (typep mp 'cl-mpm/particle::particle-damage)
           (setf (cl-mpm/particle::mp-enable-damage mp) nil))
         (when (typep mp 'cl-mpm/particle::particle-plastic)
           (setf (cl-mpm/particle::mp-enable-plasticity mp) nil))))

      (setf (cl-mpm::sim-dt-scale sim) dt-scale)
      (setf (cl-mpm:sim-mass-scale sim) 1d0)

      (setf (cl-mpm:sim-enable-damage sim) nil)
      (let ((vel-algo (cl-mpm::sim-velocity-algorithm sim))
            (temp-sim (copy-instance sim)))
        (change-class temp-sim elastic-solver)
        (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep temp-sim) 0d0)
        (setf (cl-mpm::sim-velocity-algorithm temp-sim) :QUASI-STATIC)
        ;; find initial quasi-static formation
        (dotimes (i initial-lstps)
          (cl-mpm/dynamic-relaxation:converge-quasi-static
           temp-sim
           :energy-crit conv-criteria
           :oobf-crit conv-criteria
           :dt-scale 0.9d0
           :conv-steps 10000
           :substeps substeps
           :damping-factor (sqrt 2d0)
           :post-iter-step
           (lambda (i e o)
             (save-conv-step sim output-dir *total-iter* 0 0d0 o e)
             (incf *total-iter* substeps)
             (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_conv_~5,'0d.vtk" i)) temp-sim)
             (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_conv_nodes__~5,'0d.vtk" i)) temp-sim)
             (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_conv_cells__~5,'0d.vtk" i)) temp-sim)
             (funcall plotter temp-sim)))
          (cl-mpm::finalise-loadstep temp-sim))
        (cl-mpm::reset-grid (cl-mpm:sim-mesh temp-sim))
        (setf (cl-mpm::sim-time temp-sim) 0d0)
        (cl-mpm/dynamic-relaxation::reset-mp-velocity temp-sim)
        (setf (cl-mpm::sim-velocity-algorithm temp-sim) vel-algo)
        ;; (change-class sim sim-type)
        ))

    (setf (cl-mpm:sim-mass-scale sim) mass-scale)
    (setf (cl-mpm:sim-damping-factor sim) (*
                                           (sqrt mass-scale)
                                           damping (cl-mpm/setup:estimate-critical-damping sim)))

    (setf (cl-mpm:sim-enable-damage sim) enable-damage)
    (cl-mpm:iterate-over-mps
     (cl-mpm:sim-mps sim)
     (lambda (mp)
       (when (typep mp 'cl-mpm/particle::particle-damage)
         (setf (cl-mpm/particle::mp-enable-damage mp) enable-damage))
       (when (typep mp 'cl-mpm/particle::particle-plastic)
         (setf (cl-mpm/particle::mp-enable-plasticity mp) enable-plastic))))

    (funcall post-conv-step sim)
    (setf (cl-mpm::sim-dt-scale sim) dt-scale)
    (setf (cl-mpm:sim-dt sim) (* (cl-mpm/setup:estimate-elastic-dt sim :dt-scale dt-scale)))
    (let* ((sim-time 0d0))
      (let ((total-iter 0)
            (dt-accumulator 0d0))
        (let ((state :accelerate)
              (oobf 0d0)
              (energy 0d0)
              (work 0d0)
              (hist adaptive-hist)
              )
          (when adaptive-mass-enabled
            (setf (cl-mpm:sim-mass-scale sim) adaptive-mass-scale))
          (time (loop for steps from 0 to (round total-time dt)
                      while (and (cl-mpm::sim-run-sim sim)
                                 (< sim-time total-time))
                      do
                         (let ((substeps 0))
                           (cl-mpm:sim-format sim t "Step ~d ~%" steps)
                           (cl-mpm/dynamic-relaxation::save-vtks sim output-dir steps)
                           (save-timestep sim output-dir total-iter :DYNAMIC)
                           (incf dt-accumulator dt)
                           (cl-mpm:sim-format sim t "Estimated substeps ~d ~%" (round dt (cl-mpm::sim-dt sim)))
                           (time
                            (loop while (> dt-accumulator 0d0)
                                  do
                                     (progn
                                       (cl-mpm::update-sim sim)

                                       (when adaptive-mass-enabled 
                                         (incf oobf (estimate-static-oobf sim))
                                         (incf energy (cl-mpm::sim-stats-energy sim))
                                         (incf work (estimate-strain-energy sim))
                                         (incf substeps)
                                         )

                                       (incf total-iter)
                                       (decf dt-accumulator (cl-mpm::sim-dt sim))
                                       (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm::calculate-min-dt sim))))))
                           (funcall post-iter-step sim)
                           (funcall plotter sim)
                           (when save-vtk-loadstep
                             (save-vtks sim output-dir steps))
                           (when adaptive-mass-enabled 
                             (setf oobf (/ oobf substeps)
                                   energy (/ energy work))
                             (format t "~A - OOBF ~E - energy ~E~%" state oobf energy)
                             (case state
                               (:accelerate
                                (when (or (> energy (* hist adaptive-crit))
                                          (> oobf (* hist adaptive-crit)))
                                  (format t "Switched to dynamic timestep~%")
                                  (setf
                                   state :dynamic
                                   (cl-mpm:sim-mass-scale sim) 1d0
                                   (cl-mpm:sim-dt sim) (/ (cl-mpm:sim-dt sim) (sqrt adaptive-mass-scale)))))
                               (:dynamic
                                (when (and (< energy (/ adaptive-crit hist))
                                           (< oobf (/ adaptive-crit hist)))
                                  (format t "Switched to accelerate timestep~%")
                                  (cl-mpm/dynamic-relaxation::reset-mp-velocity sim)
                                  (setf
                                   state :accelerate
                                   (cl-mpm:sim-mass-scale sim) adaptive-mass-scale
                                   (cl-mpm:sim-dt sim) (* (cl-mpm:sim-dt sim) (sqrt adaptive-mass-scale)))))))
                           (swank.live:update-swank)))))))
    result))
