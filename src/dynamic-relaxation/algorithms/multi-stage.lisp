(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)

(declaim (notinline step-real-time))
(defun step-real-time (sim
                       global-step
                       &key
                         (target-time 1d0)
                         (step-time 100d0)
                         (plotter (lambda (sim)))
                         (dt-scale 0.5)
                         (output-dir "./output/")
                         (max-steps 1000)
                         (damping 1d-4)
                         (damping-0 nil)
                         (enable-damage t)
                         (criteria 1d-3)
                         (enable-plastic t)
                         (enable-mass-scaling t)
                         (mass-scaler 1d1))

  (let ((state :dynamic))
    (setf (cl-mpm:sim-mass-scale sim) 1d0)
    (when enable-mass-scaling
        (setf
         state :accelerate
         (cl-mpm:sim-mass-scale sim) mass-scaler))
    (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm/setup:estimate-elastic-dt sim)))
    (setf (cl-mpm:sim-damping-factor sim) (* damping
                                             (if damping-0
                                                 damping-0
                                                 (cl-mpm/setup:estimate-critical-damping sim))))
    (set-mp-plastic-damage sim :enable-damage enable-damage :enable-plastic enable-plastic)
    (cl-mpm::update-stiffness-mps sim)
    (cl-mpm::reset-grid (cl-mpm:sim-mesh sim))
    (reset-mp-velocity sim)
    (setf (cl-mpm:sim-enable-damage sim) enable-damage)
    (let* ((e-crit    (* 0.5d0 criteria))
           (oobf-crit (* 0.5d0 criteria))
           (energy e-crit)
           (oobf oobf-crit)
           (work 0d0)
           (intertial-passed t)
           (dt-0 (* dt-scale (cl-mpm/setup:estimate-elastic-dt sim)))
           (substeps (max 1 (round target-time (cl-mpm:sim-dt sim)))))
      (cl-mpm:sim-format sim t "Substeps ~D~%" substeps)
      (cl-mpm:sim-format sim t "E crit ~E - OOBF crit ~E~%" e-crit oobf-crit)
      (time (loop for step from 0 to max-steps
                  while (and (cl-mpm::sim-run-sim sim)
                             (or (>= energy e-crit)
                                 (>= oobf oobf-crit)
                                 (not intertial-passed)
                                 (< step (/ step-time target-time))))
                  do
                     (let ((substeps (max 1 (round target-time (cl-mpm:sim-dt sim)))))
                       (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_real_~5,'0d_~5,'0d.vtk" global-step step)) sim)
                       (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_real_nodes_~5,'0d_~5,'0d.vtk" global-step step)) sim)
                       (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_real_cells_~5,'0d_~5,'0d.vtk" global-step step)) sim)
                       (setf energy 0d0)
                       (setf oobf 0d0)
                       (setf work 0d0)
                       (cl-mpm:sim-format sim t "Real time step ~d - substeps ~d - time ~E - dt ~E~%" step substeps target-time (cl-mpm:sim-dt sim))
                       (time
                        (dotimes (i substeps)
                          (cl-mpm::update-sim sim)
                          ;; (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm::calculate-min-dt sim)))
                          (incf oobf (estimate-static-oobf sim))
                          (incf energy (cl-mpm::sim-stats-energy sim))
                          (incf work (estimate-strain-energy sim))
                          ;; (incf work (cl-mpm::sim-stats-power sim))
                          ))

                       (setf
                        energy (/ energy substeps)
                        oobf (/ oobf substeps)
                        work (/ work substeps)
                        )
                       (if (= work 0d0)
                           (setf energy 0d0)
                           (setf energy (abs (/ energy work))))

                       (let* ((hist 1d0)
                              (hist-power 0.5d0)
                              (hist-energy (* hist (expt e-crit hist-power)))
                              (hist-oobf (* hist (expt e-crit hist-power))))
                         (when (or (> energy hist-energy)
                                   (> oobf hist-oobf))
                           (cl-mpm:sim-format sim t "Inertia passed~%")
                           (setf intertial-passed t)))
                       (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm/setup::estimate-elastic-dt sim)))
                       (when enable-mass-scaling
                         (let* ((hist 2d0)
                                (hist-power 0.6d0)
                                (hist-energy (expt e-crit hist-power))
                                (hist-oobf (expt e-crit hist-power)))
                           (format t "Current state ~A - ~E ~E~%" state hist-energy hist-oobf)
                           (case state
                             (:accelerate
                              (when (or ;; (> energy (* hist-energy hist))
                                     (> oobf (* hist-oobf hist)))
                                (format t "Switched to dynamic timestep~%")
                                (setf
                                 state :dynamic
                                 (cl-mpm:sim-mass-scale sim) 1d0
                                 (cl-mpm:sim-dt sim) (/ (cl-mpm:sim-dt sim) (sqrt mass-scaler) ))))
                             (:dynamic
                              (when (and ;; (< energy (/ hist-energy hist))
                                     (< oobf (/ hist-oobf hist)))
                                (format t "Switched to accelerate timestep~%")
                                (reset-mp-velocity sim)
                                (setf
                                 state :accelerate
                                 (cl-mpm:sim-mass-scale sim) mass-scaler
                                 (cl-mpm:sim-dt sim) (* (cl-mpm:sim-dt sim) (sqrt mass-scaler))))))))
                       (cl-mpm:sim-format sim t "Residuals ~E ~E ~%" energy oobf)
                       (setf (cl-mpm::sim-stats-oobf sim) oobf)
                       (save-timestep sim output-dir global-step state)
                       (funcall plotter sim)
                       (swank.live:update-swank)))))))

(declaim (notinline run-multi-stage))
(defun run-multi-stage (sim
                        &key (output-dir "./output/")
                          (conv-load-steps 1)
                          (post-conv-step (lambda (sim)))
                          (plotter (lambda (sim)))
                          (dt-scale 1d0)
                          (conv-dt-scale nil)
                          (damping-factor (sqrt 2d0))
                          (explicit-dt-scale 0.9d0)
                          (explicit-damping-factor 1d-3)
                          (explicit-mass-scaling t)
                          (elastic-dt-margin 1000)
                          (elastic-dt nil)
                          (dt 1d0)
                          (steps nil)
                          (total-time 1d0)
                          (max-adaptive-steps 5)
                          (min-adaptive-steps -1)
                          (adaption-constant 2)
                          (conv-criteria 1d-3)
                          (explicit-conv-criteria nil)
                          (substeps 50)
                          (sub-conv-steps 50)
                          (enable-damage t)
                          (enable-plastic t)
                          (save-vtk-dr t)
                          (save-vtk-loadstep t)
                          (max-plastic-inc 1d0)
                          (max-deformation-gradient 10d0)
                          (max-damage-inc 0.6d0)
                          (min-damage-inc 0d0)
                          (setup-quasi-static (lambda (sim)))
                          (setup-dynamic (lambda (sim)))
                          (elastic-solver 'mpm-sim-quasi-static)
                          (explicit-dynamic-solver 'cl-mpm::mpm-sim-usf))
  (let ((damping-0 (cl-mpm/setup::estimate-critical-damping sim)))
    (when steps
      (setf total-time (* dt steps)))
    (unless explicit-conv-criteria
      (setf explicit-conv-criteria conv-criteria))
    (uiop:ensure-all-directories-exist (list output-dir))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
    (cl-mpm/output::save-simulation-parameters (merge-pathnames output-dir "settings.json")
                                               sim
                                               (list :dt dt
                                                     :criteria-energy conv-criteria
                                                     :criteria-oobf conv-criteria
                                                     :criteria-hist 1d0))
    (save-timestep-preamble sim output-dir)
    (save-conv-preamble sim output-dir)
    (setf (cl-mpm::sim-dt sim) 0d0)
    ;; (cl-mpm::update-sim sim)
    (format t "Estimated dt ~E~%" (* dt-scale (cl-mpm/setup:estimate-elastic-dt sim)))
    (format t "Computed dt ~E~%" (* dt-scale (cl-mpm::calculate-min-dt sim)))
    (setf (cl-mpm::sim-dt-scale sim) dt-scale)
    (setf (cl-mpm:sim-mass-scale sim) 1d0)
    (setf (cl-mpm:sim-damping-factor sim) 0d0)
    (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm::calculate-min-dt sim)))
    ;(setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm/setup:estimate-elastic-dt sim)))
    (defparameter *total-iter* 0)
    (let ((quasi-static-solver (class-of sim))
          (load-steps conv-load-steps)
          (i 0))
      (let ((temp-sim (copy-instance sim)))
        (change-class temp-sim elastic-solver)

        (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep temp-sim) 0d0)
        (setf (cl-mpm::sim-dt temp-sim) 0d0)
        (setf (cl-mpm::sim-velocity-algorithm temp-sim) :QUASI-STATIC)
        (set-mp-plastic-damage temp-sim :enable-plastic nil :enable-damage nil)
        ;; find initial quasi-static formation
        (cl-mpm/dynamic-relaxation:converge-quasi-static
         temp-sim
         :energy-crit conv-criteria
         :oobf-crit conv-criteria
         :kinetic-damping nil
         :dt-scale dt-scale
         :conv-steps 10000
         :substeps substeps
         :damping-factor (sqrt 2d0)
         :post-iter-step
         (lambda (i e o)
           (save-conv-step temp-sim output-dir *total-iter* 0 0d0 o e)
           (incf *total-iter* substeps)
           (save-vtks temp-sim output-dir i "conv")
           (funcall plotter temp-sim)))
        (cl-mpm::finalise-loadstep temp-sim)
        (cl-mpm::reset-grid (cl-mpm:sim-mesh temp-sim) :reset-displacement t))

      (cl-mpm:iterate-over-mps
       (cl-mpm:sim-mps sim)
       (lambda (mp)
         (when (typep mp 'cl-mpm/particle::particle-damage)
           (setf (cl-mpm/particle::mp-enable-damage mp) enable-damage))
         (when (typep mp 'cl-mpm/particle::particle-plastic)
           (setf (cl-mpm/particle::mp-enable-plasticity mp) enable-plastic))))

      (funcall post-conv-step sim)
      (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim) dt)
      (funcall setup-quasi-static sim)
      (let* ((current-adaptivity 0)
             (easy-step-counter 0)
             (elastic-dt (if elastic-dt elastic-dt (cl-mpm/setup::estimate-elastic-dt sim)))
             (max-steps (floor total-time (/ dt (expt adaption-constant max-adaptive-steps))))
             )
        (cl-mpm:sim-format sim t "Elastic dt ~E, override quasi-static at ~E~%" elastic-dt (* elastic-dt elastic-dt-margin))
        (loop for step from 1 to (+ 1 max-steps)
              while (and (cl-mpm::sim-run-sim sim)
                         (< (cl-mpm::sim-time sim) total-time))
              do
                 (let ((quasi-conv nil)
                       )
                   (cl-mpm:sim-format sim t "Quasi-timestep ~D, dt refine ~D - dt ~E~%" step current-adaptivity
                                      (/ dt (expt adaption-constant current-adaptivity)))
                   (defparameter *trial-iter* 0)
                   (let ((stagger-iters 0))
                     (loop for i from 0 to (- max-adaptive-steps min-adaptive-steps)
                           while (not quasi-conv)
                           do (progn
                                (let* ((adapted-dt (/ dt (expt adaption-constant current-adaptivity))))
                                  (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim)
                                        (+ (min (- total-time (cl-mpm::sim-time sim)) adapted-dt) 1d-15))
                                  (cl-mpm:sim-format sim t "Trial step ~D, dt refine ~D - dt ~E ~%" i current-adaptivity (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim)))
                                (when (<= (sim-dt-loadstep sim) (* elastic-dt-margin elastic-dt))
                                  (cl-mpm:sim-format sim t "Quasi-time terminated as we got within ~E of the elastic dt~%" elastic-dt-margin)
                                  (setf quasi-conv nil)
                                  (loop-finish))

                                (setf *trial-iter* i)
                                (multiple-value-bind (conv inc-steps)
                                    (step-quasi-time sim step
                                                     :total-steps *total-iter*
                                                     :output-dir output-dir
                                                     :dt-scale dt-scale
                                                     :substeps substeps
                                                     :sub-conv-steps sub-conv-steps
                                                     :enable-damage enable-damage
                                                     :enable-plastic enable-plastic
                                                     :conv-criteria conv-criteria
                                                     :conv-criteria-damage conv-criteria
                                                     :max-damage-inc max-damage-inc
                                                     :max-plastic-inc max-plastic-inc
                                                     :max-deformation-gradient max-deformation-gradient
                                                     :stagger-damage nil
                                                     :save-vtk-dr save-vtk-dr)
                                  (setf quasi-conv conv
                                        stagger-iters inc-steps)
                                  (cl-mpm:sim-format sim t "Quasi-conv? ~A~%" quasi-conv)
                                  (unless quasi-conv
                                    (setf easy-step-counter 0)
                                    (when (= current-adaptivity max-adaptive-steps)
                                      (cl-mpm:sim-format sim t "Quasi-time terminated as too many dt refinemets are required~%"))
                                    (when (<= (sim-dt-loadstep sim) (* elastic-dt-margin elastic-dt))
                                      (cl-mpm:sim-format sim t "Quasi-time terminated as we got within ~E of the elastic dt~%" elastic-dt-margin))
                                    (when (or (= current-adaptivity max-adaptive-steps)
                                              (<= (sim-dt-loadstep sim) (* elastic-dt-margin elastic-dt)))
                                      (loop-finish))
                                    (incf current-adaptivity)
                                    (when (> i 4)
                                      (cl-mpm:sim-format sim t "We've adaptive multiple sucessive times, double jump dt adaption~%")
                                      (incf i 1)))))
                           finally (progn
                                     (cl-mpm:sim-format sim t "Finished with ~D dt adaptions - stagger iters ~D - conv ~A~%" (- i 1) stagger-iters quasi-conv)
                                     (when (> min-damage-inc 0d0)
                                       (format t "End of step max damage inc ~E~%" (damage-increment-criteria sim)))
                                     (if (and (= i 1)
                                              quasi-conv
                                              (if (> min-damage-inc 0d0)
                                                  (< (damage-increment-criteria sim) min-damage-inc)
                                                  t))
                                         (progn
                                           (incf easy-step-counter)
                                           (cl-mpm:sim-format sim t "Potential adaption easy steps ~A~%" easy-step-counter)
                                           (when (>= easy-step-counter 8)
                                             (setf current-adaptivity
                                                   (max min-adaptive-steps
                                                        (- current-adaptivity 1))))
                                           ;;   (setf current-adaptivity
                                           ;;         (max min-adaptive-steps
                                           ;;              (- current-adaptivity 1))))
                                           )
                                         (progn
                                           (setf easy-step-counter 0)
                                           )))))
                   (unless quasi-conv
                     ;;We've adapted down to a min
                     (cl-mpm:sim-format sim t "Start real-timestepping~%")
                     (let ((dt-loadstep (* 1d0 (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim))))
                       (cl-mpm/dynamic-relaxation::reset-mp-velocity sim)
                       (change-class sim explicit-dynamic-solver)
                       (setf (cl-mpm::sim-velocity-algorithm sim) :FLIP)
                       (funcall setup-dynamic sim)
                       (step-real-time sim step
                                       :output-dir output-dir
                                       :plotter plotter
                                       :dt-scale explicit-dt-scale
                                       :criteria explicit-conv-criteria
                                       :damping-0 damping-0
                                       :damping explicit-damping-factor
                                       :target-time (* 1d0 dt-loadstep)
                                       ;:target-time (* 0.1d0 elastic-dt explicit-dt-scale)
                                       :step-time (* dt-loadstep adaption-constant)
                                       :enable-mass-scaling explicit-mass-scaling
                                       :enable-damage enable-damage
                                       :enable-plastic enable-plastic)
                       (change-class sim quasi-static-solver)
                       (cl-mpm:sim-format sim t "Finished real-timestepping~%")
                       (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
                       ;;De-refine our timestepping
                       (decf current-adaptivity)
                       (cl-mpm/dynamic-relaxation::reset-mp-velocity sim))
                     (funcall setup-quasi-static sim))

                   (funcall plotter sim)
                   ;; (plot sim)
                   ;; (vgplot:title (format nil "Step ~D" step))
                   (save-vtks sim output-dir step)
                   ;; (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" step)) sim)
                   ;; (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) sim)
                   ;; (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d.vtk" step)) sim)
                   (swank.live:update-swank)))))
    )
  (cl-mpm:sim-format sim t "Finished algorithm~%")
  )
