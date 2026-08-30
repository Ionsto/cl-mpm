(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)

(defun step-quasi-time (sim
                        global-step
                        &key (dt-scale 1d0)
                          (substeps 50)
                          (total-steps 0)
                          (damping 1d0)
                          (sub-conv-steps 100)
                          (conv-criteria 1d-3)
                          (conv-criteria-damage 1d-3)
                          (output-dir "./output/")
                          (save-vtk-dr t)
                          (enable-damage t)
                          (enable-plastic t)
                          (max-damage-inc 0.6d0)
                          (max-plastic-inc 10d0)
                          (max-deformation-gradient 10d0)
                          (stagger-damage nil)
                          (plotter (lambda (sim))))
  (let ((total-i 0))
    (handler-bind
        ((cl-mpm/errors:error-simulation
           (lambda (c)
             (format t "Handled error~%")
             (incf *total-iter* substeps)
             (save-conv-step sim output-dir *total-iter* global-step 0d0 (cl-mpm::sim-stats-oobf sim) 0d0)
             (when save-vtk-dr
               (save-vtks-dr-step sim output-dir global-step *trial-iter* total-i))
             (incf *total-iter* substeps)
             (cl-mpm/utils::kill-errors)
             (princ c)
             (cl-mpm::reset-loadstep sim)
             (return-from step-quasi-time (values nil 0)))))
        (progn
          (let* ((oobf-crit   conv-criteria)
                 (energy-crit 1d0)
                 ;; (crit (cl-mpm/dynamic-relaxation::sim-convergence-critera sim))
                 (damage-crit conv-criteria-damage)
                 (dconv (if enable-damage damage-crit 0d0))
                 (stagger-iters 0))
            (when (equal (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
              (reset-mp-velocity sim))
            (set-mp-plastic-damage sim :enable-damage enable-damage :enable-plastic enable-plastic)
            (setf (cl-mpm:sim-enable-damage sim) nil)
            (format t "Checkpoint passed~%")
            (setf (cl-mpm::sim-stats-oobf sim) oobf-crit)
            (let ((alt-conv-crit nil))
              (loop for ac from 0 to 10
                    while (not alt-conv-crit)
                    do (progn
                         (loop for stagger-i from 0 to 100
                               while (or
                                      ;; (<= stagger-i 1)
                                      (>= dconv damage-crit)
                                      (>= (cl-mpm::sim-stats-oobf sim) oobf-crit))
                               do
                                  (progn
                                    (format t "Stagger iter ~D ~A~%" stagger-i (cl-mpm:sim-enable-damage sim))
                                    (refine-mesh sim)
                                    (cl-mpm/dynamic-relaxation:converge-quasi-static
                                     sim
                                     :oobf-crit oobf-crit
                                     :energy-crit energy-crit
                                     :dt-scale dt-scale
                                     :substeps substeps
                                     :conv-steps sub-conv-steps
                                     :convergance-criteria
                                     (lambda (sim f o)
                                       (setf dconv (compute-damage-delta sim))
                                       (and
                                        (convergence-criteria sim)
                                        (<= o (cl-mpm/dynamic-relaxation::sim-convergence-critera sim))
                                        (< dconv damage-crit)))
                                     :damping-factor damping
                                     :post-iter-step
                                     (lambda (i e o)
                                       (format t "Updated damage inside of non-stagger ~A~%" (cl-mpm:sim-enable-damage sim))
                                       (funcall plotter sim)
                                       (convergence-check sim)
                                       (check-deformation-gradient sim :max-deformation-gradient max-deformation-gradient)
                                       (check-damage-increment sim :max-damage-inc max-damage-inc)
                                       (check-plastic-increment sim :max-plastic-inc max-plastic-inc)
                                       (incf total-i)
                                       (save-conv-step sim output-dir *total-iter* global-step 0d0 (cl-mpm::sim-stats-oobf sim) 0d0)
                                       (when save-vtk-dr
                                         (save-vtks-dr-step sim output-dir global-step *trial-iter* total-i))
                                       (incf *total-iter* substeps)
                                       (when (cl-mpm::sim-enable-damage sim)
                                         (setf dconv (compute-damage-delta sim))
                                         (format t "d-conv subiter ~E~%" dconv))

                                       ;; (when (and enable-damage
                                       ;;            (not stagger-damage)
                                       ;;            (> stagger-i 1))
                                       ;;   (let ((dconv-1 (compute-damage-delta sim)))
                                       ;;     (setf dconv dconv-1)
                                       ;;     (cl-mpm:sim-format sim t "subiter-step ~D/~D - d-conv ~E ~A~%" stagger-i 0 dconv stagger-damage)
                                       ;;     (loop for d from 1 to 100
                                       ;;           when (>= dconv damage-crit)
                                       ;;             do (progn
                                       ;;                  (cl-mpm/damage::calculate-damage
                                       ;;                   sim
                                       ;;                   (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim))
                                       ;;                  (setf dconv (compute-damage-delta sim))
                                       ;;                  (cl-mpm:sim-format sim t "step ~D/~D - d-conv ~E~%" stagger-i d dconv)
                                       ;;                  (save-conv-step sim output-dir *total-iter* global-step 0d0 (cl-mpm::sim-stats-oobf sim) 0d0)
                                       ;;                  (incf *total-iter*)
                                       ;;                  (check-damage-increment sim :max-damage-inc max-damage-inc)
                                       ;;                  ))
                                       ;;     (setf dconv dconv-1)))
                                       (refine-mesh sim)))
                                    (when enable-damage
                                      (setf (cl-mpm:sim-enable-damage sim) enable-damage)
                                      (cl-mpm/damage::calculate-damage
                                       sim
                                       (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim))
                                      (let ((dconv-1 (compute-damage-delta sim)))
                                        (setf dconv dconv-1)
                                        (cl-mpm:sim-format sim t "step ~D/~D - d-conv ~E ~A~%" stagger-i 0 dconv stagger-damage)
                                        (when stagger-damage
                                          (format t "Staggering damage~%")
                                          (loop for d from 1 to 100
                                                when (>= dconv damage-crit)
                                                  do (progn
                                                       (cl-mpm/damage::calculate-damage
                                                        sim
                                                        (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim))
                                                       (setf dconv (compute-damage-delta sim))
                                                       (cl-mpm:sim-format sim t "step ~D/~D - d-conv ~E~%" stagger-i d dconv)
                                                       (save-conv-step sim output-dir *total-iter* global-step 0d0 (cl-mpm::sim-stats-oobf sim) 0d0)
                                                       (incf *total-iter*)
                                                       (check-damage-increment sim :max-damage-inc max-damage-inc)
                                                       ))
                                          (setf (cl-mpm:sim-enable-damage sim) nil))
                                        (cl-mpm:update-sim sim)
                                        (cl-mpm::update-dynamic-stats sim)
                                        (setf dconv dconv-1)))
                                    (when save-vtk-dr
                                      (save-vtks-dr-step sim output-dir global-step *trial-iter* total-i))
                                    (incf stagger-iters)))
                         (setf alt-conv-crit (convergence-check sim))
                         (unless alt-conv-crit
                           (setf (cl-mpm::sim-stats-oobf sim) oobf-crit)))))
            (when (or (> (cl-mpm::sim-stats-oobf sim) oobf-crit)
                      (> dconv damage-crit))
              (cl-mpm:sim-format sim t "Staggered solve didn't converge ~E ~E~%" dconv (cl-mpm::sim-stats-oobf sim))
              (error (make-instance 'non-convergence-error
                                    :text "Staggered solve didn't converge ~E ~E"
                                    :ke-norm  dconv
                                    :oobf-norm  (cl-mpm::sim-stats-oobf sim))))
            (when save-vtk-dr
              (save-vtks-dr-step sim output-dir global-step *trial-iter* total-i))
            (cl-mpm::finalise-loadstep sim)
            (save-timestep sim output-dir global-step :QUASI-STATIC)
            (values t stagger-iters)))
      ))
  )

(defun run-quasi-time (sim
                       &key (output-dir "./output/")
                         (post-conv-step (lambda (sim)))
                         (post-load-step (lambda (sim)))
                         (plotter (lambda (sim)))
                         (dt-scale 0.5d0)
                         (dt 1d0)
                         (total-time 1d0)
                         (substeps 50)
                         (sub-conv-steps 50)
                         (enable-damage t)
                         (enable-plastic t)
                         (max-adaptive-steps 5)
                         (min-adaptive-steps -1)
                         (adaption-constant 2)
                         (adaption-easy-steps 2)
                         (save-vtk-loadstep t)
                         (save-vtk-conv t)
                         (save-vtk-dr t)
                         (damping (sqrt 2d0))
                         (max-damage-inc 0.5d0)
                         (max-plastic-inc nil)
                         (max-deformation-gradient 10d0)
                         (elastic-solver 'mpm-sim-quasi-static)
                         (initial-quasi-static t)
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
    (setf (cl-mpm::sim-dt-scale sim) dt-scale)
    (setf (cl-mpm:sim-mass-scale sim) 1d0)
    (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm/setup:estimate-elastic-dt sim)))
    (setf (cl-mpm:sim-damping-factor sim) 0d0)
    (setf (cl-mpm:sim-enable-damage sim) nil)
    (defparameter *total-iter* 0)
    (when initial-quasi-static

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
         :damping-factor damping
         :post-iter-step
         (lambda (i e o)
           (save-conv-step temp-sim output-dir *total-iter* 0 0d0 o e)
           (incf *total-iter* substeps)
           (when save-vtk-conv
             (save-vtks temp-sim output-dir i "conv"))
           (funcall plotter temp-sim)))
        (cl-mpm::finalise-loadstep temp-sim)
        (cl-mpm::reset-grid (cl-mpm:sim-mesh temp-sim) :reset-displacement t)))
    (let ()
      (set-mp-plastic-damage sim :enable-plastic enable-plastic :enable-damage enable-damage)
      (cl-mpm:sim-format sim t "Call post-conv~%")
      (funcall post-conv-step sim)
      (cl-mpm:sim-format sim t "Start iter~%")
      (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim) dt)
      (let* ((current-adaptivity 0)
             (max-steps (floor total-time (/ dt (expt adaption-constant max-adaptive-steps))))
             (easy-step-counter 0)
             (sim-time 0d0))
        (loop for step from 1 to (+ 1 max-steps)
              while (and (cl-mpm::sim-run-sim sim)
                         (< sim-time total-time))
              do
                 (let ((quasi-conv nil))
                   (cl-mpm:sim-format sim t "quasi-timestep ~d - time ~E, dt refine ~d~%" step sim-time current-adaptivity)
                   (defparameter *trial-iter* 0)
                   (let ((stagger-iters 0))
                     (loop for i from 0 to (- max-adaptive-steps current-adaptivity)
                           while (not quasi-conv)
                           do (progn
                                (let* ((adapted-dt (/ dt (expt adaption-constant current-adaptivity))))
                                  (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim)
                                        (+ (min (- total-time sim-time) adapted-dt) 1d-15))
                                  (cl-mpm:sim-format sim t "trial step ~d, dt refine ~d - dt ~E~%" i current-adaptivity (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim))
                                  (setf *trial-iter* i)
                                  ;; (let (;; (conv-crit 1d-3)
                                  ;;       (conv-crit 1d-3)
                                  ;;       (residual-normaliser nil)
                                  ;;       (substeps 20))
                                  ;;   (setf (sim-dt-loadstep sim) (* 1d0 dt))
                                  ;;   (generalised-staggered-solve
                                  ;;    sim
                                  ;;    :crit conv-crit
                                  ;;    :substeps substeps
                                  ;;    :sub-conv-steps 100
                                  ;;    :dt-scale 0.9d0
                                  ;;    :damping (sqrt 2d0)
                                  ;;    :max-damage-inc 0.1d0
                                  ;;    :max-plastic-inc nil;
                                  ;;    :post-iter-step (lambda (i e o)
                                  ;;                      (format t "Dynamic substep ~D~%" i)
                                  ;;                      ;; (save-conv-step sim "./output/" *total-iter* *total-step* 0d0 o e)
                                  ;;                      ))


                                  ;;   (setf (cl-mpm::sim-dt-scale sim) dt-scale)
                                  ;;   (setf (cl-mpm::sim-dt sim) adapted-dt)
                                  ;;   (setf (cl-mpm::sim-damping-factor sim) 0d0)
                                  ;;   (cl-mpm::finalise-loadstep sim)
                                  ;;   (save-timestep sim output-dir step :QUASI-STATIC)
                                  ;;   (setf quasi-conv t))
                                  )
                                (multiple-value-bind (conv inc-steps)
                                    (step-quasi-time sim step
                                                     :total-steps *total-iter*
                                                     :plotter plotter
                                                     :output-dir output-dir
                                                     :dt-scale dt-scale
                                                     :substeps substeps
                                                     :sub-conv-steps sub-conv-steps
                                                     :enable-damage enable-damage
                                                     :damping damping
                                                     :save-vtk-dr save-vtk-dr
                                                     :conv-criteria conv-criteria
                                                     :conv-criteria-damage conv-criteria
                                                     :max-damage-inc max-damage-inc
                                                     :max-plastic-inc max-plastic-inc
                                                     :max-deformation-gradient max-deformation-gradient
                                                     :enable-plastic enable-plastic
                                                     ;; :stagger-damage t
                                                     )
                                  (setf quasi-conv conv
                                        stagger-iters inc-steps)
                                  (format t "Step ~A current adaptivity ~D~%" quasi-conv current-adaptivity)
                                  (unless quasi-conv
                                    (setf easy-step-counter 0)
                                    (if (= current-adaptivity max-adaptive-steps)
                                        (progn
                                          (cl-mpm::reset-node-displacement sim)
                                          (cl-mpm:sim-format sim t "quasi-time terminated as too many dt refinemets are required~%")
                                          (loop-finish))
                                        (incf current-adaptivity))))
                                )
                           finally (progn
                                     (cl-mpm:sim-format sim t "finished with ~d dt adaptions - stagger iters ~d~%" (- i 1) stagger-iters)
                                     (when (= i 1)
                                       (incf easy-step-counter))
                                     (cl-mpm:sim-format sim t "easy steps ~D~%" easy-step-counter)
                                     (when (and (= i 1)
                                                (> easy-step-counter adaption-easy-steps)
                                                ;; (< stagger-iters 4)
                                                )
                                       (setf current-adaptivity
                                             (max min-adaptive-steps
                                                  (- current-adaptivity 1)))))))
                   (unless quasi-conv
                     (setf (cl-mpm::sim-run-sim sim) nil
                           result nil))
                   (incf sim-time (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim))
                   (funcall plotter sim)
                   (funcall post-load-step sim)

                   (when save-vtk-loadstep
                     (save-vtks sim output-dir step))
                   ;; (when save-vtk-loadstep
                   ;;   (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" step)) sim)
                   ;;   (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) sim)
                   ;;   (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d.vtk" step)) sim))
                   (swank.live:update-swank)))))
    result))
