(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)
(defparameter *total-iter* 0)

(defun set-mp-plastic-damage (sim &key (enable-damage t) (enable-plastic t))
  (cl-mpm:iterate-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (when (typep mp 'cl-mpm/particle::particle-damage)
       (setf (cl-mpm/particle::mp-enable-damage mp) enable-damage))
     (when (typep mp 'cl-mpm/particle::particle-plastic)
       (setf (cl-mpm/particle::mp-enable-plasticity mp) enable-plastic))))
  )

(defun save-timestep-preamble (output-dir)
  (with-open-file (stream (merge-pathnames output-dir "./timesteps.csv") :direction :output :if-exists :supersede)
    (format stream "steps,time,damage,plastic,energy,oobf,work,step-type,mass~%")))

(defun save-timestep (sim output-dir step type)
  (when (uiop:file-exists-p (merge-pathnames output-dir "timesteps.csv"))
    (with-open-file (stream (merge-pathnames "timesteps.csv" output-dir) :direction :output :if-exists :append)
      (format stream "~D,~f,~f,~f,~f,~f,~f,~A,~f~%"
              step
              (cl-mpm::sim-time sim)
              (get-damage sim)
              (get-plastic sim)
              (cl-mpm::sim-stats-energy sim)
              (cl-mpm::sim-stats-oobf sim)
              (cl-mpm::sim-stats-work sim)
              type
              0d0))))

(defun save-conv-preamble (output-dir)
  (with-open-file (stream (merge-pathnames output-dir "conv.csv") :direction :output :if-exists :supersede)
    (format stream "iter,step,real-time,plastic,damage,oobf,energy~%")))

(defun save-conv-step (sim output-dir total-iter step real-time oobf energy)
  (when (uiop:file-exists-p (merge-pathnames output-dir "conv.csv"))
    (with-open-file (stream (merge-pathnames output-dir "conv.csv") :direction :output :if-exists :append)
      (format stream "~D,~D,~f,~f,~f,~f,~f~%" total-iter step real-time (get-plastic sim) (get-damage sim)
              (if (sb-ext:float-infinity-p oobf) 0d0 oobf) energy))))

(defun save-conv (sim output-dir iter)
  (let ((oobf (cl-mpm::sim-stats-oobf sim))
        (energy (/ (cl-mpm::sim-stats-energy sim) (cl-mpm::sim-stats-work sim)))
        (real-time (cl-mpm::sim-time sim)))
    (with-open-file (stream (merge-pathnames output-dir "conv.csv") :direction :output :if-exists :append)
      (format stream "~D,~D,~f,~f,~f,~f,~f~%" iter 0 real-time (get-plastic sim) (get-damage sim)
              oobf energy))))

(defun get-damage (sim)
  (lparallel:pmap-reduce
   (lambda (mp)
     (if (typep mp 'cl-mpm/particle::particle-damage)
         (*
          (cl-mpm/particle::mp-mass mp)
          (cl-mpm/particle:mp-damage mp))
         0d0))
   #'+
   (cl-mpm:sim-mps sim)))
(defun get-plastic (sim)
  (lparallel:pmap-reduce
   (lambda (mp)
     (if (typep mp 'cl-mpm/particle::particle-plastic)
         (*
          (cl-mpm/particle::mp-mass mp)
          (cl-mpm/particle::mp-strain-plastic-vm mp))
         0d0))
   #'+
   (cl-mpm:sim-mps sim)))

(defun true-intertial-criteria (sim loadstep-dt)
  (let ((energy 0d0)
        (mass 0d0)
        (lock (sb-thread:make-mutex)))
    (cl-mpm:iterate-over-nodes
     (cl-mpm:sim-mesh sim)
     (lambda (n)
       (when (cl-mpm/mesh:node-active n)
         (sb-thread:with-mutex (lock)
           (incf mass (cl-mpm/mesh:node-mass n))
           (incf energy
                 (*
                  ;; (cl-mpm/mesh::node-true-mass n)
                  (cl-mpm/mesh::node-true-mass n)
                  (cl-mpm/fastmaths::mag-squared
                   (cl-mpm/fastmaths:fast-scale-vector
                    (cl-mpm/mesh::node-displacment n) (expt (/ 1d0 loadstep-dt) 1)))))))))
    (let ((power (cl-mpm::sim-stats-power sim)))
      (if (= power 0d0)
          0d0
          (/
           energy
           ;; (if (= mass 0d0)
           ;;     0d0
           ;;     (/ energy mass))
           power)))))

(defun damage-increment-criteria (sim &key (criteria 0.5d0))
  (let ((lock (sb-thread:make-mutex))
        (result nil))
    (cl-mpm:iterate-over-mps
     (cl-mpm:sim-mps sim)
     (lambda (mp)
       (when (typep mp 'cl-mpm/particle::particle-damage)
         (when (> (cl-mpm/particle::mp-damage-increment mp) criteria)
           (setf result t)))))
    result))

(defun generalised-staggered-solve (sim &key
                                          (crit 1d-3)
                                          (dt-scale 1d0)
                                          (substeps 10)
                                          (post-iter-step (lambda (i e o)))
                                          (enable-damage t)
                                          (enable-plastic t)
                                          )
  (let* ((damage-prev (get-damage sim))
         (damage damage-prev)
         (oobf-crit   crit)
         (energy-crit crit)
         (damage-crit crit)
         (dconv damage-crit)
         (total-i 0)
         (stagger-iters 0)
         )
    ;; (set-mp-plastic-damage sim :enable-damage enable-damage :enable-plastic enable-plastic)
    (loop for stagger-i from 0 to 10
                while (or (>= dconv damage-crit))
                do
                   (progn
                     (setf (cl-mpm:sim-enable-damage sim) nil)
                     (cl-mpm/dynamic-relaxation:converge-quasi-static
                      sim
                      :oobf-crit oobf-crit
                      :energy-crit energy-crit
                      :dt-scale dt-scale
                      :substeps substeps
                      :conv-steps 200
                      :damping-factor 1d0
                      :post-iter-step post-iter-step)
                     (if enable-damage
                         (progn
                           (setf (cl-mpm:sim-enable-damage sim) t)
                           (when (typep sim 'cl-mpm/damage::mpm-sim-damage)
                             (cl-mpm/damage::calculate-damage sim (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim)))
                           (setf damage (get-damage sim))
                           (setf dconv (if (> damage 0d0)
                                           (if (> damage-prev 0d0)
                                               (/ (- damage damage-prev) damage-prev)
                                               sb-ext:double-float-positive-infinity)
                                           0d0))
                           (setf damage-prev damage))
                       (setf dconv 0d0))))))

(declaim (notinline step-quasi-time))
(defun step-quasi-time (sim
                        global-step
                        &key (dt-scale 1d0)
                          (substeps 50)
                          (total-steps 0)
                          (damping 1d0)
                          (conv-criteria 1d-3)
                          (conv-criteria-damage 1d-3)
                          (output-dir "./output/")
                          (enable-damage t)
                          (enable-plastic t)
                          (plotter (lambda (sim))))
  (handler-case
      (progn
        (let* ((damage-prev (get-damage sim))
               (damage damage-prev)
               (oobf-crit   conv-criteria)
               (energy-crit conv-criteria)
               (damage-crit conv-criteria-damage)
               (dconv damage-crit)
               (inertia 0d0)
               (intertia-crit 1d-3)
               (total-i 0)
               (stagger-iters 0)
               )
          (set-mp-plastic-damage sim :enable-damage enable-damage :enable-plastic enable-plastic)
          (loop for stagger-i from 0 to 10
                while (or (>= dconv damage-crit))
                do
                   (progn
                     ;; (setf (cl-mpm:sim-enable-damage sim) t)
                     (cl-mpm/dynamic-relaxation:converge-quasi-static
                      sim
                      :oobf-crit oobf-crit
                      :energy-crit energy-crit
                      :dt-scale dt-scale
                      :substeps substeps
                      :conv-steps 200
                      :damping-factor 1d0
                      :post-iter-step
                      (lambda (i e o)
                        (funcall plotter sim)
                        (vgplot:title (format nil "substep ~D - oobf: ~E" i o))
                        (incf total-i)
                        ;; (when (= (mod i 50) 0)
                        ;;   (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_step_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
                        ;;   (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_step_nodes_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
                        ;;   )
                        (format t "Def crit ~E~%" (compute-max-deformation sim))
                        (when (criteria-deformation-gradient sim :criteria 10d0)
                          (format t "Deformation gradient criteria exceeded~%")
                          (error (make-instance 'non-convergence-error
                                                :text "Deformation gradient J exceeded"
                                                :ke-norm 0d0
                                                :oobf-norm 0d0)))
                        (let ((true-intertia (true-intertial-criteria sim (sim-dt-loadstep sim))))
                          (format t "True intertia ~E~%" true-intertia)
                          (save-conv-step sim output-dir *total-iter* global-step 0d0 o true-intertia)
                          (setf inertia true-intertia))
                        ;; (save-conv-step sim output-dir *total-iter* global-step 0d0 o 0d0)
                        (incf *total-iter* substeps)))

                     (setf (cl-mpm:sim-enable-damage sim) t)
                     (when (typep sim 'cl-mpm/damage::mpm-sim-damage)
                       (cl-mpm/damage::calculate-damage sim (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim)))
                     (setf (cl-mpm:sim-enable-damage sim) nil)
                     (setf damage (get-damage sim))
                     (setf dconv (if (> damage 0d0)
                                     (if (> damage-prev 0d0)
                                         (/ (- damage damage-prev) damage-prev)
                                         sb-ext:double-float-positive-infinity)
                                     0d0))
                     (format t "Damage ~E - prev damage ~E ~%" damage damage-prev)
                     (format t "step ~D - d-conv ~E~%" stagger-i dconv)
                     (setf damage-prev damage)
                     (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_step_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
                     (incf stagger-iters)
                     (when (damage-increment-criteria sim :criteria 0.8d0)
                       (format t "Damage criteria failed~%")
                       (error (make-instance 'non-convergence-error
                                             :text "Damage criteria exeeded"
                                             :ke-norm 0d0
                                             :oobf-norm 0d0)))
                     ;; (when (> inertia intertia-crit)
                     ;;   (error (make-instance 'non-convergence-error
                     ;;                         :text "Quasi-time inertia was too large"
                     ;;                         :ke-norm inertia
                     ;;                         :oobf-norm oobf-crit)))
                     ))
          (when (or (> (cl-mpm::sim-stats-oobf sim) oobf-crit)
                    (> dconv damage-crit))
            (error (make-instance 'non-convergence-error
                                  :text "Staggered solve didn't converge ~E ~E"
                                  :ke-norm  dconv
                                  :oobf-norm  (cl-mpm::sim-stats-oobf sim))))
          (cl-mpm::finalise-loadstep sim)
          (save-timestep sim output-dir global-step :QUASI-STATIC)
          (values t stagger-iters)))
    (error (c)
      (princ c)
      ;; (princ (text c))
      (cl-mpm::reset-loadstep sim)
      (values nil 0)))
  )
(declaim (notinline step-real-time))
(defun step-real-time (sim
                       global-step
                       &key
                         (target-time 1d0)
                         (plotter (lambda (sim)))
                         (dt-scale 0.5)
                         (output-dir "./output/")
                         (max-steps 1000)
                         (mass-scale 1d0)
                         (damping 1d-4)
                         (enable-mass-scaling nil)
                         (enable-damage t)
                         (criteria 1d-3)
                         (enable-plastic t))
  (setf (cl-mpm:sim-mass-scale sim) 1d0)
  (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm/setup:estimate-elastic-dt sim)))
  (setf (cl-mpm:sim-damping-factor sim) (* damping (cl-mpm/setup:estimate-critical-damping sim)))
  (set-mp-plastic-damage sim :enable-damage enable-damage :enable-plastic enable-plastic)
  (setf (cl-mpm:sim-enable-damage sim) enable-damage)
  (let* ((e-crit criteria)
         (oobf-crit criteria)
         (energy e-crit)
         (oobf oobf-crit)
         (work 0d0)
         (intertial-passed nil)
         (dt-0 (* dt-scale (cl-mpm/setup:estimate-elastic-dt sim)))
         (substeps (round target-time (cl-mpm:sim-dt sim))))
    (format t "Substeps ~D~%" substeps)
    (time (loop for step from 0 to max-steps
                while (and (cl-mpm::sim-run-sim sim)
                           (or (>= energy e-crit)
                               (>= oobf oobf-crit)
                               (< step 4)
                               (not intertial-passed)
                               ))
                do
                   (let ((substeps (round target-time (cl-mpm:sim-dt sim))))
                     (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_real_~5,'0d_~5,'0d.vtk" global-step step)) sim)
                     (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_real_nodes_~5,'0d_~5,'0d.vtk" global-step step)) sim)
                     (setf energy 0d0)
                     (setf oobf 0d0)
                     (format t "Real time step ~d - substeps ~d - time ~E ~%" step substeps target-time)
                     (time
                      (dotimes (i substeps)
                        (cl-mpm::update-sim sim)
                        (incf oobf (cl-mpm::sim-stats-oobf sim))
                        (incf energy (cl-mpm::sim-stats-energy sim))
                        (incf work (cl-mpm::sim-stats-power sim))))

                     (setf
                      energy (/ energy substeps)
                      oobf (/ oobf substeps))
                     (if (= work 0d0)
                         (setf energy 0d0)
                         (setf energy (abs (/ energy work))))

                     (when enable-mass-scaling
                       (let ((res (and (< energy 1d-2)
                                       (< oobf 1d-2))))
                         (setf (cl-mpm:sim-mass-scale sim)
                               (if res
                                   1d2
                                   1d0))
                         (format t "Accelerate real time ~A~%" res)))
                     ;; (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm/setup:estimate-elastic-dt sim)))
                     ;; (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time sim target-time :dt-scale dt-scale)
                     ;;   (format t "CFL dt estimate: ~f~%" dt-e)
                     ;;   (format t "CFL step count estimate: ~D~%" substeps-e)
                     ;;   (setf (cl-mpm:sim-dt sim) dt-e)
                     ;;   (setf substeps substeps-e))
                     (when (or (>= energy e-crit)
                               (>= oobf oobf-crit))
                       (setf intertial-passed t))
                     (format t "Residuals ~E ~E ~%" energy oobf)
                     (save-timestep sim output-dir global-step :DYNAMIC)
                     (funcall plotter sim)
                     (swank.live:update-swank))))))

(declaim (notinline run-multi-stage))
(defun run-multi-stage (sim
                        &key (output-dir "./output/")
                          (post-conv-step (lambda (sim)))
                          (plotter (lambda (sim)))
                          (dt-scale 1d0)
                          (explicit-dt-scale 0.9d0)
                          (dt 1d0)
                          (steps 1d0)
                          (max-adaptive-steps 5)
                          (min-adaptive-steps -1)
                          (conv-criteria 1d-3)
                          (substeps 50)
                          (enable-damage t)
                          (enable-plastic t)
                          (save-vtk-dr t)
                          (explicit-dynamic-solver 'cl-mpm::mpm-sim-usf))
  (let ()
    (uiop:ensure-all-directories-exist (list output-dir))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
    (cl-mpm/output::save-simulation-parameters (merge-pathnames output-dir "settings.json")
                                               sim
                                               (list :dt dt
                                                     :criteria-energy conv-criteria
                                                     :criteria-oobf conv-criteria
                                                     :criteria-hist 1d0
                                                     ))
    (save-timestep-preamble output-dir)
    (save-conv-preamble output-dir)
    (setf (cl-mpm::sim-dt-scale sim) dt-scale)
    (setf (cl-mpm:sim-mass-scale sim) 1d0)
    (setf (cl-mpm:sim-damping-factor sim) 0d0)
    (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm/setup:estimate-elastic-dt sim)))
    (defparameter *total-iter* 0)
    (let ((quasi-static-solver (class-of sim))
          (vel-algo :BLEND)
          ;; (vel-algo (cl-mpm::sim-velocity-algorithm sim))
          )
      ;; (change-class sim explicit-dynamic-solver)
      ;; (change-class sim quasi-static-solver)
      (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim) 0d0)
      (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
      (set-mp-plastic-damage sim :enable-plastic nil :enable-damage nil)
      (cl-mpm/dynamic-relaxation:converge-quasi-static
       sim
       :energy-crit conv-criteria
       :oobf-crit conv-criteria
       :kinetic-damping nil
       :dt-scale dt-scale
       :conv-steps 1000
       :substeps substeps
       :damping-factor 1d0;1(sqrt 2)
       :post-iter-step
       (lambda (i e o)
         (save-conv-step sim output-dir *total-iter* 0 0d0 o e)
         (incf *total-iter* substeps)
         (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_conv_~5,'0d.vtk" i)) sim)
         ;; (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_conv_nodes__~5,'0d.vtk" i)) sim)
         ;; (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_conv_cells__~5,'0d.vtk" i)) sim)
         ))
      (cl-mpm::finalise-loadstep sim)
      (setf (cl-mpm::sim-time sim) 0d0)
      (cl-mpm/dynamic-relaxation::reset-mp-velocity sim)
      (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)

      (cl-mpm:iterate-over-mps
       (cl-mpm:sim-mps sim)
       (lambda (mp)
         (when (typep mp 'cl-mpm/particle::particle-damage)
           (setf (cl-mpm/particle::mp-enable-damage mp) enable-damage))
         (when (typep mp 'cl-mpm/particle::particle-plastic)
           (setf (cl-mpm/particle::mp-enable-plasticity mp) enable-plastic))))

      (funcall post-conv-step sim)
      (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim) dt)
      (let* ((current-adaptivity 0)
             (elastic-dt (cl-mpm/setup::estimate-elastic-dt sim))
             (elastic-dt-margin 1000))
        (format t "Elastic dt ~E, override quasi-static at ~E~%" elastic-dt (* elastic-dt elastic-dt-margin))
        (loop for step from 1 to steps
              while (cl-mpm::sim-run-sim sim)
              do
                 (let ((quasi-conv nil))
                   (format t "Quasi-timestep ~D, dt refine ~D - dt ~E~%" step current-adaptivity (/ dt (expt 2 current-adaptivity)))
                   (defparameter *trial-iter* 0)
                   (let ((stagger-iters 0))
                     (loop for i from 0 to max-adaptive-steps
                           while (not quasi-conv)
                           do (progn
                                (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim) (/ dt (expt 2 current-adaptivity)))
                                (format t "Trial step ~D, dt refine ~D~%" i current-adaptivity)
                                (setf *trial-iter* i)
                                (multiple-value-bind (conv inc-steps)
                                    (step-quasi-time sim step
                                                     :total-steps *total-iter*
                                                     :plotter plotter
                                                     :output-dir output-dir
                                                     :dt-scale dt-scale
                                                     :substeps substeps
                                                     :enable-damage enable-damage
                                                     :enable-plastic enable-plastic
                                                     :conv-criteria conv-criteria
                                                     :conv-criteria-damage conv-criteria)
                                  (format t "Quasi-conv? ~A~%" quasi-conv)
                                  (setf quasi-conv conv
                                        stagger-iters inc-steps)
                                  (unless quasi-conv
                                    (when (= current-adaptivity max-adaptive-steps)
                                      (format t "Quasi-time terminated as too many dt refinemets are required~%"))
                                    (when (<= (sim-dt-loadstep sim) (* elastic-dt-margin elastic-dt))
                                      (format t "Quasi-time terminated as we got within ~E of the elastic dt~%" elastic-dt-margin))
                                    (if (or (= current-adaptivity max-adaptive-steps)
                                            (<= (sim-dt-loadstep sim) (* elastic-dt-margin elastic-dt)))
                                        (progn
                                          ;; (format t "Quasi-time terminated as too many dt refinemets are required~%")
                                          (loop-finish))
                                        (incf current-adaptivity)))))
                           finally (progn
                                     (format t "Finished with ~D dt adaptions - stagger iters ~D - conv ~A~%" (- i 1) stagger-iters quasi-conv)
                                     (when (and (= i 1) (<= stagger-iters 3))
                                       (setf current-adaptivity
                                             (max min-adaptive-steps
                                                  (- current-adaptivity 1)))))))
                   (unless quasi-conv
                     ;;We've adapted down to a min
                     (format t "Start real-timestepping~%")
                     (let ((dt-loadstep (* 1d0 (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim))))
                       (cl-mpm/dynamic-relaxation::reset-mp-velocity sim)
                       (change-class sim explicit-dynamic-solver)
                       (setf (cl-mpm::sim-velocity-algorithm sim) :FLIP)
                       ;; (setf (cl-mpm/aggregate::sim-enable-aggregate sim) t)
                       (step-real-time sim step
                                       :output-dir output-dir
                                       :plotter plotter
                                       :dt-scale explicit-dt-scale
                                       :criteria 1d-3
                                       :damping 1d-3
                                       :target-time (* 0.1d0 dt-loadstep)
                                       :enable-damage enable-damage
                                       :enable-plastic enable-plastic)
                       (change-class sim quasi-static-solver)
                       ;; (setf (cl-mpm/aggregate::sim-enable-aggregate sim) t)
                       (format t "Finished real-timestepping~%")
                       (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
                       (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim) dt-loadstep))
                     )

                   (funcall plotter sim)
                   ;; (plot sim)
                   ;; (vgplot:title (format nil "Step ~D" step))
                   (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" step)) sim)
                   (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) sim)
                   (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d.vtk" step)) sim)
                   (swank.live:update-swank)))))
    )
  (format t "Finished algorithm~%")
  )



;; (defun test-multi-step ()
;;   (setup :mps 2 :refine 0.5)
;;   (run-multi-stage
;;    :output-dir "./output/"
;;    :dt 0.5d0
;;    ;; :dt-scale (/ 0.8d0 (sqrt 1d0))
;;    ;; :load-steps 20
;;    )
;;   )


(defun save-test-vtks (sim &key (output-dir "./output/"))
  (cl-mpm/output:save-vtk (merge-pathnames "test.vtk" output-dir) sim)
  (cl-mpm/output:save-vtk-nodes (merge-pathnames "test_nodes_0.vtk" output-dir) sim)
  (cl-mpm/output:save-vtk-nodes (merge-pathnames "test_nodes_1.vtk" output-dir) sim)
  (cl-mpm/output:save-vtk-cells (merge-pathnames "test_cells.vtk" output-dir) sim)
  )


(defun converge-staggered (sim
                           criteria
                           damage-criteria
                           kinetic-damping
                           adaptive-damping
                           dt-scale
                           substeps
                           damping
                           enable-damage
                           enable-plasticity
                           step
                           output-dir
                           &key (save-vtk-dr nil)
                             (plotter (lambda (sim)))
                             )
  (let* ((damage-prev (get-damage sim))
         (damage-eval t)
         (damage damage-prev)
         (max-iter 1000)
         (dconv damage-criteria))
    (setf (cl-mpm::sim-enable-damage sim) nil)
    (cl-mpm/dynamic-relaxation:converge-quasi-static
     sim
     :oobf-crit criteria
     :energy-crit criteria
     :kinetic-damping kinetic-damping
     :damping-factor (if adaptive-damping damping nil)
     :dt-scale dt-scale
     :substeps substeps
     :conv-steps 1000
     :post-iter-step
     (lambda (i energy oobf)
       (incf *total-iter* substeps)
       (save-conv-step sim output-dir *total-iter* step 0d0 oobf energy)))

    (break)
    (loop for i from 0 to max-iter
          while (>= dconv damage-criteria)
          do (progn
               (when (typep sim 'cl-mpm/damage::mpm-sim-damage)
                 (setf (cl-mpm::sim-enable-damage sim) damage-eval)
                 (cl-mpm/damage::calculate-damage sim (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim))
                 (setf damage (get-damage sim)))
               (if (= damage damage-prev)
                   (setf dconv 0d0)
                   (when (> damage damage-prev)
                     (setf dconv (if (> damage-prev 0d0) (/ (- damage damage-prev) damage-prev) sb-ext:double-float-positive-infinity))))
               (setf damage-prev damage)
               (setf (cl-mpm::sim-enable-damage sim) nil)
               (cl-mpm/dynamic-relaxation:converge-quasi-static
                sim
                :oobf-crit criteria
                :energy-crit criteria
                :kinetic-damping kinetic-damping
                :damping-factor (if adaptive-damping damping nil)
                :dt-scale dt-scale
                :substeps substeps
                :conv-steps 1000
                :post-iter-step
                (lambda (i energy oobf)
                  (incf *total-iter* substeps)
                  (save-conv-step sim output-dir *total-iter* step 0d0 oobf energy))
                )
               (funcall plotter sim)
               (when save-vtk-dr
                 (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d_~5,'0d.vtk" step i)) sim)
                 (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d_~5,'0d.vtk" step i)) sim)
                 (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d_~5,'0d.vtk" step i)) sim))

               ))))


(defun run-load-control (sim
                         &key (output-dir "./output/")
                           (loading-function nil)
                           (load-steps 10)
                           (substeps 50)
                           (damping 1d0)
                           (kinetic-damping nil)
                           (adaptive-damping t)
                           (criteria 1d-3)
                           (post-conv-step (lambda (sim)))
                           (post-iter-step (lambda (sub-iter oobf energy)))
                           (pre-step (lambda ()))
                           (plotter (lambda (sim)))
                           (save-vtk-dr t)
                           (save-vtk-loadstep t)
                           (dt-scale 1d0))
  (uiop:ensure-all-directories-exist (list output-dir))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (save-conv-preamble output-dir)
  (cl-mpm/output::save-simulation-parameters (merge-pathnames output-dir "settings.json") sim)
  (funcall pre-step)
  ;; (save-conv-step sim output-dir 0 0 0d0 0d0)
  (defparameter *total-iter* 0)
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (let* ((load (cl-mpm:sim-gravity sim)))
      (unless loading-function
        (setf loading-function (lambda (percent)
                                 (format t "Loading factor ~E~%" percent)
                                 (setf (cl-mpm:sim-gravity sim)
                                       (* load percent)))))
      (setf (cl-mpm::sim-dt-scale sim) dt-scale)
      (let ((t0 (get-internal-real-time))
            (tunits internal-time-units-per-second)
            )
        (loop for step from 1 to load-steps
              while (cl-mpm::sim-run-sim sim)
              do
                 (progn
                   (format t "Load step ~D~%" step)
                   (defparameter *ke-last* 0d0)
                   ;; (pprint step)
                   (funcall loading-function (/ (float step) load-steps))
                   (let ((conv-steps 0)
                         (i 0))
                     (when damping
                       (setf (cl-mpm::sim-damping-factor sim)
                             (* damping (cl-mpm/setup:estimate-critical-damping sim))))
                     (time
                      (cl-mpm/dynamic-relaxation:converge-quasi-static
                       sim
                       :oobf-crit criteria
                       :energy-crit criteria
                       :kinetic-damping kinetic-damping
                       :damping-factor (if adaptive-damping damping nil)
                       :dt-scale dt-scale
                       :substeps substeps
                       :conv-steps 1000
                       :post-iter-step
                       (lambda (i energy oobf)
                         (funcall post-iter-step i energy oobf)
                         (setf conv-steps (* substeps i))
                         (vgplot:title (format nil "Step ~D - substep ~D - KE ~,3E - OOBF ~,3E - damp ~,3E"  step (* i substeps) energy oobf
                                               (/ (cl-mpm::sim-damping-factor sim)
                                                  (cl-mpm/setup:estimate-critical-damping sim))))
                         (format t "Substep ~D~%" i)
                         (let ((i (+ 0 i)))
                           (when save-vtk-dr
                             (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_step_~5,'0d_~5,'0d.vtk" step i)) sim)
                             (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_step_nodes_~5,'0d_~5,'0d.vtk" step i)) sim)
                             (cl-mpm/penalty:save-vtk-penalties (merge-pathnames output-dir (format nil "sim_step_p_~5,'0d_~5,'0d.vtk" step i)) sim)

                             ;; (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_step_cells_~5,'0d_~5,'0d.vtk" step i)) sim)
                             ))

                         (incf *total-iter* substeps)
                         (save-conv-step sim output-dir *total-iter* step 0d0 oobf energy)
                         ;; (incf total-iter substeps)
                         ;; (save-conv-step sim output-dir total-iter step (/ (- (get-internal-real-time) t0) tunits) oobf energy)
                         ;; (setf t0 (get-internal-real-time))
                         (funcall plotter sim)
                         )))
                     (vgplot:title (format nil "Step ~D - ~D" step conv-steps))
                     ))
                 (funcall post-conv-step sim)
                 (when save-vtk-loadstep
                   (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) sim)
                   (cl-mpm/output::save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d.vtk" step)) sim))
                 (cl-mpm::finalise-loadstep sim)
                 (funcall plotter sim)
                 (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step))
                                    :terminal "png size 1920,1080"
                                    )
                 (when save-vtk-loadstep
                   (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" step)) sim)
                   (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim_p_~5,'0d.vtk" step)) sim ))
                 (sleep 0.1d0)
                 (swank.live:update-swank)
              )))))

;; (defun run-strength-reduction (sim
;;                                &key (output-dir "./output/")
;;                                  (loading-function nil)
;;                                  (load-steps 10)
;;                                  (substeps 50)
;;                                  (damping 1d-1)
;;                                  (kinetic-damping t)
;;                                  (adaptive-damping t)
;;                                  (criteria 1d-3)
;;                                  (post-conv-step (lambda (sim)))
;;                                  (plotter (lambda (sim)))
;;                                  (save-vtk-dr t)
;;                                  (save-vtk-loadstep t)
;;                                  (dt-scale 0.5d0))
;;   (uiop:ensure-all-directories-exist (list output-dir))
;;   (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
;;   (save-conv-preamble output-dir)
;;   ;; (save-conv-step sim output-dir 0 0 0d0 0d0)
;;   (with-accessors ((mps cl-mpm:sim-mps))
;;       sim
;;     (let* ((load (cl-mpm/particle::mp-gravity (aref (cl-mpm:sim-mps sim) 0)))
;;            (total-iter 0))
;;       (unless loading-function
;;         (setf loading-function (lambda (percent)
;;                                  (format t "Loading factor ~E~%" percent)
;;                                  (cl-mpm:iterate-over-mps
;;                                   mps
;;                                   (lambda (mp)
;;                                     (setf (cl-mpm/particle:mp-gravity mp) (* load percent)))))))
;;       (setf (cl-mpm::sim-dt-scale sim) dt-scale)
;;       (let ((t0 (get-internal-real-time))
;;             (tunits internal-time-units-per-second)
;;             )
;;         (loop for step from 1 to load-steps
;;               while (cl-mpm::sim-run-sim sim)
;;               do
;;                  (progn
;;                    (format t "Load step ~D~%" step)
;;                    (defparameter *ke-last* 0d0)
;;                    ;; (pprint step)
;;                    (funcall loading-function (/ (float step) load-steps))
;;                    (let ((conv-steps 0)
;;                          (i 0))
;;                      (when damping
;;                        (setf (cl-mpm::sim-damping-factor sim)
;;                              (* damping (cl-mpm/setup:estimate-critical-damping sim))))
;;                      (time
;;                       (cl-mpm/dynamic-relaxation:converge-quasi-static
;;                        sim
;;                        :oobf-crit criteria
;;                        :energy-crit criteria
;;                        :kinetic-damping kinetic-damping
;;                        :damping-factor (if adaptive-damping damping nil)
;;                        :dt-scale dt-scale
;;                        :substeps substeps
;;                        :conv-steps 1000
;;                        :post-iter-step
;;                        (lambda (i energy oobf)
;;                          (setf conv-steps (* substeps i))
;;                          (vgplot:title (format nil "Step ~D - substep ~D - KE ~,3E - OOBF ~,3E - damp ~,3E"  step (* i substeps) energy oobf
;;                                                (/ (cl-mpm::sim-damping-factor sim)
;;                                                   (cl-mpm/setup:estimate-critical-damping sim))))
;;                          (format t "Substep ~D~%" i)
;;                          (let ((i (+ 0 i)))
;;                            (when save-vtk-dr
;;                              (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d_~5,'0d.vtk" step i)) sim)
;;                              (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d_~5,'0d.vtk" step i)) sim)
;;                              (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d_~5,'0d.vtk" step i)) sim)))
;;                          (incf total-iter substeps)
;;                          (save-conv-step sim output-dir total-iter step (/ (- (get-internal-real-time) t0) tunits) oobf energy)
;;                          (setf t0 (get-internal-real-time))
;;                          (funcall plotter sim)
;;                          )))
;;                      (vgplot:title (format nil "Step ~D - ~D" step conv-steps))
;;                      ))
;;                  (funcall post-conv-step sim)
;;                  (when save-vtk-loadstep
;;                    (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) sim)
;;                    (cl-mpm/output::save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d.vtk" step)) sim))
;;                  (cl-mpm::finalise-loadstep sim)
;;                  (funcall plotter sim)
;;                  (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step))
;;                                     :terminal "png size 1920,1080"
;;                                     )
;;                  (when save-vtk-loadstep
;;                    (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" step)) sim)
;;                    (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim_p_~5,'0d.vtk" step)) sim ))
;;                  (sleep 0.1d0)
;;                  (swank.live:update-swank)
;;               )))))

(defun run-quasi-time (sim
                        &key (output-dir "./output/")
                          (post-conv-step (lambda (sim)))
                          (plotter (lambda (sim)))
                          (dt-scale 0.5d0)
                          (dt 1d0)
                          (steps 1d0)
                          (substeps 50)
                          (enable-damage t)
                          (enable-plastic t)
                          (max-adaptive-steps 5)
                          (min-adaptive-steps -1)
                          (save-vtk-dr t)
                          (save-vtk-loadstep t)
                          (conv-criteria 1d-3))
  (let ()
    (uiop:ensure-all-directories-exist (list output-dir))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
    (cl-mpm/output::save-simulation-parameters (merge-pathnames output-dir "settings.json")
                                               sim
                                               (list :dt dt
                                                     :criteria-energy conv-criteria
                                                     :criteria-oobf conv-criteria
                                                     :criteria-hist 1d0
                                                     ))
    (save-timestep-preamble output-dir)
    (save-conv-preamble output-dir)
    (cl-mpm:iterate-over-mps
     (cl-mpm:sim-mps sim)
     (lambda (mp)
       (when (typep mp 'cl-mpm/particle::particle-damage)
         (setf (cl-mpm/particle::mp-enable-damage mp) nil))
       (when (typep mp 'cl-mpm/particle::particle-plastic)
         (setf (cl-mpm/particle::mp-enable-plasticity mp) nil))))

    (setf (cl-mpm::sim-dt-scale sim) dt-scale)
    (setf (cl-mpm:sim-mass-scale sim) 1d0)

    (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm/setup:estimate-elastic-dt sim)))
    (setf (cl-mpm:sim-damping-factor sim) 0d0)

    (setf (cl-mpm:sim-enable-damage sim) nil)
    (defparameter *total-iter* 0)
    (let (;(substeps 50)
          (vel-algo (cl-mpm::sim-velocity-algorithm sim))
          (sim-type (class-of sim)))
      (change-class sim 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul)
      (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim) 0d0)
      (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
      ;; find initial quasi-static formation
      (cl-mpm/dynamic-relaxation:converge-quasi-static
       sim
       :energy-crit conv-criteria
       :oobf-crit conv-criteria
       :kinetic-damping nil
       :dt-scale dt-scale
       :conv-steps 10000
       :substeps substeps
       :damping-factor 1d0
       :post-iter-step
       (lambda (i e o)
         (save-conv-step sim output-dir *total-iter* 0 0d0 o e)
         (incf *total-iter* substeps)
         (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_conv_~5,'0d.vtk" i)) sim)
         (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_conv_nodes__~5,'0d.vtk" i)) sim)
         (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_conv_cells__~5,'0d.vtk" i)) sim)
         ))
      (cl-mpm::finalise-loadstep sim)
      (setf (cl-mpm::sim-time sim) 0d0)
      (cl-mpm/dynamic-relaxation::reset-mp-velocity sim)
      (setf (cl-mpm::sim-velocity-algorithm sim) vel-algo)

      (change-class sim sim-type)
      ;; (setf (cl-mpm:sim-enable-damage sim) t)
      (cl-mpm:iterate-over-mps
       (cl-mpm:sim-mps sim)
       (lambda (mp)
         (when (typep mp 'cl-mpm/particle::particle-damage)
           (setf (cl-mpm/particle::mp-enable-damage mp) enable-damage))
         (when (typep mp 'cl-mpm/particle::particle-plastic)
           (setf (cl-mpm/particle::mp-enable-plasticity mp) enable-plastic))))
      (funcall post-conv-step sim)

      (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim) dt)
      (let* ((current-adaptivity 0))
        (loop for step from 1 to steps
              while (cl-mpm::sim-run-sim sim)
              do
                 (let ((quasi-conv nil))
                   (format t "quasi-timestep ~d, dt refine ~d~%" step current-adaptivity)
                   (defparameter *trial-iter* 0)
                   (let ((stagger-iters 0))
                     (loop for i from 0 to max-adaptive-steps
                           while (not quasi-conv)
                           do (progn
                                (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim) (/ dt (expt 2 current-adaptivity)))
                                (format t "trial step ~d, dt refine ~d~%" i current-adaptivity)
                                (setf *trial-iter* i)
                                (multiple-value-bind (conv inc-steps)
                                    (step-quasi-time sim step
                                                     :total-steps *total-iter*
                                                     :plotter plotter
                                                     :output-dir output-dir
                                                     :dt-scale dt-scale
                                                     :substeps substeps
                                                     :enable-damage enable-damage
                                                     :conv-criteria conv-criteria
                                                     :conv-criteria-damage conv-criteria
                                                     :enable-plastic enable-plastic)
                                  (setf quasi-conv conv
                                        stagger-iters inc-steps)
                                  (unless quasi-conv
                                    (if (= current-adaptivity max-adaptive-steps)
                                        (progn
                                          (format t "quasi-time terminated as too many dt refinemets are required")
                                          (loop-finish))
                                        (incf current-adaptivity)))))
                           finally (progn
                                     (format t "finished with ~d dt adaptions - stagger iters ~d" i stagger-iters)
                                     (when (and (= i 1) (< stagger-iters 2))
                                       (setf current-adaptivity
                                             (max min-adaptive-steps
                                                  (- current-adaptivity 1)))))))
                   (unless quasi-conv
                     (setf (cl-mpm::sim-run-sim sim) nil))
                   (funcall plotter sim)
                   (when save-vtk-loadstep
                     (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" step)) sim)
                     (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) sim)
                     (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d.vtk" step)) sim))
                   (swank.live:update-swank)))))
    ))



(defun elastic-static-solution (sim &key (crit 1d-3)
                                      (substeps 10))
  "Take a sim, switch it to implicit quasi-static, disable elastic and plastic and converge switch the class back to the original then enable damage and plastic"
  (let ((vel-algo (cl-mpm::sim-velocity-algorithm sim))
        (sim-type (class-of sim)))
    (change-class sim 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul)
    (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim) 0d0)
    (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
    (set-mp-plastic-damage sim :enable-damage nil :enable-plastic nil)
    ;; Find initial quasi-static formation
    (cl-mpm/dynamic-relaxation:converge-quasi-static
     sim
     :energy-crit crit
     :oobf-crit crit
     :kinetic-damping nil
     :dt-scale 1d0
     :conv-steps 100
     :substeps substeps
     :damping-factor 1d0
     :post-iter-step
     (lambda (i e o)
       (format t "Static step ~D - ~E~%" i o)))
    (cl-mpm::finalise-loadstep sim)
    (set-mp-plastic-damage sim :enable-damage t :enable-plastic t)
    (setf (cl-mpm::sim-time sim) 0d0)
    (cl-mpm/dynamic-relaxation::reset-mp-velocity sim)
    (setf (cl-mpm::sim-velocity-algorithm sim) vel-algo)
    (change-class sim sim-type)))
