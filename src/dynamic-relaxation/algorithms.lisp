(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)
(defparameter *total-iter* 0)


(defgeneric convergence-criteria (sim))

(defmethod convergence-criteria ((sim cl-mpm::mpm-sim))
  t)


(defgeneric convergence-check (sim)
  (:documentation "A check to see if our converged configuration has failed some extra criteria"))

(defmethod convergence-check ((sim cl-mpm::mpm-sim)))




(deftype plastic-damage-type () `(and
                                  cl-mpm/particle::particle-damage
                                  cl-mpm/particle::particle-plastic))

(defun set-mp-plastic-damage (sim &key (enable-damage t) (enable-plastic t))
  (cl-mpm:iterate-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (typecase mp
       (plastic-damage-type
        (setf (cl-mpm/particle::mp-enable-damage mp) enable-damage
              (cl-mpm/particle::mp-enable-plasticity mp) enable-plastic))
       (cl-mpm/particle::particle-damage
        (setf (cl-mpm/particle::mp-enable-damage mp) enable-damage))
       (cl-mpm/particle::particle-plastic
        (setf (cl-mpm/particle::mp-enable-plasticity mp) enable-plastic)))
     ;; (when (typep mp 'cl-mpm/particle::particle-damage)
     ;;   (setf (cl-mpm/particle::mp-enable-damage mp) enable-damage))
     ;; (when (typep mp 'cl-mpm/particle::particle-plastic)
     ;;   (setf (cl-mpm/particle::mp-enable-plasticity mp) enable-plastic))
     )))

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


(defgeneric save-vtks (sim output-dir step))
(defmethod save-vtks (sim output-dir step)
  (format t "Save vtks ~D~%" step)
  (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" step)) sim)
  (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) sim)
  ;; (cl-mpm/output::save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d.vtk" step)) sim)
  (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim_p_~5,'0d.vtk" step)) sim ))

(defgeneric save-vtks-dr-step (sim output-dir step iter))
(defmethod save-vtks-dr-step (sim output-dir step iter)
  (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_step_~5,'0d_~5,'0d.vtk" step iter)) sim)
  (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_step_nodes_~5,'0d_~5,'0d.vtk" step iter)) sim)
  (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_step_cells_~5,'0d_~5,'0d.vtk" step iter)) sim)
  (cl-mpm/penalty:save-vtk-penalties (merge-pathnames output-dir (format nil "sim_step_p_~5,'0d_~5,'0d.vtk" step iter)) sim))

(defun save-conv-preamble (output-dir)
  (with-open-file (stream (merge-pathnames output-dir "conv.csv") :direction :output :if-exists :supersede)
    (format stream "iter,step,real-time,plastic,damage,oobf,energy~%")))

(defun save-conv-step (sim output-dir total-iter step real-time oobf energy)
  (unless (uiop:file-exists-p (merge-pathnames output-dir "conv.csv"))
    (save-conv-preamble output-dir))
  (with-open-file (stream (merge-pathnames output-dir "conv.csv") :direction :output :if-exists :append)
    (format stream "~D,~D,~f,~f,~f,~f,~f~%" total-iter step real-time (get-plastic sim) (get-damage sim)
            (if (sb-ext:float-infinity-p oobf) 0d0 oobf) energy)))

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





(defun generalised-staggered-solve (sim &key
                                          (output-dir "./output/")
                                          (crit 1d-3)
                                          (dt-scale 1d0)
                                          (substeps 10)
                                          (post-iter-step (lambda (i e o)))
                                          (enable-damage t)
                                          (enable-plastic t)
                                          (max-damage-inc 0.6d0)
                                          (damping 1d0)
                                          (staggered-steps 100)
                                          (true-stagger t)
                                          )
  (let* ((damage-prev (get-damage sim))
         (damage damage-prev)
         (oobf-crit   crit)
         (energy-crit 1d0)
         (damage-crit crit)
         (dconv damage-crit)
         (total-i 0)
         (stagger-iters 0))
    (set-mp-plastic-damage sim :enable-damage enable-damage :enable-plastic enable-plastic)

    (setf (cl-mpm:sim-enable-damage sim) nil)
    (loop for stagger-i from 0 to staggered-steps
          while (or (>= dconv damage-crit))
          do
             (progn
               (let ((iv 0))
                 (when true-stagger
                   (setf (cl-mpm:sim-enable-damage sim) nil))
                 (cl-mpm/dynamic-relaxation:converge-quasi-static
                  sim
                  :kinetic-damping nil
                  :oobf-crit oobf-crit
                  :energy-crit energy-crit
                  :dt-scale dt-scale
                  :substeps substeps
                  :convergance-criteria (lambda (sim f o)
                                          (and
                                           (<= o oobf-crit)
                                           (convergence-criteria sim)))
                  :conv-steps 1000
                  :damping-factor damping
                  :post-iter-step
                  (lambda (i e o)
                    (incf iv)
                    (funcall post-iter-step i e o)
                    (unless true-stagger
                      (let ((damage-inc (damage-increment-criteria sim)))
                        (when (> damage-inc max-damage-inc)
                          (format t "Damage criteria failed~%")
                          (error (make-instance 'non-convergence-error
                                                :text "Damage criteria exeeded"
                                                :ke-norm 0d0
                                                :oobf-norm 0d0)))))
                    ))
                 (if (and
                      enable-damage
                      (typep sim 'cl-mpm/damage::mpm-sim-damage))
                     (progn
                       (let ((fast-trial-conv oobf-crit)
                             (damage-iter t)
                             (save-update nil))
                         (loop for d from 0 to 1000
                               while (and (<= fast-trial-conv oobf-crit)
                                          damage-iter)
                               do
                                  (setf (cl-mpm:sim-enable-damage sim) t)
                                  (cl-mpm/damage::calculate-damage sim (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim))
                                  (setf (cl-mpm:sim-enable-damage sim) nil)

                                  (setf damage (get-damage sim))
                                  (setf dconv
                                        (compute-damage-delta sim)
                                        ;; (if (> damage 0d0)
                                        ;;     (if (> damage-prev 0d0)
                                        ;;         (/ (compute-damage-delta sim) damage-prev)
                                        ;;         sb-ext:double-float-positive-infinity)
                                        ;;     0d0)
                                        )
                                  (setf damage-prev damage)
                                  (unless (>= dconv damage-crit)
                                    (setf damage-iter nil))
                                  (format t "Damage ~E - prev damage ~E ~%" damage damage-prev)
                                        (format t "step ~D/~D - d-conv ~E~%" stagger-i d dconv)
                                        (let ((damage-inc (damage-increment-criteria sim)))
                                          (when (> damage-inc max-damage-inc)
                                            (format t "Damage criteria failed~%")
                                            (error (make-instance 'non-convergence-error
                                                                  :text "Damage criteria exeeded"
                                                                  :ke-norm 0d0
                                                                  :oobf-norm 0d0))))
                                        (format t "Def crit ~E~%" (compute-max-deformation sim))
                                        (convergence-check sim)

                                        (setf damage-prev damage)
                                        (when damage-iter
                                          (dotimes (i 2)
                                            (cl-mpm:update-sim sim))
                                          (setf fast-trial-conv (cl-mpm::sim-stats-oobf sim))
                                          (format t "Fast trial oobf ~E~%" fast-trial-conv)
                                          (setf save-update t)))
                               (when save-update
                                 (incf iv)
                                 (funcall post-iter-step iv (cl-mpm::sim-stats-energy sim) (cl-mpm::sim-stats-oobf sim))))
                             ;; (setf (cl-mpm:sim-enable-damage sim) t)
                             ;; (when (typep sim 'cl-mpm/damage::mpm-sim-damage)
                             ;;   (cl-mpm/damage::calculate-damage sim (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim)))
                             ;; (setf damage (get-damage sim))
                             ;; (setf dconv (if (> damage 0d0)
                             ;;                 (if (> damage-prev 0d0)
                             ;;                     (/ (- damage damage-prev) damage-prev)
                             ;;                     sb-ext:double-float-positive-infinity)
                             ;;                 0d0))
                             ;; (setf damage-prev damage)
                             )
                           (setf dconv 0d0)
                           ))))
    (when (>= dconv damage-crit)
      (format t "Failed to converge damage during stagger~%")
      (error (make-instance 'non-convergence-error
                            :text "Failed to converge damage during stagger"
                            :ke-norm 0d0
                            :oobf-norm 0d0)))))

(declaim (notinline step-quasi-time))
(defun step-quasi-time (sim
                        global-step
                        &key (dt-scale 1d0)
                          (substeps 50)
                          (total-steps 0)
                          (damping 1d0)
                          (sub-conv-steps 200)
                          (conv-criteria 1d-3)
                          (conv-criteria-damage 1d-3)
                          (output-dir "./output/")
                          (save-vtk-dr t)
                          (enable-damage t)
                          (enable-plastic t)
                          (max-damage-inc 0.6d0)
                          (true-stagger t)
                          (plotter (lambda (sim))))
  (let ((total-i 0))
    (handler-case
        (progn
          (let* ((damage-prev (get-damage sim))
                 (damage damage-prev)
                 (oobf-crit   conv-criteria)
                 (energy-crit
                   1d0
                   ;conv-criteria
                              )
                 (damage-crit conv-criteria-damage)
                 (dconv damage-crit)
                 (inertia 0d0)
                 (intertia-crit 1d-3)
                 (stagger-iters 0)
                 )
            (reset-mp-velocity sim)
            (set-mp-plastic-damage sim :enable-damage enable-damage :enable-plastic enable-plastic)
            (setf (cl-mpm:sim-enable-damage sim) nil)
            (loop for stagger-i from 0 to 100
                  while (or (>= dconv damage-crit)
                            (>= (cl-mpm::sim-stats-oobf sim) oobf-crit))
                  do
                     (progn
                       (cl-mpm/dynamic-relaxation:converge-quasi-static
                        sim
                        :oobf-crit oobf-crit
                        :energy-crit energy-crit
                        :dt-scale dt-scale
                        :substeps substeps
                        :conv-steps sub-conv-steps
                        :damping-factor damping
                        :post-iter-step
                        (lambda (i e o)
                          (funcall plotter sim)
                                        ;(vgplot:title (format nil "substep ~D - oobf: ~E" i o))
                          (incf total-i)
                          (when save-vtk-dr
                            (when (= (mod i 1) 0)
                              (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_step_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
                              (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_step_nodes_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
                              (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_step_cells_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
                              ))
                          (format t "Def crit ~E~%" (compute-max-deformation sim))
                          (let ((max-def (compute-max-deformation sim)))
                            (when (> max-def 2d0)
                              (format t "Deformation gradient criteria exceeded~%")
                              (error (make-instance 'non-convergence-error
                                                    :text "Deformation gradient J exceeded"
                                                    :ke-norm 0d0
                                                    :oobf-norm 0d0)))
                            (let ((true-intertia (true-intertial-criteria sim (sim-dt-loadstep sim))))
                              (format t "True intertia ~E~%" true-intertia)
                              (save-conv-step sim output-dir *total-iter* global-step
                                              0d0
                                              o
                                              ;; true-intertia
                                              max-def
                                              )
                              (setf inertia true-intertia)
                              (when (> true-intertia 1d-4)
                                (format t "Inertia criteria exceeded~%")
                                (error (make-instance 'error-inertia-criteria
                                                      :text "True inertia exceeded"
                                                      :inertia-norm true-intertia)))))
                          (unless true-stagger
                            (let ((damage-inc (damage-increment-criteria sim)))
                              (format t "Damage ~E - prev damage ~E - inc ~E~%" damage damage-prev damage-inc)
                              (when (> damage-inc max-damage-inc)
                                (format t "Damage criteria failed ~E~%" damage-inc)
                                (error (make-instance 'error-damage-criteria
                                                      :text "Damage criteria exeeded"
                                                      :max-damage-inc 0d0)))))
                          (convergence-check sim)
                          ;; (save-conv-step sim output-dir *total-iter* global-step 0d0 o 0d0)
                          (incf *total-iter* substeps)))

                       (let ((fast-trial-conv oobf-crit)
                             (damage-iter t))
                         (unless (typep sim 'cl-mpm/damage::mpm-sim-damage)
                           (setf dconv 0d0))
                         (when (typep sim 'cl-mpm/damage::mpm-sim-damage)
                           (loop for d from 0 to 100
                                 while (and (<= fast-trial-conv (* 1d1 oobf-crit))
                                            damage-iter)
                                 do
                                    (when (typep sim 'cl-mpm/damage::mpm-sim-damage)
                                      (setf (cl-mpm:sim-enable-damage sim) t)
                                      (cl-mpm/damage::calculate-damage sim (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim))
                                      (when true-stagger
                                        (setf (cl-mpm:sim-enable-damage sim) nil)))
                                    (setf damage (get-damage sim))
                                    (setf dconv
                                          (compute-damage-delta sim)
                                          ;; (if (> damage 0d0)
                                          ;;     (if (> damage-prev 0d0)
                                          ;;         (/ (compute-damage-delta sim) damage-prev)
                                          ;;         ;; (abs (/ (- damage damage-prev) damage-prev))
                                          ;;         sb-ext:double-float-positive-infinity)
                                          ;;     0d0)
                                          )
                                    (setf damage-prev damage)
                                    (when (< dconv damage-crit)
                                      (setf damage-iter nil))
                                    (format t "step ~D/~D - d-conv ~E - ~E~%" stagger-i d dconv fast-trial-conv)
                                    (let ((damage-inc (damage-increment-criteria sim)
                                            ;; (compute-max-damage-energy-crit sim)
                                            ))
                                      (format t "Damage ~E - prev damage ~E - inc ~E~%" damage damage-prev damage-inc)
                                      ;; (format t "Damage energy criterion ~E~%" (compute-max-damage-energy-crit sim))
                                      (when (> damage-inc max-damage-inc)
                                        (format t "Damage criteria failed ~E~%" damage-inc)
                                        (when save-vtk-dr
                                          (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_step_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
                                          (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_step_nodes_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
                                          (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_step_cells_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
                                          )
                                        
                                        (error (make-instance 'error-damage-criteria
                                                              :text "Damage criteria exeeded"
                                                              :max-damage-inc 0d0))))
                                    (when t;damage-iter
                                      (dotimes (i 1)
                                        (cl-mpm:update-sim sim))
                                      (setf fast-trial-conv (cl-mpm::sim-stats-oobf sim))
                                      (format t "fast trial ~E~%" fast-trial-conv))
                                    ;; (incf *total-iter* 1)
                                    (incf total-i)
                                    (when t;(= (mod d 5) 0)
                                      (when save-vtk-dr
                                        (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_step_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
                                        (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_step_nodes_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
                                        (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_step_cells_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
                                        )
                                      )
                                 )))
                       (setf damage-prev damage)
                       (when save-vtk-dr
                         (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_step_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
                         (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_step_nodes_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
                         (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_step_cells_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
                         )
                       (incf stagger-iters)
                       ;; (when (> inertia intertia-crit)
                       ;;   (error (make-instance 'non-convergence-error
                       ;;                         :text "Quasi-time inertia was too large"
                       ;;                         :ke-norm inertia
                       ;;                         :oobf-norm oobf-crit)))
                       ))
            (when (or (> (cl-mpm::sim-stats-oobf sim) oobf-crit)
                      (> dconv damage-crit))
              (format t "Staggered solve didn't converge ~E ~E~%" dconv (cl-mpm::sim-stats-oobf sim))
              (error (make-instance 'non-convergence-error
                                    :text "Staggered solve didn't converge ~E ~E"
                                    :ke-norm  dconv
                                    :oobf-norm  (cl-mpm::sim-stats-oobf sim))))

            (when save-vtk-dr
                (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_step_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
              (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_step_nodes_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim))
            (cl-mpm::finalise-loadstep sim)
            (save-timestep sim output-dir global-step :QUASI-STATIC)
            (values t stagger-iters)))
      (cl-mpm/errors:error-simulation (c)
        (princ c)
        (cl-mpm::reset-loadstep sim)
        (values nil 0))))
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
    (format t "E crit ~E - OOBF crit ~E~%" e-crit oobf-crit)
    (time (loop for step from 0 to max-steps
                while (and (cl-mpm::sim-run-sim sim)
                           (or (>= energy e-crit)
                               (>= oobf oobf-crit)
                               ;; (< step 4)
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
                        (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm::calculate-min-dt sim)))
                        (incf oobf (estimate-static-oobf sim))
                        (incf energy (cl-mpm::sim-stats-energy sim))
                        (incf work (cl-mpm::sim-stats-power sim))))

                     (setf
                      energy (/ energy substeps)
                      oobf (/ oobf substeps))
                     (if (= work 0d0)
                         (setf energy 0d0)
                         (setf energy (abs (/ energy work))))

                     ;; (setf (cl-mpm/aggregate::sim-enable-aggregate sim) (and (< energy 1d-2) (< oobf 1d-2)))

                     ;; (when enable-mass-scaling
                     ;;   (let ((res (and (< energy 1d-2)
                     ;;                   (< oobf 1d-2))))
                     ;;     (setf (cl-mpm:sim-mass-scale sim)
                     ;;           (if res
                     ;;               1d2
                     ;;               1d0))
                     ;;     (format t "Accelerate real time ~A~%" res)))
                     ;; (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm/setup:estimate-elastic-dt sim)))
                     ;; (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time sim target-time :dt-scale dt-scale)
                     ;;   (format t "CFL dt estimate: ~f~%" dt-e)
                     ;;   (format t "CFL step count estimate: ~D~%" substeps-e)
                     ;;   (setf (cl-mpm:sim-dt sim) dt-e)
                     ;;   (setf substeps substeps-e)u
                     (when (or (>= energy e-crit)
                               (>= oobf oobf-crit))
                       (format t "Inertia passed~%")
                       (setf intertial-passed t))
                     (format t "Residuals ~E ~E ~%" energy oobf)
                     (setf (cl-mpm::sim-stats-oobf sim) oobf)
                     (save-timestep sim output-dir global-step :DYNAMIC)
                     (funcall plotter sim)
                     (swank.live:update-swank))))))

(declaim (notinline run-multi-stage))
(defun run-multi-stage (sim
                        &key (output-dir "./output/")
                          (conv-load-steps 1)
                          (post-conv-step (lambda (sim)))
                          (plotter (lambda (sim)))
                          (dt-scale 1d0)
                          (damping-factor 1d0)
                          (explicit-dt-scale 0.9d0)
                          (explicit-damping-factor 1d-3)
                          (elastic-dt-margin 1000)
                          (dt 1d0)
                          (steps nil)
                          (total-time 1d0)
                          (max-adaptive-steps 5)
                          (min-adaptive-steps -1)
                          (conv-criteria 1d-3)
                          (substeps 50)
                          (enable-damage t)
                          (enable-plastic t)
                          (save-vtk-dr t)
                          (save-vtk-loadstep t)
                          (max-damage-inc 0.6d0)
                          (setup-quasi-static (lambda (sim)))
                          (setup-dynamic (lambda (sim)))
                          (explicit-dynamic-solver 'cl-mpm::mpm-sim-usf))
  (let ()
    (when steps
      (setf total-time (* dt steps)))
    (uiop:ensure-all-directories-exist (list output-dir))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
    (cl-mpm/output::save-simulation-parameters (merge-pathnames output-dir "settings.json")
                                               sim
                                               (list :dt dt
                                                     :criteria-energy conv-criteria
                                                     :criteria-oobf conv-criteria
                                                     :criteria-hist 1d0))
    (save-timestep-preamble output-dir)
    (save-conv-preamble output-dir)
    (setf (cl-mpm::sim-dt-scale sim) dt-scale)
    (setf (cl-mpm:sim-mass-scale sim) 1d0)
    (setf (cl-mpm:sim-damping-factor sim) 0d0)
    (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm/setup:estimate-elastic-dt sim)))
    (defparameter *total-iter* 0)
    (let ((quasi-static-solver (class-of sim))
          (vel-algo :BLEND)
          (load-steps conv-load-steps)
          (i 0)
          (grav (cl-mpm:sim-gravity sim)))

      (reset-mp-velocity sim)
      (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim) 0d0)
      (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
      (set-mp-plastic-damage sim :enable-plastic nil :enable-damage nil)
      (funcall setup-quasi-static sim)
      (loop for lstp from 1 to load-steps
            do
               (progn
                 (cl-mpm/dynamic-relaxation:converge-quasi-static
                  sim
                  :energy-crit 1d0;conv-criteria
                  :oobf-crit conv-criteria
                  :kinetic-damping nil
                  :dt-scale dt-scale
                  :conv-steps 1000
                  :substeps substeps
                  :damping-factor damping-factor
                  :post-iter-step
                  (lambda (j e o)
                    (incf i)
                    (save-conv-step sim output-dir *total-iter* 0 0d0 o e)
                    (incf *total-iter* substeps)
                    (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_conv_~5,'0d.vtk" i)) sim)
                    (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_conv_nodes__~5,'0d.vtk" i)) sim)
                    (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_conv_cells__~5,'0d.vtk" i)) sim)
                    ;; (cl-mpm/penalty::save-vtk-penalties (merge-pathnames output-dir (format nil "sim_conv_p_~5,'0d.vtk" i)) sim)
                    ))
                 (cl-mpm::finalise-loadstep sim)
                 ))
      (setf (cl-mpm::sim-time sim) 0d0)
      (cl-mpm/dynamic-relaxation::reset-mp-velocity sim)
      (reset-mp-velocity sim)
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
      (funcall setup-quasi-static sim)
      (let* ((current-adaptivity 0)
             (prev-steps-easy (list t t))
             (prev-step-iter 0)
             (elastic-dt (cl-mpm/setup::estimate-elastic-dt sim))
             (max-steps (floor total-time (/ dt (expt 2 max-adaptive-steps))))
             )
        (format t "Elastic dt ~E, override quasi-static at ~E~%" elastic-dt (* elastic-dt elastic-dt-margin))
        (loop for step from 1 to (+ 1 max-steps)
              while (and (cl-mpm::sim-run-sim sim)
                         (< (cl-mpm::sim-time sim) total-time))
              do
                 (let ((quasi-conv nil))
                   (format t "Quasi-timestep ~D, dt refine ~D - dt ~E~%" step current-adaptivity (/ dt (expt 2 current-adaptivity)))
                   (defparameter *trial-iter* 0)
                   (let ((stagger-iters 0))
                     (loop for i from 0 to max-adaptive-steps
                           while (not quasi-conv)
                           do (progn
                                (let* ((adapted-dt (/ dt (expt 2 current-adaptivity))))
                                  (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim)
                                        (+ (min (- total-time (cl-mpm::sim-time sim)) adapted-dt) 1d-15))
                                  (format t "Trial step ~D, dt refine ~D - dt ~E ~%" i current-adaptivity (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim)))
                                (setf *trial-iter* i)
                                (multiple-value-bind (conv inc-steps)
                                    (step-quasi-time sim step
                                                     :total-steps *total-iter*
                                                     ;; :plotter plotter
                                                     :output-dir output-dir
                                                     :dt-scale dt-scale
                                                     :substeps substeps
                                                     :enable-damage enable-damage
                                                     :enable-plastic enable-plastic
                                                     :conv-criteria conv-criteria
                                                     :conv-criteria-damage conv-criteria
                                                     :max-damage-inc max-damage-inc
                                                     :save-vtk-dr save-vtk-dr
                                                     )
                                  (setf quasi-conv conv
                                        stagger-iters inc-steps)
                                  (format t "Quasi-conv? ~A~%" quasi-conv)
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
                                        (incf current-adaptivity))
                                    (when (> i 4)
                                      (format t "We've adaptive multiple sucessive times, double jump dt adaption~%")
                                      (incf i 1)))))
                           finally (progn
                                     (format t "Finished with ~D dt adaptions - stagger iters ~D - conv ~A~%" (- i 1) stagger-iters quasi-conv)
                                     (incf prev-step-iter)
                                     (if (and (= i 1) (<= stagger-iters 3))
                                         (progn
                                           (setf (nth (mod prev-step-iter (length prev-steps-easy)) prev-steps-easy) t)
                                           (format t "Potential adaption easy steps ~A~%" prev-steps-easy)
                                           (when (every #'identity prev-steps-easy)
                                             (setf current-adaptivity
                                                   (max min-adaptive-steps
                                                        (- current-adaptivity 1)))))
                                         (setf (nth (mod prev-step-iter 2) prev-steps-easy) nil)
                                         ))))
                   (unless quasi-conv
                     ;;We've adapted down to a min
                     (format t "Start real-timestepping~%")
                     (let ((dt-loadstep (* 1d0 (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim))))
                       (cl-mpm/dynamic-relaxation::reset-mp-velocity sim)
                       (change-class sim explicit-dynamic-solver)
                       (setf (cl-mpm::sim-velocity-algorithm sim) :FLIP)
                       (funcall setup-dynamic sim)
                       (step-real-time sim step
                                       :output-dir output-dir
                                       :plotter plotter
                                       :dt-scale explicit-dt-scale
                                       :criteria (* conv-criteria 1d0)
                                       :damping explicit-damping-factor
                                       :target-time (* 0.1d0 dt-loadstep)
                                       :enable-damage enable-damage
                                       :enable-plastic enable-plastic)
                       (change-class sim quasi-static-solver)
                       (format t "Finished real-timestepping~%")
                       (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
                       (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim) dt-loadstep)
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






;; (defun converge-staggered (sim
;;                            criteria
;;                            damage-criteria
;;                            kinetic-damping
;;                            adaptive-damping
;;                            dt-scale
;;                            substeps
;;                            damping
;;                            enable-damage
;;                            enable-plasticity
;;                            step
;;                            output-dir
;;                            &key (save-vtk-dr nil)
;;                              (plotter (lambda (sim)))
;;                              )
;;   (let* ((damage-prev (get-damage sim))
;;          (damage-eval t)
;;          (damage damage-prev)
;;          (max-iter 1000)
;;          (dconv damage-criteria))
;;     (setf (cl-mpm::sim-enable-damage sim) nil)
;;     (cl-mpm/dynamic-relaxation:converge-quasi-static
;;      sim
;;      :oobf-crit criteria
;;      :energy-crit criteria
;;      :kinetic-damping kinetic-damping
;;      :damping-factor (if adaptive-damping damping nil)
;;      :dt-scale dt-scale
;;      :substeps substeps
;;      :conv-steps 1000
;;      :post-iter-step
;;      (lambda (i energy oobf)
;;        (incf *total-iter* substeps)
;;        (save-conv-step sim output-dir *total-iter* step 0d0 oobf energy)))
;;     (loop for i from 0 to max-iter
;;           while (>= dconv damage-criteria)
;;           do (progn
;;                (when (typep sim 'cl-mpm/damage::mpm-sim-damage)
;;                  (setf (cl-mpm::sim-enable-damage sim) damage-eval)
;;                  (cl-mpm/damage::calculate-damage sim (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim))
;;                  (setf damage (get-damage sim)))
;;                (if (= damage damage-prev)
;;                    (setf dconv 0d0)
;;                    (when (> damage damage-prev)
;;                      (setf dconv (if (> damage-prev 0d0) (/ (- damage damage-prev) damage-prev) sb-ext:double-float-positive-infinity))))
;;                (setf damage-prev damage)
;;                (setf (cl-mpm::sim-enable-damage sim) nil)
;;                (cl-mpm/dynamic-relaxation:converge-quasi-static
;;                 sim
;;                 :oobf-crit criteria
;;                 :energy-crit criteria
;;                 :kinetic-damping kinetic-damping
;;                 :damping-factor (if adaptive-damping damping nil)
;;                 :dt-scale dt-scale
;;                 :substeps substeps
;;                 :conv-steps 1000
;;                 :post-iter-step
;;                 (lambda (i energy oobf)
;;                   (incf *total-iter* substeps)
;;                   (save-conv-step sim output-dir *total-iter* step 0d0 oobf energy))
;;                 )
;;                (funcall plotter sim)
;;                (when save-vtk-dr
;;                  (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d_~5,'0d.vtk" step i)) sim)
;;                  (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d_~5,'0d.vtk" step i)) sim)
;;                  (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d_~5,'0d.vtk" step i)) sim))

;;                ))))


(defun run-load-control (sim
                         &key (output-dir "./output/")
                           (loading-function nil)
                           (load-steps 10)
                           (initial-load 0d0)
                           (substeps 50)
                           (damping 1d0)
                           (kinetic-damping nil)
                           (criteria 1d-3)
                           (post-conv-step (lambda (sim)))
                           (post-iter-step (lambda (sub-iter oobf energy)))
                           (pre-step (lambda ()))
                           (plotter (lambda (sim)))
                           (enable-damage t)
                           (enable-plastic t)
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
            (tunits internal-time-units-per-second))
        (loop for step from 1 to load-steps
              while (cl-mpm::sim-run-sim sim)
              do
                 (progn
                   (format t "Load step ~D~%" step)
                   (defparameter *ke-last* 0d0)
                   ;; (pprint step)
                   (let* ((percent (/ (float step) load-steps)))
                     ;; (declare (double-float percent initial-load))
                     (funcall loading-function (+ initial-load (* (- 1d0 initial-load) percent))))
                   (let ((conv-steps 0)
                         (i 0))
                     ;; (when damping
                     ;;   (setf (cl-mpm::sim-damping-factor sim)
                     ;;         (* damping (cl-mpm/setup:estimate-critical-damping sim))))
                     (setf (cl-mpm/dynamic-relaxation::sim-damping-scale sim) damping)
                     (time
                      (;;cl-mpm/dynamic-relaxation:converge-quasi-static
                       generalised-staggered-solve
                       sim
                       :crit criteria
                       :enable-damage enable-damage
                       :enable-plastic enable-plastic
                       :dt-scale dt-scale
                       :substeps substeps
                       :damping damping
                       ;; :conv-steps 1000
                       :post-iter-step
                       (lambda (i energy oobf)
                         (funcall post-iter-step i energy oobf)
                         (setf conv-steps (* substeps i))
                         (format t "Substep ~D~%" i)
                         (let ((i (+ 0 i)))
                           (when save-vtk-dr
                             (save-vtks-dr-step sim output-dir step i)))
                         (incf *total-iter* substeps)
                         (save-conv-step sim output-dir *total-iter* step 0d0 oobf energy)
                         (funcall plotter sim))))))
                 ;; (when save-vtk-loadstep
                 ;;   (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) sim)
                 ;;   (cl-mpm/output::save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d.vtk" step)) sim))
                 (cl-mpm::finalise-loadstep sim)
                 (funcall post-conv-step sim)
                 (save-vtks sim output-dir step)
                 (funcall plotter sim)
                 ;; (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step))
                 ;;                    :terminal "png size 1920,1080"
                 ;;                    )
                 ;; (when save-vtk-loadstep
                 ;;   (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" step)) sim)
                 ;;   (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim_p_~5,'0d.vtk" step)) sim ))
                 ;; (sleep 0.1d0)
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
                         (post-load-step (lambda (sim)))
                         (plotter (lambda (sim)))
                         (dt-scale 0.5d0)
                         (dt 1d0)
                         (total-time 1d0)
                         (substeps 50)
                         (enable-damage t)
                         (enable-plastic t)
                         (max-adaptive-steps 5)
                         (min-adaptive-steps -1)
                         (save-vtk-loadstep t)
                         (save-vtk-conv t)
                         (save-vtk-dr t)
                         (damping 1d0)
                         (elastic-solver 'mpm-sim-dr-ul)
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
    (save-timestep-preamble output-dir)
    (save-conv-preamble output-dir)
    (setf (cl-mpm::sim-dt-scale sim) dt-scale)
    (setf (cl-mpm:sim-mass-scale sim) 1d0)
    (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm/setup:estimate-elastic-dt sim)))
    (setf (cl-mpm:sim-damping-factor sim) 0d0)
    (setf (cl-mpm:sim-enable-damage sim) nil)
    (defparameter *total-iter* 0)
    (let (;(substeps 50)
          (vel-algo (cl-mpm::sim-velocity-algorithm sim))
          (sim-type (class-of sim)))
      (when (not (subtypep (type-of sim) elastic-solver))
        (change-class sim elastic-solver))
      (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim) 0d0)
      (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
      (set-mp-plastic-damage sim :enable-plastic nil :enable-damage nil)
      ;; find initial quasi-static formation
      (cl-mpm/dynamic-relaxation:converge-quasi-static
       sim
       :energy-crit conv-criteria
       :oobf-crit conv-criteria
       :kinetic-damping nil
       :dt-scale dt-scale
       :conv-steps 10000
       :substeps substeps
       :damping-factor damping
       :post-iter-step
       (lambda (i e o)
         (save-conv-step sim output-dir *total-iter* 0 0d0 o e)
         (incf *total-iter* substeps)
         (when save-vtk-conv
           (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_conv_~5,'0d.vtk" i)) sim)
           (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_conv_nodes__~5,'0d.vtk" i)) sim)
           (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_conv_cells__~5,'0d.vtk" i)) sim))
         (funcall plotter sim)
         ))
      (cl-mpm::finalise-loadstep sim)
      (cl-mpm::reset-grid (cl-mpm:sim-mesh sim))
      (format t "Reset grid, finalise loadstep~%")
      (setf (cl-mpm::sim-time sim) 0d0)
      (cl-mpm/dynamic-relaxation::reset-mp-velocity sim)
      (setf (cl-mpm::sim-velocity-algorithm sim) vel-algo)
      (format t "Set plastic-damage~%")
      (set-mp-plastic-damage sim :enable-plastic enable-plastic :enable-damage enable-damage)
      (format t "Change class~%")
      (change-class sim sim-type)
      (format t "Call post-conv~%")
      (funcall post-conv-step sim)
      (format t "Start iter~%")
      (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim) dt)
      (let* ((current-adaptivity 0)
             (max-steps (floor total-time (/ dt (expt 2 max-adaptive-steps))))
             (sim-time 0d0))
        (loop for step from 1 to (+ 1 max-steps)
              while (and (cl-mpm::sim-run-sim sim)
                         (< sim-time total-time))
              do
                 (let ((quasi-conv nil))
                   (format t "quasi-timestep ~d - time ~E, dt refine ~d~%" step sim-time current-adaptivity)
                   (defparameter *trial-iter* 0)
                   (let ((stagger-iters 0))
                     (loop for i from 0 to max-adaptive-steps
                           while (not quasi-conv)
                           do (progn
                                (let* ((adapted-dt (/ dt (expt 2 current-adaptivity))))
                                  (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim)
                                        (+ (min (- total-time sim-time) adapted-dt) 1d-15))
                                  (format t "trial step ~d, dt refine ~d - dt ~E~%" i current-adaptivity (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim)))
                                (setf *trial-iter* i)
                                (multiple-value-bind (conv inc-steps)
                                    (step-quasi-time sim step
                                                     :total-steps *total-iter*
                                                     :plotter plotter
                                                     :output-dir output-dir
                                                     :dt-scale dt-scale
                                                     :substeps substeps
                                                     :enable-damage enable-damage
                                                     :damping damping
                                                     :save-vtk-dr save-vtk-dr
                                                     :conv-criteria conv-criteria
                                                     :conv-criteria-damage conv-criteria
                                                     :max-damage-inc 0.5d0
                                                     :enable-plastic enable-plastic)
                                  (setf quasi-conv conv
                                        stagger-iters inc-steps)
                                  (unless quasi-conv
                                    (if (= current-adaptivity max-adaptive-steps)
                                        (progn
                                          (cl-mpm::reset-node-displacement sim)
                                          (format t "quasi-time terminated as too many dt refinemets are required~%")
                                          (loop-finish))
                                        (incf current-adaptivity)))))
                           finally (progn
                                     (format t "finished with ~d dt adaptions - stagger iters ~d~%" (- i 1) stagger-iters)
                                     (when (and (= i 1) (< stagger-iters 4))
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


(defun reset-mp-velocity (sim)
  (cl-mpm:iterate-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (cl-mpm/fastmaths:fast-zero (cl-mpm/particle::mp-velocity mp)))))

(defun reset-nominal-displacement (sim)
  (cl-mpm:iterate-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (cl-mpm/fastmaths:fast-zero (cl-mpm/particle::mp-displacement mp)))))

(defun elastic-static-solution (sim &key (crit 1d-3)
                                      (substeps 10)
                                      (dt-scale 1d0)
                                      (post-iter-step
                                       (lambda (i e o)
                                         (format t "Static step ~D - ~E~%" i o))))
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
     :dt-scale dt-scale
     :conv-steps 100
     :substeps substeps
     :damping-factor 1d0
     :post-iter-step
     post-iter-step
     )
    (cl-mpm::finalise-loadstep sim)
    (set-mp-plastic-damage sim :enable-damage t :enable-plastic t)
    (setf (cl-mpm::sim-time sim) 0d0)
    (cl-mpm/dynamic-relaxation::reset-mp-velocity sim)
    (setf (cl-mpm::sim-velocity-algorithm sim) vel-algo)
    (change-class sim sim-type)))

(defun run-adaptive-load-control (sim
                         &key (output-dir "./output/")
                           (loading-function nil)
                           (load-steps 10)
                           (initial-load 0d0)
                           (substeps 50)
                           (damping 1d0)
                           (kinetic-damping nil)
                           (criteria 1d-3)
                           (post-conv-step (lambda (sim)))
                           (post-iter-step (lambda (sub-iter oobf energy)))
                           (pre-step (lambda ()))
                           (plotter (lambda (sim)))
                           (enable-damage t)
                           (enable-plastic t)
                           (save-vtk-dr t)
                           (save-vtk-loadstep t)
                           (max-adaptive-steps 5)
                           (min-adaptive-steps 0)
                           (max-damage-inc 0.3d0)
                           (dt-scale 1d0))
  (uiop:ensure-all-directories-exist (list output-dir))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (save-conv-preamble output-dir)
  (cl-mpm/output::save-simulation-parameters (merge-pathnames output-dir "settings.json") sim)
  (funcall pre-step)
  ;; (save-conv-step sim output-dir 0 0 0d0 0d0)
  (setf lparallel:*debug-tasks-p* nil)
  (defparameter *total-iter* 0)
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (let* ((load (cl-mpm:sim-gravity sim))
           (finished t))
      (unless loading-function
        (setf loading-function (lambda (percent)
                                 (format t "Loading factor ~E~%" percent)
                                 (setf (cl-mpm:sim-gravity sim)
                                       (* load percent)))))
      (setf (cl-mpm::sim-dt-scale sim) dt-scale)
      (let* ((initial-load-step-size (/ (- 1d0 initial-load) load-steps))
             (load-step-size initial-load-step-size)
             (current-load initial-load)
             (current-adaptivity 0)
             )
        (loop for step from 1; to load-steps
              while (and (cl-mpm::sim-run-sim sim)
                         (<= current-load 1d0))
              do
                 (progn
                   (format t "Load step ~D~%" step)
                   (defparameter *ke-last* 0d0)
                   (let ((conv-steps 0)
                         (i-total 0))
                     (flet ((trial-solve ()
                              (handler-case
                                  (progn
                                    (format t "Load applied ~E~%" (min 1d0 (+ current-load load-step-size)))
                                    (funcall loading-function (min 1d0 (+ current-load load-step-size)))
                                    (time
                                     (;cl-mpm/dynamic-relaxation:converge-quasi-static
                                      generalised-staggered-solve
                                      sim
                                      :crit criteria
                                      :damping damping
                                      :output-dir output-dir
                                      :dt-scale dt-scale
                                      :substeps substeps
                                      :enable-damage enable-damage
                                      :enable-plastic enable-plastic
                                      :max-damage-inc max-damage-inc
                                      :post-iter-step
                                      (lambda (i-g energy oobf)
                                        (incf i-total)
                                        (let ((i i-total))
                                          (funcall post-iter-step i energy oobf)
                                          (setf conv-steps (* substeps i))
                                          (format t "Substep ~D~%" i)
                                          (let ((i (+ 0 i)))
                                            (when save-vtk-dr
                                              (save-vtks-dr-step sim output-dir step i)))

                                          (incf *total-iter* substeps)
                                          (save-conv-step sim output-dir *total-iter* step 0d0 oobf energy)
                                          (funcall plotter sim)))))
                                    (format t "Generalised solve passed~%")
                                    t)
                                (cl-mpm/errors:error-simulation (c)
                                  (princ c)
                                  (cl-mpm::reset-loadstep sim)
                                  (values nil 0))
                                (error (c)
                                  (format t "A non simulation error was thrown!")
                                  (princ c)
                                  (cl-mpm::reset-loadstep sim)
                                  (cl-mpm::sim-run-sim sim) nil
                                  (values nil 0)))))
                       (let ((quasi-conv nil))
                         (loop for i from 0 to max-adaptive-steps
                               while (not quasi-conv)
                               do (progn
                                    (setf load-step-size (* initial-load-step-size (expt 2 (- current-adaptivity))))
                                    (format t "Load step size ~E adapted by ~D~%" load-step-size current-adaptivity)
                                    (setf quasi-conv (trial-solve))
                                    (unless quasi-conv
                                      (format t "Failed conv, adapting~%" )
                                      (incf current-adaptivity))
                                    (if (or (= current-adaptivity max-adaptive-steps))
                                        (progn
                                          (loop-finish))))
                               finally (progn
                                         (format t "Finished with ~D adaptions- conv ~A~%" (- i 1) quasi-conv)
                                         (when (and (= i 1))
                                           (setf current-adaptivity
                                                 (max min-adaptive-steps
                                                      (- current-adaptivity 1))))
                                         ))
                         (unless quasi-conv;(>= current-adaptivity max-adaptive-steps)
                           (format t "Solve failed completly~%" )
                           (setf (cl-mpm::sim-run-sim sim) nil
                                 finished nil))
                         ))))
                 (incf current-load load-step-size)
                 (funcall post-conv-step sim)
                 (cl-mpm::finalise-loadstep sim)
                 (save-vtks sim output-dir step)
                 (funcall plotter sim)
                 (swank.live:update-swank)))
      finished)))

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
                   (elastic-solver 'mpm-sim-dr-ul)
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
    (save-timestep-preamble output-dir)
    (save-conv-preamble output-dir)
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

      (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm/setup:estimate-elastic-dt sim)))
      (setf (cl-mpm:sim-damping-factor sim) 0d0)

      (setf (cl-mpm:sim-enable-damage sim) nil)
      (let (;(substeps 50)
            (vel-algo (cl-mpm::sim-velocity-algorithm sim))
            (sim-type (class-of sim)))
        (change-class sim elastic-solver)
        (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim) 0d0)
        (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
        ;; find initial quasi-static formation
        (cl-mpm/dynamic-relaxation:converge-quasi-static
         sim
         :energy-crit conv-criteria
         :oobf-crit conv-criteria
         :dt-scale 1d0;dt-scale
         :conv-steps 10000
         :substeps substeps
         :damping-factor 1d0;damping
         :post-iter-step
         (lambda (i e o)
           (save-conv-step sim output-dir *total-iter* 0 0d0 o e)
           (incf *total-iter* substeps)
           (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_conv_~5,'0d.vtk" i)) sim)
           (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_conv_nodes__~5,'0d.vtk" i)) sim)
           (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_conv_cells__~5,'0d.vtk" i)) sim)
           (funcall plotter sim)))
        (cl-mpm::finalise-loadstep sim)
        (cl-mpm::reset-grid (cl-mpm:sim-mesh sim))
        (setf (cl-mpm::sim-time sim) 0d0)
        (cl-mpm/dynamic-relaxation::reset-mp-velocity sim)
        (setf (cl-mpm::sim-velocity-algorithm sim) vel-algo)
        (change-class sim sim-type)))

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
    (setf (cl-mpm:sim-dt sim)
          (*
           (cl-mpm/setup:estimate-elastic-dt sim :dt-scale dt-scale)))
    (let* ((sim-time 0d0))
      (let ((total-iter 0)
            (dt-accumulator 0d0))
        (time (loop for steps from 0 to (round total-time dt)
                    while (and (cl-mpm::sim-run-sim sim)
                               (< sim-time total-time))
                    do
                       (let ()
                         (format t "Step ~d ~%" steps)
                         (cl-mpm/dynamic-relaxation::save-vtks sim output-dir steps)
                         (save-timestep sim output-dir total-iter :DYNAMIC)
                         (incf dt-accumulator dt)
                         (format t "Estimated substeps ~d ~%" (round dt (cl-mpm::sim-dt sim)))
                         (time
                          (loop while (> dt-accumulator 0d0)
                                do
                                   (progn
                                     (cl-mpm::update-sim sim)
                                     (incf total-iter)
                                     (decf dt-accumulator (cl-mpm::sim-dt sim))
                                     (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm::calculate-min-dt sim))))))
                         (funcall post-iter-step sim)
                         (funcall plotter sim)
                         (when save-vtk-loadstep
                           (save-vtks sim output-dir steps))
                         (swank.live:update-swank))))))
    result))
