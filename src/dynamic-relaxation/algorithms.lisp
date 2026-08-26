(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)
(defparameter *total-iter* 0)


(defgeneric convergence-criteria (sim))

(defmethod convergence-criteria ((sim cl-mpm::mpm-sim))
  t)


(defgeneric convergence-check (sim)
  (:documentation "A check to see if our converged configuration has failed some extra criteria"))

(defmethod convergence-check ((sim cl-mpm::mpm-sim))
  t)




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
        (setf (cl-mpm/particle::mp-enable-plasticity mp) enable-plastic))))))

(defgeneric save-timestep-preamble (sim output-dir))
(defmethod save-timestep-preamble (sim output-dir)
  (with-open-file (stream (merge-pathnames output-dir "./timesteps.csv") :direction :output :if-exists :supersede)
    (cl-mpm:sim-format sim stream "steps,time,damage,plastic,energy,oobf,work,step-type,mass~%")))

(defmethod save-timestep ((sim cl-mpm:mpm-sim) output-dir step type)
  (when (uiop:file-exists-p (merge-pathnames output-dir "timesteps.csv"))
    (with-open-file (stream (merge-pathnames "timesteps.csv" output-dir) :direction :output :if-exists :append)
      (cl-mpm:sim-format sim stream "~D,~f,~f,~f,~f,~f,~f,~A,~f~%"
              step
              (cl-mpm::sim-time sim)
              (get-damage sim)
              (get-plastic sim)
              (cl-mpm::sim-stats-energy sim)
              (cl-mpm::sim-stats-oobf sim)
              (cl-mpm::sim-stats-work sim)
              type
              0d0))))


(defgeneric save-vtks (sim output-dir step &optional prefix))
(defmethod save-vtks ((sim cl-mpm:mpm-sim) output-dir step &optional prefix)
  (let ((pre (if prefix (format nil "_~A" prefix) "")))
    (cl-mpm:sim-format sim t "Save vtks ~D~%" step)
    (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim~A_~5,'0d.vtk" pre step)) sim)
    (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim~A_nodes_~5,'0d.vtk" pre step)) sim)
    (cl-mpm/output::save-vtk-cells (merge-pathnames output-dir (format nil "sim~A_cells_~5,'0d.vtk" pre step)) sim)
    (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim~A_p_~5,'0d.vtk" pre step)) sim)))

(defgeneric save-vtks-dr-step (sim output-dir trial-solve step iter))
(defmethod save-vtks-dr-step ((sim cl-mpm:mpm-sim) output-dir trial-solve step iter)
  (let ((post (format nil "~5,'0d_~5,'0d_~5,'0d.vtk" trial-solve step iter)))
    (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_step_~A.vtk" post)) sim)
    (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_step_nodes_~A.vtk" post)) sim)
    (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_step_cells_~A.vtk" post)) sim)
    (cl-mpm/penalty:save-vtk-penalties (merge-pathnames output-dir (format nil "sim_step_p_~A.vtk" post)) sim)))

(defgeneric save-conv-preamble (sim output-dir))
(defmethod save-conv-preamble ((sim cl-mpm:mpm-sim) output-dir)
  (with-open-file (stream (merge-pathnames output-dir "conv.csv") :direction :output :if-exists :supersede)
    (cl-mpm:sim-format sim stream "iter,step,real-time,plastic,damage,oobf,energy~%")))

(defmethod save-conv-step ((sim cl-mpm:mpm-sim) output-dir total-iter step real-time oobf energy)
  (unless (uiop:file-exists-p (merge-pathnames output-dir "conv.csv"))
    (save-conv-preamble sim output-dir))
  (with-open-file (stream (merge-pathnames output-dir "conv.csv") :direction :output :if-exists :append)
    (cl-mpm:sim-format sim stream "~D,~D,~f,~f,~f,~f,~f~%" total-iter step real-time (get-plastic sim) (get-damage sim)
            (if (sb-ext:float-infinity-p oobf) 0d0 oobf) energy)))

(defgeneric save-conv (sim output-dir iter))

(defmethod save-conv (sim output-dir iter)
  (let ((oobf (cl-mpm::sim-stats-oobf sim))
        (energy (/ (cl-mpm::sim-stats-energy sim) (cl-mpm::sim-stats-work sim)))
        (real-time (cl-mpm::sim-time sim)))
    (with-open-file (stream (merge-pathnames output-dir "conv.csv") :direction :output :if-exists :append)
      (cl-mpm:sim-format sim stream "~D,~D,~f,~f,~f,~f,~f~%" iter 0 real-time (get-plastic sim) (get-damage sim)
              oobf energy))))

(defun get-damage (sim)
  (cl-mpm::reduce-over-global-mps-sum
   sim
   (lambda (mp)
     (if (typep mp 'cl-mpm/particle::particle-damage)
         (*
          (cl-mpm/particle::mp-mass mp)
          (cl-mpm/particle:mp-damage mp))
         0d0))))
(defun get-plastic (sim)
  (cl-mpm::reduce-over-global-mps-sum
   sim
   (lambda (mp)
     (if (typep mp 'cl-mpm/particle::particle-plastic)
         (*
          (cl-mpm/particle::mp-mass mp)
          (cl-mpm/particle::mp-strain-plastic-vm mp))
         0d0))))




(defparameter *max-deformation-gradient* 10d0)



(defun generalised-staggered-solve (sim
                                    &key
                                      (output-dir "./output/")
                                      (crit 1d-3)
                                      (dt-scale 1d0)
                                      (substeps 10)
                                      (post-iter-step (lambda (i e o)))
                                      (enable-damage t)
                                      (enable-plastic t)
                                      (max-damage-inc 0.6d0)
                                      (max-plastic-inc 1d0)
                                      (max-deformation-gradient 10d0)
                                      (damping 1d0)
                                      (staggered-steps 500)
                                      (sub-conv-steps 50)
                                      (save-vtk-dr t)
                                      (convergence-criteria nil)
                                      (stagger-damage nil))
  (let* ((damage-prev (get-damage sim))
         (damage damage-prev)
         ;; (crit (cl-mpm/dynamic-relaxation::sim-convergence-critera sim))
         ;; (crit (cl-mpm/dynamic-relaxation::sim-convergence-critera sim))
         (oobf-crit   crit)
         (energy-crit 1d0)
         (damage-crit crit)
         (dconv damage-crit)
         (total-i 0))
    (setf
     (cl-mpm/dynamic-relaxation::sim-convergence-critera sim) crit
     (cl-mpm:sim-enable-damage sim) nil)
    (set-mp-plastic-damage
     sim
     :enable-damage enable-damage
     :enable-plastic enable-plastic)
    (unless stagger-damage
      (setf (cl-mpm:sim-enable-damage sim) stagger-damage))
    (let ((additional-conv nil))
      (loop for i from 0 to 100
            while (not additional-conv)
            do
               (progn
                 (loop for stagger-i from 0 to staggered-steps
                       while
                       (or
                        (>= dconv damage-crit)
                        (>= (cl-mpm::sim-stats-oobf sim) crit))
                       do (progn
                            (let ((iv 0))
                              (refine-mesh sim)
                              (cl-mpm/dynamic-relaxation:converge-quasi-static
                               sim
                               :kinetic-damping nil
                               :oobf-crit oobf-crit
                               :energy-crit energy-crit
                               :dt-scale dt-scale
                               :substeps substeps
                               :convergance-criteria
                               (lambda (sim f o)
                                 (if convergence-criteria
                                     (funcall convergence-criteria sim)
                                     (and
                                      (<= o (cl-mpm/dynamic-relaxation::sim-convergence-critera sim))
                                      (convergence-criteria sim))))
                               :conv-steps sub-conv-steps
                               :damping-factor damping
                               :post-iter-step
                               (lambda (i e o)
                                 (incf iv)
                                 ;; (when nil
                                 ;;   (check-deformation-gradient sim)
                                 ;;   (check-true-inertia sim))
                                 (check-deformation-gradient sim :max-deformation-gradient max-deformation-gradient)
                                 (check-damage-increment sim :max-damage-inc max-damage-inc)
                                 (check-plastic-increment sim :max-plastic-inc max-plastic-inc)
                                 (funcall post-iter-step i e o)))
                              (if enable-damage
                                  (progn
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
                                                     ;; (save-conv-step sim output-dir *total-iter* global-step 0d0 (cl-mpm::sim-stats-oobf sim) 0d0)
                                                     (incf *total-iter*)
                                                     (check-damage-increment sim :max-damage-inc max-damage-inc)
                                                     ))
                                        (setf (cl-mpm:sim-enable-damage sim) nil))
                                      (cl-mpm:update-sim sim)
                                      (cl-mpm::update-dynamic-stats sim)
                                      (setf dconv dconv-1)))
                                  (setf dconv 0d0))
                              ;; (if (and
                              ;;      enable-damage
                              ;;      (typep sim 'cl-mpm/damage::mpm-sim-damage))
                              ;;     (progn
                              ;;       (let ((fast-trial-conv oobf-crit)
                              ;;             (damage-iter t)
                              ;;             (save-update nil))
                              ;;         (unless (cl-mpm:sim-enable-damage sim)
                              ;;           (setf (cl-mpm:sim-enable-damage sim) t)
                              ;;           (cl-mpm/damage::calculate-damage
                              ;;            sim
                              ;;            (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim)))
                              ;;         (let ((dconv-1 (compute-damage-delta sim)))
                              ;;           (setf dconv dconv-1)
                              ;;           (cl-mpm:sim-format sim t "step ~D/~D - d-conv ~E ~A~%" stagger-i 0 dconv stagger-damage)
                              ;;           (when stagger-damage
                              ;;             (format t "Staggering damage~%")
                              ;;             (loop for d from 1 to 100
                              ;;                   while damage-iter
                              ;;                   do (setf (cl-mpm:sim-enable-damage sim) t)
                              ;;                      (if (<= d 5)
                              ;;                          (cl-mpm/damage::calculate-damage sim (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim))
                              ;;                          (dotimes (i 5)
                              ;;                            (cl-mpm/damage::calculate-damage sim (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim))))
                              ;;                      (when stagger-damage
                              ;;                        (setf (cl-mpm:sim-enable-damage sim) nil))
                              ;;                      (setf damage (get-damage sim))
                              ;;                      (setf dconv (compute-damage-delta sim))
                              ;;                      (setf damage-prev damage)
                              ;;                      (unless (>= dconv damage-crit)
                              ;;                        (setf damage-iter nil))
                              ;;                      (cl-mpm:sim-format sim t "step ~D/~D - d-conv ~E - ~E~%" stagger-i d dconv fast-trial-conv)
                              ;;                      (cl-mpm:sim-format sim t "Damage ~E - prev damage ~E ~%" damage damage-prev)
                              ;;                      ;; (cl-mpm:sim-format sim t "step ~D/~D - d-conv ~E~%" stagger-i d dconv)
                              ;;                      (check-damage-increment sim :max-damage-inc max-damage-inc)
                              ;;                      (setf damage-prev damage))
                              ;;             (cl-mpm:update-sim sim)
                              ;;             (cl-mpm::update-dynamic-stats sim)
                              ;;             (setf dconv dconv-1)))
                              ;;         (when save-update
                              ;;           (incf iv)
                              ;;           (funcall post-iter-step iv (cl-mpm::sim-stats-energy sim) (cl-mpm::sim-stats-oobf sim)))))
                              ;;     ;; No damage evolution -> instantly satisfy dconv
                              ;;     (setf dconv 0d0))
                              )))
                 (setf additional-conv (convergence-check sim)))))
    (when (>= dconv damage-crit)
      (cl-mpm:sim-format sim t "Failed to converge damage during stagger~%")
      (error (make-instance 'non-convergence-error
                            :text "Failed to converge damage during stagger"
                            :ke-norm 0d0
                            :oobf-norm 0d0)))))

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
;; (defun step-quasi-time (sim
;;                         global-step
;;                         &key (dt-scale 1d0)
;;                           (substeps 50)
;;                           (total-steps 0)
;;                           (damping 1d0)
;;                           (sub-conv-steps 50)
;;                           (conv-criteria 1d-3)
;;                           (conv-criteria-damage 1d-3)
;;                           (output-dir "./output/")
;;                           (save-vtk-dr t)
;;                           (enable-damage t)
;;                           (enable-plastic t)
;;                           (max-damage-inc 0.6d0)
;;                           (max-plastic-inc 10d0)
;;                           (stagger-damage nil)
;;                           (plotter (lambda (sim))))
;;   (let ((total-i 0))
;;     (handler-bind
;;         ((cl-mpm/errors:error-simulation
;;            (lambda (c)
;;              (format t "Handled error~%")
;;              (cl-mpm/utils::kill-errors)
;;              (princ c)
;;              (cl-mpm::reset-loadstep sim)
;;              (return-from step-quasi-time (values nil 0)))))
;;         (progn
;;           (let* ((oobf-crit   conv-criteria)
;;                  (energy-crit 1d0)
;;                  (damage-crit conv-criteria-damage)
;;                  (dconv (if enable-damage damage-crit 0d0))
;;                  (stagger-iters 0))
;;             (reset-mp-velocity sim)
;;             (set-mp-plastic-damage sim :enable-damage enable-damage :enable-plastic enable-plastic)
;;             (setf (cl-mpm:sim-enable-damage sim) nil)
;;             ;; (setf (cl-mpm:sim-enable-damage sim) enable-damage)
;;             (let ((alt-conv-crit nil))
;;               (loop for ac from 0 to 10
;;                     while (not alt-conv-crit)
;;                     do
;;                        (progn
;;                          (loop for stagger-i from 0 to 100
;;                                while (or
;;                                       (< stagger-i 2)
;;                                       (>= dconv damage-crit)
;;                                       (>= (cl-mpm::sim-stats-oobf sim) oobf-crit))
;;                                do
;;                                   (progn
;;                                     (format t "Stagger iter ~D~%" stagger-i)
;;                                     (refine-mesh sim)
;;                                     (cl-mpm/dynamic-relaxation:converge-quasi-static
;;                                      sim
;;                                      :oobf-crit oobf-crit
;;                                      :energy-crit energy-crit
;;                                      :dt-scale dt-scale
;;                                      :substeps substeps
;;                                      :conv-steps sub-conv-steps
;;                                      :convergance-criteria
;;                                      (lambda (sim f o)
;;                                        (setf dconv (compute-damage-delta sim))
;;                                        (and
;;                                         (<= o (cl-mpm/dynamic-relaxation::sim-convergence-critera sim))
;;                                         ;; (convergence-criteria sim)
;;                                         (< dconv damage-crit)))
;;                                      :damping-factor damping
;;                                      :post-iter-step
;;                                      (lambda (i e o)
;;                                        (funcall plotter sim)
;;                                        (incf total-i)
;;                                        (when save-vtk-dr
;;                                          (when (= (mod i 1) 0)
;;                                            (save-vtks-dr-step sim output-dir global-step *trial-iter* total-i)))
;;                                        ;; (cl-mpm:sim-format sim t "Def crit ~E~%" (compute-max-deformation sim))
;;                                        ;; (let ((max-def (compute-max-deformation sim)))
;;                                        ;;   (when (> max-def 10d0)
;;                                        ;;     (cl-mpm:sim-format sim t "Deformation gradient criteria exceeded~%")
;;                                        ;;     (error (make-instance 'non-convergence-error
;;                                        ;;                           :text "Deformation gradient J exceeded"
;;                                        ;;                           :ke-norm 0d0
;;                                        ;;                           :oobf-norm 0d0)))
;;                                        ;;   (let ((true-intertia (true-intertial-criteria sim (sim-dt-loadstep sim))))
;;                                        ;;     (cl-mpm:sim-format sim t "True intertia ~E~%" true-intertia)
;;                                        ;;     (save-conv-step sim output-dir *total-iter* global-step
;;                                        ;;                     0d0
;;                                        ;;                     o
;;                                        ;;                     max-def)
;;                                        ;;     (when (sim-inertia-criteria sim)
;;                                        ;;         (when (> true-intertia (sim-inertia-criteria sim))
;;                                        ;;           (cl-mpm:sim-format sim t "Inertia criteria exceeded~%")
;;                                        ;;           (error (make-instance 'error-inertia-criteria
;;                                        ;;                                 :text "True inertia exceeded"
;;                                        ;;                                 :inertia-norm true-intertia))))))
;;                                        ;; (when max-plastic-inc
;;                                        ;;   (let ((plastic-inc (plastic-increment-criteria sim)))
;;                                        ;;     (format t "Plastic inc criteria ~E~%" plastic-inc)
;;                                        ;;     (when (> plastic-inc max-plastic-inc)
;;                                        ;;       (cl-mpm:sim-format sim t "Damage criteria failed~%")
;;                                        ;;       (error (make-instance 'error-plastic-criteria
;;                                        ;;                             :max-plastic-inc plastic-inc)))))
;;                                        ;; (unless stagger-damage
;;                                        ;;   (let ((damage-inc (damage-increment-criteria sim))
;;                                        ;;         (damage (get-damage sim)))
;;                                        ;;     (cl-mpm:sim-format sim t "Damage ~E - inc ~E~%" damage damage-inc)
;;                                        ;;     (when (> damage-inc max-damage-inc)
;;                                        ;;       (cl-mpm:sim-format sim t "Damage criteria failed ~E~%" damage-inc)
;;                                        ;;       (error (make-instance 'error-damage-criteria
;;                                        ;;                             :text "Damage criteria exeeded"
;;                                        ;;                             :max-damage-inc 0d0)))))
;;                                        (incf *total-iter* substeps)
;;                                        (refine-mesh sim)))
;;                                     (setf (cl-mpm:sim-enable-damage sim) enable-damage)
;;                                     (when save-vtk-dr
;;                                       (save-vtks-dr-step sim output-dir global-step *trial-iter* total-i))
;;                                     (incf stagger-iters)))
;;                          (setf alt-conv-crit (convergence-check sim))
;;                          (unless alt-conv-crit
;;                            (setf (cl-mpm::sim-stats-oobf sim) oobf-crit)))))
;;             (when (or (> (cl-mpm::sim-stats-oobf sim) oobf-crit)
;;                       (> dconv damage-crit))
;;               (cl-mpm:sim-format sim t "Staggered solve didn't converge ~E ~E~%" dconv (cl-mpm::sim-stats-oobf sim))
;;               (error (make-instance 'non-convergence-error
;;                                     :text "Staggered solve didn't converge ~E ~E"
;;                                     :ke-norm  dconv
;;                                     :oobf-norm  (cl-mpm::sim-stats-oobf sim))))

;;             (when save-vtk-dr
;;               (save-vtks-dr-step sim output-dir global-step *trial-iter* total-i)
;;               ;;   (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_step_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
;;               ;; (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_step_nodes_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
;;               )
;;             (cl-mpm::finalise-loadstep sim)
;;             (save-timestep sim output-dir global-step :QUASI-STATIC)
;;             (values t stagger-iters)))
;;       ))
;;   )

;; (declaim (notinline step-quasi-time))
;; (defun step-quasi-time (sim
;;                         global-step
;;                         &key (dt-scale 1d0)
;;                           (substeps 50)
;;                           (total-steps 0)
;;                           (damping 1d0)
;;                           (sub-conv-steps 50)
;;                           (conv-criteria 1d-3)
;;                           (conv-criteria-damage 1d-3)
;;                           (output-dir "./output/")
;;                           (save-vtk-dr t)
;;                           (enable-damage t)
;;                           (enable-plastic t)
;;                           (max-damage-inc 0.6d0)
;;                           (max-plastic-inc 10d0)
;;                           (stagger-damage nil)
;;                           (plotter (lambda (sim))))
;;   (let ((total-i 0))
;;     (handler-bind
;;         ((cl-mpm/errors:error-simulation
;;            (lambda (c)
;;              (format t "Handled error~%")
;;              (cl-mpm/utils::kill-errors)
;;              (princ c)
;;              (cl-mpm::reset-loadstep sim)
;;              (return-from step-quasi-time (values nil 0)))))
;;         (progn
;;           (let* ((damage-prev (get-damage sim))
;;                  (damage damage-prev)
;;                  (oobf-crit   conv-criteria)
;;                  (energy-crit
;;                    1d0
;;                    ;conv-criteria
;;                               )
;;                  (damage-crit conv-criteria-damage)
;;                  (dconv (if enable-damage damage-crit 0d0))
;;                  (inertia 0d0)
;;                  (intertia-crit 1d-3)
;;                  (stagger-iters 0))
;;             (reset-mp-velocity sim)
;;             (set-mp-plastic-damage sim :enable-damage enable-damage :enable-plastic enable-plastic)
;;             (setf (cl-mpm:sim-enable-damage sim) nil)
;;             ;; (setf (cl-mpm:sim-enable-damage sim) enable-damage)
;;             (let ((alt-conv-crit nil)
;;                   (staggered stagger-damage)
;;                   )
;;               (loop for ac from 0 to 10
;;                     while (not alt-conv-crit)
;;                     do
;;                        (progn
;;                          (loop for stagger-i from 0 to 100
;;                                while (or
;;                                       (< stagger-i 2)
;;                                       (>= dconv damage-crit)
;;                                       (>= (cl-mpm::sim-stats-oobf sim) oobf-crit))
;;                                do
;;                                   (progn
;;                                     (format t "Stagger iter ~D~%" stagger-i)
;;                                     (when (and (> stagger-i 1)
;;                                                (not staggered))
;;                                       (setf stagger-damage nil))
;;                                     (refine-mesh sim)
;;                                     (cl-mpm/dynamic-relaxation:converge-quasi-static
;;                                      sim
;;                                      :oobf-crit oobf-crit
;;                                      :energy-crit energy-crit
;;                                      :dt-scale dt-scale
;;                                      :substeps substeps
;;                                      :conv-steps sub-conv-steps
;;                                      :damping-factor damping
;;                                      :post-iter-step
;;                                      (lambda (i e o)
;;                                        (funcall plotter sim)
;;                                        (incf total-i)
;;                                        (when save-vtk-dr
;;                                          (when (= (mod i 1) 0)
;;                                            (save-vtks-dr-step sim output-dir global-step *trial-iter* total-i)))
;;                                        (cl-mpm:sim-format sim t "Def crit ~E~%" (compute-max-deformation sim))
;;                                        (let ((max-def (compute-max-deformation sim)))
;;                                          (when (> max-def 10d0)
;;                                            (cl-mpm:sim-format sim t "Deformation gradient criteria exceeded~%")
;;                                            (error (make-instance 'non-convergence-error
;;                                                                  :text "Deformation gradient J exceeded"
;;                                                                  :ke-norm 0d0
;;                                                                  :oobf-norm 0d0)))
;;                                          (let ((true-intertia (true-intertial-criteria sim (sim-dt-loadstep sim))))
;;                                            (cl-mpm:sim-format sim t "True intertia ~E~%" true-intertia)
;;                                            (save-conv-step sim output-dir *total-iter* global-step
;;                                                            0d0
;;                                                            o
;;                                                            ;; true-intertia
;;                                                            max-def
;;                                                            )
;;                                            (setf inertia true-intertia)
;;                                            (when (sim-inertia-criteria sim)
;;                                                (when (> true-intertia (sim-inertia-criteria sim))
;;                                                  (cl-mpm:sim-format sim t "Inertia criteria exceeded~%")
;;                                                  (error (make-instance 'error-inertia-criteria
;;                                                                        :text "True inertia exceeded"
;;                                                                        :inertia-norm true-intertia))))))
;;                                        (when max-plastic-inc
;;                                          (let ((plastic-inc (plastic-increment-criteria sim)))
;;                                            (format t "Plastic inc criteria ~E~%" plastic-inc)
;;                                            (when (> plastic-inc max-plastic-inc)
;;                                              (cl-mpm:sim-format sim t "Damage criteria failed~%")
;;                                              (error (make-instance 'error-plastic-criteria
;;                                                                    :max-plastic-inc plastic-inc)))))
;;                                        (unless stagger-damage
;;                                          (let ((damage-inc (damage-increment-criteria sim)))
;;                                            (cl-mpm:sim-format sim t "Damage ~E - prev damage ~E - inc ~E~%" damage damage-prev damage-inc)
;;                                            (when (> damage-inc max-damage-inc)
;;                                              (cl-mpm:sim-format sim t "Damage criteria failed ~E~%" damage-inc)
;;                                              (error (make-instance 'error-damage-criteria
;;                                                                    :text "Damage criteria exeeded"
;;                                                                    :max-damage-inc 0d0)))))
;;                                        ;; (save-conv-step sim output-dir *total-iter* global-step 0d0 o 0d0)
;;                                        (incf *total-iter* substeps)
;;                                        (refine-mesh sim)))

;;                                     (let ((fast-trial-conv oobf-crit)
;;                                           (damage-iter t))
;;                                       (unless (typep sim 'cl-mpm/damage::mpm-sim-damage)
;;                                         (setf dconv 0d0))
;;                                       (when
;;                                           (and
;;                                            enable-damage
;;                                            (typep sim 'cl-mpm/damage::mpm-sim-damage))
;;                                         ;; (loop for k from 0 to 10
;;                                         ;;       while (and
;;                                         ;;              (<= fast-trial-conv (sqrt oobf-crit))
;;                                         ;;              ;damage-iter
;;                                         ;;              )
;;                                         ;;       do)
;;                                         (loop for d from 0 to 100
;;                                               while (and
;;                                                      ;; (<= fast-trial-conv (sqrt oobf-crit))
;;                                                      damage-iter)
;;                                               do
;;                                                  (setf (cl-mpm:sim-enable-damage sim) t)
;;                                                  (if (<= d 5)
;;                                                      (cl-mpm/damage::calculate-damage sim (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim))
;;                                                      (dotimes (i 5)
;;                                                        (cl-mpm/damage::calculate-damage sim (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim))))
;;                                                  (when stagger-damage
;;                                                    (setf (cl-mpm:sim-enable-damage sim) nil))
;;                                                  (setf damage (get-damage sim))
;;                                                  (setf dconv (compute-damage-delta sim))
;;                                                  (setf damage-prev damage)
;;                                                  (when (< dconv damage-crit)
;;                                                    (setf damage-iter nil))
;;                                                  (cl-mpm:sim-format sim t "step ~D/~D - d-conv ~E - ~E~%" stagger-i d dconv fast-trial-conv)
;;                                                  (let ((damage-inc (damage-increment-criteria sim)
;;                                                                    ;; (compute-max-damage-energy-crit sim)
;;                                                                    ))
;;                                                    (cl-mpm:sim-format sim t "Damage ~E - prev damage ~E - inc ~E~%" damage damage-prev damage-inc)
;;                                                    ;; (cl-mpm:sim-format sim t "Damage energy criterion ~E~%" (compute-max-damage-energy-crit sim))
;;                                                    (when (> damage-inc max-damage-inc)
;;                                                      (cl-mpm:sim-format sim t "Damage criteria failed ~E~%" damage-inc)
;;                                                      (when save-vtk-dr
;;                                                        (save-vtks-dr-step sim output-dir global-step *trial-iter* total-i)
;;                                                        ;; (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_step_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
;;                                                        ;; (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_step_nodes_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
;;                                                        ;; (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_step_cells_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
;;                                                        )
;;                                                      (error (make-instance 'error-damage-criteria
;;                                                                            :text "Damage criteria exeeded"
;;                                                                            :max-damage-inc 0d0))))
;;                                                  (save-conv-step sim output-dir *total-iter* global-step
;;                                                                  0d0
;;                                                                  fast-trial-conv
;;                                                                  0d0)
;;                                                  (incf *total-iter*)
;;                                                  (incf total-i))
;;                                         (when t
;;                                           (dotimes (i 1)
;;                                             (cl-mpm:update-sim sim))
;;                                           (cl-mpm::update-dynamic-stats sim)
;;                                           (setf fast-trial-conv (cl-mpm::sim-stats-oobf sim))
;;                                           (cl-mpm:sim-format sim t "fast trial ~E~%" fast-trial-conv))))

;;                        (setf damage-prev damage)
;;                        (when save-vtk-dr
;;                          (save-vtks-dr-step sim output-dir global-step *trial-iter* total-i)
;;                          ;; (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_step_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
;;                          ;; (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_step_nodes_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
;;                          ;; (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_step_cells_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
;;                          )
;;                        (incf stagger-iters)
;;                        ;; (when (> inertia intertia-crit)
;;                        ;;   (error (make-instance 'non-convergence-error
;;                        ;;                         :text "Quasi-time inertia was too large"
;;                        ;;                         :ke-norm inertia
;;                        ;;                         :oobf-norm oobf-crit)))
;;                        ;; (when (= (mod stagger-i 10) 0)
;;                        ;;   (convergence-check sim)
;;                        ;;   )
;;                        ))
;;                          (setf alt-conv-crit (convergence-check sim))
;;                          (unless alt-conv-crit
;;                            (setf (cl-mpm::sim-stats-oobf sim) oobf-crit)))))
;;             (when (or (> (cl-mpm::sim-stats-oobf sim) oobf-crit)
;;                       (> dconv damage-crit))
;;               (cl-mpm:sim-format sim t "Staggered solve didn't converge ~E ~E~%" dconv (cl-mpm::sim-stats-oobf sim))
;;               (error (make-instance 'non-convergence-error
;;                                     :text "Staggered solve didn't converge ~E ~E"
;;                                     :ke-norm  dconv
;;                                     :oobf-norm  (cl-mpm::sim-stats-oobf sim))))

;;             (when save-vtk-dr
;;               (save-vtks-dr-step sim output-dir global-step *trial-iter* total-i)
;;               ;;   (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_step_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
;;               ;; (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_step_nodes_~5,'0d_~5,'0d_~5,'0d.vtk" global-step *trial-iter* total-i)) sim)
;;               )
;;             (cl-mpm::finalise-loadstep sim)
;;             (save-timestep sim output-dir global-step :QUASI-STATIC)
;;             (values t stagger-iters)))
;;       ))
;;   )
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
                                       :target-time (* 0.1d0 dt-loadstep)
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
                           (max-plastic-inc nil)
                           (conv-steps 50)
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
  (save-conv-preamble sim output-dir)
  (cl-mpm/output::save-simulation-parameters (merge-pathnames output-dir "settings.json") sim)
  (funcall pre-step)
  ;; (save-conv-step sim output-dir 0 0 0d0 0d0)
  (defparameter *total-iter* 0)
  (save-vtks sim output-dir 0)
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (let* ((load (cl-mpm:sim-gravity sim)))
      (unless loading-function
        (setf loading-function (lambda (percent)
                                 (cl-mpm:sim-format sim t "Loading factor ~E~%" percent)
                                 (setf (cl-mpm:sim-gravity sim)
                                       (* load percent)))))
      (setf (cl-mpm::sim-dt-scale sim) dt-scale)
      (let ((t0 (get-internal-real-time))
            (tunits internal-time-units-per-second))
        (loop for step from 1 to load-steps
              while (cl-mpm::sim-run-sim sim)
              do
                 (progn
                   (cl-mpm:sim-format sim t "Load step ~D~%" step)
                   (defparameter *ke-last* 0d0)
                   ;; (pprint step)
                   (let* ((percent (/ (float step) load-steps)))
                     ;; (declare (double-float percent initial-load))
                     (funcall loading-function (+ initial-load (* (- 1d0 initial-load) percent))))
                   (let ((iter-conv-steps 0)
                         (i 0))
                     ;; (when damping
                     ;;   (setf (cl-mpm::sim-damping-factor sim)
                     ;;         (* damping (cl-mpm/setup:estimate-critical-damping sim))))
                     (setf (cl-mpm/dynamic-relaxation::sim-damping-scale sim) damping)
                     (time
                      (generalised-staggered-solve
                       sim
                       :crit criteria
                       :enable-damage enable-damage
                       :enable-plastic enable-plastic
                       :dt-scale dt-scale
                       :substeps substeps
                       :damping damping
                       :sub-conv-steps conv-steps
                       :max-plastic-inc max-plastic-inc
                       :post-iter-step
                       (lambda (i energy oobf)
                         (funcall post-iter-step i energy oobf)
                         (setf iter-conv-steps (* substeps i))
                         (cl-mpm:sim-format sim t "Substep ~D~%" i)
                         (let ((i (+ 0 i)))
                           (when save-vtk-dr
                             (save-vtks-dr-step sim output-dir 0 step i)))
                         (incf *total-iter* substeps)
                         (save-conv-step sim output-dir *total-iter* step 0d0 oobf energy)
                         (funcall plotter sim))))))

                 (cl-mpm::finalise-loadstep sim)
                 (funcall post-conv-step sim)
                 (when save-vtk-loadstep
                   (save-vtks sim output-dir step))
                 (funcall plotter sim)
                 (swank.live:update-swank)
              )))))


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
           (funcall plotter temp-sim)
           ))
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


(defun reset-mp-velocity (sim)
  (cl-mpm:iterate-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (cl-mpm/fastmaths::fast-zero (cl-mpm/particle::mp-velocity-gradient mp))
     (cl-mpm/fastmaths:fast-zero (cl-mpm/particle::mp-velocity mp)))))

(defun reset-nominal-displacement (sim)
  (cl-mpm:iterate-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (cl-mpm/fastmaths:fast-zero (cl-mpm/particle::mp-displacement mp)))))

(defun elastic-static-solution (sim &key (crit 1d-3)
                                      (elastic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static)
                                      (substeps 10)
                                      (dt-scale 1d0)
                                      (post-iter-step
                                       (lambda (i e o)
                                         (cl-mpm:sim-format sim t "Static step ~D - ~E~%" i o))))
  "Take a sim, switch it to implicit quasi-static, disable elastic and plastic and converge switch the class back to the original then enable damage and plastic"
  (let ((vel-algo (cl-mpm::sim-velocity-algorithm sim))
        (sim-type (class-of sim)))
    (unless (subtypep (type-of sim) elastic-solver)
      (format t "Changed class from ~A to ~A~%" (type-of sim) elastic-solver)
      (change-class sim elastic-solver))
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
                           (sub-conv-steps 50)
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
                           (stagger-damage nil)
                           (save-vtk-loadstep t)
                           (max-adaptive-steps 5)
                           (min-adaptive-steps 0)
                           (compute-zero-loadstep nil)
                           (adaption-constant 2)
                           (max-damage-inc 0.3d0)
                           (max-plastic-inc nil)
                           (max-deformation-gradient 10d0)
                           (min-damage-inc 0d0)
                           (dt-scale 1d0))
  (uiop:ensure-all-directories-exist (list output-dir))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (save-conv-preamble sim output-dir)
  (cl-mpm/output::save-simulation-parameters (merge-pathnames output-dir "settings.json") sim)
  (funcall pre-step)
  ;; (save-conv-step sim output-dir 0 0 0d0 0d0)
  (setf lparallel:*debug-tasks-p* nil)
  (cl-mpm:iterate-over-mps (cl-mpm::sim-mps sim) (lambda (mp)))
  (defparameter *total-iter* 0)
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (let* ((load (cl-mpm:sim-gravity sim))
           (finished t))
      (unless loading-function
        (setf loading-function (lambda (percent)
                                 (cl-mpm:sim-format sim t "Loading factor ~E~%" percent)
                                 (setf (cl-mpm:sim-gravity sim)
                                       (* load percent)))))
      (setf (cl-mpm::sim-dt-scale sim) dt-scale)
      (let* ((initial-load-step-size (/ (- 1d0 initial-load) load-steps))
             (prev-steps-easy (list t t))
             (prev-step-iter 0)
             (load-step-size initial-load-step-size)
             (current-load initial-load)
             (current-adaptivity 0))
        (loop for step from (if compute-zero-loadstep 0 1)
              while (and (cl-mpm::sim-run-sim sim)
                         (< (+ 1d-15 current-load) 1d0))
              do
                 (progn
                   (cl-mpm:sim-format sim t "Load step ~D~%" step)
                   (defparameter *ke-last* 0d0)
                   (let ((conv-steps 0)
                         (i-total 0))
                     (flet ((trial-solve (trial-step)
                              (handler-case
                                  (progn
                                    (cl-mpm:sim-format sim t "Load applied ~E~%" (min 1d0 (+ current-load load-step-size)))
                                    (funcall loading-function (+ current-load load-step-size))
                                    (time
                                     (;cl-mpm/dynamic-relaxation:converge-quasi-static
                                      generalised-staggered-solve
                                      sim
                                      :crit criteria
                                      :damping damping
                                      :output-dir output-dir
                                      :dt-scale dt-scale
                                      :substeps substeps
                                      :sub-conv-steps sub-conv-steps
                                      :enable-damage enable-damage
                                      :enable-plastic enable-plastic
                                      :max-damage-inc max-damage-inc
                                      :max-plastic-inc max-plastic-inc
                                      :max-deformation-gradient max-deformation-gradient
                                      :stagger-damage stagger-damage
                                      :post-iter-step
                                      (lambda (i-g energy oobf)
                                        (incf i-total)
                                        (let ((i i-total))
                                          (funcall post-iter-step i energy oobf)
                                          (setf conv-steps (* substeps i))
                                          (cl-mpm:sim-format sim t "Substep ~D~%" i)
                                          (let ((i (+ 0 i)))
                                            (when save-vtk-dr
                                              (save-vtks-dr-step sim output-dir step trial-step i)))

                                          (incf *total-iter* substeps)
                                          (save-conv-step sim output-dir *total-iter* step 0d0 oobf energy)
                                          (funcall plotter sim)))))
                                    (cl-mpm:sim-format sim t "Generalised solve passed~%")
                                    t)
                                (cl-mpm/errors:error-simulation (c)
                                  (princ c)
                                  (save-conv-step sim output-dir *total-iter* step 0d0 (cl-mpm::sim-stats-oobf sim) 0d0)
                                  (when save-vtk-dr
                                    (save-vtks-dr-step sim output-dir step trial-step i-total))
                                  (incf *total-iter* substeps)
                                  (cl-mpm::reset-loadstep sim)
                                  (values nil 0)
                                  )
                                ;; (error (c)
                                ;;   (cl-mpm:sim-format sim t "A non simulation error was thrown!")
                                ;;   (princ c)
                                ;;   (cl-mpm::reset-loadstep sim)
                                ;;   (cl-mpm::sim-run-sim sim) nil
                                ;;   (return-from run-adaptive-load-control)
                                ;;   (values nil 0))
                                )))
                       (let ((quasi-conv nil))
                         (loop for i from 0 to max-adaptive-steps
                               while (not quasi-conv)
                               do (progn
                                    (setf load-step-size
                                          (min (- 1d0 current-load)
                                               (* initial-load-step-size (expt adaption-constant (- current-adaptivity)))))

                                    (when (= step 0)
                                      (setf load-step-size 0d0))
                                    (cl-mpm:sim-format sim t "Load step size ~E adapted by ~D~%" load-step-size current-adaptivity)
                                    (setf quasi-conv (trial-solve i))
                                    (unless quasi-conv
                                      (cl-mpm:sim-format sim t "Failed conv, adapting~%" )
                                      (incf current-adaptivity))
                                    (if (or (= current-adaptivity max-adaptive-steps))
                                        (progn
                                          (loop-finish))))
                               finally (progn
                                         (cl-mpm:sim-format sim t "Finished with ~D adaptions- conv ~A~%" (- i 1) quasi-conv)
                                         (when (> min-damage-inc 0d0)
                                           (format t "Easy damage inc criteria ~E true ~E~%" min-damage-inc (damage-increment-criteria sim)))
                                         (when (and (= i 1)
                                                    (if (> min-damage-inc 0d0)
                                                        (< (damage-increment-criteria sim) min-damage-inc)
                                                        t))
                                           (setf (nth (mod prev-step-iter (length prev-steps-easy)) prev-steps-easy) t)
                                           (cl-mpm:sim-format sim t "Potential adaption easy steps ~A~%" prev-steps-easy)
                                           (when (every #'identity prev-steps-easy)
                                             (setf current-adaptivity
                                                   (max min-adaptive-steps
                                                        (- current-adaptivity 1)))))
                                         ))
                         (unless quasi-conv;(>= current-adaptivity max-adaptive-steps)
                           (cl-mpm:sim-format sim t "Solve failed completly~%" )
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

(defgeneric copy-instance (object &rest initargs &key &allow-other-keys)
  (:documentation "Makes and returns a shallow copy of OBJECT.

  An uninitialized object of the same class as OBJECT is allocated by
  calling ALLOCATE-INSTANCE.  For all slots returned by
  CLASS-SLOTS, the returned object has the
  same slot values and slot-unbound status as OBJECT.

  REINITIALIZE-INSTANCE is called to update the copy with INITARGS.")
  (:method ((object standard-object) &rest initargs &key &allow-other-keys)
    (let* ((class (class-of object))
           (copy (allocate-instance class)))
      (dolist (slot-name (mapcar #'sb-mop:slot-definition-name (sb-mop:class-slots class)))
        (when (slot-boundp object slot-name)
          (setf (slot-value copy slot-name)
                (slot-value object slot-name))))
      (apply #'reinitialize-instance copy initargs))))

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

(defun run-load-control-timed (sim
                               &key (output-dir "./output/")
                                 (loading-function nil)
                                 (load-steps 10)
                                 (initial-load 0d0)
                                 (substeps 50)
                                 (damping 1d0)
                                 (kinetic-damping nil)
                                 (criteria 1d-3)
                                 (conv-steps 50)
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
  (save-conv-preamble sim output-dir)
  (cl-mpm/output::save-simulation-parameters (merge-pathnames output-dir "settings.json") sim)
  (funcall pre-step)
  ;; (save-conv-step sim output-dir 0 0 0d0 0d0)
  (defparameter *total-iter* 0)
  (when save-vtk-loadstep
    (save-vtks sim output-dir 0))
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (let* ((load (cl-mpm:sim-gravity sim)))
      (unless loading-function
        (setf loading-function (lambda (percent)
                                 (cl-mpm:sim-format sim t "Loading factor ~E~%" percent)
                                 (setf (cl-mpm:sim-gravity sim)
                                       (* load percent)))))
      (setf (cl-mpm::sim-dt-scale sim) dt-scale)
      (let ((t0 (get-internal-real-time))
            (tunits internal-time-units-per-second))
        (declare (function loading-function))
        (loop for step from 1 to load-steps
              while (cl-mpm::sim-run-sim sim)
              do
                 (progn
                   (cl-mpm:sim-format sim t "Load step ~D~%" step)
                   (defparameter *ke-last* 0d0)
                   (let* ((percent (/ (float step) load-steps)))
                     (funcall loading-function (+ initial-load (* (- 1d0 initial-load) percent))))
                   (let ()
                     (setf (cl-mpm/dynamic-relaxation::sim-damping-scale sim) damping)
                     (time
                      (cl-mpm/dynamic-relaxation:converge-quasi-static
                       sim
                       :oobf-crit criteria
                       :energy-crit 1d0
                       :dt-scale dt-scale
                       :substeps substeps
                       :kinetic-damping nil
                       :damping-factor (sqrt 2d0)
                       :conv-steps conv-steps
                       :post-iter-step
                       (lambda (i energy oobf)
                         (funcall post-iter-step i energy oobf))))))
                 (cl-mpm::finalise-loadstep sim)
                 ;; (funcall post-conv-step sim)
                 ;; (funcall plotter sim)
                 ;; (swank.live:update-swank)
              )
        (let ((dt (/ (- (get-internal-run-time) t0) tunits)))
          dt
          (when save-vtk-loadstep
            (save-vtks sim output-dir 1))
          dt)))))
