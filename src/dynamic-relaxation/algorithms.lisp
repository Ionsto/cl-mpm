(in-package :cl-mpm/dynamic-relaxation)



(defun save-timestep (sim output-dir step type)
  (with-open-file (stream (merge-pathnames "timesteps.csv" output-dir) :direction :output :if-exists :append)
    (format stream "~D,~f,~f,~f,~f,~f,~A,~f~%"
            step
            (cl-mpm::sim-time sim)
            (get-damage sim)
            (get-plastic sim)
            0d0
            0d0
            type
            0d0)))

(defun save-conv-preamble (output-dir)
  (with-open-file (stream (merge-pathnames output-dir "conv.csv") :direction :output :if-exists :supersede)
    (format stream "iter,step,plastic,damage,oobf,energy~%")))
(defun save-conv-step (sim output-dir total-iter step oobf energy)
  (with-open-file (stream (merge-pathnames output-dir "conv.csv") :direction :output :if-exists :append)
    (format stream "~D,~D,~f,~f,~f,~f~%" total-iter step (get-plastic sim) (get-damage sim)
            oobf energy)))

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
                  (/ 1d0 loadstep-dt)
                  (cl-mpm/mesh:node-mass n)
                  (cl-mpm/mesh::node-mass n)
                  (cl-mpm/fastmaths::mag-squared (cl-mpm/mesh::node-displacment n))
                  ))))))
    (if (= mass 0d0)
        0d0
        (/ energy mass))))


(declaim (notinline step-quasi-time))
(defun step-quasi-time (sim
                        global-step
                        &key (dt-scale 0.5d0)
                          (substeps 5)
                          (total-steps 0)
                          (damping 1d-1)
                          (output-dir "./output/")
                          (enable-damage t)
                          (enable-plastic t)
                          (plotter (lambda (sim)))
                          )
  (handler-case
      (let ((total-iter 0))
        (cl-mpm:iterate-over-mps
         (cl-mpm:sim-mps sim)
         (lambda (mp)
           (when (typep mp 'cl-mpm/particle::particle-damage)
             (setf (cl-mpm/particle::mp-enable-damage mp) nil))
           (when (typep mp 'cl-mpm/particle::particle-plastic)
             (setf (cl-mpm/particle::mp-enable-plasticity mp) nil))))
        (setf (cl-mpm:sim-enable-damage sim) nil)
        (cl-mpm/dynamic-relaxation:converge-quasi-static
         sim
         :oobf-crit 1d-2
         :energy-crit 1d-2
         :kinetic-damping t
         :dt-scale dt-scale
         :substeps substeps
         :conv-steps 1000
         :damping-factor damping
         :post-iter-step
         (lambda (i energy oobf)
           (save-conv-step sim output-dir (+ total-iter total-steps) global-step oobf energy)
           (incf total-iter substeps)))
        (cl-mpm:iterate-over-mps
         (cl-mpm:sim-mps sim)
         (lambda (mp)
           (when (typep mp 'cl-mpm/particle::particle-damage)
             (setf (cl-mpm/particle::mp-enable-damage mp) enable-damage))
           (when (typep mp 'cl-mpm/particle::particle-plastic)
             (setf (cl-mpm/particle::mp-enable-plasticity mp) enable-plastic))))
        (when enable-damage
          (cl-mpm/damage:calculate-damage sim (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim)))
        (setf (cl-mpm:sim-enable-damage sim) t)
        (let ((damage-prev (get-damage sim))
              (damage-eval t)
              (oobf-crit 1d-2)
              (energy-crit 1d-2)
              ;; (substeps 1)
              )
          (cl-mpm/dynamic-relaxation:converge-quasi-static
           sim
           :oobf-crit 1d-2
           :energy-crit 1d-2
           :convergance-criteria
           (lambda (sim fnorm oobf)
             (let* ((dn (get-damage sim))
                    (dconv (if (> dn 0d0)
                               (if (> damage-prev 0d0) (/ (- dn damage-prev) damage-prev) sb-ext:double-float-positive-infinity)
                               0d0
                               )))
               (format t "d-conv ~E~%" dconv)
               (setf damage-prev dn)
               (and
                (< fnorm energy-crit)
                (< oobf oobf-crit)
                (< dconv 1d-3)
                damage-eval
                ))
             )

           :kinetic-damping t
           :dt-scale dt-scale
           :substeps substeps
           :conv-steps 200
           :damping-factor 1d0
           :post-iter-step
           (lambda (i e o)
             (funcall plotter sim)
             (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d_~5,'0d.vtk" global-step i)) sim)
             (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d_~5,'0d.vtk" global-step i)) sim)
             (let ((true-intertia (true-intertial-criteria sim (sim-dt-loadstep sim))))
               ;; (if (= power 0))
               (format t "True intertia ~E~%" true-intertia))
             (let ((dev (and (< e energy-crit)
                             (< o oobf-crit))))
               (setf (cl-mpm:sim-enable-damage sim) dev)
               (setf damage-eval dev))
             (save-conv-step sim output-dir (+ total-iter total-steps) global-step o e)
             (incf total-iter substeps)
             )))
        (cl-mpm::finalise-loadstep sim)
        (save-timestep sim output-dir global-step :QUASI-STATIC)
        (values t total-iter))
    (error (c)
      (princ c)
      (cl-mpm::reset-loadstep sim)
      (values nil 0))))
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
                         (damping 1d-1)
                         (enable-damage t)
                         (enable-plastic t))
  (change-class sim 'cl-mpm/damage::mpm-sim-usl-damage)
  (setf (cl-mpm:sim-mass-scale sim) 1d0)
  (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm/setup:estimate-elastic-dt sim)))
  (setf (cl-mpm:sim-damping-factor sim)
        (* 1d-2 (cl-mpm/setup:estimate-critical-damping sim)))
  (cl-mpm:iterate-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (when (typep mp 'cl-mpm/particle::particle-damage)
       (setf (cl-mpm/particle::mp-enable-damage mp) enable-damage))
     (when (typep mp 'cl-mpm/particle::particle-plastic)
       (setf (cl-mpm/particle::mp-enable-plasticity mp) enable-damage))))
  (setf (cl-mpm:sim-enable-damage sim) enable-damage)
  (let* ((e-crit 1d-2)
         (oobf-crit 1d-2)
         (energy e-crit)
         (oobf oobf-crit)
         (work 0d0)
         (substeps (round target-time (cl-mpm:sim-dt sim))))
    (format t "Substeps ~D~%" substeps)
    (time (loop for step from 0 to max-steps
                while (and (cl-mpm::sim-run-sim sim)
                           (or (>= energy e-crit)
                               (>= oobf oobf-crit)))
                do
                   (let ((substeps (round target-time (cl-mpm:sim-dt sim))))
                     (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d_~5,'0d_real.vtk" global-step step)) sim)
                     (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d_~5,'0d_real.vtk" global-step step)) sim)
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

                     (let ((res (and (< energy 1d-1)
                                     (< oobf 1d-1))))
                       (setf (cl-mpm:sim-mass-scale sim)
                             (if res
                                 1d2
                                 1d0))
                       (format t "Accelerate real time ~A~%" res))
                     (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm/setup:estimate-elastic-dt sim)))
                     (format t "Residuals ~E ~E ~%" energy oobf)
                     (save-timestep sim output-dir global-step :DYNAMIC)
                     (funcall plotter sim)
                     (swank.live:update-swank)))))

  (change-class sim 'cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul-usl))

(declaim (notinline run-multi-stage))
(defun run-multi-stage (sim
                        &key (output-dir "./output/")
                          (post-conv-step (lambda (sim)))
                          (plotter (lambda (sim)))
                          (dt-scale 0.5d0)
                          (dt 1d0)
                          (steps 1d0)
                          (enable-damage t)
                          (enable-plastic t))
  (let ()
    (uiop:ensure-all-directories-exist (list output-dir))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
    (with-open-file (stream (merge-pathnames output-dir "./timesteps.csv") :direction :output :if-exists :supersede)
      (format stream "steps,time,damage,plastic,energy,oobf,step-type,mass~%"))
    (save-conv-preamble output-dir)
    (cl-mpm:iterate-over-mps
     (cl-mpm:sim-mps sim)
     (lambda (mp)
       (when (typep mp 'cl-mpm/particle::particle-damage)
         (setf (cl-mpm/particle::mp-enable-damage mp) nil))
       (when (typep mp 'cl-mpm/particle::particle-plastic)
         (setf (cl-mpm/particle::mp-enable-plasticity mp) nil))))

    (setf (cl-mpm:sim-mass-scale sim) 1d0)
    (setf (cl-mpm:sim-damping-factor sim)
          (* 1d-2 (cl-mpm/setup:estimate-critical-damping sim)))
    (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm/setup:estimate-elastic-dt sim)))
    (let ((total-iter 0)
          (substeps 50))
      (cl-mpm/dynamic-relaxation:converge-quasi-static
       sim
       :oobf-crit 1d-2
       :energy-crit 1d-2
       :kinetic-damping t
       :dt-scale dt-scale
       :conv-steps 1000
       :substeps substeps
       :damping-factor 1d-2
       :post-iter-step
       (lambda (i e o)
         ;; (plot sim)
         (save-conv-step sim output-dir total-iter 0 o e)
         (incf total-iter substeps)
         (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_conv_~5,'0d.vtk" i)) sim)
         (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_conv_nodes__~5,'0d.vtk" i)) sim)
         (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_conv_cells__~5,'0d.vtk" i)) sim)
         ))
      (cl-mpm::finalise-loadstep sim)
      (setf (cl-mpm::sim-time sim) 0d0)
      (cl-mpm/dynamic-relaxation::reset-mp-velocity sim)

      (cl-mpm:iterate-over-mps
       (cl-mpm:sim-mps sim)
       (lambda (mp)
         (when (typep mp 'cl-mpm/particle::particle-damage)
           (setf (cl-mpm/particle::mp-enable-damage mp) enable-damage))
         (when (typep mp 'cl-mpm/particle::particle-plastic)
           (setf (cl-mpm/particle::mp-enable-plasticity mp) enable-plastic))))
      (funcall post-conv-step sim)
      (loop for step from 1 to steps
            while (cl-mpm::sim-run-sim sim)
            do
               (let ((quasi-conv nil))
                 (setf (cl-mpm::sim-ghost-factor sim)
                       nil;; (* 1d9 1d-4)
                       )
                 (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim) dt)
                 (multiple-value-bind (conv inc-steps)
                     (step-quasi-time sim step
                                      :total-steps total-iter
                                      :plotter plotter
                                      :output-dir output-dir
                                      :dt-scale dt-scale
                                      :enable-damage enable-damage
                                      :enable-plastic enable-plastic)
                   (setf quasi-conv conv)
                   (incf total-iter inc-steps))

                 (unless quasi-conv
                   (cl-mpm/dynamic-relaxation::reset-mp-velocity sim)
                   (step-real-time sim step
                                   :output-dir output-dir
                                   :plotter plotter
                                   :dt-scale dt-scale
                                   :target-time (* 0.1d0 dt)
                                   :enable-damage enable-damage
                                   :enable-plastic enable-plastic))

                 (funcall plotter sim)
                 ;; (plot sim)
                 ;; (vgplot:title (format nil "Step ~D" step))
                 (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" step)) sim)
                 (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) sim)
                 (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d.vtk" step)) sim)
                 (swank.live:update-swank))))
    ))



;; (defun test-multi-step ()
;;   (setup :mps 2 :refine 0.5)
;;   (run-multi-stage
;;    :output-dir "./output/"
;;    :dt 0.5d0
;;    ;; :dt-scale (/ 0.8d0 (sqrt 1d0))
;;    ;; :load-steps 20
;;    )
;;   )




(defun run-load-control (sim
                         &key (output-dir "./output/")
                           (loading-function nil)
                           (load-steps 10)
                           (substeps 50)
                           (damping 1d-1)
                           (kinetic-damping t)
                           (adaptive-damping t)
                           (criteria 1d-3)
                           (post-conv-step (lambda (sim)))
                           (plotter (lambda (sim)))
                           (save-vtk-dr t)
                           (save-vtk-loadstep t)
                           (dt-scale 0.5d0))
  (uiop:ensure-all-directories-exist (list output-dir))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (save-conv-preamble output-dir)
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (let* ((load (cl-mpm/particle::mp-gravity (aref (cl-mpm:sim-mps sim) 0)))
           (total-iter 0))
      (unless loading-function
        (setf loading-function (lambda (percent)
                                 (format t "Loading factor ~E~%" percent)
                                 (cl-mpm:iterate-over-mps
                                  mps
                                  (lambda (mp)
                                    (setf (cl-mpm/particle:mp-gravity mp) (* load percent)))))))
      (when damping
        (setf (cl-mpm::sim-damping-factor sim)
              (* damping (cl-mpm/setup:estimate-critical-damping sim))))
      (setf (cl-mpm::sim-dt-scale sim) dt-scale)
      (loop for step from 1 to load-steps
            while (cl-mpm::sim-run-sim sim)
            do
               (progn
                 (defparameter *ke-last* 0d0)
                 ;; (pprint step)
                 (funcall loading-function (/ (float step) load-steps))
                 (let ((conv-steps 0)
                       (i 0))
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
                       (setf conv-steps (* substeps i))
                       (vgplot:title (format nil "Step ~D - substep ~D - KE ~E - OOBF ~E"  step (* i substeps) energy oobf))
                       (format t "Substep ~D~%" i)
                       (let ((i (+ 0 i)))
                         (when save-vtk-dr
                           (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d_~5,'0d.vtk" step i)) sim)
                           (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d_~5,'0d.vtk" step i)) sim)
                           (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d_~5,'0d.vtk" step i)) sim)))
                       (save-conv-step sim output-dir total-iter step oobf energy)
                       (funcall plotter sim)
                       (incf total-iter substeps)
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

            ))))
