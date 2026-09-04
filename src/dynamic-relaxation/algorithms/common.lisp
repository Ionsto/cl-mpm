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
  "Setup the mp-specific flags for plasticity and damage evolution"
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
  (let ((post (format nil "~5,'0d_~5,'0d_~5,'0d" trial-solve step iter)))
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
                                  (setf dconv 0d0)))))
                 (setf additional-conv (convergence-check sim)))))
    (when (>= dconv damage-crit)
      (cl-mpm:sim-format sim t "Failed to converge damage during stagger~%")
      (error (make-instance 'non-convergence-error
                            :text "Failed to converge damage during stagger"
                            :ke-norm 0d0
                            :oobf-norm 0d0)))))

;; (defun assemble-global-vec (sim accessor dim &optional (res nil))
;;   (declare (function accessor))
;;   (let* ((active-nodes (sim-agg-nodes-fd sim))
;;          (ndof (length active-nodes))
;;          (v (if res res (cl-mpm/utils::arb-matrix ndof 1))))
;;     (cl-mpm::iterate-over-nodes-array
;;      active-nodes
;;      (lambda (n)
;;        (setf (varef v (cl-mpm/mesh::node-agg-fd n)) (varef (funcall accessor n) dim))))
;;     v))

;; (defun anderson-acceleration (sim)
;;   (with-accessors ((mps (cl-mpm:sim-mps sim))))
;;   (let (
;;         (size (length ))
;;         )
;;     ))
