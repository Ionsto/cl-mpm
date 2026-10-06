(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)

(defmethod %converge-quasi-static ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-mpi)
                                   energy-crit
                                   oobf-crit
                                   live-plot
                                   dt-scale
                                   substeps
                                   conv-steps
                                   post-iter-step
                                   convergance-criteria
                                   kinetic-damping
                                   damping-factor)
  (declare (double-float dt-scale))
  (setf *run-convergance* t)
  (cl-mpm/mpi::with-mpi-errors
      (with-accessors ((mps cl-mpm:sim-mps))
          sim
        (let* ((fnorm 0d0)
               (energy-total 0d0)
               (oobf 0d0)
               (load 0d0)
               (converged nil)
               (rank (cl-mpi:mpi-comm-rank)))
          (setf (cl-mpm::sim-dt-scale sim) dt-scale)
          (loop for i from 0 to conv-steps
                while (and *run-convergance*
                           (cl-mpm::sim-run-sim sim)
                           (not converged))
                do
                   (progn
                     (setf fnorm 0d0)
                     (optional-time
                      t
                      (dotimes (j substeps)
                        (cl-mpm:update-sim sim)))
                     (cl-mpm::update-dynamic-stats sim)
                     (setf energy-total (cl-mpm::sim-stats-energy sim))
                     (setf oobf (cl-mpm::sim-stats-oobf sim))
                     (when (= 0 rank)
                       (format t "Conv step ~D - KE norm: ~E - Work: ~E - OOBF: ~E~%" i fnorm *work* oobf))
                     (when (if convergance-criteria
                               (funcall convergance-criteria sim fnorm oobf)
                               (and
                                (< oobf oobf-crit)))
                       (when (= 0 rank)
                         (format t "Took ~D steps to converge~%" i))
                       (setf converged t))
                     (cl-mpm/mpi::with-mpi-errors
                         (when post-iter-step
                           (funcall post-iter-step i fnorm oobf)))
                     (swank.live:update-swank)))
          (when (not converged)
            (when (= 0 rank)
              (format t "System didn't converge~%"))
            (error (make-instance 'non-convergence-error
                                  :text "System failed to converge"
                                  :ke-norm fnorm
                                  :oobf-norm oobf)))
          (values load fnorm oobf)))))




(defmethod save-timestep ((sim cl-mpm/mpi:mpm-sim-mpi) output-dir step type)
  (when (uiop:file-exists-p (merge-pathnames output-dir "timesteps.csv"))
    (let ((str 
            (with-output-to-string (stream)
              (cl-mpm:sim-format sim stream
                                 "~D,~f,~f,~f,~f,~f,~f,~A,~f~%"
                                 step
                                 (cl-mpm::sim-time sim)
                                 (get-damage sim)
                                 (get-plastic sim)
                                 (cl-mpm::sim-stats-energy sim)
                                 (cl-mpm::sim-stats-oobf sim)
                                 (cl-mpm::sim-stats-work sim)
                                 type
                                 0d0))))
      (when (= (cl-mpi:mpi-comm-rank) 0)
        (with-open-file (stream (merge-pathnames "timesteps.csv" output-dir) :direction :output :if-exists :append)
          (format stream str))))))

(defmethod save-conv-step ((sim cl-mpm/mpi:mpm-sim-mpi) output-dir total-iter step real-time oobf energy)
  (when (uiop:file-exists-p (merge-pathnames output-dir "timesteps.csv"))
    (let ((str
            (with-output-to-string (stream)
              (cl-mpm:sim-format sim stream "~D,~D,~f,~f,~f,~f,~f~%" total-iter step real-time (get-plastic sim) (get-damage sim)
                                 (if (sb-ext:float-infinity-p oobf) 0d0 oobf) energy))))
      (when (= (cl-mpi:mpi-comm-rank) 0)
        (with-open-file (stream (merge-pathnames output-dir "conv.csv") :direction :output :if-exists :append)
          (format stream str))))))
