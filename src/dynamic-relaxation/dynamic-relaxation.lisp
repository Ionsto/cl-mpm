(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)


(define-condition non-convergence-error (cl-mpm/errors:error-simulation)
  ((text :initarg :text :reader text)
   (ke-norm :initarg :ke-norm :reader ke-norm)
   (oobf-norm :initarg :oobf-norm :reader oobf-norm)))

(defparameter *run-convergance* t)
(declaim (notinline converge-quasi-static))
(defun converge-quasi-static (sim &key
                                    (energy-crit 1d-8)
                                    (oobf-crit 1d-8)
                                    (live-plot nil)
                                    (dt-scale 0.5d0)
                                    (substeps 50)
                                    (conv-steps 50)
                                    (post-iter-step nil)
                                    (convergance-criteria nil)
                                    (pic-update nil)
                                    (kinetic-damping nil)
                                    (damping-factor 1d-1))
  "Converge a simulation to a quasi-static solution via dynamic relaxation, whatever simulation is updated until it converges"
  (let ((current-vel (cl-mpm::sim-velocity-algorithm sim)))
    (when pic-update
      (setf (cl-mpm::sim-velocity-algorithm sim) :BLEND))
    (when (typep sim 'mpm-sim-dr)
      (setf (cl-mpm/dynamic-relaxation::sim-damping-scale sim)
            (if damping-factor
                damping-factor
                0d0)))
    (restart-case
        (%converge-quasi-static sim energy-crit oobf-crit live-plot dt-scale substeps conv-steps post-iter-step convergance-criteria kinetic-damping damping-factor)
      (continue ())
      (retry-convergence ()
        (%converge-quasi-static sim energy-crit oobf-crit live-plot dt-scale substeps conv-steps post-iter-step convergance-criteria kinetic-damping damping-factor)
        ))
    (setf (cl-mpm::sim-velocity-algorithm sim) current-vel)))

(defgeneric %converge-quasi-static (sim
                                    energy-crit
                                    oobf-crit
                                    live-plot
                                    dt-scale
                                    substeps
                                    conv-steps
                                    post-iter-step
                                    convergance-criteria
                                    kinetic-damping
                                    damping-factor
                                    )
 )
(defparameter *work* 0d0)
(defmethod %converge-quasi-static (sim
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
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (let* ((fnorm 0d0)
           (energy-total 0d0)
           (oobf 0d0)
           (target-time 1d-4)
           (estimated-t 1d-5)
           (total-step 0)
           (load 0d0)
           (converged nil))

      (setf (cl-mpm::sim-dt-scale sim) dt-scale)
      (setf (cl-mpm:sim-dt sim) (cl-mpm/setup::estimate-elastic-dt sim :dt-scale dt-scale))
      (if damping-factor
          (setf (cl-mpm:sim-damping-factor sim)
              (* damping-factor (cl-mpm/setup::estimate-critical-damping sim)))
          (setf (cl-mpm:sim-damping-factor sim) 0d0))

      (format t "Substeps ~D~%" substeps)
      (let ((full-load (list))
            (full-step (list))
            (full-enegy (list))
            (energy-list (list))
            (power-last 0d0)
            (power-current 0d0)
            (energy-first 0d0)
            (energy-last 0d0)
            ;; (dt-0 (cl-mpm/setup::estimate-elastic-dt sim :dt-scale dt-scale))
            )
        (setf *work* 0d0)
        (setf (cl-mpm::sim-stats-work sim) 0d0)
        (loop for i from 0 to conv-steps
              while (and *run-convergance*
                         (cl-mpm::sim-run-sim sim)
                         (not converged))
              do
                 (progn
                   (setf fnorm 0d0
                         load 0d0)
                   (time
                    (dotimes (j substeps)
                      (setf cl-mpm/penalty::*debug-force* 0d0)
                      (cl-mpm:update-sim sim)
                      ;; (when damping-factor
                      ;;   (setf (cl-mpm:sim-damping-factor sim) (* damping-factor (dr-estimate-damping sim))))
                      (let ((power (cl-mpm::sim-stats-power sim))
                            (energy (cl-mpm::sim-stats-energy sim)))
                        (incf *work* power)
                        (when kinetic-damping
                          (if (and
                               ;; (< (* power power-last) 0d0)
                               (> energy-last energy-first)
                               (> energy-last energy))
                              (progn
                                (format t "Peak found resetting KE - ~E ~E ~E~%" energy-first energy-last energy)
                                (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
                                (cl-mpm:iterate-over-mps
                                 mps
                                 (lambda (mp)
                                   (cl-mpm/fastmaths:fast-zero (cl-mpm/particle:mp-velocity mp))))
                                (setf power-last 0d0
                                      energy-first 0d0
                                      energy-last 0d0
                                      energy 0d0))
                              (progn
                                (setf energy-first energy-last)
                                (setf power-last power
                                      energy-last energy)))))))
                   (setf load cl-mpm/penalty::*debug-force*)
                   (setf energy-total (cl-mpm::sim-stats-energy sim))
                   (if (= *work* 0d0)
                       (setf fnorm 0d0)
                       (setf fnorm (abs (/ energy-total *work*))))
                   (setf oobf (cl-mpm::sim-stats-oobf sim))
                   (format t "Estimated dt ~E~%" (cl-mpm:sim-dt sim))
                   (format t "Conv step ~D - KE norm: ~E - Work: ~E - OOBF: ~E - Load: ~E~%" i fnorm *work* oobf
                           load)
                   (when (if convergance-criteria
                             (funcall convergance-criteria sim fnorm oobf)
                             (and
                              (< fnorm energy-crit)
                              (< oobf oobf-crit)))
                     (format t "Took ~D steps to converge~%" i)
                     (setf converged t))
                   (when post-iter-step
                     (funcall post-iter-step i fnorm oobf))
                   (swank.live:update-swank))))
      (when (not converged)
        (error (make-instance 'non-convergence-error
                              :text "System failed to converge"
                              :ke-norm fnorm
                              :oobf-norm 0d0)))
      (values load fnorm oobf))))

(defmethod %converge-quasi-static ((sim cl-mpm/mpi:mpm-sim-mpi)
                                   energy-crit
                                   oobf-crit
                                   live-plot
                                   dt-scale
                                   substeps
                                   conv-steps
                                   post-iter-step
                                   convergance-criteria
                                   kinetic-damping
                                    damping-factor
                                   )
  (setf *run-convergance* t)
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (let* ((fnorm 0d0)
           (oobf 0d0)
           (rank (cl-mpi:mpi-comm-rank))
           (load 0d0)
           (target-time 1d-4)
           (converged nil))
      (setf (cl-mpm:sim-dt sim)
            (cl-mpm/setup::estimate-elastic-dt sim :dt-scale dt-scale))
      (when (= rank 0)
        (format t "Substeps ~D~%" substeps))
      (setf *work* 0d0)
      (let ((full-load (list))
            (full-step (list))
            (energy-list (list))
            (energy-total 0d0)
            (power-last 0d0))
        (loop for i from 0 to conv-steps
              while (and *run-convergance*
                         (not converged))
              do
                 (progn
                   (dotimes (j substeps)
                     (setf cl-mpm/penalty::*debug-force* 0d0)
                     (cl-mpm:update-sim sim)
                     (let ((power (cl-mpm::sim-stats-power sim)))
                       (incf *work* power)
                       (when kinetic-damping
                         (if (< (* power-last power) 0d0)
                             (progn
                               (when (= rank 0)
                                 (format t "Peak found resetting KE~%"))
                               (cl-mpm:iterate-over-mps
                                mps
                                (lambda (mp)
                                  (cl-mpm/fastmaths:fast-zero (cl-mpm/particle:mp-velocity mp))))
                               (setf power-last 0d0))
                             (setf power-last power)))))

                   (setf load (cl-mpm/mpi::mpi-sum cl-mpm/penalty::*debug-force*))
                   (setf energy-total (cl-mpm::sim-stats-energy sim))
                   (if (= *work* 0d0)
                       (setf fnorm 0d0)
                       (setf fnorm (abs (/ energy-total *work*))))

                   (setf oobf (cl-mpm::sim-stats-oobf sim))
                   (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm::calculate-min-dt sim)))

                   (when (= 0 rank)
                     (format t "Estimated dt ~E~%" (cl-mpm:sim-dt sim)))


                   (when (= 0 rank)
                     (format t "Conv step ~D - KE norm: ~E - Work: ~E - OOBF: ~E - Load: ~E~%" i fnorm *work* oobf load))

                   (push energy-total energy-list)
                   (when (and (< fnorm energy-crit)
                              (< oobf oobf-crit)
                              (if convergance-criteria (funcall convergance-criteria sim) t))
                     (when (= 0 rank)
                       (format t "Took ~D steps to converge~%" i))
                     (setf converged t))
                   (when post-iter-step
                     (funcall post-iter-step i fnorm oobf))
                   (swank.live:update-swank))))
      (when (not converged)
        (error (make-instance 'non-convergence-error
                              :text "System failed to converge"
                              :ke-norm fnorm
                              :oobf-norm 0d0))
        (when (= 0 rank)
          (format t "System didn't converge~%"))
        )
      (values load fnorm oobf))))


(defgeneric dr-estimate-damping (sim))

(defmethod dr-estimate-damping (sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (dt cl-mpm:sim-dt))
      sim
    (let ((num 0d0)
          (denom 0d0))
      (setf
       num
       (cl-mpm::reduce-over-nodes
        mesh
        (lambda (node)
          (if (and (cl-mpm/mesh:node-active node)
                   (not (cl-mpm/mesh::node-agg node)))
              (cl-mpm/fastmaths:dot
               (cl-mpm/mesh:node-velocity node)
               (cl-mpm/fastmaths:fast-.-
                (cl-mpm/mesh::node-residual-prev node)
                (cl-mpm/mesh::node-residual node)))
              0d0))
        #'+))
      (setf
       denom
       (* dt
          (cl-mpm::reduce-over-nodes
           mesh
           (lambda (node)
             (if (and (cl-mpm/mesh:node-active node)
                      (not (cl-mpm/mesh::node-agg node)))
                 (* (cl-mpm/mesh:node-mass node)
                    (cl-mpm/fastmaths::dot
                     (cl-mpm/mesh::node-velocity node)
                     (cl-mpm/mesh::node-velocity node)))
                 0d0))
           #'+)))
      (if (> num 0d0)
          (if (= denom 0d0)
              0d0
              (* (sqrt 2) (sqrt (/ num denom))))
          0d0))))
(defmethod dr-estimate-damping ((sim cl-mpm/aggregate::mpm-sim-aggregated))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (dt cl-mpm:sim-dt))
      sim
    (let ((num 0d0)
          (denom 0d0))
      (setf
       num
       (cl-mpm::reduce-over-nodes
        mesh
        (lambda (node)
          (if (and (cl-mpm/mesh:node-active node)
                   (not (cl-mpm/mesh::node-agg node)))
              (cl-mpm/fastmaths:dot
               (cl-mpm/mesh:node-velocity node)
               (cl-mpm/fastmaths:fast-.-
                (cl-mpm/mesh::node-residual-prev node)
                (cl-mpm/mesh::node-residual node)))
              0d0))
        #'+))
      (setf
       denom
       (* dt
          (cl-mpm::reduce-over-nodes
           mesh
           (lambda (node)
             (if (and (cl-mpm/mesh:node-active node)
                      (not (cl-mpm/mesh::node-agg node)))
                 (* (cl-mpm/mesh:node-mass node)
                    (cl-mpm/fastmaths::dot
                     (cl-mpm/mesh::node-velocity node)
                     (cl-mpm/mesh::node-velocity node)))
                 0d0))
           #'+)))

      (when (cl-mpm/aggregate::sim-enable-aggregate sim)
        (loop for d from 0 below (cl-mpm/mesh::mesh-nd mesh)
              do
                 (let* ((res (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-residual d))
                        (res-prev (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-residual-prev d))
                        (E (cl-mpm/aggregate::sim-global-e sim))
                        (ma (cl-mpm/aggregate::sim-global-ma sim))
                        (vi (cl-mpm/aggregate::assemble-internal-vec sim #'cl-mpm/mesh::node-velocity d)))
                   (incf num (cl-mpm/fastmaths:dot
                              vi
                              (magicl:@ (magicl:transpose E)
                                        (cl-mpm/fastmaths::fast-.-
                                         res-prev
                                         res))))
                   (incf denom (* 0.5d0 (cl-mpm/utils::mtref (magicl:@ (magicl:transpose vi) ma vi) 0 0))))))
      (if (> num 0d0)
          (if (= denom 0d0)
              0d0
              (* (sqrt 2) (sqrt (/ num denom))))
          0d0))))

(defun compute-condition (sim)
  (/ 1d0 (expt (/ (dr-estimate-damping sim) 2) 2)))


(defun reset-mp-velocity (sim)
  (cl-mpm:iterate-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (cl-mpm/fastmaths::fast-zero (cl-mpm/particle:mp-velocity mp)))))


(defmethod cl-mpm::finalise-loadstep ((sim mpm-sim-dr))
  ;;DR algorithm requires that finalisation is called once
  (incf (cl-mpm::sim-time sim) (sim-dt-loadstep sim))
  ;; (cl-mpm::reset-node-displacement sim)
  (cl-mpm::new-loadstep sim))


(defun estimate-ul-enhancement (particle)
  (with-accessors ((df cl-mpm/particle::mp-deformation-gradient-increment))
      particle
    (let* ((grads (cl-mpm::gradient-push-forwards (list 1d0 1d0 1d0) df))
           (peak-grad (abs (reduce #'max (mapcar (lambda (x) (expt x 2)) grads)))))
      peak-grad)))





