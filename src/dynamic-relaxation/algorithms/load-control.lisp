(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)
(defun run-load-control (sim
                         &key (output-dir "./output/")
                           (loading-function nil)
                           (load-steps 10)
                           (initial-load 0d0)
                           (substeps 50)
                           (damping (sqrt 2d0))
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
