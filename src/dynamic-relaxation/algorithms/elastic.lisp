(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)
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
