(defpackage :cl-mpm/examples/ice/cliff-stability
  (:use :cl
   :cl-mpm/example))
(in-package :cl-mpm/examples/ice/cliff-stability)

;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)

(defparameter *angle* 38d0)
(defparameter *angle-r* 20d0)
(defparameter *angle-psi* 0d0)
(defparameter *rt* 1d0)
(defparameter *rc* 0.9d0)
(defparameter *enable-plastic-damage* nil)
(defparameter *delay-time* 1d4)
(defparameter *delay-exponent* 4d0)
(defparameter *enable-viscosity* nil)
(defparameter *viscosity* 1d13)
(defparameter *gf* 10000d0)
(defparameter *pd-oversize* 1d-3)

(defparameter *length-scaler* 2d0)
(defparameter *length-scale* nil)

(defparameter *biot-coefficent* 0d0)

(defparameter *ductility* 10d0)
;; (defparameter *tensile-strength* 0.1185d6)
(defparameter *tensile-strength* 0.1d6)


(defun est-shear-from-angle (angle angle-r rc)
  (let* ((angle-plastic (* angle (/ pi 180)))
         (angle-plastic-residual (* angle-r (/ pi 180))))
    (- 1d0
       (* (- 1d0 rc)
          (/ (tan angle-plastic-residual)
             (tan angle-plastic))))))

(declaim (notinline plot-domain))
(defun plot-domain ()
  (when *sim*
    (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
           (h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
           (ms-x (first ms))
           (ms-y (second ms)))
      (vgplot:format-plot t "set object 1 rect from 0,0 to ~f,~f fc rgb 'blue' fs transparent solid 0.5 noborder behind" ms-x *water-height*))
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :trial nil
     ;; :colour-func (lambda (mp) (sqrt (cl-mpm/constitutive::voigt-j2 (cl-mpm/utils:deviatoric-voigt (cl-mpm/particle:mp-stress mp)))))
     ;; :colour-func #'cl-mpm/particle::mp-strain-plastic-vm
     ;; :colour-func (lambda (mp) (/ (cl-mpm/particle:mp-mass mp)
     ;;                              (cl-mpm/particle:mp-volume mp)))
     :colour-func #'cl-mpm/particle::mp-damage
     ;; :colour-func (lambda (mp) (cl-mpm/utils::varef  (cl-mpm/particle::mp-stress mp ) 1))
     ))
  )


(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-ice-brittle) dt)
  (with-accessors ((undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                   (y cl-mpm/particle::mp-damage-y-local)
                   (strain cl-mpm/particle::mp-strain)
                   (trial-strain cl-mpm/particle::mp-trial-strain)
                   (plastic-strain cl-mpm/particle::mp-strain-plastic)
                   (ps-vm cl-mpm/particle::mp-strain-plastic-vm)
                   (damage cl-mpm/particle:mp-damage)
                   (pressure cl-mpm/particle::mp-pressure)
                   (init-stress cl-mpm/particle::mp-initiation-stress)
                   (ybar cl-mpm/particle::mp-damage-ybar)
                   (angle cl-mpm/particle::mp-friction-angle)
                   (de cl-mpm/particle::mp-elastic-matrix)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (E cl-mpm/particle::mp-e)
                   (nu cl-mpm/particle::mp-nu)
                   (j cl-mpm/particle::mp-deformation-jacobian-strain)
                   (pd-inc cl-mpm/particle::mp-plastic-damage-evolution))
      mp
    (declare (double-float E ps-vm angle pressure))
    (progn
      (let* ((ps-y (sqrt (* E (expt ps-vm 2))))
             ;; (stress (cl-mpm/constitutive:linear-elastic-mat strain de))
             (stress-pressure
               (cl-mpm/fastmaths:fast-.+
                undamaged-stress
                ;; stress
                (cl-mpm/utils:voigt-eye (*
                                         0d0
                                         j
                                         (- pressure)
                                         )))))
        (setf y
              (*
               (+
                (if pd-inc ps-y 0d0)
                ;; (cl-mpm/damage::tensile-energy-norm strain e de)
                (cl-mpm/damage::criterion-mohr-coloumb-rankine-stress-tensile stress-pressure angle)
                ;; (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile stress-pressure angle)
                ;; (cl-mpm/damage::criterion-mohr-coloumb-rankine-stress-tensile stress angle)
                )))))))



(defparameter *penalty-epsilon-scale* 1d0)

(defun setup (&key (refine 1) (mps 2)
                (pressure-condition t)
                (cryo-static t)
                (hydro-static nil)
                (elastic-static nil)
                (friction 0d0)
                (ice-height 800d0)
                (bench-length 0d0)
                (aspect 1)
                (floatation-ratio 0.9)
                (bench-extra-cut 0d0)
                (slope 0d0)
                (use-penalty t)
                (stick-base t)
                (multigrid-refines 0))
  (let* ((density 918d0)
         (water-density 1028d0)
         (mesh-resolution (/ 10d0 refine))
         (h-fine mesh-resolution)
         (offset (* mesh-resolution (if use-penalty 2 0)))
         (end-height ice-height)
         (ice-length (* ice-height aspect))
         (start-height (+ ice-height (* slope ice-length)))
         (ice-height end-height)
         (floating-point (* ice-height (/ density water-density)))
         (water-level (* floating-point floatation-ratio))
         ;; (water-level (+ floating-point extra-cliff-height))
         (datum (+ water-level offset))
         ;; (datum (* (round datum mesh-resolution) mesh-resolution))
         (domain-size (list (+ ice-length (* 2 ice-height)) (* start-height 2)))
         (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
         (block-size (list ice-length (max start-height end-height)))
         (E 1d9)
         )
    (pprint element-count)
    (defparameter *sim* nil)
    (defparameter *water-height* datum)
    (defparameter *offset* offset)
    (defparameter *ice-length* ice-length)
    (setf *sim* (cl-mpm/setup::make-simple-sim mesh-resolution element-count
                                               :sim-type
                                               ;; Took [TODO] seconds
                                               ;; 'cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul
                                               ;; Took 57 seconds
                                               'cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic
                                               ;; 'cl-mpm/dynamic-relaxation::mpm-sim-octree-damage-quasi-static
                                               :args-list
                                               (list
                                                :enable-fbar t
                                                :enable-aggregate t
                                                :split-factor (* 1.2d0 (sqrt 2) (/ 1d0 mps))
                                                :enable-split nil
                                                :max-split-depth 6
                                                ;; :refinement multigrid-refines
                                                )))

    (setf mesh-resolution (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
    (setf h-fine (* mesh-resolution (expt 2 (- multigrid-refines))))
    (let* ((angle *angle*)
           (E 1d9)
           ;; (init-stress (* 0.1185d6 1d0))
           (init-stress *tensile-strength*)
           (length-scale
             (if *length-scale*
                 *length-scale*
                 (* h-fine *length-scaler*)))
           (gf (/ (* 10 length-scale (expt init-stress 2)) E))
           ;; (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress E))
           ;; (ductility 10d0)
           (ductility *ductility*)
           (oversize-ratio (cl-mpm/damage::compute-oversize-factor (- 1d0 *pd-oversize*) ductility))
           )
      (format t "Ice length ~F~%" ice-length)
      (format t "Water height ~F~%" water-level)
      (format t "True Water height ~F~%" (- datum offset))
      (format t "Cliff height ~F~%" (- (+ offset ice-height) datum))
      (format t "Init stress ~E - ductilty ~E ~%" init-stress ductility)
      (format t "Mesh size ~F~%" mesh-resolution)
      (format t "Estimated lc ~E~%" length-scale)
      (format t "Estimated ductility ~E~%" ductility)
      (format t "Estimated oversize factor ~E~%" oversize-ratio)
      (format t "Init stress ~E~%" init-stress)
      (let* ((rt *rt*)
             (rc *rc*))
        (cl-mpm:add-mps
         *sim*
         (cl-mpm/setup:make-block-mps
          (list 0 offset)
          block-size
          (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
          density
          'cl-mpm/particle::particle-ice-delayed
          ;; 'cl-mpm/particle::particle-ice-brittle
          :E E
          :nu 0.24d0
          :kt-res-ratio rt
          :kc-res-ratio rc
          :initiation-stress init-stress
          :friction-angle (cl-mpm/utils:deg-to-rad angle)
          :residual-friction (cl-mpm/utils:deg-to-rad *angle-r*)

          :psi (cl-mpm/utils::deg-to-rad *angle-psi*)

          :ductility ductility

          :local-length length-scale
          :delay-time *delay-time*
          :delay-exponent *delay-exponent*
          :oversize (- 1d0 *pd-oversize*)
          :enable-plasticity t
          :enable-damage t
          :residual-strength 1d0
          :biot-coeff *biot-coefficent*

          :enable-viscosity nil
          :viscosity *viscosity*
          )))
      (when hydro-static
        (cl-mpm/setup::initialise-stress-self-weight-vardatum
         *sim*
         (lambda (pos) datum)
         :k-x 1d0
         :k-z 1d0
         :scaler (lambda (pos) (/ water-density density))))
      (when cryo-static
        (cl-mpm/setup::initialise-stress-self-weight-vardatum
         *sim*
         (lambda (pos)
           (let ((alpha (- 1d0 (/ (abs (-
                                       ice-length
                                       (cl-mpm/utils::varef pos 0))) ice-length))))
             (+ offset
                (* alpha end-height)
                (* (- 1d0 alpha) start-height))))
         :k-x 1d0
         :k-z 1d0
         ))
      (when elastic-static
        (cl-mpm/setup::initialise-stress-self-weight-vardatum
         *sim*
         (lambda (pos)
           (let ((alpha (- 1d0 (/ (abs (-
                                        ice-length
                                        (cl-mpm/utils::varef pos 0))) ice-length))))
             (+ offset
                (* alpha end-height)
                (* (- 1d0 alpha) start-height))))
         :index 0))

      (unless (= start-height end-height)
        (cl-mpm/setup::remove-sdf *sim*
                                  (lambda (p)
                                    (cl-mpm/setup::plane-point-point-sdf
                                     p
                                     (cl-mpm/utils:vector-from-list (list 0d0 (+ offset start-height) 0d0))
                                     (cl-mpm/utils:vector-from-list (list ice-length (+ offset end-height) 0d0))))
                                  :refine 2))


      (let ((cutout (+ (- ice-height water-level) bench-extra-cut))
            (cutback bench-length))
        (when (> cutback 0d0)
          (cl-mpm/setup:remove-sdf
           *sim*
           (cl-mpm/setup::rectangle-sdf (list (first block-size) (+ offset ice-height ice-height))
                                        (list cutback (+ cutout ice-height)))))))
    (cl-mpm/setup:setup-bcs
     *sim*
     :left '(0 nil 0)
     :bottom '(nil 0 0))

    (setf (cl-mpm:sim-mass-scale *sim*) 1d0)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 0.1d0
             (sqrt 1d0)
             (cl-mpm/setup:estimate-critical-damping *sim*)))
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)
    (setf (cl-mpm::sim-allow-mp-split *sim*) t)
    (setf (cl-mpm::sim-allow-mp-damage-removal *sim*) nil)
    (setf (cl-mpm::sim-mp-damage-removal-instant *sim*) nil)
    (setf (cl-mpm::sim-mp-damage-removal-criteria *sim*) 0.99d0)
    (setf (cl-mpm::sim-ghost-factor *sim*) nil)
    (setf (cl-mpm::sim-velocity-algorithm *sim*) :BLEND)
    (setf (cl-mpm:sim-dt *sim*) (* 0.5d0 (cl-mpm/setup:estimate-elastic-dt *sim*)))
    (setf (cl-mpm::sim-enable-damage *sim*) nil)
    (setf *run-sim* t)
    (defparameter *water-bc*
      (if pressure-condition
          (cl-mpm/buoyancy::make-bc-buoyancy-clip
           *sim*
           datum
           water-density
           (lambda (pos datum)
             (>= (cl-mpm/utils:varef pos 1) (* mesh-resolution 0)))
           :visc-damping 1d0)
          (cl-mpm/buoyancy::make-bc-buoyancy-body
           *sim*
           datum
           water-density
           (lambda (pos) t))))

    (cl-mpm:add-bcs-force-list
     *sim*
     *water-bc*)

    (let ((domain-half (* 0.5d0 (first domain-size)))
          (friction friction))
      (defparameter *floor-bc*
        (cl-mpm/penalty::make-bc-penalty-distance-point
         *sim*
         (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
         (cl-mpm/utils:vector-from-list (list
                                         domain-half
                                         offset
                                         0d0))
         (* domain-half 1.1d0)
         (* E *penalty-epsilon-scale*)
         friction
         ;; 0.1d0
         0d0
         )))

    (when use-penalty
      (cl-mpm:add-bcs-force-list
       *sim*
       *floor-bc*))
    (unless use-penalty
      (if stick-base
          (cl-mpm/setup:setup-bcs
           *sim*
           :left '(0 nil 0)
           :bottom '(0 0 0))
          (cl-mpm/setup:setup-bcs
           *sim*
           :left '(0 nil 0)
           :bottom '(nil 0 0))))
    ;; (cl-mpm:add-bcs-force-list
    ;;  *sim*
    ;;  *bc-erode*
    ;;  )
    (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
    )
  (cl-mpm/output:add-mp-output *sim* :SCALAR "water-pressure" #'cl-mpm/particle::mp-pressure)
  (cl-mpm/output::add-mp-output
   *sim* :SCALAR "damage-increment-error"
   (lambda (mp)
     (if (typep mp 'cl-mpm/particle::particle-damage)
         (with-accessors ((damage cl-mpm/particle::mp-damage)
                          (damage-prev cl-mpm/particle::mp-damage-prev-trial)
                          (inc cl-mpm/particle::mp-damage-increment)
                          (mass cl-mpm/particle::mp-mass))
             mp
           (- damage damage-prev))
         0d0)))

  (cl-mpm/output:add-mp-output
     *sim*
     :SCALAR
     "damage-tcs-c"
     #'cl-mpm/particle::mp-damage-compression)
    (cl-mpm/output:add-mp-output
     *sim*
     :SCALAR
     "damage-tcs-t"
     #'cl-mpm/particle::mp-damage-tension)
    (cl-mpm/output:add-mp-output
     *sim*
     :SCALAR
     "damage-tcs-s"
     #'cl-mpm/particle::mp-damage-shear)
    (cl-mpm/output:add-mp-output
     *sim*
     :SCALAR
     "current-effective-angle"
     (lambda (mp)
       (if (typep mp 'cl-mpm/particle::particle-ice-brittle)
           (* (/ 180 pi) (atan (* (/ (- 1d0 (cl-mpm/particle::mp-damage-shear mp))
                                     (- 1d0 (cl-mpm/particle::mp-damage-compression mp)))
                                  (tan (cl-mpm/particle::mp-phi mp)))))
           0d0)))
    (cl-mpm/output:add-mp-output
     *sim*
     :SCALAR
     "current-cohesion"
     (lambda (mp)
       (if (typep mp 'cl-mpm/particle::particle-ice-brittle)
           (*
            (if (> (cl-mpm/constitutive::voight-trace (cl-mpm/particle::mp-stress mp)) 0d0)
                (- 1d0 (cl-mpm/particle::mp-damage-tension mp))
                (- 1d0 (cl-mpm/particle::mp-damage-compression mp)))
            (max 0d0 (cl-mpm/particle::mp-c mp)))
           0d0)))
  )

(defun save-stabilty-data (output-dir sim stable height floatation notch)
  (let ((filename (merge-pathnames (format nil "data_~A_~A_~A.json" height floatation notch) output-dir)))
    (str:to-file
     filename
     (jonathan:to-json
      (list
       :stable stable
       :height height
       :time (cl-mpm::sim-time sim)
       :floatation floatation
       :notch notch
       )))))

(defparameter *early-exit-flag* nil)

(define-condition early-exit-condition (error)
  ((stable
    :initform t
    :accessor stable
    :initarg :stable)))

(defmethod cl-mpm/dynamic-relaxation::convergence-check ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-ul))
  (if (> (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim) 0d0)
      (progn
        ;; (pprint "Check velocity")
        ;; (cl-mpm/dynamic-relaxation::check-max-velocity sim :max-velocity 1d0)
        (let ((v (cl-mpm/dynamic-relaxation::max-velocity-criteria sim (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim))))
          (format t "Max velocity ~E~%" v)
          (when (> v 1d0)
            (format t "Maximum velocity error~%")
            (error (make-instance 'early-exit-condition
                                  :stable nil))))
        (let ((d (cl-mpm::reduce-over-global-mps-max
                  sim
                  (lambda (mp)
                    (cl-mpm/fastmaths:mag
                      (cl-mpm/particle::mp-displacement mp))))))
          (format t  "Max disp: ~E~%" d)
          (when (> d 1d0)
            (format t  "Maximum global displacement reached ~E~%" d)
            (error (make-instance 'early-exit-condition
                                  :stable nil))))
        (let ((y (cl-mpm::reduce-over-global-mps-max
                  sim
                  #'cl-mpm/particle::mp-damage-ybar)))
          (when (> y 0d0)
            (when (< y *tensile-strength*)
              (defparameter *early-exit-flag* t)
              (format t  "Maximum global value of y-bar is ~E - yield strength is ~E~%" y *tensile-strength*)
              (error (make-instance 'early-exit-condition
                                    :stable t)))))
        t)
      t))

(defmethod cl-mpm/dynamic-relaxation::damage-increment-criteria ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-ul))
  (cl-mpm/dynamic-relaxation::damage-increment-criteria-mp sim)
  ;; (cl-mpm/dynamic-relaxation::compute-max-damage-energy-crit sim)
  ;; (cl-mpm/dynamic-relaxation::compute-max-damage-energy-crit-mp sim)
  )
(defun compute-energy-release (sim mp)
  (if (typep mp 'cl-mpm/particle::particle-damage)
      (with-accessors ((volume cl-mpm/particle::mp-volume)
                       (stress cl-mpm/particle::mp-stress)
                       (stress-u cl-mpm/particle::mp-undamaged-stress)
                       (strain cl-mpm/particle::mp-strain)
                       (damage-inc cl-mpm/particle::mp-damage-increment)
                       (j cl-mpm/particle::mp-deformation-jacobian-strain)
                       (k cl-mpm/particle::mp-history-stress)
                       (k-n cl-mpm/particle::mp-history-stress-n))
          mp
        (let ((k-temp k)
              (e-d-n 0d0)
              (e-d 0d0)
              (damage-n 0d0)
              (damage-n1 0d0)
              )
          (declare (double-float volume j))
          (setf k k-n)
          (cl-mpm/damage::compute-damage mp)
          (cl-mpm/particle::post-damage-step mp (cl-mpm:sim-dt sim))
          (setf e-d-n (* 0.5d0 volume j (cl-mpm/fastmaths:dot stress strain)))
          (setf damage-n (cl-mpm/particle::mp-damage mp))
          (setf k k-temp)
          (cl-mpm/damage::compute-damage mp)
          (cl-mpm/particle::post-damage-step mp (cl-mpm:sim-dt sim))
          (setf e-d (* 0.5d0 volume j (cl-mpm/fastmaths:dot stress strain)))
          (setf damage-n1 (cl-mpm/particle::mp-damage mp))
          ;; (format t "~E ~E - ~E ~E~%" damage-n damage-n1 e-d e-d-n)
          (if (> e-d 0d0)
              (abs (/ (abs (- e-d e-d-n)) e-d))
              0d0)))
      0d0))

(defun stability-qt-test ()
  (cl-mpm/utils:set-workers 16)
  (let* ((heights (list
                   500d0
                   ;; 200d0
                   ;; 300d0
                   ;; ;; 400d0
                   ;; 500d0
                   ;; 600d0
                   ))
         (floatations (list
                       0.5d0
                       ;; 0d0
                       ;; 0.5d0
                       ;; 0.9d0
                       ;; 1d0
                       ))
         (cliff-step 20d0)
         (density 918d0)
         (water-density 1028d0)
         )
    (defparameter *stability* (make-array (list (length heights) (length floatations)) :initial-element nil
                                                                                       :element-type t))
    (let ((stability-dir (merge-pathnames (format nil "./analysis_scripts/ice/ice-cliff-stability/data-cliff-stability/"))))
      (ensure-directories-exist stability-dir)
      (defparameter *heights* heights)
      (defparameter *floatations* floatations)
      (loop for hi from 0
            for height in heights
            do
               (let* ((res t)
                      (floating-point (/ density water-density))
                      (floating-cliff (- height (* height floating-point))))
                 (loop for fi from 0
                       for flotation in floatations
                       ;; for i from 0 to (ceiling (- height floating-cliff) cliff-step)
                       do
                          (let* (;(cliff-size (min height (+ floating-cliff (* i cliff-step))))
                                 ;; (cliff-size (min height 100d0))
                                 ;; (flotation (/ (- height cliff-size) (* floating-point height)))
                                 (mps 3)
                                 (output-dir (format nil "./output-~f-~f/" height flotation)))
                            (format t "Problem ~f ~f~%" height flotation)
                            (defparameter *length-scaler* 1d0)
                            ;; (defparameter *length-scale* 20d0)
                            (defparameter *early-exit-flag* nil)
                            (setup :refine 0.125d0
                                   :multigrid-refines 0
                                   :friction 0.5d0
                                   :ice-height height
                                   :mps mps
                                   :hydro-static nil
                                   :cryo-static nil
                                   :elastic-static t
                                   :aspect 6d0
                                   :slope 0d0
                                   :bench-length (* 1d0 height)
                                   :floatation-ratio flotation
                                   :use-penalty t
                                   :stick-base nil)
                            (cl-mpm::domain-sort-mps *sim*)

                            ;; (cl-mpm/output:add-mp-output *sim* :SCALAR "energy-release" (lambda (mp)
                            ;;                                                               (compute-energy-release *sim* mp)))
                            (when (typep *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-octree)
                              (setf (cl-mpm/dynamic-relaxation::sim-intra-mesh-aggregation *sim*) t)
                              (setf (cl-mpm/dynamic-relaxation::sim-octree-refinement-criteria *sim*)
                                    (lambda (sim mesh c)
                                      (multiple-value-bind (damage damage-y length)
                                          (cl-mpm/dynamic-relaxation::damage-refinement-criteria sim mesh c)
                                        (> damage 0d0)))))
                            (plot-domain)
                            (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0)
                            (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
                            (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
                                  (cl-mpm::sim-ghost-factor *sim*) nil
                                  (cl-mpm::sim-velocity-algorithm *sim*) :TBLEND
                                  ;; (cl-mpm::sim-velocity-algorithm *sim*) :QUASI-STATIC
                                  )
                            (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
                            (time
                             (let* ((dt-0 1d3)
                                    (dt-start 1d3)
                                    (max-time 1d8)
                                    (adaptive-constant 2)
                                    (max-adaptive (round (log (/ dt-start dt-0) adaptive-constant)))
                                    (res t))
                               (handler-case
                                   (setf res (cl-mpm/dynamic-relaxation::run-quasi-time
                                              *sim*
                                              :output-dir output-dir
                                              :dt (* dt-0 (expt adaptive-constant max-adaptive))
                                              :total-time max-time
                                              ;; :steps 1000
                                              :dt-scale 0.9d0
                                              :conv-criteria 1d-3
                                              :substeps (max 5 (* 1 (round (/ height (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*))) 1)))
                                              :sub-conv-steps 50
                                              :min-adaptive-steps -20
                                              :max-adaptive-steps max-adaptive
                                              :adaption-constant adaptive-constant
                                              :max-damage-inc 0.90d0
                                              :max-plastic-inc nil
                                              :max-deformation-gradient 4d0
                                              :save-vtk-dr t
                                              :save-vtk-loadstep t
                                              :enable-damage t
                                              :enable-plastic t
                                              ;; :elastic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
                                              ;; :elastic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-dr-multigrid
                                              :plotter (lambda (sim) (plot-domain))
                                              :post-conv-step (lambda (sim)
                                                                (cl-mpm:iterate-over-mps
                                                                 (cl-mpm:sim-mps *sim*)
                                                                 (lambda (mp)
                                                                   (cl-mpm/fastmaths:fast-zero (cl-mpm/particle::mp-displacement mp))))
                                                                (plot-domain)
                                                                )))
                                   (early-exit-condition (c)
                                     (setf res (stable c))
                                     (when res
                                       (setf (cl-mpm:sim-time *sim*) max-time))
                                     (format t "Early exit~%")))
                               (cl-mpm/dynamic-relaxation::save-vtks *sim* output-dir 1)
                               ;; (setf (aref *stability* hi fi) (if res 1 0))
                               ;; (unless res
                               ;;   (loop for j from fi below (length floatations)
                               ;;         do (setf (aref *stability* hi j) 0)))
                               (format t "Stability: ~A~%" res)
                               ;; (print-stab)
                               (save-stabilty-data stability-dir *sim* res height flotation 0d0)
                               (unless res
                                 (loop-finish))))
                            )))))))

(defun stability-ssr-test ()
  (cl-mpm/utils:set-workers 8)
  (let* ((heights (list 500d0))
         (floatations (list
                       0.8d0
                            )))
    (defparameter *stability* (make-array (list (length heights) (length floatations)) :initial-element nil
                                                                                       :element-type t))
    (let ((stability-dir (merge-pathnames (format nil "./analysis_scripts/ice/ice-cliff-stability/data-cliff-stability/"))))
      (ensure-directories-exist stability-dir)
      (defparameter *heights* heights)
      (defparameter *floatations* floatations)
      (loop for hi from 0
            for height in heights
            do
               (let ((res t))
                 (loop for fi from 0
                       for flotation in floatations
                       do
                          (let* ((mps 3)
                                 (output-dir (format nil "./output-~f-~f/" height flotation)))
                            (format t "Problem ~f ~f~%" height flotation)
                            (defparameter *length-scaler* 2d0)
                            (defparameter *length-scale* 20d0)
                            (setup :refine 0.5d0
                                   :multigrid-refines 0
                                   :friction 0d0
                                   :ice-height height
                                   :mps mps
                                   :hydro-static nil
                                   :cryo-static nil
                                   :aspect 6d0
                                   :slope 0d0
                                   ;; :bench-length (* 1d0 height)
                                   :floatation-ratio flotation
                                   :use-penalty nil
                                   :stick-base nil)
                            (cl-mpm::domain-sort-mps *sim*)
                            (plot-domain)
                            (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0)
                            (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
                            (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
                                  (cl-mpm::sim-ghost-factor *sim*) nil)
                            (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep *sim*) 0d0)
                            (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
                            (let ((tensile-max 1d6)
                                  (tensile-min 0.1d6))
                              (let ((res nil))
                                (handler-case
                                    (setf res
                                          (cl-mpm/dynamic-relaxation::run-adaptive-load-control
                                           *sim*
                                           :output-dir output-dir
                                           :dt-scale 0.9d0
                                           :criteria 1d-3
                                           :substeps (* 5 (round height 100))
                                           :sub-conv-steps 200
                                           :loading-function
                                           (lambda (load)
                                             (let ((strength (- tensile-max
                                                                (* (- tensile-max tensile-min) load))))
                                               (format t "Trial strength ~E~%" strength)
                                               (cl-mpm:iterate-over-mps
                                                (cl-mpm:sim-mps *sim*)
                                                (lambda (mp)
                                                  (setf (cl-mpm/particle::mp-initiation-stress mp)
                                                        strength)
                                                  (cl-mpm/particle::update-tcs mp)
                                                  (cl-mpm/damage::compute-damage mp)))))
                                           ;; :min-adaptive-steps -8
                                           ;; :max-adaptive-steps 8
                                           ;; :adaption-constant 4
                                           :max-damage-inc 1.10d0
                                           :max-plastic-inc nil
                                           :max-deformation-gradient 4d0
                                           :save-vtk-dr t
                                           :save-vtk-loadstep t
                                           :enable-damage t
                                           :enable-plastic t
                                           :stagger-damage nil
                                           :load-steps 10
                                           :plotter (lambda (sim) (plot-domain))
                                           :post-conv-step (lambda (sim) (plot-domain))))
                                  (early-exit-condition (c)
                                    (setf res t)
                                    (format t "Early exit~%")))
                                (cl-mpm/dynamic-relaxation::save-vtks *sim* output-dir 1)
                                (setf (aref *stability* hi fi) (if res 1 0))
                                (unless res
                                  (loop for j from fi below (length floatations)
                                        do (setf (aref *stability* hi j) 0)))
                                (format t "Stability:~%")
                                (print-stab)
                                (save-stabilty-data stability-dir *sim* res height flotation 0d0)
                                (unless res
                                  (loop-finish))))
                            )))))))

(declaim (notinline print-stab))
(defun print-stab ()
  (format t "    - " )
  (loop for fi from 0
        for flotation in *floatations*
        do
           (format t "~A - " flotation))
  (format t "~%")
  (loop for hi from 0
        for height in *heights*
        do
           (format t "~A - " height)
           (loop for fi from 0
                 for flotation in *floatations*
                 do
                    (format t "~A - " (aref *stability* hi fi)))
           (format t "~%")
        ))


