(defpackage :cl-mpm/examples/ice/cliff-stress
  (:use :cl
   :cl-mpm/example))
(in-package :cl-mpm/examples/ice/cliff-stress)

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
      (let* ((ps-y ;(sqrt (* E (expt ps-vm 2)))
               (sqrt (* E (expt ps-vm 2)))
               )
             ;; (stress (cl-mpm/constitutive:linear-elastic-mat trial-strain de))
             (stress
               ;; undamaged-stress
               ;; (cl-mpm/fastmaths:fast-scale!
                (cl-mpm/constitutive:linear-elastic-mat strain de)
               ;;  undamaged-stress
               ;;  1d0
               ;;  ;; (/ 1d0 (magicl:det def))
               ;;  )
               )
             (stress-pressure
               (cl-mpm/fastmaths:fast-.+
                stress
                (cl-mpm/utils:voigt-eye (*
                                         ;; 1/3
                                         j
                                         ;; (- 1d0 damage)
                                         ;; damage
                                         ;; (- pressure)
                                         (/ (- pressure) 3)
                                         ))))
             )
        (setf y
              (*
               (+
                (if pd-inc ps-y 0d0)
                ;; (cl-mpm/damage::criterion-max-principal-stress stress-pressure)
                ;; (cl-mpm/damage::criterion-j2 stress)
                ;; (cl-mpm/damage::criterion-j2 stress-pressure)
                ;; (cl-mpm/damage::tensile-energy-norm strain e de)
                (cl-mpm/damage::criterion-mohr-coloumb-rankine-stress-tensile stress-pressure angle)
                ;; (cl-mpm/damage::criterion-mohr-coloumb-rankine-stress-tensile stress angle-rad)
                ;; (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile stress-pressure angle-rad)
                ;; (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile stress angle-rad)
                ;; (cl-mpm/damage::drucker-prager-criterion stress angle-rad)
                ;; (cl-mpm/damage::drucker-prager-criterion stress-pressure angle-rad)
                )))))))


(declaim (notinline plot-domain))
(defun plot-domain (&key (trial t))
  (when *sim*
    (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
           (h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
           (ms-x (first ms))
           (ms-y (second ms)))
      (vgplot:format-plot t "set object 1 rect from 0,0 to ~f,~f fc rgb 'blue' fs transparent solid 0.5 noborder behind" ms-x *water-height*)
      (vgplot:format-plot t "set object 2 rect from 0,0 to ~f,~f fc rgb 'black' fs transparent solid 1 noborder behind" ms-x *offset*))
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :trial trial
     ;; :colour-func (lambda (mp) (sqrt (cl-mpm/constitutive::voigt-j2 (cl-mpm/utils:deviatoric-voigt (cl-mpm/particle:mp-stress mp)))))
     ;; :colour-func #'cl-mpm/particle::mp-strain-plastic-vm
     ;; :colour-func (lambda (mp) (/ (cl-mpm/particle:mp-mass mp)
     ;;                              (cl-mpm/particle:mp-volume mp)))
     :colour-func #'cl-mpm/particle::mp-damage
     ;; :colour-func (lambda (mp) (cl-mpm/utils::varef  (cl-mpm/particle::mp-stress mp ) 1))
     ))
  )

(defparameter *angle* 40d0)
(defparameter *angle-r* 10d0)
(defparameter *angle-psi* 0d0)
(defparameter *rt* 1d0)
(defparameter *rc* 0d0)
(defparameter *rs* 1d0)
(defparameter *enable-plastic-damage* nil)
(defparameter *delay-time* 1d6)
(defparameter *delay-exponent* 4d0)
(defparameter *enable-viscosity* nil)
(defparameter *length-scaler* 1d0)
(defparameter *gf* 10000d0)
(defparameter *pd-oversize* 1d-4)

(declaim (notinline setup))
(defun setup (&key (refine 1) (mps 2)
                (pressure-condition t)
                (cryo-static t)
                (hydro-static nil)
                (melange nil)
                (friction 0d0)
                (ice-height 400d0)
                (bench-length 0d0)
                (bench-extra-cut 0d0)
                (aspect 1)
                (floatation-ratio 0.9)
                (slope 0.05d0)
                (multigrid-refines 0)
                (extra-offset 0)
                (use-penalty t)
                (stick-base t))
  (let* ((density 918d0)
         (water-density 1028d0)
         (mesh-resolution (/ 10d0 refine))
         (h-fine mesh-resolution)
         (offset (* mesh-resolution (+ (if use-penalty 2 0) extra-offset)))
         (end-height ice-height)
         (ice-length (* ice-height aspect))
         (start-height (+ ice-height (* slope ice-length)))
         (ice-height end-height)
         (floating-point (* ice-height (/ density water-density)))
         (water-level (* floating-point floatation-ratio))
         (datum (+ water-level offset))
         (domain-size (list (+ ice-length (* 2 ice-height))
                            (+ (* 1d0 offset)
                               (* start-height 2))
                            ;; ice-length
                            )
                      )
         (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
         (block-size (list ice-length (max start-height end-height)
                           ;; ice-length
                           ))
         (E 1d9))
    (defparameter *water-height* datum)
    (defparameter *offset* offset)
    (defparameter *ice-length* ice-length)
    (setf *sim* (cl-mpm/setup::make-simple-sim mesh-resolution element-count
                                               :sim-type
                                               'cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul
                                               ;; 'cl-mpm/dynamic-relaxation::mpm-sim-dr-multigrid
                                               ;; 'cl-mpm/dynamic-relaxation::mpm-sim-octree-damage-quasi-static
                                               :args-list
                                               (list
                                                :enable-fbar nil
                                                :enable-aggregate t
                                                :mass-update-count 1
                                                :split-factor (* 1.2d0 (/ 1d0 mps))
                                                ;; :refinement multigrid-refines
                                                :max-split-depth 6)))
    (setf mesh-resolution (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
    ;; (unless (typep *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-octree)
    ;;         (setf multigrid-refines 0))
    (setf h-fine (* mesh-resolution (expt 2 (- multigrid-refines))))
    (let* ((angle *angle*)
           (init-stress (* 0.1185d6 1d0))
           (init-c (cl-mpm/damage::mohr-coloumb-tensile-to-coheasion init-stress (* angle (/ pi 180))))
           (gf *gf*)
           (length-scale
             ;; 10d0
             (* h-fine *length-scaler*)
                         )
           (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress E))
           (oversize (cl-mpm/damage::compute-oversize-factor (- 1d0 *pd-oversize*) ductility)))
      (defparameter *length-scale* length-scale)

      (format t "Ice length ~F~%" ice-length)
      (format t "Water height ~F~%" water-level)
      (format t "True Water height ~F~%" (- datum offset))
      (format t "Cliff height ~F~%" (- (+ offset ice-height) datum))
      (format t "Mesh size ~F~%" mesh-resolution)
      (format t "Fine mesh size ~F~%" h-fine)
      (format t "Estimated oversize ~F~%" oversize)
      (format t "Estimated lc ~E~%" length-scale)
      (format t "Estimated ductility ~E~%" ductility)
      (format t "Init stress ~E~%" init-stress)
      (format t "Init c ~E~%" init-c)
      (let* ((rt *rt*)
             (rc *rc*)
             (rs (est-shear-from-angle angle *angle-r* rc))
             )
        (format t "Strengths: Tension ~E - Compression ~E - shear ~E~%" rt rc rs)
        (cl-mpm:add-mps
         *sim*
         (cl-mpm/setup:make-block-mps
          (list 0 offset 0)
          block-size
          (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
          density
          'cl-mpm/particle::particle-ice-erodable
          ;; 'cl-mpm/particle::particle-ice-delayed
          :E 1d9
          :nu 0.24d0

          :kt-res-ratio rt
          :kc-res-ratio rc
          :residual-strength 1d0
          :initiation-stress init-stress
          :friction-angle (cl-mpm/utils:deg-to-rad angle)
          :residual-friction (cl-mpm/utils:deg-to-rad *angle-r*)

          :psi (cl-mpm/utils::deg-to-rad 0d0)
          :oversize (- 1d0 *pd-oversize*)

          :ductility ductility
          :local-length length-scale
          :delay-time *delay-time*
          :delay-exponent *delay-exponent*
          :enable-plasticity t
          :enable-damage t
          :enable-viscosity *enable-viscosity*
          :viscosity 1d13
          :plastic-damage-evolution *enable-plastic-damage*

          :index 0

          ;; 'cl-mpm/particle::particle-finite-viscoelastic-ice
          ;; :E 1d9
          ;; :nu 0.325d0
          ;; :visc-factor 111d6
          ;; :visc-power 3d0
          ;; :enable-viscosity nil


          ;; 'cl-mpm/particle::particle-finite-viscoelastic-ice
          ;; :E 1d9
          ;; :nu 0.325d0
          ;; :visc-factor 11.1d6
          ;; :visc-power 3d0
          ;; :enable-viscosity nil

          ;; :gravity -9.8d0
          ))
        (est-angle angle rs rc)
        (when melange
          (let* ((melange-depth (* ice-height 0.1d0))
                 (melange-length (- (first domain-size) (first block-size)))
                 (block-size (list melange-length melange-depth))
                 (offset (- datum (* melange-depth (/ density water-density))))
                 )
            (defparameter *bc-melange*
              (cl-mpm/buoyancy:make-bc-pressure
               *sim*
               -1d3
               0d0
               :clip-func
               (lambda (pos)
                 (and
                  (>= (cl-mpm/utils::varef pos 1)
                      (- datum (* melange-depth (/ density water-density))))
                  (<= (cl-mpm/utils::varef pos 1) datum)))))
            (cl-mpm::add-bcs-force-list
             *sim*
             *bc-melange*
             )
            ;; (cl-mpm:add-mps
            ;;  *sim*
            ;;  (cl-mpm/setup:make-block-mps
            ;;   (list ice-length offset)
            ;;   block-size
            ;;   (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
            ;;   density
            ;;   'cl-mpm/particle::particle-ice-erodable
            ;;   ;; 'cl-mpm/particle::particle-ice-delayed
            ;;   :E 1d9
            ;;   :nu 0.24d0

            ;;   :kt-res-ratio rt
            ;;   :kc-res-ratio rc
            ;;   :g-res-ratio rs

            ;;   :initiation-stress init-stress;18d3
            ;;   :friction-angle angle
            ;;   :psi (* 0d0 (/ pi 180))
            ;;   :phi (* angle (/ pi 180))
            ;;   :c (* init-c oversize)
            ;;   :softening 0d0
            ;;   :ductility ductility
            ;;   :local-length length-scale
            ;;   :delay-time 1d3
            ;;   :delay-exponent 1d0
            ;;   :enable-plasticity t
            ;;   :enable-damage t
            ;;   :enable-viscosity nil
            ;;   :plastic-damage-evolution *enable-plastic-damage*
            ;;   :index 1
            ;;   ))
            ;; (cl-mpm::iterate-over-mps
            ;;  (cl-mpm::sim-mps *sim*)
            ;;  (lambda (mp)
            ;;    (when (= (cl-mpm/particle::mp-index mp) 1)
            ;;      (cl-mpm/damage::set-mp-damage mp 0.99d0))))
            )
          ))

      ;; (cl-mpm/setup::remove-sdf *sim*
      ;;                           (lambda (p)
      ;;                             (cl-mpm/setup::plane-point-point-sdf
      ;;                              p
      ;;                              (cl-mpm/utils:vector-from-list (list 0d0 (+ offset start-height) 0d0))
      ;;                              (cl-mpm/utils:vector-from-list (list ice-length (+ offset end-height) 0d0))))
      ;;                           :refine 1
      ;;                           )
      
      (unless (= start-height end-height)
        (cl-mpm/setup::remove-sdf *sim*
                                  (lambda (p)
                                    (cl-mpm/setup::plane-point-point-sdf
                                     p
                                     (cl-mpm/utils:vector-from-list (list 0d0 (+ offset start-height) 0d0))
                                     (cl-mpm/utils:vector-from-list (list ice-length (+ offset end-height) 0d0))))
                                  :refine 2))
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
         :index 0
         ))


      (let ((cutout (+ (- ice-height water-level) bench-extra-cut))
            (cutback bench-length))
        (when (> cutback 0d0)
          (cl-mpm/setup:remove-sdf
           *sim*
           (cl-mpm/setup::rectangle-sdf (list (first block-size) (+ offset ice-height ice-height))
                                        (list cutback (+ cutout ice-height)))))))


    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)
    (when (typep *sim* 'cl-mpm/damage::mpm-sim-damage)
      (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t))
    (setf (cl-mpm::sim-allow-mp-split *sim*) t)
    (setf (cl-mpm::sim-ghost-factor *sim*) nil)
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
           :visc-damping 0d0)
          (cl-mpm/buoyancy::make-bc-buoyancy-body
           *sim*
           datum
           water-density
           (lambda (pos) t))))

    (cl-mpm:add-bcs-force-list
     *sim*
     *water-bc*)

    (let ((domain-half (* 0.5d0 (first domain-size)))
          (friction friction)
          (epsilon-scale 1d-2))
      (defparameter *floor-bc*
        (cl-mpm/penalty::make-bc-penalty-distance-point
         *sim*
         (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
         (cl-mpm/utils:vector-from-list (list
                                         domain-half
                                         offset
                                         0d0))
         (* domain-half 1.1d0 100)
         (* E epsilon-scale)
         friction
         0d0)))

    (defparameter *bc-erode*
      (cl-mpm/erosion::make-bc-erode
       *sim*
       :enable nil
       :rate 1d-3
       :scalar-func (lambda (pos)
                      1d0
                      ;; (min 1d0 (exp (* 0.5d0 (- (cl-mpm/utils:varef pos 1) datum))))
                      )
       :clip-func (lambda (pos)
                    (and
                     (>= datum (cl-mpm/utils:varef pos 1))
                     (>= (cl-mpm/utils:varef pos 1) (+ offset mesh-resolution) )))))
    (when use-penalty
      (cl-mpm:add-bcs-force-list
       *sim*
       *floor-bc*
       )
      )
    (unless use-penalty
      (if stick-base
          (cl-mpm/setup:setup-bcs
           *sim*
           :left '(0 nil 0)
           :bottom '(0 0 0)
           ;; :front '(0 0 0)
           ;; :back '(0 0 0)
           )
          (cl-mpm/setup:setup-bcs
           *sim*
           :left '(0 nil 0)
           :bottom '(nil 0 0)
           ;; :front '(0 0 0)
           ;; :back '(0 0 0)
           )
          ))
    ;; (cl-mpm:add-bcs-force-list
    ;;  *sim*
    ;;  *bc-erode*
    ;;  )
    (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
    ))







(defun elastic-solution ()
  (cl-mpm/utils::set-workers 8)
  (dolist (cryo (list nil t))
    (dolist (friction (list 0d0 0.25 0.5d0))
      (let* ((mps 3)
             (dt 1d3)
             (total-time 1d8)
             (H 400d0)
             (ice-aspect 4d0)
             (density 918d0)
             (explicit-dt-scale 0.45d0)
             (water-damping 5d0)
             (floatation-ratio 0.95d0)
             (output-dir (format nil "./output-H_~F_float_~E_cryo_~A_friction_~F/" H floatation-ratio cryo friction)))
        (defparameter *length-scaler* 2d0)
        (setup
         :refine 0.5
         :multigrid-refines 0
         :friction friction
         :bench-length (* 0d0 H)
         :bench-extra-cut 00d0
         :ice-height H
         :mps mps
         :hydro-static nil
         :cryo-static cryo
         :melange nil
         :aspect ice-aspect
         :slope 0d0
         :floatation-ratio floatation-ratio
         :use-penalty nil
         :extra-offset 0
         :stick-base nil)
        ;; (cl-mpm/dynamic-relaxation::refine-mesh *sim*)
        (cl-mpm::domain-sort-mps *sim*)

        (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0)
        (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
              (cl-mpm::sim-ghost-factor *sim*) nil)

        (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
        (setf (cl-mpm:sim-settings *sim*)
              (list :OCEAN-HEIGHT *water-height*
                    :EXPLICIT-DT-SCALE explicit-dt-scale
                    :EKL  (cl-mpm/damage::sim-enable-ekl *sim*)
                    :LENGTH-LOCALISATION  (cl-mpm/damage::sim-enable-length-localisation *sim*)
                    :PLASTIC-DAMAGE-DRIVING *enable-plastic-damage*
                    :PLASTIC-DAMAGE-OVERSIZE *pd-oversize*
                    :DELAY-TIME *delay-time*
                    :DELAY-EXP *delay-exponent*
                    :ANGLE *angle*
                    :ANGLE-R *angle-r*
                    :ANGLE-PSI *angle-psi*
                    :WATER-DAMPING water-damping
                    :R-C *rc*
                    :GF *gf*
                    :LENGTH-SCALER *length-scaler*))
        (uiop:ensure-all-directories-exist (list output-dir))
        (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
        (cl-mpm/dynamic-relaxation::save-conv-preamble *sim* output-dir)
        (cl-mpm/output::save-simulation-parameters
         (merge-pathnames output-dir "settings.json")
         *sim* (list) )
        (let ( (total-iter 0)
               )
          (dotimes (i 2)
            (cl-mpm/dynamic-relaxation::elastic-static-solution
             *sim*
             :crit 1d-6
             :dt-scale 0.9d0
             :substeps 50
             :output-dir output-dir
             :post-iter-step
             (lambda (step energy oobf)
               (cl-mpm/dynamic-relaxation::save-vtks-dr-step *sim* output-dir i step)
               (cl-mpm/dynamic-relaxation::save-conv-step *sim* output-dir total-iter i 0d0 oobf energy)
               (incf total-iter)
               )
             )
            (setf (cl-mpm:sim-enable-damage *sim*) t)
            (cl-mpm/damage::calculate-damage *sim* 1d-9)
            (setf (cl-mpm:sim-enable-damage *sim*) nil)
            (cl-mpm/dynamic-relaxation::save-vtks *sim* output-dir i)))
        ;; (cl-mpm/dynamic-relaxation::save-vtks *sim* output-dir 1)
        ))))
