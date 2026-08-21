(defpackage :cl-mpm/examples/ice-buoyancy
  (:use :cl
   :cl-mpm/example))
(in-package :cl-mpm/examples/ice-buoyancy)

(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)



(defclass cl-mpm/particle::particle-ice-erodable (cl-mpm/particle::particle-ice-delayed
                                                    cl-mpm/particle::particle-erosion)
  ())


(defmethod cl-mpm/erosion::mp-erosion-enhancment ((mp cl-mpm/particle::particle-ice-erodable))
  (+ 1d0 (* 10 (cl-mpm/particle::mp-damage mp)))
  ;; (+ 1d0 (* 10 (cl-mpm/particle::mp-strain-plastic-vm mp)))
  )

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-ice-brittle) dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-damaged mesh mp dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-dynamic-relaxation mesh mp dt fbar)
  (cl-mpm::update-stress-kirchoff mesh mp dt fbar)
  )

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-ice-brittle) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  ;; (cl-mpm::update-domain-det mesh mp dt)
  ;; (cl-mpm::co-domain-corner-2d mesh mp dt)
  ;; (cl-mpm::update-domain-polar-2d mesh mp dt)
  (cl-mpm::update-domain-polar mesh mp dt)
  ;; (cl-mpm::update-domain-midpoint mesh mp dt)
  ;; (cl-mpm::update-domain-deformation mesh mp dt)
  ;; (cl-mpm::scale-domain-size mesh mp)
  ;; (when (= (cl-mpm/particle::mp-split-depth mp) (- (cl-mpm::sim-max-split-depth *sim*) 0))
  ;;   (cl-mpm::clamp-domains mesh mp))
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
    (declare (double-float E ps-vm angle pressure j))
    (progn
      (let* ((ps-y (the double-float (sqrt (* E (expt ps-vm 2)))))
             (stress (cl-mpm/constitutive:linear-elastic-mat strain de))
             (stress-pressure
               (cl-mpm/fastmaths:fast-.+
                stress
                (cl-mpm/utils:voigt-eye
                 (* j (* (- pressure) 0.3d0))))))
        (setf
         y
         (*
          (+
           (if pd-inc ps-y 0d0)
           ;; (cl-mpm/damage::criterion-max-principal-stress stress-pressure)
           ;; (cl-mpm/damage::criterion-j2 stress)
           ;; (cl-mpm/damage::criterion-j2 stress-pressure)
           ;; (cl-mpm/damage::tensile-energy-norm strain e de)
           (cl-mpm/damage::criterion-mohr-coloumb-rankine-stress-tensile stress-pressure angle)
           ;; (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile stress-pressure angle)
           ;; (cl-mpm/damage::criterion-mohr-coloumb-rankine-stress-tensile stress angle)
           ;; (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile stress-pressure angle)
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

(defparameter *angle* 50d0)
(defparameter *angle-r* 30d0)
(defparameter *angle-psi* 0d0)
(defparameter *rt* (- 1d0 1d-9))
(defparameter *rc* 0d0)
(defparameter *rs* 1d0)
(defparameter *enable-plastic-damage* nil)
(defparameter *delay-time* 1d6)
(defparameter *delay-exponent* 2d0)
(defparameter *enable-viscosity* nil)
(defparameter *length-scaler* 4d0)
(defparameter *gf* 10000d0)
(defparameter *pd-oversize* 1d-4)
(defparameter *ductility* 2d0)
;; (defparameter *tensile-strength* 0.1185d6)
(defparameter *tensile-strength* 0.1d6)

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
                                               ;; 'cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul
                                               'cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic
                                               ;; 'cl-mpm/dynamic-relaxation::mpm-sim-damage-quasi-static-mpi
                                               ;; 'cl-mpm/dynamic-relaxation::mpm-sim-dr-multigrid
                                               ;; 'cl-mpm/dynamic-relaxation::mpm-sim-octree-damage-quasi-static
                                               :args-list
                                               (list
                                                :enable-fbar nil
                                                :enable-aggregate t
                                                :mass-update-count 1
                                                :split-factor (* 1.2d0 (/ 1d0 mps))
                                                ;; :refinement multigrid-refines
                                                :max-split-depth 3)))
    (setf mesh-resolution (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
    ;; (unless (typep *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-octree)
    ;;         (setf multigrid-refines 0))
    (setf h-fine (* mesh-resolution (expt 2 (- multigrid-refines))))
    (let* ((angle *angle*)
           (init-stress *tensile-strength*)
           (init-c (cl-mpm/damage::mohr-coloumb-tensile-to-coheasion init-stress (* angle (/ pi 180))))
           (gf *gf*)
           (length-scale
             ;; 10d0
             (* h-fine *length-scaler*)
                         )
           ;; (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress E))
           (ductility *ductility*)
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
          :E E
          :nu 0.24d0

          :residual-strength rt
          :kt-res-ratio rt
          :kc-res-ratio rc
          :residual-strength (- 1d0 1d-9)
          :initiation-stress init-stress
          :friction-angle (cl-mpm/utils:deg-to-rad angle)
          :residual-friction (cl-mpm/utils:deg-to-rad *angle-r*)

          :psi (cl-mpm/utils::deg-to-rad *angle-psi*)
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
          :material-damping 0d-2
          :index 0))
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
             ))))

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
         :index 0))

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
          (epsilon-scale 1d-2)
          (penalty-damping 0d0)
          )
      (defparameter *floor-bc*
        (cl-mpm/penalty::make-bc-penalty-distance-point
         *sim*
         (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
         (cl-mpm/utils:vector-from-list (list
                                         domain-half
                                         offset
                                         0d0))
         (* domain-half 1.1d0)
         (* E epsilon-scale)
         friction
         penalty-damping)))
    ;; (setf (cl-mpm/penalty::bc-penalty-stiffness-scale *floor-bc*) 1d0)

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
    (cl-mpm/output:add-mp-output
     *sim*
     :SCALAR
     "j"
     #'cl-mpm/particle::mp-deformation-jacobian-strain)
    ))




(defun test-buoy ()

  (let ((step 0)
        (output-dir "./output/"))
    (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir (format nil "sim_~5,'0d.vtk" step)) *sim* )
    (cl-mpm/output:save-vtk-nodes (uiop:merge-pathnames* output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) *sim* )
    (cl-mpm:update-sim *sim*)
    (incf step)
    (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir (format nil "sim_~5,'0d.vtk" step)) *sim* )
    (cl-mpm/output:save-vtk-nodes (uiop:merge-pathnames* output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) *sim* )))

(defun get-damage ()
  (lparallel:pmap-reduce
   (lambda (mp)
     (if (typep mp 'cl-mpm/particle::particle-damage)
         (*
          (cl-mpm/particle::mp-mass mp)
          (cl-mpm/particle:mp-damage mp))
         0d0))
   #'+
   (cl-mpm:sim-mps *sim*)))


(defun stop ()
  (setf *run-sim* nil)
  (setf cl-mpm/dynamic-relaxation::*run-convergance* nil))


(defun damage-refinement-criteria (sim mesh c)
  ;; (let ((damage 0d0)
  ;;       (damage-ybar 0d0)
  ;;       )
  ;;   (cl-mpm/damage::iterate-over-point-neighbour-mps
  ;;    (aref (cl-mpm::sim-mesh-list sim) 0)
  ;;    (cl-mpm/mesh::cell-centroid c)
  ;;    ;; (* 2 *length-scale*)
  ;;    (* 1d0 (cl-mpm/mesh::cell-h c))
  ;;    (lambda (mesh mp dist)
  ;;      (declare (ignore mesh dist))
  ;;      (with-accessors ((d-ybar cl-mpm/particle::mp-damage-ybar)
  ;;                       (d cl-mpm/particle::mp-damage)
  ;;                       (initiation-stress cl-mpm/particle::mp-initiation-stress))
  ;;          mp
  ;;        (declare (double-float damage-ybar initiation-stress damage))
  ;;        (setf damage-ybar (max (* ;; (- 1d0 damage)
  ;;                                  (/ d-ybar initiation-stress)) damage-ybar))
  ;;        (setf damage (max d damage))
  ;;        )))
  ;;   (case (cl-mpm/dynamic-relaxation::cell-mesh-index c)
  ;;     (0  (or (> damage-ybar 2d0)
  ;;             (> damage 0d0)))
  ;;     (1  (or (> damage 0.1d0)))
  ;;     (2  (> damage 0.2d0))
  ;;     (3  (> damage 0.3d0))
  ;;     ;; (2  (> damage 0.85d0))
  ;;     ;; (3  (> damage 0.95d0))
  ;;     (t nil))
  ;;   )
  (multiple-value-bind (damage damage-ybar) (cl-mpm/dynamic-relaxation::damage-refinement-criteria sim mesh c)
    (> damage 0d0)
    ;; (> damage-ybar (* (cl-mpm/dynamic-relaxation::cell-mesh-index c) 1d0))
    )
  )

(defmethod cl-mpm/dynamic-relaxation::damage-increment-criteria ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-ul))
  ;; (cl-mpm/dynamic-relaxation::compute-max-damage-energy-crit sim)
  (cl-mpm/dynamic-relaxation::damage-increment-criteria-mp sim)
  ;; (damage-increment-criteria-mesh sim)
  )

(defun calving-test ()
  (cl-mpm/utils::set-workers 16)
  (let* ((mps 3)
         (dt 1d3)
         (total-time 1d8)
         (H 900d0)
         (ice-aspect 2d0)
         (density 918d0)
         (explicit-dt-scale 0.45d0)
         (water-damping 1d0)
         (floatation-ratio 0.75d0)
         (output-dir "./output/"))
    (defparameter *length-scaler* 2d0)
    (setup
     :refine 0.25
     :multigrid-refines 0
     :friction 0.5d0
     :bench-length (* 0d0 H)
     :bench-extra-cut (* 0d0 (* H 0.25d0))
     :ice-height H
     :mps mps
     :hydro-static nil
     :cryo-static nil
     :melange nil
     :aspect ice-aspect
     :slope 0d0
     :floatation-ratio floatation-ratio
     :use-penalty t
     :extra-offset 0
     :stick-base nil)

    (push (list :SCALAR "water-pressure" #'cl-mpm/particle::mp-pressure) (cl-mpm::sim-output-list *sim*))

    (push
     (list
      :SCALAR
      "damage-delta"
      (lambda (mp)
        (if (typep mp 'cl-mpm/particle::particle-damage)
            (with-accessors ((damage cl-mpm/particle::mp-damage)
                             (damage-prev cl-mpm/particle::mp-damage-prev-trial))
                mp
              (expt (- damage damage-prev) 2))
            0d0)))
     (cl-mpm::sim-output-list *sim*))

    (cl-mpm/output:add-mp-output
     *sim*
     :SCALAR
     "effective-pressure"
     (lambda (mp)
       (if (typep mp 'cl-mpm/particle::particle-ice-brittle)
           (-
            (* 1/3 (cl-mpm/utils::trace-voigt (cl-mpm/particle::mp-undamaged-stress mp)))
            (* 1d0 (cl-mpm/particle::mp-pressure mp)))
           0d0)))
    (cl-mpm/output:add-mp-output
     *sim*
     :SCALAR
     "j"
     #'cl-mpm/particle::mp-deformation-jacobian-strain)

    ;; (cl-mpm/dynamic-relaxation::refine-mesh *sim*)
    (cl-mpm::domain-sort-mps *sim*)
    ;; (setf (cl-mpm/buoyancy::bc-enable *water-bc*) nil)

    (when (typep *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-octree)
      (setf (cl-mpm/dynamic-relaxation::sim-intra-mesh-aggregation *sim*) t)
      (setf (cl-mpm/dynamic-relaxation::sim-octree-refinement-criteria *sim*)
            (lambda (sim mesh c)
              (or
               (damage-refinement-criteria sim mesh c)))))


    (plot-domain)

    (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 1d0)
    (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
          (cl-mpm::sim-ghost-factor *sim*) nil)

    (setf lparallel:*debug-tasks-p* nil)
    (setf (cl-mpm::sim-allow-mp-damage-removal *sim*) nil)
    (setf (cl-mpm::sim-mp-damage-removal-instant *sim*) nil)

    (setf (cl-mpm/damage::sim-enable-stress-based-length *sim*) nil)
    (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
    ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
    ;; (setf (cl-mpm/damage::sim-enable-stress-based-length *sim*) nil)

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

    (setf (cl-mpm:sim-enable-damage *sim*) nil)

    (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)

    ;; (break)
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (let ((step 0))
      (cl-mpm/dynamic-relaxation::run-multi-stage
       *sim*
       :output-dir output-dir
       :dt dt
       :conv-dt-scale 0.9d0
       :dt-scale 0.9d0
       :damping-factor (sqrt 2d0)
       :conv-criteria 1d-3
       :conv-load-steps 1
       ;; :min-adaptive-steps -4
       ;; :max-adaptive-steps 10
       :min-adaptive-steps -14
       :max-adaptive-steps 14
       :adaption-constant 2
       :max-damage-inc 0.9d0
       :max-plastic-inc nil
       ;; :min-damage-inc 0.005d0
       :elastic-dt-margin 1d3
       :substeps 20
       :sub-conv-steps 100
       :total-time total-time
       :save-vtk-loadstep t
       :save-vtk-dr t
       :enable-plastic t
       :enable-damage t
       :plotter (lambda (sim)
                  (plot-domain)
                  (vgplot:title (format nil "Step ~D - Time ~F - ~A"
                                        step
                                        (cl-mpm::sim-time sim)
                                        (if (equal (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
                                            "Quasi-Static"
                                            "Dynamic")))
                  (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step)) :terminal "png size 1920,1080")
                  (incf step))
       ;; :explicit-conv-criteria 1d-2
       :explicit-mass-scaling t
       :explicit-damping-factor 1d-4
       :explicit-dt-scale 0.45d0
       :explicit-dynamic-solver 'cl-mpm/damage::mpm-sim-agg-damage
       ;; :explicit-dynamic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-octree-damage-usf
       ;; :elastic-dt-margin 1d4
       ;; :explicit-mass-scaling nil
       ;; :explicit-dt-scale 10d0
       ;; :explicit-damping-factor 0d-3
       ;; ;; :explicit-dynamic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-octree-implicit-dynamic
       ;; :explicit-dynamic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic
       ;; :elastic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static
       ;; :initial-quasi-static t
       :post-conv-step
       (lambda (sim)
         ;; (cl-mpm::domain-sort-mps *sim*)
         (setf (cl-mpm/buoyancy::bc-enable *bc-erode*) t))
       :setup-quasi-static
       (lambda (sim)
         (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
         (setf (cl-mpm/dynamic-relaxation::sim-true-damping *sim*) (* 1d-3 (cl-mpm/setup::estimate-critical-damping *sim*)))
         (setf
          (cl-mpm/aggregate::sim-enable-aggregate sim) t
          (cl-mpm::sim-ghost-factor *sim*) nil
          (cl-mpm::sim-velocity-algorithm sim) :TPIC
          ;; (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0
          (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) water-damping
          ))
       :setup-dynamic
       (lambda (sim)
         (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
         (setf
          (cl-mpm/damage::sim-damage-delocal-counter-max sim) 10
          (cl-mpm/aggregate::sim-enable-aggregate sim) t
          (cl-mpm::sim-ghost-factor *sim*) nil
          (cl-mpm::sim-velocity-algorithm sim) :TBLEND
          (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) water-damping))
       ))
    ))

(defun calving-real-test ()
  (cl-mpm/utils:set-workers 16)
  (let* ((mps 3)
         (H 900d0)
         (density 918d0))
    (setf *delay-time* 1d4)
    (defparameter *length-scaler* 2d0)
    (setup :refine 0.25
           ;; :multigrid-refines 1
           :friction 0.2d0
           :bench-length (* 0d0 H)
           :bench-extra-cut 100d0
           :ice-height H
           :mps mps
           :hydro-static nil
           :cryo-static t
           :melange nil
           :aspect 4d0
           :slope 0d0
           :floatation-ratio 0.75d0
           :use-penalty t
           ;; :stick-base t
           )
    ;; (change-class *sim* 'cl-mpm/damage::mpm-sim-agg-damage)
    (change-class *sim* 'cl-mpm/damage::mpm-sim-agg-damage)
    ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic)
       ;; ;; :explicit-dynamic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-octree-implicit-dynamic
    ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-octree-damage-usf)
    ;; (setf (cl-mpm/dynamic-relaxation::sim-intra-mesh-aggregation *sim*) t)
    ;; (change-class *sim* 'cl-mpm/damage::mpm-sim-agg-damage)
    (plot-domain)
    (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 1d0)
    (setf
     (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
     (cl-mpm::sim-ghost-factor *sim*) nil)
    ;; (let ((d (- 1d0 1d-4)))
    ;;   (cl-mpm:iterate-over-mps
    ;;    (cl-mpm:sim-mps *sim*)
    ;;    (lambda (mp)
    ;;      (let ((k (cl-mpm/damage::find-k-damage-mp mp d)))
    ;;        (setf (cl-mpm/particle::mp-history-stress mp)
    ;;              k
    ;;              (cl-mpm/particle::mp-history-stress-n mp) k))
    ;;      (cl-mpm/damage::compute-damage mp))))
    (setf (cl-mpm::sim-allow-mp-damage-removal *sim*) nil)
    (setf (cl-mpm::sim-mp-damage-removal-instant *sim*) nil)
    (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
    (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) nil)
    (setf (cl-mpm:sim-enable-damage *sim*) nil)
    (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (setf (cl-mpm::sim-velocity-algorithm *sim*) :TBLEND)
    ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic)
    (let ((step 0))
      (cl-mpm/dynamic-relaxation::run-time
       *sim*
       :output-dir "./output/"
       :dt 1d2
       :total-time 1d5
       ;; :dt-scale 1000d0
       :dt-scale 0.9d0
       :mass-scale 1d4
       :damping 1d-3
       :enable-plastic t
       :enable-damage t
       :conv-criteria 1d-3
       :save-vtk-loadstep t
       :initial-quasi-static t
       ;; :elastic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-octree-damage-quasi-static
       :elastic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul
       :post-conv-step
       (lambda (sim)
         (setf
          (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
          (cl-mpm::sim-ghost-factor *sim*) nil
          ;; (cl-mpm/aggregate::sim-enable-aggregate *sim*) nil
          ;; (cl-mpm::sim-ghost-factor *sim*) (* density 1d-4)
          ;; (cl-mpm::sim-ghost-factor *sim*) (* 1d9 1d-9)
          ))

       :plotter
       (lambda (sim)
         (plot-domain)
         (vgplot:title (format nil "Step ~D - Time ~F - ~A"
                               step
                               (cl-mpm::sim-time sim)
                               (if (equal (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
                                   "Quasi-Static"
                                   "Dynamic")))
         (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step)) :terminal "png size 1920,1080")
         (incf step))
       ))))

(defmethod cl-mpm/dynamic-relaxation::convergence-check ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic))
  (if (> (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim) 0d0)
      (progn
        (format t "Check inertia~%")
        ;; (let* ((elastic-dt (cl-mpm/setup::estimate-elastic-dt sim))
        ;;        (current-dt (cl-mpm/dynamic-relaxation::sim-dt-loadstep sim))
        ;;        ;; (current-intertia (cl-mpm/dynamic-relaxation::true-intertial-criteria sim current-dt))
        ;;        (current-intertia (cl-mpm/dynamic-relaxation::true-intertial-criteria sim elastic-dt))
        ;;        (thresh 1d-1))
        ;;   (cl-mpm:sim-format sim t "Current inertia ~E - ~E - ~E~%" current-intertia (/ current-intertia thresh) (/ current-dt elastic-dt))
        ;;   (when (> (/ current-intertia thresh) (/ current-dt elastic-dt))
        ;;     (when (> (/ current-dt elastic-dt) 1d0)
        ;;       (cl-mpm:sim-format sim t "Inertia criteria exceeded~%")
        ;;       (error (make-instance 'cl-mpm/dynamic-relaxation::error-inertia-criteria
        ;;                             :text "Ratio of true inertia to elastic inertia exceeded"
        ;;                             :inertia-norm current-intertia))))
        ;;   )
        t)
      t))

(defmethod cl-mpm/dynamic-relaxation::damage-increment-criteria ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic))
  (cl-mpm/dynamic-relaxation::damage-increment-criteria-mp sim)
  ;; (cl-mpm/dynamic-relaxation::compute-max-damage-energy-crit sim)
  )

(defun calving-qs-test ()
  (cl-mpm/utils:set-workers 16)
  (let* ((mps 3)
         (H 900d0)
         (density 918d0))
    (setf *delay-time* 1d5)
    (defparameter *length-scaler* 2d0)
    (setup :refine 0.125
           ;; :multigrid-refines 1
           :friction 0.2d0
           :bench-length (* 0d0 H)
           :bench-extra-cut 00d0
           :ice-height H
           :mps mps
           :hydro-static nil
           :cryo-static t
           :melange nil
           :aspect 2d0
           :slope 0.0d0
           :floatation-ratio 0.5d0
           :use-penalty t
           ;; :stick-base t
           )
    ;; (cl-mpm/output:add-mp-output
    ;;  *sim*
    ;;  :VECTOR
    ;;  "bf"
    ;;  #'cl-mpm/particle::mp-body-force)
    (push (list :SCALAR "water-pressure" #'cl-mpm/particle::mp-pressure) (cl-mpm::sim-output-list *sim*))
    (push (list :SCALAR "water-pressure-damage"
                (lambda (mp) (*
                              (cl-mpm/particle::mp-damage mp)
                              (cl-mpm/particle::mp-pressure mp)))) (cl-mpm::sim-output-list *sim*))
    (cl-mpm/output:add-node-output
     *sim*
     :SCALAR
     "true-mass"
     #'cl-mpm/mesh::node-true-mass)
    (cl-mpm/output:add-node-output
     *sim*
     :VECTOR
     "inertia"
     #'cl-mpm/mesh::node-inertia-force)
    ;; (change-class *sim* 'cl-mpm/damage::mpm-sim-agg-damage)
    (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic)
    ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul)
    (plot-domain)
    (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0)
    (setf (cl-mpm::sim-allow-mp-damage-removal *sim*) nil)
    (setf (cl-mpm::sim-mp-damage-removal-instant *sim*) nil)
    (setf (cl-mpm/damage::sim-enable-ekl *sim*) nil)
    (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
    (setf (cl-mpm:sim-enable-damage *sim*) nil)
    (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (setf (cl-mpm::sim-velocity-algorithm *sim*) :TBLEND)
    ;; (setf (cl-mpm/dynamic-relaxation::sim-enable-dynamics *sim*) nil)
    (let ((step 0))
      (cl-mpm/dynamic-relaxation::run-quasi-time
       *sim*
       :output-dir "./output/"
       :dt 1d3
       :total-time 1d6
       :dt-scale 0.9d0
       :enable-plastic t
       :enable-damage t
       :substeps 50
       :sub-conv-steps 100
       ;; :min-adaptive-steps 0
       ;; :max-adaptive-steps 0
       :min-adaptive-steps -12
       :max-adaptive-steps 12
       :adaption-constant 2
       :adaption-easy-steps 8
       :max-damage-inc 0.9d0
       :max-plastic-inc nil
       :max-deformation-gradient 50d0

       :conv-criteria 1d-3
       :save-vtk-loadstep t
       :elastic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static
       :initial-quasi-static t
       :post-conv-step
       (lambda (sim)
         ;; (setf (cl-mpm/dynamic-relaxation::sim-true-damping *sim*) (* 1d-4 (cl-mpm/setup::estimate-critical-damping *sim*)))
         ;; (setf (cl-mpm::sim-mass-scale *sim*) 1d0)
         (setf
          (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
          (cl-mpm::sim-ghost-factor *sim*) nil))
       :plotter
       (lambda (sim)
         (plot-domain))
       :post-load-step
       (lambda (sim)
         (plot-domain)
         (vgplot:title (format nil "Step ~D - Time ~F - ~A"
                               step
                               (cl-mpm::sim-time sim)
                               (if (equal (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
                                   "Quasi-Static"
                                   "Dynamic")))
         (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step)) :terminal "png size 1920,1080")
         (incf step))))))





(defun save-test-vtks (&key (output-dir "./output/"))
  (cl-mpm::finalise-loadstep *sim*)
  (cl-mpm/output:save-vtk (merge-pathnames "test.vtk" output-dir) *sim*)
  (cl-mpm/output:save-vtk-nodes (merge-pathnames "test_nodes.vtk" output-dir) *sim*)
  (cl-mpm/output:save-vtk-cells (merge-pathnames "test_cells.vtk" output-dir) *sim*)
  )


(defun est-shear-from-angle (angle angle-r rc)
  (let* ((angle-plastic (* angle (/ pi 180)))
         (angle-plastic-residual (* angle-r (/ pi 180))))
    (- 1d0
       (* (- 1d0 rc)
          (/ (tan angle-plastic-residual)
             (tan angle-plastic))))))

(defun est-angle (angle rs rc)
  (let* ((ratio (if (< rc 1d0) (/ (- 1d0 rs) (- 1d0 rc)) 0d0))
         (angle-plastic (* angle (/ pi 180)))
         (angle-plastic-damaged (atan (* ratio (tan angle-plastic))))
         )
    (format t "Plastic virgin angle: ~F~%"
            (* (/ 180 pi) angle-plastic))
    (format t "Plastic residual angle: ~F~%"
            (* (/ 180 pi) angle-plastic-damaged))))





(defun test-dist ()
  (let* ((mp-a (find-mp *sim* (cl-mpm/utils:vector-from-list (list 270d0 100d0 0d0))))
         (mp-b (find-mp *sim* (cl-mpm/utils:vector-from-list (list 280d0 100d0 0d0))))
         (mesh (cl-mpm:sim-mesh *sim*)))
    (pprint (cl-mpm/damage::diff-squared mp-a mp-b))
    (pprint (cl-mpm/damage::diff-damaged mesh mp-a mp-b))))

(defun plot-nonlocal-inter ()
  (let* ((find-pos (cl-mpm/utils:vector-from-list (list 207d0 470d0 0d0)))
         (mp (cl-mpm/setup::find-mp *sim* find-pos)))
    (multiple-value-bind (pos weights) (cl-mpm/damage::get-nonlocal-interactions-stress-based *sim* mp)
      (let ((x (loop for p in pos collect (cl-mpm/utils:varef p 0)))
            (y (loop for p in pos collect (cl-mpm/utils:varef p 1))))
        (loop for p in pos
              do (format t "~A ~A~%" (cl-mpm/utils:varef p 0) (cl-mpm/utils:varef p 1)))
        (vgplot:format-plot t "set xrange [~f:~f]" -20d0 20d0)
        (vgplot:format-plot t "set yrange [~f:~f]" -20d0 20d0)
        (vgplot:3d-plot x y weights ";;with points lc palette")
        (vgplot:xlabel "x")
        (vgplot:ylabel "y")
        ))))



(defun elastic-solution ()
  (cl-mpm/utils::set-workers 8)
  (let* ((mps 3)
         (dt 1d3)
         (total-time 1d8)
         (H 400d0)
         (ice-aspect 4d0)
         (density 918d0)
         (explicit-dt-scale 0.45d0)
         (water-damping 5d0)
         (floatation-ratio 0.95d0)
         (output-dir "./output-noslip-isostress/")
         (post-iter-step (lambda (sim)))
         )
    (defparameter *length-scaler* 2d0)
    (setup
     :refine 0.5
     :multigrid-refines 0
     :friction 0.5d0
     :bench-length (* 0d0 H)
     :bench-extra-cut 00d0
     :ice-height H
     :mps mps
     :hydro-static nil
     :cryo-static t
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
         :substeps 5
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
    ))


;; (let* ((strain (cl-mpm/utils:voigt-from-list (list -2d0 -2d0 -2d0 0d0 0d0 5d0)))
;;        (de (cl-mpm/constitutive::linear-elastic-matrix 1d0 0d0))
;;        (stress (cl-mpm/constitutive:linear-elastic-mat strain de))
;;        (j 1d0)
;;        (pressure (- 1d0))
;;        (stress-pressure
;;          (cl-mpm/fastmaths:fast-.+
;;           stress
;;           (cl-mpm/utils:voigt-eye (*
;;                                    j
;;                                    (/ (- pressure) 1)
;;                                    ))))
;;        (angle (cl-mpm/utils:deg-to-rad 40d0))
;;        )
;;   (pprint (* 1/3 (cl-mpm/utils:trace-voigt stress)))
;;   (pprint (* 1/3 (cl-mpm/utils:trace-voigt stress-pressure)))
;;   (pprint (cl-mpm/damage::criterion-mohr-coloumb-rankine-stress-tensile stress angle))
;;   (pprint (cl-mpm/damage::criterion-mohr-coloumb-rankine-stress-tensile stress-pressure angle)))
