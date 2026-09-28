(defpackage :cl-mpm/examples/ice/elastic-stress
  (:use :cl
   :cl-mpm/example))
(in-package :cl-mpm/examples/ice/elastic-stress)

(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
;; (setf sb-ext::*block-compile-default* t)
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)


(defparameter *angle* 38d0)
(defparameter *angle-r* 38d0)
(defparameter *angle-psi* 0d0)
(defparameter *rt* (- 1d0 1d-3))
(defparameter *rc* 0d0)
(defparameter *enable-plastic-damage* nil)
(defparameter *delay-time* 1d4)
(defparameter *delay-exponent* 4d0)
(defparameter *enable-viscosity* nil)
(defparameter *length-scaler* 2d0)
(defparameter *gf* 10000d0)
(defparameter *pd-oversize* 1d-6)
(defparameter *ductility* 10d0)
(defparameter *tensile-strength* 0.1d6)
(defparameter *biot-coefficent* 1d0)
(defparameter *alpha* 0.3d0)

;; (defparameter *alpha* 0.4d0)


(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-ice-brittle) dt fbar)
  (cl-mpm::update-stress-kirchoff mesh mp dt fbar))

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-ice-brittle) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  (cl-mpm::update-domain-polar mesh mp dt))

(cl-mpm/utils::with-voigt-pool
    (defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-ice-brittle) dt)
      (with-accessors ((undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                       (y cl-mpm/particle::mp-damage-y-local)
                       (strain cl-mpm/particle::mp-strain)
                       (trial-strain cl-mpm/particle::mp-trial-strain)
                       (plastic-strain cl-mpm/particle::mp-strain-plastic)
                       (ps-vm cl-mpm/particle::mp-strain-plastic-vm)
                       (ps-vm-inc cl-mpm/particle::mp-strain-plastic-vm-inc)
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
          (let* ((ps-y (the double-float
                            (* E (max 0d0 ps-vm-inc))
                            ;; (sqrt (* E (expt ps-vm-inc 2)))
                            ))
                 ;; (undamaged-stress (cl-mpm/constitutive:linear-elastic-mat strain de))
                 ;; (undamaged)
                 (stress-pressure
                   (cl-mpm/fastmaths:fast-.+
                    undamaged-stress
                    ;; stress
                    (cl-mpm/utils:voigt-eye
                     (*
                      (the double-float *alpha*)
                      (the double-float (cl-mpm/particle::mp-biot-coefficent mp))
                      j
                      (- pressure))
                     (grab-new-voigt)
                     )
                    (grab-new-voigt))))
            (setf
             y
             (*
              ;; (- 1d0 damage)
              (+
               (if pd-inc ps-y 0d0)
               ;; (cl-mpm/damage::criterion-max-principal-stress stress-pressure)
               ;; (cl-mpm/damage::criterion-j2 stress)
               ;; (cl-mpm/damage::criterion-j2 stress-pressure)
               ;; (cl-mpm/damage::tensile-energy-norm strain e de)
               ;; (cl-mpm/damage::criterion-mohr-coloumb-rankine-stress-tensile stress-pressure angle)
               ;; (cl-mpm/damage::criterion-mohr-coloumb-rankine-stress-tensile stress angle)
               (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile stress-pressure angle)
               ;; (cl-mpm/damage::criterion-mohr-coloumb-rankine-stress-tensile stress-pressure angle)
               ;; (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile stress-pressure angle)
               ;; (cl-mpm/damage::drucker-prager-criterion stress angle)
               ;; (cl-mpm/damage::drucker-prager-criterion stress-pressure angle)
               ))))))))

(defparameter *water-height* 0d0)
(defparameter *offset* 0d0)
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
     ;; :colour-func #'cl-mpm/particle::mp-damage-ybar
     ;; :colour-func #'cl-mpm/particle::mp-pressure
     ;; :colour-func (lambda (mp) (cl-mpm/utils::varef (cl-mpm/particle::mp-av-damage-gradient mp) 0))

     ;; :colour-func (lambda (mp) (if (> (cl-mpm/particle::mp-damage-ybar mp) 0.1d6) 1d0 0d0))
     )))

(defparameter *bc-melange* nil)
(defparameter *water-bc* nil)
(defparameter *floor-bc* nil)

(declaim (notinline setup))
(defun setup (&key (refine 1) (mps 2)
                (pressure-condition t)
                (cryo-static t)
                (elastic-static nil)
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
                (stick-base t)
                (epsilon-scale 1d-2)
                (E 1d9)
                )
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
         (domain-size (list (+ ice-length (* 4 ice-height))
                            (+ (* 1d0 offset)
                               (* start-height 2))
                            ;; ice-length
                            )
                      )
         (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
         (block-size (list ice-length (max start-height end-height)
                           ;; ice-length
                           )))
    (defparameter *water-height* datum)
    (defparameter *offset* offset)
    (defparameter *ice-length* ice-length)
    (setf *sim* (cl-mpm/setup::make-simple-sim mesh-resolution element-count
                                               :sim-type
                                               'cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul
                                               ;; 'cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic
                                               ;; 'cl-mpm/dynamic-relaxation::mpm-sim-damage-quasi-static-mpi
                                               ;; 'cl-mpm/dynamic-relaxation::mpm-sim-dr-multigrid
                                               ;; 'cl-mpm/dynamic-relaxation::mpm-sim-octree-damage-quasi-static
                                               :args-list
                                               (list
                                                :enable-fbar nil
                                                :enable-aggregate t
                                                :ghost-factor nil
                                                :remove-on-oversplit t
                                                :mass-update-count 1
                                                :split-factor (* 1.2d0 (/ 1d0 mps))
                                                ;; :refinement multigrid-refines
                                                :max-split-depth 3
                                                :enable-split t
                                                ;; :mp-removal-size 0.5d0
                                                ;; :damage-removal t
                                                :damage-removal nil
                                                :damage-removal-crit (- 1d0 *pd-oversize*)
                                                :damage-removal-instant nil
                                                )))
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
          'cl-mpm/particle::particle-ice-delayed
          :E E
          :nu 0.3d0

          :kt-res-ratio rt
          :kc-res-ratio rc
          :residual-strength 1d0;(- 1d0 1d-3)
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
          :viscosity 1d10
          :plastic-damage-evolution *enable-plastic-damage*
          :material-damping 0d-2
          :density-degredation-max 0d0
          :biot-coeff *biot-coefficent*
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
             *bc-melange*))))

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

      (let ((cutout (+ (- ice-height water-level) bench-extra-cut))
            (cutback bench-length))
        (format t "Total cutdown ~F~%" (- ice-height cutout))
        (when (> cutback 0d0)
          (cl-mpm/setup:remove-sdf
           *sim*
           (cl-mpm/setup::rectangle-sdf (list (first block-size) (+ offset ice-height ice-height))
                                        (list cutback (+ cutout ice-height)))))))


    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)
    (when (typep *sim* 'cl-mpm/damage::mpm-sim-damage)
      (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) nil))
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
    (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
    (cl-mpm/output:add-mp-output
     *sim*
     :VECTOR
     "damage-grads"
     #'cl-mpm/particle::mp-av-damage-gradient)
    (cl-mpm/output:add-mp-output
     *sim*
     :SCALAR
     "damage-tcs-c"
     #'cl-mpm/particle::mp-damage-compression)
    (cl-mpm/output:add-mp-output
     *sim*
     :SCALAR
     "i1-undamaged"
     (lambda (mp) (cl-mpm/utils::trace-voigt (cl-mpm/particle::mp-undamaged-stress mp))))
    (cl-mpm/output:add-mp-output
     *sim*
     :SCALAR
     "density"
     (lambda (mp) (/ (cl-mpm/particle::mp-mass mp) (cl-mpm/particle::mp-volume mp))))
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
     "damage-nl-delta"
     (lambda (mp) (- (cl-mpm/particle::mp-damage-ybar mp) (cl-mpm/particle::mp-damage-y-local mp))))
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
    (cl-mpm/output:add-mp-output *sim* :SCALAR "j" #'cl-mpm/particle::mp-deformation-jacobian-strain)
    (cl-mpm/output:add-mp-output *sim* :SCALAR "water-pressure" #'cl-mpm/particle::mp-pressure)
    ;; (cl-mpm/output:add-mp-output *sim* :SCALAR "effective-total-pressure" (lambda (mp)
    ;;                                                     (-
    ;;                                                      (/ (cl-mpm/utils:trace-voigt (cl-mpm/particle::mp-undamaged-stress mp))
    ;;                                                         (* 3d0 (cl-mpm/particle::mp-deformation-jacobian-strain mp)))
    ;;                                                      (cl-mpm/particle::mp-pressure mp))))
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

;; (let ((x (loop for x from -2d0 to 2d0 by 0.01d0 collect x)))
;;   (vgplot:plot x (mapcar (lambda (x) (cl-mpm/damage::weight-func (* x x) 1d0)) x)))




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
  (cl-mpm/utils::set-workers 16)
  (let ((output-dir "./output/"))
    (ensure-directories-exist output-dir)
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
    (dolist (alpha (list 0d0))
      (dolist (friction (list 0d0))
        (dolist (notch (list 0d0))
          (dolist (aspect (list 1d0 2d0 4d0 6d0))
            (let* ((mps 3)
                   (H 900d0)
                   (ice-aspect aspect)
                   (floatation-ratio 0.0d0))
              (defparameter *alpha* alpha)
              (defparameter *length-scaler* 1d0)
              (setup
               :refine 0.125
               ;; :multigrid-refines 0
               :friction friction
               :bench-length (* notch H)
               :bench-extra-cut (* 0d0 (* H 1d0))
               :ice-height H
               :mps mps
               :hydro-static nil
               :cryo-static t
               :elastic-static nil
               :melange nil
               :aspect ice-aspect
               :slope 0d0
               :floatation-ratio floatation-ratio
               :use-penalty nil
               ;; :extra-offset 2
               :stick-base t)
              (cl-mpm/dynamic-relaxation::elastic-static-solution
               *sim*)
              (setf (cl-mpm::sim-enable-damage *sim*) t)
              (cl-mpm/damage:calculate-damage *sim* 1d0)
              (cl-mpm/output:save-vtk (uiop:merge-pathnames*
                                       (format nil "./sim_stress_~A_notch_~F_friction_~F_alpha_~F.vtk" aspect notch friction alpha)
                                       output-dir) *sim*)
              (plot-domain))))))))


(defun initial-stress ()
  (cl-mpm/utils::set-workers 12)
  (let ((output-dir "./output/"))
    (ensure-directories-exist output-dir)
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
    (dolist (initial-stress (list :CRYO-STATIC
                                  :ELASTIC-STATIC
                                  :NIL))
      (dolist (stick (list t nil))
        (let* ((mps 3)
               (output-dir (format nil "./output-stress_~A_stick_~A/" initial-stress stick))
               (H 300d0)
               (ice-aspect 4d0)
               (floatation-ratio 0.7d0))
          (defparameter *alpha* 0d0)
          (defparameter *length-scaler* 1d0)
          (format t "Running ~A~%" output-dir)
          (setup
           :refine 0.5
           :friction 0d0
           :bench-length (* 0d0 H)
           :bench-extra-cut (* 0d0 (* H 1d0))
           :ice-height H
           :mps mps
           :cryo-static (equal initial-stress :CRYO-STATIC)
           :elastic-static (equal initial-stress :ELASTIC-STATIC)
           :melange nil
           :aspect ice-aspect
           :slope 0d0
           :floatation-ratio floatation-ratio
           :use-penalty nil
           :stick-base stick)
          (setf
           (cl-mpm:sim-settings *sim*)
           (list :OCEAN-HEIGHT *water-height*))
          ;; (setf (cl-mpm/dynamic-relaxation::dt-loadstep *sim* 0d0))
          (cl-mpm/dynamic-relaxation::run-load-control
           *sim*
           :crit 1d-3
           :output-dir output-dir
           :load-steps 2
           :loading-function (lambda (p))
           :enable-damage nil
           :enable-plastic nil
           :dt-scale 0.9d0
            :post-iter-step
            (lambda (i o e)
              (setf (cl-mpm::sim-enable-damage *sim*) t)
              ;; (cl-mpm/damage:calculate-damage *sim* 1d-15)
              (setf (cl-mpm::sim-enable-damage *sim*) nil)))

          (setf (cl-mpm::sim-enable-damage *sim*) t)
          (cl-mpm/damage:calculate-damage *sim* 1d-15)
          (cl-mpm/dynamic-relaxation::save-vtks *sim* output-dir 1)
          (plot-domain)
          (vgplot:title output-dir)
          )))))

(defun test-angle-alpha (&key (angle 40d0) (alpha 0d0))
  (defparameter *alpha* alpha)
  (cl-mpm:iterate-over-mps
   (cl-mpm:sim-mps *sim*)
   (lambda (mp)
     (setf (cl-mpm/particle::mp-friction-angle mp) (cl-mpm/utils:deg-to-rad angle))))
  (setf (cl-mpm::sim-enable-damage *sim*) t)
  (cl-mpm/damage:calculate-damage *sim* 1d0)
  (plot-domain)
  )

(defun notch-stress ()
  (cl-mpm/utils::set-workers 16)
  (vgplot:close-all-plots)
  (dolist (initial-stress (list
                           ;; :CRYO-STATIC
                           :ELASTIC-STATIC
                           ;; :NIL
                           ))
    (dolist (height (list 400d0))
      (dolist (float (list 0.9d0))
        (dolist (notch-ratio (list 0.5d0))
          (dolist (friction (list 0.5d0))
            (let* ((mps 3)
                   ;; (name (format nil "stress_~A_angle_~F_friction_~F_alpha_~F" initial-stress angle friction alpha))
                   (name (format nil "height_~F_stress_~A_friction_~F_notch_~F_floatation_~F" height initial-stress friction notch-ratio float))
                   (output-dir (format nil "./output-~A/" name))
                   (H height)
                   (ice-aspect 4d0)
                   (floatation-ratio float))
              (defparameter *length-scaler* 1d0)
              (defparameter *alpha* 0.3d0)
              (defparameter *angle* 38d0)
              (format t "Running ~A~%" output-dir)
              (setup
               :refine 0.5
               :friction friction
               :bench-length (* notch-ratio H)
               :bench-extra-cut (* 0d0 (* H 1d0))
               :ice-height H
               :mps mps
               :cryo-static (equal initial-stress :CRYO-STATIC)
               :elastic-static (equal initial-stress :ELASTIC-STATIC)
               :melange nil
               :aspect ice-aspect
               :slope 0d0
               :floatation-ratio floatation-ratio
               :use-penalty t
               :stick-base nil)
              (setf
               (cl-mpm:sim-settings *sim*)
               (list :OCEAN-HEIGHT *water-height*
                     :OFFSET 2))
              (cl-mpm/output:add-mp-output
               *sim*
               :SCALAR
               "sigma-1"
               (lambda (mp)
                 (multiple-value-bind (s1 s2 s3) (cl-mpm/utils::principal-stresses-3d (cl-mpm/particle::mp-stress mp))
                   s1)))
              (cl-mpm/output:add-mp-output
               *sim*
               :SCALAR
               "sigma-3"
               (lambda (mp)
                 (multiple-value-bind (s1 s2 s3) (cl-mpm/utils::principal-stresses-3d (cl-mpm/particle::mp-stress mp))
                   s3)))
              (cl-mpm/dynamic-relaxation::run-load-control
               *sim*
               :crit 1d-3
               :output-dir output-dir
               :load-steps 1
               :loading-function (lambda (p))
               :enable-damage nil
               :enable-plastic nil
               :dt-scale 0.9d0
               :plotter
               (lambda (sim)
                 (setf (cl-mpm::sim-enable-damage *sim*) t)
                 (cl-mpm/damage:calculate-damage *sim* 1d-15)
                 (setf (cl-mpm::sim-enable-damage *sim*) nil)
                 (plot-domain)
                 )
               :post-iter-step
               (lambda (i o e)
                 ))
              (setf (cl-mpm::sim-enable-damage *sim*) t)
              (cl-mpm/damage:calculate-damage *sim* 1d-15)
              (cl-mpm/dynamic-relaxation::save-vtks *sim* output-dir 1)
              (plot-domain)

              (vgplot:title output-dir)
              (vgplot:print-plot (merge-pathnames (format nil "frame_~A.png" name)) :terminal "png size 1920,1080")
              )))))))

(defun notch-angle-stress ()
  (cl-mpm/utils::set-workers 16)
  (vgplot:close-all-plots)
  (dolist (initial-stress (list
                           ;; :CRYO-STATIC
                           :ELASTIC-STATIC
                           ;; :NIL
                           ))
    (dolist (friction (list 0.5d0))
      (let* ((mps 3)
             ;; (name (format nil "stress_~A_angle_~F_friction_~F_alpha_~F" initial-stress angle friction alpha))
             (name (format nil "stress_~A_friction_~F" initial-stress friction))
             (output-dir (format nil "./output-~A/" name))
             (H 400d0)
             (ice-aspect 4d0)
             (floatation-ratio 0.90d0))
        (defparameter *length-scaler* 1d0)
        (format t "Running ~A~%" output-dir)
        (setup
         :refine 0.5
         :friction friction
         :bench-length (* 0.5d0 H)
         :bench-extra-cut (* 0d0 (* H 1d0))
         :ice-height H
         :mps mps
         :cryo-static (equal initial-stress :CRYO-STATIC)
         :elastic-static (equal initial-stress :ELASTIC-STATIC)
         :melange nil
         :aspect ice-aspect
         :slope 0d0
         :floatation-ratio floatation-ratio
         :use-penalty t
         :stick-base nil)
        (setf
         (cl-mpm:sim-settings *sim*)
         (list :OCEAN-HEIGHT *water-height*
               :OFFSET 2))
        (cl-mpm/output:add-mp-output
         *sim*
         :SCALAR
         "sigma-1"
         (lambda (mp)
           (multiple-value-bind (s1 s2 s3) (cl-mpm/utils::principal-stresses-3d (cl-mpm/particle::mp-stress mp))
             s1)))
        (cl-mpm/output:add-mp-output
         *sim*
         :SCALAR
         "sigma-3"
         (lambda (mp)
           (multiple-value-bind (s1 s2 s3) (cl-mpm/utils::principal-stresses-3d (cl-mpm/particle::mp-stress mp))
             s3)))
        (cl-mpm/dynamic-relaxation::run-load-control
         *sim*
         :crit 1d-3
         :output-dir output-dir
         :load-steps 1
         :loading-function (lambda (p))
         :enable-damage nil
         :enable-plastic nil
         :dt-scale 0.9d0
         :plotter
         (lambda (sim)
           (setf (cl-mpm::sim-enable-damage *sim*) t)
           (cl-mpm/damage:calculate-damage *sim* 1d-15)
           (setf (cl-mpm::sim-enable-damage *sim*) nil)
           (plot-domain)
           )
         :post-iter-step
         (lambda (i o e)
           ))
        (setf (cl-mpm::sim-enable-damage *sim*) t)
        (cl-mpm/damage:calculate-damage *sim* 1d-15)
        (cl-mpm/dynamic-relaxation::save-vtks *sim* output-dir 1)
        (plot-domain)

        (vgplot:title output-dir)
        (vgplot:print-plot (merge-pathnames (format nil "frame_~A.png" name)) :terminal "png size 1920,1080")

        (dolist (angle (list 40d0 30d0 20d0))
          (dolist (alpha (list 0d0 0.25d0 0.5d0 0.75d0 1d0))
            (defparameter *angle* angle)
            (defparameter *alpha* alpha)
            (let* ((name (format nil "stress_~A_angle_~F_friction_~F_alpha_~F" initial-stress angle friction alpha))
                   (output-dir (merge-pathnames (format nil "./output-~A/" name))))
              (uiop:ensure-all-directories-exist (list output-dir))
              (cl-mpm/output::save-simulation-parameters (merge-pathnames output-dir "settings.json")
                                                         *sim*
                                                         (list))
              (cl-mpm:iterate-over-mps
               (cl-mpm:sim-mps *sim*)
               (lambda (mp)
                 (setf (cl-mpm/particle::mp-friction-angle mp) (cl-mpm/utils:deg-to-rad angle))))
              (setf (cl-mpm::sim-enable-damage *sim*) t)
              (cl-mpm/damage:calculate-damage *sim* 1d0)
              (cl-mpm/dynamic-relaxation::save-vtks *sim* output-dir 0)
              (plot-domain)
              (vgplot:title output-dir)
              (vgplot:print-plot (merge-pathnames (format nil "frame_~A.png" name)) :terminal "png size 1920,1080")
              )))))))



;; (let ((i 0))
;;   (cl-mpm::iterate-over-mps-serial
;;    (cl-mpm:sim-mps *sim*)
;;    (lambda (mp)
;;      (when (cl-mpm/particle::mp-penalty-contact-point mp)
;;        (format t "~D - ~A~%" i (cl-mpm/particle::mp-penalty-frictional-force mp))
;;        (cl-mpm/particle::iterate-over-mp-corners
;;         (cl-mpm:sim-mesh *sim*)
;;         mp
;;         (lambda (c)
;;           (when (cl-mpm/particle::corner-contact c)
;;             (format t "corner ~A - ~A ~%"
;;                     (cl-mpm/particle::corner-contact c)
;;                     (cl-mpm/particle::corner-penalty-frictional-force c))))))
;;      (incf i))))



;; (let ((iter 1000))
;;   (time (dotimes (i iter) (make-instance 'cl-mpm/particle::particle-damage)))
;;   (time (dotimes (i iter) (allocate-instance (find-class 'cl-mpm/particle::particle-damage)))))





