(defpackage :cl-mpm/examples/ice-buoyancy
  (:use :cl
   :cl-mpm/example))
(in-package :cl-mpm/examples/ice-buoyancy)

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
(defparameter *rt* (- 1d0 1d0))

(defparameter *pd-oversize* 1d-2)
(defparameter *rc*
  (- 1d0 *pd-oversize*))

(defparameter *enable-plastic-damage* nil)
(defparameter *delay-time* 1d4)
(defparameter *delay-exponent* 2d0)
(defparameter *enable-viscosity* nil)
(defparameter *length-scaler* 2d0)
(defparameter *gf* 10000d0)
(defparameter *ductility* 10d0)
(defparameter *tensile-strength* 0.2d6)
(defparameter *biot-coefficent* 1d0)
(defparameter *alpha* 1d0)
;; (defparameter *alpha* 0.5d0)
(defparameter *material-damping* 0d-3)

;; (defparameter *alpha* 0.4d0)


(defclass cl-mpm/particle::particle-ice-erodable (cl-mpm/particle::particle-ice-delayed
                                                    cl-mpm/particle::particle-erosion)
  ())


(defmethod cl-mpm/erosion::mp-erosion-enhancment ((mp cl-mpm/particle::particle-ice-erodable))
  ;; (+ 1d0 (* 10 (cl-mpm/particle::mp-damage mp)))
  (expt (the double-float (cl-mpm/particle::mp-damage mp)) 2)
  ;; (+ 1d0 (* 10 (cl-mpm/particle::mp-strain-plastic-vm mp)))
  )

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-ice-brittle) dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-damaged mesh mp dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-dynamic-relaxation mesh mp dt fbar)
  (cl-mpm::update-stress-kirchoff mesh mp dt fbar))

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
        (declare (double-float E ps-vm angle pressure j damage))
        (progn
          (let* ((ps-y (the double-float
                            (* E (max 0d0 ps-vm-inc))
                            ;; (sqrt (* E (expt ps-vm-inc 2)))
                            ;; (sqrt (* E ps-vm))
                            ;; (sqrt (* E ps-vm-inc))
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
                      (- pressure)
                      ;; (the double-float (- 1d0 (* 2d0 damage)))
                      )
                     (grab-new-voigt))
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
               ;; (cl-mpm/damage::criterion-mohr-coloumb-rankine-stress-tensile stress-pressure angle)
               (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile stress-pressure angle)
               ;; (cl-mpm/damage::criterion-mohr-coloumb-gill strain angle e nu)
               ;; (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile stress-pressure angle)
               ;; (cl-mpm/damage::drucker-prager-criterion stress angle)
               ;; (cl-mpm/damage::drucker-prager-criterion stress-pressure angle)
               ))))))))

(defparameter *water-height* 0d0)
(defparameter *offset* 0d0)
(declaim (notinline plot-domain))
(defun plot-domain (&key (trial t))
  ;; (plot-pq)
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
     ))
  )

(defparameter *bc-melange* nil)
(defparameter *bc-erode* nil)
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
                (undercut-length 0d0)
                (aspect 1)
                (floatation-ratio 0.9)
                (slope 0d0)
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
                                                :enable-fbar t
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
                                                :damage-removal-instant nil)))
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
             (rs (est-shear-from-angle angle *angle-r* rc)))

        (let* ((pd (- 1d0 *pd-oversize*))
               (k (cl-mpm/damage::find-k-damage E init-stress ductility pd))
               (ds (cl-mpm/damage::damage-response-exponential-peerlings-residual k E init-stress ductility rs))
               (dc (cl-mpm/damage::damage-response-exponential-peerlings-residual k E init-stress ductility rc)))
          (format t "Damage pd ~E~%" pd)
          (format t "Damage t ~E~%" (cl-mpm/damage::damage-response-exponential-peerlings-residual k E init-stress ductility rt))
          (format t "Damage c ~E~%" dc)
          (format t "Damage s ~E~%" ds)
          (format t "Real residual angle ~E~%" (cl-mpm/utils::rad-to-deg (atan (* (/ (- 1d0 ds) (- 1d0 dc)) (tan (cl-mpm/utils::deg-to-rad *angle*))))))
          (format t "0 dc residual angle ~E~%" (cl-mpm/utils::rad-to-deg (atan (* (- 1d0 ds)  (tan (cl-mpm/utils::deg-to-rad *angle*))))))
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
          ;; 'cl-mpm/particle::particle-ice-brittle
          :E E
          :nu 0.3d0

          :kt-res-ratio rt
          :kc-res-ratio rc
          :residual-strength 1d0;(- 1d0 1d-3)
          :initiation-stress init-stress
          :friction-angle (cl-mpm/utils:deg-to-rad angle)
          :residual-friction (cl-mpm/utils:deg-to-rad *angle-r*)

          :psi (cl-mpm/utils::deg-to-rad *angle-psi*)
          :softening 0d0
          :oversize (- 1d0 *pd-oversize*)
          :ductility ductility
          :local-length length-scale
          :delay-time *delay-time*
          :delay-exponent *delay-exponent*
          :enable-plasticity t
          :enable-damage t
          :enable-viscosity *enable-viscosity*
          :viscosity 1d14
          :plastic-damage-evolution *enable-plastic-damage*
          :material-damping *material-damping*
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
         :force-3d t
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
    (let ((cutout water-level)
          (cutback undercut-length))
      ;; (format t "Total undercut ~F~%" (- ice-height cutout))
      (when (> cutback 0d0)
        (cl-mpm/setup:remove-sdf
         *sim*
         (cl-mpm/setup::rectangle-sdf (list (first block-size) offset)
                                      (list cutback (+ cutout))))))


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
        (cl-mpm/penalty::make-bc-penalty-point-normal
         *sim*
         (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
         (cl-mpm/utils:vector-from-list (list
                                         domain-half
                                         offset
                                         0d0))
         ;; (* domain-half 1.1d0)
         (* E epsilon-scale)
         friction
         penalty-damping)))
    ;; (setf (cl-mpm/penalty::bc-penalty-stiffness-scale *floor-bc*) 1d0)

    (defparameter *bc-erode*
      (cl-mpm/erosion::make-bc-erode-uniform
       *sim*
       :enable nil
       :rate 1d-4)
      ;; (cl-mpm/erosion::make-bc-erode
      ;;  *sim*
      ;;  :enable nil
      ;;  :rate 1d0
      ;;  :scalar-func (lambda (pos)
      ;;                 1d0
      ;;                 ;; (min 1d0 (exp (* 0.5d0 (- (cl-mpm/utils:varef pos 1) datum))))
      ;;                 )
      ;;  :clip-func (lambda (pos)
      ;;               (and
      ;;                (>= datum (cl-mpm/utils:varef pos 1))
      ;;                (<= (- datum (* 0.25d0 start-height)) (cl-mpm/utils:varef pos 1))
      ;;                ;; (>= (cl-mpm/utils:varef pos 1) (+ offset mesh-resolution) )
      ;;                )))
      )
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
    ;;  *bc-erode*)
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
           (if (> (cl-mpm/constitutive::voight-trace (cl-mpm/particle::mp-undamaged-stress mp)) 0d0)
               (* (/ 180 pi) (atan (* (/ (- 1d0 (cl-mpm/particle::mp-damage-shear mp))
                                         1d0
                                         ;; (- 1d0 (cl-mpm/particle::mp-damage-tension mp))
                                         )
                                      (tan (cl-mpm/particle::mp-phi mp)))))
               (* (/ 180 pi) (atan (* (/ (- 1d0 (cl-mpm/particle::mp-damage-shear mp))
                                         1d0
                                         ;; (- 1d0 (cl-mpm/particle::mp-damage-compression mp))
                                         )
                                      (tan (cl-mpm/particle::mp-phi mp))))))
           0d0)))
    (cl-mpm/output:add-mp-output
     *sim*
     :SCALAR
     "current-cohesion"
     (lambda (mp)
       (if (typep mp 'cl-mpm/particle::particle-ice-brittle)
           (*
            (if (> (cl-mpm/constitutive::voight-trace (cl-mpm/particle::mp-undamaged-stress mp)) 0d0)
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

(defmethod cl-mpm/dynamic-relaxation::damage-increment-criteria ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-ul))
  ;; (cl-mpm/dynamic-relaxation::compute-max-damage-energy-crit sim)
  ;; (cl-mpm/dynamic-relaxation::compute-max-damage-energy-crit-mp sim)
  (cl-mpm/dynamic-relaxation::damage-increment-criteria-mp sim)
  ;; (damage-increment-criteria-mesh sim)
  )
(defmethod cl-mpm/dynamic-relaxation::damage-increment-criteria ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic))
  (cl-mpm/dynamic-relaxation::damage-increment-criteria-mp sim)
  ;; (cl-mpm/dynamic-relaxation::compute-max-damage-energy-crit-mp sim)
  ;; (cl-mpm/dynamic-relaxation::compute-max-damage-energy-crit sim)
  )

(defun reduced-output-data (sim)
  (cl-mpm/output::reset-mp-output sim)
  (cl-mpm/output::add-mp-output sim :SCALAR "split-depth" #'cl-mpm/particle::mp-split-depth)
  (cl-mpm/output::add-mp-output sim :SCALAR "yield-function" (lambda (mp) (if (slot-exists-p mp 'cl-mpm/particle::yield-func) (cl-mpm/particle::mp-yield-func mp) 0d0)))
  (cl-mpm/output::add-mp-output sim :SCALAR "plastic-strain" (lambda (mp) (if (slot-exists-p mp 'cl-mpm/particle::strain-plastic-vm) (cl-mpm/particle::mp-strain-plastic-vm mp) 0d0)))
  (cl-mpm/output::add-mp-output sim :VECTOR "disp" #'cl-mpm/particle::mp-displacement)
  (cl-mpm/output::add-mp-output sim :SCALAR "fric-normal" #'cl-mpm/particle::compute-normal-force)
  (cl-mpm/output::add-mp-output sim :VOIGT "sig" #'cl-mpm/particle::mp-stress)

  (macrolet ((damage-val (mp &body body)
               `(if (typep ,mp 'cl-mpm/particle::particle-damage)
                    ,@body
                    0d0))
             (has-slot-val (mp slot &body body)
               `(if (typep ,mp 'cl-mpm/particle::particle-damage)
                    ,@body
                    0d0))
             )
    (cl-mpm/output::add-mp-output sim :SCALAR "damage-inc" (lambda (mp) (damage-val mp (cl-mpm/particle::mp-damage-increment mp))))
    (cl-mpm/output::add-mp-output sim :SCALAR "damage" (lambda (mp) (damage-val mp (cl-mpm/particle::mp-damage mp))))
    (cl-mpm/output::add-mp-output sim :SCALAR "damage-ybar" (lambda (mp) (damage-val mp (cl-mpm/particle::mp-damage-ybar mp)))))
  )

(defun calving-test (&key (output-dir "./output/"))
  (cl-mpm/utils::set-workers 16)
  (let* ((mps 3)
         (dt 1d3)
         (total-time 1d10)
         (H 400d0)
         (ice-aspect 4d0)
         (density 918d0)
         (explicit-dt-scale 0.50d0)
         (water-damping 100d0)
         (friction 0.5d0)
         (floatation-ratio 0.76d0)
         )
    (defparameter *length-scaler* 3d0)
    (setup
     :refine 0.25
     :friction friction
     :bench-length (* 0d0 H)
     :bench-extra-cut (* 0d0 (* H 1d0))
     :undercut-length (* 0d0 H)
     :ice-height H
     :mps mps
     :hydro-static nil
     :cryo-static t
     :elastic-static nil
     :melange nil
     :aspect ice-aspect
     :slope 0d0
     :floatation-ratio floatation-ratio
     :use-penalty t
     :stick-base nil)

    (cl-mpm/output:add-mp-output *sim* :SCALAR "eroded"
                                 (lambda (mp)
                                   (/ (cl-mpm/particle::mp-eroded-volume mp) (cl-mpm/particle::mp-mass mp))))

    (cl-mpm/output:add-mp-output *sim* :SCALAR "boundary" #'cl-mpm/particle::mp-boundary)
    (cl-mpm/output:add-mp-output *sim* :SCALAR "plastic-c" #'cl-mpm/particle::mp-c)
    (cl-mpm/output:add-mp-output *sim* :VECTOR "body-force" #'cl-mpm/particle::mp-body-force)
    (cl-mpm/output:add-node-output *sim* :VECTOR "boundary-vec" #'cl-mpm/mesh::node-boundary-vec)

    (reduced-output-data *sim*)
    (cl-mpm::domain-sort-mps *sim*)

    (when (typep *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-octree)
      (setf (cl-mpm/dynamic-relaxation::sim-intra-mesh-aggregation *sim*) t)
      (setf (cl-mpm/dynamic-relaxation::sim-octree-refinement-criteria *sim*)
            (lambda (sim mesh c)
              (or
               (damage-refinement-criteria sim mesh c)))))


    (plot-domain)
    ;; (break)

    (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0)

    (setf (cl-mpm/damage::sim-enable-stress-based-length *sim*) nil)
    (setf (cl-mpm/damage::sim-enable-ekl *sim*) nil)
    (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
    ;; (setf (cl-mpm::sim-nonlocal-damage *sim*) nil)

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
                :LENGTH-SCALER *length-scaler*
                :TENSILE-STRENGTH *tensile-strength*))

    (setf (cl-mpm:sim-enable-damage *sim*) nil)

    (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)

    ;; (break)
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (let ((step 0)
          (substeps (ceiling (* 10 (floor H 100) (/ 10d0 (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*))))))
          )
      (format t "Substeps ~D~%" substeps)
      ;; (setf (cl-mpm/penalty::bc-penalty-friction *floor-bc*) 0d0)

      (cl-mpm/dynamic-relaxation::run-multi-stage
       *sim*
       :output-dir output-dir
       :dt dt
       :conv-dt-scale 0.9d0
       :dt-scale 0.9d0
       :damping-factor (sqrt 2d0)
       :conv-criteria 1d-6
       :conv-load-steps 1
       ;; :min-adaptive-steps -4
       ;; :max-adaptive-steps 10
       :min-adaptive-steps -14
       :max-adaptive-steps 14
       :adaption-constant 4
       :easy-step-adaption-constant 4
       :max-damage-inc 1.9d0
       :min-tangent-ratio 1d-1
       :max-deformation-gradient 2d0
       :max-plastic-inc nil
       :max-inertia-norm 1d-4
       :stagger-damage :HYBRID-FULL
       ;; :MONOLITH-QS
       ;; :stagger-damage :MONOLITH-QS
       ;; :stagger-damage :FULL
       ;; :stagger-damage :HYBRID-FULL
       ;; :min-damage-inc 0.005d0
       :substeps substeps
       :sub-conv-steps 500
       :total-time total-time
       :save-vtk-loadstep t
       :save-vtk-dr t
       :enable-plastic t
       :enable-damage t
       :plotter (lambda (sim)
                  (plot-domain)
                  (vgplot:title (format nil "Step ~D - Time ~F - oobf ~E - ~A"
                                        step
                                        (cl-mpm::sim-time sim)
                                        (cl-mpm::sim-stats-oobf sim)
                                        (if (typep *sim* 'cl-mpm::mpm-sim-usf)
                                            "Explicit"
                                            "Implicit")
                                        ;; (if (equal (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
                                        ;;     "Quasi-Static"
                                        ;;     "Dynamic")
                                        ))
                  ;; (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step)) :terminal "png size 1920,1080")
                  (incf step))
       :explicit-conv-criteria 1d-3
       :elastic-dt-margin 1d2
       :explicit-mass-scaling nil
       :explicit-damping-factor 1d-6
       :explicit-dt-scale explicit-dt-scale
       :explicit-dynamic-solver 'cl-mpm/damage::mpm-sim-agg-damage

       ;; :explicit-dynamic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-octree-damage-usf
       ;; :elastic-dt-margin 1d4
       ;; :explicit-mass-scaling nil
       ;; :explicit-dt-scale 10d0
       ;; :explicit-damping-factor 0d-3
       ;; ;; :explicit-dynamic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-octree-implicit-dynamic
       ;; :explicit-dynamic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic
       ;; :elastic-solver ;'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static
       ;; 'cl-mpm/dynamic-relaxation::mpm-sim-octree-quasi-static
       ;; :initial-quasi-static t
       :post-conv-step
       (lambda (sim)
         ;; (cl-mpm::domain-sort-mps *sim*)
         (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) water-damping
         (setf (cl-mpm/bc::bc-enable *bc-erode*) nil))
       :setup-quasi-static
       (lambda (sim)
         ;(cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
         (setf
          (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC
          )
         (when (typep *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic)

           (setf (cl-mpm/dynamic-relaxation::sim-true-damping *sim*) (* 0d-6 (cl-mpm/setup::estimate-critical-damping *sim*))
                 (cl-mpm::sim-velocity-algorithm sim) :TBLEND
                 ))
         ;; (cl-mpm::remove-mps-func
         ;;  *sim*
         ;;  (lambda (mp)
         ;;    (and
         ;;     (typep mp 'cl-mpm/particle::particle-damage)
         ;;     (> (cl-mpm/particle::mp-damage mp) 0.9d0))))
         ;; (cl-mpm::reset-grid (cl-mpm:sim-mesh *sim*) :reset-displacement t)
         ;; (cl-mpm/dynamic-relaxation::pre-step *sim*)
         (cl-mpm::check-mps *sim*)
         (setf
          (cl-mpm/penalty::bc-penalty-friction *floor-bc*) friction
          (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) water-damping))
       :setup-dynamic
       (lambda (sim)
         ;(cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
         (setf
          (cl-mpm/damage::sim-damage-delocal-counter-max sim) 2
          (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) water-damping
          (cl-mpm::sim-velocity-algorithm sim) :TBLEND))))))

(defun calving-real-test ()
  (vgplot:close-all-plots)
  (cl-mpm/utils:set-workers 16)
  (let* ((mps 3)
         (H 600d0)
         (density 918d0)
         (water-damping 50d0)
         )
    (setf *delay-time* 1d2)
    (defparameter *length-scaler* 2d0)
    (setup :refine 0.125
           :friction 0.5d0
           :bench-length (* 0.5d0 H)
           :bench-extra-cut (* 0d0 (* H 1d0))
           :ice-height H
           :mps mps
           :hydro-static nil
           :cryo-static nil
           :elastic-static t
           :melange nil
           :aspect 4d0
           :slope 0d0
           :floatation-ratio 0.9d0
           ;; :floatation-ratio 1.00d0
           :use-penalty t
           :extra-offset 0
           :stick-base t
           )
    (change-class *sim* 'cl-mpm/damage::mpm-sim-agg-damage)
    (setf (cl-mpm::sim-max-split-depth *sim*) 0)
    ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic)
       ;; ;; :explicit-dynamic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-octree-implicit-dynamic
    ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-octree-damage-usf)
    ;; (setf (cl-mpm/dynamic-relaxation::sim-intra-mesh-aggregation *sim*) t)
    ;; (change-class *sim* 'cl-mpm/damage::mpm-sim-agg-damage)
    (plot-domain)
    (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) water-damping)
    (setf
     (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
     (cl-mpm::sim-ghost-factor *sim*) nil)

    (setf (cl-mpm/damage::sim-enable-ekl *sim*) nil)
    (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
    (setf (cl-mpm:sim-enable-damage *sim*) nil)

    (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (setf (cl-mpm::sim-velocity-algorithm *sim*) :TBLEND)
    ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic)
    (let ((step 0))
      (cl-mpm/dynamic-relaxation::run-time
       *sim*
       :output-dir "./output/"
       :dt 1d0
       :total-time 1d5
       ;; :dt-scale 1000d0
       :dt-scale 0.5d0
       :mass-scale 1d0
       :damping 1d-4
       :enable-plastic t
       :enable-damage t
       :conv-criteria 1d-3
       :save-vtk-loadstep t
       :initial-quasi-static t
       ;; :elastic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-octree-damage-quasi-static
       :elastic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul
       :post-conv-step
       (lambda (sim)
         (cl-mpm:iterate-over-mps
          (cl-mpm:sim-mps *sim*)
          (lambda (mp)
            (cl-mpm/fastmaths::fast-zero (cl-mpm/particle::mp-displacement mp))))
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



(defun calving-qs-test ()
  (cl-mpm/utils:set-workers 16)
  (let* ((mps 3)
         (H 900d0)
         (density 918d0))
    (setf *delay-time* 1d5)
    (defparameter *length-scaler* 2d0)
    (setup :refine 0.125
           ;; :multigrid-refines 1
           :friction 0.5d0
           :bench-length (* 0.5d0 H)
           :bench-extra-cut 00d0
           :ice-height H
           :mps mps
           :hydro-static nil
           :cryo-static nil
           :elastic-static t
           :melange nil
           :aspect 4d0
           :slope 0.0d0
           :floatation-ratio 0.60d0
           :use-penalty t
           ;; :stick-base t
           )
    ;; (cl-mpm/output:add-mp-output
    ;;  *sim*
    ;;  :VECTOR
    ;;  "bf"
    ;;  #'cl-mpm/particle::mp-body-force)
    ;; (push (list :SCALAR "water-pressure" #'cl-mpm/particle::mp-pressure) (cl-mpm::sim-output-list *sim*))
    ;; (push (list :SCALAR "water-pressure-damage"
    ;;             (lambda (mp) (*
    ;;                           (cl-mpm/particle::mp-damage mp)
    ;;                           (cl-mpm/particle::mp-pressure mp)))) (cl-mpm::sim-output-list *sim*))
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
    (setf (cl-mpm/damage::sim-enable-ekl *sim*) nil)
    (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
    (setf (cl-mpm:sim-enable-damage *sim*) t)
    (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (setf (cl-mpm::sim-velocity-algorithm *sim*) :TPIC)
    ;; (setf (cl-mpm/dynamic-relaxation::sim-enable-dynamics *sim*) nil)
    (let ((step 0))
      (cl-mpm/dynamic-relaxation::run-quasi-time
       *sim*
       :output-dir "./output/"
       :dt 1d4
       :total-time 1d7
       :dt-scale 0.9d0
       :enable-plastic t
       :enable-damage t
       :substeps 50
       :sub-conv-steps 10
       ;; :min-adaptive-steps 0
       ;; :max-adaptive-steps 0
       :min-adaptive-steps -12
       :max-adaptive-steps 12
       :adaption-constant 4
       :adaption-easy-steps 8
       :max-damage-inc 1.1d0
       :max-plastic-inc nil
       :max-deformation-gradient 10d0

       :conv-criteria 1d-3
       :save-vtk-loadstep t
       :elastic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static
       :initial-quasi-static t
       :post-conv-step
       (lambda (sim)
         (setf (cl-mpm/dynamic-relaxation::sim-true-damping *sim*) (* 1d-4 (cl-mpm/setup::estimate-critical-damping *sim*)))
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
          (dolist (aspect (list 2d0))
            (let* ((mps 3)
                   (H 600d0)
                   (ice-aspect aspect)
                   (floatation-ratio 1d0))
              (defparameter *alpha* alpha)
              (defparameter *length-scaler* 1d0)
              (setup
               :refine 0.25
               ;; :multigrid-refines 0
               :friction friction
               :bench-length (* notch H)
               :bench-extra-cut (* 0d0 (* H 1d0))
               :ice-height H
               :mps mps
               :hydro-static nil
               :cryo-static nil
               :elastic-static t
               :melange nil
               :aspect ice-aspect
               :slope 0d0
               :floatation-ratio floatation-ratio
               :use-penalty nil
               :extra-offset 2
               :stick-base nil)
              ;; (cl-mpm/dynamic-relaxation::elastic-static-solution
              ;;  *sim*)

              (cl-mpm/dynamic-relaxation::run-elastic
               *sim*
               :crit 1d-9
               :dt-scale 0.9d0
               )
              ;; (setf (cl-mpm::sim-enable-damage *sim*) t)
              ;; (cl-mpm/damage:calculate-damage *sim* 1d0)
              ;; (cl-mpm/output:save-vtk (uiop:merge-pathnames*
              ;;                          (format nil "./sim_stress_~A_notch_~F_friction_~F_alpha_~F.vtk" aspect notch friction alpha)
              ;;                          output-dir) *sim*)
              (plot-domain))))))))


(defun initial-stress ()
  (cl-mpm/utils::set-workers 16)
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
              (cl-mpm/damage:calculate-damage *sim* 1d-15)
              (setf (cl-mpm::sim-enable-damage *sim*) nil)))
          ;; (cl-mpm/dynamic-relaxation::run-elastic
          ;;  *sim*
          ;;  :conv 1d-3
          ;;  :output-dir output-dir
          ;;  :post-iter-step
          ;;  (lambda (i o e)
          ;;    (setf (cl-mpm::sim-enable-damage *sim*) t)
          ;;    (cl-mpm/damage:calculate-damage *sim* 1d-15)
          ;;    (setf (cl-mpm::sim-enable-damage *sim*) nil))
          ;; )
          (setf (cl-mpm::sim-enable-damage *sim*) t)
          (cl-mpm/damage:calculate-damage *sim* 1d0)
          (cl-mpm/dynamic-relaxation::save-vtks *sim* output-dir 1)
          (plot-domain)
          (break)
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
                           :CRYO-STATIC
                           ;; :ELASTIC-STATIC
                           ;; :NIL
                           ))
    (dolist (height (list 300d0 500d0 700d0 900d0))
      (dolist (float (list 0.5d0 0.75d0 0.9d0))
        (dolist (notch-ratio (list 0.5d0 1d0))
          (dolist (friction (list 1d0))
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
                           :CRYO-STATIC
                           ;; :ELASTIC-STATIC
                           ;; :NIL
                           ))
    (dolist (friction (list 0.5d0))
      (let* ((mps 3)
             ;; (name (format nil "stress_~A_angle_~F_friction_~F_alpha_~F" initial-stress angle friction alpha))
             (name (format nil "stress_~A_friction_~F" initial-stress friction))
             (output-dir (format nil "./output-~A/" name))
             (H 900d0)
             (ice-aspect 4d0)
             (floatation-ratio 0.76d0))
        (defparameter *length-scaler* 1d0)
        (format t "Running ~A~%" output-dir)
        (setup
         :refine 0.125
         :friction friction
         :bench-length (* 1d0 H)
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

        (dolist (angle (list 50d0 40d0 30d0 20d0))
          (dolist (alpha (list 0d0 0.2d0 0.5d0 0.75d0 1d0))
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

(defun calving-test-sweep ()
  (dolist (f (list 0.6d0 0.7d0 0.8d0))
    (cl-mpm/utils::set-workers 16)
    (let* ((mps 4)
           (dt 1d3)
           (total-time 1d9)
           (H 400d0)
           (ice-aspect 2d0)
           (density 918d0)
           (explicit-dt-scale 0.5d0)
           (water-damping 1d0)
           (friction 0.5d0)
           (floatation-ratio f)
           (output-dir (format nil "./output-~F/" f)))
      (defparameter *length-scaler* 2d0)
      (setup
       :refine 0.25
       ;; :multigrid-refines 0
       :friction friction
       :bench-length (* 0d0 H)
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
       :use-penalty t
       ;; :extra-offset 2
       :stick-base nil)

      (cl-mpm/output:add-mp-output *sim* :SCALAR "def-aspect" #'cl-mpm/dynamic-relaxation::compute-deformation-aspect-2d)
      (cl-mpm/output:add-mp-output *sim* :SCALAR "eroded"
                                   (lambda (mp)
                                     (/ (cl-mpm/particle::mp-eroded-volume mp) (cl-mpm/particle::mp-mass mp))))
      (cl-mpm/output:add-mp-output *sim* :SCALAR "boundary" #'cl-mpm/particle::mp-boundary)
      (cl-mpm/output:add-mp-output *sim* :SCALAR "vm-plastic-inc" #'cl-mpm/particle::mp-strain-plastic-vm-inc)

      (cl-mpm/output:add-node-output
       *sim*
       :SCALAR
       "volume"
       (lambda (mp) (cl-mpm/mesh::node-volume mp)))
      (cl-mpm/output:add-node-output
       *sim*
       :SCALAR
       "volume-ratio"
       (lambda (mp) (/
                     (cl-mpm/mesh::node-volume mp)
                     (cl-mpm/mesh::node-volume-true mp))))
      (cl-mpm/output:add-node-output
       *sim*
       :VECTOR
       "inertia"
       #'cl-mpm/mesh::node-inertia-force)
      (cl-mpm::domain-sort-mps *sim*)

      (when (typep *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-octree)
        (setf (cl-mpm/dynamic-relaxation::sim-intra-mesh-aggregation *sim*) t)
        (setf (cl-mpm/dynamic-relaxation::sim-octree-refinement-criteria *sim*)
              (lambda (sim mesh c)
                (or
                 (damage-refinement-criteria sim mesh c)))))


      (plot-domain)

      (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0)
      (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
            (cl-mpm::sim-ghost-factor *sim*) nil)

      (setf (cl-mpm/damage::sim-enable-stress-based-length *sim*) nil)
      (setf (cl-mpm/damage::sim-enable-ekl *sim*) nil)
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

      (setf (cl-mpm:sim-enable-damage *sim*) nil)

      (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-9)

      ;; (break)
      (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
      (let ((step 0))
        ;; (setf (cl-mpm/penalty::bc-penalty-friction *floor-bc*) 0d0)
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
         :adaption-constant 4
         :max-damage-inc 0.9d0
         :max-deformation-gradient 2d0
         :max-plastic-inc nil
         ;; :min-damage-inc 0.005d0
         :substeps (* 2 (floor H 200) (floor (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)) 10d0))
         :sub-conv-steps 50
         :total-time total-time
         :save-vtk-loadstep t
         :save-vtk-dr t
         :enable-plastic t
         :enable-damage t
         :plotter (lambda (sim)
                    ;; (format t "Agg CFL ~E - ~E~%" (cl-mpm/aggregate::estimate-aggregated-cfl *sim*) (cl-mpm/setup::estimate-elastic-dt *sim*))
                    (plot-domain)
                    (vgplot:title (format nil "Step ~D - Time ~F - oobf ~E - ~A"
                                          step
                                          (cl-mpm::sim-time sim)
                                          (cl-mpm::sim-stats-oobf sim)
                                          (if (typep *sim* 'cl-mpm::mpm-sim-usf)
                                              "Explicit"
                                              "Implicit")
                                          ;; (if (equal (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
                                          ;;     "Quasi-Static"
                                          ;;     "Dynamic")
                                          ))
                    (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step)) :terminal "png size 1920,1080")
                    (incf step))
         ;; :explicit-conv-criteria 1d-2

         :elastic-dt-margin 1d2
         :explicit-mass-scaling nil
         :explicit-damping-factor 1d-4
         :explicit-dt-scale explicit-dt-scale
         :explicit-dynamic-solver 'cl-mpm/damage::mpm-sim-agg-damage

         ;; :explicit-dynamic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-octree-damage-usf
         ;; :elastic-dt-margin 1d4
         ;; :explicit-mass-scaling nil
         ;; :explicit-dt-scale 10d0
         ;; :explicit-damping-factor 0d-3
         ;; ;; :explicit-dynamic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-octree-implicit-dynamic
         ;; :explicit-dynamic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic
         ;; :elastic-solver ;'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static
         ;; 'cl-mpm/dynamic-relaxation::mpm-sim-octree-quasi-static
         ;; :initial-quasi-static t
         :post-conv-step
         (lambda (sim)
           ;; (cl-mpm::domain-sort-mps *sim*)
           (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) water-damping
           (setf (cl-mpm/bc::bc-enable *bc-erode*) nil))
         :setup-quasi-static
         (lambda (sim)
           (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-9)
           (when (typep *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic)
             (setf (cl-mpm/dynamic-relaxation::sim-true-damping *sim*) (* 1d-4 (cl-mpm/setup::estimate-critical-damping *sim*))))
           (cl-mpm::remove-mps-func
            *sim*
            (lambda (mp)
              (and
               (typep mp 'cl-mpm/particle::particle-damage)
               (> (cl-mpm/particle::mp-damage mp) 0.999d0))))
           (cl-mpm::reset-grid (cl-mpm:sim-mesh *sim*) :reset-displacement t)
           (cl-mpm/dynamic-relaxation::pre-step *sim*)
           (cl-mpm::check-mps *sim*)
           (setf
            (cl-mpm/penalty::bc-penalty-friction *floor-bc*) friction
            (cl-mpm::sim-velocity-algorithm sim) :TBLEND
            (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) water-damping))
         :setup-dynamic
         (lambda (sim)
           (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
           (setf
            (cl-mpm/damage::sim-damage-delocal-counter-max sim) 10
            (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) water-damping
            (cl-mpm::sim-velocity-algorithm sim) :TBLEND)))))))

(defmacro time-form (it form)
  `(progn
     (declaim (optimize speed))
     (let* ((iterations ,it)
            (start (get-internal-real-time))
            (gc-start sb-ext:*gc-real-time*)
            )
       (time
        (dotimes (i ,it)
          ,form))

       (let* ((end (get-internal-real-time))
              (gc-end sb-ext:*gc-real-time*)

              (units internal-time-units-per-second)
              (dt (/ (- end start) (* iterations units))))
         (format t "Total time: ~f ~%" (/ (- end start) units))
         (format t "Time per iteration: ~f~%" (/ (- end start) (* iterations units)))
         (format t "Total gc time: ~f ~%" (/ (- gc-end gc-start) units))
         (format t "Throughput: ~f~%" (/ 1 dt))
         (format t "Time per MP: ~E~%" (/ dt (length (cl-mpm:sim-mps *sim*))))
         (format t "MP Throughput: ~E~%" (/ (length (cl-mpm:sim-mps *sim*)) dt))
         dt))))

(defun profile ()
  ;; (setup :refine 16)
  (sb-profile:unprofile)
  (sb-profile:profile "CL-MPM")
  (sb-profile:profile "CL-MPM/PARTICLE")
  (sb-profile:profile "CL-MPM/MESH")
  (sb-profile:profile "CL-MPM/SHAPE-FUNCTION")
  (sb-profile:reset)
  (time-form
   10
   (progn
     (cl-mpm::update-sim *sim*)))
  (format t "MPS ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (sb-profile:report))




;; (progn
;;   (time
;;    (dotimes (i 100000000)
;;      (cl-mpm/mesh::get-node-values (cl-mpm:sim-mesh *sim*) 0 0 0)))
;;   (time
;;    (dotimes (i 100000000)
;;      (cl-mpm/mesh::get-node (cl-mpm:sim-mesh *sim*) (list 0 0 0))
;;      )))






(let ((mp-pos (list 567d0 340d0 0d0)))
  (defun plot-nonlocal-inter ()
    (vgplot:close-all-plots)
    (let* ((find-pos (cl-mpm/utils:vector-from-list mp-pos))
           (mp (cl-mpm/setup::find-mp *sim* find-pos)))
      (multiple-value-bind (pos weights) (cl-mpm/damage::get-nonlocal-interactions *sim* mp)
        (let ((x (loop for p in pos collect (cl-mpm/utils:varef p 0)))
              (y (loop for p in pos collect (cl-mpm/utils:varef p 1))))
          (vgplot:3d-plot x y weights ";;with points lc palette")
          (vgplot:xlabel "x")
          (vgplot:ylabel "y")))))

  (defun plot-damage-domain ()
    (when *sim*
      (let* ((find-pos (cl-mpm/utils:vector-from-list mp-pos))
             (found-mp (cl-mpm/setup::find-mp *sim* find-pos)))
        ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
        ;; (setf (cl-mpm/damage::sim-enable-stress-based-length *sim*) nil)
        ;; (setf (cl-mpm/damage::sim-enable-ekl *sim*) nil)
        ;; (cl-mpm::iterate-over-mps
        ;;  (cl-mpm::sim-mps *sim*)
        ;;  (lambda (mp)
        ;;    (setf (cl-mpm/particle::mp-local-length mp) (* 4d0 (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*))))))
        ;; (cl-mpm/damage::setup-mp-local-list *sim*)
        (multiple-value-bind (pos weights) (cl-mpm/damage::get-nonlocal-interactions *sim* found-mp))

        (cl-mpm/plotter:simple-plot
         *sim*
         :plot :deformed
         :trial nil
         :colour-func (lambda (mp) (if (eq mp found-mp) 1d0 0d0))))
      (sleep 0.5)
      (cl-mpm/plotter:simple-plot
       *sim*
       :plot :deformed
       :trial nil
       :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp)))
      (sleep 0.5)
      (cl-mpm/plotter:simple-plot
       *sim*
       :plot :deformed
       :trial nil
       :colour-func (lambda (mp) (cl-mpm/particle::mp-debug-j mp)))
      )))

(defun brittle-cracks ()
  (cl-mpm/utils::set-workers 12)
  (vgplot:close-all-plots)
  (dolist (initial-stress (list
                           ;; :CRYO-STATIC
                           ;; :ELASTIC-STATIC
                           :NIL
                           ))
    (dolist (height (list 150d0))
      (dolist (float (list 0d0))
        (dolist (notch-ratio (list 0d0))
          (dolist (friction (list 0d0))
            (let* ((mps 3)
                   ;; (name (format nil "height_~F_stress_~A_friction_~F_notch_~F_floatation_~F" height initial-stress friction notch-ratio float))
                   (name (format nil "height_~F" height))
                   (output-dir (format nil "./output-~A/" name))
                   (H height)
                   (ice-aspect 2d0)
                   (floatation-ratio float))

              (defparameter *length-scaler* 2d0)
              (defparameter *alpha* 0d0)
              ;; (defparameter *angle* 38d0)
              ;; (defparameter *angle-r* 10d0)
              ;; (defparameter *rt* 1d0)
              ;; (defparameter *rc* 0d0)
              (defparameter *ductility* 10d0)
              (format t "Running ~A~%" output-dir)
              (setup
               :refine 1
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
              (cl-mpm:iterate-over-mps
               (cl-mpm:sim-mps *sim*)
               (lambda (mp)
                 (change-class mp 'cl-mpm/particle::particle-ice-brittle)))
              (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) nil)
              (setf (cl-mpm/damage::sim-enable-stress-based-length *sim*) nil)
              (setf (cl-mpm/damage::sim-enable-ekl *sim*) nil)
              (setf
               (cl-mpm:sim-settings *sim*)
               (list :OCEAN-HEIGHT *water-height*
                     :OFFSET 2))
              (cl-mpm/output:add-node-output *sim* :SCALAR "damage" #'cl-mpm/mesh::node-damage)
              (cl-mpm/dynamic-relaxation::run-adaptive-load-control
               *sim*
               :criteria 1d-6
               :output-dir output-dir
               :max-adaptive-steps 0
               :load-steps 1
               :substeps 10
               :sub-conv-steps 500
               ;; :loading-function (lambda (p))
               :enable-damage t
               :enable-plastic t
               ;; :stagger-damage :MONOLITH-QS
               ;; :stagger-damage :FULL
               :stagger-damage :HYBRID-FULl
               :max-damage-inc 100d0
               :dt-scale 0.9d0
               :plotter
               (lambda (sim)
                 ;; (setf (cl-mpm::sim-enable-damage *sim*) t)
                 ;; (cl-mpm/damage:calculate-damage *sim* 1d-15)
                 ;; (setf (cl-mpm::sim-enable-damage *sim*) nil)
                 (plot-domain)
                 )
               :post-iter-step
               (lambda (i o e)
                 ))

              (vgplot:title output-dir)
              (vgplot:print-plot (merge-pathnames (format nil "frame_~A.png" name)) :terminal "png size 1920,1080")
              )))))))




(defun test-grads ()
  (cl-mpm::iterate-over-mps
   (cl-mpm::sim-mps *sim*)
   (lambda (mp)
     (setf (cl-mpm/particle::mp-damage mp)
           (cl-mpm/utils:varef (cl-mpm/particle::mp-position mp) 0))))
  (cl-mpm::iterate-over-mps
   (cl-mpm::sim-mps *sim*)
   (lambda (mp)
     (cl-mpm/damage::calculate-average-damage-grads (cl-mpm::sim-mesh *sim*) mp (cl-mpm/particle::mp-local-length mp))))
  (plot-domain))


(defun plot-shape-functions ()
  (let* ((R 0.1d0)
         (h 1d0)
         (x (loop for x from (* -2d0 h) to (* 2d0 h) by 0.01d0 collect x)))
    (vgplot:plot x (mapcar (lambda (x) (cl-mpm/shape-function::shape-gimp-fast x R h)) x) "fast"
                 x (mapcar (lambda (x) (cl-mpm/shape-function::shape-gimp-fbar x R h)) x) "fbar")))

(let* ((angle 48d0)
       (angle-r 30d0)
       (rt 1d0)
       (rc 0d0)
       (rs (est-shear-from-angle angle angle-r rc))
       (pd-oversize 1d-3)
       (E 1d9)
       (init-stress 100d3)
       (ductility 10d0)
       )
  (let* ((pd (- 1d0 pd-oversize))
         (k (cl-mpm/damage::find-k-damage E init-stress ductility pd))
         (ds (cl-mpm/damage::damage-response-exponential-peerlings-residual k E init-stress ductility rs))
         (dc (cl-mpm/damage::damage-response-exponential-peerlings-residual k E init-stress (* ductility 100d0) rc)))
    (format t "Damage pd ~E~%" pd)
    (format t "Damage t ~E~%" (cl-mpm/damage::damage-response-exponential-peerlings-residual k E init-stress ductility rt))
    (format t "Damage c ~E~%" dc)
    (format t "Damage s ~E~%" ds)
    (format t "Real residual angle ~E~%" (cl-mpm/utils::rad-to-deg (atan (* (/ (- 1d0 ds) (- 1d0 dc)) (tan (cl-mpm/utils::deg-to-rad angle))))))
    (format t "0 dc residual angle ~E~%" (cl-mpm/utils::rad-to-deg (atan (* (- 1d0 ds)  (tan (cl-mpm/utils::deg-to-rad angle))))))
    )
  )

(let* ((density 918d0)
       (water-density 1028d0)
       (h 900d0)
       ;; (f 0.76d0)
       (f 0.96d0)
       )
  (format t "Cliff ~F~%" (- h (* h (* f (/ density water-density))))))


(defun plot-regular-pressure ()
  (let ((x (loop for x from -2d0 to 2d0 by 0.01d0 collect x))
        (g -9.8d0)
        (rho 1d3)
        (datum 0d0)
        (l 0.5d0)
        )
    (vgplot:close-all-plots)
    (vgplot:plot
     ;; x (mapcar (lambda (x) (cl-mpm/buoyancy::pressure-at-depth x datum rho g)) x) "fast"
     ;; x (mapcar (lambda (x) (cl-mpm/buoyancy::pressure-at-depth-regular x datum rho g l)) x) "fast"

     ;; x (mapcar (lambda (x) (cl-mpm/utils::varef (cl-mpm/buoyancy::buoyancy-virtual-stress x datum rho g) 1)) x) "fast"
     ;; x (mapcar (lambda (x) (cl-mpm/utils::varef (cl-mpm/buoyancy::buoyancy-virtual-stress-regular x datum rho g l) 1)) x) "fast"
     x (mapcar (lambda (x) (cl-mpm/utils::varef (cl-mpm/buoyancy::buoyancy-virtual-div x datum rho g) 1)) x) "fast"
     x (mapcar (lambda (x) (cl-mpm/utils::varef (cl-mpm/buoyancy::buoyancy-virtual-div-regular x datum rho g l) 1)) x) "fast"
     )))


(defun plot-pq ()
  (let* ((mp-count (length (cl-mpm:sim-mps *sim*)))
         (p (make-array mp-count))
         (q (make-array mp-count)))
    (vgplot:close-all-plots)
    (let ((i 0 ))
      (cl-mpm::iterate-over-mps-serial
       (cl-mpm:sim-mps *sim*)
       (lambda (mp)
         (setf (aref p i) (* 1/3 (cl-mpm/utils::trace-voigt (cl-mpm/particle::mp-stress mp))))
         (setf (aref q i) (cl-mpm/fastmaths::voigt-j2 (cl-mpm/utils::deviatoric-voigt (cl-mpm/particle::mp-stress mp))))
         (incf i)))
      (vgplot:plot p q ";;with points"))))

;; (pprint
;;  (cl-mpm/buoyancy::buoyancy-virtual-div-regular
;;   -1d0
;;   0d0
;;   (* -1d0 (- 1000d0 900d0))
;;   (cl-mpm::sim-gravity *sim*)
;;   0d0))



(let ((x 10d0)
      (density 918d0)
      (water-density 1028d0)
      (g 9.81d0)
      (nu 0.3d0))
  (let* ((p-water (* water-density g x))
         (p-ice (* density g x))
         (e-elastic (* density g x))
         (k (/ nu (- 1d0 nu)))
         (p-elastic (* 1/3 (+ e-elastic  (* k e-elastic) (* k e-elastic))))
         )
    (pprint p-water)
    (pprint p-ice)
    (pprint p-elastic)
    (pprint (/ p-ice p-water))
    (pprint (/ p-elastic p-water))
    (pprint (/ p-elastic p-ice))

    )
  )

(let ((a :CRYO))
  (case a
    (:CRYO (print "Hello"))
    )
  )
(defun flotation-ratio (h w)
  (let* ((density 918d0)
         (water-density 1028d0)
         (f (/ density water-density)))
    (pprint (/ w h))
    (pprint (/ (/ w h) f))))
