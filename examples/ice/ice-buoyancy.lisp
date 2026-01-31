(defpackage :cl-mpm/examples/ice-buoyancy
  (:use :cl
   :cl-mpm/example))
(in-package :cl-mpm/examples/ice-buoyancy)

;; (sb-ext:restrict-compiler-policy 'speed  3 3)
;; (sb-ext:restrict-compiler-policy 'debug  0 0)
;; (sb-ext:restrict-compiler-policy 'safety 0 0)
(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)


(defclass cl-mpm/particle::particle-ice-erodable (cl-mpm/particle::particle-ice-delayed
                                                    cl-mpm/particle::particle-erosion)
  ())


(defmethod cl-mpm/erosion::mp-erosion-enhancment ((mp cl-mpm/particle::particle-ice-erodable))
  (+ 1d0 (* 10 (cl-mpm/particle::mp-damage mp)))
  ;; (+ 1d0 (* 10 (cl-mpm/particle::mp-strain-plastic-vm mp)))
  )

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-elastic) dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-damaged mesh mp dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-dynamic-relaxation mesh mp dt fbar)
  (cl-mpm::update-stress-kirchoff-dynamic-relaxation mesh mp dt fbar)
  )
(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-elastic) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  ;; (cl-mpm::update-domain-det mesh mp dt)
  ;; (cl-mpm::co-domain-corner-2d mesh mp dt)
  (cl-mpm::update-domain-polar-2d mesh mp dt)
  ;; (cl-mpm::update-domain-midpoint mesh mp dt)
  ;; (cl-mpm::update-domain-deformation mesh mp dt)
  ;; (cl-mpm::scale-domain-size mesh mp)
  )


(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-ice-brittle) dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-damaged mesh mp dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-dynamic-relaxation mesh mp dt fbar)
  (cl-mpm::update-stress-kirchoff-dynamic-relaxation mesh mp dt fbar)
  )

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-ice-brittle) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  ;; (cl-mpm::update-domain-det mesh mp dt)
  ;; (cl-mpm::co-domain-corner-2d mesh mp dt)
  (cl-mpm::update-domain-polar-2d mesh mp dt)
  ;; (cl-mpm::update-domain-midpoint mesh mp dt)
  ;; (cl-mpm::update-domain-deformation mesh mp dt)
  ;; (cl-mpm::scale-domain-size mesh mp)
  )

(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-ice-brittle) dt)
  (with-accessors (;(stress cl-mpm/particle::mp-undamaged-stress)
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
                   (pd-inc cl-mpm/particle::mp-plastic-damage-evolution))
      mp
    (progn
      (let ((ps-y (sqrt (* E (expt ps-vm 2))))
            (stress (cl-mpm/constitutive:linear-elastic-mat strain de))
            )
        (setf y
              (*
               ;; (- 1d0 damage)
               (+
                ;; ps-y
                ;; (cl-mpm/damage::tensile-energy-norm
                ;;  strain
                ;;  E
                ;;  de)
                ;; (cl-mpm/damage::tensile-energy-norm-pressure
                ;;  strain
                ;;  E
                ;;  nu
                ;;  de
                ;;  (* 1d0 (* 1d0 (magicl:det def)) (/ (- pressure) 3)))
                ;; (cl-mpm/damage::criterion-mohr-coloumb-rankine-stress-tensile
                ;;  ;; cl-mpm/damage::criterion-mohr-coloumb-stress-tensile
                ;;  (cl-mpm/fastmaths:fast-.+
                ;;   stress
                ;;   (cl-mpm/utils:voigt-eye (* 0d0
                ;;                              (* 1d0 (magicl:det def))
                ;;                              (/ (- pressure) 3))))
                ;;  (* angle (/ pi 180d0)))
                (if pd-inc ps-y 0d0)
                (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile
                 (cl-mpm/fastmaths:fast-.+
                  (cl-mpm/fastmaths::fast-scale-voigt stress
                                                      (/ 1d0 (magicl:det def)))
                  (cl-mpm/utils:voigt-eye (* 1d0
                                             ;; (* 1d0 (magicl:det def))
                                             (/ (- pressure) 3))))
                 (* angle (/ pi 180d0)))
                ))))
      )))
(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-elastic-damage) dt)
  (with-accessors ((strain cl-mpm/particle::mp-strain)
                   (stress cl-mpm/particle::mp-undamaged-stress)
                   (E cl-mpm/particle::mp-e)
                   (nu cl-mpm/particle::mp-nu)
                   (y cl-mpm/particle::mp-damage-y-local)
                   (de cl-mpm/particle::mp-elastic-matrix)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (pressure cl-mpm/particle::mp-pressure)
                   ) mp
    (let ((angle 60d0))
      (setf y
            (+
             (cl-mpm/damage::criterion-mohr-coloumb-rankine-stress-tensile
              stress
              ;; (cl-mpm/fastmaths:fast-.+
              ;;  (cl-mpm/utils:voigt-eye (* 1d0 (magicl:det def) (/ (- pressure) 3)))
              ;;  )
              (* angle (/ pi 180d0)))
             )))
    ;; (setf y (cl-mpm/damage::tensile-energy-norm-pressure strain E nu de
    ;;                                                      0d0
    ;;                                                      ;; (* 0.3d0 ;(magicl:det def)
    ;;                                                      ;;    pressure)
    ;;                                                      ))
    ))


(declaim (notinline plot-domain))
(defun plot-domain (&key (trial t))
  (when *sim*
    (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
           (h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
           (ms-x (first ms))
           (ms-y (second ms)))
      (vgplot:format-plot t "set object 1 rect from 0,0 to ~f,~f fc rgb 'blue' fs transparent solid 0.5 noborder behind" ms-x *water-height*)
      (vgplot:format-plot t "set object 2 rect from 0,0 to ~f,~f fc rgb 'black' fs transparent solid 1 noborder behind" ms-x *offset*)
      )
    ;; (cl-mpm::g2p (cl-mpm:sim-mesh *sim*)
    ;;              (cl-mpm:sim-mps *sim*)
    ;;              (cl-mpm:sim-dt *sim*)
    ;;              0d0
    ;;              :TRIAL)
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

(defparameter *rs* 0.3d0)
(defparameter *enable-plastic-damage* t)

(defun setup (&key (refine 1) (mps 2)
                (pressure-condition t)
                (cryo-static t)
                (hydro-static nil)
                (melange nil)
                (friction 0d0)
                (ice-height 800d0)
                (bench-length 0d0)
                (aspect 1)
                (floatation-ratio 0.9)
                (slope 0.1d0)
                )
  (let* ((density 918d0)
         (water-density 1028d0)
         ;; (density 900d0)
         ;; (water-density 1000d0)
         (mesh-resolution (/ 10d0 refine))
         ;; (refines 0)
         ;; (mps (* mps (expt 2 (- refines 1))))
         ;; (h-fine (/ mesh-resolution (expt 2 (- refines 1))))
         (h-fine mesh-resolution)
         (offset (* mesh-resolution 2))
         ;; (end-height ice-height)
         ;; (start-height ice-height)

         (end-height ice-height)
         (ice-length (* ice-height aspect))
         (start-height (+ ice-height (* slope ice-length)))
         (ice-height end-height)
         (floating-point (* ice-height (/ density water-density)))
         (water-level (* floating-point floatation-ratio))
         (datum (+ water-level offset))
         ;; (datum (* (round datum mesh-resolution) mesh-resolution))
         (domain-size (list (+ ice-length (* 3 ice-height)) (* start-height 2)))
         (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
         (block-size (list ice-length (max start-height end-height)))
         (E 1d9)
         )
    (defparameter *water-height* datum)
    (defparameter *offset* offset)
    (defparameter *ice-length* ice-length)
    (setf *sim* (cl-mpm/setup::make-simple-sim mesh-resolution element-count
                                               :sim-type
                                               ;; 'cl-mpm/damage::mpm-sim-usl-damage
                                               ;; 'cl-mpm/damage::mpm-sim-agg-damage
                                               ;; 'cl-mpm::mpm-sim-usf
                                               'cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul
                                               ;; 'cl-mpm/dynamic-relaxation::mpm-sim-dr-multigrid
                                               ;; 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul-usl
                                               ;; 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic
                                               :args-list
                                               (list
                                                :enable-fbar t
                                                :enable-aggregate t
                                                :split-factor (* 1.2d0 (sqrt 2) 0.33d0)
                                                :max-split-depth 4
                                                ;; :refinement refines
                                                )))
    (let* ((angle 40d0)
           (init-stress (* 0.1185d6 1d0))
           (init-c (cl-mpm/damage::mohr-coloumb-tensile-to-coheasion init-stress (* angle (/ pi 180))))
           (gf 10d0)
           ;; (length-scale (* h-fine 2d0))
           (length-scale (* mesh-resolution 2d0))
           (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress E))
           ;; (ductility (/ (/ 5d0 0.125d0) refine))
           (ductility 20d0)
           ;; (ductility 2d0)
           (oversize (cl-mpm/damage::compute-oversize-factor (- 1d0 1d-3) ductility)))
      (format t "Ice length ~F~%" ice-length)
      (format t "Water height ~F~%" water-level)
      (format t "True Water height ~F~%" (- datum offset))
      (format t "Cliff height ~F~%" (- (+ offset ice-height) datum))
      (format t "Mesh size ~F~%" mesh-resolution)
      (format t "Estimated oversize ~F~%" oversize)
      (format t "Estimated lc ~E~%" length-scale)
      (format t "Estimated ductility ~E~%" ductility)
      (format t "Init stress ~E~%" init-stress)
      (format t "Init c ~E~%" init-c)
      (let ((rs *rs*))
        (cl-mpm:add-mps
         *sim*
         (cl-mpm/setup:make-block-mps
          (list 0 offset)
          block-size
          (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
          density
          'cl-mpm/particle::particle-ice-erodable
          ;; 'cl-mpm/particle::particle-ice-delayed
          :E 1d9
          :nu 0.24d0

          :kt-res-ratio 1d0
          :kc-res-ratio 0d0
          :g-res-ratio rs

          :initiation-stress init-stress;18d3
          :friction-angle angle
          :psi (* 0d0 (/ pi 180))
          :phi (* angle (/ pi 180))
          :c (* init-c oversize)
          :softening 0d0
          :ductility ductility
          :local-length length-scale
          :delay-time 1d4
          :delay-exponent 2
          :enable-plasticity t
          :enable-damage t
          :plastic-damage-evolution *enable-plastic-damage*


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
        (est-angle angle rs)
        )
      (when melange
        (let* ((melange-depth 50d0)
               (melange-pressure -1000d3))
          (cl-mpm::add-bcs-force-list
           *sim*
           (cl-mpm/buoyancy::make-bc-pressure
            *sim*
            melange-pressure
            0d0
            :clip-func
            (lambda (pos)
              (and
               (<= (cl-mpm/utils:varef pos 1) datum)
               (>= (cl-mpm/utils:varef pos 1) (- datum melange-depth)))))))
        ;; (let* ((melange-depth 50d0)
        ;;        (melange-length (- (first domain-size) (first block-size)))
        ;;        (block-size (list melange-length melange-depth))
        ;;        (offset (- datum (* melange-depth (/ density water-density))))
        ;;        )
        ;;   (cl-mpm:add-mps
        ;;    *sim*
        ;;    (cl-mpm/setup:make-block-mps
        ;;     (list ice-length offset)
        ;;     block-size
        ;;     (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
        ;;     density
        ;;     'cl-mpm/particle::particle-ice-delayed
        ;;     :E 1d9
        ;;     :nu 0.24d0

        ;;     :kt-res-ratio 1d0
        ;;     :kc-res-ratio 0d0
        ;;     :g-res-ratio 0.95d0

        ;;     :index 2
        ;;     :initiation-stress init-stress;18d3
        ;;     :friction-angle angle
        ;;     :psi (* 0d0 (/ pi 180))
        ;;     :phi (* angle (/ pi 180))
        ;;     :c (* init-c oversize)
        ;;     :softening 0d0
        ;;     :ductility ductility
        ;;     :local-length length-scale
        ;;     :delay-time 1d5
        ;;     :delay-exponent 2
        ;;     :enable-plasticity t
        ;;     :enable-damage t
        ;;     ))
        ;;   (cl-mpm::iterate-over-mps
        ;;    (cl-mpm::sim-mps *sim*)
        ;;    (lambda (mp)
        ;;      (when (= (cl-mpm/particle::mp-index mp) 2)
        ;;        (cl-mpm/damage::set-mp-damage mp 0.99d0))))
        ;;   )
        )

      ;; (cl-mpm/setup::remove-sdf *sim*
      ;;                           (lambda (p)
      ;;                             (cl-mpm/setup::plane-point-point-sdf
      ;;                              p
      ;;                              (cl-mpm/utils:vector-from-list (list 0d0 (+ offset start-height) 0d0))
      ;;                              (cl-mpm/utils:vector-from-list (list ice-length (+ offset end-height) 0d0))))
      ;;                           :refine 1
      ;;                           )
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
         :k-z 1d0))
      (unless (= start-height end-height)
        (cl-mpm/setup::remove-sdf *sim*
                                  (lambda (p)
                                    (cl-mpm/setup::plane-point-point-sdf
                                     p
                                     (cl-mpm/utils:vector-from-list (list 0d0 (+ offset start-height) 0d0))
                                     (cl-mpm/utils:vector-from-list (list ice-length (+ offset end-height) 0d0))))
                                  :refine 2))


      (let ((cutout (+ (- ice-height water-level) 50d0))
            (cutback bench-length))
        (when (> cutback 0d0)
          (cl-mpm/setup:remove-sdf
           *sim*
           (cl-mpm/setup::rectangle-sdf (list (first block-size) (+ offset ice-height))
                                        (list cutback cutout))
           )
          ;; (cl-mpm/setup::remove-sdf *sim*
          ;;                           (lambda (p)
          ;;                             (cl-mpm/setup::plane-point-point-sdf
          ;;                              p
          ;;                              (cl-mpm/utils:vector-from-list (list ice-length datum 0d0))
          ;;                              (cl-mpm/utils:vector-from-list (list (- ice-length cutback) offset 0d0))))
          ;;                           :refine 3
          ;;                           )
          )
        ))
    (setf
     (cl-mpm:sim-bcs *sim*)
     (cl-mpm/bc::make-outside-bc-varfix
      (cl-mpm:sim-mesh *sim*)
      '(0 nil 0)
      '(0 nil 0)
      '(nil 0 nil)
      '(nil 0 nil)
      '(nil nil 0)
      '(nil nil 0)))
    (cl-mpm/setup:setup-bcs
     *sim*
     :bottom '(0 0 0))

    (setf (cl-mpm:sim-mass-scale *sim*) 1d0)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 0.1d0
             (sqrt 1d0)
             (cl-mpm/setup:estimate-critical-damping *sim*)))
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)
    (when (typep *sim* 'cl-mpm/damage::mpm-sim-damage)
      (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t))
    (setf (cl-mpm::sim-allow-mp-split *sim*) t)
    (setf (cl-mpm::sim-ghost-factor *sim*)
          nil
          ;; (* E 1d-4)
          )
    ;; (setf (cl-mpm::sim-velocity-algorithm *sim*) :PIC)
    (setf (cl-mpm::sim-velocity-algorithm *sim*) :FLIP)
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
         (* E 0.01d0)
         friction
         ;; 0.1d0
         0d0
         )))

    (defparameter *bc-erode*
      (cl-mpm/erosion::make-bc-erode
       *sim*
       :enable nil
       :rate 1d-4
       :scalar-func (lambda (pos)
                      1d0
                      ;; (min 1d0 (exp (* 0.5d0 (- (cl-mpm/utils:varef pos 1) datum))))
                      )
       :clip-func (lambda (pos)
                    (and
                     (>= datum (cl-mpm/utils:varef pos 1))
                     (>= (cl-mpm/utils:varef pos 1) (+ offset mesh-resolution) )))))
    (when (> offset 0d0)
      (cl-mpm:add-bcs-force-list
       *sim*
       *floor-bc*
       ))
    (cl-mpm:add-bcs-force-list
     *sim*
     *bc-erode*
     )
    (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
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


(defun calving-test ()
  (loop for dt in (list 1d3)
        do
           (let* ((mps 2)
                  (H 600d0))
             (setup :refine 0.25
                    :friction 0.5d0
                    :bench-length (* 0d0 H)
                    :ice-height H
                    :mps mps
                    :hydro-static nil
                    :cryo-static t
                    :melange nil
                    :aspect 2d0
                    :slope 0d0
                    :floatation-ratio 0.7d0)
             (plot-domain)
             (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0)
             (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
             (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) nil
                   ;; (cl-mpm::sim-ghost-factor *sim*) (* 1d9 1d-3)
                   (cl-mpm::sim-ghost-factor *sim*) nil
                   ;; (cl-mpm/dynamic-relaxation::sim-convergence-critera *sim*) 1d-3
                   )

             (setf lparallel:*debug-tasks-p* nil)
             (setf (cl-mpm::sim-allow-mp-damage-removal *sim*) nil)
             (setf (cl-mpm::sim-mp-damage-removal-instant *sim*) nil)
             (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
             (setf (cl-mpm:sim-enable-damage *sim*) nil)
             (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
             (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
             (let ((step 0))
               (cl-mpm/dynamic-relaxation::run-multi-stage
                *sim*
                :output-dir "./output/"
                :dt dt
                :dt-scale 1d0
                :damping-factor 1d0;(float (sqrt 2) 0d0)
                :conv-criteria 1d-3
                :conv-load-steps 1
                ;; :min-adaptive-steps -4
                ;; :max-adaptive-steps 10
                :min-adaptive-steps -8
                :max-adaptive-steps 14
                :substeps 20
                :total-time 1d7
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
                ;; :explicit-dt-scale 0.25d0
                ;; :explicit-damping-factor 1d-3
                ;; :explicit-dynamic-solver 'cl-mpm/damage::mpm-sim-agg-damage
                :explicit-damping-factor 0d-4
                :explicit-dt-scale 1d0
                :explicit-dynamic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic
                :post-conv-step (lambda (sim)
                                  (setf (cl-mpm/buoyancy::bc-enable *bc-erode*) nil))
                :setup-quasi-static
                (lambda (sim)
                  (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
                  (setf
                   (cl-mpm/aggregate::sim-enable-aggregate sim) t
                   (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC
                   (cl-mpm::sim-ghost-factor sim) nil
                   ;; (cl-mpm::sim-ghost-factor sim) (* 1d9 1d-3)
                   (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0))
                :setup-dynamic
                (lambda (sim)
                  (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
                  (setf (cl-mpm/aggregate::sim-enable-aggregate sim) t
                        (cl-mpm::sim-velocity-algorithm sim) :BLEND
                        (cl-mpm::sim-ghost-factor sim) nil;(* 1d9 1d-4)
                        (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 2d0
                        ))
                ))
             )))






(defun save-test-vtks (&key (output-dir "./output/"))
  (cl-mpm::finalise-loadstep *sim*)
  (cl-mpm/output:save-vtk (merge-pathnames "test.vtk" output-dir) *sim*)
  (cl-mpm/output:save-vtk-nodes (merge-pathnames "test_nodes.vtk" output-dir) *sim*)
  (cl-mpm/output:save-vtk-cells (merge-pathnames "test_cells.vtk" output-dir) *sim*)
  )



(defun est-angle (angle rs)
  (let* ((rc 0d0)
         (ratio (/ (- 1d0 rs) (- 1d0 rc)))
         (angle-plastic (* angle (/ pi 180)))
         (angle-plastic-damaged (atan (* ratio (tan angle-plastic))))
         )
    (format t "Plastic virgin angle: ~F~%"
            (* (/ 180 pi) angle-plastic))
    (format t "Plastic residual angle: ~F~%"
            (* (/ 180 pi) angle-plastic-damaged))))
