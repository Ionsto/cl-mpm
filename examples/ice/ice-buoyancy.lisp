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


(defmethod cl-mpm::update-dynamic-stats (sim)
  (with-accessors ((stats-energy cl-mpm::sim-stats-energy)
                   (stats-oobf cl-mpm::sim-stats-oobf)
                   (stats-power cl-mpm::sim-stats-power))
      sim
      (setf stats-energy (cl-mpm/dynamic-relaxation:estimate-energy-norm sim)
            stats-oobf (cl-mpm/dynamic-relaxation:estimate-oobf sim)
            stats-power (cl-mpm/dynamic-relaxation:estimate-power-norm sim))))

(defmethod cl-mpm/erosion::mp-erosion-enhancment ((mp cl-mpm/particle::particle-ice-erodable))
  (+ 1d0 (* 10 (cl-mpm/particle::mp-damage mp)))
  ;; (+ 1d0 (* 10 (cl-mpm/particle::mp-strain-plastic-vm mp)))
  )


(defun get-damage ()
  (lparallel:pmap-reduce
   (lambda (mp)
     (if (typep mp 'cl-mpm/particle::particle-damage)
         (*
          1d2
          (cl-mpm/particle::mp-mass mp)
          (cl-mpm/particle:mp-damage mp))
         0d0))
   #'+
   (cl-mpm:sim-mps *sim*)))
(defun get-plastic ()
  (lparallel:pmap-reduce
   (lambda (mp)
     (if (typep mp 'cl-mpm/particle::particle-plastic)
         (*
          (cl-mpm/particle::mp-mass mp)
          (cl-mpm/particle::mp-strain-plastic-vm mp))
         0d0))
   #'+
   (cl-mpm:sim-mps *sim*)))

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
  (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
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
                   )
      mp
    (progn
      (let ((ps-y (sqrt (* E (expt ps-vm 2)))))
        (setf y
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
               ;;  (* 0d0 -1d0 (* 1d0 (magicl:det def)) (/ (- pressure) 1)))
               ;; (cl-mpm/damage::criterion-mohr-coloumb-rankine-stress-tensile
               ;;  ;; cl-mpm/damage::criterion-mohr-coloumb-stress-tensile
               ;;  (cl-mpm/fastmaths:fast-.+
               ;;   stress
               ;;   (cl-mpm/utils:voigt-eye (* 0d0
               ;;                              (* 1d0 (magicl:det def))
               ;;                              (/ (- pressure) 3))))
               ;;  (* angle (/ pi 180d0)))
               (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile
                (cl-mpm/fastmaths:fast-.+
                 (cl-mpm/fastmaths::fast-scale-voigt stress (/ 1d0 (magicl:det def)))
                 (cl-mpm/utils:voigt-eye (* 1d0
                                            ;; (* 1d0 (magicl:det def))
                                            (/ (- pressure) 3))))
                (* angle (/ pi 180d0)))
               )))
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
(defun plot-domain ()
  (when *sim*
    (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
           (h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
           (ms-x (first ms))
           (ms-y (second ms)))
      (vgplot:format-plot t "set object 1 rect from 0,0 to ~f,~f fc rgb 'blue' fs transparent solid 0.5 noborder behind" ms-x *water-height*))
    ;; (cl-mpm::g2p (cl-mpm:sim-mesh *sim*)
    ;;              (cl-mpm:sim-mps *sim*)
    ;;              (cl-mpm:sim-dt *sim*)
    ;;              0d0
    ;;              :TRIAL)
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :trial t
     ;; :colour-func (lambda (mp) (sqrt (cl-mpm/constitutive::voigt-j2 (cl-mpm/utils:deviatoric-voigt (cl-mpm/particle:mp-stress mp)))))
     ;; :colour-func #'cl-mpm/particle::mp-strain-plastic-vm
     ;; :colour-func (lambda (mp) (/ (cl-mpm/particle:mp-mass mp)
     ;;                              (cl-mpm/particle:mp-volume mp)))
     :colour-func #'cl-mpm/particle::mp-damage
     ;; :colour-func (lambda (mp) (cl-mpm/utils::varef  (cl-mpm/particle::mp-stress mp ) 1))
     ))
  )

(defun setup (&key (refine 1) (mps 2)
                (pressure-condition t)
                (cryo-static t)
                (hydro-static nil)
                (melange t)
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
         (domain-size (list (+ ice-length (* 1 ice-height)) (* start-height 2)))
         (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
         (block-size (list ice-length (max start-height end-height)))
         (E 1d9)
         )
    (defparameter *water-height* datum)
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
                                                ;; :refinement refines
                                                )))
    (let* ((angle 40d0)
           (init-stress (* 0.1185d6 1d0))
           (init-c (cl-mpm/damage::mohr-coloumb-tensile-to-coheasion init-stress (* angle (/ pi 180))))
           (gf 1000d0)
           (length-scale (* h-fine 1d0))
           (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress E))
           (ductility 5d0)
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
      (cl-mpm:add-mps
       *sim*
       (cl-mpm/setup:make-block-mps
        (list 0 offset)
        block-size
        (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
        density
        ;'cl-mpm/particle::particle-ice-erodable
        'cl-mpm/particle::particle-ice-delayed
        :E 1d9
        :nu 0.24d0

        :kt-res-ratio 1d0
        :kc-res-ratio 0d0
        :g-res-ratio 0.95d0

        :initiation-stress init-stress;18d3
        :friction-angle angle
        :psi (* 0d0 (/ pi 180))
        :phi (* angle (/ pi 180))
        :c (* init-c oversize)
        :softening 0d0
        :ductility ductility
        :local-length length-scale
        :delay-time 1d4
        :delay-exponent 1
        :enable-plasticity t
        :enable-damage t


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
         :k-z 1d0)
        ;; (break)
        ;; (cl-mpm/setup::initialise-stress-self-weight *sim* (+ offset start-height))
        )
      ;; (break)
      ;; (cl-mpm/setup::initialise-stress-self-weight *sim* (+ offset ice-height))
      (unless (= start-height end-height)
        (cl-mpm/setup::remove-sdf *sim*
                                  (lambda (p)
                                    (cl-mpm/setup::plane-point-point-sdf
                                     p
                                     (cl-mpm/utils:vector-from-list (list 0d0 (+ offset start-height) 0d0))
                                     (cl-mpm/utils:vector-from-list (list ice-length (+ offset end-height) 0d0))))
                                  :refine 2
                                  ))


      (let ((cutout (+ (- ice-height water-level) 00d0))
            (cutback bench-length)
            )
        ;; (pprint cutout)
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

    (setf (cl-mpm:sim-mass-scale *sim*) 1d0)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 0.1d0
             (sqrt 1d0)
             (cl-mpm/setup:estimate-critical-damping *sim*)))
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)
    (when (typep *sim* 'cl-mpm/damage::mpm-sim-damage)
      (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t))
    (setf (cl-mpm::sim-allow-mp-split *sim*) t)
    (setf (cl-mpm::sim-max-split-depth *sim*) 3)
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
       :rate 1d0
       :scalar-func (lambda (pos)
                      (min 1d0 (exp (* 0.5d0 (- (cl-mpm/utils:varef pos 1) datum)))))
       :clip-func (lambda (pos)
                    (>= datum (cl-mpm/utils:varef pos 1))
                    ;; (and (>= (+ offset mesh-resolution) (cl-mpm/utils:varef pos 1)))
                    )))
    (when (> offset 0d0)
      (cl-mpm:add-bcs-force-list
       *sim*
       *floor-bc*
       ))
    ;; (cl-mpm:add-bcs-force-list
    ;;  *sim*
    ;;  *bc-erode*
    ;;  )
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

(defun run-static (&key (output-dir "./output/"))
  (uiop:ensure-all-directories-exist (list (uiop:merge-pathnames* output-dir)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  ;; (cl-mpm/dynamic-relaxation:converge-quasi-static
  ;;  *sim*
  ;;  :energy-crit 1d-2
  ;;  :oobf-crit 1d-2
  ;;  :substeps 50
  ;;  :conv-steps 100
  ;;  :dt-scale 0.5d0)
  ;; (setf (cl-mpm::sim-enable-damage *sim*) t)
  (setf (cl-mpm:sim-mass-scale *sim*) 1d0)
  (setf (cl-mpm:sim-damping-factor *sim*)
        (* 0.8d0
           (sqrt (cl-mpm:sim-mass-scale *sim*))
           (cl-mpm/setup:estimate-critical-damping *sim*)))

  (setf (cl-mpm/buoyancy::bc-enable *bc-erode*) t)
  (setf (cl-mpm:sim-dt *sim*)
        (* 0.5d0 (cl-mpm/setup:estimate-elastic-dt *sim*)))
  (let* ((dt-scale .5d0)
         (target-time 1d0)
         (substeps (ceiling target-time (cl-mpm:sim-dt *sim*)))
         (work 0d0)
         (oobf 0d0)
         (energy 0d0)
         (step-count 100)
         (load -9.8d0)
         (elastic-load 0d0)
         (load-0 (* elastic-load load))
         (load-inc (* (- 1d0 elastic-load) (/ load step-count))))
    (format t "Substeps ~D~%" substeps)
    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (setf (cl-mpm/particle::mp-enable-plasticity mp) t)))
    (loop for step from 0 below step-count
          while *run-sim*
          do
             (format t "Step ~D~%" step)
             ;; (setf work 0d0)
             (cl-mpm::iterate-over-mps
              (cl-mpm:sim-mps *sim*)
              (lambda (mp)
                (setf (cl-mpm/particle::mp-gravity mp) (+ load-0 (* step load-inc)))))

             (cl-mpm/dynamic-relaxation:converge-quasi-static
              *sim*
              :energy-crit 1d-1
              :oobf-crit 1d-1
              :substeps 50
              :conv-steps 1000
              :dt-scale dt-scale
              :post-iter-step
              (lambda (&rest args)
                (plot-domain)))

             ;; (let ((damage-0 0d0)
             ;;       (damage-1 0d0))
             ;;   (loop for i from 0 to 10
             ;;         while (or (= i 0) (> damage-1 (* (+ 1d0 1d-6) damage-0)))
             ;;         do
             ;;            (progn
             ;;              (setf damage-0 damage-1)
             ;;              (format t "Staggered solve: ~D~%" i)
             ;;              (cl-mpm/dynamic-relaxation:converge-quasi-static
             ;;               *sim*
             ;;               :energy-crit 1d-2
             ;;               :oobf-crit 1d-2
             ;;               :substeps 50
             ;;               :conv-steps 1000
             ;;               :dt-scale dt-scale
             ;;               :post-iter-step
             ;;               (lambda (&rest args)
             ;;                 (plot-domain)))
             ;;              (setf (cl-mpm/damage::sim-damage-delocal-counter-max *sim*) 1)
             ;;              (setf (cl-mpm::sim-enable-damage *sim*) t)
             ;;              (cl-mpm/damage::calculate-damage *sim*)
             ;;              (setf (cl-mpm::sim-enable-damage *sim*) nil)
             ;;              (setf damage-1 (get-damage))
             ;;              (format t "Staggered damage diff: ~E~%" (- damage-1 damage-0))
             ;;              ))
             ;;   (when (> damage-1 (* damage-0 1.1d0))
             ;;     (format t "Staggered error~%"))
             ;;   )
             ;; (cl-mpm/dynamic-relaxation:converge-quasi-static
             ;;  *sim*
             ;;  :oobf-crit 1d-1
             ;;  :energy-crit 1d-1
             ;;  :dt-scale 0.5d0
             ;;  ;; :post-iter-step
             ;;  ;; (lambda (&rest args)
             ;;  ;;   (plot-domain))
             ;;  )


             (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir (format nil "sim_~5,'0d.vtk" step)) *sim* )
             (cl-mpm/output:save-vtk-nodes (uiop:merge-pathnames* output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) *sim* )
             (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim_p_~5,'0d.vtk" step)) *sim* )
             (plot-domain)
             (swank.live:update-swank)
             (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step))
                                :terminal "png size 1920,1080"
                                )
             )))

(defun relax-elastic (output-dir dt-scale step crit-damp)
  (setf (cl-mpm::sim-enable-damage *sim*) nil)
  (setf (cl-mpm:sim-mass-scale *sim*) 1d0)
  (cl-mpm::iterate-over-mps
   (cl-mpm:sim-mps *sim*)
   (lambda (mp)
     (setf (cl-mpm/particle::mp-enable-plasticity mp) nil)))
  (let ((damping-0 (cl-mpm::sim-damping-factor *sim*)))
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 0.01d0
             (sqrt (cl-mpm:sim-mass-scale *sim*))
             crit-damp))
    (cl-mpm/dynamic-relaxation:converge-quasi-static
     *sim*
     :oobf-crit 1d-2
     :energy-crit 1d-2
     :dt-scale dt-scale
     :substeps 50
     :conv-steps 1000
     :pic-update t
     :post-iter-step
     (lambda (i oobf energy)
       (format t "Sub converge ~D~%" i)
       ;; (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir (format nil "sim_subconv_~5,'0d_~5,'0d.vtk" step i)) *sim*)
       ;; (cl-mpm/output::save-vtk-nodes (uiop:merge-pathnames* output-dir (format nil "sim_subconv_nodes_~5,'0d_~5,'0d.vtk" step i)) *sim*)
       ;; (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir (format nil "sim_conv_~5,'0d.vtk" i)) *sim*)
       ;; (cl-mpm/output::save-vtk-nodes (uiop:merge-pathnames* output-dir (format nil "sim_conv_nodes_~5,'0d.vtk" i)) *sim*)
       ;; (plot-domain)
       ))
    (setf (cl-mpm:sim-damping-factor *sim*) damping-0)
    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (setf (cl-mpm/particle::mp-enable-plasticity mp) t))))
  (setf (cl-mpm::sim-enable-damage *sim*) t)

  )

(defun run (&key (output-dir "./output/"))
  (uiop:ensure-all-directories-exist (list (uiop:merge-pathnames* output-dir)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (with-open-file (stream (merge-pathnames output-dir "./timesteps.csv") :direction :output :if-exists :supersede)
    (format stream "steps,time,damage,plastic,energy,oobf,step-type,mass~%"))


  (defparameter *damping-water* 0d0)
  (let ((dt-scale 0.5d0)
        (visc-damping (cl-mpm/buoyancy::bc-viscous-damping *water-bc*))
        (fbar (cl-mpm:sim-enable-fbar *sim*))
        )
    (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0)
    (setf (cl-mpm:sim-mass-scale *sim*) 1d0)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 1d-2
             (sqrt (cl-mpm:sim-mass-scale *sim*))
             (cl-mpm/setup:estimate-critical-damping *sim*)))
    (setf (cl-mpm:sim-enable-fbar *sim*) nil)
    (cl-mpm/dynamic-relaxation:converge-quasi-static
     *sim*
     :oobf-crit 1d-3
     :energy-crit 1d-3
     :dt-scale dt-scale
     :substeps 50
     :conv-steps 1000
     :post-iter-step
     (lambda (i oobf energy)
       (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir (format nil "sim_conv_~5,'0d.vtk" i)) *sim*)
       (cl-mpm/output::save-vtk-nodes (uiop:merge-pathnames* output-dir (format nil "sim_conv_nodes_~5,'0d.vtk" i)) *sim*)
       (cl-mpm/output::save-vtk-cells (uiop:merge-pathnames* output-dir (format nil "sim_conv_cells_~5,'0d.vtk" i)) *sim*)
       (plot-domain)))
    (setf (cl-mpm:sim-enable-fbar *sim*) fbar)
    (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) visc-damping)
    )

  (cl-mpm::iterate-over-mps
   (cl-mpm:sim-mps *sim*)
   (lambda (mp)
     (setf (cl-mpm/particle::mp-enable-plasticity mp) t
           (cl-mpm/particle::mp-enable-damage mp) t)))

  (setf (cl-mpm::sim-enable-damage *sim*) t)
  (setf (cl-mpm:sim-mass-scale *sim*) 1d0)
  (setf (cl-mpm/buoyancy::bc-enable *bc-erode*) nil)
  (let* ((dt-scale 0.90d0)
         (substeps 0d0)
         (work 0d0)
         (oobf 0d0)
         (energy 0d0)
         (sim-state :accelerate)
         (accelerate-target-time 1d1)
         (accelerate-mass-scale 1d4)
         (collapse-target-time 1d0)
         (collapse-mass-scale 1d0)
         (criteria-energy 5d-2)
         (criteria-oobf 5d-2)
         (criteria-hist 1.5d0)
         (target-time 1d0)
         (time 0d0)
         (damage-est 0d0)
         (dt-0 (/ (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale)
                  (sqrt (cl-mpm:sim-mass-scale *sim*))))
         (crit-damp (cl-mpm/setup:estimate-critical-damping *sim*))
         )
    (cl-mpm/output::save-simulation-parameters (merge-pathnames output-dir "settings.json")
                                               *sim*
                                               (list :dt target-time
                                                     :criteria-energy criteria-energy
                                                     :criteria-oobf criteria-oobf
                                                     :criteria-hist criteria-hist
                                                     :ocean-height *water-height*
                                                     :ice-length *ice-length*
                                                     :delay-time delay-time
                                                     ))
    (setf (cl-mpm::sim-mass-scale *sim*) accelerate-mass-scale
          target-time accelerate-target-time)
    (setf (cl-mpm:sim-dt *sim*) (* dt-scale (cl-mpm/setup:estimate-elastic-dt *sim*)))
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 1d-3
             (sqrt (cl-mpm:sim-mass-scale *sim*))
             (cl-mpm/setup:estimate-critical-damping *sim*)))
    (cl-mpm:iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (cl-mpm/fastmaths::fast-zero (cl-mpm/particle:mp-velocity mp))
       ))
    (setf substeps (ceiling target-time (cl-mpm:sim-dt *sim*)))
    (format t "Substeps ~D~%" substeps)
    (loop for step from 0 below 1000
          while *run-sim*
          do
             (format t "Step ~D~%" step)
             (setf damage-est
                   (lparallel:pmap-reduce (lambda (mp)
                                            (*
                                             (cl-mpm/particle::mp-mass mp)
                                             (cl-mpm/particle::mp-damage mp)))
                                          #'+ (cl-mpm:sim-mps *sim*)
                                          :initial-value 0d0)
                   )
             (with-open-file (stream (merge-pathnames "timesteps.csv" output-dir) :direction :output :if-exists :append)
               (format stream "~D,~f,~f,~f,~f,~f,~A,~f~%"
                       step
                       (cl-mpm::sim-time *sim*)
                       damage-est
                       0d0 ;; *data-plastic*
                       energy
                       oobf
                       sim-state
                       0d0
                       ;; *data-mass*
                       ))
             (setf work 0d0)
             (setf energy 0d0
                   oobf 0d0
                   )
             (time
              (dotimes (i substeps)
                (cl-mpm:update-sim *sim*)
                (incf time (cl-mpm:sim-dt *sim*))
                ;; (cl-mpm::remove-mps-func
                ;;  *sim*
                ;;  (lambda (mp)
                ;;    ;; (> (cl-mpm/fastmaths:det (cl-mpm/particle:mp-deformation-gradient mp)) 1.5d0)
                ;;    (multiple-value-bind (s1 s2 s3) (cl-mpm/utils::principal-stresses-3d (cl-mpm/particle::mp-undamaged-stress mp) )
                ;;      (> (max s1 s2 s3) *removal-strain*))
                ;;    ))
                ;; (cl-mpm::split-mps-eigenvalue *sim*)
                (incf oobf (cl-mpm::sim-stats-oobf *sim*))
                (incf energy (cl-mpm::sim-stats-energy *sim*))
                (incf work (abs (cl-mpm::sim-stats-power *sim*)))
                ))
             (setf
              energy (/ energy substeps)
              oobf (/ oobf substeps))

             (if (= work 0d0)
                 (setf energy 0d0)
                 (setf energy (abs (/ energy work))))

             (let ((hist-factor criteria-hist))
               (format t "Hist-correct factors ~A ~A ~%" (* criteria-energy (+ 1d0 hist-factor)) (* criteria-oobf   (+ 1d0 hist-factor)))
               (when (or (>= energy (* criteria-energy hist-factor))
                         (>= oobf   (* criteria-oobf   hist-factor)))
                 (when (not (eq sim-state :collapse))
                   (setf sim-state :collapse)
                   (format t "Changed to collapse~%")
                   (setf work 0d0)))
               (when (and (< energy (/ criteria-energy criteria-hist))
                          (< oobf   (/ criteria-oobf   criteria-hist)))
                 (when (not (eq sim-state :accelerate))
                   (format t "Changed to accelerate~%")
                   (setf work 0d0)
                   (setf sim-state :accelerate)
                   (cl-mpm:iterate-over-mps
                    (cl-mpm:sim-mps *sim*)
                    (lambda (mp)
                      (cl-mpm/fastmaths::fast-zero (cl-mpm/particle:mp-velocity mp))
                      ))))
               (case sim-state
                 (:accelerate
                  (format t "Accelerate timestep~%")
                  (setf
                   target-time accelerate-target-time
                   (cl-mpm::sim-mass-scale *sim*) accelerate-mass-scale))
                 (:collapse
                  (format t "Collapse timestep~%")
                  (setf
                   ;; work 0d0
                   target-time collapse-target-time
                   (cl-mpm::sim-mass-scale *sim*) collapse-mass-scale))))


             (format t "OOBF ~E - Energy ~E~%" oobf energy)
             (format t "State ~A~%" sim-state)

             (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir (format nil "sim_~5,'0d.vtk" step)) *sim* )
             (cl-mpm/output:save-vtk-nodes (uiop:merge-pathnames* output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) *sim* )
             (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim_p_~5,'0d.vtk" step)) *sim* )
             (plot-domain)
             (vgplot:title (format nil "Time:~F - KE ~E - OOBF ~E - Work ~E - ~A"  time energy oobf work sim-state))
             ;; (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
             ;;   (format t "CFL dt estimate: ~f~%" dt-e)
             ;;   (format t "CFL step count estimate: ~D~%" substeps-e)
             ;;   (setf substeps substeps-e))

             (setf (cl-mpm:sim-dt *sim*) (* (sqrt (cl-mpm:sim-mass-scale *sim*)) dt-0))
             (setf substeps (round target-time (cl-mpm:sim-dt *sim*)))

             (setf (cl-mpm:sim-damping-factor *sim*)
                   (* 1d-4
                      (sqrt (cl-mpm:sim-mass-scale *sim*))
                      crit-damp))
             (swank.live:update-swank)
             (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step))
                                :terminal "png size 1920,1080"))))



(defun elastic-sim (&key (output-dir "./output/")
                      (damping 0.1d0)
                      (dt-scale 0.5d0)
                      )
  (uiop:ensure-all-directories-exist (list (uiop:merge-pathnames* output-dir)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (let ((ed (cl-mpm::sim-enable-damage *sim*)))
    (setf (cl-mpm::sim-enable-damage *sim*) nil)
    (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0)
    (setf (cl-mpm:sim-mass-scale *sim*) 1d0)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* damping (cl-mpm/setup:estimate-critical-damping *sim*)))
    (setf (cl-mpm:sim-dt *sim*) (* dt-scale (cl-mpm/setup:estimate-elastic-dt *sim*)))
    (cl-mpm/output::save-simulation-parameters (merge-pathnames output-dir "settings.json")
                                               *sim*
                                               (list :ocean-height *water-height*))
    (format t "Estimated dt ~E~%" (cl-mpm:sim-dt *sim*))
    (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir (format nil "sim_conv_~5,'0d.vtk" 0)) *sim*)
    (cl-mpm/output:save-vtk-nodes (uiop:merge-pathnames* output-dir (format nil "sim_conv_nodes_~5,'0d.vtk" 0)) *sim*)
    (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim_p_~5,'0d.vtk" 0)) *sim* )
    ;; (setf (cl-mpm::sim-velocity-algorithm *sim*) :PIC)
    (let ((current-step 0)
          (step-list (list))
          (energy-list (list))
          (oobf-list (list))
          (superstep 0)
          )
      (dotimes (i 50)
        (setf (cl-mpm:sim-enable-damage *sim*) nil)
        (cl-mpm/dynamic-relaxation:converge-quasi-static
         *sim*
         :oobf-crit 1d-3
         :energy-crit 1d-3
         :dt-scale dt-scale
         :kinetic-damping nil
         :conv-steps 1000
         :substeps 50
         :post-iter-step
         (lambda (i energy oobf)
           (when (and (< energy 1d-1) (< oobf 1d-1))
             (setf (cl-mpm:sim-enable-damage *sim*) t))
           (incf current-step)
           (push current-step step-list)
           (push energy energy-list)
           (push oobf oobf-list)
           (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir (format nil "sim_conv_~5,'0d_~5,'0d.vtk" superstep (1+ i))) *sim*)
          (cl-mpm/output:save-vtk-nodes (uiop:merge-pathnames* output-dir (format nil "sim_conv_nodes_~5,'0d_~5,'0d.vtk" superstep (1+ i))) *sim*)
          (cl-mpm/output:save-vtk-cells (uiop:merge-pathnames* output-dir (format nil "sim_conv_cells_~5,'0d_~5,'0d.vtk" superstep (1+ i))) *sim*)
          ;; (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim_p_~5,'0d.vtk" (1+ i))) *sim* )
           (vgplot:semilogy
            (reverse step-list)
            (reverse energy-list)
            "Energy"
            (reverse step-list)
            (reverse oobf-list)
            "OOBF")
           (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_conv_~5,'0d.png" i)) :terminal "png size 3840,2160")))
        (cl-mpm::finalise-loadstep *sim*)
        (let ((i superstep))
          (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir (format nil "sim_conv_~5,'0d.vtk" (1+ i))) *sim*)
          (cl-mpm/output:save-vtk-nodes (uiop:merge-pathnames* output-dir (format nil "sim_conv_nodes_~5,'0d.vtk" (1+ i))) *sim*)
          (cl-mpm/output:save-vtk-cells (uiop:merge-pathnames* output-dir (format nil "sim_conv_cells_~5,'0d.vtk" (1+ i))) *sim*)
          (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim_p_~5,'0d.vtk" (1+ i))) *sim* )
          (incf superstep)
          ))
      ;; (setf (cl-mpm::sim-enable-damage *sim*) t)
      ;; (cl-mpm/damage::calculate-damage *sim*)
      ;; (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir (format nil "sim_conv_~5,'0d.vtk" (1+ current-step))) *sim*)
      ;; (cl-mpm/output:save-vtk-nodes (uiop:merge-pathnames* output-dir (format nil "sim_conv_nodes_~5,'0d.vtk" (1+ current-step))) *sim*)
      )

    ;; (setf (cl-mpm::sim-enable-damage *sim*) ed)
    ;; (cl-mpm::iterate-over-mps
    ;;  (cl-mpm:sim-mps *sim*)
    ;;  (lambda (mp)
    ;;    (setf (cl-mpm/particle::mp-enable-plasticity mp) t)))
    ))

(defun est-angle ()
  (let* ((rc 0d0)
         (rs 0.5d0)
         (ratio (/ (- 1d0 rs) (- 1d0 rc)))
         (angle-plastic (* 40d0 (/ pi 180)))
         (angle-plastic-damaged (atan (* ratio (tan angle-plastic))))
         )
    (format t "Chalk plastic virgin angle: ~F~%"
            (* (/ 180 pi) angle-plastic))
    (format t "Chalk plastic residual angle: ~F~%"
            (* (/ 180 pi) angle-plastic-damaged))))


(defun conv-test (&key (output-dir "./output/"))
  (setup :refine 0.5 :mps 3)
  (defparameter *data-energy* (list))
  (defparameter *data-oobf* (list))
  (defparameter *data-t* (list))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (vgplot:close-all-plots)
  (vgplot:figure)
  (let ((ed (cl-mpm::sim-enable-damage *sim*))
        (dt-scale 0.5d0)
        (steps 1000)
        )
    (setf (cl-mpm::sim-enable-damage *sim*) nil)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 0.1d0
             (sqrt (cl-mpm:sim-mass-scale *sim*))
             (cl-mpm/setup:estimate-critical-damping *sim*))
          )
    (format t "Estimated dt ~E~%" (cl-mpm:sim-dt *sim*))
    (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir (format nil "sim_conv_~5,'0d.vtk" 0)) *sim*)
    (cl-mpm/output:save-vtk-nodes (uiop:merge-pathnames* output-dir (format nil "sim_conv_nodes_~5,'0d.vtk" 0)) *sim*)
    ;; (setf (cl-mpm::sim-velocity-algorithm *sim*) :PIC)
    (let ((substeps 10))
      (loop for i from 0 to steps
            do
               (let ((power 0d0))
                 (format t "Step ~D~%" i)
                 (dotimes (j substeps)
                   (cl-mpm:update-sim *sim*)
                   (incf power (cl-mpm::sim-stats-power *sim*))
                   )
                 (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir (format nil "sim_conv_~5,'0d.vtk" (1+ i))) *sim*)
                 (cl-mpm/output:save-vtk-nodes (uiop:merge-pathnames* output-dir (format nil "sim_conv_nodes_~5,'0d.vtk" (1+ i))) *sim*)
                 (let ((time (cl-mpm::sim-time *sim*))
                       (oobf (cl-mpm::sim-stats-oobf *sim*))
                       (eps (if (= power 0d0) 0d0 (abs (/ (cl-mpm::sim-stats-energy *sim*) power)))))
                   (push time *data-t*)
                   (push eps *data-energy*)
                   (push oobf *data-oobf*)
                   (format t "Energy ~E - OOBF ~E~%" eps oobf))
                 (setf (cl-mpm:sim-dt *sim*) (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
                 ;; (plot-domain)
                 (plot-conv-data)
                 (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_conv_~5,'0d.png" i)) :terminal "png size 3840,2160")
                 (swank.live:update-swank)
                 )))))

(defun plot-conv-data ()
  ;; (vgplot:figure)
  (vgplot:semilogy
   *data-t* *data-energy* "Energy"
   *data-t* *data-oobf* "OOBF"
   ))


(defun save-data (name)
  (with-open-file (stream (merge-pathnames name) :direction :output :if-exists :supersede)
    (format stream "time,oobf,energy~%")
    (loop for time in (reverse *data-t*)
          for oobf in (reverse *data-oobf*)
          for energy in (reverse *data-energy*)
          do (format stream "~f,~f,~f~%" time oobf energy))))



(defun stop ()
  (setf *run-sim* nil)
  (setf cl-mpm/dynamic-relaxation::*run-convergance* nil))

(defun body-pressure-test ()
  (let ((mps 3)
        (refine 0.25))
    (setup :refine refine
           :mps mps
           :pressure-condition t
           :friction 0.5d0
           ::bench-length 800d0
           :aspect 4
           )
    (elastic-sim :output-dir "./output-cryo-pressure/")
    ;; (setup :refine refine :mps mps :pressure-condition nil)
    ;; (elastic-sim :output-dir "./output-cryo-body/")
    ;; (setup :refine refine :mps mps :pressure-condition t :cryo-static nil)
    ;; (elastic-sim :output-dir "./output-zero-pressure/")
    ;; (setup :refine refine :mps mps :pressure-condition nil :cryo-static nil)
    ;; (elastic-sim :output-dir "./output-zero-body/")
    )
  )

(defun test-friction ()
  (let ((mps 2)
        (refine 0.5))
    (setup :refine refine :mps mps)
    (elastic-sim :output-dir "./output-0.0/")
    (setup :refine refine :mps mps :friction 0.5d0)
    (elastic-sim :output-dir "./output-0.5/")))

(defun test-bench ()
  (let ((mps 3)
        (refine 0.25)
        (friction 0.5d0))
    ;; (setup :refine refine :mps mps :friction friction)
    ;; (elastic-sim :output-dir "./output-nobench/")
    (setup :refine refine :mps mps :bench-length 1200d0  :friction friction
           :aspect 6
           )
    (elastic-sim :output-dir "./output-bench/"))
  )

;; (defun zero-body-pressure-test ()
;;   )


(defun calving-test ()
  (loop for dt in (list 1d4)
        do
           (let* ((mps 3)
                  (H 400d0)
                  )
             (setup :refine 0.25
                    :friction 0.2d0
                    :bench-length (* 0d0 H)
                    :ice-height H
                    :mps mps
                    :hydro-static nil
                    :cryo-static t
                    :melange nil
                    :aspect 4d0
                    :slope 0d0
                    :floatation-ratio 0.8d0)
             (plot-domain)
             (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0)
             (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
             (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
                   ;; (cl-mpm::sim-ghost-factor *sim*) (* 1d9 1d-3)
                   (cl-mpm::sim-ghost-factor *sim*) nil
                   ;; (cl-mpm/dynamic-relaxation::sim-convergence-critera *sim*) 1d-3
                   )
             (setf (cl-mpm:sim-enable-damage *sim*) nil)
             (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
             ;; (cl-mpm/dynamic-relaxation::run-time
             ;;  *sim*
             ;;  :output-dir "./output/"
             ;;  :dt 1d3
             ;;  :total-time 1d7
             ;;  :dt-scale 0.9d0
             ;;  ;; :dt-scale 1d4
             ;;  :damping 1d-4
             ;;  :mass-scale 1d7
             ;;  :enable-damage t
             ;;  :enable-plastic t
             ;;  :plotter (lambda (sim) (plot-domain))
             ;;  :post-conv-step
             ;;  (lambda (sim)
             ;;    (setf
             ;;     (cl-mpm/aggregate::sim-enable-aggregate sim) t
             ;;     (cl-mpm/damage::sim-enable-length-localisation sim) t
             ;;     (cl-mpm::sim-ghost-factor sim) nil
             ;;     ;; (cl-mpm/dynamic-relaxation::sim-convergence-critera *sim*) 1d-3
             ;;     (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0))

             ;;  )
             (cl-mpm/dynamic-relaxation::run-multi-stage
              *sim*
              :output-dir "./output/"
              :dt dt
              :dt-scale 1d0
              :damping-factor 1d0;(sqrt 2)
              :conv-criteria 1d-2
              :conv-load-steps 2
              :min-adaptive-steps -4
              :max-adaptive-steps 10
              :substeps 10
              :steps 1000
              :enable-plastic t
              :enable-damage t
              ;; :post-conv-step (lambda (sim)
              ;;                   (setf (cl-mpm/penalty::bc-penalty-friction *floor-bc*) 0.9d0))
              :plotter (lambda (sim) (plot-domain))
              ;; :explicit-dt-scale 100d0
              :explicit-dt-scale 0.25d0
              :explicit-damping-factor 1d-4
              :explicit-dynamic-solver
              ;; 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic
              'cl-mpm/damage::mpm-sim-agg-damage
              :setup-quasi-static
              (lambda (sim)
                ;; (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
                (setf
                 (cl-mpm/aggregate::sim-enable-aggregate sim) t
                 (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC
                 (cl-mpm::sim-ghost-factor sim) nil
                 ;; (cl-mpm::sim-ghost-factor sim) (* 1d9 1d-3)
                 (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0))
              :setup-dynamic
              (lambda (sim)
                ;; (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
                (setf (cl-mpm/aggregate::sim-enable-aggregate sim) nil
                      (cl-mpm::sim-velocity-algorithm sim) :BLEND
                      (cl-mpm::sim-ghost-factor sim) nil;(* 1d9 1d-4)
                      (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 1d0
                      ))
              )
             )))



(defun test-square ()
  (defparameter *bc-square*
    (cl-mpm/penalty::make-bc-penalty-square
     nil
     (cl-mpm/utils:vector-from-list (list 0d0 0d0 1d0));;vector
     (cl-mpm/utils:vector-from-list (list 0d0 0d0 0d0));;vector
     (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0));;vector
     (cl-mpm/utils:vector-from-list (list 0d0 1d0 0d0));;vector
     1d0 0d0 0d0
     ))
  )

(defun test ()
  (setup :refine 1;0.25
         :mps 2
         :aspect 1
         :bench-length 000d0
         )
  (elastic-sim
   :damping 1d-2
   :dt-scale 0.8d0
   )
  )

(defun test-refined-elastic ()

  (let ((mps 2)
        (refine-0 0.125)
        (mesh-size 10)
        (last-mps nil)
        (refine-steps 2)
        )
    (format t "Multi-grid~%")
    (time
     (dotimes (i (1+ refine-steps))
       (print i)
       (let* ((refine (* (expt 2 i) refine-0))
              (h (/ mesh-size refine))
              (offset (* 0 h))
              (mps (/ (expt 2 (1+ refine-steps)) (expt 2 i)))
              )
         (setup :refine refine :mps mps
                :ice-height 800d0
                :aspect 4
                :mps mps
                )
         (format t "MPs ~D~%" mps)
         (unless (= i 0)
           (cl-mpm::remove-mps-func
            *sim*
            (lambda (mp)
              t))
           (cl-mpm:add-mps *sim* last-mps)
           )
         (elastic-sim :output-dir (format nil "./output-~E/" refine)
                      :dt-scale 0.25d0
                      )
         (cl-mpm:iterate-over-mps
          (cl-mpm:sim-mps *sim*)
          (lambda (mp)
            (cl-mpm/fastmaths::fast-zero (cl-mpm/particle:mp-velocity mp))
            ))
         (setf last-mps (cl-mpm:sim-mps *sim*)))
       ))
    (format t "Single-grid~%")
    (time
     (let ((i refine-steps))
       (print i)
       (let* ((refine (* (expt 2 i) refine-0))
              (h (/ mesh-size refine))
              (offset (* 0 h))
              (mps (/ (expt 2 (1+ refine-steps)) (expt 2 i)))
              )
         (setup :refine refine :mps mps
                :ice-height 800d0
                :aspect 4
                :mps mps
                )
         (format t "MPs ~D~%" mps)
         (elastic-sim :output-dir (format nil "./output-single-~E/" refine)
                      :dt-scale 0.25d0))))))




(declaim (notinline quasi-time))
(defun quasi-time (&key (output-dir "./output/")
                     (time-steps 10)
                     (time 1d0)
                     (enable-damage t)
                     (enable-plastic t)
                     (dt-scale 0.5d0))
  (uiop:ensure-all-directories-exist (list output-dir))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (setf (cl-mpm:sim-mass-scale *sim*) 1d0)
  (setf (cl-mpm:sim-damping-factor *sim*)
        (* 1d-1
           (cl-mpm/setup:estimate-critical-damping *sim*)))
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (with-open-file (stream (merge-pathnames output-dir "conv.csv") :direction :output :if-exists :supersede)
      (format stream "iter,step,damage,plastic,oobf,energy~%"))
    (with-open-file (stream (merge-pathnames output-dir "./timesteps.csv") :direction :output :if-exists :supersede)
      (format stream "steps,time,damage,plastic,energy,oobf,step-type,mass~%"))
    (with-open-file (stream (merge-pathnames "timesteps.csv" output-dir) :direction :output :if-exists :append)
      (format stream "~D,~f,~f,~f,~f,~f,~A,~f~%"
              0
              (cl-mpm::sim-time *sim*)
              (get-damage)
              (get-plastic)
              0d0
              0d0
              :QUASI-STATIC
              0d0
              ))
    (let ((step 0))
      (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) *sim*)
      (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" step)) *sim*)
      (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d.vtk" step)) *sim*)
      (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim_p_~5,'0d.vtk" step)) *sim* ))
    (let* ((total-iter 0))
      (loop for step from 1 to time-steps
            while *run-sim*
            do
               (let ((crit (cl-mpm/setup:estimate-critical-damping *sim*)))
                 (setf (cl-mpm:sim-damping-factor *sim*)
                       (* 0.1d0 crit))
                 ;; (setf (cl-mpm::sim-velocity-algorithm *sim*) :PIC)
                 (setf (cl-mpm/dynamic-relaxation::sim-dt-loadstep *sim*) (/ time time-steps))
                 (defparameter *ke-last* 0d0)
                 (let ((conv-steps 0)
                       (substeps 10)
                       (damping 1d-1)
                       (damage-prev 0d0))
                   (time
                    (progn
                      (cl-mpm:iterate-over-mps
                       (cl-mpm:sim-mps *sim*)
                       (lambda (mp)
                         (when (typep mp 'cl-mpm/particle::particle-damage)
                           (setf (cl-mpm/particle::mp-enable-damage mp) nil))
                         (when (typep mp 'cl-mpm/particle::particle-plastic)
                           (setf (cl-mpm/particle::mp-enable-plasticity mp) nil))))
                      (setf (cl-mpm:sim-enable-damage *sim*) nil)
                      (cl-mpm/dynamic-relaxation:converge-quasi-static
                       *sim*
                       :oobf-crit 1d-2
                       :energy-crit 1d-2
                       :kinetic-damping t
                       :dt-scale dt-scale
                       :substeps substeps
                       :conv-steps 1000
                       :damping-factor 1d-2
                       :post-iter-step
                       (lambda (i e o)
                         (plot-domain)
                         (incf conv-steps)
                         (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d_~5,'0d.vtk" step conv-steps)) *sim*)
                         (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d_~5,'0d.vtk" step conv-steps)) *sim*)
                         (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d_~5,'0d.vtk" step conv-steps)) *sim*) (incf total-iter substeps)
                         (with-open-file (stream (merge-pathnames output-dir "conv.csv") :direction :output :if-exists :append)
                           (format stream "~D,~D,~f,~f,~f,~f~%" total-iter step (get-plastic) (get-damage) 
                                   o e)))
                       )
                      (cl-mpm:iterate-over-mps
                       (cl-mpm:sim-mps *sim*)
                       (lambda (mp)
                         (when (typep mp 'cl-mpm/particle::particle-damage)
                           (setf (cl-mpm/particle::mp-enable-damage mp) enable-damage))
                         (when (typep mp 'cl-mpm/particle::particle-plastic)
                           (setf (cl-mpm/particle::mp-enable-plasticity mp) enable-plastic))))
                      (when enable-damage
                        (cl-mpm/damage:calculate-damage *sim* (cl-mpm/dynamic-relaxation::sim-dt-loadstep *sim*)))
                      (setf (cl-mpm:sim-enable-damage *sim*) enable-damage)
                      (let ((damage-prev (get-damage)))
                        (cl-mpm/dynamic-relaxation:converge-quasi-static
                         *sim*
                         :oobf-crit 1d-2
                         :energy-crit 1d-2
                         :convergance-criteria
                         (lambda (sim fnorm oobf)
                           (let* ((dn (get-damage))
                                  (dconv (if (> dn 0d0)
                                             (if (> damage-prev 0d0) (/ (- dn damage-prev) damage-prev) sb-ext:double-float-positive-infinity)
                                             0d0
                                             )))
                             (format t "d-conv ~E~%" dconv)
                             (setf damage-prev dn)
                             (and
                              (< fnorm 1d-2)
                              (< oobf 1d-2)
                              (< dconv 1d-3)
                              )))
                         :kinetic-damping t
                         :dt-scale dt-scale
                         :substeps substeps
                         :conv-steps 100
                         :damping-factor 1d-2
                         :post-iter-step
                         (lambda (i e o)
                           (plot-domain)
                           (incf conv-steps)
                           (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d_~5,'0d.vtk" step conv-steps)) *sim*)
                           (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d_~5,'0d.vtk" step conv-steps)) *sim*)
                           (incf total-iter substeps)
                           (with-open-file (stream (merge-pathnames output-dir "conv.csv") :direction :output :if-exists :append)
                             (format stream "~D,~D,~f,~f,~f,~f~%" total-iter step (get-plastic) (get-damage) 
                                     o e)))))))
                   (vgplot:title (format nil "Step ~D - ~D" step conv-steps))
                   ))
               ;; (cl-mpm/penalty::bc-increment-center *bc-squish* (cl-mpm/utils:vector-from-list (list 0d0 load-inc 0d0)))
               (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) *sim*)
               (cl-mpm::finalise-loadstep *sim*)
               (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step))
                                  :terminal "png size 1920,1080"
                                  )
               (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" step)) *sim*)
               (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d.vtk" step)) *sim*)
               (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim_p_~5,'0d.vtk" step)) *sim* )
               (with-open-file (stream (merge-pathnames "timesteps.csv" output-dir) :direction :output :if-exists :append)
                 (format stream "~D,~f,~f,~f,~f,~f,~A,~f~%"
                         step
                         (cl-mpm::sim-time *sim*)
                         (get-damage)
                         (get-plastic)
                         0d0
                         0d0
                         :QUASI-STATIC
                         0d0
                         ))
               (sleep 0.1d0)
               (swank.live:update-swank)

            ))))

(defun save-test-vtks (&key (output-dir "./output/"))
  (cl-mpm::finalise-loadstep *sim*)
  (cl-mpm/output:save-vtk (merge-pathnames "test.vtk" output-dir) *sim*)
  (cl-mpm/output:save-vtk-nodes (merge-pathnames "test_nodes.vtk" output-dir) *sim*)
  (cl-mpm/output:save-vtk-cells (merge-pathnames "test_cells.vtk" output-dir) *sim*)
  )


(defun agg-test ()
  (let* ((mesh (cl-mpm:sim-mesh *sim*))
         (cell (cl-mpm/mesh::get-cell mesh (list 0 0 0)))
         (pos (cl-mpm/utils:vector-from-list (list 10d0 0d0 0d0)))
         )
    (pprint cell)
    (pprint (cl-mpm/mesh:mesh-resolution mesh))
    (cl-mpm/aggregate::find-local-coords-agg mesh cell pos)
    ))

(defun calving-quasi-time-test ()
  (loop for dt in (list 1d4)
        do
           (let* ((mps 2))
             (setup :refine 0.125
                    :friction 0d0
                    :bench-length 000d0
                    :ice-height 400d0
                    :mps mps
                    :cryo-static t
                    :aspect 1d0)
             (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) nil)
             (plot-domain)
             (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0)
             (setf (cl-mpm::sim-dt-scale *sim*) 1d0)
             (cl-mpm/dynamic-relaxation::run-quasi-time
              *sim*
              :output-dir "./output/"
              :dt dt
              :dt-scale 1d0
              :enable-plastic t
              :enable-damage t
              :steps 1000
              :plotter (lambda (sim) (plot-domain))
              :post-conv-step
              (lambda (sim)
                (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0))))))

(defun calving-multi-stage ()
  (loop for dt in (list 1d4)
        do
           (let* ((mps 2))
             (setup :refine 1
                    :friction 0d0
                    :bench-length 000d0
                    :ice-height 400d0
                    :floatation-ratio 0.7
                    :mps mps
                    :cryo-static t
                    :aspect 1d0)
             (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
             (plot-domain)
             (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0)
             (setf (cl-mpm::sim-dt-scale *sim*) 1d0)
             (cl-mpm/dynamic-relaxation::run-multi-stage
              *sim*
              :output-dir "./output/"
              :dt dt
              :dt-scale 1d0
              :enable-plastic nil
              :enable-damage t
              :steps 1000
              :conv-criteria 1d-3
              ;; :post-conv-step (lambda (sim)
              ;;                   (setf (cl-mpm/penalty::bc-penalty-friction *floor-bc*) 0.5d0))
              :plotter (lambda (sim) (plot-domain))
              :max-adaptive-steps 10
              :explicit-dt-scale 0.5d0
              :elastic-dt-margin 1d3
              :explicit-damping-factor 1d-3
              :explicit-dynamic-solver 'cl-mpm/damage::mpm-sim-agg-damage
              ;; :explicit-dynamic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic
              :setup-quasi-static
              (lambda (sim)
                (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
                (setf
                 (cl-mpm/aggregate::sim-enable-aggregate sim) t
                 (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC
                 (cl-mpm::sim-ghost-factor sim) nil;; (* 1d9 1d-4)
                 (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0))
              :setup-dynamic
              (lambda (sim)
                (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-5)
                (setf (cl-mpm/aggregate::sim-enable-aggregate sim) nil
                      (cl-mpm::sim-velocity-algorithm sim) :FLIP
                      (cl-mpm::sim-ghost-factor sim) nil;(* 1d9 1d-4)
                      (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0))
              ))))
(defun plot-domain-sxy ()
  (when *sim*
    (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
           (h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
           (ms-x (first ms))
           (ms-y (second ms)))
      (vgplot:format-plot t "set object 1 rect from 0,0 to ~f,~f fc rgb 'blue' fs transparent solid 0.5 noborder behind" ms-x *water-height*))
    (cl-mpm::g2p (cl-mpm:sim-mesh *sim*)
                 (cl-mpm:sim-mps *sim*)
                 (cl-mpm:sim-dt *sim*)
                 0d0
                 :TRIAL)
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :trial t
     :colour-func (lambda (mp) (cl-mpm/utils:varef (cl-mpm/particle::mp-stress mp) 5))
     ))
  )
(defun plot-domain-ybar ()
  (when *sim*
    (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
           (h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
           (ms-x (first ms))
           (ms-y (second ms)))
      (vgplot:format-plot t "set object 1 rect from 0,0 to ~f,~f fc rgb 'blue' fs transparent solid 0.5 noborder behind" ms-x *water-height*))
    (cl-mpm::g2p (cl-mpm:sim-mesh *sim*)
                 (cl-mpm:sim-mps *sim*)
                 (cl-mpm:sim-dt *sim*)
                 0d0
                 :TRIAL)
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :trial t
     :colour-func #'cl-mpm/particle::mp-damage-ybar
     ))
  )

(defun test-static-damage-ybar ()
  (let ((output-dir "./output/"))
    (uiop:ensure-all-directories-exist (list output-dir))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
    (setup :refine 0.125
           :friction 0.9d0
           :bench-length 100d0
           :ice-height 1000d0
           :mps 2
           :cryo-static t
           :aspect 4d0
           :floatation-ratio 0.8d0)
    (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t)
    (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
    (setf (cl-mpm:sim-mass-scale *sim*) 1d0)
    (cl-mpm/dynamic-relaxation::save-conv-preamble output-dir)
    (cl-mpm/dynamic-relaxation::elastic-static-solution
     *sim*
     :criteria 1d-3
     :post-iter-step
     (lambda (i e o)
       (plot-domain-sxy)
       (cl-mpm/dynamic-relaxation::save-conv-step *sim* output-dir i 0 0d0 o e)
       (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_conv_~5,'0d.vtk" i)) *sim*)
       (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_conv_nodes__~5,'0d.vtk" i)) *sim*)
       (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_conv_cells__~5,'0d.vtk" i)) *sim*)))
    (setf (cl-mpm:sim-enable-damage *sim*) t)
    (cl-mpm/dynamic-relaxation::set-mp-plastic-damage *sim* :damage t)
    (cl-mpm/damage::calculate-damage *sim* 1d-15)
    (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir "sim.vtk") *sim*)
    (cl-mpm/output:save-vtk-nodes (uiop:merge-pathnames* output-dir "sim_nodes.vtk") *sim*))
  (plot-domain-ybar)
  )
  ;)


(defun test-visco ()
  (let* ((mps 4)
         (dt 1d12))
    (setup :refine 0.125
           :friction 0d0
           :bench-length 000d0
           :ice-height 400d0
           :mps mps
           :cryo-static t
           :aspect 4d0
           :floatation-ratio 0.5d0)
    (plot-domain)
    (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0)
    (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
    (cl-mpm/setup::set-mass-filter *sim* 900d0 :proportion 1d-9)
    (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t)
    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (setf (cl-mpm/particle::mp-enable-viscosity mp) nil)))
    (cl-mpm/dynamic-relaxation::run-quasi-time
     *sim*
     :output-dir "./output/"
     :dt dt
     :dt-scale 0.25d0
     :substeps 10
     :steps 1000
     :max-adaptive-steps 0
     :enable-plastic t
     :plotter (lambda (sim) (plot-domain))
     :post-conv-step
     (lambda (sim)
       (cl-mpm::iterate-over-mps
        (cl-mpm:sim-mps *sim*)
        (lambda (mp)
          (setf (cl-mpm/particle::mp-enable-viscosity mp) t)))))
    )
  )

(defun calving-lc-test ()
  (let* ((mps 16))
    (setup :refine 0.125
           :friction 0d0
           :bench-length 0d0
           :ice-height 600d0
           :mps mps
           :cryo-static nil
           :aspect 2d0
           :slope 0d0
           :floatation-ratio 0.9d0)
    (plot-domain)
    (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0)
    (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
    (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
          ;; (cl-mpm::sim-ghost-factor *sim*) (* 1d9 1d-3)
          (cl-mpm::sim-ghost-factor *sim*) nil
          )
    (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-9)
    (cl-mpm/dynamic-relaxation::run-adaptive-load-control
     *sim*
     :output-dir "./output/"
     :load-steps 10
     :initial-load 0.5d0
     ;; :loading-function
     ;; (lambda (percent)
     ;;   (format t "Loading factor ~E~%" percent)
     ;;   (setf (cl-mpm:sim-gravity *sim*) (* -9.8d0 percent)))
     :dt-scale 1d0
     :conv-criteria 1d-3
     :substeps 10
     :steps 1000
     :enable-damage t
     :enable-plastic t
     :save-vtk-dr nil
     :max-adaptive-steps 10
     :post-conv-step (lambda (sim) (plot-domain))
     )))

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

(defun stability-test ()
  (let* ((heights (reverse (list 200d0 300d0 400d0 500d0 600d0))
                  )
         (floatations (list 0.95d0 0.9d0 0.85d0 0.8d0 0.75d0)))
    (defparameter *stability* (make-array (list (length heights) (length floatations)) :initial-element nil
                                                                                       :element-type t))
    (defparameter *heights* heights)
    (defparameter *floatations* floatations)
    (loop for hi from 0
          for height in heights
          do
             (let ((res t))
               (loop for fi from 0
                     for flotation in floatations
                     do
                        (let* ((mps 2)
                               (output-dir (format nil "./output-~f-~f/" height flotation)))
                          (setup :refine 0.5
                                 :friction 0.5d0
                                 :bench-length 0d0
                                 :ice-height height
                                 :mps mps
                                 :cryo-static nil
                                 :aspect 2d0
                                 :slope 0.1d0
                                 :floatation-ratio flotation)
                          (plot-domain)
                          (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0)
                          (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) nil)
                          (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
                                ;; (cl-mpm::sim-ghost-factor *sim*) (* 1d9 1d-3)
                                (cl-mpm::sim-ghost-factor *sim*) nil
                                )
                          (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-9)
                          (let ((res (cl-mpm/dynamic-relaxation::run-adaptive-load-control
                                      *sim*
                                      :output-dir output-dir
                                      :load-steps 10
                                      :initial-load 0.5d0
                                      :loading-function
                                      (lambda (percent)
                                        (format t "Loading factor ~E~%" percent)
                                        (setf (cl-mpm:sim-gravity *sim*) (* -9.8d0 percent)))
                                      :dt-scale 1d0
                                      :conv-criteria 1d-3
                                      :substeps 10
                                      :steps 1000
                                      :enable-damage t
                                      :enable-plastic t
                                      :max-adaptive-steps 6
                                      :post-conv-step (lambda (sim) (plot-domain)))))
                            (setf (aref *stability* hi fi) (if res 1 0))
                            (unless res
                              (loop for j from fi below (length floatations)
                                    do (setf (aref *stability* hi j) 0)))
                            (format t "Stability:~%")
                            (print-stab)
                            (unless res
                              (loop-finish)))
                          ))))
    ))

(defun stability-qt-test ()
  (let* ((heights
           (list 1000)
           ;(reverse (list 200d0 300d0 400d0 500d0 600d0))
                  )
         (floatations (list 0.95d0 0.9d0 0.85d0 0.8d0 0.75d0)))
    (defparameter *stability* (make-array (list (length heights) (length floatations)) :initial-element nil
                                                                                       :element-type t))
    (defparameter *heights* heights)
    (defparameter *floatations* floatations)
    (loop for hi from 0
          for height in heights
          do
             (let ((res t))
               (loop for fi from 0
                     for flotation in floatations
                     do
                        (let* ((mps 2)
                               (output-dir (format nil "./output-~f-~f/" height flotation)))
                          (format t "Problem ~f ~f~%" height flotation)
                          (setup :refine (expt 2 -4)
                                 :friction 0.5d0
                                 :bench-length 0d0
                                 :ice-height height
                                 :mps mps
                                 :cryo-static nil
                                 :aspect 1d0
                                 :slope 0d0
                                 :floatation-ratio flotation)
                          (plot-domain)
                          (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0)
                          (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) nil)
                          (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
                                ;; (cl-mpm::sim-ghost-factor *sim*) (* 1d9 1d-3)
                                (cl-mpm::sim-ghost-factor *sim*) nil
                                )
                          (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-9)
                          (let ((res (cl-mpm/dynamic-relaxation::run-quasi-time
                                      *sim*
                                      :output-dir output-dir
                                      :dt 1d3
                                      :total-time 1d6
                                      :dt-scale 1d0
                                      :conv-criteria 1d-3
                                      :substeps 50
                                      :enable-damage t
                                      :enable-plastic t
                                      :min-adaptive-steps -4
                                      :max-adaptive-steps 6
                                      :save-vtk-dr nil
                                      :elastic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-dr-multigrid
                                      :plotter (lambda (sim) (plot-domain))
                                      :post-conv-step (lambda (sim) (plot-domain)))))
                            (setf (aref *stability* hi fi) (if res 1 0))
                            (unless res
                              (loop for j from fi below (length floatations)
                                    do (setf (aref *stability* hi j) 0)))
                            (format t "Stability:~%")
                            (print-stab)
                            (unless res
                              (loop-finish)))
                          ))))
    ))


