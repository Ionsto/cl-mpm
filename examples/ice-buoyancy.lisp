(defpackage :cl-mpm/examples/ice-buoyancy
  (:use :cl
   :cl-mpm/example))
(in-package :cl-mpm/examples/ice-buoyancy)

(defclass cl-mpm/particle::particle-mc-erodable (cl-mpm/particle::particle-mc
                                                    cl-mpm/particle::particle-erosion)
  ())

(defmethod cl-mpm/erosion::mp-erosion-enhancment ((mp cl-mpm/particle::particle-mc-erodable))
  ;; (+ 1d0 (* 10 (cl-mpm/particle::mp-damage mp)))
  (+ 1d0 (* 10 (cl-mpm/particle::mp-strain-plastic-vm mp)))
  )

(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 0 0)

(defun cl-mpm/damage::length-localisation (local-length local-length-damaged damage)
  (declare (double-float local-length damage))
  (* local-length (max (- 1d0 damage) 1d-10)))

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-chalk-delayed) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  ;; (cl-mpm::update-domain-det mesh mp dt)
  ;; (cl-mpm::co-domain-corner-2d mesh mp dt)
  (cl-mpm::update-domain-polar-2d mesh mp dt)
  ;; ;; (cl-mpm::update-domain-midpoint mesh mp dt)
  (cl-mpm::scale-domain-size mesh mp)
  )
(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-chalk-delayed) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
                     (strain cl-mpm/particle::mp-strain)
                     (damage cl-mpm/particle:mp-damage)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (angle cl-mpm/particle::mp-friction-angle)
                     (de cl-mpm/particle::mp-elastic-matrix)
                     (E cl-mpm/particle::mp-e)
                     )
        mp
      (progn
        (setf damage-increment
              (min
               ;; (cl-mpm/damage::tensile-energy-norm strain E de)
               (max 0d0
                    (cl-mpm/damage::criterion-mohr-coloumb-rankine-stress-tensile
                     stress
                     (* angle (/ pi 180d0))))
               ))
        (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
        (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
        ))))


(defun plot-domain ()
  (when *sim*
    ;; (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
    ;;        (h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
    ;;        (floor-datum (* 2 h))
    ;;        (ms-x (first ms))
    ;;        (ms-y (second ms)))
    ;;   (vgplot:format-plot t "set object 1 rect from 0,0 to ~f,~f fc rgb 'blue' fs transparent solid 0.5 noborder behind" ms-x *water-height*)
    ;;   (vgplot:format-plot t "set style fill solid")
    ;;   (vgplot:format-plot t "set yrange [~f:~f]" floor-datum ms-y)
    ;;   (vgplot:format-plot t "set size ratio ~f" (/ (- ms-y floor-datum) ms-x)))
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     ;; :colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle:mp-stress mp) :xy))
     ;; :colour-func #'cl-mpm/particle::mp-strain-plastic-vm
     ;; :colour-func (lambda (mp) (/ (cl-mpm/particle:mp-mass mp)
     ;;                              (cl-mpm/particle:mp-volume mp)))
     :colour-func #'cl-mpm/particle::mp-damage
     )))
(defun setup (&key (refine 1) (mps 2))
  (let* ((density 916.7d0)
         (mesh-resolution (/ 10d0 refine))
         (offset (* mesh-resolution 2))
         (datum (+ 100d0 offset))
         (domain-size (list 800d0 400d0))
         (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
         (block-size (list 400d0 200d0)))
    (setf *sim* (cl-mpm/setup::make-simple-sim mesh-resolution element-count
                                               :sim-type
                                               'cl-mpm/damage::mpm-sim-damage
                                               ;; 'cl-mpm::mpm-sim-usf
                                               ))

    (let* ((E 1d9)
           (angle 60d0)
           (init-c 50d3)
           (init-stress (cl-mpm/damage::mohr-coloumb-coheasion-to-tensile init-c (* angle (/ pi 180))))
           (gf 40d0)
           (length-scale (* mesh-resolution 2d0))
           (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress E))
           (oversize (cl-mpm/damage::compute-oversize-factor 0.9999d0 ductility)))
      (format t "Mesh size ~F~%" mesh-resolution)
      (format t "Estimated oversize ~F~%" oversize)
      (format t "Estimated lc ~E~%" length-scale)
      (format t "Estimated ductility ~E~%" ductility)
      (format t "Init stress ~E~%" init-stress)
      (defparameter *removal-strain* (* 100d0 (/ init-stress 1)))
      (format t "Removal strain ~E~%" *removal-strain*)
      (cl-mpm:add-mps
       *sim*
       (cl-mpm/setup:make-block-mps
        (list 0 offset)
        block-size
        (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
        density
        'cl-mpm/particle::particle-chalk-delayed
        :E 1d9
        :nu 0.24d0

        :kt-res-ratio 1d0
        :kc-res-ratio 0d0
        :g-res-ratio 0.9d0

        :initiation-stress init-stress;18d3
        :friction-angle angle
        :psi (* 0d0 (/ pi 180))
        :phi (* angle (/ pi 180))
        :c (* init-c oversize)
        :softening 0d0
        :ductility ductility
        :local-length length-scale
        :delay-time 1d1
        :delay-exponent 1d0
        :enable-plasticity nil
        :enable-damage t
        ;; 'cl-mpm/particle::particle-finite-viscoelastic-ice
        ;; :E 1d9
        ;; :nu 0.325d0
        ;; :visc-factor 11.1d6
        ;; :visc-power 3d0
        ;; :enable-viscosity nil

        :gravity -9.8d0
        )))
    (setf
     (cl-mpm:sim-bcs *sim*)
     (cl-mpm/bc::make-outside-bc-varfix
      (cl-mpm:sim-mesh *sim*)
      '(0 nil 0)
      '(0 nil 0)
      '(nil 0 0)
      '(nil 0 0)
      '(nil nil 0)
      '(nil nil 0)))

    (setf (cl-mpm:sim-mass-scale *sim*) 1d4)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 0.1d0
             (sqrt 1d4)
             (cl-mpm/setup:estimate-critical-damping *sim*)))
    ;; (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-4)
    (setf (cl-mpm::sim-enable-fbar *sim*) t)
    (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
    (setf (cl-mpm::sim-allow-mp-split *sim*) t)
    ;; (setf (cl-mpm::sim-velocity-algorithm *sim*) :PIC)
    (setf (cl-mpm::sim-velocity-algorithm *sim*) :BLEND)



    (setf (cl-mpm:sim-dt *sim*)
          (* 0.5d0 (cl-mpm/setup:estimate-elastic-dt *sim*)))
    (setf (cl-mpm::sim-enable-damage *sim*) nil)
    (setf *run-sim* t)
    (if t
        (cl-mpm:add-bcs-force-list
         *sim*
         (cl-mpm/buoyancy::make-bc-buoyancy-clip
          *sim*
          datum
          1000d0
          (lambda (pos datum) t)))
        (cl-mpm:add-bcs-force-list
         *sim*
         (cl-mpm/buoyancy::make-bc-buoyancy-body
          *sim*
          datum
          1000d0
          (lambda (pos) t))))
    (let ((domain-half (* 0.5d0 (first domain-size)))
          (friction 0.9d0))
      (defparameter *floor-bc*
        (cl-mpm/penalty::make-bc-penalty-distance-point
         *sim*
         (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
         (cl-mpm/utils:vector-from-list (list
                                         domain-half
                                         offset
                                         0d0))
         (* domain-half 1.1d0)
         (* 1d9 0.1d0)
         friction
         0d0)))

    ;; (defparameter *bc-erode*
    ;;   (cl-mpm/erosion::make-bc-erode
    ;;    *sim*
    ;;    :enable nil
    ;;    :rate 1d-1
    ;;    :scalar-func (lambda (pos)
    ;;                   (min 1d0 (exp (* 0.5d0 (- (cl-mpm/utils:varef pos 1) datum)))))
    ;;    :clip-func (lambda (pos)
    ;;                 (>= datum (cl-mpm/utils:varef pos 1))
    ;;                 ;; (and (>= (+ offset mesh-resolution) (cl-mpm/utils:varef pos 1)))
    ;;                 )))
    (cl-mpm:add-bcs-force-list
     *sim*
     *floor-bc*
     )
    ;; (cl-mpm:add-bcs-force-list
    ;;  *sim*
    ;;  *bc-erode*
    ;;  )
    (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
    ))



(defun run (&key (output-dir "./output/"))
  (uiop:ensure-all-directories-exist (list (uiop:merge-pathnames* output-dir)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (let ((dt-scale 0.5d0))
    (cl-mpm/dynamic-relaxation:converge-quasi-static
     *sim*
     :oobf-crit 1d-1
     :energy-crit 1d-1
     :dt-scale dt-scale
     :post-iter-step
     (lambda (&rest args)
       (plot-domain))))

  (cl-mpm::iterate-over-mps
   (cl-mpm:sim-mps *sim*)
   (lambda (mp)
     (setf (cl-mpm/particle::mp-enable-plasticity mp) nil)))

  (setf (cl-mpm::sim-enable-damage *sim*) t)
  (setf (cl-mpm:sim-mass-scale *sim*) 1d0)
  (setf (cl-mpm:sim-damping-factor *sim*)
        (* 0d-4
           ;; (sqrt (cl-mpm:sim-mass-scale *sim*))
           (cl-mpm/setup:estimate-critical-damping *sim*)))

  ;; (setf (cl-mpm/buoyancy::bc-enable *bc-erode*) t)
  (let* ((dt-scale 0.5d0)
         (target-time 1d0)
         (substeps 0d0)
         (work 0d0)
         (oobf 0d0)
         (energy 0d0)
         )
    (setf (cl-mpm:sim-dt *sim*) (* dt-scale (cl-mpm/setup:estimate-elastic-dt *sim*)))
    (setf substeps (ceiling target-time (cl-mpm:sim-dt *sim*)))
    (format t "Substeps ~D~%" substeps)
    (loop for step from 0 below 1000
          while *run-sim*
          do
             (format t "Step ~D~%" step)
             (setf work 0d0)
             (time
              (dotimes (i substeps)
                (cl-mpm:update-sim *sim*)

                ;; (cl-mpm::remove-mps-func
                ;;  *sim*
                ;;  (lambda (mp)
                ;;    ;; (> (cl-mpm/fastmaths:det (cl-mpm/particle:mp-deformation-gradient mp)) 1.5d0)
                ;;    (multiple-value-bind (s1 s2 s3) (cl-mpm/utils::principal-stresses-3d (cl-mpm/particle::mp-undamaged-stress mp) )
                ;;      (> (max s1 s2 s3) *removal-strain*))
                ;;    ))
                ;; (cl-mpm::split-mps-eigenvalue *sim*)
                (incf work (cl-mpm/dynamic-relaxation::estimate-power-norm *sim*))))
             (setf oobf (cl-mpm/dynamic-relaxation::estimate-oobf *sim*))
             (setf energy (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
             (setf
              energy (/ energy substeps)
              oobf (/ oobf substeps))
             (if (= work 0d0)
                 (setf energy 0d0)
                 (setf energy (abs (/ energy work))))

             (format t "OOBF ~E - Energy ~E~%" oobf energy)

             (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir (format nil "sim_~5,'0d.vtk" step)) *sim* )
             (cl-mpm/output:save-vtk-nodes (uiop:merge-pathnames* output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) *sim* )
             (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim_p_~5,'0d.vtk" step)) *sim* )
             (plot-domain)
             (swank.live:update-swank)
             (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step))
                                :terminal "png size 1920,1080"
                                )
             )))

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
         (load-inc (* (- 1d0 elastic-load) (/ load step-count)))
         )
    (format t "Substeps ~D~%" substeps)
    ;; (cl-mpm::iterate-over-mps
    ;;  (cl-mpm:sim-mps *sim*)
    ;;  (lambda (mp)
    ;;    (setf (cl-mpm/particle::mp-gravity mp) load-0)))
    ;; (cl-mpm/dynamic-relaxation:converge-quasi-static
    ;;  *sim*
    ;;  :energy-crit 1d-2
    ;;  :oobf-crit 1d-2
    ;;  :substeps 50
    ;;  :conv-steps 1000
    ;;  :dt-scale dt-scale)
    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (setf (cl-mpm/particle::mp-enable-plasticity mp) t)))
    (loop for step from 0 below step-count
          while *run-sim*
          do
             (format t "Step ~D~%" step)
             (setf work 0d0)
             (cl-mpm::iterate-over-mps
              (cl-mpm:sim-mps *sim*)
              (lambda (mp)
                (setf (cl-mpm/particle::mp-gravity mp) (+ load-0 (* step load-inc)))))

             (cl-mpm/dynamic-relaxation:converge-quasi-static
              *sim*
              :energy-crit 1d-2
              :oobf-crit 1d-2
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
