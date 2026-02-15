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
    ;; (cl-mpm::g2p (cl-mpm:sim-mesh *sim*)
    ;;              (cl-mpm:sim-mps *sim*)
    ;;              (cl-mpm:sim-dt *sim*)
    ;;              0d0
    ;;              :TRIAL)
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
      (let ((ps-y (sqrt (* E ps-vm))))
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
               (cl-mpm/damage::criterion-mohr-coloumb-rankine-stress-tensile
                ;;cl-mpm/damage::criterion-mohr-coloumb-stress-tensile
                (cl-mpm/fastmaths:fast-.+
                 (cl-mpm/fastmaths::fast-scale-voigt stress (/ 1d0 (magicl:det def)))
                 (cl-mpm/utils:voigt-eye (* 1d0
                                            ;; (* 1d0 (magicl:det def))
                                            (/ (- pressure) 3))))
                (* angle (/ pi 180d0)))
               )))
      ;; (setf y (cl-mpm/damage::tensile-energy-norm strain E de))
      ;; (setf (cl-mpm/particle::mp-damage-y-local mp) )
      )))


(defparameter *angle* 40d0)
(defparameter *angle-r* 10d0)
(defparameter *rc* 0d0)
(defparameter *rs* 1d0)
(defparameter *viscosity* 1d13)
(defparameter *delay-time* 1d3)
(defparameter *delay-exponent* 1d0)
(defparameter *gf* 1000d0)
(defparameter *length-scaler* 2d0)
(defparameter *enable-plastic-damage* nil)

(defun setup (&key (refine 1) (mps 2)
                (pressure-condition t)
                (cryo-static t)
                (hydro-static nil)
                (friction 0d0)
                (ice-height 800d0)
                (bench-length 0d0)
                (aspect 1)
                (floatation-ratio 0.9)
                ;; (extra-cliff-height 0)
                (slope 0d0)
                (use-penalty t)
                (stick-base t)
                (multigrid-refines 0)
                )
  (let* ((density 918d0)
         (water-density 1028d0)
         ;; (density 900d0)
         ;; (water-density 1000d0)
         ;; (refine (/ refine (expt 2 (- refines 1))))
         (mesh-resolution (/ 10d0 refine))
         ;; (mps (* mps (expt 2 (- refines 1))))
         (h-fine mesh-resolution)
         ;; (mesh-resolution h-fine)
         ;; (h-fine mesh-resolution)
         (offset (* mesh-resolution (if use-penalty 2 0)))
         ;; (end-height ice-height)
         ;; (start-height ice-height)

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
                                               ;; 'cl-mpm/damage::mpm-sim-usl-damage
                                               ;; 'cl-mpm/damage::mpm-sim-damage
                                               ;; 'cl-mpm::mpm-sim-usf
                                               'cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul
                                               ;; 'cl-mpm/dynamic-relaxation::mpm-sim-dr-multigrid
                                               ;; 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul-usl
                                               :args-list
                                               (list
                                                :enable-fbar t
                                                :enable-aggregate t
                                                :split-factor (* 1.2d0 (sqrt 2) (/ 1d0 mps))
                                                :enable-split nil
                                                ;; :refinement multigrid-refines
                                                )))
    (setf mesh-resolution (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
    (setf h-fine mesh-resolution)
    (let* ((angle *angle*)
           (init-stress (* 0.1185d6 1d0))
           (init-c (cl-mpm/damage::mohr-coloumb-tensile-to-coheasion init-stress (* angle (/ pi 180))))
           (gf *gf*)
           ;; (length-scale 10d0)
           (length-scale (* mesh-resolution *length-scaler*))
           (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress E))
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
      (let* ((rt 1d0)
             (rc *rc*)
             (rs (est-shear-from-angle angle *angle-r* rc)))
        (cl-mpm:add-mps
         *sim*
         (cl-mpm/setup:make-block-mps
          (list 0 offset)
          block-size
          (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
          density
          'cl-mpm/particle::particle-ice-delayed
          :E 1d9
          :nu 0.24d0

          :kt-res-ratio rt
          :kc-res-ratio rc
          :g-res-ratio rs

          :initiation-stress init-stress;18d3
          :friction-angle angle

          :psi (* 0d0 (/ pi 180))
          :phi (* angle (/ pi 180))
          :c (* init-c oversize)

          :softening 0d0
          :ductility ductility
          :local-length length-scale
          :delay-time *delay-time*
          :delay-exponent *delay-exponent*
          :enable-plasticity t
          :enable-damage t
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
         :k-z 1d0))

      (unless (= start-height end-height)
        (cl-mpm/setup::remove-sdf *sim*
                                  (lambda (p)
                                    (cl-mpm/setup::plane-point-point-sdf
                                     p
                                     (cl-mpm/utils:vector-from-list (list 0d0 (+ offset start-height) 0d0))
                                     (cl-mpm/utils:vector-from-list (list ice-length (+ offset end-height) 0d0))))
                                  :refine 2))


      (let* ((domain-height (second domain-size))
             (midpoint (/ (+ domain-height (+ water-level offset)) 2))
             (dist (- domain-height midpoint))
             (cutout (/ dist 1))
             ;; (cutout (+ (- end-height water-level) 00d0))
             (cutback bench-length)
             )
        ;; (pprint cutout)
        (when (> cutback 0d0)
          (cl-mpm/setup:remove-sdf
           *sim*
           (cl-mpm/setup::rectangle-sdf (list (first block-size)
                                              ;; (+ offset ice-height)
                                              midpoint
                                              )
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
    ;; (setf
    ;;  (cl-mpm:sim-bcs *sim*)
    ;;  (cl-mpm/bc::make-outside-bc-varfix
    ;;   (cl-mpm:sim-mesh *sim*)
    ;;   '(0 nil 0)
    ;;   '(0 nil 0)
    ;;   '(nil 0 nil)
    ;;   '(nil 0 nil)
    ;;   '(nil nil 0)
    ;;   '(nil nil 0)))
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
    (when (typep *sim* 'cl-mpm/damage::mpm-sim-damage)
      (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t))
    (setf (cl-mpm::sim-allow-mp-split *sim*) t)
    (setf (cl-mpm::sim-allow-mp-damage-removal *sim*) nil)
    (setf (cl-mpm::sim-mp-damage-removal-instant *sim*) nil)
    (setf (cl-mpm::sim-mp-damage-removal-criteria *sim*) 0.99d0)
    (setf (cl-mpm::sim-max-split-depth *sim*) 3)
    (setf (cl-mpm::sim-ghost-factor *sim*)
          ;; nil
          (* E 1d-4)
          )
    ;; (setf (cl-mpm::sim-velocity-algorithm *sim*) :PIC)
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
         (* E 0.01d0)
         friction
         ;; 0.1d0
         0d0
         )))

    (when use-penalty
      (cl-mpm:add-bcs-force-list
       *sim*
       *floor-bc*
       ))
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
    ))

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

(defun stability-qt-test ()
  (let* ((heights (list 300d0))
         (floatations (list 0.7d0)))
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
                          (let* ((mps 2)
                                 (output-dir (format nil "./output-~f-~f/" height flotation)))
                            (format t "Problem ~f ~f~%" height flotation)
                            (setup :refine 0.25d0
                                   :multigrid-refines 2
                                   :friction 0d0
                                   :ice-height height
                                   :mps mps
                                   :hydro-static nil
                                   :cryo-static t
                                   :aspect 2d0
                                   :slope 0d0
                                   :bench-length 0d0
                                   :floatation-ratio flotation
                                   :use-penalty nil
                                   :stick-base t
                                   )
                            (plot-domain)
                            (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0)
                            (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) nil)
                            (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
                                  ;; (cl-mpm::sim-ghost-factor *sim*) (* 1d9 1d-3)
                                  (cl-mpm::sim-ghost-factor *sim*) nil
                                  )
                            (cl-mpm/setup::set-mass-filter *sim* 918d0 :proportion 1d-15)
                            (let ((res (cl-mpm/dynamic-relaxation::run-quasi-time
                                        *sim*
                                        :output-dir output-dir
                                        :dt 1d4
                                        :total-time 1d6
                                        ;; :steps 1000
                                        :dt-scale 1d0
                                        :conv-criteria 1d-3
                                        :substeps 10
                                        :enable-damage t
                                        :enable-plastic t
                                        :min-adaptive-steps -4
                                        :max-adaptive-steps 9
                                        :save-vtk-dr t
                                        :save-vtk-loadstep t
                                        ;; :elastic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
                                        ;; :elastic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-dr-multigrid
                                        :plotter (lambda (sim) (plot-domain))
                                        :post-conv-step (lambda (sim) (plot-domain)
                                                          (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)))))
                              (cl-mpm/dynamic-relaxation::save-vtks *sim* output-dir 1)
                              (setf (aref *stability* hi fi) (if res 1 0))
                              (unless res
                                (loop for j from fi below (length floatations)
                                      do (setf (aref *stability* hi j) 0)))
                              (format t "Stability:~%")
                              (print-stab)
                              (save-stabilty-data stability-dir *sim* res height flotation 0d0)
                              (unless res
                                (loop-finish)))
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
