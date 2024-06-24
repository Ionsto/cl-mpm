(defpackage :cl-mpm/examples/shear-box
  (:use :cl))
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)

(in-package :cl-mpm/examples/shear-box)
(declaim (optimize (debug 3) (safety 3) (speed 0)))
(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-chalk-brittle) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
                     (strain cl-mpm/particle::mp-strain)
                     (damage cl-mpm/particle:mp-damage)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (pressure cl-mpm/particle::mp-pressure)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     (angle cl-mpm/particle::mp-friction-angle)
                     (c cl-mpm/particle::mp-coheasion)
                     (nu cl-mpm/particle::mp-nu)
                     (ft cl-mpm/particle::mp-ft)
                     (fc cl-mpm/particle::mp-fc)
                     (E cl-mpm/particle::mp-e)
                     (de cl-mpm/particle::mp-elastic-matrix)
                     (kc-r cl-mpm/particle::mp-k-compressive-residual-ratio)
                     (kt-r cl-mpm/particle::mp-k-tensile-residual-ratio)
                     (g-r cl-mpm/particle::mp-shear-residual-ratio)
                     ) mp
      (declare (double-float pressure damage))
      (progn
        (when (< damage 1d0)
          ;; (setf damage-increment (cl-mpm/damage::tensile-energy-norm strain E de))
          ;; (setf damage-increment
          ;;       (sqrt (* 3/2
          ;;                (cl-mpm/constitutive::voigt-j2
          ;;                 (cl-mpm/utils::deviatoric-voigt (magicl:scale stress (/ 1d0 (magicl:det def))))))))
          ;; (setf damage-increment (* E (cl-mpm/damage::modified-vm-criterion strain nu (/ fc ft))))

          ;; (multiple-value-bind (s_1 s_2 s_3) (cl-mpm/damage::principal-stresses-3d (magicl:scale stress (/ 1d0 (magicl:det def))))
          ;;   (setf damage-increment (max damage-increment
          ;;                               (max s_1
          ;;                                    s_2
          ;;                                    s_3
          ;;                                    0d0
          ;;                                    )))
          ;;   )
          ;(setf damage-increment (criterion-mc strain (* angle (/ pi 180d0)) E nu))
          ;; (setf damage-increment (criterion-dp strain (* angle (/ pi 180d0)) E nu))
          ;; (setf damage-increment (cl-mpm/damage::smooth-rankine-criterion (magicl:scale stress (/ 1d0 (magicl:det def)))))
          (setf damage-increment
                (max 0d0
                     (cl-mpm/damage::drucker-prager-criterion
                      (magicl:scale stress (/ 1d0 (magicl:det def))) (* angle (/ pi 180d0)))))
          )
        (when (>= damage 1d0)
          (setf damage-increment 0d0))
        ;;Delocalisation switch
        (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
        (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
        ))))

(defun plot (sim)
  (plot-load-disp))
(defun plot-domain (sim)

  ;; (vgplot:plot (mapcar (lambda (x) (* x 1d3)) *data-disp*) *data-v*)
  (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms)))
    (vgplot:format-plot t "set object 1 rect from 0,~f to ~f,~f fc rgb 'black' fs transparent solid 0.8 noborder behind"
                        (* 0.5d0 *box-size*)
                        (+ *box-size* *displacement-increment*)
                        ms-y)
    (vgplot:format-plot t "set object 2 rect from 0,0 to ~f,~f fc rgb 'red' fs transparent solid 0.8 noborder behind"
                        (* 1d0 *box-size*)
                        (* 0.5d0 *box-size*)
                        )
    (vgplot:format-plot t "set object 3 rect from ~f,0 to ~f,~f fc rgb 'red' fs transparent solid 0.8 noborder behind"
                        (* 2d0 *box-size*)
                        ms-x
                        (* 0.5d0 *box-size*))
    (vgplot:format-plot t "set object 4 rect from ~f,~f to ~f,~f fc rgb 'black' fs transparent solid 0.8 noborder behind"
                        (+ (* 2 *box-size*) *displacement-increment*)
                        (* 0.5d0 *box-size*)
                        ms-x
                        ms-y)
    )
  (vgplot:format-plot t "set style fill solid")
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :deformed
   ;:colour-func (lambda (mp) (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0))
   ;; :colour-func (lambda (mp) (magicl:tref (cl-mpm/particle::mp-penalty-frictional-force mp) 0 0))
   ;; :colour-func #'cl-mpm/particle::mp-index
   ;:colour-func #'cl-mpm/particle::mp-damage
   :colour-func 
   (lambda (mp)
     (if (= 0 (cl-mpm/particle::mp-index mp))
         (cl-mpm/particle::mp-strain-plastic-vm mp)
         0d0)))
  )
(defun plot-load-disp ()

  ;; (vgplot:figure)
  (let ((width 0.06d0)
        (depth (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
    (vgplot:plot (mapcar (lambda (x) (* x 1d3)) *data-disp*) 
                 (mapcar (lambda (x) (* depth (/ x width) 1d-3)) *data-v*)))
  (vgplot:xlabel "Displacement (mm)")
  (vgplot:ylabel "Shear stress (kN/m^2)")
  )
 
(defun apply-penalty-box (left-x right-x height friction)
  (let* ((left-normal (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0)))
         (right-normal (cl-mpm/utils:vector-from-list (list -1d0 0d0 0d0)))
         (plane-normal (cl-mpm/utils:vector-from-list (list 0d0 1d0 0d0)))
         (plane-normal-left (cl-mpm/utils:vector-from-list (list 0d0 -1d0 0d0)))
         (epsilon (* 1d2 1d9))
         ;; (friction 0.0d0)
         )
    (with-accessors ((mesh cl-mpm:sim-mesh)
                     (mps cl-mpm:sim-mps)
                     (dt cl-mpm:sim-dt))
        *sim*
      (cl-mpm/penalty::apply-force-mps
       mesh mps dt
       plane-normal
       height
       epsilon
       friction
       (lambda (mp) (and
                     ;; (>= (magicl:tref mp 1 0) (- height 1d-3))
                     (>= (magicl:tref mp 0 0) right-x)))
       )
      (cl-mpm/penalty::apply-force-mps
       mesh mps dt
       plane-normal-left
       (- height)
       epsilon
       friction
       (lambda (mp) (and
                     ;; (<= (magicl:tref mp 1 0) (+ height 1d-3))
                     (<= (magicl:tref mp 0 0) (+ left-x *displacement-increment*))))
       )
      (cl-mpm/penalty::apply-force-mps
       mesh mps dt
       right-normal
       (- right-x)
       epsilon
       friction
       (lambda (mp) (<= (magicl:tref mp 1 0) height)))
      (cl-mpm/penalty::apply-force-mps
       mesh mps dt
       left-normal
       left-x
       epsilon
       friction
       (lambda (mp) (<= (magicl:tref mp 1 0) height)))
      (setf cl-mpm/penalty::*debug-force* 0d0)
      (cl-mpm/penalty::apply-force-mps
       mesh mps dt
       right-normal
       (- (+ right-x *displacement-increment*))
       epsilon
       friction
       (lambda (mp) (> (magicl:tref mp 1 0) height)))
      (setf cl-mpm/penalty::*debug-force* 0d0)
      (cl-mpm/penalty::apply-force-mps
       mesh mps dt
       left-normal
       (+ left-x *displacement-increment*)
       epsilon
       friction
       (lambda (mp) (> (magicl:tref mp 1 0) height))))

    ))
(defparameter *displacement-increment* 0d0)
(defparameter *box-size* 0d0)
(defun apply-loading-force (left-x right-x height)
  (let* ((left-normal (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0)))
         (right-normal (cl-mpm/utils:vector-from-list (list -1d0 0d0 0d0)))
         (epsilon 1d9)
         (friction 0.0d0))
    (with-accessors ((mesh cl-mpm:sim-mesh)
                     (mps cl-mpm:sim-mps)
                     (dt cl-mpm:sim-dt))
        *sim*
      (cl-mpm/penalty::apply-force-mps
       mesh mps dt
       right-normal
       (- (+ right-x *displacement-increment*))
       epsilon
       friction
       (lambda (mp) (> (magicl:tref mp 1 0) height)))
      (setf cl-mpm/penalty::*debug-force* 0d0)
      (cl-mpm/penalty::apply-force-mps
       mesh mps dt
       left-normal
       (+ left-x *displacement-increment*)
       epsilon
       friction
       (lambda (mp) (> (magicl:tref mp 1 0) height))))
    ))

(declaim (notinline setup-test-column))
(defun setup-test-column (size offset block-size &optional (e-scale 1) (mp-scale 1) &key (angle 0d0) (friction 0.1d0) (surcharge-load 72.5d3))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               :sim-type 'cl-mpm::mpm-sim-usf
               ;; :sim-type 'cl-mpm/damage::mpm-sim-damage
               ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         ;; (floor-offset (* h-y 2))
         (floor-offset 2d0)
         (density 1.7d3)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
         )
    (declare (double-float h density))
    (progn
      (let* ((angle-rad (* angle (/ pi 180)))
             (init-stress 60d3)
             (downscale (/ 1d0 1d0))
                                        ;(gf (/ (expt (/ init-stress 6.88d0) 2) 1d9))
             (gf 4d0)
             (length-scale h)
             ;; (length-scale (/ (* 1d9 gf) (expt init-stress 2)))
             (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress 1d9))
             )
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-mps-from-list
               (cl-mpm/setup::make-block-mps-list
                offset
                block-size
                (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                density
                ;; 'cl-mpm/particle::particle-elastic
                'cl-mpm/particle::particle-mc
                :E 1d9
                :nu 0.24d0
                ;; :rho 20d3

                ;; 'cl-mpm/particle::particle-chalk-delayed
                ;; :E 1d9
                ;; :nu 0.24d0

                ;; :enable-plasticity t

                ;; :kt-res-ratio 1d-10
                ;; :kc-res-ratio 1d0
                ;; :g-res-ratio 1d-2

                ;; :friction-angle 43d0
                ;; :initiation-stress init-stress;18d3
                ;; :delay-time 1d-3
                ;; :delay-exponent 1d0
                ;; ;; :ductility 5d0
                ;; :ductility ductility
                ;; :local-length length-scale
                ;; :local-length-damaged 10d-10
                :enable-plasticity t
                :psi 0d0
                ;; :phi (* 42d0 (/ pi 180))
                ;; :c 131d3
                :phi (* 30d0 (/ pi 180))
                :c 0d3

                :index 0
                :gravity 0.0d0
                )))
        )
      (let* ((sur-height (* 1 h-x))
             (sur-size (list 0.06d0 sur-height))
             ;(load 72.5d3)
             (load surcharge-load)
             ;; (gravity 10d0)
             ;; (density (/ load (* gravity sur-height)))
             (gravity (/ load (* density sur-height)))
             )
        (format t "Gravity ~F~%" gravity)
        (cl-mpm::add-mps
         sim
         (cl-mpm/setup::make-mps-from-list
          (cl-mpm/setup::make-block-mps-list
           (mapcar #'+ offset (list 0d0 (second block-size)))
           sur-size
           (mapcar (lambda (e) (* e e-scale mp-scale)) sur-size)
           density
           'cl-mpm/particle::particle-elastic-damage
           :E 1d9
           :nu 0.24d0
           :initiation-stress 1d20
           :index 1
           :gravity (- gravity)))))
      (setf (cl-mpm:sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-enable-fbar sim) t)
      ;; (setf (cl-mpm::sim-mass-filter sim) 1d0)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (let ((ms 1d0))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        (setf (cl-mpm:sim-damping-factor sim)
              (* 0.1d0 (cl-mpm/setup::estimate-critical-damping sim))))

      (let ((dt-scale 1d0))
        (setf
         (cl-mpm:sim-dt sim)
         (* dt-scale h
            (sqrt (cl-mpm::sim-mass-scale sim))
            (sqrt (/ density (cl-mpm/particle::mp-p-modulus (aref (cl-mpm:sim-mps sim) 0)))))))

      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))
      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-var
             (cl-mpm:sim-mesh sim)
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil 0 nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 0 nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
            ))
      ;; (defparameter *floor-bc*
      ;;   (cl-mpm/penalty::make-bc-penalty-point-normal
      ;;    sim
      ;;    (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
      ;;    (cl-mpm/utils:vector-from-list (list 0d0 floor-offset 0d0))
      ;;    (* 1d9 1d0)
      ;;    friction)
      ;;   )
      ;; (setf (cl-mpm::sim-bcs-force-list sim)
      ;;       (list
      ;;        (cl-mpm/bc:make-bcs-from-list
      ;;         (list
      ;;          *floor-bc*))))
      sim)))


(defparameter *sim* nil)
(defparameter *run-sim* t)
(defparameter *t* 0)
(defparameter *sim-step* 0)
(defparameter *refine* (/ 1d0 1d0))
(let ((refine (uiop:getenv "REFINE")))
  (when refine
    (setf *refine* (parse-integer (uiop:getenv "REFINE")))  
    ))
(defmethod cl-mpm::post-stress-step (mesh (mp cl-mpm/particle::particle-vm) dt)
  (with-accessors ((ps cl-mpm/particle::mp-strain-plastic-vm)
                   (rho cl-mpm/particle::mp-rho))
      mp
    (let* ((rho-0 200d3)
           (rho-1 200d2)
           (soft 100d0)
           )
      (setf rho (+ rho-1 (* (- rho-0 rho-1) (exp (- (* soft ps)))))))))

(defmethod cl-mpm::post-stress-step (mesh (mp cl-mpm/particle::particle-mc) dt)
  (with-accessors ((ps cl-mpm/particle::mp-strain-plastic-vm)
                   (c cl-mpm/particle::mp-c)
                   (phi cl-mpm/particle::mp-phi)
                   )
      mp
    (let ((phi_0 (* 42d0 (/ pi 180)))
          (phi_1 (* 30d0 (/ pi 180)))
          (c_0 131d3)
          (soft 20d0))
      ;; (setf
      ;;  c (* c_0 (exp (- (* soft ps))))
      ;;  phi (+ phi_1 (* (- phi_0 phi_1) (exp (- (* soft ps))))))
      )))

(defun setup (&key (refine 1d0) (mps 2) (friction 0.0d0) (surcharge-load 72.5d3))
  (defparameter *displacement-increment* 0d0)
  (let* ((mps-per-dim mps)
         (mesh-size (/ 0.03d0 refine))
         (sunk-size 0.03d0)
         (box-size (* 2d0 sunk-size))
         (domain-size (* 3d0 box-size))
         (offset (list box-size 0))
         )
    (setf *box-size* box-size)
    (defparameter *sim* (setup-test-column
                         (list domain-size (* 2 box-size))
                         offset
                         (list box-size box-size)
                         (/ 1d0 mesh-size)
                         mps-per-dim
                         :friction friction
                         :surcharge-load surcharge-load))
    (setf (cl-mpm::sim-bcs-force-list *sim*)
          (list
           (cl-mpm/bc:make-bcs-from-list
            (list
             (cl-mpm/bc::make-bc-closure
              nil
              (lambda ()
                (apply-penalty-box box-size (* 2d0 box-size) sunk-size friction)
                ;; (apply-loading-force box-size (* 2d0 box-size) sunk-size)
                )))))))
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (format t "Mesh-size: ~E~%" (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defun stop ()
  (setf *run-sim* nil
        cl-mpm/dynamic-relaxation::*run-convergance* nil))

(defun run (&optional (output-directory "./output/"))
  (cl-mpm/output:save-vtk-mesh (merge-pathnames output-directory "mesh.vtk")
                          *sim*)

  (defparameter *data-t* (list))
  (defparameter *data-disp* (list))
  (defparameter *data-v* (list))
  (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :supersede)
    (format stream "disp,load~%"))
  (vgplot:close-all-plots)
  (let* ((target-time 0.0001d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         ;; (substeps 10)
         (dt-scale 0.8d0)
         (load-steps 100)
         (displacment 1d-3)
         ;; (displacment 1d-5)
         (disp-inc (/ displacment load-steps)))

    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (cl-mpm::update-sim *sim*)
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                    (format t "CFL dt estimate: ~f~%" dt-e)
                    (format t "CFL step count estimate: ~D~%" substeps-e)
                    (setf substeps substeps-e))
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 10d0 (cl-mpm/setup::estimate-critical-damping *sim*)))

    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (when (= (cl-mpm/particle::mp-index mp) 0)
               (setf (cl-mpm/particle::mp-enable-plasticity mp) nil)))

    (cl-mpm/dynamic-relaxation:converge-quasi-static
     *sim*
     :energy-crit 1d-2
     :oobf-crit 1d-2
     :substeps 20
     :conv-steps 200
     :post-iter-step
     (lambda (i energy oobf)
       (cl-mpm/output:save-vtk (merge-pathnames output-directory (format nil "sim_conv_~5,'0d.vtk" i)) *sim*)))

    ;; (loop for steps from 0 below load-steps
    ;;       while *run-sim*
    ;;       do
    ;;          (progn
    ;;            (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_grav_~5,'0d.vtk" steps)) *sim*)
    ;;           (dotimes (i 10)
    ;;             (cl-mpm::update-sim *sim*))
    ;;           (format t "Energy ~F~%" (cl-mpm/dynamic-relaxation:estimate-energy-norm *sim*))
    ;;           (swank.live:update-swank)))

    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (when (= (cl-mpm/particle::mp-index mp) 0)
               (setf (cl-mpm/particle::mp-enable-plasticity mp) nil)))

    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 1d-1 (cl-mpm/setup::estimate-critical-damping *sim*))
          ;; (cl-mpm::sim-enable-damage *sim*) t
          )
    (defparameter *displacement-increment* 0d0)
    (format t "Substeps ~D~%" substeps)
    (vgplot:figure)
    (setf cl-mpm/penalty::*debug-force* 0)
    (time (loop for steps from 0 below load-steps
                while *run-sim*
                do
                   (progn
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames output-directory (format nil "sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-nodes (merge-pathnames output-directory (format nil "sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     (let ((load-av 0d0)
                           (disp-av 0d0))
                       (time
                        (dotimes (i substeps)
                          (cl-mpm::update-sim *sim*)
                          (incf load-av (/ cl-mpm/penalty::*debug-force* substeps))
                          (incf disp-av (/ *displacement-increment* substeps))
                          (incf *displacement-increment* (/ disp-inc substeps))
                          (incf *t* (cl-mpm::sim-dt *sim*))))

                       (setf load-av cl-mpm/penalty::*debug-force*)
                       (setf disp-av *displacement-increment*)

                       (push *t* *data-t*)
                       (push disp-av *data-disp*)
                       (push load-av *data-v*)
                       (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :supersede)
                         (format stream "~f,~f~%" disp-av load-av)))
                     (incf *sim-step*)
                     (plot *sim*)
                     (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                        :terminal "png size 1920,1080"
                                        )
                     (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       (format t "CFL dt estimate: ~f~%" dt-e)
                       (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf substeps substeps-e))
                     (swank.live:update-swank))))
    )
  (vgplot:figure)
  (plot-load-disp)
  )
