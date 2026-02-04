(defpackage :cl-mpm/examples/collapse
  (:use :cl))

(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)
;; (sb-ext:restrict-compiler-policy 'speed  3 3)
;; (sb-ext:restrict-compiler-policy 'debug  0 0)
;; (sb-ext:restrict-compiler-policy 'safety 0 0)
;; (setf *block-compile-default* nil)
;(sb-int:set-floating-point-modes :traps '(:overflow :invalid :inexact :divide-by-zero :underflow))
;; (sb-int:set-floating-point-modes :traps '(:overflow :divide-by-zero :underflow))

(in-package :cl-mpm/examples/collapse)
;; (declaim (optimize (debug 0) (safety 0) (speed 3)))
(declaim (optimize (debug 3) (safety 3) (speed 0)))

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-elastic) dt fbar)
  (cl-mpm::update-stress-kirchoff-dynamic-relaxation mesh mp dt fbar))

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-elastic) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  (cl-mpm::update-domain-polar-2d mesh mp dt))

(defun plot-load-disp ()
  (vgplot:semilogy *data-steps* *data-energy*))
(declaim (notinline plot))
(defun plot (sim)
  ;; (cl-mpm::g2p (cl-mpm:sim-mesh *sim*)
  ;;              (cl-mpm:sim-mps *sim*)
  ;;              (cl-mpm:sim-dt *sim*)
  ;;              0d0
  ;;              :TRIAL)
  ;; (plot-load-disp)
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
         (h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
         (ms-x (first ms))
         (ms-y (second ms))
         (datum (cl-mpm/penalty::bc-penalty-datum *bc-squish*))
        )
    ;; (vgplot:format-plot t "set object 1 rect from ~f,~f to ~f,~f fc rgb 'black' fs transparent solid 0.5 noborder behind" 0 ms-y ms-x (- datum))
    )
  ;; (format t "~E~%" (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :deformed
   :trial t
   :colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :xy))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-damage-ybar mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-strain-plastic-vm mp))
   )
  )

(defparameter *eta* 1d0)

(defun stop ()
  (setf *run-sim* nil))

(declaim (notinline setup-test-column))
(defun setup-test-column (size block-size sim-type &optional (e-scale 1) (mp-scale 1) (multigrid-refinement 0))
  (let* (;; (multigrid-refinement 3)
         (multigrid-enabled (> multigrid-refinement 0))
         (e-scale (if multigrid-enabled (/ e-scale (expt 2 (- multigrid-refinement 1)))
                      e-scale))
         (mp-scale (if multigrid-enabled (* mp-scale (expt 2 (- multigrid-refinement 1)))
                       mp-scale))
         (sim (cl-mpm/setup::make-simple-sim
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               ;; :sim-type
               ;; 'cl-mpm::mpm-sim-sd
               ;; :sim-type sim-type
               ;; :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-usf
               ;; :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
               :sim-type (if multigrid-enabled
                             'cl-mpm/dynamic-relaxation::mpm-sim-dr-multigrid
                             'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul)
               ;; :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
               ;; :sim-type 'cl-mpm/aggregate:mpm-sim-agg-usf
               ;; :sim-type 'cl-mpm/aggregate:mpm-sim-agg-usf
               ;; :sim-type 'cl-mpm/damage::mpm-sim-agg-damage
               ;; :sim-type 'cl-mpm/damage::mpm-sim-agg-damage
               :args-list
               (append
                (list
                 :split-factor 0.8d0
                 :enable-fbar t
                 :enable-aggregate t
                 :vel-algo :QUASI-STATIC
                 ;; :vel-algo :BLEND
                 :max-split-depth 6
                 ;; :enable-length-localisation nil
                 :enable-split nil
                 :gravity -10d0
                 )
                (when multigrid-enabled
                  (list :refinement multigrid-refinement))
                )))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 1d3)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
         (offset (* h 0)))
    (declare (double-float h density))
    (progn
      (let* ((angle 40d0)
             (E 1d9))
        (cl-mpm::add-mps
         sim
         (cl-mpm/setup::make-block-mps
          (mapcar #'+
                  (mapcar (lambda (x) 0) size)
                  (list 0 offset 0))
          block-size
          (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
          density
          ;; 'cl-mpm/particle::particle-elastic
          ;; :E 1d6
          ;; :nu 0.30d0
          ;; 'cl-mpm/particle::particle-vm
          ;; :E 1d6
          ;; :nu 0.3d0
          ;; :rho 20d3
          'cl-mpm/particle::particle-mc
          :E 1d6
          :nu 0.3d0
          :phi (* 30d0 (/ 3.14d0 180d0))
          :psi (* 00d0 (/ 3.14d0 180d0))
          :c 10d3
          ;; :enable-plasticity t
          ;; 'cl-mpm/particle::particle-elastic-damage-delayed
          ;; :E 1d9
          ;; :nu 0.3d0
          ;; :initiation-stress 1d4
          ;; :local-length (* 2 h)
          ;; :ductility 50d0
          ;; :delay-time 100d0
          ;; :delay-exponent 2d0

          ;; 'cl-mpm/particle::particle-chalk-erodable
          ;; :E E
          ;; :nu 0.24d0
          ;; :enable-damage t
          ;; :enable-plasticity t
          ;; :friction-angle angle

          ;; :kt-res-ratio 1d0
          ;; :kc-res-ratio 0d0
          ;; :g-res-ratio 0.51d0

          ;; :initiation-stress 1d3
          ;; :delay-time 1d0
          ;; :delay-exponent 1d0
          ;; :ductility ductility
          ;; :local-length h
          ;; :psi (* 0d0 (/ pi 180))
          ;; :phi (* angle (/ pi 180))
          ;; :c (* init-c oversize)
          ;; :softening 0d0
          ;; :gravity -10d0
          ;; :gravity-axis (cl-mpm/utils:vector-from-list '(-1d0 1d0 0d0))
          )))
      ;; (format t "Charictoristic time ~E~%" (/ ))
      ;; (setf (cl-mpm:sim-allow-mp-split sim) t)
      (setf (cl-mpm::sim-split-factor sim) 0.4d0)
      ;; (setf (cl-mpm::sim-max-split-depth sim) 6)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      ;; (setf (cl-mpm::sim-enable-fbar sim) t)
      (defparameter *density* density)
      (cl-mpm/setup::set-mass-filter sim density :proportion 1d-9)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) t)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)

      (when (typep sim 'cl-mpm/damage::mpm-sim-damage)
        (setf (cl-mpm/damage::sim-enable-length-localisation sim) t))

      (let ((ms 1d0))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        (setf (cl-mpm:sim-damping-factor sim)
              (* 0d-1
                 (cl-mpm/setup:estimate-critical-damping sim))))

      (setf *eta* (* 0.5d0 (cl-mpm/setup::estimate-stiffness-critical-damping sim 1d6 density)))

      (let ((dt-scale 1d0))
        (setf
         (cl-mpm:sim-dt sim)
         (* dt-scale h
            (sqrt (cl-mpm::sim-mass-scale sim))
            (sqrt (/ density (cl-mpm/particle::mp-p-modulus (aref (cl-mpm:sim-mps sim) 0)))))))

      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))

      (cl-mpm/setup::setup-bcs
       sim
       :top '(nil nil nil)
       :bottom '(nil 0 nil)
       ;; :left '(0 nil nil)
       ;; :front '(nil nil 0)
       )


      ;; (let ((water-density 1000d0)
      ;;       (datum 3d0))
      ;;   (defparameter *water-bc*
      ;;     (cl-mpm/buoyancy::make-bc-buoyancy-clip
      ;;      sim
      ;;      datum
      ;;      water-density
      ;;      (lambda (pos datum)
      ;;        (>= (cl-mpm/utils:varef pos 1) (* h 0)))
      ;;      :visc-damping 0d0))
      ;;   (cl-mpm:add-bcs-force-list
      ;;    sim
      ;;    *water-bc*))
      (defparameter *bc-squish*
        (cl-mpm/penalty::make-bc-penalty-distance-point
         sim
         (cl-mpm/utils:vector-from-list (list 0d0 -1d0 0d0))
         (cl-mpm/utils:vector-from-list (list (* (first size) 0.5d0) (float (second block-size) 0d0) 0d0))
         (* (first size) 0.5d0)
         (* 1d6 1d0)
         0d0
         0d0
         ))
      (defparameter *bc-pressure*
        (cl-mpm/buoyancy::make-bc-pressure
         sim
         0d0
         0d0)
        )

      (let ((domain-half (* 0.5d0 (first size)))
            (friction 0.5d0)
            (E 1d6))
        (defparameter *floor-bc*
          (cl-mpm/penalty::make-bc-penalty-distance-point
           sim
           (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
           (cl-mpm/utils:vector-from-list (list
                                           domain-half
                                           offset
                                           0d0))
           (* domain-half 1d0)
           (* E 1d0)
           friction
           1d0)))
      (when (> offset 0d0)
        (cl-mpm:add-bcs-force-list sim *floor-bc*))
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

(declaim (notinline setup))
(defun setup (&key (refine 1) (mps 2)
                (sim-type 'cl-mpm:mpm-sim-usf)
                (multigrid-refine 0))
  (defparameter *sim* nil)
  (let ((mps-per-dim mps))
    (defparameter *sim* (setup-test-column '(32 16) '(8 8) sim-type refine mps-per-dim multigrid-refine)))
  ;; (cl-mpm/setup::initialise-stress-self-weight
  ;;  *sim*
  ;;  15d0
  ;;  )
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defun save-csv (output-file)
  (let* ((mp-list
           (loop for mp across (cl-mpm:sim-mps *sim*)
                 collect mp))
         (y (loop for mp in mp-list collect (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)))
         (y-ref (loop for pos in *original-configuration*
                      collect (float pos 1e0)))
         (syy (loop for mp in mp-list collect (float (magicl:tref (cl-mpm/particle::mp-stress mp) 1 0) 1e0)))
         (rho 80d0)
         (E 1d5)
         (g 10d0)
         (vp-0-list (loop for size in *original-size*
                          collect (float (* (cl-mpm/utils:varef size 0) (cl-mpm/utils:varef size 1)) 1e0)))
         (pressure -1e4)
         (max-y 50)
         (syy-ref (mapcar (lambda (x) pressure) y-ref))
         (df (lisp-stat:make-df '(:y :syy :syy-ref :vp)
                                (list
                                 (coerce y-ref '(vector single-float))
                                 (coerce syy 'vector)
                                 (coerce syy-ref 'vector)
                                 (coerce vp-0-list 'vector)))))
    (lisp-stat:write-csv df output-file :add-first-row t)))

(defun run (&key (output-dir "./output/")
              (ms 1d0)
              (dt-scale 0.5d0)
              )
  (uiop:ensure-all-directories-exist (list output-dir))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (cl-mpm/output:save-vtk-mesh (merge-pathnames output-dir "mesh.vtk")
                          *sim*)

  (defparameter *data-steps* (list))
  (defparameter *data-oobf* (list))
  (defparameter *data-energy* (list))
  (cl-mpm/dynamic-relaxation::elastic-static-solution
   *sim*
   :crit 1d-3
   )
  ;; (setf (cl-mpm:sim-damping-factor *sim*)
  ;;       (* 1d0
  ;;          (cl-mpm/setup:estimate-critical-damping *sim*)))
  ;; (cl-mpm:iterate-over-mps
  ;;  (cl-mpm:sim-mps *sim*)
  ;;  (lambda (mp)
  ;;    (when (typep mp 'cl-mpm/particle::particle-plastic)
  ;;      (setf (cl-mpm/particle::mp-enable-plasticity mp) nil))))
  ;; (cl-mpm/dynamic-relaxation:converge-quasi-static
  ;;  *sim*
  ;;  :dt-scale dt-scale
  ;;  :kinetic-damping t
  ;;  :conv-steps 1000
  ;;  :oobf-crit 1d-3
  ;;  :energy-crit 1d-3
  ;;  :damping-factor nil
  ;;  :post-iter-step
  ;;  (lambda (i energy oobf)
  ;;    (plot *sim*)
  ;;    (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_conv_~5,'0d_~5,'0d.vtk" 0 i)) *sim*)
  ;;    (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_conv_nodes_~5,'0d_~5,'0d.vtk" 0 i)) *sim*)))
  (cl-mpm:iterate-over-mps
   (cl-mpm:sim-mps *sim*)
   (lambda (mp)
     (when (typep mp 'cl-mpm/particle::particle-plastic)
       (setf (cl-mpm/particle::mp-enable-plasticity mp) t))))
  (let* ((target-time 0.1d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-min (cl-mpm:sim-dt *sim*)))
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 0.1d0
             (cl-mpm/setup:estimate-critical-damping *sim*)))
    (setf (cl-mpm:sim-mass-scale *sim*) ms)
    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup:estimate-elastic-dt *sim*))
    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (setf substeps (round target-time (cl-mpm:sim-dt *sim*)))
    (cl-mpm/dynamic-relaxation::save-timestep-preamble output-dir)
    (cl-mpm/output::save-simulation-parameters (merge-pathnames "settings.json" output-dir)
                                               *sim*
                                               (list :dt target-time))
    (format t "Substeps ~D~%" substeps)
    (setf (cl-mpm::sim-stats-work *sim*) 0d0)
    (time (loop for steps from 0 to 1000
                while *run-sim*
                do
                   (let ((energy 0d0)
                         (work 0d0)
                         (oobf 0d0))
                     (cl-mpm/dynamic-relaxation::save-timestep *sim* output-dir steps :DYNAMIC)
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d.vtk" *sim-step*)) *sim*)
                     ;; (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim_p_~5,'0d.vtk" *sim-step*)) *sim* )
                     ;; (cl-mpm/output::save-vtk-mesh-nodes (merge-pathnames output-dir (format nil "sim_nodes_p_~5,'0d.vtk" *sim-step*)) (cl-mpm::sim-mesh-p *sim*))
                     (setf dt-min (cl-mpm::calculate-min-dt *sim*))
                     (time
                      (dotimes (i substeps)
                        (cl-mpm::update-sim *sim*)
                        ;; (setf dt-min (min (cl-mpm::calculate-min-dt *sim*) dt-min))
                        ;; (setf (cl-mpm:sim-dt *sim*) (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
                        (incf oobf (cl-mpm::sim-stats-oobf *sim*))
                        (incf energy (cl-mpm::sim-stats-energy *sim*))
                        (incf work (cl-mpm::sim-stats-power *sim*))
                        (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))

                     (setf
                      energy (/ energy substeps)
                      oobf (/ oobf substeps))

                     (if (= work 0d0)
                         (setf energy 0d0)
                         (setf energy (abs (/ energy work))))

                     ;; (format t "Min dt est ~f~%" dt-min)
                     (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       (format t "CFL dt estimate: ~f~%" dt-e)
                       (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf substeps substeps-e))
                     ;; (setf (cl-mpm:sim-dt *sim*) (* dt-min dt-scale))
                     ;; (setf substeps (floor target-time (cl-mpm:sim-dt *sim*)))
                     (format t "Substeps ~D~%" substeps)
                     ;; (let ((energy 0d0))
                     ;;   (setf energy (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
                     ;;   (push energy *data-energy*)
                     ;;   (push steps *data-steps*)
                     ;;   (format t "Energy ~F~%" energy))
                     (format t "~E ~E~%" energy oobf)

                     (incf *sim-step*)
                     (plot *sim*)
                     (vgplot:title (format nil "Time:~F - KE ~E - OOBF ~E - Work ~E"  (cl-mpm::sim-time *sim*) energy oobf work))
                     ;; (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                     ;;                    :terminal "png size 1920,1080"
                     ;;                    )

                     (swank.live:update-swank)
                     ))))
  (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" *sim-step*)) *sim*))


;; (setf lparallel:*kernel* (lparallel:make-kernel 8 :name "custom-kernel"))
;; (push (lambda ()
;;         (format t "Closing kernel~%")
;;         (lparallel:end-kernel))
;;       sb-ext:*exit-hooks*)
;; (setup)
;; (run)

;; (time
;;  (dotimes (i 1)
;;    (cl-mpm::update-stress (cl-mpm:sim-mesh *sim*)
;;                           (cl-mpm:sim-mps *sim*)
;;                           (cl-mpm:sim-dt *sim*))))

(defun save-data (name)
  (with-open-file (stream (merge-pathnames name) :direction :output :if-exists :supersede)
    (format stream "time,oobf,energy,damage~%")
    (loop for time in (reverse *data-t*)
          for oobf in (reverse *data-oobf*)
          for energy in (reverse *data-energy*)
          for damage in (reverse *data-damage*)
          do (format stream "~f,~f,~f,~f~%" time oobf energy damage))))


(defmacro time-form (it form)
  `(progn
     (declaim (optimize speed))
     (let* ((iterations ,it)
            (start (get-internal-real-time)))
       (time
        (dotimes (i ,it)
              ,form))
       (let* ((end (get-internal-real-time))
              (units internal-time-units-per-second)
              (dt (/ (- end start) (* iterations units)))
              )
         (format t "Total time: ~f ~%" (/ (- end start) units)) (format t "Time per iteration: ~f~%" (/ (- end start) (* iterations units)))
         (format t "Throughput: ~f~%" (/ 1 dt))
         (format t "Time per MP: ~E~%" (/ dt (length (cl-mpm:sim-mps *sim*))))
         dt))))
(defun test-mult ()
  (time
   (cl-mpm:iterate-over-mps
    (cl-mpm:sim-mps *sim*)
    (lambda (mp)
      (with-accessors ((strain cl-mpm/particle:mp-strain)
                       (de cl-mpm/particle::mp-elastic-matrix)
                       (stress cl-mpm/particle::mp-stress-kirchoff))
          mp
        (cl-mpm/constitutive::linear-elastic-mat strain de stress))))))
(defun profile ()
  (setup :refine 16)
  ;; (sb-profile:unprofile)
  ;; (sb-profile:profile "CL-MPM")
  ;; ;; (sb-profile:profile "CL-MPM/PARTICLE")
  ;; ;; (sb-profile:profile "CL-MPM/MESH")
  ;; ;; (sb-profile:profile "CL-MPM/SHAPE-FUNCTION")
  ;; (sb-profile:reset)
  (time-form 100
             (progn
               (format t "~D~%" i)
                  (cl-mpm::update-sim *sim*)))
  ;; (time-form 1000000
  ;;            (progn
  ;;              ;; (cl-mpm:iterate-over-neighbours (cl-mpm:sim-mesh *sim*) (aref (cl-mpm:sim-mps *sim*) 0) (lambda (mesh mp &rest args) (setf (fill-pointer (cl-mpm/particle::mp-cached-nodes mp)) 0)))
  ;;              (cl-mpm::iterate-over-neighbours-shape-gimp-simd (cl-mpm:sim-mesh *sim*) (aref (cl-mpm:sim-mps *sim*) 0) (lambda (&rest args)))
  ;;              ))
  ;; (time
  ;;  (dotimes (i 100)
  ;;    (format t "~D~%" i)
  ;;    (cl-mpm::update-sim *sim*)))
  (format t "MPS ~D~%" (length (cl-mpm:sim-mps *sim*)))
  ;; (sb-profile:report)
  )

(defun profile ()
  (setup :refine 16)
  (time-form 100
             (progn
               (format t "~D~%" i)
               (cl-mpm::update-sim *sim*)))
  ;; (time-form
  ;;  100
  ;;  (progn
  ;;    (cl-mpm::check-mps *sim*)))
  ;; (time
  ;;  (dotimes (i 100)
  ;;    (format t "~D~%" i)
  ;;    (cl-mpm::update-sim *sim*)))
  (format t "MPS ~D~%" (length (cl-mpm:sim-mps *sim*)))
  ;; (sb-profile:report)
  )

;; (lparallel:end-kernel)
;; (sb-ext::exit)
;; (uiop:quit)



(defun find-nans ()
  (cl-mpm::iterate-over-nodes
   (cl-mpm:sim-mesh *sim*)
   (lambda (mp)
     (loop for v across (magicl::storage (cl-mpm/mesh::node-velocity mp))
           do
              ;; (when (> v 0d0)
              ;;   (format t "~E~%" v))
              (when (or (sb-ext::float-nan-p v)
                        ;; (> v 1d10)
                        ;; (< v -1d10)
                        )
                (format t "Found nan vel~%")
                (pprint mp)
                )
           ))))



;;Old
;; Evaluation took:
;; 7.969 seconds of real time
;; 71.424494 seconds of total run time (53.687638 user, 17.736856 system)
;; [ Real times consist of 2.410 seconds GC time, and 5.559 seconds non-GC time. ]
;; [ Run times consist of 2.467 seconds GC time, and 68.958 seconds non-GC time. ]
;; 896.27% CPU
;; 33,467,579,605 processor cycles
;; 22,937,955,888 bytes consed

;; Total time: 7.969999 
;; Time per iteration: 0.07969999
;; Throughput: 12.547053
;; Time per MP: 4.8645017e-6
;; MPS 16384

;;New


(defun stop ()
  (setf *run-sim* nil)
  (setf cl-mpm/dynamic-relaxation::*run-convergance* nil))
(defun test-conv ()
  (defparameter *data-steps* (list))
  (defparameter *data-oobf* (list))
  (defparameter *data-energy* (list))
  (defparameter *data-name* (list))

  (vgplot:figure)
  (let ((path (merge-pathnames "./analysis_scripts/vel_algo/data/")))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* path)) do (uiop:delete-file-if-exists f))
    (dolist (algo (list :FLIP))
      (dolist (dt-scale (list 0.5d0))
        (dolist (refine (list 1 2 3))
          (dolist (fbar (list nil))
            (dolist (sim-type (list
                               'cl-mpm::mpm-sim-usf
                               'cl-mpm::mpm-sim-usl
                               ))
              (let* ((steps (list))
                     (energy (list))
                     (oobf (list))
                     ;; (dt-scale 0.5d0)
                     (substeps (round 50 dt-scale))
                     )
                (setup :refine refine :mps 3
                       :sim-type sim-type
                       )
                ;; (change-class *sim* sim-type)
                (setf (cl-mpm:sim-damping-factor *sim*)
                      (* 1d-1 (cl-mpm/setup::estimate-critical-damping *sim*)))
                (setf
                 (cl-mpm::sim-velocity-algorithm *sim*)
                 algo
                 (cl-mpm::sim-enable-fbar *sim*) fbar)
                (let ((name (format nil "~A-~A-~E-~A-~A"
                                    sim-type
                                    (cl-mpm::sim-velocity-algorithm *sim*)
                                    dt-scale
                                    refine
                                    fbar
                                    )))
                  (vgplot:title name)
                  (handler-case
                      (cl-mpm/dynamic-relaxation:converge-quasi-static
                       *sim*
                       :dt-scale dt-scale
                       :energy-crit 1d-9
                       :oobf-crit 1d-9
                       :substeps substeps
                       :conv-steps 100
                       :dt-scale dt-scale
                       :pic-update nil
                       :post-iter-step
                       (lambda (i e o)
                         ;; (plot *sim*)
                         (push i steps)
                         (push e energy)
                         (push o oobf)
                         (apply #'vgplot:semilogy
                                (reduce #'append
                                        (append
                                         (mapcar #'list
                                                 (append *data-steps* (list steps))
                                                 (append *data-oobf*  (list oobf))
                                                 (append *data-name*  (list name))
                                                 )))))
                       )
                    (cl-mpm/dynamic-relaxation::non-convergence-error (c)
                      (format t "Sim failed to converge~%"))
                    )
                  (with-open-file (stream (merge-pathnames path (format nil "~A.csv" name)) :direction :output :if-exists :supersede)
                    (format stream "step,energy,oobf~%")
                    (loop for s in steps
                          for e in energy
                          for o in oobf
                          do (format stream "~D,~F,~F~%" s e o)))
                  (cl-mpm/output:save-vtk (merge-pathnames path (format nil "sim_~A.vtk" name)) *sim*)
                  (cl-mpm/output:save-vtk-nodes (merge-pathnames path (format nil "sim_nodes_~A.vtk" name)) *sim*)

                  (push steps *data-steps*)
                  (push energy *data-energy*)
                  (push oobf *data-oobf*)
                  (push name *data-name*)
                  )))))))))

(defun plot-all ()
  (apply #'vgplot:semilogy
         (reduce #'append
                 (mapcar #'list
                         *data-steps*
                         *data-oobf*
                         *data-name*
                         ))))

(defun plot-some ()
  (let* ((indexes (list 0 1 2
                        9 10 11
                        ))
         (steps (loop for i in indexes collect (nth i *data-steps*)))
         (energy (loop for i in indexes collect (nth i *data-energy*)))
         (oobf (loop for i in indexes collect (nth i *data-oobf*)))
         (name (loop for i in indexes collect (nth i *data-name*)))
        )
    (apply #'vgplot:semilogy
           (reduce #'append
                   (mapcar #'list
                           steps
                           oobf
                           name
                           )))))

;; (let ((mesh (cl-mpm:sim-mesh *sim*))
;;       (mp0 (aref (cl-mpm:sim-mps *sim*) 0)))
;;   (pprint (cl-mpm/particle:mp-position mp0))
;;   (cl-mpm::iterate-over-neighbours
;;    mesh mp0
;;    (lambda (mesh mp node svp grads fsvp fgrads)
;;      (with-accessors ((node-active cl-mpm/mesh:node-active)
;;                       (node-j-inc cl-mpm/mesh::node-jacobian-inc)
;;                       (node-volume cl-mpm/mesh::node-volume))
;;          node
;;        (pprint svp)))))

(defun test-ms ()
  (setf *run-sim* t)
  (loop for  ms in (list 1d0 1d1 1d2 1d3)
        while *run-sim*
        do
           (progn
             (setup :refine 2)
             (run :output-dir (format nil "./output-~F/" ms) :ms ms))))

(defun test-usl-usf ()
  (setf *run-sim* t)
  (loop for sim-type in (list
                         'cl-mpm::mpm-sim-usf
                         'cl-mpm::mpm-sim-usl)
        while *run-sim*
        do
           (progn
             (setup :refine 8 :mps 2 :sim-type sim-type)
             (setf (cl-mpm:sim-enable-fbar *sim*) t)
             (run :output-dir (format nil "./output-~A/" sim-type)))))


(defun test-max-conv ()
  (defparameter *data-steps* (list))
  (defparameter *data-oobf* (list))
  (defparameter *data-energy* (list))
  (defparameter *data-name* (list))

  (vgplot:figure)
  (let ((path (merge-pathnames "./analysis_scripts/vel_algo/data/")))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* path)) do (uiop:delete-file-if-exists f))
    (dolist (algo (list :FLIP))
      (dolist (dt-scale (list 1d0 0.5d0 0.1d0))
        (dolist (refine (list 1 2 3))
          (dolist (fbar (list nil))
            (dolist (sim-type (list
                               'cl-mpm::mpm-sim-usf
                               'cl-mpm::mpm-sim-usl
                               ))
              (let* ((steps (list))
                     (energy (list))
                     (oobf (list))
                     ;; (dt-scale 0.5d0)
                     (substeps (round (* refine 50) dt-scale))
                     )
                (setup :refine refine :mps 3
                       :sim-type sim-type
                       )
                ;; (change-class *sim* sim-type)
                (setf (cl-mpm:sim-damping-factor *sim*)
                      (* 1d-1 (cl-mpm/setup::estimate-critical-damping *sim*)))
                (setf
                 (cl-mpm::sim-velocity-algorithm *sim*)
                 algo
                 (cl-mpm::sim-enable-fbar *sim*) fbar)
                (let ((name (format nil "~A-~A-~E-~A-~A"
                                    sim-type
                                    (cl-mpm::sim-velocity-algorithm *sim*)
                                    dt-scale
                                    refine
                                    fbar
                                    )))
                  (vgplot:title name)
                  (let ((work 0d0)
                        (e 0d0)
                        (o 0d0)
                        (sim *sim*)
                        (max-steps 200))
                      (loop for i from 0 to max-steps
                            while *run-sim*
                            do
                               (progn
                                 (dotimes (j substeps)
                                   (cl-mpm:update-sim *sim*))
                                 (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm::calculate-min-dt sim)))
                                 (setf e (cl-mpm::sim-stats-energy sim))
                                 (if (= work 0d0)
                                     (setf e 0d0)
                                     (setf e (abs (/ e work))))
                                 (setf o (cl-mpm::sim-stats-oobf sim))
                                 (push i steps)
                                 (push e energy)
                                 (push o oobf)
                                 (apply #'vgplot:semilogy
                                        (reduce #'append
                                                (append
                                                 (mapcar #'list
                                                         (append *data-steps* (list steps))
                                                         (append *data-oobf*  (list oobf))
                                                         (append *data-name*  (list name))
                                                         )))))
                               (swank.live:update-swank)
                       ))
                  (with-open-file (stream (merge-pathnames path (format nil "~A.csv" name)) :direction :output :if-exists :supersede)
                    (format stream "step,energy,oobf~%")
                    (loop for s in steps
                          for e in energy
                          for o in oobf
                          do (format stream "~D,~F,~F~%" s e o)))
                  (cl-mpm/output:save-vtk (merge-pathnames path (format nil "sim_~A.vtk" name)) *sim*)
                  (cl-mpm/output:save-vtk-nodes (merge-pathnames path (format nil "sim_nodes_~A.vtk" name)) *sim*)

                  (push steps *data-steps*)
                  (push energy *data-energy*)
                  (push oobf *data-oobf*)
                  (push name *data-name*)
                  )))))))))


(defun get-displacement-norm (sim)
  (cl-mpm::reduce-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (* (cl-mpm/particle::mp-mass mp)
        (cl-mpm/fastmaths::mag
         (cl-mpm/particle::mp-displacement-increment mp))))
   #'+))

(defun test-arc ()
  (setup)
  (run-static-arc)
  )

(defun run-static-arc (&key
                         (output-dir "./output/")
                     )
  (uiop:ensure-all-directories-exist (list output-dir))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (setf (cl-mpm:sim-damping-factor *sim*)
        (* 1d-1
           (cl-mpm/setup:estimate-critical-damping *sim*)))
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (let* ((load-steps 50)
           (load-0 0.1d0)
           (load (cl-mpm/particle::mp-gravity (aref (cl-mpm:sim-mps *sim*) 0)))
           (load-inc (/ (- 1d0 load-0) load-steps))
           (current-load load-0)
           (arc-length 0d0)
           (target-arc-length nil)
           )
      (defparameter *loading-path* (list))
      (defparameter *data-disp* (list))
      (loop for step from 0 below load-steps
            while *run-sim*
            do
               (let ((load-factor (+ current-load load-inc))
                     (dl 0d0))
                 (pprint load-factor)
                 (cl-mpm:iterate-over-mps
                  mps
                  (lambda (mp)
                    ;; (new-loadstep mp)
                    (setf (cl-mpm/particle:mp-gravity mp) (* load load-factor))))
                 (cl-mpm/dynamic-relaxation:converge-quasi-static
                  *sim*
                  :oobf-crit 1d-2
                  :energy-crit 1d-2
                  :substeps 100
                  :conv-steps 1000
                  :convergance-criteria
                  (lambda (sim energy oobf)
                    (setf arc-length (get-displacement-norm sim))
                    (and
                     (if target-arc-length
                         (< (/ (abs (- arc-length target-arc-length)) target-arc-length) 1d-2)
                         t)
                     (< oobf 1d-2)
                     (< energy 1d-2)
                     ))
                  :post-iter-step
                  (lambda (i energy oobf)
                    (setf arc-length (get-displacement-norm *sim*))
                    (plot)
                    (vgplot:title (format nil "Step ~D - Time ~F - KE ~E - OOBF ~E - Target ~E - Arc ~A - LF ~E"  step (cl-mpm::sim-time *sim*) energy oobf target-arc-length
                                          (if target-arc-length
                                              (/ arc-length target-arc-length)
                                              0d0) load-factor))
                    (format t "Substep ~D~%" i)
                    (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d_~5,'0d.vtk" step i)) *sim*)
                    (when target-arc-length
                      (let ((delta (/ (- target-arc-length arc-length ) target-arc-length)))
                        (format t "Delta ~E~%" delta)
                        (incf load-factor
                              (*
                               1d0
                               ;; (cl-mpm:sim-dt *sim*)
                               ;; (max (- 1d0 oobf) 0d0)
                               oobf
                               delta
                               ))
                        (cl-mpm:iterate-over-mps
                         mps
                         (lambda (mp)
                           (setf (cl-mpm/particle:mp-gravity mp) (* load-factor load))))))
                    ))
                 (unless target-arc-length
                   (setf target-arc-length (get-displacement-norm *sim*))
                   (format t "Arc-length set at ~E~%" target-arc-length)
                   )
                 (setf current-load load-factor)
                 (push load-factor *loading-path*)
                 (push (get-displacement-norm *sim*) *data-disp*)
                 )
               (cl-mpm::new-loadstep *sim*)
               ;; (plot)

               ;; (vgplot:title (format nil "Time:~F - KE ~E - OOBF ~E"  (cl-mpm::sim-time *sim*) energy oobf))
               (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step))
                                  :terminal "png size 1920,1080"
                                  )
               (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" step)) *sim*)
               (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) *sim*)
               (sleep 0.1d0)
               (swank.live:update-swank)

            )
      (vgplot:figure)
      (vgplot:plot data-disp loading-path "Load")
      )))


(defun test-load-control ()
  (setup :mps 3 :refine 0.5 :multigrid-refine 0)
  (cl-mpm/setup::set-mass-filter *sim* *density* :proportion 1d-15)
  ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul)
  (setf (cl-mpm::sim-gravity *sim*) -10d0)
  ;; (setf (cl-mpm::sim-gravity *sim*) -1000d0)
  (cl-mpm/dynamic-relaxation::run-load-control
   *sim*
   :output-dir (format nil "./output/")
   :plotter #'plot
   :load-steps 10
   :damping 1d0;(sqrt 2)
   :substeps 10
   :criteria 1d-4
   :save-vtk-dr t
   :save-vtk-loadstep t
   :dt-scale 1d0))


(defun test-3d ()
  (setup :mps 3 :refine 1)
  (cl-mpm/setup::set-mass-filter *sim* *density* :proportion 1d-9)

  (setf (cl-mpm::sim-ghost-factor *sim*) (* 1d6 1d-4))
  (setf (cl-mpm::sim-gravity *sim*) -1000d0)
  (vgplot:close-all-plots)
  (let ((step (list))
        (res (list))
        (total-step 0))
    (cl-mpm/dynamic-relaxation::run-load-control
     *sim*
     :output-dir (format nil "./output/")
     :load-steps 10
     :damping (sqrt 2)
     :substeps 10
     :criteria 1d-9
     :save-vtk-dr t
     :save-vtk-loadstep nil

     :plotter (lambda (sim)
                           (vgplot:semilogy (reverse step) (reverse res))
                           ); #'plot-sigma-yy
                :post-iter-step (lambda (i o e)
                                  (push total-step step)
                                  (incf total-step)
                                  (push e res)
                                  )

     :dt-scale 0.5d0)))

(defun test-pressure ()
  (setup :mps 2 :refine 2)
  (cl-mpm/setup::set-mass-filter *sim* *density* :proportion 1d-15)
  (cl-mpm:add-bcs-force-list
   *sim*
   *bc-pressure*)
  (cl-mpm/setup::remove-sdf *sim*
                            (lambda (pos)
                              (-
                               (funcall
                                (cl-mpm/setup::circle-sdf
                                 (list 0d0 0d0 0d0)
                                 8d0
                                 )
                                pos
                                )))
                            
                            :refine 1)
  ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul)
  (vgplot:close-all-plots)
  (setf (cl-mpm::sim-gravity *sim*) 0d0)
  (let ((pressure -1d7)
        (step (list))
        (res (list))
        (total-step 0)
        )
    (cl-mpm/dynamic-relaxation::run-load-control
     *sim*
     :output-dir (format nil "./output/")
     :loading-function (lambda (percent)
                         (setf
                          (cl-mpm/buoyancy::bc-pressure-pressures *bc-pressure*)
                          (list
                           (* pressure percent)
                           (* pressure percent))))
     ;; :plotter #'plot
     :load-steps 10
     :damping (sqrt 2)
     :substeps 10
     :criteria 1d-9
     :save-vtk-dr t
     :save-vtk-loadstep nil

     :plotter (lambda (sim)
                           (vgplot:semilogy (reverse step) (reverse res))
                           ); #'plot-sigma-yy
                :post-iter-step (lambda (i o e)
                                  (push total-step step)
                                  (incf total-step)
                                  (push e res)
                                  )

     :dt-scale 0.5d0)))
(defun test-static ()
  (setup :mps 2 :refine 0.5)
  (let ((output-dir "./output/")
        (dt-scale 0.5d0))
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 0.1d0
             (cl-mpm/setup:estimate-critical-damping *sim*)))
    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (cl-mpm/setup::set-mass-filter *sim* *density* :proportion 1d-9)
    (loop for i from 0 to 1000
          while *run-sim*
          do
             (progn
               (dotimes (j 10)
                 (cl-mpm:update-sim *sim*))
               (cl-mpm/output:save-vtk (merge-pathnames (format nil "sim_~5,'0d.vtk" i) output-dir) *sim*)
               (cl-mpm/output:save-vtk-nodes (merge-pathnames (format nil "sim_nodes_~5,'0d.vtk" i) output-dir) *sim*)
               (cl-mpm/output:save-vtk-cells (merge-pathnames (format nil "sim_cells_~5,'0d.vtk" i) output-dir) *sim*)
               (plot *sim*)
               (vgplot:title (format nil "~D" i))
               (swank.live:update-swank)))))


(defun run-static (&key (output-dir "./output/")
                     (load-steps 10)
                     (dt-scale 0.5d0))
  (uiop:ensure-all-directories-exist (list output-dir))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (setf (cl-mpm:sim-damping-factor *sim*)
        (* 1d-1
           (cl-mpm/setup:estimate-critical-damping *sim*)))
  (with-open-file (stream (merge-pathnames output-dir "conv.csv") :direction :output :if-exists :supersede)
    (format stream "iter,step,damage,plastic,oobf,energy~%"))
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (let* ((load-0 0d0)
           ;; (load -5d0)
           (load (cl-mpm/particle::mp-gravity (aref (cl-mpm:sim-mps *sim*) 0)))
           (load-inc (/ (- load load-0) load-steps))
           (current-load load-0)
           (total-iter 0)
           )
      (loop for step from 0 below load-steps
            while *run-sim*
            do
               (incf current-load load-inc)
               (cl-mpm:iterate-over-mps
                mps
                (lambda (mp)
                  (setf (cl-mpm/particle:mp-gravity mp) current-load)))
               (let ((crit (cl-mpm/setup:estimate-critical-damping *sim*)))
                 ;; (setf (cl-mpm:sim-damping-factor *sim*)
                 ;;       (* 0.1d0 crit))
                 (defparameter *ke-last* 0d0)
                 ;; (cl-mpm:update-sim *sim*)
                 (let ((conv-steps 0)
                       (substeps 50)
                       (i 0))
                   ;; (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d_~5,'0d.vtk" step i)) *sim*)
                   ;; (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d_~5,'0d.vtk" step i)) *sim*)
                   ;; (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d_~5,'0d.vtk" step i)) *sim*)
                   (time
                    (cl-mpm/dynamic-relaxation:converge-quasi-static
                     *sim*
                     :oobf-crit 1d-3
                     :energy-crit 1d-3
                     :kinetic-damping t
                     :damping-factor 1d-2
                     :dt-scale dt-scale
                     :substeps substeps
                     :conv-steps 1000
                     :post-iter-step
                     (lambda (i energy oobf)
                       (setf conv-steps (* substeps i))
                       (plot *sim*)
                       (vgplot:title (format nil "Step ~D - substep ~D - KE ~E - OOBF ~E"  step i energy oobf))
                       (format t "Substep ~D~%" i)
                       (let ((i (+ 0 i)))
                         (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d_~5,'0d.vtk" step i)) *sim*)
                         (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d_~5,'0d.vtk" step i)) *sim*)
                         (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d_~5,'0d.vtk" step i)) *sim*))
                       (with-open-file (stream (merge-pathnames output-dir "conv.csv") :direction :output :if-exists :append)
                         (format stream "~D,~D,~f,~f,~f,~f~%" total-iter step (get-plastic) (get-damage)
                                 oobf energy))
                       (incf total-iter)
                       )))
                   (plot *sim*)
                   (vgplot:title (format nil "Step ~D - ~D" step conv-steps))
                   ))
               ;; (cl-mpm/penalty::bc-increment-center *bc-squish* (cl-mpm/utils:vector-from-list (list 0d0 load-inc 0d0)))
               (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) *sim*)
               (cl-mpm/output::save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d.vtk" step)) *sim*)
               (cl-mpm::finalise-loadstep *sim*)
               (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step))
                                  :terminal "png size 1920,1080"
                                  )
               (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" step)) *sim*)
               (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim_p_~5,'0d.vtk" step)) *sim* )
               (sleep 0.1d0)
               (swank.live:update-swank)

            ))))



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

(defun test-quasi-time ()
  (setup :mps 2 :refine 1)
  (let ((dt 0.01d0)
        (total-time 100d0))
    (cl-mpm/dynamic-relaxation::run-quasi-time
     *sim*
     :output-dir "./output/"
     :dt-scale 1d0
     :dt dt
     :steps (round total-time dt)
     ;; :time total-time
     )))


(defun save-test-vtks (&key (output-dir "./output/"))
  (cl-mpm/output:save-vtk (merge-pathnames "test.vtk" output-dir) *sim*)
  (cl-mpm/output:save-vtk-nodes (merge-pathnames "test_nodes_0.vtk" output-dir) *sim*)
  (cl-mpm/output:save-vtk-nodes (merge-pathnames "test_nodes_1.vtk" output-dir) *sim*)
  (cl-mpm/output:save-vtk-cells (merge-pathnames "test_cells.vtk" output-dir) *sim*))





(defun test-multi-step ()
  (setup :mps 2 :refine 0.5)
  (run-multi-stage
   :output-dir "./output/"
   :dt 0.5d0
   ;; :dt-scale (/ 0.8d0 (sqrt 1d0))
   ;; :load-steps 20
   )
  )


(require 'sb-sprof)

(defun spy (mat)
  (let ((pos-x (list))
        (pos-y (list)))
    (vgplot:close-all-plots)
    (vgplot:figure)
    (destructuring-bind (lx ly) (magicl:shape mat)
      (loop for x from 0 below lx
            do (loop for y from 0 below ly
                     do (when (> (abs (cl-mpm/utils::mtref mat x y)) 0d0)
                          (push x pos-x)
                          (push y pos-y)
                          ;; (push (- ly (+ 1 y)) pos-y)
                          )))
      (let ((ad 0.25d0))
        (vgplot:format-plot t "set xrange [~f:~f]" (- ad) (- (+ lx ad) 1d0))
        (vgplot:format-plot t "set yrange [~f:~f]" (- ad) (- (+ ly ad) 1d0)))
      (vgplot:plot pos-x pos-y ";;with points"))))


(defun test-e ()
  (setup :mps 2 :refine 0.125)
  (cl-mpm:update-sim *sim*)
  (let* ((sim *sim*)
         (E (cl-mpm/aggregate::assemble-global-e sim))
         (mii (cl-mpm/aggregate::assemble-global-mass sim))
         (f (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-force))
         (ma
           (magicl:@ (magicl:transpose E)
                     mii
                     E))
         )
    (let ((elem (aref (cl-mpm/aggregate::sim-agg-elems *sim*) 0)))
      (format t "~A~%" (cl-mpm/aggregate::agg-shape-functions elem)))
    ;; (spy ma)
    ;; (pprint E)
    ;; (pprint (magicl:inv ma))
    ;; (spy ma)
    ;; (spy E)
    ;; (spy (cl-mpm/aggregate::assemble-global-mass *sim*))
    ;; (spy E)
    ;; (cl-mpm/aggregate::project-global-vec sim acc cl-mpm/mesh::node-acceleration)
    (save-test-vtks)
    ))


(defun test-ssr ()
  (let ((phi 50d0)
        (phi-r 0d0)
        (c 30d3)
        (c-r 0d0))
    (setup :refine 1 :mps 2)
    (cl-mpm/dynamic-relaxation::run-load-control
     *sim*
     :output-dir (format nil "./output/")
     :plotter #'plot
     :load-steps 50
     :damping 1d-3
     :substeps 50
     ;; :conv-steps 1000
     :criteria 1d-2
     :adaptive-damping nil
     :kinetic-damping t
     :save-vtk-dr t
     :save-vtk-loadstep t
     :dt-scale 0.5d0
     :loading-function
     (lambda (percent)
       (cl-mpm:iterate-over-mps
        (cl-mpm:sim-mps *sim*)
        (lambda (mp)
          (setf (cl-mpm/particle::mp-c mp)
                (+ c-r (* (- 1d0 percent) (- c c-r))))))))))



(defun test-static ()
  (dolist (gf-factor (list 1d-3))
    (let ((dt-scale 0.50d0)
          (output-dir (format nil "./output-~E/" gf-factor)))
      (ensure-directories-exist output-dir)
      (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
      (setup :mps 2 :refine 0.5)
      (defparameter *data-steps* (list))
      (defparameter *data-oobf* (list))
      (defparameter *data-energy* (list))
      (setf (cl-mpm::sim-ghost-factor *sim*) (* 1d6 gf-factor))
      (cl-mpm/dynamic-relaxation::save-conv-preamble output-dir)
      (handler-case
          (cl-mpm/dynamic-relaxation:converge-quasi-static
           *sim*
           :dt-scale dt-scale
           :damping-factor 0.1d0
           :kinetic-damping nil
           :conv-steps 100
           :substeps 10
           :oobf-crit 1d-5
           :energy-crit 1d-5
           :post-iter-step
           (lambda (i energy oobf)
             (push i *data-steps*)
             (push oobf *data-energy*)
             (plot-load-disp)
             (cl-mpm/dynamic-relaxation::save-conv-step *sim* output-dir i 0 0d0 oobf energy)
             (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" i)) *sim*)
             (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" i)) *sim*)))
          (cl-mpm/dynamic-relaxation::non-convergence-error (c)
            (progn
              (format t "Sim failed to converge~%")
              t))
          (t () t)
        )
      ))
  )

(defun test-implicit-dynamic ()
  (setup :mps 3 :refine 2)
  (let* ((dt 0.1d0)
         (total-time 100d0)
         (output-dir "./output/")
         (target-time 0.1d0)
         (substeps (floor target-time dt)))
    (uiop:ensure-all-directories-exist (list output-dir))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic)
    (cl-mpm/dynamic-relaxation::save-timestep-preamble output-dir)
    (cl-mpm/dynamic-relaxation::save-conv-preamble output-dir)
    (cl-mpm/dynamic-relaxation::elastic-static-solution *sim* :substeps 100)
    (setf (cl-mpm:sim-dt *sim*) dt)
    (loop for steps from 0 to 1000
     while *run-sim*
     do
     (let ((energy 0d0)
           (work 0d0)
           (oobf 0d0))
       (cl-mpm/dynamic-relaxation::save-timestep *sim* output-dir steps :DYNAMIC)
       (format t "Step ~d ~%" steps)
       (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" *sim-step*)) *sim*)
       (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
       (time
        (dotimes (i substeps)
          (cl-mpm::update-sim *sim*)))
       (incf *sim-step*)
       (plot *sim*)
       (vgplot:title (format nil "Time:~F"  (cl-mpm::sim-time *sim*)))
       (swank.live:update-swank)
       ))))




;; (let* ((h 1d0)
;;        (g (list (/ 1d0 h) (/ 1d0 h) (/ 1d0 h)))
;;        (df (cl-mpm/utils:matrix-from-list (list 0.1d0 1d0 0d0
;;                                                 0d0 0.1d0 0d0
;;                                                 0d0 0d0 1d0))))
;;   (pprint g)
;;   (pprint (cl-mpm::gradient-push-forwards g df))
;;   ;; (multiple-value-bind (l v) (magicl:eig df)
;;   ;;   (pprint l)
;;   ;;   (pprint (reduce #'max (mapcar (lambda (x) (/ 1d0 (abs x))) l)))
;;   ;;   )
;;   )

(let* ((E 1d0)
       (nu 0.1d0)
       (phi 0.1d0)
       (psi 0.05d0)
       (c 1d0)
       (strain (cl-mpm/utils::voigt-from-list (list 1d0 1d0 1d0 0d0 0d0 0d0)))
       (de (cl-mpm/constitutive:linear-elastic-matrix E nu))
       (stress (magicl:@ de strain))
       )
  ;; (pprint strain)
  ;; (pprint stress)
  (multiple-value-bind (sig eps f psinc pmod)
      ;; (cl-mpm/constitutive::mc-plastic stress de strain E nu phi psi c)
    (cl-mpm/ext::constitutive-mohr-coulomb stress de strain E nu phi psi c)
    ;; (pprint sig)
    (pprint eps)
    ;; (pprint de)
    (pprint psinc)
    (pprint pmod)
    ))
