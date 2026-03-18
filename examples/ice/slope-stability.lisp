(defpackage :cl-mpm/examples/ice/slope-stability
  (:use :cl
   :cl-mpm/example))
(in-package :cl-mpm/examples/ice/slope-stability)


(declaim (notinline plot-domain))
(defun plot-domain ()
  (when *sim*
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :trial t
     :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp)))))

(defun get-load ()
  (cl-mpm/penalty::resolve-load *penalty*))

(declaim (notinline setup))
(defun setup (&key
                (refine 1)
                (mps 3)
                (enable-fbar t)
                (multigrid-refines 0)
                (angle 20d0)
                (angle-r 0d0)
                (rt 1d0)
                (rc 0d0)
                (l-scale 1d0)
                (gf 0.1d3)
                )

  (let* ((refine (* refine (expt 2 (- multigrid-refines))))
         (height 10d0)
         (width 20d0)
         (h (/ 1d0 refine))
         (density 2.04d3)
         (E 10d6)
         (nu 0.4d0)
         (domain-width 30d0)
         (domain-height (+ domain-width (* 2 h)))
         (domain-size (list domain-width domain-height))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list width height)))
    (setf
     *sim*
     (cl-mpm/setup::make-simple-sim
      h
      element-count
      :sim-type
       ;; 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
       'cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul
       ;; 'cl-mpm/dynamic-relaxation::mpm-sim-dr-multigrid
       :args-list
       (list
        :enable-aggregate nil
        :ghost-factor (* E 1d-4)
        ;; :enable-aggregate t
        ;; :ghost-factor nil
        :enable-split nil
        :enable-fbar enable-fbar
        ;; :refinement multigrid-refines
        )))
    (setf h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
    (let* ((c 40d3)
           (init-stress (cl-mpm/damage::mohr-coloumb-coheasion-to-tensile c angle))
           (rs (cl-mpm/damage::est-shear-from-angle angle angle-r rc))
           (L (* h l-scale))
           (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf L init-stress E))
           ;; (rt 1d0)
           ;; (rc 1d0)
           ;; (rs 1d0)
           )
      (when (< ductility 1d0)
        (error "Ductility less than 1"))
      (format t "Init stress ~E~%" init-stress)
      (format t "Ductility ~E~%" ductility)
      (format t "rc ~E~%rt ~E~%rs ~E~%" rc rt rs)
      (format t "Angle ~E residual ~E~%" angle angle-r)
      (cl-mpm:add-mps
       *sim*
       (cl-mpm/setup:make-block-mps
        (list 0 0)
        block-size
        (mapcar (lambda (e) (* (/ e h) mps)) block-size)
        density
        'cl-mpm/particle::particle-damage-frictional
        :E E
        :nu nu
        :local-length L
        :ductility ductility
        :friction-angle (cl-mpm/utils:deg-to-rad angle)
        :initiation-stress init-stress
        :friction-model :DP
        :kt-res-ratio rt
        :kc-res-ratio rc
        :g-res-ratio rs
        ;; :delay-time 100d0
        ;; :delay-exponent 2d0

        ;; 'cl-mpm/particle::particle-elastic
        ;; :E E
        ;; :nu nu
        ;; 'cl-mpm/particle::particle-vm
        ;; :E E
        ;; :nu nu
        ;; :rho (* 1d4 (sqrt 3/2))
        ;; :rho 1d6
        ;; :rho 1d6
        ;; 'cl-mpm/particle::particle-mc
        ;; :E E
        ;; :nu 0.25d0
        ;; :psi (* 30d0 (/ pi 180))
        ;; :phi (* 30d0 (/ pi 180))
        ;; :c 1d3
        )))
    (cl-mpm/setup:remove-sdf
     *sim*
     (lambda (p)
       (cl-mpm/setup::plane-point-point-sdf
        p
        (cl-mpm/utils:vector-from-list (list 10d0 10d0 0d0))
        (cl-mpm/utils:vector-from-list (list 20d0 0d0 0d0))))
     :refine 1)

    (setf (cl-mpm::sim-gravity *sim*) -9.8d0)
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)
    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 nil nil)
     :right '(0 nil nil)
     :bottom '(0 0 nil))


    (let* ((friction 0d0)
           (epsilon-scale 1d0)
           (epsilon (* (cl-mpm/particle::calculate-p-wave-modulus E nu) epsilon-scale))
           (width 2d0)
           (offset 8d0)
           )
      (format t "Penalty parameter ~E~%" epsilon)
      (defparameter *penalty-down*
        (cl-mpm/penalty::make-bc-penalty-distance-point
         *sim*
         (cl-mpm/utils:vector-from-list '(0d0 -1d0 0d0))
         (cl-mpm/utils:vector-from-list (list
                                         offset
                                         height
                                         0d0))
         width
         epsilon
         friction
         0d0))
      (defparameter *penalty-right*
        (cl-mpm/penalty::make-bc-penalty-distance-point
         *sim*
         (cl-mpm/utils:vector-from-list '(1d0 0d0 0d0))
         (cl-mpm/utils:vector-from-list (list
                                         (+ offset (* width))
                                         (+ height (/ width 2))
                                         0d0))
         (/ width 2)
         epsilon
         friction
         0d0))
      (defparameter *penalty-left*
        (cl-mpm/penalty::make-bc-penalty-distance-point
         *sim*
         (cl-mpm/utils:vector-from-list '(-1d0 0d0 0d0))
         (cl-mpm/utils:vector-from-list (list
                                         (- offset (* width))
                                         (+ height (/ width 2))
                                         0d0))
         (/ width 2)
         epsilon
         friction
         0d0))
      (defparameter *penalty*
        (cl-mpm/penalty::make-bc-penalty-structure
         *sim*
         epsilon
         friction
         0d0
         (list
          *penalty-left*
          *penalty-down*
          *penalty-right*
          ))))
    (cl-mpm:add-bcs-force-list
     *sim*
     *penalty*)
    )
  (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
  )

(defun plot-load-disp ()
  (vgplot:plot
   (mapcar (lambda (x) (* x -1d3)) *data-disp*)
   (mapcar (lambda (x) x) *data-load*)
   "Solution"
   )
  (vgplot:xlabel "Displacement (mm)")
  )

(defun save-csv (output-dir filename data-disp data-load)
  (with-open-file (stream (merge-pathnames filename output-dir) :direction :output :if-exists :supersede)
    (format stream "disp,load~%")
    (loop for disp in data-disp
          for load in data-load
          do (format stream "~E,~E~%" (float disp 0e0) (float load 0e0)))))
(declaim (notinline run))
(defun run (&key (output-dir (format nil "./output/"))
              (csv-dir (format nil "./output/"))
              (csv-filename (format nil "load-disp.csv")))
  (ensure-directories-exist output-dir)
  (ensure-directories-exist csv-dir)
  (cl-mpm/dynamic-relaxation::elastic-static-solution
   *sim*)
  (let* ((lstps 50)
         (total-disp -1d0)
         (current-disp 0d0)
         (disp-0 (cl-mpm::reduce-over-mps (cl-mpm:sim-mps *sim*)
                                          (lambda (mp)
                                            (cl-mpm/utils:get-vector (cl-mpm/particle::mp-displacement mp) :y))
                                          #'min))
         (step 0))
    (defparameter *data-disp* (list 0d0))
    (defparameter *data-load* (list 0d0))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))

    (vgplot:close-all-plots)
    (time
     (cl-mpm/dynamic-relaxation::run-adaptive-load-control
      *sim*
      :output-dir output-dir
      ;:plotter (lambda (sim) (plot-load-disp))
      :plotter (lambda (sim) (plot-domain))
      :loading-function
      (lambda (i)
        (setf current-disp (+ (* i total-disp) disp-0))
        (cl-mpm/penalty::bc-set-displacement
         *penalty*
         (cl-mpm/utils:vector-from-list (list 0d0 current-disp 0d0))))
      :post-conv-step
      (lambda (sim)
        (push current-disp *data-disp*)
        (let ((load (get-load)))
          (format t "Load ~E~%" load)
          (push load *data-load*))
        (plot-load-disp)
        (save-csv csv-dir csv-filename *data-disp* *data-load*)
        (incf step))
      :load-steps lstps
      :max-adaptive-steps 20
      :enable-plastic t
      :enable-damage t
      :damping 1d0;(sqrt 2d0)
      :max-damage-inc 0.5d0
      :substeps 10
      :criteria 1d-3
      :save-vtk-dr t
      :save-vtk-loadstep t
      :dt-scale 1d0))))

(defun test ()
  (setup :mps 3 :refine 1 :enable-fbar t :multigrid-refines 0 :gf 10000d0)
  ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
  (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
  (run)
  ;; (save-csv "./examples/fbar/rigid-footing/" (format nil "data_fbaradjust_~A.csv" t) *data-disp* *data-load*)
  ;; (dolist (fbar (list t nil))
  ;;   (save-csv "./examples/fbar/rigid-footing/" (format nil "data_fbar_~A.csv" fbar) *data-disp* *data-load*)
  ;;   )
  )



