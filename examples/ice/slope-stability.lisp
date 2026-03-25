(defpackage :cl-mpm/examples/ice/slope-stability
  (:use :cl
   :cl-mpm/example))
(in-package :cl-mpm/examples/ice/slope-stability)



(defclass cl-mpm/particle::particle-fpd-isotropic (cl-mpm/particle::particle-plastic-damage-frictional)
  ())
(defclass cl-mpm/particle::particle-fpd-tcs (cl-mpm/particle::particle-plastic-damage-frictional)
  ())
(defclass cl-mpm/particle::particle-fpd-tcs-strain (cl-mpm/particle::particle-plastic-damage-frictional)
  ())
(defclass cl-mpm/particle::particle-fpd-spectral (cl-mpm/particle::particle-plastic-damage-frictional)
  ())

(defclass cl-mpm/particle::particle-fpd-spectral-strain (cl-mpm/particle::particle-plastic-damage-frictional)
  ())

(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-fpd-tcs) dt)
  (cl-mpm/damage::apply-tcs-degredation mp))

(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-fpd-tcs-strain) dt)
  (cl-mpm/damage::apply-tcs-strain-degredation mp))

(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-fpd-spectral) dt)
  ;; (setf (cl-mpm/particle::mp-damage mp) (cl-mpm/particle::mp-damage-tension mp))
  (cl-mpm/damage::apply-tensile-stress-degredation mp))

(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-fpd-spectral-strain) dt)
  (setf (cl-mpm/particle::mp-damage mp) (cl-mpm/particle::mp-damage-tension mp))
  (cl-mpm/damage::apply-tensile-strain-degredation mp))

(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-fpd-isotropic) dt)
  (setf (cl-mpm/particle::mp-damage mp) (cl-mpm/particle::mp-damage-tension mp))
  (cl-mpm/damage::apply-isotropic-degredation mp))

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
                (angle 15d0)
                (angle-r 0d0)
                (rt 1d0)
                (rc 0d0)
                (l-scale 1d0)
                (gf 0.1d3)
                (oversize-factor (- 1d0 1d-1))
                )

  (let* ((refine (* refine (expt 2 (- multigrid-refines))))
         (height 10d0)
         (width (* height 2))
         (h (/ 1d0 refine))
         (density 20d3)
         (E 10d6)
         (nu 0.4d0)
         (domain-width (+ width height))
         (domain-height (+ height (* 2 h)))
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
        :ghost-factor (* E 1d-2)
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
           (oversize (cl-mpm/damage::compute-oversize-factor oversize-factor ductility))
           ;; (rt 1d0)
           ;; (rc 1d0)
           ;; (rs 1d0)
           )
      (when (< ductility 1d0)
        (error "Ductility less than 1"))
      (format t "Init stress ~E~%" init-stress)
      (format t "Ductility ~E~%" ductility)
      (format t "Oversize ~E~%" oversize)
      (format t "rc ~E~%rt ~E~%rs ~E~%" rc rt rs)
      (format t "Angle ~E residual ~E~%" angle angle-r)
      (cl-mpm:add-mps
       *sim*
       (cl-mpm/setup:make-block-mps
        (list 0 0)
        block-size
        (mapcar (lambda (e) (* (/ e h) mps)) block-size)
        density
        ;; 'cl-mpm/particle::particle-damage-frictional
        ;; :E E
        ;; :nu nu
        ;; :local-length L
        ;; :ductility ductility
        ;; :friction-angle (cl-mpm/utils:deg-to-rad angle)
        ;; :initiation-stress init-stress
        ;; :friction-model :MC
        ;; :kt-res-ratio rt
        ;; :kc-res-ratio rc
        ;; :g-res-ratio rs

        ;; 'cl-mpm/particle::particle-ice-brittle
        ;; :E E
        ;; :nu nu

        ;; :kt-res-ratio rt
        ;; :kc-res-ratio rc
        ;; :g-res-ratio rs

        ;; :initiation-stress init-stress;18d3
        ;; :friction-angle angle
        ;; :psi (cl-mpm/utils:deg-to-rad angle)
        ;; :phi (cl-mpm/utils:deg-to-rad angle)
        ;; :c (* c oversize)

        ;; :ductility ductility
        ;; :local-length L
        ;; :enable-plasticity nil
        ;; :enable-damage t

        'cl-mpm/particle::particle-plastic-damage-frictional
        :E E
        :nu nu
        :local-length L
        :ductility ductility
        :friction-angle (cl-mpm/utils:deg-to-rad angle)
        :residual-friction (cl-mpm/utils:deg-to-rad angle-r)
        :initiation-stress init-stress
        :friction-model :MC
        :oversize oversize-factor
        :kt-res-ratio rt
        :kc-res-ratio rc
        :psi 0d0
        :plastic-damage-evolution nil

        ;; 'cl-mpm/particle::particle-mc
        ;; :E E
        ;; :nu nu
        ;; :psi 0d0;(cl-mpm/utils:deg-to-rad angle)
        ;; :phi (cl-mpm/utils:deg-to-rad angle)
        ;; :c c
        )
       ))
    (cl-mpm/setup:remove-sdf
     *sim*
     (lambda (p)
       (cl-mpm/setup::plane-point-point-sdf
        p
        (cl-mpm/utils:vector-from-list (list height height 0d0))
        (cl-mpm/utils:vector-from-list (list width 0d0 0d0))))
     :refine 1)

    (setf (cl-mpm::sim-gravity *sim*) -1d0)
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)
    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 nil nil)
     :right '(0 nil nil)
     :bottom '(0 0 nil))


    (let* ((friction 1d0)
           (epsilon-scale 1d2)
           (epsilon (* (cl-mpm/particle::calculate-p-wave-modulus E nu) epsilon-scale))
           (penalty-width 2d0)
           (offset (- height penalty-width))
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
         penalty-width
         epsilon
         friction
         0d0))
      (defparameter *penalty-right*
        (cl-mpm/penalty::make-bc-penalty-distance-point
         *sim*
         (cl-mpm/utils:vector-from-list '(1d0 0d0 0d0))
         (cl-mpm/utils:vector-from-list (list
                                         (+ offset (* penalty-width))
                                         (+ height (/ penalty-width 2))
                                         0d0))
         (/ penalty-width 2)
         epsilon
         friction
         0d0))
      (defparameter *penalty-left*
        (cl-mpm/penalty::make-bc-penalty-distance-point
         *sim*
         (cl-mpm/utils:vector-from-list '(-1d0 0d0 0d0))
         (cl-mpm/utils:vector-from-list (list
                                         (- offset (* penalty-width))
                                         (+ height (/ penalty-width 2))
                                         0d0))
         (/ penalty-width 2)
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
          )))
      ;; (let ((normal (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))))
      ;;   (defparameter *penalty*
      ;;     (cl-mpm/penalty::make-bc-penalty-displacment
      ;;      *sim*
      ;;      normal
      ;;      epsilon
      ;;      :clip-function (lambda (p) (cl-mpm/penalty::clip-radial
      ;;                                  p
      ;;                                  normal
      ;;                                  (cl-mpm/utils:vector-from-list (list offset height 0d0))
      ;;                                  penalty-width)))))

      )
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
              (csv-dir nil)
              (enable-plastic t)
              (csv-filename (format nil "load-disp.csv")))
  (unless csv-dir
    (setf csv-dir output-dir))
  (ensure-directories-exist output-dir)
  (ensure-directories-exist csv-dir)
  (let ((eps (cl-mpm/penalty::bc-penalty-epsilon *penalty*)))
    (setf (cl-mpm/penalty::bc-penalty-epsilon *penalty*) 1d-15)
    (cl-mpm/dynamic-relaxation::elastic-static-solution
     *sim*
     :elastic-solver (type-of *sim*))
    (setf (cl-mpm/penalty::bc-penalty-epsilon *penalty*) eps))
  (let* ((lstps 10)
         (total-disp -0.4d0)
         (current-disp 0d0)
         (disp-0
           (cl-mpm::reduce-over-mps (cl-mpm:sim-mps *sim*)
                                          (lambda (mp)
                                            (cl-mpm/utils:get-vector (cl-mpm/particle::mp-displacement mp) :y))
                                          #'min)
                 )
         (step 0))
    (defparameter *data-disp* (list))
    (defparameter *data-load* (list))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (push disp-0 *data-disp*)
    (push (get-load) *data-load*)

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
      :enable-plastic enable-plastic
      :enable-damage t
      :damping 1d0;(sort 2d0)
      :max-damage-inc 1.1d0
      :substeps 50
      :criteria 1d-3
      :save-vtk-dr t
      :save-vtk-loadstep t
      :dt-scale 1d0))))
;; (defun run-gravity (&key (output-dir (format nil "./output/"))
;;               (csv-dir nil)
;;               (enable-plastic t)
;;               (csv-filename (format nil "load-disp.csv")))
;;   (unless csv-dir
;;     (setf csv-dir output-dir))
;;   (ensure-directories-exist output-dir)
;;   (ensure-directories-exist csv-dir)
;;   (let* ((lstps 50)
;;          (gravity -5d0)
;;          (step 0))
;;     (setf (cl-mpm::sim-gravity *sim*) gravity)
;;     (defparameter *data-disp* (list))
;;     (defparameter *data-load* (list))
;;     (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))

;;     (vgplot:close-all-plots)
;;     (time
;;      (cl-mpm/dynamic-relaxation::run-adaptive-load-control
;;       *sim*
;;       :output-dir output-dir
;;       ;:plotter (lambda (sim) (plot-load-disp))
;;       :plotter (lambda (sim) (plot-domain))
;;       :load-steps lstps
;;       :max-adaptive-steps 10
;;       :adaption-constant 2
;;       :enable-plastic enable-plastic
;;       :enable-damage t
;;       :damping 1d0;(sort 2d0)
;;       :max-damage-inc 1.1d0
;;       :substeps 20
;;       :criteria 1d-3
;;       :post-conv-step
;;       (lambda (sim)
;;         (push
;;          (/
;;           (cl-mpm::reduce-over-mps
;;            (cl-mpm:sim-mps *sim*)
;;            (lambda (mp)
;;              (cl-mpm/fastmaths:mag (cl-mpm/particle::mp-displacement mp)))
;;            #'+)
;;           (length (cl-mpm::sim-mps *sim*)))
;;          *data-disp*)
;;         (let ((load (cl-mpm::sim-gravity *sim*)))
;;           (format t "Load ~E~%" load)
;;           (push load *data-load*))
;;         (plot-load-disp)
;;         (save-csv csv-dir csv-filename *data-disp* *data-load*)
;;         (incf step))
;;       :save-vtk-dr t
;;       :save-vtk-loadstep t
;;       :dt-scale 1d0))))

(defun test ()
  (setup :mps 3
         :refine 3
         :enable-fbar t
         :multigrid-refines 0
         :gf 10000d0
         :l-scale 4d0
         :angle 16.7d0
         :angle-r 10d0
         :rt 1d0
         :rc 0d0)
  ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
  ;; (setf (cl-mpm::sim-nonlocal-damage *sim*) t)
  (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
  (run)
  ;; (save-csv "./examples/fbar/rigid-footing/" (format nil "data_fbaradjust_~A.csv" t) *data-disp* *data-load*)
  ;; (dolist (fbar (list t nil))
  ;;   (save-csv "./examples/fbar/rigid-footing/" (format nil "data_fbar_~A.csv" fbar) *data-disp* *data-load*)
  ;;   )
  )




(defun test-degredation ()
  (let ((refine 2)
        (mps 3))
    (dolist (particle (list
                       'cl-mpm/particle::particle-fpd-tcs-strain
                       'cl-mpm/particle::particle-fpd-tcs
                       ;; 'cl-mpm/particle::particle-fpd-spectral
                       ;; 'cl-mpm/particle::particle-fpd-spectral-strain
                       ;; 'cl-mpm/particle::particle-fpd-isotropic
                       ))
      (setup :mps mps
             :refine refine
             :oversize-factor (- 1d0 1d-2)
             :gf 10000d0
             :l-scale 2d0
             :angle 16.7d0
             :angle-r 10d0
             :rt 1d0)
      (cl-mpm::iterate-over-mps
       (cl-mpm:sim-mps *sim*)
       (lambda (mp)
         (change-class mp particle)))
      ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
      (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
      (ignore-errors
       (run :output-dir (format nil "./output-ekl-~D-~A/" refine particle)
                    :enable-plastic nil)))))
