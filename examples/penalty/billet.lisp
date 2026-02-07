(defpackage :cl-mpm/examples/billet
  (:use :cl
   :cl-mpm/example
   :cl-mpm/utils)
  )
(in-package :cl-mpm/examples/billet)
(declaim (notinline plot-domain))
(defun plot-domain ()
  (when *sim*
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :trial t
     :colour-func (lambda (mp) (cl-mpm/particle::mp-index mp)))))

(defun setup (&key (refine 1) (mps 3))
  (let* ((L 10d0)
         (d 1d0)
         (domain-width 30d0)
         (h (/ 1d0 refine))
         (height 10d0)
         (domain-height (+ height (* 2 h)))
         (density 7750d3)
         (E 1d6)
         (domain-size (list domain-width domain-height))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list height height)))

    (setf *sim* (cl-mpm/setup::make-simple-sim
                 h
                 element-count
                 :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
                 :args-list (list :enable-aggregate t
                                  :enable-fbar t
                                  :enable-split t
                                  :vel-algo :QUASI-STATIC
                                  )))
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0 0)
      block-size
      (mapcar (lambda (e) (* (/ e h) mps)) block-size)
      density
      'cl-mpm/particle::particle-vm
      :E E
      :nu 0.3d0
      :rho 20d3
      ;; 'cl-mpm/particle::particle-elastic
      ;; :E E
      ;; :nu 0.2d0
      ;; :index 0
      :gravity-axis (cl-mpm/utils:vector-zeros)))
    (cl-mpm/setup::set-mass-filter *sim* 1d0 :proportion 1d-15)
    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 nil nil))

    (defparameter *penalty*
      (cl-mpm/penalty::make-bc-penalty-distance-point
       *sim*
       (cl-mpm/utils:vector-from-list '(0d0 -1d0 0d0))
       (cl-mpm/utils:vector-from-list (list (* 0.5d0 domain-width)
                                            height
                                            0d0))
       (/ domain-width 2)
       (* E 1d2)
       0.5d0
       0d0))
    (defparameter *current-inc* 0d0)
    (cl-mpm:add-bcs-force-list
     *sim*
     *penalty*)
    (setf (cl-mpm:sim-dt *sim*)
          (* 0.5d0 (cl-mpm/setup:estimate-elastic-dt *sim*)))
    (setf *run-sim* t))
  (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
  )

(defun run ()
  (let* ((lstps 20)
         (total-disp -5d0)
         (delta (/ total-disp lstps)))
    (cl-mpm/dynamic-relaxation::run-load-control
     *sim*
     :output-dir (format nil "./output/")
     :plotter (lambda (sim) (plot-domain))
     :loading-function (lambda (i)
                         (cl-mpm/penalty::bc-increment-center
                          *penalty*
                          (cl-mpm/utils:vector-from-list (list 0d0 delta 0d0))))
     :load-steps lstps
     :kinetic-damping nil
     :damping 1d0
     :substeps 50
     :criteria 1d-3
     :save-vtk-dr t
     :save-vtk-loadstep t
     :dt-scale 0.9d0
     )))

(defun test ()
  (setup :mps 2 :refine 1)
  (run))
