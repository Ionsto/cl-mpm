(defpackage :cl-mpm/examples/beam
  (:use :cl
   :cl-mpm/example
   :cl-mpm/utils)
  )
(in-package :cl-mpm/examples/beam)
(defun plot-domain ()
  (when *sim*
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :colour-func (lambda (mp) (cl-mpm/particle::mp-index mp)))))

(defun setup (&key (refine 1) (mps 3))
  (let* ((L 10d0)
         (d 1d0)
         (dsize 11d0)
         (density 7750d3)
         (mesh-resolution (/ 0.25d0 refine))
         (domain-size (list dsize dsize))
         (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
         (depth 1d0)
         (offset-y (- 11 (* 2 d)))
         (center-y (+ offset-y (/ d 2)))
         (block-size (list L d)))
    (setf *sim* (cl-mpm/setup::make-simple-sim mesh-resolution element-count
                                               :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
                                               :args-list (list :enable-aggregate nil)))
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0 offset-y)
      block-size
      (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
      density
      'cl-mpm/particle::particle-elastic
      :E 12d6
      :nu 0.2d0
      :index 0
      :gravity-axis (cl-mpm/utils:vector-zeros)))
    (let ((mass 0d0))
      (cl-mpm:iterate-over-mps
       (cl-mpm:sim-mps *sim*)
       (lambda (mp)
         (when
             (and
              (<= (abs (- (+ (varef (cl-mpm/particle:mp-position mp) 0)
                             (* 0.5d0 (varef (cl-mpm/particle::mp-domain-size mp) 0)))
                          L)) 1d-3)
              (< (abs (-  (varef (cl-mpm/particle:mp-position mp) 1) center-y)) (+ 1d-3 (* 1d0 (varef (cl-mpm/particle::mp-domain-size mp) 1)))))
           (setf mass (cl-mpm/particle::mp-mass mp))
           (setf
            ;; (cl-mpm/particle::mp-gravity mp) (/ -100d3 (* 2d0 (cl-mpm/particle::mp-mass mp)))
            (cl-mpm/particle::mp-index mp) 1)
           (setf (cl-mpm/particle::mp-gravity-axis mp) (cl-mpm/utils:vector-from-list (list 0d0 1d0 0d0)))
           )))
      (setf (cl-mpm:sim-gravity *sim*) (/ -100d3 (* 2d0 mass))))

    (cl-mpm/setup::set-mass-filter *sim* 1d0 :proportion 1d-15)
    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 0 nil))

    (setf (cl-mpm:sim-dt *sim*)
          (* 0.5d0 (cl-mpm/setup:estimate-elastic-dt *sim*)))
    (setf *run-sim* t)))

(defun run ()
  (cl-mpm/dynamic-relaxation::run-load-control
   *sim*
   :output-dir (format nil "./output/")
   :plotter (lambda (sim) (plot-domain))
   :load-steps 25
   :kinetic-damping nil
   :damping 1d0
   :substeps 100
   :criteria 1d-3
   :save-vtk-dr t
   :save-vtk-loadstep t
   :dt-scale 1d0))

(defun test ()
  (setup :mps 2)
  (run))
