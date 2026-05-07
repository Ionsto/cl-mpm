(defpackage :cl-mpm/examples/beam
  (:use :cl
   :cl-mpm/example
   :cl-mpm/utils
   ))
(in-package :cl-mpm/examples/beam)
(declaim (notinline plot-domain))

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-elastic) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  (cl-mpm::update-domain-stretch mesh mp dt)
  )

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
         (dsize 11d0)
         (density 7750d3)
         (mesh-resolution (/ 0.25d0 refine))
         (domain-size (list dsize dsize))
         (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
         (depth 1d0)
         (offset-y (- 11 (* 2 d)))
         (center-y (+ offset-y (/ d 2)))
         (block-size (list L d)))
    (setf *sim* (cl-mpm/setup::make-simple-sim
                 mesh-resolution
                 element-count
                 ;; :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-multi grid
                 :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static
                 :args-list (list :enable-aggregate nil
                                  :mass-update-count 1
                                  :mp-removal-size nil
                                  ;; :refinement 2
                                  ;; :ghost-factor (* 12d6 1d-4)
                                  :ghost-factor nil
                                  )))
    (setf mesh-resolution (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
    (format t "Mesh resolution ~E~%" mesh-resolution)
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

    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-9)
    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 0 nil))
    (setf *run-sim* t))
  (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
  )

(defun run (&key (output-dir "./output/"))
  (cl-mpm/dynamic-relaxation::run-load-control
   *sim*
   :output-dir output-dir
   :plotter (lambda (sim) (plot-domain))
   :load-steps 50
   :damping (sqrt 2d0)
   :dt-scale 1d0
   :substeps 1000
   :criteria 1d-3
   :sub-conv-steps 1000
   :save-vtk-dr nil
   :save-vtk-loadstep t))

(defun test ()
  (setup :mps 3 :refine 1)
  (run))
