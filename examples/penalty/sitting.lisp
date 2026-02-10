(defpackage :cl-mpm/examples/penalty/sitting
  (:use :cl :cl-mpm/example :cl-mpm/utils))
(in-package :cl-mpm/examples/penalty/sitting)
(declaim (optimize (debug 3) (safety 3) (speed 0)))

(defun setup (&key (refine 1) (mps 2) (mu 0.25d0)
                (eps-scale 1d1)
                )
  (let* ((d 1d0)
         (dsize 2d0)
         (density 1d3)
         (mesh-resolution (/ 0.25d0 refine))
         (offset (* 2d0 mesh-resolution))
         (domain-size (list dsize (+ offset (* 2 d))))
         (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
         (block-size (list d d))
         (E 1d6))
    (setf *sim* (cl-mpm/setup::make-simple-sim
                 mesh-resolution
                 element-count
                 :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
                 :args-list (list :enable-aggregate nil
                                  :enable-fbar nil)))
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0d0 offset)
      block-size
      (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
      density
      'cl-mpm/particle::particle-elastic
      :E E
      :nu 0.2d0
      :index 0))

    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)

    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 nil nil))

    (let ()
      (defparameter *floor*
        (cl-mpm/penalty::make-bc-penalty-distance-point  
         *sim*
         (cl-mpm/utils:vector-from-list (list 0d0 1d0 0d0))
         (cl-mpm/utils:vector-from-list (list (* dsize 0.5 ) offset 0d0))
         (* dsize 0.5d0)
         (* E eps-scale)
         mu
         0d0
         )))
    (cl-mpm::add-bcs-force-list
     *sim*
     *floor*)
    (setf *run-sim* t))
  (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))

  )

(defun run (&key (output-dir (format nil "./output/")))
  (vgplot:close-all-plots)
  (cl-mpm/dynamic-relaxation::run-load-control
   *sim*
   :output-dir output-dir
   :plotter (lambda (sim)
              (plot-domain))
   :load-steps 1
   :kinetic-damping nil
   :damping 1d0;(sqrt 2d0)
   :substeps 50
   :criteria 1d-2
   :save-vtk-dr t
   :save-vtk-loadstep t
   :dt-scale 1d0)
  )

(defun test ()
  (let ((mu 0.1d0))
    (dolist (r (list 1 2 4))
      (dolist (mps (list 2 4))
        (setup :mps mps :refine r :mu mu
               :eps-scale 1d2)
        (run :output-dir (format nil "./output-~D-~D/" r mps))))))
