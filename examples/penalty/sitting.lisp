(defpackage :cl-mpm/examples/penalty/sitting
  (:use :cl :cl-mpm/example :cl-mpm/utils))
(in-package :cl-mpm/examples/penalty/sitting)
(declaim (optimize (debug 3) (safety 3) (speed 0)))

(defun setup (&key (refine 1) (mps 2) (mu 0.25d0)
                (eps-scale 1d1)
                )
  (let* ((d 1d0)
         (dsize 2d0)
         (density 1d4)
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
                 :args-list (list :enable-aggregate t
                                  :enable-fbar t
                                  :enable-split t
                                  :split-factor (* 2d0 (/ 1d0 mps))
                                  )))
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
      ;; 'cl-mpm/particle::particle-vm
      ;; :E E
      ;; :nu 0.2d0
      ;; :rho 1d3
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

(defun run-time (&key (output-dir (format nil "./output/")))
  (setup :refine 1 :eps-scale 1d2)
  (change-class *sim* 'cl-mpm/aggregate::mpm-sim-agg-usf)
  (setf (cl-mpm::sim-velocity-algorithm *sim*) :BLEND)
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (let ((step 0))
    (cl-mpm/dynamic-relaxation::run-time
     *sim*
     :output-dir output-dir
     :plotter (lambda (sim)
                (plot-domain)
                (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step)) :terminal "png size 1920,1080")
                (incf step)
                )
     :total-time 30d0
     :damping 1d-3
     :dt 0.1d0
     :initial-quasi-static nil
     :dt-scale 0.5d0))
  )

(defun run (&key (output-dir (format nil "./output/")))
  (vgplot:close-all-plots)
  (cl-mpm/dynamic-relaxation::run-load-control
   *sim*
   :output-dir output-dir
   :plotter (lambda (sim)
              (plot-domain))
   :load-steps 10
   :kinetic-damping nil
   :damping 1d0;(sqrt 2d0)
   :substeps 10
   :criteria 1d-6
   :save-vtk-dr t
   :save-vtk-loadstep t
   :dt-scale 1d0)
  )

(defun test ()
  (let ((mu 0.9d0))
    (dolist (r (list 8))
      (dolist (mps (list 2))
        (setup :mps mps :refine r :mu mu
               :eps-scale 1d2)
        (run :output-dir (format nil "./output-~D-~D/" r mps))))))


