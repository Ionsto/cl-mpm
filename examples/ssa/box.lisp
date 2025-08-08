(defpackage :cl-mpm/examples/ssa/box
  (:use :cl
   :cl-mpm/example))
(in-package :cl-mpm/examples/ssa/box)
(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)


(defun setup (&key (refine 1) (mps 2))
  (let* ((density 916.7d0)
         (mesh-resolution (/ 10d0 refine))
         (domain-size (list 1000d0 200d0))
         (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
         (block-size (list 500d0 100d0)))
    (defparameter *sim* (cl-mpm/setup::make-simple-sim mesh-resolution element-count
                                                       :sim-type
                                                       'cl-mpm/ssa::mpm-sim-ssa
                                                       :args-list (list
                                                                   :enable-fbar t)
                                                       ))
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0 0)
      block-size
      (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
      (* density 1d3)
      'cl-mpm/particle::particle-ssa-glen
      :E 1d9
      :nu 0.325d0
      :visc-factor 1d6
      ;; :height 1d3
      ))
    (cl-mpm/setup::setup-bcs
     *sim*
     :top '(nil 0 nil)
     :bottom '(nil 0 nil))
    (setf (cl-mpm:sim-mass-scale *sim*) 1d2)
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-9)
    (setf (cl-mpm::sim-velocity-algorithm *sim*) :PIC)
    (setf *run-sim* t)
    (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
   ))

(defun run ()
  (let ((dt-scale 0d0)
        (target-time 1d0))
    )
  )
(defun run (&key (output-dir "./output/"))
  (uiop:ensure-all-directories-exist (list (uiop:merge-pathnames* output-dir)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))

  (let ((dt-scale 0.5d0))
    ;; (setf (cl-mpm:sim-mass-scale *sim*) 1d0)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 1d-1
             (sqrt (cl-mpm:sim-mass-scale *sim*))
             (cl-mpm/setup:estimate-critical-damping *sim*)))
    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup:estimate-elastic-dt *sim* :dt-scale dt-scale))
    (let* ((target-time 1d0)
           (substeps (ceiling target-time (cl-mpm:sim-dt *sim*))))
      (format t "Substeps ~D~%" substeps)
      (loop for step from 0 below 1000
            while *run-sim*
            do
               (format t "Step ~D~%" step)
               (time
                (dotimes (i substeps)
                  (cl-mpm:update-sim *sim*)
                  ))
               (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir (format nil "sim_~5,'0d.vtk" step)) *sim* )
               (cl-mpm/output:save-vtk-nodes (uiop:merge-pathnames* output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) *sim* )
               (plot-domain)
               (swank.live:update-swank)
               (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step))
                                  :terminal "png size 1920,1080"
                                  )
            ))))
