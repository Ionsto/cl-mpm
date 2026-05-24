(defpackage :cl-mpm/examples/octree/float
  (:use :cl
   :cl-mpm/example
   :cl-mpm/utils)
  )
(in-package :cl-mpm/examples/octree/float)

(defparameter *datum* 0)
(defparameter *water-bc* nil)
(declaim (notinline plot-domain))
(defun plot-domain ()
  (when *sim*
    (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
           (h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
           (ms-x (first ms))
           (ms-y (second ms)))
      (vgplot:format-plot t "set object 1 rect from 0,0 to ~f,~f fc rgb 'blue' fs transparent solid 0.5 noborder behind" ms-x *datum*)
      (cl-mpm/plotter:simple-plot
       *sim*
       :plot :deformed
       :trial t
       :colour-func (lambda (mp) (cl-mpm/utils:varef (cl-mpm/particle::mp-stress mp) 1))))))

(defun setup (&key (refine 1) (mps 3) (multigrid-refine 0))
  (let* ((L 10d0)
         (d 1d0)
         (domain-width 30d0)
         (h (/ 1d0 refine))
         (offset l)
         (height 10d0)
         (domain-height (* 3 L))
         (density 1000d0)
         (E 1d9)
         (domain-size (list domain-width domain-height))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list height height)))

    (setf *sim* (cl-mpm/setup::make-simple-sim
                 h
                 element-count
                 ;; :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static
                 :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-octree-quasi-static
                 ;:sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-octree-usf
                 ;; :sim-type 'cl-mpm/aggregate::mpm-sim-agg-usf
                 :args-list (list :enable-aggregate t
                                  :enable-fbar nil
                                  :enable-split nil
                                  :mp-removal-size nil
                                  :intra-mesh-agg t
                                  :refinement multigrid-refine
                                  )))
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0 (+ 1d0 offset))
      block-size
      (mapcar (lambda (e) (* (/ e h) mps)) block-size)
      density
      'cl-mpm/particle::particle-elastic
      :E E
      :nu 0.2d0))

    (setf (cl-mpm::sim-gravity *sim*) -1d1)
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)
    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 nil nil))

    (let ((datum (+ offset (* 0.5d0 height)))
          (water-density (* density 2d0))
          (pressure-condition t))
      (defparameter *datum* datum)
      (defparameter *water-bc*
        (if pressure-condition
            (cl-mpm/buoyancy::make-bc-buoyancy-clip
             *sim*
             datum
             water-density
             (lambda (pos datum) t)
             :visc-damping 0d0)
            (cl-mpm/buoyancy::make-bc-buoyancy-body
             *sim*
             datum
             water-density
             (lambda (pos) t)))))
    (cl-mpm:add-bcs-force-list
     *sim*
     *water-bc*)
    (defparameter *current-inc* 0d0)
    ;; (setf
    ;;  (cl-mpm:sim-dt *sim*)
    ;;  (* 0.5d0 (cl-mpm/setup:estimate-elastic-dt *sim*)))
    (setf *run-sim* t))
  (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*))))

(defun run ()
  (let* ((lstps 1))
    (cl-mpm/dynamic-relaxation::run-load-control
     *sim*
     :output-dir (format nil "./output/")
     :plotter (lambda (sim) (plot-domain)
                (let ((top-pos
                        (cl-mpm::reduce-over-mps
                         (cl-mpm::sim-mps *sim*)
                         (lambda (mp)
                           (+
                            (cl-mpm/utils::varef (cl-mpm/particle::mp-domain-size mp) 1)
                            (cl-mpm/utils::varef (cl-mpm/particle::mp-position-trial mp) 1)))
                         #'max)))
                  (format t "Top pos ~F - target ~F~%" top-pos (+ *datum* (* 0.5d0 10d0)))))
     :load-steps lstps
     :damping (sqrt 1d0)
     :substeps 100
     :conv-steps 1000
     :criteria 1d-9
     :save-vtk-dr t
     :save-vtk-loadstep t
     :dt-scale 0.5d0
     ))
  )
(defun run-time (&key (output-dir "./output/"))
  (let ((step 0))
    (setf (cl-mpm::sim-velocity-algorithm *sim*) :TBLEND)
    (cl-mpm/dynamic-relaxation::run-time
     *sim*
     :output-dir output-dir
     :plotter (lambda (sim)
                (plot-domain)
                (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step)) :terminal "png size 1920,1080")
                (incf step))
     :total-time 100d0
     :damping 0.1d0
     :dt 0.1d0
     :initial-quasi-static nil
     :dt-scale 0.5d0)))

(defun test ()
  (cl-mpm/utils:set-workers 8)
  (setup :mps 4 :refine 1 :multigrid-refine 0)
  (setf
   (cl-mpm/dynamic-relaxation::sim-octree-refinement-criteria *sim*)
   (lambda (sim mesh c)
     ;; nil
     (and
      (> (cl-mpm/utils::varef (cl-mpm/mesh::cell-centroid c) 0) 4d0))
     ))
  (run)
  ;; (run-time)
  )
