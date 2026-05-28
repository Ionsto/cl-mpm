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
  (let* ((L 2d0)
         (d 1d0)
         (domain-width 5d0)
         (h (/ 0.125d0 refine))
         (offset l)
         (height 1d0)
         (domain-height (* 2 L))
         (density 1000d0)
         (E 1d6)
         (domain-size (list domain-width domain-height))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list height height)))

    (setf *sim* (cl-mpm/setup::make-simple-sim
                 h
                 element-count
                 ;; :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static
                 ;; :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-octree-quasi-static
                 :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-octree-usf
                 ;; :sim-type 'cl-mpm/aggregate::mpm-sim-agg-usf
                 :args-list (list :enable-aggregate t
                                  :enable-fbar nil
                                  :enable-split nil
                                  :mp-removal-size nil
                                  :intra-mesh-agg nil
                                  :refinement multigrid-refine
                                  )))
    (defparameter *block-size* height)
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0 (+ 0.05d0 offset))
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

    (let ((datum L)
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
    (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-octree-usf)
    (setf (cl-mpm::sim-velocity-algorithm *sim*) :TBLEND)

    (vgplot:close-all-plots)
    (let ((time (list))
          (top-height (list)))
      (cl-mpm/dynamic-relaxation::run-time
       *sim*
       :output-dir output-dir
       :plotter
       (lambda (sim)
         ;; (plot-domain)
         ;; (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step)) :terminal "png size 1920,1080")
         (format t "Top level ~E~%" (/
                                     (- (compute-top-level *sim*) (+ *datum* (* 0.5d0 *block-size*)))
                                     *block-size*))
         (push (cl-mpm::sim-time *sim*) time)
         (push (compute-top-level *sim*) top-height)
         (vgplot:plot time top-height)
         (incf step))
       :total-time 1000d0
       :damping 1d-2
       :dt 0.1d0
       :initial-quasi-static nil
       :dt-scale 0.9d0))))

(defun compute-top-level (sim)
  (cl-mpm::reduce-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (cl-mpm/utils:varef
      (cl-mpm/fastmaths:fast-.+
       (cl-mpm/particle::mp-position mp)
       (cl-mpm/particle::mp-domain-size mp))
      1))
   #'max))

(defun test ()
  (cl-mpm/utils:set-workers 8)
  (setup :mps 4 :refine 1 :multigrid-refine 1)
  (setf
   (cl-mpm/dynamic-relaxation::sim-octree-refinement-criteria *sim*)
   (lambda (sim mesh c)
     (declare (ignore mesh sim))
     (> (cl-mpm/utils::varef (cl-mpm/mesh::cell-centroid c) 0) 0.5d0)))
  ;; (run)
  (run-time)
  )
