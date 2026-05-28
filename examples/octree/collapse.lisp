(defpackage :cl-mpm/examples/octree/collapse
  (:use
   :cl
   :cl-mpm/example
   :cl-mpm/utils))
(in-package :cl-mpm/examples/octree/collapse)


(defparameter *sim* nil)
(defun setup (&key (refine 1)
                (mps 2)
                (multigrid-refine 1))
  (let* ((E 1d6)
         (density 1d3)
         (block-size (list 8d0 8d0 8d0))
         (size (list 16d0 16d0 16d0))
         (sim (cl-mpm/setup::make-simple-sim
               (/ 1d0 refine)
               (mapcar (lambda (x) (* x refine)) size)
               ;; :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-octree-quasi-static
               :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-octree-implicit-dynamic
               :args-list
               (list
                :split-factor (* (sqrt 2) (/ 1d0 mps))
                :enable-fbar nil
                :enable-aggregate t
                :mass-update-count 1
                :damping-update-count 1
                :max-split-depth 6
                :enable-split t
                :mp-removal-size nil
                :refinement multigrid-refine
                :gravity -10d0)))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim))))
    (declare (double-float h density))
    (progn
      (cl-mpm::add-mps
       sim
       (cl-mpm/setup::make-block-mps
        (list 0d0 0d0 0d0)
        block-size
        (mapcar (lambda (e) (* (/ e h) mps)) block-size)
        density
        ;; 'cl-mpm/particle::particle-elastic
        ;; :E E
        ;; :nu 0d0
        'cl-mpm/particle::particle-vm-implicit
        :E E
        :nu 0.3d0
        :rho 20d3))
      (cl-mpm::domain-sort-mps sim)
      (cl-mpm/setup::set-mass-filter sim density :proportion 1d-15)
      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))
      )
    (defparameter *sim* sim))
  (cl-mpm/setup::setup-bcs
   *sim*
   :left (list 0 nil nil)
   :bottom (list nil 0 nil)
   :top (list nil nil nil))
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defun test ()
  (cl-mpm/utils:set-workers 16)
  (setup :mps 2 :refine 0.5
         :multigrid-refine 0)

  (setf
   (cl-mpm/dynamic-relaxation::sim-octree-refinement-criteria *sim*)
   (lambda (sim mesh c)
     (and
      (> (cl-mpm/utils::varef (cl-mpm/mesh::cell-centroid c) 1) 4))))

  (cl-mpm/dynamic-relaxation::run-load-control
   *sim*
   :output-dir (format nil "./output-1/")
   :plotter (lambda (sim) (plot-domain))
   :load-steps 10
   :damping (sqrt 2d0)
   :substeps 20
   :criteria 1d-3
   :conv-steps 5000
   ;; :sub-conv-steps 1000
   :post-iter-step
   (lambda (i o e))
   :save-vtk-dr t
   :save-vtk-loadstep t
   :dt-scale 1d0)
  )

(defun test-real-time ()
  (cl-mpm/utils:set-workers 8)
  (setup :mps 3 :refine 1 :multigrid-refine 1)
  (let ((output-dir (format nil "./output/")))
    (vgplot:close-all-plots)
    (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-octree-implicit-dynamic)

    (when (typep *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-octree)
      (setf (cl-mpm/dynamic-relaxation::sim-octree-refinement-criteria *sim*)
            (lambda (sim mesh c)
              (and
               ;; (= (cl-mpm/dynamic-relaxation::cell-mesh-index c) 0)
               (< (cl-mpm/utils::varef (cl-mpm/mesh::cell-centroid c) 1)
                  4d0)))))

    (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t)
    (setf (cl-mpm::sim-velocity-algorithm *sim*) :TBLEND)
    (cl-mpm/setup::set-mass-filter *sim* 1d3 :proportion 1d-9)
    ;; (cl-mpm:update-sim *sim*)
    ;; (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (let ((step 0))
      (cl-mpm/dynamic-relaxation::run-time
       *sim*
       :output-dir output-dir
       :plotter (lambda (sim)
                  (plot-domain)
                  (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step)) :terminal "png size 1920,1080")
                  (incf step))
       :total-time 100d0
       :damping 1d-2
       :dt 0.1d0
       :initial-quasi-static nil
       :dt-scale 10d0))))
