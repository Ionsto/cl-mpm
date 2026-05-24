(defpackage :cl-mpm/examples/octree/column
  (:use
   :cl
   :cl-mpm/example
   :cl-mpm/utils))
(in-package :cl-mpm/examples/octree/column)


(defparameter *sim* nil)
(defun setup (&key
                (refine 1)
                (mps 2)
                (multigrid-refine 1))
  (let* ((E 10d3)
         (L 50d0)
         (h (/ 1 refine))
         (density 80d0)
         (size (list h (+ L (* 2 h))))
         (block-size (list h L))
         (sim (cl-mpm/setup::make-simple-sim
               h
               (mapcar (lambda (x) (* x refine)) size)
               :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-octree-quasi-static
               :args-list
               (list
                :split-factor (* (sqrt 2) (/ 1d0 mps))
                :enable-fbar nil
                :enable-aggregate t
                :max-split-depth 6
                :enable-split nil
                :mp-removal-size nil
                :refinement multigrid-refine
                :gravity -10d0
                )))
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
        'cl-mpm/particle::particle-elastic
        :E E
        :nu 0d0))
      ;; (cl-mpm::domain-sort-mps sim)
      (cl-mpm/setup::set-mass-filter sim density :proportion 1d-3)
      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))
      )
    (defparameter *sim* sim))
  (cl-mpm/setup::setup-bcs
   *sim*
   :left (list 0 nil nil)
   :right (list 0 nil nil)
   :bottom (list 0 0 nil)
   :top (list nil nil nil))
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defun test ()
  (cl-mpm/utils:set-workers 8)
  (setup :mps 4 :refine 0.25
         :multigrid-refine 1)
  (cl-mpm/dynamic-relaxation::set-mesh-default *sim*)
  (let ((h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*))))
    (setf (cl-mpm/dynamic-relaxation::sim-octree-refinement-criteria
           *sim*)
          (lambda (sim mesh c)
            (= (mod (nth-value 0 (round (cl-mpm/utils::varef (cl-mpm/mesh::cell-centroid c) 1) h)) 5)
               0))))
  ;; (cl-mpm::update-sim *sim*)
  ;; (cl-mpm/dynamic-relaxation::save-vtks *sim* "./output/" 0)
  ;; (cl-mpm::update-sim *sim*)
  ;; (cl-mpm/dynamic-relaxation::save-vtks *sim* "./output/" 1)
  ;; (cl-mpm::update-sim *sim*)
  ;; (cl-mpm/dynamic-relaxation::save-vtks *sim* "./output/" 2)
  (cl-mpm/dynamic-relaxation::run-load-control
   *sim*
   :output-dir (format nil "./output/")
   :plotter (lambda (sim) (plot-domain))
   :load-steps 10
   :damping (sqrt 2d0)
   :substeps 50
   :criteria 1d-9
   :conv-steps 5000
   ;; :sub-conv-steps 1000
   :post-iter-step
   (lambda (i o e))
   :save-vtk-dr t
   :save-vtk-loadstep t
   :dt-scale 1d0)
  )
