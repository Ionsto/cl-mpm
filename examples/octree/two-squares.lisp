(defpackage :cl-mpm/examples/octree/two-squares
  (:use :cl
   :cl-mpm/example
   :cl-mpm/utils)
  )
(in-package :cl-mpm/examples/octree/two-squares)
(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)
(setf cl-mpm/settings:*optimise-setting* cl-mpm/settings::*optimise-debug*)
(declaim (notinline plot-domain))
(defun plot-domain ()
  (when *sim*
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :fill nil
     :trial t
     ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-index mp))
     :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp))
     )))

(declaim (notinline setup))
(defun setup (&key
                (refine 1)
                (mps 3)
                (enable-fbar nil)
                (multigrid-refines 0)
                )
  (let* ((domain-width 16d0)
         (height 4d0)
         (h (/ 1d0 refine))
         (domain-height (+ domain-width))
         (density 1d3)
         (E 1d6)
         (nu 0.3d0)
         (domain-size (list 16 8))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list 16 4)))
    (setf
     *sim*
     (cl-mpm/setup::make-simple-sim
      h
      element-count
      :sim-type
      ;; 'cl-mpm/dynamic-relaxation::mpm-sim-octree-usf
      ;; 'cl-mpm/dynamic-relaxation::mpm-sim-octree-damage-quasi-static
      'cl-mpm/dynamic-relaxation::mpm-sim-octree-quasi-static
      :args-list
      (list
       :enable-aggregate t
       :intra-mesh-agg t
       :enable-split t
       :split-factor (* 1.01d0 (sqrt 2) (/ 1d0 mps))
       :max-split-depth 10
       :refinement multigrid-refines)))
    (setf h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
    (let ((h h))
      (cl-mpm:add-mps
       *sim*
       (cl-mpm/setup:make-block-mps
        (list 0d0 0d0)
        block-size
        (mapcar (lambda (e) (* (/ e h) mps)) block-size)
        density
        'cl-mpm/particle::particle-elastic
        :E E
        :nu nu
        ;; 'cl-mpm/particle::particle-vm-implicit
        ;; :E E
        ;; :nu nu
        ;; :rho 10d3
        )))

    (setf (cl-mpm::sim-gravity *sim*) -10d0)
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)
    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 nil nil)
     :right '(0 nil nil)
     :bottom '(nil 0 nil)))
  (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*))))



(defun save-csv (output-dir filename data-disp data-load)
  (with-open-file (stream (merge-pathnames filename output-dir) :direction :output :if-exists :supersede)
    (format stream "disp,load~%")
    (loop for disp in data-disp
          for load in data-load
          do (format stream "~E,~E~%" (float disp 0e0) (float load 0e0)))))
(declaim (notinline run))

(defmethod cl-mpm/dynamic-relaxation::save-vtks ((sim cl-mpm/dynamic-relaxation::mpm-sim-octree) output-dir step)
  (cl-mpm:sim-format sim t "Save vtks ~D~%" step)
  (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" step)) sim)
  ;; (break)
  (dotimes (i (1+ (cl-mpm::sim-multigrid-refinement sim)))
    (setf (cl-mpm::sim-mesh sim) (aref (cl-mpm::sim-mesh-list sim) i))
    (cl-mpm/output::save-vtk-mesh-nodes (merge-pathnames output-dir (format nil "sim_nodes_~D_~5,'0d.vtk" i step)) (aref (cl-mpm::sim-mesh-list sim) i))
    (cl-mpm/output::save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~D_~5,'0d.vtk" i step)) sim))
  (setf (cl-mpm::sim-mesh sim) (aref (cl-mpm::sim-mesh-list sim) 0))
  )

(defun refinement-criteria (sim mesh c)
  (let ((d-damage (cl-mpm/utils:vector-zeros))
        (damage 0d0)
        (volume 0d0))
    (cl-mpm/damage::iterate-over-point-neighbour-mps
     (aref (cl-mpm::sim-mesh-list sim) 0)
     (cl-mpm/mesh::cell-centroid c)
     (lambda (mesh mp dist)
       (setf damage (max (/ (cl-mpm/particle::mp-damage-ybar mp)
                            (cl-mpm/particle::mp-initiation-stress mp)) damage))))
    (> damage 0.95d0)))

(defun test (&key (output-dir "./output/"))
  (cl-mpm/utils::set-workers 8)
  (setup :mps 4 :refine 1 :multigrid-refines 2)
  (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-octree-quasi-static)

  (setf (cl-mpm/dynamic-relaxation::sim-octree-refinement-criteria
         *sim*)
        (lambda (sim mesh c)
          (and
           (>
            (cl-mpm/utils::varef (cl-mpm/mesh::cell-centroid c) 1)
               2d0)
           (<
            (abs (- (cl-mpm/utils::varef (cl-mpm/mesh::cell-centroid c) 0)
                   8d0))
            2))
          ;; (and
          ;;  (>
          ;;   (cl-mpm/utils::varef (cl-mpm/mesh::cell-centroid c) 0)
          ;;   (* ;; (1+ (cl-mpm/dynamic-relaxation::cell-mesh-index c))
          ;;        8d0)))
          ))

  ;; (cl-mpm:update-sim *sim*)
  ;; (cl-mpm/dynamic-relaxation::save-vtks *sim* "./output/" 0)
  ;; (cl-mpm/dynamic-relaxation::save-vtks *sim* "./output/" 1)
  (cl-mpm/setup::set-mass-filter *sim* 1d3 :proportion 0d-15)
  (cl-mpm/dynamic-relaxation::run-load-control
   *sim*
   :output-dir (format nil "./output/")
   :plotter (lambda (s) (plot-domain))
   :load-steps 1
   :damping (sqrt 2d0)
   :substeps 50
   :conv-steps 1000
   :criteria 1d-9
   :save-vtk-dr t
   :save-vtk-loadstep nil
   :dt-scale 0.9d0)
  ;; (cl-mpm/dynamic-relaxation::run-time
  ;;  *sim*
  ;;  :output-dir output-dir
  ;;  :plotter (lambda (sim)
  ;;             (plot-domain)
  ;;             )
  ;;  :total-time 100d0
  ;;  :damping 0.1d0
  ;;  :dt 0.01d0
  ;;  :initial-quasi-static nil
  ;;  :dt-scale 0.5d0)
  )





(defun disp-test ()
  (let* ((d 1)
         (et (cl-mpm/aggregate::sim-global-sparse-et *sim*))
         (e (cl-mpm/aggregate::sim-global-sparse-e *sim*))
         (m (cl-mpm/aggregate::sim-global-sparse-ma *sim*))
         (bcs-int (aref (cl-mpm/aggregate::sim-global-bcs-int *sim*) d))
         (bcs (aref (cl-mpm/aggregate::sim-global-bcs *sim*) d))
         (disp (cl-mpm/utils::arb-matrix 4 1))
         )
    ;; (pprint (magicl:nrows bcs))
    (cl-mpm::reset-node-displacement *sim*)
    (loop for i from 0 to 10
          do
             (progn
               (cl-mpm/aggregate::project-global-vec
                *sim*
                (cl-mpm/aggregate::extend-vec
                 *sim*
                 disp
                 d)
                #'cl-mpm/mesh::node-displacment
                d
                )
               (cl-mpm/dynamic-relaxation::save-vtks *sim* "./output/" i)
               ;; (incf (varef disp 1) 0.1d0)
               ;; (incf (varef disp 5) 0.1d0)
               ;; (incf (varef disp 15) 0.1d0)
               (incf (varef disp 1) 0.1d0)
               ;; (incf (varef disp 19) 0.1d0)
               ;; (incf (varef disp 17) 0.1d0)
               ;; (incf (varef disp 3) 0.1d0)
               ))
    ;; (pprint
    ;;  ; (cl-mpm/utils::diagonal-to-arb m)
    ;;  (cl-mpm/fastmaths::fast-@-arb-arb
    ;;   (cl-mpm/fastmaths::fast-@-arb-arb
    ;;    (cl-mpm/implicit::reduce-with-bcs
    ;;     (cl-mpm/utils::sparse-to-mat et)
    ;;     bcs-int
    ;;     bcs
    ;;     )
    ;;    (cl-mpm/implicit::reduce-with-bcs
    ;;     (cl-mpm/utils::diagonal-to-arb m)
    ;;     bcs
    ;;     bcs)
    ;;    ;; :multithreaded t
    ;;    )
    ;;   (cl-mpm/implicit::reduce-with-bcs
    ;;    (cl-mpm/utils::sparse-to-mat e)
    ;;    bcs
    ;;    bcs-int
    ;;    )
    ;;   :multithreaded t)
    ;;  )
    ))
