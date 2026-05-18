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
         (block-size (list 8 height)))
    (setf
     *sim*
     (cl-mpm/setup::make-simple-sim
      h
      element-count
      :sim-type
      'cl-mpm/dynamic-relaxation::mpm-sim-octree-usf
      ;; 'cl-mpm/dynamic-relaxation::mpm-sim-octree-damage-quasi-static
      :args-list
      (list
       :enable-aggregate t
       :intra-mesh-agg t
       :enable-split t
       :enable-fbar enable-fbar
       :refinement 1
       :vel-algo :TBLEND
       )))
    (setf h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
    (let ((h h))
      (cl-mpm:add-mps
       *sim*
       (cl-mpm/setup:make-block-mps
        (list 0d0 0d0)
        block-size
        (mapcar (lambda (e) (* (/ e h) mps)) block-size)
        density
        'cl-mpm/particle::particle-vm-implicit
        :E E
        :nu nu
        :rho 10d3

        ;; 'cl-mpm/particle::particle-elastic-damage-delayed
        ;; :E 1d6
        ;; :nu 0.3d0
        ;; :initiation-stress 1d4
        ;; :local-length (* h)
        ;; :ductility 10d0
        ;; :delay-time 1000d0
        ;; :delay-exponent 1d0
        ))

      ;; (cl-mpm:add-mps
      ;;  *sim*
      ;;  (cl-mpm/setup:make-block-mps
      ;;   (list 8d0 0d0)
      ;;   block-size
      ;;   (mapcar (lambda (e) (* (/ e h) mps)) block-size)
      ;;   density
      ;;   'cl-mpm/particle::particle-elastic
      ;;   :E E
      ;;   :nu nu))
      )

    (setf (cl-mpm::sim-gravity *sim*) -10d0)
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-3)
    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 nil nil)
     :right '(0 nil nil)
     :bottom '(nil 0 nil)))
  (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
  )



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
  (dotimes (i (1+ (cl-mpm::sim-multigrid-refinement sim)))
    (setf (cl-mpm::sim-mesh sim) (aref (cl-mpm::sim-mesh-list sim) i))
    (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~D_~5,'0d.vtk" i step)) sim)
    (cl-mpm/output::save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~D_~5,'0d.vtk" i step)) sim))

  (setf (cl-mpm::sim-mesh sim) (aref (cl-mpm::sim-mesh-list sim) 0))
  
  )

(defun refinement-criteria (mesh c)
  (let ((d-damage (cl-mpm/utils:vector-zeros))
        (damage 0d0)
        (volume 0d0)
        )
    (cl-mpm/aggregate::iterate-over-cell-shape-local
     mesh
     c
     (cl-mpm/mesh::cell-centroid c)
     (lambda (node weight grads)
       ;; (let ((dsvp (cl-mpm/utils::arb-matrix-from-list (list (cl-mpm/utils::gradient-dx grads) 0d0 0d0
       ;;                                                       0d0 (cl-mpm/utils::gradients-dy grads) 0d0
       ;;                                                       0d0 0d0 (cl-mpm/utils::gradients-dz grads))
       ;;                                                 3 3))))
       ;; (cl-mpm/fastmaths::fast-.+
       ;;  (cl-mpm/fastmaths::fast-@-arb-vector
       ;;   )
       ;;  )
       (incf damage (* weight (cl-mpm/mesh::node-damage node)))
       (incf volume (* weight (cl-mpm/mesh::node-volume node)))
       ))
    (if (> volume 0d0)
        (> (/ damage volume) 0.5d0)
        nil)))

(defun test (&key (output-dir "./output/"))
  (cl-mpm/utils::set-workers 8)
  (setup :refine 2)
  ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-octree-damage-quasi-static)
  (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-octree-quasi-static)
  ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-octree-damage-usf)
  ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-octree-usf)
  ;; (setf (cl-mpm::sim-velocity-algorithm *sim*) :QUASI-STATIC)
  (let* ((h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
         (transition-pos 4d0)
         (upper-pos (+ 4d0 h)))
    (setf (cl-mpm/dynamic-relaxation::sim-octree-refinement-criteria
           *sim*)
          ;; #'refinement-criteria
          (lambda (mesh c)
            (and
             (> (cl-mpm/utils::varef (cl-mpm/mesh::cell-centroid c) 0) 4)
             ))
          ))

  ;; (cl-mpm:update-sim *sim*)
  ;; (cl-mpm/dynamic-relaxation::save-vtks *sim* "./output/" 0)
  ;; (cl-mpm/dynamic-relaxation::save-vtks *sim* "./output/" 1)
  ;; (cl-mpm/dynamic-relaxation::run-quasi-time
  ;;  *sim*
  ;;  :output-dir "./output/"
  ;;  :dt-scale 1d0
  ;;  :plotter (lambda (sim) (plot-domain))
  ;;  :conv-criteria 1d-3
  ;;  ;; :min-adaptive-steps -8
  ;;  ;; :max-adaptive-steps 8
  ;;  ;; :adaption-constant 2
  ;;  :substeps 20
  ;;  :dt 100d0
  ;;  :total-time 100000d0
  ;;  )
  ;; (cl-mpm/dynamic-relaxation::run-time
  ;;  *sim*
  ;;  :output-dir output-dir
  ;;  :plotter (lambda (sim)
  ;;             (plot-domain))
  ;;  :total-time 100d0
  ;;  :damping 0.1d0
  ;;  :dt 0.1d0
  ;;  :initial-quasi-static nil
  ;;  :dt-scale 0.25d0)
  (cl-mpm/dynamic-relaxation::run-load-control
   *sim*
   :output-dir (format nil "./output/")
   :plotter (lambda (s) (plot-domain))
   :load-steps 40
   :damping 1d0
   :substeps 100
   :criteria 1d-3
   :save-vtk-dr t
   :save-vtk-loadstep t
   :dt-scale 0.9d0)
  )




