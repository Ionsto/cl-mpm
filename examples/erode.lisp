(defpackage :cl-mpm/examples/erode
  (:use :cl
        :cl-mpm/example))
(in-package :cl-mpm/examples/erode)

(defun plot-domain ()
  (when *sim*
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :colour-func #'cl-mpm/particle::mp-damage
     )))

(defclass bc-water-damage (cl-mpm/buoyancy::bc-scalar)
  ((damage-rate
    :initform 1d0
    :initarg :damage-rate
    :accessor bc-water-damage-damage-rate)))

(defun make-bc-water-damage (sim datum rate)
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (with-accessors ((h cl-mpm/mesh:mesh-resolution))
        mesh
      (make-instance 'bc-water-damage
                     :index nil
                     :damage-rate rate
                     :damage-volume t
                     :sim sim))))
(defmethod cl-mpm/bc::apply-bc ((bc bc-water-damage) node mesh dt)
  (call-next-method)
  ;; (break)
  (when (cl-mpm::sim-enable-damage *sim*)
    (with-accessors ((sim cl-mpm/buoyancy::bc-buoyancy-sim))
        bc
      (when (cl-mpm::sim-enable-damage sim)
        (loop for mp across (cl-mpm:sim-mps sim)
              do
                 (let ((weathering 0d0))
                   (cl-mpm:iterate-over-neighbours
                    mesh
                    mp
                    (lambda (mesh mp node svp grads fsvp fgrads)
                      (incf weathering (* svp (cl-mpm/mesh::node-boundary-scalar node)))))
                   (setf (cl-mpm/particle::mp-boundary mp)
                         weathering)
                   (incf
                    ;(cl-mpm/particle::mp-history-stress mp)
                    (cl-mpm/particle::mp-damage mp)
                    (abs (*
                          (bc-water-damage-damage-rate bc)
                          weathering dt)))
                   (setf
                    (cl-mpm/particle::mp-damage mp)
                    (min
                      (cl-mpm/particle::mp-damage mp)
                      1d0))))))))

(defclass bc-erode (cl-mpm/buoyancy::bc-scalar)
  ((damage-rate
    :initform 1d0
    :initarg :damage-rate
    :accessor bc-water-damage-damage-rate)))

(defun make-bc-erode (sim datum rate)
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (with-accessors ((h cl-mpm/mesh:mesh-resolution))
        mesh
      (make-instance 'bc-erode
                     :index nil
                     :damage-rate rate
                     :damage-volume t
                     :clip-func (lambda (pos) (and (> (cl-mpm/utils:varef pos 1) datum)))
                     :sim sim))))
(defmethod cl-mpm/bc::apply-bc ((bc bc-erode) node mesh dt)
  (call-next-method)
  ;; (break)
  (with-accessors ((sim cl-mpm/buoyancy::bc-buoyancy-sim))
      bc
    (when (cl-mpm::sim-enable-damage sim)
      (loop for mp across (cl-mpm:sim-mps sim)
            do
               (let ((weathering 0d0))
                 (cl-mpm:iterate-over-neighbours
                  mesh
                  mp
                  (lambda (mesh mp node svp grads fsvp fgrads)
                    (incf weathering (* svp (cl-mpm/mesh::node-boundary-scalar node)))))

                 (setf (cl-mpm/particle::mp-boundary mp)
                       weathering)
                 (let ((density (/ (cl-mpm/particle::mp-mass mp) (cl-mpm/particle::mp-volume mp))))
                   (setf
                    (cl-mpm/particle::mp-volume mp)
                    (max
                     0d0
                     (-
                      (cl-mpm/particle::mp-volume mp)
                      (abs (*
                            (bc-water-damage-damage-rate bc)
                            weathering dt)))))
                   (setf (cl-mpm/particle::mp-mass mp) (* density (cl-mpm/particle::mp-volume mp))))))
      (cl-mpm::remove-mps-func sim (lambda (mp) (= 0d0 (cl-mpm/particle::mp-mass mp)))))))

(defun setup (&key (refine 1) (mps 2))

  (let* ((density 1d3)
         (mesh-resolution 1d0)
         (domain-size (list 20d0 10d0))
         (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
         (block-size (list 10d0 10d0)))
    (setf *sim* (cl-mpm/setup::make-simple-sim mesh-resolution element-count
                                               :sim-type 'cl-mpm/damage:mpm-sim-damage
                                               ))
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0 0)
      block-size
      (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
      density
      'cl-mpm/particle::particle-elastic-damage
      :E 1d6
      :nu 0.3d0))
    (setf (cl-mpm:sim-dt *sim*)
          (* 0.5d0 (cl-mpm/setup:estimate-elastic-dt *sim*)))
    (setf (cl-mpm::sim-enable-damage *sim*) t)
    (cl-mpm:add-bcs-force-list
     *sim*
     (make-bc-erode *sim* 5d0 0.1d0))
    ))


(defun run ()
  (let* ((dt-scale 0.5d0)
         (target-time 1d0)
         (substeps (ceiling target-time (cl-mpm:sim-dt *sim*))))
    (format t "Substeps ~D~%" substeps)
    (loop for step from 0 below 100
          do
          (format t "Step ~D~%" step)
          (dotimes (i substeps)
            (cl-mpm:update-sim *sim*))
          (cl-mpm/output:save-vtk (merge-pathnames (format nil "./output/sim_~5,'0d.vtk" step)) *sim* )
          (cl-mpm/output:save-vtk-nodes (merge-pathnames (format nil "./output/sim_nodes_~5,'0d.vtk" step)) *sim* )
          (plot-domain)
          (swank.live:update-swank))))
