(defpackage :cl-mpm/examples/erode
  (:use :cl
        :cl-mpm/example))
(in-package :cl-mpm/particle)


(defclass particle-erosion (particle-elastic-damage)
  ((eroded-volume
    :initform 0d0
    :accessor mp-eroded-volume)))
(in-package :cl-mpm/examples/erode)

(defmethod cl-mpm/output::save-vtk (filename (sim cl-mpm::mpm-sim-usf))
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)) sim
    (with-open-file (fs filename :direction :output :if-exists :supersede)
      (format fs "# vtk DataFile Version 2.0~%")
      (format fs "Lisp generated vtk file, SJVS~%")
      (format fs "ASCII~%")
      (format fs "DATASET UNSTRUCTURED_GRID~%")
      (format fs "POINTS ~d double~%" (length mps))
      (loop for mp across mps
            do (format fs "~E ~E ~E ~%"
                       (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 1 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 2 0) 'single-float)
                       ))
      (format fs "~%")
      (let ((id 1)
            (nd (cl-mpm/mesh:mesh-nd mesh)))
        (declare (special id))
        (format fs "POINT_DATA ~d~%" (length mps))
        (cl-mpm/output::save-parameter "mass" (cl-mpm/particle:mp-mass mp))

        (cl-mpm/output::save-parameter "density" (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)))
        (cl-mpm/output::save-parameter "volume" (cl-mpm/particle::mp-volume mp))

        (cl-mpm/output::save-parameter "boundary" (cl-mpm/particle::mp-boundary mp))
        (cl-mpm/output::save-parameter "erode" (cl-mpm/particle::mp-eroded-volume mp))
        (cl-mpm/output::save-parameter "size_x" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 0 0))
        (cl-mpm/output::save-parameter "size_y" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0))
        (when (= 3 nd)
          (cl-mpm/output::save-parameter "size_z" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 2 0)))

        (cl-mpm/output::save-parameter "damage"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage)
                                           (cl-mpm/particle:mp-damage mp)
                                           0d0))

        )
      )))

(defun plot-domain ()
  (when *sim*
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :colour-func #'cl-mpm/particle::mp-eroded-volume
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
                     :damage-volume nil
                     :clip-func (lambda (pos)
                                  t
                                  ;; (and (> (cl-mpm/utils:varef pos 1) datum))
                                  )
                     :sim sim))))
(defmethod cl-mpm/bc::apply-bc ((bc bc-erode) node mesh dt)
  (call-next-method)
  (with-accessors ((sim cl-mpm/buoyancy::bc-buoyancy-sim)
                   (rate bc-water-damage-damage-rate))
      bc
    (when t;(cl-mpm::sim-enable-damage sim)
      (cl-mpm:iterate-over-mps
       (cl-mpm:sim-mps sim)
       (lambda (mp)
         (with-accessors ((volume cl-mpm/particle::mp-volume)
                          (mass cl-mpm/particle::mp-mass)
                          (erode cl-mpm/particle::mp-eroded-volume)
                          )
             mp
             (let ((weathering 0d0))
               (cl-mpm:iterate-over-neighbours
                mesh
                mp
                (lambda (mesh mp node svp grads fsvp fgrads)
                  (with-accessors ((node-boundary-scalar cl-mpm/mesh::node-boundary-scalar)
                                   (node-volume cl-mpm/mesh::node-volume))
                      node
                    (incf weathering (* svp (/ node-boundary-scalar node-volume))))))

               (setf weathering (min weathering 0d0))
               (setf weathering (* rate (min weathering 0d0) volume))
               ;; (setf weathering (- (sqrt (abs weathering))))
               ;; (setf weathering (* weathering (+ 1d0 (* 8 (cl-mpm/particle:mp-damage mp)))))
               (setf (cl-mpm/particle::mp-boundary mp) weathering)
               (incf erode (- weathering))
               ;; (let ((density (/ mass volume)))
               ;;   (setf
               ;;    mass
               ;;    (max
               ;;     0d0
               ;;     (-
               ;;      mass
               ;;      (abs (*
               ;;            (bc-water-damage-damage-rate bc)
               ;;            weathering dt)))))
               ;;   ;; (setf mass (* density volume))
               ;;   )
               ))))

      (cl-mpm::remove-mps-func sim (lambda (mp) (>= (cl-mpm/particle::mp-eroded-volume mp)
                                                    (cl-mpm/particle::mp-mass mp)))))))

(defun setup (&key (refine 1) (mps 2))
  (let* ((density 1d3)
         (mesh-resolution (/ 1d0 refine))
         (domain-size (list 20d0 10d0))
         (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
         (block-size (list 10d0 10d0)))
    (setf *sim* (cl-mpm/setup::make-simple-sim mesh-resolution element-count
                                               :sim-type 'cl-mpm::mpm-sim-usf
                                               ;:sim-type 'cl-mpm/damage:mpm-sim-damage
                                               ))
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0 0)
      block-size
      (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
      density
      'cl-mpm/particle::particle-erosion
      :E 1d6
      :nu 0.3d0))
    (setf (cl-mpm:sim-dt *sim*)
          (* 0.5d0 (cl-mpm/setup:estimate-elastic-dt *sim*)))
    ;; (setf (cl-mpm::sim-enable-damage *sim*) t)
    (setf *run-sim* t)
    (cl-mpm:add-bcs-force-list
     *sim*
     (make-bc-erode *sim* 100d0 10d0))
    ))


(defun run (&key (output-dir "./output/"))
  (uiop:ensure-all-directories-exist (list (uiop:merge-pathnames* output-dir)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (let* ((dt-scale 1d0)
         (target-time 1d0)
         (substeps (ceiling target-time (cl-mpm:sim-dt *sim*))))
    (format t "Substeps ~D~%" substeps)
    (loop for step from 0 below 25
          while *run-sim*
          do
          (format t "Step ~D~%" step)
          (time 
           (dotimes (i substeps)
             (cl-mpm:update-sim *sim*)))
          (cl-mpm/output:save-vtk (uiop:merge-pathnames* output-dir (format nil "sim_~5,'0d.vtk" step)) *sim* )
          (cl-mpm/output:save-vtk-nodes (uiop:merge-pathnames* output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) *sim* )
          (plot-domain)
          (swank.live:update-swank))))

(defun conv-test ()
  (dolist (refine (list 1 2))
    (dolist (mps (list 2 4))
      (setup :refine refine :mps mps)
      (run :output-dir (format nil "./output-~D-~D/" refine mps)))))
