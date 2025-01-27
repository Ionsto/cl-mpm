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
     (cl-mpm/erosion::make-bc-erode *sim* :rate 100d0))))


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


(defun p2g-erosion-mp (mesh mp)
  "P2G transfer for one MP"
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (with-accessors ((mp-vel  cl-mpm/particle:mp-velocity)
                   (mp-mass cl-mpm/particle:mp-mass)
                   (mp-volume cl-mpm/particle:mp-volume)
                   (mp-pmod cl-mpm/particle::mp-p-modulus)
                   (mp-damage cl-mpm/particle::mp-damage))
      mp
    (let ((mp-mass mp-mass)
          (mp-vel mp-vel)
          (mp-volume mp-volume)
          (mp-pmod mp-pmod)
          (mp-damage mp-damage))
      (declare (type double-float mp-mass mp-volume))
      (cl-mpm:iterate-over-neighbours
       mesh mp
       (lambda (mesh mp node svp grads fsvp fgrads)
         (declare
          (cl-mpm/particle:particle mp)
          (cl-mpm/mesh::node node)
          (double-float svp)
          )
         (with-accessors ((node-vel   cl-mpm/mesh:node-velocity)
                          (node-active  cl-mpm/mesh::node-active)
                          (node-mass  cl-mpm/mesh:node-mass)
                          (node-volume  cl-mpm/mesh::node-volume)
                          (node-volume-true  cl-mpm/mesh::node-volume-true)
                          (node-svp-sum  cl-mpm/mesh::node-svp-sum)
                          (node-force cl-mpm/mesh:node-force)
                          (node-p-wave cl-mpm/mesh::node-pwave)
                          (node-damage cl-mpm/mesh::node-damage)
                          (node- cl-mpm/mesh::node-damage)
                          (node-lock  cl-mpm/mesh:node-lock)) node
           (declare (type double-float node-mass node-volume mp-volume mp-pmod mp-damage node-svp-sum svp node-p-wave node-damage)
                    (type sb-thread:mutex node-lock))
           (sb-thread:with-mutex (node-lock)
             (setf node-active t)
             (incf node-mass
                   (* mp-mass svp))
             (incf node-volume
                   (* mp-volume svp))
             (incf node-p-wave
                   (* mp-pmod svp))
             (incf node-damage
                   (* mp-damage svp))
             (incf node-svp-sum svp)
             (cl-mpm/fastmaths::fast-fmacc node-vel mp-vel (* mp-mass svp))
             )
           ;;Ideally we include these generic functions for special mapping operations, however they are slow
           ;; (special-p2g mp node svp dsvp)
           )))
      ))
  (values))

(declaim (notinline p2g))
(defun p2g-erosion (mesh mps)
  "Map particle momentum to the grid"
  (declare (type (vector cl-mpm/particle:particle) mps) (cl-mpm/mesh::mesh mesh))
  (cl-mpm:iterate-over-mps
   mps
   (lambda (mp)
     (p2g-erosion-mp mesh mp))))
