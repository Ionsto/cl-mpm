(defpackage :cl-mpm/dynamic-relaxation-mpi
  (:use :cl
        :cl-mpm/dynamic-relaxation)
  (:export)
  )
(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)

(defclass mpm-sim-dr-mpi (cl-mpm/dynamic-relaxation::mpm-sim-dr
                          cl-mpm/mpi::mpm-sim-mpi-nodes)
  ()
  (:default-initargs
   :vel-algo :QUASI-STATIC)
  (:documentation "DR psudo-linear step with update stress last update"))

(defclass mpm-sim-quasi-static-mpi (mpm-sim-dr-mpi cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul)
  ()
  (:default-initargs
   :vel-algo :QUASI-STATIC)
  (:documentation "DR psudo-linear step with update stress last update"))

(defmethod cl-mpm/dynamic-relaxation::map-stiffness :after ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-mpi))
  (cl-mpm/mpi::mpi-sync-mass sim))


(defmethod save-vtks ((sim mpm-sim-dr-mpi) output-dir step)
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (when (= rank 0)
      (format t "Save vtks ~D~%" step))
    (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d_~5,'0d.vtk" rank step)) sim)
    (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d_~5,'0d.vtk" rank step)) sim)
    (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim_p_~5,'0d_~5,'0d.vtk" rank step)) sim )))

(defmethod save-vtks-dr-step ((sim mpm-sim-dr-mpi) output-dir step iter)
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_step_~5,'0d_~5,'0d_~5,'0d.vtk" rank step iter)) sim)
    (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_step_nodes_~5,'0d_~5,'0d_~5,'0d.vtk" rank step iter)) sim)
    (cl-mpm/penalty:save-vtk-penalties (merge-pathnames output-dir (format nil "sim_step_p_~5,'0d_~5,'0d_~5,'0d.vtk" rank step iter)) sim)))

(defun pre-step-mpi (sim)
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (dt cl-mpm::dt)
               (dt-loadstep dt-loadstep)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (nonlocal-damage cl-mpm::nonlocal-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (ghost-factor cl-mpm::ghost-factor)
               (initial-setup initial-setup)
               (enable-aggregate cl-mpm/aggregate::enable-aggregate)
               (damping cl-mpm::damping-factor)
               (vel-algo cl-mpm::velocity-algorithm))
      sim
    (cl-mpm::reset-grid mesh)
    (cl-mpm::reset-node-displacement sim)
    (cl-mpm::p2g mesh mps)
    (cl-mpm/mpi::mpi-sync-momentum sim)
    (when (> mass-filter 0d0)
      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
    (cl-mpm::filter-cells sim)
    (cl-mpm::update-node-kinematics sim)
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (n)
       (when (cl-mpm/mpi::node-in-computational-domain sim n)
         (setf
          (cl-mpm/mesh::node-true-mass n) (cl-mpm/mesh:node-mass n))
         (cl-mpm/utils:vector-copy-into (cl-mpm/mesh::node-velocity n) (cl-mpm/mesh::node-true-velocity n)))))
    (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
    (update-node-fictious-mass sim)
    (cl-mpm::filter-cells sim)

    (cl-mpm::apply-bcs mesh bcs dt)
    (cl-mpm::update-cells sim)
    (when enable-aggregate
      (cl-mpm/aggregate::update-aggregate-elements sim))
    (cl-mpm::apply-bcs mesh bcs dt)
    (midpoint-starter sim)
    (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
    (setf initial-setup t)))


(defmethod cl-mpm::update-sim ((sim mpm-sim-quasi-static-mpi))
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (dt cl-mpm::dt)
               (dt-loadstep dt-loadstep)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (nonlocal-damage cl-mpm::nonlocal-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (ghost-factor cl-mpm::ghost-factor)
               (initial-setup initial-setup)
               (enable-aggregate cl-mpm/aggregate::enable-aggregate)
               (damping cl-mpm::damping-factor)
               (damping-scale cl-mpm/dynamic-relaxation::damping-scale)
               (vel-algo cl-mpm::velocity-algorithm))
      sim
    (declare (double-float damping-scale damping))
    (unless initial-setup
      (pre-step sim)
      (setf (cl-mpm/damage::sim-damage-delocal-counter-max sim) -1)
      (cl-mpm/damage::update-delocalisation-list mesh mps))
    (cl-mpm/penalty::reset-penalty sim)
    (setf dt 1d0)
    (cl-mpm::update-nodes sim)
    (cl-mpm::update-cells sim)
    (cl-mpm::reset-nodes-force sim)

    (cl-mpm::iterate-over-nodes
     (cl-mpm::sim-mesh sim)
     (lambda (n)
       (when (and (cl-mpm/mesh:node-active n)
                  (not (cl-mpm/mpi::node-in-computational-domain sim n)))
         (cl-mpm/fastmaths:fast-zero (cl-mpm/mesh::node-displacment n))
         (cl-mpm/fastmaths:fast-zero (cl-mpm/mesh::node-velocity n)))))

    (cl-mpm/mpi::mpi-sync-displacement sim)

    (cl-mpm::update-stress mesh mps dt-loadstep fbar)
    ;; (cl-mpm/damage::calculate-damage sim dt-loadstep)
    (cl-mpm::p2g-force-fs sim)
    (cl-mpm::apply-bcs mesh bcs-force dt)

    (loop for bcs-f in bcs-force-list
          do (cl-mpm::apply-bcs mesh bcs-f dt))

    ;; (cl-mpm::iterate-over-nodes
    ;;  (cl-mpm::sim-mesh sim)
    ;;  (lambda (n)
    ;;    (setf (cl-mpm/mesh::node-mass n) 0d0)))

    ;; (update-node-fictious-mass sim)
    ;; ;; (cl-mpm/mpi::mpi-sync-mass sim)

    (cl-mpm/mpi::mpi-sync-force sim)
    ;; ;; ;;Update our nodes after force mapping
    (cl-mpm::update-node-forces sim)
    (cl-mpm::apply-bcs mesh bcs dt)
    (cl-mpm::update-dynamic-stats sim)
    (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
    )
  )
(defmethod cl-mpm::finalise-loadstep :after ((sim mpm-sim-dr-mpi))
  (cl-mpm/mpi::exchange-mps sim 0d0)
  (cl-mpm/mpi::set-mp-mpi-index sim)
  (cl-mpm/mpi::clear-ghost-mps sim))

(defmethod pre-step ((sim mpm-sim-dr-mpi))
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (dt cl-mpm::dt)
               (dt-loadstep dt-loadstep)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (nonlocal-damage cl-mpm::nonlocal-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (ghost-factor cl-mpm::ghost-factor)
               (initial-setup initial-setup)
               (enable-aggregate cl-mpm/aggregate::enable-aggregate)
               (damping cl-mpm::damping-factor)
               (vel-algo cl-mpm::velocity-algorithm))
      sim
    (cl-mpm::reset-grid mesh)
    (cl-mpm::reset-node-displacement sim)
    (cl-mpm::p2g mesh mps)
    (when (> mass-filter 0d0)
      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
    (cl-mpm/mpi::mpi-sync-momentum sim)

    (cl-mpm::filter-cells sim)
    (update-node-fictious-mass sim)
    ;; (cl-mpm/mpi::mpi-sync-momentum sim)

    (cl-mpm::filter-cells sim)
    (cl-mpm::apply-bcs mesh bcs dt)
    (cl-mpm::update-cells sim)
    ;; (when enable-aggregate
    ;;   (cl-mpm/aggregate::update-aggregate-elements sim))
    (cl-mpm::apply-bcs mesh bcs dt)
    (midpoint-starter sim)
    (cl-mpm/mpi::mpi-sync-force sim)
    (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
    (setf initial-setup t)))

(defmethod dr-estimate-damping ((sim mpm-sim-dr-mpi))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (dt cl-mpm:sim-dt))
      sim
    (let ((num 0d0)
          (denom 0d0))
      (setf
       num
       (cl-mpm::reduce-over-nodes
        mesh
        (lambda (node)
          (if (and (cl-mpm/mesh:node-active node)
                   (not (cl-mpm/mesh::node-agg node))
                   (cl-mpm/mpi::node-in-computational-domain sim node))
              (cl-mpm/fastmaths:dot
               (cl-mpm/mesh:node-velocity node)
               (cl-mpm/fastmaths:fast-.-
                (cl-mpm/mesh::node-residual-prev node)
                (cl-mpm/mesh::node-residual node)))
              0d0))
        #'+))
      (setf
       denom
       (* dt
          (cl-mpm::reduce-over-nodes
           mesh
           (lambda (node)
             (if (and (cl-mpm/mesh:node-active node)
                      (not (cl-mpm/mesh::node-agg node))
                      (cl-mpm/mpi::node-in-computational-domain sim node))
                 (* (cl-mpm/mesh:node-mass node)
                    (cl-mpm/fastmaths::dot
                     (cl-mpm/mesh::node-velocity node)
                     (cl-mpm/mesh::node-velocity node)))
                 0d0))
           #'+)))
      (setf num (cl-mpm/mpi::mpi-sum num)
            denom (cl-mpm/mpi::mpi-sum denom))
      (if (> num 0d0)
          (if (= denom 0d0)
              0d0
              (* (sqrt 2) (sqrt (/ num denom))))
          0d0))))

(defun combi-stats-mpi (sim)
  (destructuring-bind (mass
                       energy
                       oobf-num
                       oobf-denom
                       power)
      (cl-mpm::reduce-over-nodes
       (cl-mpm:sim-mesh sim)
       (lambda (node)
         (if (and (cl-mpm/mesh:node-active node)
                  (cl-mpm/mpi::node-in-computational-domain sim node))
             (with-accessors ((active cl-mpm/mesh::node-active)
                              (f-ext cl-mpm/mesh::node-external-force)
                              (res cl-mpm/mesh::node-residual)
                              (node-oobf cl-mpm/mesh::node-oobf)
                              (mass cl-mpm/mesh::node-mass)
                              (volume cl-mpm/mesh::node-volume)
                              (volume-t cl-mpm/mesh::node-volume-true)
                              (vel cl-mpm/mesh::node-velocity)
                              (disp cl-mpm/mesh::node-displacment)
                              )
                 node
               (declare (double-float mass))
               (let (;(mass 1d0)
                     ;; (scale-factor (expt mass 1))
                     (scale-factor 1d0)
                     ;; (scale-factor (/ volume volume-t))
                     )
                 (list
                  scale-factor
                  (* scale-factor (* 0.5d0 mass (cl-mpm/fastmaths::mag-squared vel)))
                  (* scale-factor (cl-mpm/fastmaths::mag-squared res))
                  (* scale-factor (cl-mpm/fastmaths::mag-squared f-ext))
                  (* scale-factor
                     (cl-mpm/fastmaths:dot
                      disp f-ext)))))
             (list 0d0 0d0 0d0 0d0 0d0)))
       (lambda (a b) (mapcar (lambda (x y) (declare (double-float x y)) (+ x y)) a b)))
    (declare (double-float mass energy oobf-num oobf-denom power))
    (let ((oobf 0d0)
          (oobf-num (cl-mpm/mpi::mpi-sum oobf-num))
          (oobf-denom (cl-mpm/mpi::mpi-sum oobf-denom))
          (mass (cl-mpm/mpi::mpi-sum mass))
          (power (cl-mpm/mpi::mpi-sum power))
          (energy (cl-mpm/mpi::mpi-sum energy)))
      (if (> oobf-denom 0d0)
          (setf oobf (sqrt (/ oobf-num oobf-denom)))
          (setf oobf (if (> oobf-num 0d0) sb-ext:double-float-positive-infinity 0d0)))
      (if (> mass 0d0)
          (values energy oobf power)
          (values 0d0 0d0 0d0)))))

(defmethod cl-mpm::update-dynamic-stats ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-mpi))
  (with-accessors ((stats-energy cl-mpm::sim-stats-energy)
                   (stats-oobf cl-mpm::sim-stats-oobf)
                   (stats-power cl-mpm::sim-stats-power)
                   (stats-work cl-mpm::sim-stats-work))
      sim
    (multiple-value-bind (e o p) (combi-stats-mpi sim)
      (setf stats-energy e
            stats-oobf o
            stats-power p)
      (incf stats-work p))))
