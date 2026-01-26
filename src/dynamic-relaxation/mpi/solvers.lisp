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

(defclass mpm-sim-quasi-static-mpi (cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul
                                    mpm-sim-dr-mpi)
  ()
  (:default-initargs
   :vel-algo :QUASI-STATIC)
  (:documentation "DR psudo-linear step with update stress last update"))

(defmethod cl-mpm/dynamic-relaxation::map-stiffness :after ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-mpi))
  (cl-mpm/mpi::mpi-sync-mass sim))

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
    (cl-mpm/damage::calculate-damage sim dt-loadstep)
    (cl-mpm::p2g-force-fs sim)
    (cl-mpm::apply-bcs mesh bcs-force dt)

    (loop for bcs-f in bcs-force-list
          do (cl-mpm::apply-bcs mesh bcs-f dt))

    (cl-mpm::iterate-over-nodes
     (cl-mpm::sim-mesh sim)
     (lambda (n)
       (setf (cl-mpm/mesh::node-mass n) 0d0)))

    (update-node-fictious-mass sim)
    ;; (cl-mpm/mpi::mpi-sync-mass sim)

    ;; (when ghost-factor
    ;;   (cl-mpm/ghost::apply-ghost sim ghost-factor)
    ;;   (cl-mpm::apply-bcs mesh bcs dt))

    (cl-mpm/mpi::mpi-sync-force sim)

    ;; (setf damping-scale 0d0)
    ;; (setf damping 0d0)
    (setf damping (* damping-scale (cl-mpm/dynamic-relaxation::dr-estimate-damping sim)))
    ;; (format t "Rank ~D - ~E~%" (cl-mpi:mpi-comm-rank) damping)

    ;; ;;Update our nodes after force mapping
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
