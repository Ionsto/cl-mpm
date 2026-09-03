(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)


(defclass mpm-sim-dr-dynamic-mpi (mpm-sim-dr-dynamic
                                           mpm-sim-dr-mpi
                                           cl-mpm/mpi::mpm-sim-mpi-damage)
  ()
  (:default-initargs
   :vel-algo :TBLEND)
  (:documentation "MPI enabled implicit dynamic stepper"))

(defmethod pre-step ((sim mpm-sim-dr-dynamic-mpi))
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
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
               (vel-algo cl-mpm::velocity-algorithm)
               (solve-count cl-mpm/dynamic-relaxation::solve-count)
               )
      sim
    (setf (cl-mpm/dynamic-relaxation::sim-solve-count sim) 0)
    (cl-mpm::reset-grid mesh :reset-displacement t)
    (cl-mpm::p2g mesh mps vel-algo)
    (setf (cl-mpm::sim-dt sim) 1d0)
    (when (> mass-filter 0d0)
      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
    (cl-mpm::apply-essential-bcs sim)
    (cl-mpm::filter-cells sim)
    (cl-mpm::update-node-kinematics sim)
    (cl-mpm::apply-essential-bcs sim)
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (n)
       (setf (cl-mpm/mesh::node-true-mass n) (cl-mpm/mesh:node-mass n))
       (cl-mpm/utils:vector-copy-into (cl-mpm/mesh::node-velocity n) (cl-mpm/mesh::node-true-velocity n))))
    (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
    (cl-mpm::reset-node-displacement sim)
    (setf (cl-mpm::sim-damping-factor sim) 0d0)
    (midpoint-starter sim)
    (setf initial-setup t)))

(defmethod cl-mpm::update-sim ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic-mpi))
  "Update stress last algorithm"
  (declare (cl-mpm::mpm-sim sim))
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (dt-loadstep dt-loadstep)
               (damping-scale cl-mpm/dynamic-relaxation::damping-scale)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (ghost-factor cl-mpm::ghost-factor)
               (initial-setup initial-setup)
               (enable-aggregate cl-mpm/aggregate::enable-aggregate)
               (damping cl-mpm::damping-factor)
               (mass-update-iter cl-mpm/dynamic-relaxation::mass-update-count)
               (solve-count cl-mpm/dynamic-relaxation::solve-count)
               (vel-algo cl-mpm::velocity-algorithm))
      sim
    (unless initial-setup
      (pre-step sim))
    (setf (cl-mpm:sim-dt sim) 1d0)
    (cl-mpm/penalty::reset-penalty sim)

    (cl-mpm::reset-nodes-force sim)
    (cl-mpm::apply-essential-bcs sim)
    (cl-mpm::apply-force-bcs sim dt-loadstep)

    (cl-mpm/mpi::with-mpi-errors
        (cl-mpm::update-stress mesh mps dt-loadstep fbar))

    (cl-mpm::p2g-force-fs sim)
    (cl-mpm/mpi::with-mpi-errors
        (cl-mpm/mpi::mpi-sync-force sim))
    (cl-mpm/mpi::with-mpi-errors
        (update-node-fictious-mass sim))
    (cl-mpm::update-node-forces sim)
    (cl-mpm::apply-essential-bcs sim)
    (cl-mpm::update-nodes sim)
    (cl-mpm/mpi::with-mpi-errors
        (cl-mpm/mpi::mpi-sync-displacement sim))
    (cl-mpm::update-filtered-cells sim)
    (cl-mpm::update-dynamic-stats sim)
    (incf solve-count)))

(defmethod update-node-fictious-mass ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic-mpi))
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (dt cl-mpm::sim-dt)
                   (bcs-force-list cl-mpm::sim-bcs-force-list))
      sim

    (cl-mpm::iterate-over-nodes
     (cl-mpm:sim-mesh sim)
     (lambda (n)
       (when t
         (setf (cl-mpm/mesh::node-mass n) 0d0))))
    (map-stiffness sim)
    (loop for bcs-f in bcs-force-list
          do (loop for bc across bcs-f
                   do (cl-mpm/bc::assemble-bc-stiffness sim bc)))
    (cl-mpm/mpi::mpi-sync-mass sim)
    (cl-mpm/aggregate::update-mass-matrix sim)
    (setf dt 1d0)))
