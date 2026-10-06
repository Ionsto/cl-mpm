(in-package :cl-mpm/mpi)

(defmethod cl-mpm::update-sim ((sim mpm-sim-mpi-usf))
  (with-mpi-errors
      (with-slots ((mesh cl-mpm::mesh)
                   (mps  cl-mpm::mps)
                   (bcs  cl-mpm::bcs)
                   (bcs-force cl-mpm::bcs-force)
                   (dt cl-mpm::dt)
                   (mass-filter cl-mpm::mass-filter)
                   (split cl-mpm::allow-mp-split)
                   (enable-damage cl-mpm::enable-damage)
                   (nonlocal-damage cl-mpm::nonlocal-damage)
                   (remove-damage cl-mpm::allow-mp-damage-removal)
                   (fbar cl-mpm::enable-fbar)
                   (vel-algo cl-mpm::velocity-algorithm)
                   )
          sim
        (declare (type double-float mass-filter))
        (progn

          ;; (set-mp-mpi-index sim)
          ;; (exchange-mps sim 0d0)
          ;; (set-mp-mpi-index sim)
          ;; (clear-ghost-mps sim)
          ;; (exchange-mps sim)
          (with-mpi-errors
              (cl-mpm::reset-grid mesh)
            (cl-mpm::p2g mesh mps vel-algo))
          (mpi-sync-momentum sim)
          (with-mpi-errors
              (when (> mass-filter 0d0)
                (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
            (cl-mpm::apply-essential-bcs sim)
            (cl-mpm::filter-cells sim)
            (cl-mpm::update-node-kinematics sim)
            (cl-mpm::apply-essential-bcs sim)
            (cl-mpm::update-nodes sim)
            (cl-mpm::update-filtered-cells sim))
          (with-mpi-errors
              (cl-mpm::update-stress mesh mps dt fbar)
            (cl-mpm::update-stiffness-mps sim)
            (cl-mpm::p2g-force sim)
            (cl-mpm::apply-force-bcs sim dt))
          (mpi-sync-force sim)
          (with-mpi-errors
              (cl-mpm::update-node-forces sim)
            (cl-mpm::apply-essential-bcs sim)
            (cl-mpm::reset-node-displacement sim)
            (cl-mpm::update-nodes sim)
            (cl-mpm::apply-essential-bcs sim)
            (mpi-sync-displacement sim)
            (cl-mpm::update-dynamic-stats sim)
            (cl-mpm::g2p mesh mps dt 0d0 vel-algo)
            (cl-mpm::new-loadstep sim))
          (set-mp-mpi-index sim)
          (exchange-mps sim 0d0)
          (set-mp-mpi-index sim)
          (clear-ghost-mps sim)))))



(defmethod cl-mpm::update-particles ((sim cl-mpm/mpi::mpm-sim-mpi))
  (with-mpi-errors
    (call-next-method)))
