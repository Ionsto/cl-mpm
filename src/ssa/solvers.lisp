(in-package :cl-mpm/ssa)


(defmethod cl-mpm::update-sim ((sim mpm-sim-ssa))
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (mps cl-mpm::sim-mps)
                   (bcs cl-mpm::sim-bcs)
                   (bcs-force cl-mpm::sim-bcs-force)
                   (bcs-force-list cl-mpm::sim-bcs-force-list)
                   (dt cl-mpm::sim-dt)
                   (mass-filter cl-mpm::sim-mass-filter)
                   (split cl-mpm::sim-allow-mp-split)
                   (remove-damage cl-mpm::sim-allow-mp-damage-removal)
                   (fbar cl-mpm::sim-enable-fbar)
                   (time cl-mpm::sim-time)
                   (vel-algo cl-mpm::sim-velocity-algorithm)
                   (damping cl-mpm::sim-damping-factor)
                   )
      sim
    (declare (double-float mass-filter dt time))
    (cl-mpm::reset-grid mesh)
    (when (> (length mps) 0)
      (cl-mpm::p2g mesh mps)

      (when (> mass-filter 0d0)
        (cl-mpm::filter-grid mesh mass-filter))

      (cl-mpm::filter-cells sim)
      (cl-mpm::update-node-kinematics sim)
      (cl-mpm::apply-bcs mesh bcs dt)
      (cl-mpm::update-nodes sim)
      (cl-mpm::update-cells sim)
      (cl-mpm::update-stress mesh mps dt fbar)
      ;; Map forces onto nodes
      (p2g-force mesh mps)
      (when bcs-force-list
        (loop for bcs-f in bcs-force-list
              do (cl-mpm::apply-bcs mesh bcs-f dt)))
      (cl-mpm::update-node-forces sim)

      (cl-mpm::reset-node-displacement sim)
      (cl-mpm::update-nodes sim)

      (cl-mpm::apply-bcs mesh bcs dt)
      (cl-mpm::update-dynamic-stats sim)
      ;; Also updates mps inline
      (cl-mpm::g2p mesh mps dt damping vel-algo)
      (cl-mpm::new-loadstep sim)
      ))
  )
