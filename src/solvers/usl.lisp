(in-package :cl-mpm)
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(defclass mpm-sim-usl (mpm-sim)
  ()
  (:documentation "Explicit simulation with update stress last update"))

(defmethod update-sim ((sim mpm-sim-usl))
  "Update stress last algorithm"
  (declare (cl-mpm::mpm-sim sim))
  (with-slots ((mesh mesh)
               (mps mps)
               (bcs bcs)
               (bcs-force bcs-force)
               (dt dt)
               (mass-filter mass-filter)
               (split allow-mp-split)
               (enable-damage enable-damage)
               (nonlocal-damage nonlocal-damage)
               (remove-damage allow-mp-damage-removal)
               (fbar enable-fbar)
               (bcs-force-list bcs-force-list)
               (vel-algo velocity-algorithm)
               (damping damping-factor)
               (time time)
               )
      sim
    (declare (type double-float mass-filter))
    (progn
      (reset-grid mesh)
      (when (> (length mps) 0)
        ;; Map momentum to grid
        (p2g mesh mps)
        ;;Reset nodes below our mass-filter
        (when (> mass-filter 0d0)
          (filter-grid mesh (sim-mass-filter sim)))
        ;;Turn momentum into velocity
        (update-node-kinematics sim)
        (p2g-force mesh mps)
        (loop for bcs-f in bcs-force-list
              do (apply-bcs mesh bcs-f dt))
        ;; (apply-bcs mesh bcs-force dt)
        ;;Update our nodes after force mapping
        (update-node-forces sim)
        ;;Apply velocity bcs
        (apply-bcs mesh bcs dt)
        (reset-node-displacement sim)
        (update-nodes sim)
        ;;Grid to particle mapping
        (g2p mesh mps dt damping vel-algo)
        ;;2nd round of momentum mapping
        (reset-grid-velocity mesh)
        (p2g mesh mps)
        (when (> mass-filter 0d0)
          (filter-grid-velocity mesh (sim-mass-filter sim)))
        (update-node-kinematics sim)
        (reset-node-displacement sim)
        (update-nodes sim)
        (apply-bcs mesh bcs dt)
        ;;Update stress last
        (update-stress mesh mps dt fbar)
        (update-dynamic-stats sim)
        ;; (update-particles sim)
        (new-loadstep sim)

        (when remove-damage
          (remove-material-damaged sim))
        (when split
          (split-mps sim))
        (check-mps sim))
      (incf time dt)
      )))

