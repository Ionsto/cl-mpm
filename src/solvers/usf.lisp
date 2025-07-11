(in-package :cl-mpm)
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(defclass mpm-sim-usf (mpm-sim)
  ()
  (:documentation "Explicit simulation with update stress first update"))

(defmethod update-sim ((sim mpm-sim-usf))
  "Update stress first algorithm"
  (declare (cl-mpm::mpm-sim-usf sim))
  (with-slots ((mesh mesh)
               (mps mps)
               (bcs bcs)
               (bcs-force bcs-force)
               (bcs-force-list bcs-force-list)
               (dt dt)
               (mass-filter mass-filter)
               (split allow-mp-split)
               (enable-damage enable-damage)
               (nonlocal-damage nonlocal-damage)
               (remove-damage allow-mp-damage-removal)
               (fbar enable-fbar)
               (time time)
               (vel-algo velocity-algorithm)
               (damping damping-factor)
               (ghost-factor cl-mpm::ghost-factor)
               )
                sim
    (declare (double-float mass-filter dt time))
    (progn
      (reset-grid mesh)
      (when (> (length mps) 0)
        (p2g mesh mps)
        (when (> mass-filter 0d0)
          (filter-grid mesh (sim-mass-filter sim)))
        (filter-cells sim)
        (update-node-kinematics sim)
        (apply-bcs mesh bcs dt)
        ;;Trial update displacements
        (update-nodes sim)
        (update-cells sim)

        (update-stress mesh mps dt fbar)
        ;; Map forces onto nodes
        (p2g-force mesh mps)
        (when bcs-force-list
          (loop for bcs-f in bcs-force-list
                do (apply-bcs mesh bcs-f dt)))
        ;; (cl-mpm/ghost::apply-ghost sim ghost-factor)
        (update-node-forces sim)
        ;; (cl-mpm/ghost::apply-half-step-ghost sim)
        ;; Reapply velocity BCs
        ;; (apply-bcs mesh bcs dt)
        ;;Compute true displacements
        (reset-node-displacement sim)
        (update-nodes sim)

        (apply-bcs mesh bcs dt)
        (update-dynamic-stats sim)
        ;; Also updates mps inline
        (g2p mesh mps dt damping vel-algo)
        (new-loadstep sim)
        ;; (check-single-mps sim)
        )
      (incf time dt))))
