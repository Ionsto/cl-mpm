(in-package :cl-mpm/mpi)

(defmethod sim-add-mp ((sim mpm-sim-mpi) mp)
  (with-accessors ((uid-counter cl-mpm::sim-unique-index-counter)
                   (lock cl-mpm::sim-unique-index-lock))
      sim
    (let* ((size (cl-mpi:mpi-comm-size))
           (rank (cl-mpi:mpi-comm-rank))
           (shift (ceiling (log size 2))))
      (sb-thread:with-mutex (lock)
        (setf (cl-mpm/particle::mp-unique-index mp) (+ (ash uid-counter shift) rank))
        (incf uid-counter)))
    (call-next-method)))


(defgeneric update-min-domain-size (sim))
(defmethod update-min-domain-size ((sim mpm-sim-mpi))
  (setf (mpm-sim-mpi-min-size sim)
        (* (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) (mpm-sim-mpi-halo-depth sim))))
(defmethod update-min-domain-size ((sim mpm-sim-mpi-nodes-damage))
  (setf (mpm-sim-mpi-min-size sim)
        (max
         (* (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) (mpm-sim-mpi-halo-depth sim))
         (mpm-sim-mpi-halo-damage-size sim))))
(defmethod (setf cl-mpm::sim-mesh) :after (value (sim mpm-sim-mpi))
  (update-min-domain-size sim))
(defmethod (setf mpm-sim-mpi-halo-depth) :after (value (sim mpm-sim-mpi))
  (update-min-domain-size sim))
(defmethod (setf mpm-sim-mpi-halo-damage-size) :after (value (sim mpm-sim-mpi-nodes-damage))
  (update-min-domain-size sim))


(defmethod cl-mpm::update-sim ((sim mpm-sim-mpi-nodes))
  (with-slots ((mesh cl-mpm::mesh)
               (mps  cl-mpm::mps)
               (bcs  cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (bcs-force-list cl-mpm::bcs-force-list)
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
                    ;; (exchange-mps sim)
                  (when t;(> (length mps) 0)
                    (cl-mpm::reset-grid mesh)
                    (cl-mpm::p2g mesh mps vel-algo)
                    (mpi-sync-momentum sim)
                    (when (> mass-filter 0d0)
                      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
                    (cl-mpm::update-node-kinematics sim)
                    (cl-mpm::apply-bcs mesh bcs dt)
                    (cl-mpm::update-stress mesh mps dt fbar)
                    (cl-mpm::p2g-force sim)
                    (loop for bcs-f in bcs-force-list
                          do (cl-mpm::apply-bcs mesh bcs-f dt))
                    (mpi-sync-force sim)
                    (cl-mpm::update-node-forces sim)
                    (cl-mpm::apply-bcs mesh bcs dt)
                    (cl-mpm::g2p mesh mps dt vel-algo)
                    (cl-mpm::update-particles sim)
                    (when remove-damage
                      (cl-mpm::remove-material-damaged sim))
                    (when split
                      (cl-mpm::split-mps sim))
                    (cl-mpm::check-mps sim)
                    (set-mp-mpi-index sim)
                    )
                  (exchange-mps sim 0d0)
                  (set-mp-mpi-index sim)
                  (clear-ghost-mps sim))))



(defmethod cl-mpm::update-sim ((sim mpm-sim-mpi-nodes-damage))
  (with-slots ((mesh cl-mpm::mesh)
               (mps  cl-mpm::mps)
               (bcs  cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (bcs-force-list cl-mpm::bcs-force-list)
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
                  (cl-mpm::reset-grid mesh)
                  (cl-mpm::p2g mesh mps vel-algo)
                  (mpi-sync-momentum sim)
                  (when (> mass-filter 0d0)
                    (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
                  (cl-mpm::update-node-kinematics sim)
                  (cl-mpm::apply-bcs mesh bcs dt)
                  (cl-mpm::update-stress mesh mps dt fbar)
                  (cl-mpm/damage::calculate-damage sim)
                  (cl-mpm::p2g-force sim)
                  (loop for bcs-f in bcs-force-list
                        do (cl-mpm::apply-bcs mesh bcs-f dt))
                  (mpi-sync-force sim)
                  (cl-mpm::update-node-forces sim)
                  (cl-mpm::apply-bcs mesh bcs dt)
                  (cl-mpm::update-dynamic-stats sim)
                  (cl-mpm::g2p mesh mps dt vel-algo)
                  (cl-mpm::update-particles sim)
                  (when remove-damage
                    (cl-mpm::remove-material-damaged sim))
                  (when split
                    (cl-mpm::split-mps sim))
                  (cl-mpm::check-mps sim)
                  ;; (cl-mpm::check-single-mps sim)
                  ;; (set-mp-mpi-index sim)
                  (exchange-mps sim 0d0)
                  (set-mp-mpi-index sim)
                  (clear-ghost-mps sim))))

(defmethod cl-mpm::update-sim ((sim mpm-sim-usl-mpi-nodes-damage))
  (with-slots ((mesh cl-mpm::mesh)
               (mps  cl-mpm::mps)
               (bcs  cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (bcs-force-list cl-mpm::bcs-force-list)
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
                  (cl-mpm::reset-grid mesh)
                  (cl-mpm::p2g mesh mps vel-algo)
                  (mpi-sync-momentum sim)
                  (when (> mass-filter 0d0)
                    (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
                  (cl-mpm::update-node-kinematics sim)
                  (cl-mpm::apply-bcs mesh bcs dt)
                  (cl-mpm::p2g-force sim)
                  (loop for bcs-f in bcs-force-list
                        do
                           (cl-mpm::apply-bcs mesh bcs-f dt))
                  (mpi-sync-force sim)
                  (cl-mpm::update-node-forces sim)
                  (cl-mpm::apply-bcs mesh bcs dt)
                  (cl-mpm::g2p mesh mps dt vel-algo)


                  ;;2nd round of mapping for USL
                  (cl-mpm::reset-grid-velocity mesh)
                  (cl-mpm::p2g mesh mps vel-algo)
                  (mpi-sync-momentum sim)
                  (when (> mass-filter 0d0)
                    (cl-mpm::filter-grid-velocity mesh (cl-mpm::sim-mass-filter sim)))
                  (cl-mpm::update-node-kinematics sim)
                  (cl-mpm::apply-bcs mesh bcs dt)

                  (cl-mpm::update-stress mesh mps dt fbar)
                  (cl-mpm/damage::calculate-damage sim)
                  (cl-mpm::update-dynamic-stats sim)
                  (cl-mpm::update-particles sim)

                  (when remove-damage
                    (cl-mpm::remove-material-damaged sim))
                  (when split
                    (cl-mpm::split-mps sim))

                  ;;Must be done here
                  (set-mp-mpi-index sim)
                  (exchange-mps sim 0d0)
                  (set-mp-mpi-index sim)
                  (clear-ghost-mps sim)
                  (cl-mpm::check-mps sim)
                  (cl-mpm::check-single-mps sim)

                    )))

(defmethod cl-mpm::update-sim ((sim mpm-sim-mpi))
  (with-slots ((mesh cl-mpm::mesh)
               (mps  cl-mpm::mps)
               (bcs  cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (bcs-force-list cl-mpm::bcs-force-list)
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
                    (exchange-mps sim)
                  (when (> (length mps) 0)
                    (cl-mpm::reset-grid mesh)
                    (cl-mpm::p2g mesh mps vel-algo)
                    (when (> mass-filter 0d0)
                      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
                    (cl-mpm::update-node-kinematics sim)
                    (cl-mpm::apply-bcs mesh bcs dt)
                    ;; (cl-mpm::update-stress mesh mps dt)
                    (cl-mpm::update-stress mesh mps dt fbar)
                                        ;(exchange-mps sim)
                    (cl-mpm/damage::calculate-damage sim)
                                        ;(exchange-mps sim)
                    (cl-mpm::p2g-force sim)
                    (loop for bcs-f in bcs-force-list
                          do
                             (cl-mpm::apply-bcs mesh bcs-f dt))
                    (cl-mpm::update-node-forces sim)
                    (cl-mpm::apply-bcs mesh bcs dt)
                                        ;Also updates mps inline
                    (cl-mpm::g2p mesh mps dt vel-algo)
                    ;;MPI reduce new velocities
                                        ;(exchange-mps sim)
                    ;;Get new MPS

                    (cl-mpm::update-particles sim)
                    (when remove-damage
                      (cl-mpm::remove-material-damaged sim))
                    (when split
                      (cl-mpm::split-mps sim))
                    (cl-mpm::check-mps sim)
                    (set-mp-mpi-index sim))
                    ;; (clear-ghost-mps sim)
                    ;;Update mp list between processors
                    )))

(defmethod cl-mpm::calculate-min-dt-mps ((sim cl-mpm/mpi::mpm-sim-mpi))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mass-scale cl-mpm::sim-mass-scale))
      sim
    (let ((inner-factor 
            ;;MPI not a fan of most-positive-double-float
            1d50
            ;most-positive-double-float
                        ))
      (iterate-over-nodes-serial
       mesh
       (lambda (node)
         (with-accessors ((node-active  cl-mpm/mesh:node-active)
                          (node-pos cl-mpm/mesh::node-position)
                          (pmod cl-mpm/mesh::node-pwave)
                          (mass cl-mpm/mesh::node-mass)
                          (svp-sum cl-mpm/mesh::node-svp-sum)
                          (vol cl-mpm/mesh::node-volume)
                          ) node
           (when (and node-active
                      ;(in-computational-domain sim node-pos)
                      (> vol 0d0)
                      (> pmod 0d0)
                      (> svp-sum 0d0))
             (let ((nf (/ mass (* vol (/ pmod svp-sum)))))
                 (when (< nf inner-factor)
                   (setf inner-factor nf)))))))
      (let ((rank (cl-mpi:mpi-comm-rank))
            (size (cl-mpi:mpi-comm-size)))
        (static-vectors:with-static-vector (source 1 :element-type 'double-float :initial-element inner-factor)
          (static-vectors:with-static-vector (dest 1 :element-type 'double-float :initial-element 0d0)
            (cl-mpi:mpi-allreduce source dest cl-mpi:+mpi-min+)
            (setf inner-factor (aref dest 0))
            (if (< inner-factor most-positive-double-float)
                (progn
                  ;; (format t "Rank ~D: dt - ~F~%" rank (* (sqrt mass-scale) (sqrt inner-factor) (cl-mpm/mesh:mesh-resolution mesh)))
                  ;; (format t "global : dt - ~F~%" (* (sqrt mass-scale) (sqrt inner-factor) (cl-mpm/mesh:mesh-resolution mesh)))
                  (* (sqrt mass-scale) (sqrt inner-factor) (cl-mpm/mesh:mesh-resolution mesh)))
                (cl-mpm::sim-dt sim))))))))
