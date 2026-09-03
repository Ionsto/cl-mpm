(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)


(defun midpoint-starter-mpi (sim)
     (with-slots ((mesh cl-mpm::mesh)
                  (mps cl-mpm::mps)
                  (bcs cl-mpm::bcs)
                  (bcs-force cl-mpm::bcs-force)
                  (dt cl-mpm::dt)
                  (fbar cl-mpm::enable-fbar)
                  (dt-loadstep dt-loadstep)
                  (agg cl-mpm/aggregate::enable-aggregate)
                  (bcs-force-list cl-mpm::bcs-force-list))
         sim
       (progn
         (setf dt 1d0)
         (cl-mpm::update-stress mesh mps dt-loadstep fbar)
         (cl-mpm/damage::calculate-damage sim dt-loadstep)
         (cl-mpm::p2g-force-fs sim)
         (cl-mpm::apply-essential-bcs sim)
         (cl-mpm::apply-force-bcs sim dt-loadstep)
         (cl-mpm/mpi::mpi-sync-force sim)
         (let ((dt-scale (cl-mpm::sim-dt-scale sim)))
           (setf (cl-mpm::sim-dt-scale sim) (* 1d0 dt-scale))
           (update-node-fictious-mass sim)
           ;; (cl-mpm/aggregate::update-node-forces-agg sim (* -0.5d0 dt))
           (update-node-forces-midpoint-starter sim)
           ;; (cl-mpm::reset-node-displacement sim)
           (setf (cl-mpm::sim-dt-scale sim) dt-scale)
           ;; (cl-mpm::update-nodes sim)
           )
         ;; (cl-mpm/aggregate::project-displacement sim)
         (cl-mpm::iterate-over-nodes
          mesh
          (lambda (n)
            (when (cl-mpm/mesh:node-active n)
              (cl-mpm/fastmaths::fast-zero (cl-mpm/mesh::node-force n))
              (cl-mpm/fastmaths::fast-zero (cl-mpm/mesh::node-residual n))
              (cl-mpm/fastmaths::fast-zero (cl-mpm/mesh::node-residual-prev n))
              ;; (cl-mpm/fastmaths::fast-zero (cl-mpm/mesh::node-force n))
              ;; (cl-mpm/utils::vector-copy-into (cl-mpm/mesh::node-force n)
              ;;                                 (cl-mpm/mesh::node-residual n))
              )))
         (cl-mpm::apply-essential-bcs sim))))

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
    (setf (cl-mpm/dynamic-relaxation::sim-solve-count sim) 0)
    (cl-mpm::reset-grid mesh :reset-displacement t)
    (cl-mpm::reset-node-displacement sim)
    (cl-mpm::p2g mesh mps vel-algo)
    (setf (cl-mpm::sim-dt sim) 1d0)
    (cl-mpm/mpi::mpi-sync-momentum sim)
    (when (> mass-filter 0d0)
      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
    (cl-mpm::apply-essential-bcs sim)
    (cl-mpm::filter-cells sim)
    ;; (cl-mpm::apply-essential-bcs sim)
    ;; (cl-mpm::update-node-kinematics sim)
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (n)
       (setf
        (cl-mpm/mesh::node-true-mass n) (cl-mpm/mesh:node-mass n))
       (cl-mpm/fastmaths:fast-zero (cl-mpm/mesh::node-true-velocity n))))
    (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
    (cl-mpm::reset-node-displacement sim)
    ;; (midpoint-starter-mpi sim)
    (setf (cl-mpm::sim-damping-factor sim) 0d0)
    (setf initial-setup t)))

(defmethod update-node-fictious-mass ((sim cl-mpm/dynamic-relaxation::mpm-sim-quasi-static-mpi))
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (dt cl-mpm::sim-dt)
                   (bcs-force-list cl-mpm::sim-bcs-force-list))
      sim

    (cl-mpm::iterate-over-nodes
     (cl-mpm:sim-mesh sim)
     (lambda (n)
       ;; (not (cl-mpm/mpi::node-in-computational-domain sim n))
       (when t
         (setf (cl-mpm/mesh::node-mass n) 0d0))))
    (map-stiffness-quasi-static sim)
    (loop for bcs-f in bcs-force-list
          do (loop for bc across bcs-f
                   do (cl-mpm/bc::assemble-bc-stiffness sim bc)))
    (cl-mpm/mpi::mpi-sync-mass sim)
    (cl-mpm/aggregate::update-mass-matrix sim)
    (setf dt 1d0)))


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
      (pre-step sim))
    (cl-mpm/penalty::reset-penalty sim)
    (setf dt 1d0)
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
    ;;Update our nodes after force mapping
    (update-node-forces-quasi-static sim)
    (cl-mpm::apply-essential-bcs sim)
    (cl-mpm::update-nodes sim)
    (cl-mpm/mpi::with-mpi-errors
        (cl-mpm/mpi::mpi-sync-displacement sim))
    (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)))

(defmethod cl-mpm::update-sim ((sim mpm-sim-damage-quasi-static-mpi))
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
               ;; (solve-count cl-mpm/dynamic-relaxation::solve-count)
               (vel-algo cl-mpm::velocity-algorithm))
      sim
    (declare (double-float damping-scale damping))
    (unless initial-setup
      (pre-step sim)
      (setf (cl-mpm/damage::sim-damage-delocal-counter-max sim) -1)
      (cl-mpm/damage::update-delocalisation-list mesh mps))
    (cl-mpm/penalty::reset-penalty sim)
    (setf dt 1d0)
    (cl-mpm::reset-nodes-force sim)
    (cl-mpm::apply-essential-bcs sim)
    (cl-mpm::apply-force-bcs sim dt-loadstep)

    (cl-mpm/mpi::with-mpi-errors
        (cl-mpm::update-stress mesh mps dt-loadstep fbar))
    (cl-mpm/mpi::with-mpi-errors
        (cl-mpm/damage::calculate-damage sim dt-loadstep))

    (cl-mpm::p2g-force-fs sim)
    (cl-mpm/mpi::with-mpi-errors
        (cl-mpm/mpi::mpi-sync-force sim))

    (cl-mpm/mpi::with-mpi-errors
        (update-node-fictious-mass sim))
    (update-node-forces-quasi-static sim)
    (cl-mpm::apply-essential-bcs sim)
    (cl-mpm::update-nodes sim)
    (cl-mpm/mpi::with-mpi-errors
        (cl-mpm/mpi::mpi-sync-displacement sim))
    (cl-mpm::apply-essential-bcs sim)
    (cl-mpm::update-filtered-cells sim)
    (cl-mpm::g2p mesh mps dt damping :TRIAL)
    ;; (incf solve-count)
    (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
    (cl-mpm::update-dynamic-stats sim)))

(defmethod cl-mpm::finalise-loadstep :after ((sim mpm-sim-dr-mpi))
  (cl-mpm/mpi::set-mp-mpi-index sim)
  (cl-mpm/mpi::exchange-mps sim 0d0)
  (cl-mpm/mpi::set-mp-mpi-index sim)
  (cl-mpm/mpi::clear-ghost-mps sim))

(defmethod cl-mpm/dynamic-relaxation::pre-step ((sim mpm-sim-quasi-static-mpi))
  (pre-step-mpi sim)
  )

(let ((work-pool (cl-mpm/utils::make-object-pool
                  :constructor (lambda ()
                                 (cl-mpm/utils:vector-zeros)))))
  (defmethod dr-estimate-damping ((sim mpm-sim-dr-mpi))
    (with-accessors ((mesh cl-mpm:sim-mesh))
        sim
      (let ((num 0d0)
            (denom 0d0)
            (dt 1d0))
        (declare (double-float denom num))
        (cl-mpm/utils::object-pool-ensure-size work-pool)
        (setf
         num
         (float
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
                   (cl-mpm/mesh::node-residual node)
                   (cl-mpm/utils::object-pool-grab-unsafe work-pool)
                   ))
                 0d0))
           #'+)
          0d0))
        (setf
         denom
         (* dt
            (the double-float
                 (float
                  (cl-mpm::reduce-over-nodes
                   mesh
                   (lambda (node)
                     (the double-float
                          (if (and (cl-mpm/mesh:node-active node)
                                   (not (cl-mpm/mesh::node-agg node))
                                   (cl-mpm/mpi::node-in-computational-domain sim node))
                              (progn
                                (when (not (typep (cl-mpm/mesh:node-mass node) 'double-float))
                                  (break))
                                (*
                                 (the double-float (cl-mpm/mesh:node-mass node))
                                 (the double-float (cl-mpm/fastmaths::mag-squared
                                                    (cl-mpm/mesh::node-velocity node)))))
                              0d0)))
                   #'+)
                  0d0
                  ))))

        (setf num (cl-mpm/mpi::mpi-sum num)
              denom (cl-mpm/mpi::mpi-sum denom))
        (min 1.9d0
             (max 0d0
                  (* (the double-float (cl-mpm/dynamic-relaxation::sim-damping-scale sim))
                     (if (> num 0d0)
                         (if (= denom 0d0)
                             0d0
                             (the double-float
                                  (* (sqrt 2d0)
                                     (the double-float
                                          (sqrt (max 0d0 (the double-float (/ num denom))))))))
                         0d0))))))
    ))


(cl-mpm/utils::with-arb-pool
  (defun update-node-forces-qs-mpi (sim)
    (with-accessors ((mesh cl-mpm:sim-mesh)
                     (mass-scale cl-mpm:sim-mass-scale)
                     (damping cl-mpm:sim-damping-factor)
                     (damping-scale cl-mpm/dynamic-relaxation::sim-damping-scale)
                     (damping-algo cl-mpm::sim-damping-algorithm)
                     (agg-elems cl-mpm/aggregate::sim-agg-elems)
                     (dt cl-mpm:sim-dt)
                     (solve-count sim-solve-count)
                     (damping-update-count sim-damping-update-count)
                     (enable-aggregate cl-mpm/aggregate::sim-enable-aggregate))
        sim
      (declare (fixnum solve-count damping-update-count)
               (double-float dt damping damping-scale))
      ;; (cl-mpm::apply-essential-bcs sim)
      (cl-mpm:iterate-over-nodes
       mesh
       (lambda (node)
         (when (cl-mpm/mesh:node-active node)
           (cl-mpm::calculate-forces-midpoint node 0d0 0d0 mass-scale))))
      (cl-mpm::compute-reaction-force sim)
      ;;For each aggregated element set solve mass matrix and velocity
      (when enable-aggregate
        (let* ()
          (cl-mpm/aggregate::iterate-over-dimensions
           (cl-mpm::mesh-nd mesh)
           (lambda (d)
             (let* ((f (cl-mpm/aggregate::aggregate-vec
                        sim
                        (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-force d) d))
                    (et (cl-mpm/aggregate::sim-global-sparse-et sim))
                    (e (cl-mpm/aggregate::sim-global-sparse-e sim))
                    (sma (cl-mpm/aggregate::sim-global-sparse-ma sim))
                    (bcs-int (aref (cl-mpm/aggregate::sim-global-bcs-int sim) d))
                    (bcs (aref (cl-mpm/aggregate::sim-global-bcs sim) d))
                    (work-vec (grab-new))
                    (work-vec-agg (grab-new))
                    (fg (grab-new))
                    (f (grab-new))
                    ;; (work-vec (cl-mpm/utils::arb-vector (length (cl-mpm/aggregate::sim-agg-nodes-fd sim))))
                    ;; (work-vec-agg (cl-mpm/utils::arb-vector (length (cl-mpm/aggregate::sim-agg-nodes-fdc sim))))
                    )
               (cl-mpm/utils::resize-vector fg (cl-mpm/utils::sparse-matrix-nrows e))
               (cl-mpm/utils::resize-vector f (cl-mpm/utils::sparse-matrix-nrows et))
               (cl-mpm/aggregate::aggregate-vec
                sim (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-force d fg) d
                f)
               (cl-mpm/utils::resize-vector work-vec (cl-mpm/utils::sparse-matrix-nrows e))
               (cl-mpm/utils::resize-vector work-vec-agg (cl-mpm/utils::sparse-matrix-nrows et))
               (let* ((acc
                        (cl-mpm/linear-solver::solve-conjugant-gradients
                         (lambda (v)
                           (cl-mpm/fastmaths::fast-@-sparse-mat-dense-vec-masked-multithread
                            e
                            v
                            bcs
                            bcs-int
                            work-vec
                            )
                           (cl-mpm/fastmaths::fast-.*
                            sma
                            work-vec
                            work-vec)
                           (cl-mpm/fastmaths::fast-@-sparse-mat-dense-vec-masked-multithread
                            et
                            work-vec
                            bcs-int
                            bcs
                            work-vec-agg))
                         f
                         :tol 1d-15
                         :max-iters 10000
                         :mask bcs-int
                         )))
                 (cl-mpm/aggregate::zero-global sim #'cl-mpm/mesh::node-acceleration d)
                 (cl-mpm/aggregate::project-int-vec sim acc #'cl-mpm/mesh::node-acceleration d)))))))
      (when (= (mod solve-count damping-update-count) 0)
        (setf damping (the double-float (cl-mpm/dynamic-relaxation::dr-estimate-damping sim))))
      (cl-mpm:iterate-over-nodes
       mesh
       (lambda (node)
         (when (cl-mpm/mesh:node-active node)
           (with-accessors ((mass cl-mpm::node-mass)
                            (vel cl-mpm::node-velocity)
                            (force cl-mpm::node-force)
                            (internal cl-mpm/mesh::node-interior)
                            (agg cl-mpm/mesh::node-agg)
                            (acc cl-mpm::node-acceleration))
               node
             (when (or (not agg)
                       internal)
               (cl-mpm::integrate-vel-midpoint vel acc mass mass-scale dt damping))))))
      (cl-mpm::apply-essential-bcs sim)

      (with-accessors ((ke-prev sim-ke-prev)
                       (ke sim-ke)
                       (ke-damping sim-kinetic-damping))
          sim
        (when ke-damping
          (setf ke (calculate-ke sim))
          (setf damping 0d0)
          ;; (setf ke (cl-mpm::sim-stats-energy sim))
          (when (> ke-prev ke)
            (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
            (cl-mpm::reset-grid-acceleration (cl-mpm:sim-mesh sim))
            ;; (cl-mpm::reset-nodes-force sim)
            (setf ke 0d0))
          (setf ke-prev ke)))
      ;; (cl-mpm::apply-essential-bcs sim)
      )))
