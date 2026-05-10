(in-package :cl-mpm/dynamic-relaxation)

(defun implicit-assemble-stiffness (sim)
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (mps cl-mpm::sim-mps))
      sim
    (with-accessors ((nd cl-mpm/mesh::mesh-nd))
        mesh
      (cl-mpm::iterate-over-nodes
       mesh
       (lambda (node)
         (when (cl-mpm/mesh:node-active node)
           (setf (cl-mpm/mesh:node-mass node) 0d0))))
      (let ((mass-scale (the double-float (/ 1d0 (the (double-float ) (cl-mpm::sim-dt-scale sim))))))
        (declare (double-float mass-scale))
        (cl-mpm::iterate-over-mps
         mps
         (lambda (mp)
           (let ((stiffness (cl-mpm/implicit::assemble-mp-stiffness mesh mp))
                 (df-inv (cl-mpm/particle::mp-deformation-gradient-increment-inverse mp))
                 (mp-volume (cl-mpm/particle::mp-volume mp))
                 (g-a (cl-mpm/utils::stretch-dsvp-3d-zeros))
                 (g-b (cl-mpm/utils::stretch-dsvp-3d-zeros)))
             (declare (double-float mp-volume))
             (cl-mpm::iterate-over-neighbours
              mesh
              mp
              (lambda (node svp grads fsvp fgrads)
                (with-accessors ((node-active  cl-mpm/mesh:node-active)
                                 (node-lock  cl-mpm/mesh:node-lock))
                    node
                  (declare (boolean node-active)
                           (sb-thread:mutex node-lock))
                  (when node-active
                    (let ((g-a (cl-mpm/implicit::assemble-g-3d-prealloc
                                (cl-mpm::gradient-push-forwards-cached
                                 grads
                                 df-inv)
                                g-a)))
                      (cl-mpm::iterate-over-neighbours
                       mesh
                       mp
                       (lambda (node-b svp-b grads-b fsvp fgrads)
                         (when (cl-mpm::node-active node-b)
                           (let ((g-b (cl-mpm/implicit::assemble-g-3d-prealloc
                                       (cl-mpm::gradient-push-forwards-cached
                                        grads-b
                                        df-inv)
                                       g-b)))
                             (let ((stiff (magicl:@ (magicl:transpose g-a) stiffness g-b)))
                               (sb-thread:with-mutex (node-lock)
                                 (loop for v across (cl-mpm/utils:fast-storage stiff)
                                       do
                                          (progn
                                            (incf (the double-float (cl-mpm/mesh::node-mass node))
                                                  (the double-float
                                                       (* 0.25d0
                                                          mass-scale
                                                          (* mp-volume
                                                             (the double-float (abs v))))))))))))))))))))))))))

(defun map-stiffness-quasi-static (sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps))
      sim
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (node)
       (declare (cl-mpm/mesh::node node))
       (when (cl-mpm/mesh:node-active node)
         (setf (the double-float (cl-mpm/mesh::node-mass node)) 0d0))))

    (let* ((h (cl-mpm/mesh::mesh-resolution mesh))
           (nd (cl-mpm/mesh::mesh-nd mesh))
           (mass-scale (the double-float (/ 1d0
                                            (the double-float (cl-mpm::sim-dt-scale sim))))))
      (declare (double-float h mass-scale)
               (fixnum nd))
      (cl-mpm::iterate-over-mps
       mps
       (lambda (mp)
         ;; (setf (cl-mpm/particle::mp-p-modulus mp) (cl-mpm/particle::estimate-stiffness mp))
         (let* ((mp-volume (cl-mpm/particle::mp-volume mp))
                (mp-pmod (cl-mpm/particle::estimate-stiffness mp))
                (ul (estimate-ul-enhancement mp nd))
                (mp-factor (* mp-pmod mp-volume ul mass-scale
                              (/ 1d0  (* h h)))))
           (declare (type double-float mp-factor mp-pmod ul mp-volume))
           (cl-mpm::iterate-over-neighbours-cached
            mesh mp
            (lambda (node svp grads fsvp fgrads)
              (declare
               (cl-mpm/particle:particle mp)
               (cl-mpm/mesh::node node)
               (ignore grads fsvp fgrads)
               (double-float svp))
              (declare (type double-float mp-pmod mp-volume ul))
              (with-slots ((node-active cl-mpm/mesh::active)
                           (node-mass cl-mpm/mesh::mass))
                  node
                (declare (type double-float mp-volume mp-pmod svp)
                         (double-float node-mass))
                (when node-active
                  (sb-thread::with-mutex ((cl-mpm/mesh::node-lock node))
                    (incf node-mass
                          (the double-float
                               (* 2d0
                                  mp-factor
                                  svp))))))))))))))

(defmethod map-stiffness ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-ul))
  (map-stiffness-quasi-static sim))

(defgeneric update-node-fictious-mass (sim))

(defmethod update-node-fictious-mass ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr))
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (dt cl-mpm::sim-dt)
                   (bcs-force-list cl-mpm::sim-bcs-force-list))
      sim
    (map-stiffness-quasi-static sim)
    ;; (implicit-assemble-stiffness sim)
    (loop for bcs-f in bcs-force-list
          do (loop for bc across bcs-f
                   do (cl-mpm/bc::assemble-bc-stiffness sim bc)))
    (cl-mpm/ghost::apply-ghost-stiffness sim)
    (cl-mpm/aggregate::update-mass-matrix sim)
    (setf dt 1d0)))

(defmethod cl-mpm::finalise-loadstep ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-ul))
  ;;DR algorithm requires that finalisation is called once
  (setf (sim-initial-setup sim) nil)
  (cl-mpm::g2p (cl-mpm:sim-mesh sim)
               (cl-mpm:sim-mps sim)
               (cl-mpm:sim-dt sim)
               (cl-mpm:sim-damping-factor sim)
               (cl-mpm::sim-velocity-algorithm sim))
  (call-next-method)
  (when (cl-mpm::sim-allow-mp-damage-removal sim)
    (cl-mpm::remove-material-damaged sim))
  )
(defun midpoint-starter (sim)
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
         (cl-mpm::p2g-force-fs sim)
         (cl-mpm::apply-bcs mesh bcs-force dt-loadstep)
         (loop for bcs-f in bcs-force-list
               do (cl-mpm::apply-bcs mesh bcs-f dt-loadstep))
         (setf (cl-mpm::sim-damping-factor sim) 0d0)
         (update-node-fictious-mass sim)
         (cl-mpm/aggregate::update-node-forces-agg sim (* -0.5d0 dt))
         (cl-mpm::iterate-over-nodes
          mesh
          (lambda (n)
            (when (cl-mpm/mesh:node-active n)
              (cl-mpm/utils::vector-copy-into (cl-mpm/mesh::node-force n)
                                              (cl-mpm/mesh::node-residual n)))))
         (cl-mpm::apply-bcs mesh bcs dt-loadstep))))
(defmethod cl-mpm::reset-loadstep ((sim mpm-sim-dr-ul))
  (setf (sim-initial-setup sim) nil)
  (call-next-method))

(defgeneric pre-step (sim))
(defmethod pre-step ((sim mpm-sim-dr-ul))
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
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (initial-setup initial-setup)
               (enable-aggregate cl-mpm/aggregate::enable-aggregate)
               (damping cl-mpm::damping-factor)
               (ghost-factor cl-mpm::ghost-factor)
               (vel-algo cl-mpm::velocity-algorithm))
      sim
    (setf (cl-mpm/dynamic-relaxation::sim-solve-count sim) 0)
    (cl-mpm::apply-bcs mesh bcs dt)
    (cl-mpm::reset-grid mesh :reset-displacement t)
    (cl-mpm::reset-node-displacement sim)
    (cl-mpm::p2g mesh mps vel-algo)
    (setf (cl-mpm::sim-dt sim) 1d0)
    (when (> mass-filter 0d0)
      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
    ;; (cl-mpm::compact-mesh-active sim)
    (cl-mpm::filter-cells sim)
    (when ghost-factor
      (cl-mpm/ghost::build-ghost-cache sim))
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (n)
       (when (cl-mpm/mesh::node-active n)
         (setf
          (cl-mpm/mesh::node-true-mass n) (cl-mpm/mesh:node-mass n)))
       (cl-mpm/fastmaths:fast-zero (cl-mpm/mesh::node-true-velocity n))))
    (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
    (setf (cl-mpm::sim-damping-factor sim) 0d0)
    (midpoint-starter sim)
    (setf initial-setup t)))

(defmethod cl-mpm::update-nodes ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-ul))
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (dt cl-mpm::sim-dt)
                   (agg cl-mpm/aggregate::sim-enable-aggregate))
      sim
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (node)
       (when (and (cl-mpm/mesh:node-active node)
                  (or (not (cl-mpm/mesh::node-agg node))
                      (cl-mpm/mesh::node-interior node)))
         (cl-mpm::update-node node dt))))
    (when agg
      (cl-mpm/aggregate::project-displacement sim))))

(defmethod cl-mpm::update-sim ((sim mpm-sim-dr-ul))
  "Update stress last algorithm"
  (declare (cl-mpm::mpm-sim sim))
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
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (ghost-factor cl-mpm::ghost-factor)
               (initial-setup initial-setup)
               (enable-aggregate cl-mpm/aggregate::enable-aggregate)
               (damping cl-mpm::damping-factor)
               (damping-scale cl-mpm/dynamic-relaxation::damping-scale)
               (solve-count cl-mpm/dynamic-relaxation::solve-count)
               (mass-update-iter cl-mpm/dynamic-relaxation::mass-update-count)
               (vel-algo cl-mpm::velocity-algorithm))
      sim
    (declare (fixnum solve-count mass-update-iter))
    (unless initial-setup
      (pre-step sim)
      (with-accessors ((ke-prev sim-ke-prev)
                       (ke sim-ke))
          sim
        (setf ke 0d0
              ke-prev 0d0)))
    (setf dt 1d0)
    ;; ;; (cl-mpm/penalty::reset-penalty sim)
    (cl-mpm::reset-nodes-force sim)
    (cl-mpm::update-stress mesh mps dt-loadstep fbar)
    (cl-mpm::p2g-force-fs sim)
    (cl-mpm::apply-bcs mesh bcs-force dt)
    ;; ;; (loop for bcs-f in bcs-force-list
    ;; ;;       do (cl-mpm::apply-bcs mesh bcs-f dt-loadstep))

    ;; ;; (when ghost-factor
    ;; ;;   (cl-mpm/ghost::apply-ghost-cached sim)
    ;; ;;   (cl-mpm::apply-bcs mesh bcs dt))

    (when (= (mod solve-count mass-update-iter) 0)
      (update-node-fictious-mass sim))
    (incf solve-count)
    ;; ;; ;;Update our nodes after force mapping
    (update-node-forces-quasi-static sim)
    (cl-mpm::update-nodes sim)
    (cl-mpm::apply-bcs mesh bcs dt)
    (cl-mpm::update-filtered-cells sim)
    ;; ;; (cl-mpm::update-dynamic-stats sim)
    ;; ;; (cl-mpm::g2p mesh mps dt damping :TRIAL)
    (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC))
  )






(defmethod pre-step ((sim mpm-sim-dr-damage-ul))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps)
                   (delocal-counter cl-mpm/damage::sim-damage-delocal-counter-max))
      sim
    (call-next-method)
    (setf delocal-counter -1)
    (cl-mpm/damage::update-delocalisation-list mesh mps)))



(defmethod cl-mpm::update-sim ((sim mpm-sim-dr-damage-ul))
  "Update stress last algorithm"
  (declare (cl-mpm::mpm-sim sim))
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
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (ghost-factor cl-mpm::ghost-factor)
               (initial-setup initial-setup)
               (enable-aggregate cl-mpm/aggregate::enable-aggregate)
               (damping cl-mpm::damping-factor)
               (damping-scale cl-mpm/dynamic-relaxation::damping-scale)
               (vel-algo cl-mpm::velocity-algorithm)
               (solve-count cl-mpm/dynamic-relaxation::solve-count)
               (mass-update-iter cl-mpm/dynamic-relaxation::mass-update-count)
               )
      sim
    (declare (double-float damping-scale damping))
    (unless initial-setup
      (pre-step sim))
    (cl-mpm/penalty::reset-penalty sim)
    (setf dt 1d0)
    (cl-mpm::update-nodes sim)
    (cl-mpm::update-filtered-cells sim)
    (cl-mpm::reset-nodes-force sim)
    (cl-mpm::update-stress mesh mps dt-loadstep fbar)
    (cl-mpm/damage::calculate-damage sim dt-loadstep)
    (cl-mpm::p2g-force-fs sim)
    (cl-mpm::apply-bcs mesh bcs-force dt)

    (loop for bcs-f in bcs-force-list
          do (cl-mpm::apply-bcs mesh bcs-f dt-loadstep))

    (incf solve-count)
    (when ghost-factor
      (cl-mpm/ghost::apply-ghost-cached sim)
      (cl-mpm::apply-bcs mesh bcs dt))

    (when (= (mod solve-count mass-update-iter) 0)
      (update-node-fictious-mass sim))

    ;; ;;Update our nodes after force mapping
    (update-node-forces-quasi-static sim)
    (cl-mpm::apply-bcs mesh bcs dt)
    ;; (cl-mpm::update-dynamic-stats sim)
    ;; (cl-mpm::g2p mesh mps dt damping :TRIAL)
    (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
    ))
(defun update-node-forces-quasi-static (sim)
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
    ;; (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt)
    (cl-mpm:iterate-over-nodes
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (cl-mpm::calculate-forces-midpoint node 0d0 0d0 mass-scale))))
    (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt)
    ;; (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm::sim-active-bcs sim) dt)
    ;;For each aggregated element set solve mass matrix and velocity
    (when enable-aggregate
      (let* ()
        (cl-mpm/aggregate::iterate-over-dimensions
         (cl-mpm::mesh-nd mesh)
         (lambda (d)
           (let* ((f
                    (cl-mpm/aggregate::aggregate-vec
                     sim
                     (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-force d) d))
                  (et (cl-mpm/aggregate::sim-global-sparse-et sim))
                  (e (cl-mpm/aggregate::sim-global-sparse-e sim))
                  (sma (cl-mpm/aggregate::sim-global-sparse-ma sim))
                  (bcs-int (aref (cl-mpm/aggregate::sim-global-bcs-int sim) d))
                  (bcs (aref (cl-mpm/aggregate::sim-global-bcs sim) d))
                  (work-vec (cl-mpm/utils::arb-matrix (length (cl-mpm/aggregate::sim-agg-nodes-fd sim)) 1))
                  (work-vec-agg (cl-mpm/utils::arb-matrix (length (cl-mpm/aggregate::sim-agg-nodes-fdc sim)) 1)))
             (cl-mpm/aggregate::apply-internal-bcs sim f d)
             (let* ((acc
                      (cl-mpm/linear-solver::solve-conjugant-gradients
                       (lambda (v)
                         (cl-mpm/fastmaths::fast-@-sparse-mat-dense-vec-masked
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
                         (cl-mpm/fastmaths::fast-@-sparse-mat-dense-vec-masked
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
               (cl-mpm/aggregate::project-int-vec sim acc #'cl-mpm/mesh::node-acceleration d))))))
      )
    (when (= (mod solve-count damping-update-count) 0)
      (setf damping (the double-float (cl-mpm/dynamic-relaxation::dr-estimate-damping sim))))
    ;; ;; (with-accessors ((ke-prev sim-ke-prev)
    ;; ;;                  (ke sim-ke)
    ;; ;;                  (ke-damping sim-kinetic-damping))
    ;; ;;     sim
    ;; ;;   (when ke-damping
    ;; ;;     (setf ke (calculate-ke sim))
    ;; ;;     ;; (setf ke (cl-mpm::sim-stats-energy sim))
    ;; ;;     (when (> ke-prev ke)
    ;; ;;       (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
    ;; ;;       (cl-mpm::reset-nodes-force sim)
    ;; ;;       (setf ke 0d0))
    ;; ;;     (setf ke-prev ke)))
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
    (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm::sim-active-bcs sim) dt)
    ))
(defmethod cl-mpm::update-node-forces ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-ul))
  (update-node-forces-quasi-static sim))




(defmethod cl-mpm::finalise-loadstep ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-usf))
  ;;DR algorithm requires that finalisation is called once
  (reset-mp-velocity sim)
  (call-next-method)
  )

(defmethod cl-mpm::update-sim ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-usf))
  "update stress first algorithm"
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
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (ghost-factor cl-mpm::ghost-factor)
               (initial-setup initial-setup)
               (enable-aggregate cl-mpm/aggregate::enable-aggregate)
               (damping cl-mpm::damping-factor)
               (damping-scale cl-mpm/dynamic-relaxation::damping-scale)
               (time cl-mpm::time)
               (vel-algo cl-mpm::velocity-algorithm))
                sim
    (progn
      (cl-mpm::reset-grid mesh)
      (when (> (length mps) 0)
        (cl-mpm::p2g mesh mps vel-algo)
        (when (> mass-filter 0d0)
          (cl-mpm::filter-grid mesh mass-filter))
        (setf (cl-mpm:sim-dt sim) (* (cl-mpm::sim-dt-scale sim) (cl-mpm::calculate-min-dt sim)))
        (cl-mpm::filter-cells sim)
        (cl-mpm::apply-bcs mesh bcs dt)
        (cl-mpm::update-node-kinematics sim)
        (cl-mpm::apply-bcs mesh bcs dt)
        ;;trial update displacements
        (cl-mpm::update-nodes sim)
        (cl-mpm::update-cells sim)

        (cl-mpm::update-stress mesh mps dt fbar)
        ;; (update-node-fictious-mass sim)
        ;; (setf (cl-mpm:sim-dt sim) 1d0)
        ;; map forces onto nodes
        (cl-mpm::p2g-force sim)
        (when bcs-force-list
          (loop for bcs-f in bcs-force-list
                do (cl-mpm::apply-bcs mesh bcs-f dt-loadstep)))

        (cl-mpm::update-node-forces sim)
        (cl-mpm::reset-node-displacement sim)
        (cl-mpm::update-nodes sim)
        (cl-mpm::apply-bcs mesh bcs dt)
        ;; (cl-mpm::update-dynamic-stats sim)
        ;; Also updates mps inline
        (cl-mpm::g2p mesh mps dt damping vel-algo)
        (cl-mpm::update-particles sim)
        (cl-mpm::new-loadstep sim)
        )
      (incf time dt))))

(defmethod estimate-static-oobf ((sim mpm-sim-implict-dynamic))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (sim-agg cl-mpm/aggregate::sim-enable-aggregate))
      sim
    (destructuring-bind (energy
                         oobf-num
                         oobf-denom
                         power)
        (cl-mpm::reduce-over-nodes
         (cl-mpm:sim-mesh sim)
         (lambda (node)
           (if (and (cl-mpm/mesh:node-active node)
                    (not (cl-mpm/mesh::node-agg node)))
               (with-accessors ((active cl-mpm/mesh::node-active)
                                (f-ext cl-mpm/mesh::node-external-force)
                                (f-int cl-mpm/mesh::node-internal-force)
                                (node-oobf cl-mpm/mesh::node-oobf)
                                (mass cl-mpm/mesh::node-mass)
                                (volume cl-mpm/mesh::node-volume)
                                (volume-t cl-mpm/mesh::node-volume-true)
                                (vel cl-mpm/mesh::node-velocity)
                                (disp cl-mpm/mesh::node-displacment)
                                )
                   node
                 (declare (double-float mass))
                 (let ()
                   (list
                    (* 0.5d0 mass (cl-mpm/fastmaths::mag-squared vel))
                    (cl-mpm/fastmaths::mag-squared (cl-mpm/fastmaths::fast-.+-vector f-ext f-int))
                    (cl-mpm/fastmaths::mag-squared f-ext)
                    (cl-mpm/fastmaths:dot disp f-ext))))
               (list 0d0 0d0 0d0 0d0)))
         (lambda (a b) (mapcar (lambda (x y) (declare (double-float x y)) (+ x y)) a b)))
      (declare (double-float  energy oobf-num oobf-denom power))
      (when sim-agg
        (cl-mpm/aggregate::iterate-over-dimensions-serial
         (cl-mpm/mesh::mesh-nd mesh)
         (lambda (d)
           (let* ((f-int (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-internal-force d))
                  (f-ext (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-external-force d))
                  (E (cl-mpm/aggregate::sim-global-e sim))
                  (ma (cl-mpm/aggregate::sim-global-ma sim))
                  (vi (magicl:@ (magicl:transpose E) (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-velocity d)))
                  (disp (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-displacment d)))
             (incf oobf-num (cl-mpm/fastmaths::mag-squared
                             (cl-mpm/aggregate::apply-internal-bcs
                              sim
                              (magicl:@ (magicl:transpose E)
                                        (cl-mpm/fastmaths::fast-.+
                                         f-int
                                         f-ext))
                              d
                              )))
             (incf oobf-denom (cl-mpm/fastmaths::mag-squared
                               (cl-mpm/aggregate::apply-internal-bcs
                                sim
                                (magicl:@ (magicl:transpose E) f-ext)
                                d)))
             (incf power (cl-mpm/fastmaths:dot
                          disp f-ext))
             (incf energy (* 0.5d0 (cl-mpm/utils::mtref (magicl:@ (magicl:transpose vi) ma vi) 0 0)))
             ))))
      (let ((oobf 0d0)
            (oobf-num (sqrt oobf-num))
            (oobf-denom (sqrt oobf-denom))
            )
        (if (> oobf-denom 0d0)
            (setf oobf (/ oobf-num oobf-denom))
            (setf oobf (if (> oobf-num 0d0) sb-ext:double-float-positive-infinity 0d0)))
        oobf))))
