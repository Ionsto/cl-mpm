(in-package :cl-mpm/dynamic-relaxation)



(defun dr-calculate-forces-implicit-dynamic (node true-damping dt mass-scale real-dt)
  "Update forces and nodal velocities with viscous damping"
  (when (cl-mpm/mesh:node-active node)
    (with-accessors ((mass cl-mpm/mesh:node-mass)
                     (vel cl-mpm/mesh:node-velocity)
                     (disp cl-mpm/mesh::node-displacment)
                     (vt cl-mpm/mesh::node-true-velocity)
                     (real-mass cl-mpm/mesh::node-true-mass)
                     (force cl-mpm/mesh:node-force)
                     (force-ext cl-mpm/mesh::node-external-force)
                     (force-int cl-mpm/mesh::node-internal-force)
                     (force-damp cl-mpm/mesh::node-damping-force)
                     (force-ghost cl-mpm/mesh::node-ghost-force)
                     (residual cl-mpm/mesh::node-residual)
                     (residual-prev cl-mpm/mesh::node-residual-prev)
                     (acc cl-mpm/mesh::node-acceleration))
        node
      (declare (double-float mass dt true-damping mass-scale))
      (progn
        (cl-mpm/fastmaths:fast-zero acc)
        ;;Backwards euler acceleration
        (let ((vn (cl-mpm/utils:vector-zeros))
              (inertia-force (cl-mpm/utils:vector-zeros)))
          (unless (= real-dt 0d0)
            (cl-mpm/utils:vector-copy-into disp vn)
            (cl-mpm/fastmaths:fast-scale! vn (/ 1d0 real-dt))
            (cl-mpm/fastmaths::fast-.- vn vt inertia-force)
            (cl-mpm/fastmaths::fast-scale! inertia-force (/ (- real-mass) real-dt))
            (cl-mpm/fastmaths:fast-fmacc force-damp vn (* true-damping -1d0 real-mass))
            ;; (setf (cl-mpm/mesh::node-inertia-force node) inertia-force)
            )
          ;;Set acc to f/m
          (cl-mpm/fastmaths::fast-.+ inertia-force force force)
          (cl-mpm/fastmaths::fast-.+-vector force-int force force)
          (cl-mpm/fastmaths::fast-.+-vector force-ext force force)
          (cl-mpm/fastmaths::fast-.+-vector force-ghost force force)
          (cl-mpm/fastmaths::fast-.+-vector force-damp force force)
          (when (> mass 0d0)
            (cl-mpm/fastmaths:fast-fmacc acc force (/ 1d0 (* mass-scale mass))))
          ;;Would intergrtate but don't
          (cl-mpm/utils::vector-copy-into residual residual-prev)
          (cl-mpm/utils::vector-copy-into force residual)

          ;; (cl-mpm/fastmaths::fast-.+-vector force-ext inertia-force force-ext)
          ;; (cl-mpm/utils::vector-copy-into force-int residual)
          ;; (cl-mpm/fastmaths::fast-.+ inertia-force residual residual)
          ;; (cl-mpm/fastmaths::fast-.+-vector force-ghost residual residual)
          ))))
  (values))


(cl-mpm/utils::with-arb-pool
    (defmethod cl-mpm::update-node-forces ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic))
      (with-accessors ((mesh cl-mpm:sim-mesh)
                       (mass-scale cl-mpm::sim-mass-scale)
                       (damping cl-mpm::sim-damping-factor)
                       (true-damping cl-mpm/dynamic-relaxation::sim-true-damping)
                       (damping-algo cl-mpm::sim-damping-algorithm)
                       (agg-elems cl-mpm/aggregate::sim-agg-elems)
                       ;; (dt cl-mpm::sim-dt)
                       (dt-loadstep cl-mpm/dynamic-relaxation::sim-dt-loadstep)
                       (enable-dynamics sim-enable-dynamics)
                       (damping-scale cl-mpm/dynamic-relaxation::sim-damping-scale)
                       (solve-count cl-mpm/dynamic-relaxation::sim-solve-count)
                       (damping-update-count sim-damping-update-count)
                       (enable-aggregate cl-mpm/aggregate::sim-enable-aggregate))
          sim
        (declare (double-float damping damping-scale))
        (cl-mpm:iterate-over-nodes
         mesh
         (lambda (node)
           (when (cl-mpm/mesh:node-active node)
             (if (and
                  enable-dynamics
                  (> dt-loadstep 0d0))
                 (dr-calculate-forces-implicit-dynamic node true-damping 0d0 mass-scale dt-loadstep)
                 (cl-mpm::calculate-forces-midpoint node 0d0 0d0 mass-scale)))))
        (cl-mpm::compute-reaction-force sim)
        (cl-mpm::apply-essential-bcs sim)
        ;; (when enable-aggregate
        ;;   (let* ((ma (cl-mpm/aggregate::sim-global-ma sim)))
        ;;     (cl-mpm/aggregate::iterate-over-dimensions
        ;;      (cl-mpm/mesh:mesh-nd mesh)
        ;;      (lambda (d)
        ;;        (let ((f (cl-mpm/aggregate::aggregate-vec
        ;;                  sim
        ;;                  (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-force d)
        ;;                  d)))
        ;;          ;; (cl-mpm/aggregate::apply-internal-bcs sim f d)
        ;;          (let* ((acc (cl-mpm/aggregate::linear-solve-with-bcs sim ma f d)))
        ;;            (cl-mpm/aggregate::zero-global sim #'cl-mpm/mesh::node-acceleration d)
        ;;            (cl-mpm/aggregate::project-int-vec
        ;;             sim
        ;;             ;; (cl-mpm/aggregate::extend-vec sim acc d)
        ;;             acc
        ;;             #'cl-mpm/mesh::node-acceleration d)
        ;;            ))))))
        (when enable-aggregate
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
                 ;; (cl-mpm/fastmaths::fast-scale! acc (/ 1d0 mass-scale))
                 (cl-mpm/aggregate::zero-global sim #'cl-mpm/mesh::node-acceleration d)
                 (cl-mpm/aggregate::project-int-vec sim acc #'cl-mpm/mesh::node-acceleration d))))))

        (when (= (mod solve-count damping-update-count) 0)
          (setf damping (the double-float (cl-mpm/dynamic-relaxation::dr-estimate-damping sim))))

        (cl-mpm:iterate-over-nodes
         mesh
         (lambda (node)
           (when (cl-mpm/mesh:node-active node)
             (with-accessors ((mass cl-mpm/mesh::node-mass)
                              (vel cl-mpm/mesh::node-velocity)
                              (force cl-mpm/mesh::node-force)
                              (internal cl-mpm/mesh::node-interior)
                              (agg cl-mpm/mesh::node-agg)
                              (acc cl-mpm/mesh::node-acceleration))
                 node
               (when (or internal
                         (not agg))
                 (cl-mpm::integrate-vel-midpoint vel acc mass mass-scale 1d0 damping))))))
        (cl-mpm::apply-essential-bcs sim)
        ;; (cl-mpm/aggregate::project-displacement sim)
        ;; (cl-mpm::apply-essential-bcs sim)
        )))

(defmethod cl-mpm::update-nodes ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic))
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   ;; (dt cl-mpm/dynamic-relaxation::sim-dt-loadstep)
                   (agg cl-mpm/aggregate::sim-enable-aggregate))
      sim
    (let ((dt 1d0))
      (cl-mpm::iterate-over-nodes
       mesh
       (lambda (node)
         (when (and (cl-mpm/mesh:node-active node)
                    (or (not (cl-mpm/mesh::node-agg node))
                        (cl-mpm/mesh::node-interior node)))
           (cl-mpm::update-node node dt))))
      (when agg
        (cl-mpm/aggregate::project-displacement sim)))))

(defmethod map-stiffness ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (dt-true cl-mpm/dynamic-relaxation::sim-dt-loadstep)
                   (mps cl-mpm:sim-mps))
      sim
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (node)
       (declare (cl-mpm/mesh::node node))
       (when t;(cl-mpm/mesh:node-active node)
         (setf (the double-float (cl-mpm/mesh::node-mass node)) 0d0))))

    (let* ((h (cl-mpm/mesh::mesh-resolution mesh))
           (nd (cl-mpm/mesh::mesh-nd mesh))
           (mass-scale (the double-float (/ 1d0 (the double-float (cl-mpm::sim-dt-scale sim))))))
      (declare (double-float h mass-scale)
               (fixnum nd))
      (cl-mpm::iterate-over-mps
       mps
       (lambda (mp)
         (setf (cl-mpm/particle::mp-p-modulus mp)
               (cl-mpm/particle::estimate-stiffness mp))
         (let* ((mp-volume (cl-mpm/particle::mp-volume mp))
                (mp-mass (cl-mpm/particle::mp-mass mp))
                (mp-pmod (cl-mpm/particle::mp-p-modulus mp))
                (ul (estimate-ul-enhancement mp nd))
                (mp-factor
                  (+
                   (* 16d0
                      ul
                      mass-scale
                      (+
                       (* (/ 1d0 (expt dt-true 1))
                          (sim-true-damping sim)
                          mp-mass)
                       (* (/ 1d0 (expt dt-true 2))
                          mp-mass)))
                   (*
                    mp-pmod
                    mp-volume
                    ul
                    mass-scale
                    (/ 1d0 (* h h))))))
           (declare (type double-float mp-factor mp-pmod ul mp-volume))
           (cl-mpm::iterate-over-neighbours
            mesh mp
            (lambda (node svp grads fsvp fgrads)
              (declare
               (cl-mpm/particle:particle mp)
               (cl-mpm/mesh::node node)
               ;; (ignore grads fsvp fgrads)
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
                               (* 0.25d0
                                  mp-factor))))))))))))))

(defmethod update-node-fictious-mass ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic))
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (dt-true cl-mpm/dynamic-relaxation::sim-dt-loadstep)
                   (enable-dynamics sim-enable-dynamics)
                   (bcs-force-list cl-mpm::sim-bcs-force-list))
      sim
    (map-stiffness sim)
    (loop for bcs-f in bcs-force-list
          do (loop for bc across bcs-f
                   do (cl-mpm/bc::assemble-bc-stiffness sim bc)))

    (cl-mpm/ghost::apply-ghost-stiffness sim)

    ;; (when (and enable-dynamics
    ;;            (> dt-true 0d0))
    ;;   (let ((mass-scale (the double-float (/ 1d0 (cl-mpm::sim-dt-scale sim)))))
    ;;     (cl-mpm:iterate-over-nodes
    ;;      mesh
    ;;      (lambda (node)
    ;;        (when (cl-mpm/mesh:node-active node)
    ;;          (incf (cl-mpm/mesh:node-mass node)
    ;;                (* mass-scale
    ;;                   (+
    ;;                    (* (/ 1d0 (expt dt-true 1))
    ;;                       (sim-true-damping sim)
    ;;                       (cl-mpm/mesh::node-true-mass node))
    ;;                    (* (/ 1d0 (expt dt-true 2))
    ;;                       (cl-mpm/mesh::node-true-mass node))))))))))

    (cl-mpm/aggregate::update-mass-matrix sim)
    (setf (cl-mpm:sim-dt sim) 1d0)
    ))

(defmethod cl-mpm::update-sim ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic))
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
    ;; (setf dt 1d0)
    (cl-mpm/penalty::reset-penalty sim)
    (cl-mpm::reset-nodes-force sim)
    (cl-mpm::apply-force-bcs sim dt-loadstep)
    (cl-mpm::apply-essential-bcs sim)
    (cl-mpm::update-stress mesh mps dt-loadstep fbar)
    (cl-mpm/damage::calculate-damage sim dt-loadstep)
    (cl-mpm::p2g-force-fs sim)
    ;; (when ghost-factor
    ;;   (cl-mpm/ghost::apply-ghost sim ghost-factor)
    ;;   (cl-mpm::apply-bcs mesh bcs dt))
    (when (= (mod solve-count mass-update-iter) 0)
      (update-node-fictious-mass sim))
    ;; ;;Update our nodes after force mapping
    (cl-mpm::update-node-forces sim)
    (cl-mpm::apply-essential-bcs sim)
    (cl-mpm::update-nodes sim)
    (cl-mpm::update-filtered-cells sim)
    ;; (cl-mpm::update-dynamic-stats sim)
    ;; (cl-mpm::g2p mesh mps dt damping :TRIAL)
    (cl-mpm::update-dynamic-stats sim)
    (incf solve-count)
    ))

(defmethod cl-mpm::finalise-loadstep ((sim mpm-sim-dr-dynamic))
  ;;DR algorithm requires that finalisation is called once
  ;; (setf (cl-mpm/dynamic-relaxation::sim-sub-stepping sim) nil)
  (with-accessors ((real-dt cl-mpm/dynamic-relaxation::sim-dt-loadstep)
                   (mesh cl-mpm::sim-mesh))
      sim
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (n)

       (with-accessors ((disp cl-mpm/mesh::node-displacment)
                        (vel cl-mpm/mesh:node-velocity)
                        (vel-n cl-mpm/mesh::node-true-velocity)
                        (acc cl-mpm/mesh::node-acceleration))
           n
         (unless (= real-dt 0d0)
           (cl-mpm/utils:vector-copy-into disp vel)
           (cl-mpm/fastmaths:fast-scale! vel (/ 1d0 real-dt))
           (cl-mpm/fastmaths::fast-.- vel vel-n acc)
           (cl-mpm/fastmaths::fast-scale! acc (/ 1d0 real-dt))))
       (cl-mpm/utils:vector-copy-into (cl-mpm/mesh::node-velocity n) (cl-mpm/mesh::node-true-velocity n))
       (setf (cl-mpm/mesh:node-mass n) (cl-mpm/mesh::node-true-mass n))
       )))
  (cl-mpm/aggregate::update-mass-matrix sim)
  (cl-mpm::update-dynamic-stats sim)
  (call-next-method))


(defmethod cl-mpm::calculate-min-dt ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic))
  ;;we dont care about the stiffness from our BC constraints at all!
  (cl-mpm::calculate-min-dt-mps sim))

(defmethod cl-mpm/setup::%estimate-elastic-dt ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic))
  (with-accessors ((mps cl-mpm:sim-mps)
                   (bcs-force-list cl-mpm:sim-bcs-force-list))
      sim
    (if (> (length mps) 0)
        (progn
          (min
           (cl-mpm/setup::%estimate-elastic-dt-mps sim)))
        sb-ext:double-float-positive-infinity)))




(defmethod pre-step ((sim mpm-sim-dr-dynamic))
  ;; (break)
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
    ;; (midpoint-starter sim)
    (setf initial-setup t)))



(defun super-step (sim)
  (let ((crit (cl-mpm/dynamic-relaxation::sim-convergence-critera sim))
        (dt (cl-mpm::sim-dt sim))
        (dt-scale (cl-mpm::sim-dt-scale sim))
        (true-damping (cl-mpm::sim-damping-factor sim)))
    (setf (cl-mpm/dynamic-relaxation::sim-sub-stepping sim) t)
    (setf (sim-dt-loadstep sim) (* 1d0 dt))
    ;; (change-class sim 'cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic)
    (setf (cl-mpm/dynamic-relaxation::sim-true-damping sim) true-damping)
    (labels ((try-step ()
               (handler-case
                   (let (;; (conv-crit 1d-3)
                         (conv-crit crit)
                         (residual-normaliser nil)
                         (substeps 20))
                     (generalised-staggered-solve
                      sim
                      :crit conv-crit
                      :substeps substeps
                      :sub-conv-steps 100
                      :dt-scale 0.9d0
                      :damping (sqrt 2d0)
                      ;; :convergance-criteria
                      ;; (lambda (sim f o)
                      ;;   (let ((c ;; (cl-mpm/dynamic-relaxation::res-norm-aggregated sim)
                      ;;            (cl-mpm/dynamic-relaxation::vel-norm-aggregated sim)))
                      ;;     (if residual-normaliser
                      ;;         (setf c (/ c residual-normaliser))
                      ;;         (when (> c 0d0)
                      ;;           (setf residual-normaliser c
                      ;;                 c 1d0)))
                      ;;     (format t "Conv crit: ~E - norm ~E~%" c residual-normaliser)
                      ;;     (< c conv-crit)))
                      :max-damage-inc 0.1d0
                      :max-plastic-inc nil;10d0
                      :post-iter-step (lambda (i e o)
                                        (format t "Dynamic substep ~D~%" i)

                                        ;; (check-damage-increment sim)
                                        (when nil;(uiop:directory-exists-p "./output/")
                                          (cl-mpm/output:save-vtk (merge-pathnames "./output/" (format nil "rsim_step_~5,'0d.vtk" i)) sim)
                                          (cl-mpm/output:save-vtk-nodes (merge-pathnames "./output/" (format nil "rsim_step_nodes_~5,'0d.vtk" i)) sim)
                                          (save-conv-step sim "./output/" *total-iter* *total-step* 0d0 o e)
                                          )
                                        ;; (save-conv-step sim "./output/" *total-iter* *total-step* 0d0 o e)
                                        (incf *total-iter* substeps)))
                     t)
                 (cl-mpm/errors:error-simulation (c)
                   (princ c)
                   ;; (cl-mpm::reset-loadstep sim)
                   (cl-mpm:iterate-over-mps
                    (cl-mpm:sim-mps sim)
                    (lambda (mp)
                      (cl-mpm/particle::reset-loadstep-mp mp)))
                   (cl-mpm::reset-node-displacement sim)
                   nil))))
      (let ((conv nil)
            (dt-step (- dt 1d-15)))
        (loop while (> dt-step 0d0)
              do
                 (let ((max-dt-refines 8))
                  (loop for i from 0 to max-dt-refines
                        while (not conv)
                        do (progn
                             (setf conv (try-step))
                             (unless conv
                               (format t "Implicit dynamic refine ~D - dt step ~A~%" i (sim-dt-loadstep sim))
                               (setf (sim-dt-loadstep sim) (/ (sim-dt-loadstep sim)
                                                              2d0))))
                           finally (when (= i max-dt-refines)
                                     (error "Auto-stepping implicit dynamic scheme failed")))
                   (decf dt-step (sim-dt-loadstep sim))))))
    (incf *total-step*)
    (setf (cl-mpm::sim-dt-scale sim) dt-scale)
    (setf (cl-mpm::sim-dt sim) dt)
    (setf (cl-mpm::sim-damping-factor sim) true-damping)
    (setf (cl-mpm/dynamic-relaxation::sim-sub-stepping sim) nil)
    (cl-mpm::finalise-loadstep sim))
  )


(defparameter *total-iter* 0)
(defparameter *total-step* 0)
(defmethod cl-mpm::update-sim ((sim cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic))
  "Update stress last algorithm"
  (if (cl-mpm/dynamic-relaxation::sim-sub-stepping sim)
      (call-next-method)
      (super-step sim)))
