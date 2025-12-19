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
          ;; (cl-mpm/fastmaths:fast-fmacc force-damp vn (* damping -1d0 mass))
          (unless (= real-dt 0d0)
            (cl-mpm/utils:vector-copy-into disp vn)
            (cl-mpm/fastmaths:fast-scale! vn (/ 1d0 real-dt))
            (cl-mpm/fastmaths::fast-.- vn vt inertia-force)
            (cl-mpm/fastmaths::fast-scale! inertia-force (/ (- real-mass) real-dt))
            ;; (cl-mpm/fastmaths::fast-.+ inertia-force force force)
            )
          ;;Set acc to f/m
          (cl-mpm/fastmaths::fast-.+ inertia-force force force)
          (cl-mpm/fastmaths::fast-.+-vector force-int force force)
          (cl-mpm/fastmaths::fast-.+-vector force-ext force force)
          (cl-mpm/fastmaths::fast-.+-vector force-ghost force force)
          (cl-mpm/fastmaths:fast-fmacc acc force (/ 1d0 mass))
          ;;Would intergrtate but don't
          (cl-mpm/utils::vector-copy-into residual residual-prev)
          (cl-mpm/utils::vector-copy-into force residual)

          ;; (cl-mpm/fastmaths::fast-.+-vector force-ext inertia-force force-ext)
          ;; (cl-mpm/utils::vector-copy-into force-int residual)
          ;; (cl-mpm/fastmaths::fast-.+ inertia-force residual residual)
          ;; (cl-mpm/fastmaths::fast-.+-vector force-ghost residual residual)
          ))))
  (values))



(defmethod cl-mpm::update-node-forces ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mass-scale cl-mpm::sim-mass-scale)
                   (damping cl-mpm::sim-damping-factor)
                   (damping-algo cl-mpm::sim-damping-algorithm)
                   (agg-elems cl-mpm/aggregate::sim-agg-elems)
                   (dt cl-mpm::sim-dt)
                   (dt-loadstep cl-mpm/dynamic-relaxation::sim-dt-loadstep)
                   (enable-dynamics sim-enable-dynamics)
                   (enable-aggregate cl-mpm/aggregate::sim-enable-aggregate))
      sim
    (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt)
    (cl-mpm:iterate-over-nodes
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (if enable-dynamics
             (cl-mpm/dynamic-relaxation::dr-calculate-forces-implicit-dynamic node 0d0 0d0 mass-scale dt-loadstep)
             (cl-mpm::calculate-forces node 0d0 0d0 mass-scale)))))

    (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt)
    ;;For each aggregated element set solve mass matrix and velocity
    (when enable-aggregate
      (let* ((E (cl-mpm/aggregate::sim-global-e sim))
             (ma (cl-mpm/aggregate::sim-global-ma sim)))
        (loop for d from 0 below (cl-mpm::mesh-nd mesh)
              do (let ((f (magicl:@ (magicl:transpose E)
                                    (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-force d))))
                   (cl-mpm/aggregate::apply-internal-bcs sim f d)
                   (let* ((acc (magicl:linear-solve ma f)))
                     (cl-mpm/aggregate::apply-internal-bcs sim acc d)
                     ;; (project-global-vec sim (magicl:@ E acc) #'cl-mpm/mesh::node-acceleration d)
                     (cl-mpm/aggregate::zero-global sim #'cl-mpm/mesh::node-acceleration d)
                     (cl-mpm/aggregate::project-int-vec sim acc #'cl-mpm/mesh::node-acceleration d)))))

      (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt))
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
             (cl-mpm::integrate-vel-midpoint vel acc mass mass-scale dt damping))))))
    (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt)))

(defmethod update-node-fictious-mass ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic))
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (dt cl-mpm::sim-dt)
                   (dt-true cl-mpm/dynamic-relaxation::sim-dt-loadstep)
                   (enable-dynamics sim-enable-dynamics)
                   (bcs-force-list cl-mpm::sim-bcs-force-list)
                   )
      sim
    (map-stiffness sim)
    (cl-mpm/penalty::assemble-penalty-stiffness-matrix sim)
    (loop for bcs-f in bcs-force-list
          do (loop for bc across bcs-f
                   do (cl-mpm/bc::assemble-bc-stiffness sim bc)))
    (when enable-dynamics
      (cl-mpm:iterate-over-nodes
       mesh
       (lambda (node)
         (when (cl-mpm/mesh:node-active node)
           (incf (cl-mpm/mesh:node-mass node)
                 (* (/ 2d0 (expt dt-true 2))
                    (cl-mpm/mesh::node-true-mass node)))))))

    (cl-mpm/aggregate::update-mass-matrix sim)
    (setf dt 1d0)))

(defmethod cl-mpm::update-sim ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic))
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
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (ghost-factor cl-mpm::ghost-factor)
               (initial-setup initial-setup)
               (enable-aggregate cl-mpm/aggregate::enable-aggregate)
               (damping cl-mpm::damping-factor)
               (vel-algo cl-mpm::velocity-algorithm))
                sim
    (declare (type double-float mass-filter))
    (with-slots ((mesh cl-mpm::mesh)
                 (mps cl-mpm::mps)
                 (bcs cl-mpm::bcs)
                 (bcs-force cl-mpm::bcs-force)
                 (dt cl-mpm::dt)
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
                 (vel-algo cl-mpm::velocity-algorithm))
        sim
      (unless initial-setup
        (pre-step sim))
      (setf dt 1d0)
      (cl-mpm::update-nodes sim)
      (cl-mpm::update-cells sim)
      (cl-mpm::reset-nodes-force sim)
      (cl-mpm::update-stress mesh mps dt-loadstep fbar)
      (cl-mpm/damage::calculate-damage sim dt-loadstep)
      (cl-mpm::p2g-force-fs sim)
      (cl-mpm::apply-bcs mesh bcs-force dt)
      (loop for bcs-f in bcs-force-list
            do (cl-mpm::apply-bcs mesh bcs-f dt))
      (update-node-fictious-mass sim)
      (when ghost-factor
        (cl-mpm/ghost::apply-ghost sim ghost-factor)
        (cl-mpm::apply-bcs mesh bcs dt))
      (setf damping (* damping-scale (cl-mpm/dynamic-relaxation::dr-estimate-damping sim)))
      ;; ;;Update our nodes after force mapping
      (cl-mpm::update-node-forces sim)
      (cl-mpm::apply-bcs mesh bcs dt)
      (cl-mpm::update-dynamic-stats sim)
      (cl-mpm::g2p mesh mps dt damping :TRIAL)
      )))

(defmethod cl-mpm::finalise-loadstep ((sim mpm-sim-dr-dynamic))
  ;;DR algorithm requires that finalisation is called once
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



(defparameter *total-iter* 0)
(defparameter *total-step* 0)
(defmethod cl-mpm::update-sim ((sim cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic))
  "Update stress last algorithm"
  (let ((crit (cl-mpm/dynamic-relaxation::sim-convergence-critera sim))
        (damage-enabled (cl-mpm::sim-enable-damage sim))
        (dt (cl-mpm::sim-dt sim))
        (dt-scale (cl-mpm::sim-dt-scale sim))
        )
    (setf (sim-dt-loadstep sim) (* 1d0 dt))
    (change-class sim 'cl-mpm/dynamic-relaxation::mpm-sim-dr-dynamic)
    ;; (setf (cl-mpm::sim-dt-scale sim) 0.25d0)
    (let ((prev-res nil)
          (res 0d0)
          (conv-crit 1d-6)
          (substeps 10)
          )
      (;generalised-staggered-solve
       converge-quasi-static
       sim
       ;; :enable-damage damage-enabled
       ;; :damping 1d0
       :energy-crit 1d0;crit
       :oobf-crit crit
       :substeps substeps
       :damping-factor 1d0
       :conv-steps 10000
       :dt-scale 1d0;dt-scale

       :convergance-criteria
       (lambda (sim f o)
         (let ((c (cl-mpm/dynamic-relaxation::res-norm-aggregated sim)))
           ;; (pprint c)
           (< c crit)))

       ;;   ;; (unless prev-res
       ;;   ;;   (setf prev-res (cl-mpm/dynamic-relaxation::res-norm-aggregated sim)))
       ;; :convergance-criteria
       ;; (lambda (sim f o)
       ;;   ;; (unless prev-res
       ;;   ;;   (setf prev-res (cl-mpm/dynamic-relaxation::res-norm-aggregated sim)))
       ;;   (setf res (cl-mpm/dynamic-relaxation::res-norm-aggregated sim))
       ;;   ;; (let* ((conv (if (> prev-res 0d0)
       ;;   ;;                  (< (/ res prev-res) conv-crit)
       ;;   ;;                  nil)))
       ;;   ;;   (format t "Succesion res criteria ~E~%"
       ;;   ;;           (if (> prev-res 0d0)
       ;;   ;;               (/ res prev-res)
       ;;   ;;               0d0))
       ;;     (setf (cl-mpm::sim-stats-oobf sim) res)
       ;;     (format t "Succesion res criteria ~E~%" res)
       ;;     (< res conv-crit)
       ;;   ;;   conv)
       ;;   )
       :post-iter-step (lambda (i e o)
                         (format t "Dynamic substep ~D~%" i)
                         ;; (when (uiop:directory-exists-p "./output/")
                         ;;   (cl-mpm/output:save-vtk (merge-pathnames "./output/" (format nil "rsim_step_~5,'0d.vtk" i)) sim)
                         ;;   (cl-mpm/output:save-vtk-nodes (merge-pathnames "./output/" (format nil "rsim_step_nodes_~5,'0d.vtk" i)) sim)
                         ;;   (save-conv-step sim "./output/" *total-iter* *total-step* 0d0 o e))
                         ;; (save-conv-step sim "./output/" *total-iter* *total-step* 0d0 o e)
                         (incf *total-iter* substeps))))
    (incf *total-step*)
    (setf (cl-mpm::sim-dt-scale sim) dt-scale)
    (change-class sim 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic)
    (setf (cl-mpm::sim-dt sim) dt)
    (cl-mpm::finalise-loadstep sim))
  )



(defmethod pre-step ((sim mpm-sim-dr-dynamic))
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
    (cl-mpm::filter-cells sim)
    (cl-mpm::update-node-kinematics sim)
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (n)
       (setf
        (cl-mpm/mesh::node-true-mass n) (cl-mpm/mesh:node-mass n))
       ;;More complicated things need to happen with aggregation
       (cl-mpm/utils:vector-copy-into (cl-mpm/mesh::node-velocity n) (cl-mpm/mesh::node-true-velocity n))
       ))
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

;; (defmethod cl-mpm::update-dynamic-stats ((sim mpm-sim-implict-dynamic))
;;   (with-accessors ((stats-energy cl-mpm::sim-stats-energy)
;;                    (stats-oobf cl-mpm::sim-stats-oobf)
;;                    (stats-power cl-mpm::sim-stats-power)
;;                    (stats-work cl-mpm::sim-stats-work))
;;       sim
;;     (multiple-value-bind (e o p) (cl-mpm/dynamic-relaxation::combi-stats sim)
;;       (setf stats-energy e
;;             stats-oobf o
;;             stats-power p)
;;       (incf stats-work p))))
