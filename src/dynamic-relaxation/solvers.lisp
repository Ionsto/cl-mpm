(in-package :cl-mpm/dynamic-relaxation)



(defmethod cl-mpm::finalise-loadstep ((sim mpm-sim-dr))
  ;;DR algorithm requires that finalisation is called once
  (incf (cl-mpm::sim-time sim) (sim-dt-loadstep sim))
  (cl-mpm::new-loadstep sim))

;; (defmethod cl-mpm::update-sim ((sim mpm-sim-dr-usf))
;;   "Update stress last algorithm"
;;   (declare (cl-mpm::mpm-sim sim))
;;   (with-slots ((mesh cl-mpm::mesh)
;;                (mps cl-mpm::mps)
;;                (bcs cl-mpm::bcs)
;;                (bcs-force cl-mpm::bcs-force)
;;                (dt cl-mpm::dt)
;;                (mass-filter cl-mpm::mass-filter)
;;                (split cl-mpm::allow-mp-split)
;;                (enable-damage cl-mpm::enable-damage)
;;                (nonlocal-damage cl-mpm::nonlocal-damage)
;;                (remove-damage cl-mpm::allow-mp-damage-removal)
;;                (fbar cl-mpm::enable-fbar)
;;                (bcs-force-list cl-mpm::bcs-force-list)
;;                (vel-algo cl-mpm::velocity-algorithm))
;;                 sim
;;     (declare (type double-float mass-filter))
;;                 (progn
;;                     (cl-mpm::reset-grid mesh)
;;                     (cl-mpm::p2g mesh mps)
;;                     (when (> mass-filter 0d0)
;;                       (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
;;                     (cl-mpm::update-node-kinematics sim)
;;                     (cl-mpm::apply-bcs mesh bcs dt)
;;                     (cl-mpm::update-nodes sim)
;;                     (cl-mpm::update-stress mesh mps dt fbar)
;;                     (cl-mpm::p2g-force-fs mesh mps)
;;                     (cl-mpm::apply-bcs mesh bcs-force dt)
;;                     (loop for bcs-f in bcs-force-list
;;                           do (cl-mpm::apply-bcs mesh bcs-f dt))
;;                     ;;Update our nodes after force mapping
;;                     (cl-mpm::update-node-forces sim)
;;                     (cl-mpm::apply-bcs mesh bcs dt)
;;                     ;; (cl-mpm::update-nodes sim)
;;                     (cl-mpm::update-dynamic-stats sim)
;;                     (cl-mpm::g2p mesh mps dt vel-algo)
;;                     ;; (when remove-damage
;;                     ;;   (cl-mpm::remove-material-damaged sim))
;;                     ;; (when split
;;                     ;;   (cl-mpm::split-mps sim))
;;                     ;; (cl-mpm::check-mps sim)
;;                     )))
(defmethod cl-mpm::finalise-loadstep ((sim mpm-sim-dr-ul))
  ;;DR algorithm requires that finalisation is called once
  (setf (sim-initial-setup sim) nil)
  ;; (cl-mpm::update-nodes sim)
  (cl-mpm::g2p (cl-mpm:sim-mesh sim)
               (cl-mpm:sim-mps sim)
               (cl-mpm:sim-dt sim)
               (cl-mpm:sim-damping-factor sim)
               (cl-mpm::sim-velocity-algorithm sim))
  (call-next-method)
  )
(defmethod cl-mpm::reset-loadstep ((sim mpm-sim-dr-ul))
  (setf (sim-initial-setup sim) nil)
  (call-next-method))

(defun pre-step (sim)
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
    (setf initial-setup t)))

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
      (unless initial-setup
        (pre-step sim))
      (setf dt 1d0)

      (cl-mpm::update-nodes sim)
      (cl-mpm::update-cells sim)
      (cl-mpm::reset-nodes-force sim)
      (cl-mpm::update-stress mesh mps dt-loadstep fbar)
      (cl-mpm::p2g-force-fs mesh mps)
      (cl-mpm::apply-bcs mesh bcs-force dt)
      (loop for bcs-f in bcs-force-list
            do (cl-mpm::apply-bcs mesh bcs-f dt))
      ;; ;;Update our nodes after force mapping
      (cl-mpm::update-node-forces sim)
      (cl-mpm::apply-bcs mesh bcs dt)
      (cl-mpm::update-dynamic-stats sim)
      ;; (cl-mpm::g2p mesh mps dt damping :QUASI-STATIC)
      )))


;; (defmethod cl-mpm::update-cells ((sim mpm-sim-dr-ul))
;;   (with-accessors ((mesh cl-mpm::sim-mesh)
;;                    (dt cl-mpm::sim-dt)
;;                    (agg sim-enable-aggregate))
;;       sim
;;     (cl-mpm::iterate-over-cells
;;      mesh
;;      (lambda (cell)
;;        (cl-mpm::filter-cell mesh cell dt)
;;        (cl-mpm::update-cell mesh cell dt)))
;;     (when agg
;;       (cl-mpm/aggregate::update-aggregate-elements sim))
;;     ))

(defun midpoint-starter (sim)
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (dt cl-mpm::dt)
               (fbar cl-mpm::enable-fbar)
               (agg cl-mpm/aggregate::enable-aggregate)
               (bcs-force-list cl-mpm::bcs-force-list))
                sim
    (progn
      (cl-mpm::update-stress mesh mps dt fbar)
      (cl-mpm::p2g-force mesh mps)
      (cl-mpm::apply-bcs mesh bcs-force dt)
      (loop for bcs-f in bcs-force-list
            do (cl-mpm::apply-bcs mesh bcs-f dt))
      (when agg
        (cl-mpm/aggregate::update-node-forces-agg sim (* -0.5d0 dt)))
      (cl-mpm::apply-bcs mesh bcs dt))
    ))

;; (defmethod cl-mpm::update-cells ((sim mpm-sim-dr-ul))
;;   (with-accessors ((mesh cl-mpm::sim-mesh)
;;                    (dt cl-mpm::sim-dt)
;;                    (agg cl-mpm::sim-enable-aggregate)
;;                    )
;;       sim
;;     (cl-mpm::iterate-over-cells
;;      mesh
;;      (lambda (cell)
;;        (cl-mpm::filter-cell mesh cell dt)
;;        (cl-mpm::update-cell mesh cell dt)))))


(defun map-stiffness (sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps)
                   )
      sim
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (setf (cl-mpm/mesh:node-mass node) 0d0))))
    (let ((mass-scale (the double-float (/ 1d0 (cl-mpm::sim-dt-scale sim)))))
      (cl-mpm::iterate-over-mps
       mps
       (lambda (mp)
         (with-accessors ((mp-volume cl-mpm/particle:mp-volume)
                          (mp-pmod cl-mpm/particle::mp-p-modulus)
                          (def cl-mpm/particle::mp-deformation-gradient)
                          )
             mp
           (let ((mp-volume mp-volume)
                 (mp-pmod mp-pmod))
             (declare (type double-float mp-pmod mp-volume))
             (cl-mpm::iterate-over-neighbours
              mesh mp
              (lambda (mesh mp node svp grads fsvp fgrads)
                (declare
                 (cl-mpm/particle:particle mp)
                 (cl-mpm/mesh::node node)
                 (double-float svp)
                 (ignore mesh))
                (with-accessors ((node-active  cl-mpm/mesh::node-active)
                                 (node-mass  cl-mpm/mesh:node-mass)
                                 (node-true-v  cl-mpm/mesh::node-volume-true)
                                 (node-lock  cl-mpm/mesh:node-lock))
                    node
                  (declare (type double-float mp-volume mp-pmod svp node-mass node-true-v)
                           (type sb-thread:mutex node-lock))
                  (when node-active
                    (sb-thread:with-mutex (node-lock)
                      (incf node-mass (/ (* 8
                                            mp-pmod
                                            (/ 1d0 (expt (cl-mpm/fastmaths:det def) 2))
                                            svp mp-volume
                                            mass-scale)
                                         node-true-v)))))))))))))
  )

(defun update-node-fictious-mass (sim)
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (dt cl-mpm::sim-dt))
      sim
    (map-stiffness sim)
    (setf dt 1d0)))


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
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (ghost-factor cl-mpm::ghost-factor)
               (initial-setup initial-setup)
               (enable-aggregate cl-mpm/aggregate::enable-aggregate)
               (damping cl-mpm::damping-factor)
               (vel-algo cl-mpm::velocity-algorithm))
      sim
    (unless initial-setup
      (pre-step sim)
      (setf (cl-mpm/damage::sim-damage-delocal-counter-max sim) -1)
      (cl-mpm/damage::update-delocalisation-list mesh mps))
    (setf dt 1d0)
    (cl-mpm::update-nodes sim)
    (cl-mpm::update-cells sim)
    (cl-mpm::reset-nodes-force sim)
    (cl-mpm::update-stress mesh mps dt-loadstep fbar)
    (cl-mpm/damage::calculate-damage sim dt-loadstep)
    (cl-mpm::p2g-force-fs mesh mps)
    (cl-mpm::apply-bcs mesh bcs-force dt)
    (loop for bcs-f in bcs-force-list
          do (cl-mpm::apply-bcs mesh bcs-f dt))
    (update-node-fictious-mass sim)
    ;; ;;Update our nodes after force mapping
    (cl-mpm::update-node-forces sim)
    (cl-mpm::apply-bcs mesh bcs dt)
    (cl-mpm::update-dynamic-stats sim)
    (cl-mpm::g2p mesh mps dt damping :QUASI-STATIC)
    ))




;; (defmethod cl-mpm::finalise-loadstep ((sim mpm-sim-dr-damage-ul))
;;   ;; (incf (cl-mpm::sim-time sim) (sim-dt-loadstep sim))
;;   ;; (cl-mpm::new-loadstep sim)
;;   (call-next-method))

;; (defmethod cl-mpm::update-sim ((sim mpm-sim-dr-damage-usf))
;;   "Update stress last algorithm"
;;   (declare (cl-mpm::mpm-sim sim))
;;   (with-slots ((mesh cl-mpm::mesh)
;;                (mps cl-mpm::mps)
;;                (bcs cl-mpm::bcs)
;;                (bcs-force cl-mpm::bcs-force)
;;                (dt cl-mpm::dt)
;;                (dt-loadstep dt-loadstep)
;;                (mass-filter cl-mpm::mass-filter)
;;                (split cl-mpm::allow-mp-split)
;;                (enable-damage cl-mpm::enable-damage)
;;                (nonlocal-damage cl-mpm::nonlocal-damage)
;;                (remove-damage cl-mpm::allow-mp-damage-removal)
;;                (fbar cl-mpm::enable-fbar)
;;                (bcs-force-list cl-mpm::bcs-force-list)
;;                (vel-algo cl-mpm::velocity-algorithm))
;;       sim
;;     (declare (type double-float mass-filter))
;;     (progn
;;       (cl-mpm::reset-grid mesh)
;;       (when (> (length mps) 0)
;;         (cl-mpm::p2g mesh mps)
;;         (when (> mass-filter 0d0)
;;           (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
;;         (cl-mpm::update-node-kinematics sim)
;;         ;; (update-node-fictious-mass sim)
;;         (cl-mpm::apply-bcs mesh bcs dt)
;;         (cl-mpm::update-nodes sim)
;;         (cl-mpm::update-stress mesh mps dt-loadstep fbar)
;;         (cl-mpm/damage::calculate-damage sim dt-loadstep)

;;         (cl-mpm::p2g-force mesh mps)
;;         ;; (cl-mpm::apply-bcs mesh bcs-force dt)
;;         (loop for bcs-f in bcs-force-list
;;               do (cl-mpm::apply-bcs mesh bcs-f dt))
;;         ;;Update our nodes after force mapping
;;         (cl-mpm::update-node-forces sim)
;;         (cl-mpm::apply-bcs mesh bcs dt)
;;         (cl-mpm::update-dynamic-stats sim)
;;         (cl-mpm::g2p mesh mps dt vel-algo)
;;         ;;
;;         (cl-mpm::update-particles sim)
;;         (cl-mpm::reset-node-displacement sim)

;;         (when remove-damage
;;           (cl-mpm::remove-material-damaged sim))
;;         (when split
;;           (cl-mpm::split-mps sim))
;;         (cl-mpm::check-mps sim))
;;       )))

;; ;; (defmethod cl-mpm::finalise-loadstep ((sim mpm-sim-dr-damage-usf))
;; ;;   (incf (cl-mpm::sim-time sim) (sim-dt-loadstep sim))
;; ;;   (cl-mpm:iterate-over-mps
;; ;;    (cl-mpm:sim-mps sim)
;; ;;    (lambda (mp)
;; ;;      (cl-mpm/particle::new-loadstep-mp mp))))


;; (in-package :cl-mpm/aggregate)

(defmethod cl-mpm::update-node-forces ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-ul))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mass-scale cl-mpm:sim-mass-scale)
                   (damping cl-mpm:sim-damping-factor)
                   (damping-algo cl-mpm::sim-damping-algorithm)
                   (agg-elems cl-mpm/aggregate::sim-agg-elems)
                   (dt cl-mpm:sim-dt)
                   (enable-aggregate cl-mpm/aggregate::sim-enable-aggregate))
      sim

    (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt)
    (cl-mpm:iterate-over-nodes
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (cl-mpm::calculate-forces node 0d0 0d0 mass-scale))))

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
                     (cl-mpm/aggregate::project-int-vec sim acc #'cl-mpm/mesh::node-acceleration d)
                     ))))

      (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt))
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
           (when (or internal
                     (not agg))
             (cl-mpm::integrate-vel-midpoint vel acc mass mass-scale dt damping))))))
    (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt)))

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
            (cl-mpm/fastmaths::fast-.+ inertia-force force-int force-int))
          ;;Set acc to f/m
          (cl-mpm/fastmaths::fast-.+-vector force-int force force)
          (cl-mpm/fastmaths::fast-.+-vector force-ext force force)
          ;; (cl-mpm/fastmaths::fast-.+-vector inertia-force force force)
          ;; (cl-mpm/utils:vector-copy-into inertia-force force-ghost)
          ;;Include velocity prop damping
          ;; (cl-mpm/fastmaths::fast-.+-vector force-damp force force)
          ;; (cl-mpm/fastmaths::fast-.+-vector force-ghost force force)
          (cl-mpm/fastmaths:fast-fmacc acc force (/ 1d0 (* mass mass-scale)))
          ;; (cl-mpm::integrate-vel-euler vel acc mass mass-scale dt 0d0)
          ;;Would intergrtate but don't
          (cl-mpm/utils::vector-copy-into residual residual-prev)
          (cl-mpm/utils::vector-copy-into force-int residual)))
      ))
  (values))



(defmethod cl-mpm::update-node-forces ((sim cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mass-scale cl-mpm::sim-mass-scale)
                   (damping cl-mpm::sim-damping-factor)
                   (damping-algo cl-mpm::sim-damping-algorithm)
                   (agg-elems cl-mpm/aggregate::sim-agg-elems)
                   (dt cl-mpm::sim-dt)
                   (dt-loadstep cl-mpm/dynamic-relaxation::sim-dt-loadstep)
                   (enable-aggregate cl-mpm/aggregate::sim-enable-aggregate))
      sim
    (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt)
    (cl-mpm:iterate-over-nodes
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (cl-mpm/dynamic-relaxation::dr-calculate-forces-implicit-dynamic node 0d0 0d0 mass-scale dt-loadstep))))

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

(defmethod cl-mpm::update-sim ((sim cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic))
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
      (unless initial-setup
        (pre-step sim))
      (setf dt 1d0)
      (cl-mpm::update-nodes sim)
      (cl-mpm::update-cells sim)
      (cl-mpm::reset-nodes-force sim)
      (cl-mpm::update-stress mesh mps dt-loadstep fbar)
      (cl-mpm::p2g-force-fs mesh mps)
      (cl-mpm::apply-bcs mesh bcs-force dt)
      (loop for bcs-f in bcs-force-list
            do (cl-mpm::apply-bcs mesh bcs-f dt))
      ;; ;;Update our nodes after force mapping
      (cl-mpm::update-node-forces sim)
      (cl-mpm::apply-bcs mesh bcs dt)
      (cl-mpm::update-dynamic-stats sim)
      ;; (cl-mpm::g2p mesh mps dt damping vel-algo)
      )))

(defmethod cl-mpm::finalise-loadstep ((sim mpm-sim-implict-dynamic))
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
       ;; (cl-mpm/utils:vector-copy-into (cl-mpm/mesh::node-velocity n) (cl-mpm/mesh::node-true-velocity n))
       )))
  (call-next-method))
