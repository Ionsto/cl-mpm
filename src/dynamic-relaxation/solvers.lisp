(in-package :cl-mpm/dynamic-relaxation)



(defmethod cl-mpm::finalise-loadstep ((sim mpm-sim-dr))
  ;;DR algorithm requires that finalisation is called once
  (incf (cl-mpm::sim-time sim) (sim-dt-loadstep sim))
  (cl-mpm::new-loadstep sim))

(defmethod cl-mpm::update-sim ((sim mpm-sim-dr-usf))
  "Update stress last algorithm"
  (declare (cl-mpm::mpm-sim sim))
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (dt cl-mpm::dt)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (nonlocal-damage cl-mpm::nonlocal-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (vel-algo cl-mpm::velocity-algorithm))
                sim
    (declare (type double-float mass-filter))
                (progn
                    (cl-mpm::reset-grid mesh)
                    (cl-mpm::p2g mesh mps)
                    (when (> mass-filter 0d0)
                      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
                    (cl-mpm::update-node-kinematics sim)
                    (cl-mpm::apply-bcs mesh bcs dt)
                    (cl-mpm::update-nodes sim)
                    (cl-mpm::update-stress mesh mps dt fbar)
                    (cl-mpm::p2g-force-fs mesh mps)
                    (cl-mpm::apply-bcs mesh bcs-force dt)
                    (loop for bcs-f in bcs-force-list
                          do (cl-mpm::apply-bcs mesh bcs-f dt))
                    ;;Update our nodes after force mapping
                    (cl-mpm::update-node-forces sim)
                    (cl-mpm::apply-bcs mesh bcs dt)
                    ;; (cl-mpm::update-nodes sim)
                    (cl-mpm::update-dynamic-stats sim)
                    (cl-mpm::g2p mesh mps dt vel-algo)
                    ;; (when remove-damage
                    ;;   (cl-mpm::remove-material-damaged sim))
                    ;; (when split
                    ;;   (cl-mpm::split-mps sim))
                    ;; (cl-mpm::check-mps sim)
                    )))
(defmethod cl-mpm::finalise-loadstep ((sim mpm-sim-dr-ul))
  ;;DR algorithm requires that finalisation is called once
  (setf (sim-initial-setup sim) nil)

  (cl-mpm::update-nodes sim)
  (cl-mpm::g2p (cl-mpm:sim-mesh sim)
               (cl-mpm:sim-mps sim)
               (cl-mpm:sim-dt sim)
               (cl-mpm:sim-damping-factor sim)
               (cl-mpm::sim-velocity-algorithm sim))
  ;; (incf (cl-mpm::sim-time sim) (sim-dt-loadstep sim))
  ;; (cl-mpm::new-loadstep sim)
  (call-next-method)
  )
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
    (progn
      (unless initial-setup
        (cl-mpm::reset-grid mesh)
        (cl-mpm::reset-node-displacement sim)
        (cl-mpm::p2g mesh mps)
        (when (> mass-filter 0d0)
          (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
        (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
        (update-node-fictious-mass sim)
        ;; (set-mass sim)
        ;; (setf (cl-mpm:sim-dt sim) (cl-mpm::calculate-min-dt sim))
        ;; (setf (cl-mpm:sim-dt sim) (* (cl-mpm::sim-dt-scale sim) (cl-mpm::calculate-min-dt sim)))
        ;; (setf (cl-mpm:sim-dt sim) (* (cl-mpm::sim-dt-scale sim) (cl-mpm/setup::estimate-elastic-dt sim)))
        ;; (setf (cl-mpm:sim-dt sim) (* (cl-mpm::sim-dt-scale sim) (cl-mpm::calculate-min-dt sim)))
        (format t "Min dt~E~%" (cl-mpm:sim-dt sim))
        (cl-mpm::apply-bcs mesh bcs dt)
        (cl-mpm::update-cells sim)
        (when enable-aggregate
          (cl-mpm/aggregate::update-aggregate-elements sim))
        (cl-mpm::apply-bcs mesh bcs dt)
        (midpoint-starter sim)
       (setf initial-setup t))
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
      (cl-mpm::g2p mesh mps dt damping :QUASI-STATIC))))


(defmethod cl-mpm::update-cells ((sim mpm-sim-dr-ul))
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (dt cl-mpm::sim-dt)
                   (agg sim-enable-aggregate))
      sim
    (cl-mpm::iterate-over-cells
     mesh
     (lambda (cell)
       (cl-mpm::filter-cell mesh cell dt)
       (cl-mpm::update-cell mesh cell dt)))
    ;; (when agg
    ;;   (cl-mpm/aggregate::update-aggregate-elements sim))
    ))

(defun midpoint-starter (sim)
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (dt cl-mpm::dt)
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list))
                sim
    (progn
      (cl-mpm::update-stress mesh mps dt fbar)
      (cl-mpm::p2g-force mesh mps)
      (cl-mpm::apply-bcs mesh bcs-force dt)
      (loop for bcs-f in bcs-force-list
            do (cl-mpm::apply-bcs mesh bcs-f dt))
      (cl-mpm/aggregate::update-node-forces-agg sim (* -0.5d0 dt))
      (cl-mpm::apply-bcs mesh bcs dt))
    ))

(defmethod cl-mpm::update-cells ((sim mpm-sim-dr-ul))
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (dt cl-mpm::sim-dt)
                   (agg cl-mpm::sim-enable-aggregate)
                   )
      sim
    (cl-mpm::iterate-over-cells
     mesh
     (lambda (cell)
       (cl-mpm::filter-cell mesh cell dt)
       (cl-mpm::update-cell mesh cell dt)))))




(defun set-mass (sim)
  (let ((density 1d0))
    (cl-mpm::iterate-over-nodes
     (cl-mpm:sim-mesh sim)
     (lambda (n)
       (when (cl-mpm/mesh::node-active n)
         (setf (cl-mpm/mesh::node-mass n)
               1d0
               ;; (* 100d0 (cl-mpm/mesh::node-mass n))
               ;; (* 1d0 (cl-mpm/mesh::node-volume n))
               ;; (* alpha (/(cl-mpm/mesh::node-pwave n) (cl-mpm/mesh::node-svp-sum n)))
               )))))
  ;; (let ((max-mass 0d0))
  ;;   (setf max-mass
  ;;         (cl-mpm::reduce-over-nodes
  ;;          (cl-mpm:sim-mesh sim)
  ;;          (lambda (n)
  ;;            (if (and
  ;;                 (cl-mpm/mesh:node-active n))
  ;;                ;; (cl-mpm/mesh:node-mass n)
  ;;                ;; (/ (cl-mpm/mesh::node-pwave n) (cl-mpm/mesh::node-pwave n))
  ;;                sb-ext:double-float-negative-infinity))
  ;;          #'max))
  ;;   (cl-mpm:iterate-over-nodes
  ;;    (cl-mpm:sim-mesh sim)
  ;;    (lambda (n)
  ;;      (when (cl-mpm/mesh:node-active n)
  ;;        (setf (cl-mpm/mesh:node-mass n) max-mass)))))
  )

(defun update-node-fictious-mass (sim)
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (dt cl-mpm::sim-dt))
      sim
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (setf (cl-mpm/mesh:node-mass node)
               (* (/ 1d0 (cl-mpm::sim-dt-scale sim))
                  (cl-mpm/mesh::node-pwave node))))))
    (setf dt 1d0)
    )
  )


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
               (damaping cl-mpm::damping-factor)
               (vel-algo cl-mpm::velocity-algorithm))
      sim
    (declare (type double-float mass-filter))
    (progn
      (unless initial-setup
        (cl-mpm::reset-grid mesh)
        (cl-mpm::p2g mesh mps)
        (when (> mass-filter 0d0)
          (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
        (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
        (update-node-fictious-mass sim)
        ;; (setf (cl-mpm:sim-dt sim) (cl-mpm::calculate-min-dt sim))
        ;; (setf (cl-mpm:sim-dt sim) (* (cl-mpm::sim-dt-scale sim) (cl-mpm::calculate-min-dt sim)))
        ;; (setf (cl-mpm:sim-dt sim) (* (cl-mpm::sim-dt-scale sim) (cl-mpm/setup::estimate-elastic-dt sim)))
        ;; (setf (cl-mpm:sim-dt sim) (* (cl-mpm::sim-dt-scale sim) (cl-mpm::calculate-min-dt sim)))
        (format t "Min dt~E~%" (cl-mpm:sim-dt sim))
        (cl-mpm::apply-bcs mesh bcs dt)
        (cl-mpm::update-cells sim)
        (when enable-aggregate
          (cl-mpm/aggregate::update-aggregate-elements sim))
        (cl-mpm::apply-bcs mesh bcs dt)
        (midpoint-starter sim)
        (setf (cl-mpm/damage::sim-damage-delocal-counter sim) 0)
        (cl-mpm/damage::calculate-damage sim 0d0)
        (setf (cl-mpm/damage::sim-damage-delocal-counter sim) -1)
        (setf initial-setup t))
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
      ;; ;;Update our nodes after force mapping
      (cl-mpm::update-node-forces sim)
      (cl-mpm::apply-bcs mesh bcs dt)
      (cl-mpm::update-dynamic-stats sim)
      (cl-mpm::g2p mesh mps dt damping vel-algo)
      ))
  )




(defmethod cl-mpm::finalise-loadstep ((sim mpm-sim-dr-damage-ul))
  ;; (incf (cl-mpm::sim-time sim) (sim-dt-loadstep sim))
  ;; (cl-mpm::new-loadstep sim)
  (call-next-method))

(defmethod cl-mpm::update-sim ((sim mpm-sim-dr-damage-usf))
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
               (vel-algo cl-mpm::velocity-algorithm))
      sim
    (declare (type double-float mass-filter))
    (progn
      (cl-mpm::reset-grid mesh)
      (when (> (length mps) 0)
        (cl-mpm::p2g mesh mps)
        (when (> mass-filter 0d0)
          (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
        (cl-mpm::update-node-kinematics sim)
        ;; (update-node-fictious-mass sim)
        (cl-mpm::apply-bcs mesh bcs dt)
        (cl-mpm::update-nodes sim)
        (cl-mpm::update-stress mesh mps dt-loadstep fbar)
        (cl-mpm/damage::calculate-damage sim dt-loadstep)

        (cl-mpm::p2g-force mesh mps)
        ;; (cl-mpm::apply-bcs mesh bcs-force dt)
        (loop for bcs-f in bcs-force-list
              do (cl-mpm::apply-bcs mesh bcs-f dt))
        ;;Update our nodes after force mapping
        (cl-mpm::update-node-forces sim)
        (cl-mpm::apply-bcs mesh bcs dt)
        (cl-mpm::update-dynamic-stats sim)
        (cl-mpm::g2p mesh mps dt vel-algo)
        ;;
        (cl-mpm::update-particles sim)
        (cl-mpm::reset-node-displacement sim)

        (when remove-damage
          (cl-mpm::remove-material-damaged sim))
        (when split
          (cl-mpm::split-mps sim))
        (cl-mpm::check-mps sim))
      )))

;; (defmethod cl-mpm::finalise-loadstep ((sim mpm-sim-dr-damage-usf))
;;   (incf (cl-mpm::sim-time sim) (sim-dt-loadstep sim))
;;   (cl-mpm:iterate-over-mps
;;    (cl-mpm:sim-mps sim)
;;    (lambda (mp)
;;      (cl-mpm/particle::new-loadstep-mp mp))))


(in-package :cl-mpm/aggregate)
(defmethod cl-mpm::update-node-forces ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-ul))
  ;;For non-aggregate nodes, use simple mass matrix inversion
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mass-scale sim-mass-scale)
                   (damping sim-damping-factor)
                   (damping-algo sim-damping-algorithm)
                   (agg-elems sim-agg-elems)
                   (enable-aggregate sim-enable-aggregate)
                   (dt sim-dt))
      sim
    (declare (double-float dt damping))

    (iterate-over-nodes
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (cl-mpm::calculate-forces-midpoint node damping 0d0 mass-scale))))

    (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt)

    ;;For each aggregated element set solve mass matrix and velocity
    (when enable-aggregate
      ;; (iterate-over-agg-elem
      ;;  agg-elems
      ;;  (lambda (elem)
      ;;    (calculate-forces-agg-elem sim elem damping)))
      (let* (;; (E (cl-mpm/aggregate::assemble-global-e sim))
             ;; (mii (cl-mpm/aggregate::assemble-global-mass sim))
             (E (sim-global-e sim))
             (ma (sim-global-ma sim))
             (f (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-force))
             ;; (fa (magicl:@ E (magicl:transpose E) f))
             )

        ;; (cl-mpm/aggregate::project-global-vec sim fa cl-mpm/mesh::node-force)

        ;; (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt)
        ;; (iterate-over-nodes
        ;;  mesh
        ;;  (lambda (node)
        ;;    (with-accessors ((agg cl-mpm/mesh::node-agg)
        ;;                     (active node-active)
        ;;                     (force node-force))
        ;;        node
        ;;      (when (and active agg)
        ;;        (fast-zero force)))))

        (let* (;(f (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-force))
               (acc
                 (magicl:@
                  E
                  ;; (magicl:inv (magicl:@ (magicl:transpose E) mii E))
                  (magicl:linear-solve
                   ma
                   (magicl:@
                    (magicl:transpose E)
                    f)))))
          (cl-mpm/aggregate::project-global-vec sim acc cl-mpm/mesh::node-acceleration))

        (cl-mpm/aggregate::project-global-vec sim (magicl:@ E (magicl:transpose E) (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-internal-force)) cl-mpm/mesh::node-internal-force)
        (cl-mpm/aggregate::project-global-vec sim (magicl:@ E (magicl:transpose E) (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-external-force)) cl-mpm/mesh::node-external-force)
        (cl-mpm/aggregate::project-global-vec sim (magicl:@ E (magicl:transpose E) (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-residual)) cl-mpm/mesh::node-residual)
        )
      )

    (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt)

    ;; (setf (cl-mpm:sim-damping-factor sim) (cl-mpm/dynamic-relaxation::dr-estimate-damping sim))
    (iterate-over-nodes
     mesh
     (lambda (node)
       (when (and (cl-mpm/mesh:node-active node))
         (with-accessors ((mass node-mass)
                          (vel node-velocity)
                          (acc node-acceleration))
             node
           (cl-mpm::integrate-vel-midpoint vel acc mass mass-scale dt damping)))))

    (when enable-aggregate
      (let ((E (sim-global-e sim))
            (ma (sim-global-ma sim))
            (disp (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-displacment))
            (disp-proj (cl-mpm/aggregate::assemble-global-internal-vec sim #'cl-mpm/mesh::node-displacment))
            )
        (cl-mpm/aggregate::project-global-vec sim (magicl:@ E disp-proj) cl-mpm/mesh::node-displacment))
      (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt))
    ))
