(in-package :cl-mpm/dynamic-relaxation)
(defun map-stiffness (sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps))
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
                      (incf node-mass (/ (* 2
                                            mp-pmod
                                            (/ 1d0 (expt
                                                    (min (cl-mpm/utils::mtref def 0 0)
                                                         (cl-mpm/utils::mtref def 1 1))
                                                    ;; (cl-mpm/fastmaths:det def)
                                                    2))
                                            svp mp-volume
                                            mass-scale)
                                         node-true-v))))))))))))
    )
  )
(defun update-node-fictious-mass (sim)
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (dt cl-mpm::sim-dt))
      sim
    (map-stiffness sim)
    (cl-mpm/penalty::assemble-penalty-stiffness-matrix sim)
    (cl-mpm/aggregate::update-mass-matrix sim)
    (setf dt 1d0)))

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
         (cl-mpm::p2g-force sim)
         (cl-mpm::apply-bcs mesh bcs-force dt)
         (loop for bcs-f in bcs-force-list
               do (cl-mpm::apply-bcs mesh bcs-f dt))
         (when agg
           (cl-mpm/aggregate::update-node-forces-agg sim (* -0.5d0 dt)))
         (cl-mpm::apply-bcs mesh bcs dt))
       ))
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
    (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
    (setf initial-setup t)))

(defmethod cl-mpm::update-sim ((sim mpm-sim-dr-ul))
  "Update stress last algorithm"
  (declare (cl-mpm::mpm-sim sim))
    ;; (declare (type double-float mass-filter))
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
      (unless initial-setup
        (pre-step sim)
        (with-accessors ((ke-prev sim-ke-prev)
                         (ke sim-ke))
            sim
          (setf ke 0d0
                ke-prev 0d0)))
      (setf dt 1d0)
      (cl-mpm/penalty::reset-penalty sim)
      (cl-mpm::update-nodes sim)
      (cl-mpm::update-cells sim)
      (cl-mpm::reset-nodes-force sim)
      (cl-mpm::update-stress mesh mps dt-loadstep fbar)
      (cl-mpm::p2g-force-fs sim)
      (cl-mpm::apply-bcs mesh bcs-force dt)
      (loop for bcs-f in bcs-force-list
            do (cl-mpm::apply-bcs mesh bcs-f dt))
      (update-node-fictious-mass sim)
      (setf damping (* damping-scale (cl-mpm/dynamic-relaxation::dr-estimate-damping sim)))
      ;; ;;Update our nodes after force mapping
      (cl-mpm::update-node-forces sim)
      (cl-mpm::apply-bcs mesh bcs dt)
      ;; (with-accessors ((ke-prev sim-ke-prev)
      ;;                  (ke sim-ke))
      ;;     sim
      ;;   (setf ke (calculate-ke sim))
      ;;   ;; (setf ke (cl-mpm::sim-stats-energy sim))
      ;;   (when (> ke-prev ke)
      ;;     (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
      ;;     (cl-mpm::reset-nodes-force sim)
      ;;     (setf ke 0d0))
      ;;   (setf ke-prev ke))
      (cl-mpm::update-dynamic-stats sim)
      ;; (cl-mpm::g2p mesh mps dt damping :QUASI-STATIC)
      ))

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
    (cl-mpm/penalty::reset-penalty sim)
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
    ;; ;;Update our nodes after force mapping
    (cl-mpm::update-node-forces sim)
    (cl-mpm::apply-bcs mesh bcs dt)
    (cl-mpm::update-dynamic-stats sim)
    ;; (cl-mpm::g2p mesh mps dt damping :QUASI-STATIC)
    ))
(defmethod cl-mpm::update-node-forces ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-ul))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mass-scale cl-mpm:sim-mass-scale)
                   (damping cl-mpm:sim-damping-factor)
                   (damping-scale cl-mpm/dynamic-relaxation::sim-damping-scale)
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

    (setf damping (* damping-scale (cl-mpm/dynamic-relaxation::dr-estimate-damping sim)))
    (with-accessors ((ke-prev sim-ke-prev)
                     (ke sim-ke)
                     (ke-damping sim-kinetic-damping))
        sim
      (when ke-damping
        (setf ke (calculate-ke sim))
        ;; (setf ke (cl-mpm::sim-stats-energy sim))
        (when (> ke-prev ke)
          (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
          (cl-mpm::reset-nodes-force sim)
          (setf ke 0d0))
        (setf ke-prev ke)))
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





;;Funny switch to enable 1st order richardson (very slow)
;; (defmethod cl-mpm::update-nodes ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-ul))
;;   (with-accessors ((mesh cl-mpm::sim-mesh)
;;                    (dt cl-mpm::sim-dt)
;;                    (agg cl-mpm/aggregate::sim-enable-aggregate))
;;       sim
;;     (cl-mpm::iterate-over-nodes
;;      mesh
;;      (lambda (node)
;;        (when (and (cl-mpm/mesh:node-active node)
;;                   (or (not (cl-mpm/mesh::node-agg node))
;;                       (cl-mpm/mesh::node-interior node)))
;;          (if t
;;              (cl-mpm::update-node node dt)
;;              (with-accessors ((vel   cl-mpm::node-acceleration)
;;                               (disp   cl-mpm/mesh::node-displacment))
;;                  node
;;                (cl-mpm/fastmaths:fast-fmacc disp vel dt)))
;;          )))
;;     (when agg
;;       (loop for d from 0 below (cl-mpm/mesh::mesh-nd mesh)
;;             do
;;                (let ((E (cl-mpm/aggregate::sim-global-e sim))
;;                      (disp-proj (cl-mpm/aggregate::assemble-internal-vec sim #'cl-mpm/mesh::node-displacment d)))
;;                  (cl-mpm/aggregate::apply-internal-bcs sim disp-proj d)
;;                  (cl-mpm/aggregate::project-global-vec
;;                   sim
;;                   (magicl:@ E disp-proj)
;;                   #'cl-mpm/mesh::node-displacment
;;                   d)))
;;       (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt))
;;     ))

(defmethod cl-mpm::update-sim ((sim mpm-sim-dr-usf))
  "Update stress first algorithm"
  (declare (cl-mpm::mpm-sim-usf sim))
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
               (time cl-mpm::time)
               (vel-algo cl-mpm::velocity-algorithm))
                sim
    (progn
      (cl-mpm::reset-grid mesh)
      (when (> (length mps) 0)
        (cl-mpm::p2g mesh mps)

        (when (> mass-filter 0d0)
          (cl-mpm::filter-grid mesh mass-filter))

        (cl-mpm::filter-cells sim)

        (cl-mpm::update-node-kinematics sim)

        (cl-mpm::apply-bcs mesh bcs dt)
        ;;Trial update displacements
        (cl-mpm::update-nodes sim)
        (cl-mpm::update-cells sim)

        (cl-mpm::update-stress mesh mps dt fbar)
        (update-node-fictious-mass sim)
        (setf (cl-mpm:sim-dt sim) 1d0)
        ;; Map forces onto nodes
        (cl-mpm::p2g-force sim)
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
        (cl-mpm::update-particles sim)
        ;; (cl-mpm::new-loadstep sim)
        )
      (incf time dt))))
