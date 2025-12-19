(in-package :cl-mpm/dynamic-relaxation)
(defmethod map-stiffness ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-ul))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps))
      sim
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (setf (cl-mpm/mesh:node-mass node) 0d0))))
    (let ((h (cl-mpm/mesh::mesh-resolution mesh))
          (nd (cl-mpm/mesh::mesh-nd mesh))
          (mass-scale (the double-float (/ 1d0 (cl-mpm::sim-dt-scale sim)))))
      (cl-mpm::iterate-over-mps
       mps
       (lambda (mp)
         (with-accessors ((mp-volume cl-mpm/particle::mp-volume)
                          (mp-volume-n cl-mpm/particle::mp-volume-n)
                          (mp-pmod cl-mpm/particle::mp-p-modulus)
                          (def-n cl-mpm/particle::mp-deformation-gradient-0)
                          (def cl-mpm/particle::mp-deformation-gradient)
                          (df cl-mpm/particle::mp-deformation-gradient-increment)
                          )
             mp
           (let ((mp-volume mp-volume)
                 (mp-pmod mp-pmod)
                 (ul (estimate-ul-enhancement mp)))
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
                      (let* (;; (grads (cl-mpm::gradient-push-forwards grads df))
                             ;; (peak-grad (reduce #'max (mapcar #'abs grads)))
                             ;; (peak-grad (reduce #'max (mapcar (lambda (x) (/ 1d0 x)) (remove 0d0 grads))))
                             )
                        (declare (double-float ul h))
                        ;; (format t "~A ~%" peak-grad)
                        ;; (incf node-mass (* 1
                        ;;                    (expt (cl-mpm/fastmaths::det-3x3 def) -1)
                        ;;                    mp-pmod
                        ;;                    svp
                        ;;                    mp-volume-n
                        ;;                    (expt h -2)
                        ;;                    mass-scale))
                        (incf node-mass (* 2
                                           (expt (cl-mpm/fastmaths::det-3x3 def) -1)
                                           mp-pmod
                                           svp
                                           mp-volume
                                           ul
                                           (expt h -2)
                                           mass-scale)))
                      ;; (incf node-mass (/ (* 2
                      ;;                       mp-pmod
                      ;;                       svp
                      ;;                       mp-volume
                      ;;                       mass-scale)
                      ;;                    ;; node-true-v
                      ;;                    (expt                           h 2)
                      ;;                    ))
                      ))))))))))))

(defgeneric update-node-fictious-mass (sim))

(defmethod update-node-fictious-mass ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr))
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (dt cl-mpm::sim-dt)
                   (bcs-force-list cl-mpm::sim-bcs-force-list))
      sim
    (map-stiffness sim)
    (cl-mpm/penalty::assemble-penalty-stiffness-matrix sim)
    (loop for bcs-f in bcs-force-list
          do (loop for bc across bcs-f
                   do (cl-mpm/bc::assemble-bc-stiffness sim bc)))
    (cl-mpm/aggregate::update-mass-matrix sim)
    (setf dt 1d0)))

(defmethod cl-mpm::finalise-loadstep ((sim mpm-sim-dr-ul))
  ;;DR algorithm requires that finalisation is called once
  (setf (sim-initial-setup sim) nil)
  ;; (cl-mpm::update-nodes sim)
  ;; (pprint (cl-mpm::sim-velocity-algorithm sim))
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
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (ghost-factor cl-mpm::ghost-factor)
               (initial-setup initial-setup)
               (enable-aggregate cl-mpm/aggregate::enable-aggregate)
               (damping cl-mpm::damping-factor)
               (vel-algo cl-mpm::velocity-algorithm))
      sim

    (cl-mpm::reset-grid mesh :reset-displacement t)
    (cl-mpm::reset-node-displacement sim)
    (cl-mpm::p2g mesh mps)
    (when (> mass-filter 0d0)
      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
    (cl-mpm::filter-cells sim)
    ;; (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
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
    (when ghost-factor
      (cl-mpm/ghost::apply-ghost sim ghost-factor)
      (cl-mpm::apply-bcs mesh bcs dt))
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
    (cl-mpm::g2p mesh mps dt damping :TRIAL)
    (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)))





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
    (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
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
         (cl-mpm::calculate-forces-midpoint node 0d0 0d0 mass-scale))))

    (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt)
    ;;For each aggregated element set solve mass matrix and velocity
    (when enable-aggregate
      (let* ((E (cl-mpm/aggregate::sim-global-e sim))
             (ma (cl-mpm/aggregate::sim-global-ma sim)))
        (loop for d from 0 below (cl-mpm::mesh-nd mesh)
              do (let ((f (magicl:@ (magicl:transpose E)
                                    (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-force d))))
                   (cl-mpm/aggregate::apply-internal-bcs sim f d)
                   (let* ((acc
                            (cl-mpm/aggregate::linear-solve-with-bcs
                             ma f (cl-mpm/aggregate::assemble-internal-bcs sim d))
                            ;; (magicl:linear-solve ma f)
                               ))
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
        (loop for d from 0 below (cl-mpm/mesh::mesh-nd mesh)
              do (let* ((f-int (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-internal-force d))
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
                   )))
      (let ((oobf 0d0)
            (oobf-num (sqrt oobf-num))
            (oobf-denom (sqrt oobf-denom))
            )
        (if (> oobf-denom 0d0)
            (setf oobf (/ oobf-num oobf-denom))
            (setf oobf (if (> oobf-num 0d0) sb-ext:double-float-positive-infinity 0d0)))
        oobf))))
