(in-package :cl-mpm/dynamic-relaxation)


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



