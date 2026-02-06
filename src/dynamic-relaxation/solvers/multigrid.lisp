(in-package :cl-mpm/dynamic-relaxation)
(defmethod %converge-quasi-static ((sim mpm-sim-dr-multigrid)
                                   energy-crit
                                   oobf-crit
                                   live-plot
                                   dt-scale
                                   substeps
                                   conv-steps
                                   post-iter-step
                                   convergance-criteria
                                   kinetic-damping
                                   damping-factor)
  (let ((refine (cl-mpm::sim-multigrid-refinement sim))
        (mass-filter (/ (cl-mpm::sim-mass-filter sim) (expt (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) 2))))

    ;; (setf (cl-mpm::sim-mass-filter sim) (* mass-filter  (expt (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) 2) ))
    ;; (setf (cl-mpm:sim-mesh sim) (nth 0 (cl-mpm::sim-mesh-list sim))
    ;;       (cl-mpm:sim-bcs sim)  (nth 0 (cl-mpm::sim-bcs-list sim)))
    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps sim)
     (lambda (mp)
       (setf (fill-pointer (cl-mpm/particle::mp-cached-nodes mp)) 0)))

    (let ((damage-state (cl-mpm::sim-enable-damage sim)))
      ;; (setf (cl-mpm::sim-enable-damage sim) nil)
      (let ((total-iters 0))
        (dotimes (i refine)
          (format t "Mesh step ~D~%" i)
          (setf (cl-mpm:sim-mesh sim) (nth i (cl-mpm::sim-mesh-list sim))
                (cl-mpm:sim-bcs sim)  (nth i (cl-mpm::sim-bcs-list sim)))
          (setf (cl-mpm::sim-multigrid-current-mesh sim) i)
          ;; (when (= (+ 1 i) refine)
          ;;   (setf (cl-mpm::sim-enable-damage sim) damage-state))
          (pre-step sim)
          (format t "Solve ~%")
          ;; (save-vtks-dr-step sim "./output/" 1000 (* 2 i))
          (call-next-method
           sim
           energy-crit
           oobf-crit
           live-plot
           dt-scale
           ;; (* substeps (expt 2 (- i 1)))
           substeps
           conv-steps
           (lambda (i o e)
             (funcall post-iter-step total-iters o e)
             (incf total-iters))
           convergance-criteria
           kinetic-damping
           damping-factor)
          ;; (save-vtks-dr-step sim "./output/" 1000 (+ 1 (* 2 i)))
          ;;Remap
          (unless (= (+ 1 i) refine)
            (setf (cl-mpm:sim-mesh sim) (nth (+ i 1) (cl-mpm::sim-mesh-list sim))
                  (cl-mpm:sim-bcs sim)  (nth (+ i 1) (cl-mpm::sim-bcs-list sim)))
            (cl-mpm::iterate-over-mps
             (cl-mpm:sim-mps sim)
             (lambda (mp)
               (setf (fill-pointer (cl-mpm/particle::mp-cached-nodes mp)) 0)))
            ;; (pre-step sim)
            (let ((mesh (cl-mpm:sim-mesh sim))
                  (mps (cl-mpm:sim-mps sim)))
              (cl-mpm::reset-grid mesh :reset-displacement nil)
              (cl-mpm::p2g mesh mps)
              (when (> mass-filter 0d0)
                (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim))))
            (let ((coarse-mesh (nth i (cl-mpm::sim-mesh-list sim))))
              (cl-mpm::iterate-over-nodes
               (cl-mpm:sim-mesh sim)
               (lambda (n)
                 (when (cl-mpm/mesh::node-active n)
                   (with-accessors ((disp cl-mpm/mesh::node-displacment))
                       n
                     (cl-mpm::fast-zero disp)
                     (cl-mpm::iterate-over-neighbours-point-linear
                      coarse-mesh
                      (cl-mpm/mesh::node-position n)
                      (lambda (mc coarse-node weight grads)
                        (cl-mpm/fastmaths::fast-fmacc
                         disp
                         (cl-mpm/mesh::node-displacment coarse-node)
                         weight))))))))))))))

(defmethod cl-mpm/damage::calculate-damage :around ((sim mpm-sim-dr-multigrid) dt)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps))
      sim
    (let ((m (cl-mpm:sim-mesh sim)))
      (setf mesh (first (last (cl-mpm::sim-mesh-list sim))))
      (call-next-method)
      (setf mesh m))))

(defmethod cl-mpm::sim-add-mp :around ((sim mpm-sim-dr-multigrid) mp)
  (with-accessors ((mesh cl-mpm::sim-mesh))
      sim
    (let ((m (cl-mpm:sim-mesh sim)))
      (setf mesh (first (last (cl-mpm::sim-mesh-list sim))))
      (call-next-method)
      (setf mesh m))))

;; (defmethod pre-step ((sim mpm-sim-dr-multigrid))
;;   (with-accessors ((mesh cl-mpm:sim-mesh)
;;                    (mps cl-mpm:sim-mps)
;;                    (delocal-counter cl-mpm/damage::sim-damage-delocal-counter-max))
;;       sim
;;     (call-next-method)
;;     (setf delocal-counter -1)
;;     (progn
;;       (cl-mpm/damage::update-delocalisation-list (first (last (cl-mpm::sim-mesh-list sim))) mps))))

(defmethod pre-step ((sim mpm-sim-dr-multigrid))
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
               (vel-algo cl-mpm::velocity-algorithm)
               (enable-damage cl-mpm::enable-damage)
               (current-mesh cl-mpm::current-mesh))
      sim
    (format t "Current mesh ~D~%" current-mesh)
    (if (= current-mesh 0)
        (progn
          (cl-mpm::reset-grid mesh :reset-displacement t)
          )
        (progn
          (cl-mpm::reset-grid mesh :reset-displacement nil)))
    (cl-mpm::p2g mesh mps)
    (when (> mass-filter 0d0)
      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (n)
       (setf
        (cl-mpm/mesh::node-true-mass n) (cl-mpm/mesh:node-mass n))
       (cl-mpm/fastmaths:fast-zero (cl-mpm/mesh::node-true-velocity n))))
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
    (setf initial-setup t)
    (setf (cl-mpm/damage::sim-damage-delocal-counter-max sim) -1)
    (progn
      (cl-mpm/damage::update-delocalisation-list (first (last (cl-mpm::sim-mesh-list sim))) mps))))
