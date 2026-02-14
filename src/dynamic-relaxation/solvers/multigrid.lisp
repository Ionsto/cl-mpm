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
        (conv-crit (cl-mpm/dynamic-relaxation::sim-convergence-critera sim))
        (mass-filter (/ (cl-mpm::sim-mass-filter sim)
                        (expt (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) 2))))

    ;; (setf (cl-mpm::sim-mass-filter sim) (* mass-filter  (expt (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) 2) ))
    ;; (setf (cl-mpm:sim-mesh sim) (nth 0 (cl-mpm::sim-mesh-list sim))
    ;;       (cl-mpm:sim-bcs sim)  (nth 0 (cl-mpm::sim-bcs-list sim)))
    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps sim)
     (lambda (mp)
       (setf (fill-pointer (cl-mpm/particle::mp-cached-nodes mp)) 0)))

    (format t "Multigrid solve~%")
    (let ((damage-state (cl-mpm::sim-enable-damage sim)))
      ;; (setf (cl-mpm::sim-enable-damage sim) nil)
      (let ((total-iters 0))
        (dotimes (i (+ refine 1))
          (let ((final-step (= i refine)))
            (format t "Mesh step ~D~%" i)
            (setf (cl-mpm:sim-mesh sim) (nth i (cl-mpm::sim-mesh-list sim))
                  (cl-mpm:sim-bcs sim)  (nth i (cl-mpm::sim-bcs-list sim)))
            (setf (cl-mpm::sim-multigrid-current-mesh sim) i)
            ;; (when (= (+ 1 i) refine)
            ;;   (setf (cl-mpm::sim-enable-damage sim) damage-state))
            (pre-step sim)
            (format t "Solve final step ~A~%" final-step)
            ;; (save-vtks-dr-step sim "./output/" 1000 (* 2 i))
            (setf (cl-mpm/dynamic-relaxation::sim-convergence-critera sim)
                  (if final-step conv-crit (sqrt conv-crit)))
            (format t "Conv crit ~E~%" (cl-mpm/dynamic-relaxation::sim-convergence-critera sim))
            (call-next-method
             sim
             (cl-mpm/dynamic-relaxation::sim-convergence-critera sim)
             (cl-mpm/dynamic-relaxation::sim-convergence-critera sim)
             ;; (if final-step energy-crit (sqrt energy-crit))
             ;; (if final-step oobf-crit   (sqrt oobf-crit))
             live-plot
             dt-scale
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
            (unless final-step
              (setf (cl-mpm:sim-mesh sim) (nth (+ i 1) (cl-mpm::sim-mesh-list sim))
                    (cl-mpm:sim-bcs sim)  (nth (+ i 1) (cl-mpm::sim-bcs-list sim)))
              (cl-mpm::iterate-over-mps
               (cl-mpm:sim-mps sim)
               (lambda (mp)
                 (setf (fill-pointer (cl-mpm/particle::mp-cached-nodes mp)) 0)))

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
                           weight)))))))))))))))

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
               (current-mesh cl-mpm::current-mesh))
      sim
    (format t "Current mesh ~D~%" current-mesh)
    (if (= current-mesh 0)
        (progn
          (cl-mpm::reset-grid mesh :reset-displacement t))
        (progn
          (cl-mpm::reset-grid mesh :reset-displacement nil)))
    (setf (cl-mpm/dynamic-relaxation::sim-solve-count sim) 0)
    (cl-mpm::p2g mesh mps)
    (when (> mass-filter 0d0)
      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
    (cl-mpm::apply-bcs mesh bcs dt)
    (cl-mpm::filter-cells sim)
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (n)
       (setf
        (cl-mpm/mesh::node-true-mass n) (cl-mpm/mesh:node-mass n))
       (cl-mpm/fastmaths:fast-zero (cl-mpm/mesh::node-true-velocity n))))
    (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
    (cl-mpm::update-cells sim)
    (cl-mpm::apply-bcs mesh bcs dt)
    (midpoint-starter sim)
    (setf initial-setup t)
    (setf (cl-mpm/damage::sim-damage-delocal-counter-max sim) -1)
    (when (= current-mesh 0)
      (cl-mpm/damage::update-delocalisation-list (first (last (cl-mpm::sim-mesh-list sim))) mps))))


(defmethod cl-mpm/setup::%post-make-simple-sim ((sim cl-mpm::mpm-sim-multigrid) resolution element-count args-list)
  (let* ((size (mapcar (lambda (x) (* x resolution)) element-count)))
    (progn
      (setf (cl-mpm::sim-mesh-list sim) (list)
            (cl-mpm::sim-bcs-list sim) (list))
      (dotimes (refine (+ (cl-mpm::sim-multigrid-refinement sim) 1))
        (let* ((m (cl-mpm::make-mesh size (/ resolution (expt 2 (+ refine 0))) nil))
               (bcs (cl-mpm/bc:make-outside-bc m)))
          (push m (cl-mpm::sim-mesh-list sim))
          (push bcs (cl-mpm::sim-bcs-list sim))))
      (setf (cl-mpm::sim-mesh-list sim) (reverse (cl-mpm::sim-mesh-list sim))
            (cl-mpm::sim-bcs-list sim) (reverse (cl-mpm::sim-bcs-list sim)))
      (setf (cl-mpm:sim-mesh sim) (first (last (cl-mpm::sim-mesh-list sim))))
      (setf (cl-mpm:sim-bcs sim)  (first (last (cl-mpm::sim-bcs-list sim))))
      sim)))

(defmethod cl-mpm/setup::%setup-bcs ((sim cl-mpm::mpm-sim-multigrid)
                       left
                       right
                       top
                       bottom
                       front
                       back)
  (loop for mesh in (cl-mpm::sim-mesh-list sim)
        for i from 0
        do
           (setf (nth i (cl-mpm::sim-bcs-list sim))
                 (cl-mpm/bc::make-outside-bc-varfix
                  mesh
                  left right top bottom front back))))
