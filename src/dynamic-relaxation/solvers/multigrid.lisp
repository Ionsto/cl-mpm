(in-package :cl-mpm/dynamic-relaxation)
;; (defmethod %converge-quasi-static ((sim mpm-sim-dr-multigrid)
;;                                    energy-crit
;;                                    oobf-crit
;;                                    live-plot
;;                                    dt-scale
;;                                    substeps
;;                                    conv-steps
;;                                    post-iter-step
;;                                    convergance-criteria
;;                                    kinetic-damping
;;                                    damping-factor)
;;   (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
;;                    (bcs-list cl-mpm::sim-bcs-list)
;;                    (vel-algo cl-mpm::sim-velocity-algorithm))
;;       sim
;;     (let ((refine (cl-mpm::sim-multigrid-refinement sim))
;;           (conv-crit (cl-mpm/dynamic-relaxation::sim-convergence-critera sim))
;;           (mass-filter (/ (cl-mpm::sim-mass-filter sim)
;;                           (expt (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) 2))))
;;       (cl-mpm::iterate-over-mps
;;        (cl-mpm:sim-mps sim)
;;        (lambda (mp)
;;          (setf (fill-pointer (cl-mpm/particle::mp-cached-nodes mp)) 0)))

;;       (format t "Multigrid solve~%")
;;       (let ((damage-state (cl-mpm::sim-enable-damage sim)))
;;         ;; (setf (cl-mpm::sim-enable-damage sim) nil)
;;         (let ((total-iters 0))
;;           (dotimes (i (+ refine 1))
;;             (let ((final-step (= i refine)))
;;               (format t "Mesh step ~D~%" i)
;;               (setf (cl-mpm:sim-mesh sim) (aref mesh-list i)
;;                     (cl-mpm:sim-bcs sim)  (aref bcs-list i))
;;               (setf (cl-mpm::sim-multigrid-current-mesh sim) i)
;;               ;; (when (= (+ 1 i) refine)
;;               ;;   (setf (cl-mpm::sim-enable-damage sim) damage-state))
;;               (pre-step sim)
;;               (format t "Solve final step ~A~%" final-step)
;;               ;; (save-vtks-dr-step sim "./output/" 1000 (* 2 i))
;;               (setf (cl-mpm/dynamic-relaxation::sim-convergence-critera sim)
;;                     (if final-step conv-crit (sqrt conv-crit)))
;;               (format t "Conv crit ~E~%" (cl-mpm/dynamic-relaxation::sim-convergence-critera sim))
;;               (call-next-method
;;                sim
;;                (cl-mpm/dynamic-relaxation::sim-convergence-critera sim)
;;                (cl-mpm/dynamic-relaxation::sim-convergence-critera sim)
;;                live-plot
;;                dt-scale
;;                substeps
;;                conv-steps
;;                (lambda (i o e)
;;                  (funcall post-iter-step total-iters o e)
;;                  (incf total-iters))
;;                convergance-criteria
;;                kinetic-damping
;;                damping-factor)
;;               ;;Remap
;;               (unless final-step
;;                 (setf (cl-mpm:sim-mesh sim) (aref (cl-mpm::sim-mesh-list sim) (+ i 1))
;;                       (cl-mpm:sim-bcs sim)  (aref (cl-mpm::sim-bcs-list sim) (+ i 1)))
;;                 (cl-mpm::iterate-over-mps
;;                  (cl-mpm:sim-mps sim)
;;                  (lambda (mp)
;;                    (setf (fill-pointer (cl-mpm/particle::mp-cached-nodes mp)) 0)))

;;                 (let ((mesh (cl-mpm:sim-mesh sim))
;;                       (mps (cl-mpm:sim-mps sim)))
;;                   (cl-mpm::reset-grid mesh :reset-displacement nil)
;;                   (cl-mpm::p2g mesh mps vel-algo)
;;                   (when (> mass-filter 0d0)
;;                     (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim))))
;;                 (let ((coarse-mesh (aref (cl-mpm::sim-mesh-list sim) i)))
;;                   (cl-mpm::iterate-over-nodes
;;                    (cl-mpm:sim-mesh sim)
;;                    (lambda (n)
;;                      (when (cl-mpm/mesh::node-active n)
;;                        (with-accessors ((disp cl-mpm/mesh::node-displacment))
;;                            n
;;                          (cl-mpm::fast-zero disp)
;;                          (cl-mpm::iterate-over-neighbours-point-linear
;;                           coarse-mesh
;;                           (cl-mpm/mesh::node-position n)
;;                           (lambda (mc coarse-node weight grads)
;;                             (cl-mpm/fastmaths::fast-fmacc
;;                              disp
;;                              (cl-mpm/mesh::node-displacment coarse-node)
;;                              weight))))))))))))))))

;; (defmethod cl-mpm/damage::calculate-damage :around ((sim mpm-sim-dr-multigrid) dt)
;;   (with-accessors ((mesh cl-mpm:sim-mesh)
;;                    (mps cl-mpm:sim-mps)
;;                    (mesh-list cl-mpm::sim-mesh-list))
;;       sim
;;     (let ((m (cl-mpm:sim-mesh sim)))
;;       (setf mesh (aref mesh-list (- (length mesh-list) 1)))
;;       (call-next-method)
;;       (setf mesh m)
;;       )))

;; (defmethod cl-mpm::sim-add-mp :around ((sim mpm-sim-dr-multigrid) mp)
;;   (with-accessors ((mesh cl-mpm::sim-mesh)
;;                    (mesh-list cl-mpm::sim-mesh-list)
;;                    )
;;       sim
;;     (let ((m (cl-mpm:sim-mesh sim)))
;;       (setf mesh (aref mesh-list (- (length mesh-list) 1)))
;;       (call-next-method)
;;       (setf mesh m))))

;; (defmethod pre-step ((sim mpm-sim-dr-multigrid))
;;   (with-accessors ((mesh cl-mpm:sim-mesh)
;;                    (mps cl-mpm:sim-mps)
;;                    (delocal-counter cl-mpm/damage::sim-damage-delocal-counter-max))
;;       sim
;;     (call-next-method)
;;     (setf delocal-counter -1)
;;     (progn
;;       (cl-mpm/damage::update-delocalisation-list (first (last (cl-mpm::sim-mesh-list sim))) mps))))

;; (defmethod pre-step ((sim mpm-sim-dr-multigrid))
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
;;                (ghost-factor cl-mpm::ghost-factor)
;;                (initial-setup initial-setup)
;;                (enable-aggregate cl-mpm/aggregate::enable-aggregate)
;;                (damping cl-mpm::damping-factor)
;;                (vel-algo cl-mpm::velocity-algorithm)
;;                (current-mesh cl-mpm::current-mesh)
;;                (mesh-list cl-mpm::mesh-list)
;;                )
;;       sim
;;     (format t "Current mesh ~D~%" current-mesh)
;;     (if (= current-mesh 0)
;;         (progn
;;           (cl-mpm::reset-grid mesh :reset-displacement t))
;;         (progn
;;           (cl-mpm::reset-grid mesh :reset-displacement nil)))
;;     (setf (cl-mpm/dynamic-relaxation::sim-solve-count sim) 0)
;;     (cl-mpm::p2g mesh mps vel-algo)
;;     (when (> mass-filter 0d0)
;;       (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
;;     (cl-mpm::apply-bcs mesh bcs dt)
;;     (cl-mpm::filter-cells sim)
;;     (cl-mpm::iterate-over-nodes
;;      mesh
;;      (lambda (n)
;;        (setf
;;         (cl-mpm/mesh::node-true-mass n) (cl-mpm/mesh:node-mass n))
;;        (cl-mpm/fastmaths:fast-zero (cl-mpm/mesh::node-true-velocity n))))
;;     (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
;;     (cl-mpm::update-cells sim)
;;     (cl-mpm::apply-bcs mesh bcs dt)
;;     (midpoint-starter sim)
;;     (setf initial-setup t)
;;     (setf (cl-mpm/damage::sim-damage-delocal-counter-max sim) -1)
;;     (when (= current-mesh 0)
;;       (cl-mpm/damage::update-delocalisation-list (aref mesh-list (- (length mesh-list) 1)) mps))))

;; (defmethod cl-mpm/setup::%post-make-simple-sim ((sim cl-mpm::mpm-sim-multigrid) resolution element-count args-list)
;;   (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
;;                    (bcs-list cl-mpm::sim-bcs-list))
;;       sim
;;     (let* ((size (mapcar (lambda (x) (* x resolution)) element-count)))
;;       (let ((mesh-count (+ (cl-mpm::sim-multigrid-refinement sim) 1)))
;;         (setf (cl-mpm::sim-mesh-list sim) (make-array mesh-count)
;;               (cl-mpm::sim-bcs-list sim) (make-array mesh-count))
;;         (dotimes (refine mesh-count)
;;           (let* ((m (cl-mpm::make-mesh size (/ resolution (expt 2 (+ refine 0))) nil))
;;                  (bcs (cl-mpm/bc:make-outside-bc m)))
;;             (setf (aref mesh-list refine) m)
;;             (setf (aref bcs-list refine) bcs)))
;;         (setf (cl-mpm::sim-bcs sim) (aref bcs-list 0))
;;         (set-mesh-default sim)
;;           sim))))

;; (defmethod cl-mpm/setup::%setup-bcs ((sim cl-mpm::mpm-sim-multigrid)
;;                        left
;;                        right
;;                        top
;;                        bottom
;;                        front
;;                        back)
;;   (let ()
;;     (loop for mesh in (cl-mpm::sim-mesh-list sim)
;;           for i from 0
;;           do
;;              (progn
;;                (set-mesh sim i)
;;                (setf (aref (cl-mpm::sim-bcs-list sim) i)
;;                      (cl-mpm/bc::make-outside-bc-varfix
;;                       mesh
;;                       left right top bottom front back))
;;                ;;Avoid triggering the setf bcs active rebuild
;;                (setf (cl-mpm::sim-bcs sim)
;;                      (cl-mpm/bc::make-outside-bc-varfix
;;                       mesh
;;                       left right top bottom front back))
;;                (cl-mpm/setup::resolve-bc-nodes sim mesh (aref (cl-mpm::sim-bcs-list sim) i))
;;                )))
;;   )
(defun iterate-top-down-mesh (sim func)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                    (active-refinement sim-octree-active-refinement)
                    (refinement cl-mpm::sim-multigrid-refinement))
       sim
     (when (> refinement 0)
       (loop for mesh-index from 0 below refinement
             do (let* ((mesh (aref mesh-list mesh-index)))
                  (funcall func mesh mesh-index))))))
(defun iterate-bottom-up-mesh (sim func)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (active-refinement sim-octree-active-refinement)
                   (refinement cl-mpm::sim-multigrid-refinement))
      sim
    (when (> refinement 0)
      (loop for mesh-index from (- refinement 1) downto 0
            do (let* ((mesh (aref mesh-list mesh-index)))
                 (funcall func mesh mesh-index))))))

(defmacro project-vector-down-grid (sim accessor d)
  `(with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                    (active-refinement sim-octree-active-refinement)
                    (refinement cl-mpm::sim-multigrid-refinement))
       ,sim
     (when (> refinement 0)
       (loop for mesh-index from 0 below refinement
             do (let* ((mesh (aref mesh-list mesh-index)))
                  (cl-mpm::iterate-over-cells
                   mesh
                   (lambda (cell)
                     (when (and (= (cell-octree-refine cell) 1))
                       (iterate-over-sub-nodes
                        ,sim
                        mesh-index
                        cell
                        (lambda (n)
                          (setf (cl-mpm/utils:varef (funcall ,accessor n) ,d) 0d0)
                          (cl-mpm::iterate-over-cell-shape-local
                           mesh
                           cell
                           (cl-mpm/mesh::node-position n)
                           (lambda (nw w grads)
                             (declare (ignore grads))
                             (incf (cl-mpm/utils:varef (funcall ,accessor n) ,d)
                                   (*
                                    w
                                    (cl-mpm/utils:varef (funcall ,accessor nw) ,d)))))))))))))))

(defmacro aggregate-vector-up-grid (sim accessor d)
  `(with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                    (active-refinement sim-octree-active-refinement)
                    (refinement cl-mpm::sim-multigrid-refinement))
       ,sim
     (when (> refinement 0)
       (loop for mesh-index from (- refinement 1) downto 0
             do (let* ((mesh (aref mesh-list mesh-index)))
                  (cl-mpm::iterate-over-cells
                   mesh
                   (lambda (cell)
                     (when (and (= (cell-octree-refine cell) 1))
                       (iterate-over-sub-nodes
                        ,sim
                        mesh-index
                        cell
                        (lambda (n)
                          (cl-mpm::iterate-over-cell-shape-local
                           mesh
                           cell
                           (cl-mpm/mesh::node-position n)
                           (lambda (nw w grads)
                             (declare (ignore grads))
                             (incf
                              (cl-mpm/utils:varef (funcall ,accessor nw) ,d)
                              (*
                               w (cl-mpm/utils:varef (funcall ,accessor n) ,d)))))))))))))))

(defmacro project-scalar-down-grid (sim accessor)
  `(with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                    (active-refinement sim-octree-active-refinement)
                    (refinement cl-mpm::sim-multigrid-refinement))
       ,sim
     (when (> refinement 0)
       (loop for mesh-index from 0 below refinement
             do (let* ((mesh (aref mesh-list mesh-index)))
                  (when (< mesh-index refinement)
                    (cl-mpm::iterate-over-cells
                     mesh
                     (lambda (cell)
                       (when (and (cl-mpm/mesh::cell-active cell)
                                  (= (cell-octree-refine cell) 2))
                         (iterate-over-sub-nodes
                          ,sim
                          mesh-index
                          cell
                          (lambda (n)
                            (setf (,accessor n) 0d0)
                            (cl-mpm::iterate-over-cell-shape-local
                             mesh
                             cell
                             (cl-mpm/mesh::node-position n)
                             (lambda (nw w grads)
                               (incf (,accessor n)
                                     (*
                                      w (,accessor nw))))))))))))))))

(defmacro aggregate-scalar-up-grid (sim accessor)
  `(with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                    (active-refinement sim-octree-active-refinement)
                    (refinement cl-mpm::sim-multigrid-refinement))
       ,sim
     (when (> refinement 0)
       (loop for mesh-index from (- refinement 1) downto 0
             do (let* ((mesh (aref mesh-list mesh-index)))
                  (when (< mesh-index refinement)
                    (cl-mpm::iterate-over-cells
                     mesh
                     (lambda (cell)
                       (when (and (cl-mpm/mesh::cell-active cell)
                                  (= (cell-octree-refine cell) 2))
                         (iterate-over-sub-nodes
                          ,sim
                          mesh-index
                          cell
                          (lambda (n)
                            (cl-mpm::iterate-over-cell-shape-local
                             mesh
                             cell
                             (cl-mpm/mesh::node-position n)
                             (lambda (nw w grads)
                               (incf (,accessor nw)
                                     (*
                                      w (,accessor n))))))))))))))))


(defun set-mesh-default (sim)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (mesh cl-mpm::sim-mesh)
                   )
      sim
    (setf (cl-mpm:sim-mesh sim) (aref mesh-list 0))
    ;; (build-active-bcs sim)
    ;; (setf (cl-mpm::sim-bcs sim) (aref (cl-mpm::sim-bcs-list sim) 0))
    ;; (cl-mpm/setup::resolve-bc-nodes sim (cl-mpm::sim-mesh sim) (cl-mpm::sim-bcs sim))
    ;; (setf (cl-mpm::sim-active-bcs sim) (aref (cl-mpm::sim-bcs-list sim) 0))
    ))
(defun set-mesh (sim index)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (mesh cl-mpm::sim-mesh))
      sim
    (setf (cl-mpm:sim-mesh sim) (aref mesh-list index))
    ;; (setf (cl-mpm::sim-active-bcs sim) (aref (cl-mpm::sim-bcs-list sim) index))
    )
  )


(defclass cell-octree (cl-mpm/mesh::cell)
  ((octree-refine
    :initform 0
    :accessor cell-octree-refine
    ;;-1: transition in mesh above
    ;; 0: unrefined
    ;; 1: transition
    ;; 2: refined
    )
   (octree-refine-prev
    :initform 0
    :accessor cell-octree-refine-prev)
   (mesh-index
    :initform 0
    :accessor cell-mesh-index)
   (coupling
    :initform nil
    :accessor cell-coupling)
   (coupled-from-above
    :initform nil
    :accessor cell-coupled-from-above)
   (reaggregation-element
    :initform nil
    :accessor cell-reaggregation-element)))


(defclass node-octree (cl-mpm/mesh::node)
  ((mesh-index
    :initform 0
    :accessor node-mesh-index)
   (transition-zone
    :initform nil
    :accessor node-transition-zone)))

(defclass mesh-octree (cl-mpm/mesh::mesh)
  ())


(defclass mpm-sim-octree (cl-mpm::mpm-sim-multigrid cl-mpm/aggregate::mpm-sim-aggregated)
  ((octree-refinement-criteria
    :initform (lambda (c) nil)
    :accessor sim-octree-refinement-criteria)
   (active-refinement
    :initform 0
    :accessor sim-octree-active-refinement)
   (intra-mesh-aggregation
    :initform t
    :initarg :intra-mesh-agg
    :accessor sim-intra-mesh-aggregation)))


;;abstract class
(defclass mpm-sim-octree-damage (cl-mpm/damage::mpm-sim-damage mpm-sim-octree)
  ())

(defclass mpm-sim-octree-usf (cl-mpm/aggregate::mpm-sim-agg-usf mpm-sim-octree)
  ())

(defclass mpm-sim-octree-damage-usf (cl-mpm/damage::mpm-sim-agg-damage mpm-sim-octree-usf mpm-sim-octree-damage)
  ())

(defclass mpm-sim-octree-quasi-static (cl-mpm/dynamic-relaxation::mpm-sim-quasi-static mpm-sim-octree)
  ())

(defclass mpm-sim-octree-damage-quasi-static (cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul
                                              mpm-sim-octree-quasi-static  mpm-sim-octree-damage)
  ())

(defmethod initialize-instance :after ((sim cl-mpm/dynamic-relaxation::mpm-sim-octree) &key)
  (setf
   (cl-mpm::sim-output-list sim)
   (append
    (cl-mpm::sim-output-list sim)
    (list
     (list :SCALAR "mesh-index" #'cl-mpm/particle::mp-mesh-index)))))

(defun iterate-over-meshes (sim func)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (refinement cl-mpm::sim-multigrid-refinement))
      sim
    (loop for mesh-index from 0 to refinement
          do
             (let ((mesh (aref mesh-list mesh-index)))
               (funcall func mesh mesh-index)))))


(defmethod cl-mpm::apply-essential-bcs ((sim mpm-sim-octree))
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (bcs-list cl-mpm::sim-bcs-list)
                   (refinement cl-mpm::sim-multigrid-refinement)
                   (dt cl-mpm::sim-dt))
      sim
    (loop for mesh-index from 0 to refinement
          do
             (let ((mesh (aref mesh-list mesh-index))
                   (bcs (aref bcs-list mesh-index)))
               (with-accessors ((nodes  mesh-nodes)
                                (nD     mesh-nD))
                   mesh
                 (cl-mpm/utils::bpdotimes
                  (i (length bcs))
                  (let ((bc (aref bcs i)))
                    (when bc
                      (with-accessors ((node cl-mpm/bc::bc-node)
                                       (index cl-mpm/bc::bc-index))
                          bc
                        (if (or node
                                (not index))
                            (cl-mpm/bc:apply-bc bc node mesh dt)
                            (progn
                              (setf node (cl-mpm/mesh:get-node mesh index))
                              (if node
                                  (cl-mpm/bc:apply-bc bc node mesh dt)
                                  (error "BC attempted to get a nil node ~A ~A" bc index)))))))))))))

(defun build-active-bcs (sim)
  (with-accessors ((bcs-list cl-mpm::sim-bcs-list)
                   (refinement cl-mpm::sim-multigrid-refinement))
      sim
    (let ((bc-count 0))
      (loop for bcs across bcs-list
            do (incf bc-count (length bcs)))
      (let ((new-bcs-array (make-array bc-count)))
        (let ((i 0))
          (loop for bcs across bcs-list
                do (loop for n across bcs
                         do (progn
                              (setf (aref new-bcs-array i) n)
                              (incf i)))))
        (setf (cl-mpm::sim-bcs sim)
              new-bcs-array)
        (setf (cl-mpm::sim-active-bcs sim)
              new-bcs-array)))))


(defun resolve-all-bcs (sim)
  (loop for mesh across (cl-mpm::sim-mesh-list sim)
        for bcs across (cl-mpm::sim-bcs-list sim)
        do
           (progn
             (cl-mpm/setup::resolve-bc-nodes sim
                                             mesh
                                             bcs))))


(defmethod cl-mpm/setup::%setup-bcs ((sim cl-mpm/dynamic-relaxation::mpm-sim-octree)
                                     left
                                     right
                                     top
                                     bottom
                                     front
                                     back)
  (set-mesh-default sim)
  (with-accessors ((bcs-list cl-mpm::sim-bcs-list)
                   (refinement cl-mpm::sim-multigrid-refinement)
                   (mesh-list cl-mpm::sim-mesh-list))
      sim

    (let ((bc-count 0))
      (loop for mesh across (cl-mpm::sim-mesh-list sim)
            for i from 0
            do
               (progn
                 (set-mesh sim i)
                 (setf (aref bcs-list i)
                       (cl-mpm/bc::make-outside-bc-varfix
                        (cl-mpm::sim-mesh sim)
                        left right top bottom front back))
                 ;;Avoid triggering the setf bcs active rebuild
                 (cl-mpm/setup::resolve-bc-nodes sim
                                                 (cl-mpm::sim-mesh sim)
                                                 (aref bcs-list i))
                 (incf bc-count (length (aref bcs-list i)))))
      (let ((new-bcs-array (make-array bc-count)))
        (let ((i 0))
          (loop for bcs across bcs-list
                do (loop for n across bcs
                         do (progn
                              (setf (aref new-bcs-array i) n)
                              (incf i)))))
        (setf (cl-mpm::sim-bcs sim)
              new-bcs-array)
        (setf (cl-mpm::sim-active-bcs sim)
              new-bcs-array)))

    (set-mesh-default sim)
    ;; (setf
    ;;  (cl-mpm::sim-bcs sim)
    ;;  (aref bcs-list 0))

    ;; (setf
    ;;  (cl-mpm::sim-bcs sim)
    ;;  (aref bcs-list 0))
    ))

(defmethod update-instance-for-different-class ((prev cl-mpm::mpm-sim)
                                                (current mpm-sim-octree)
                                                &rest initargs
                                                  ;; &key (refinement 0)
                                                  )
  ;; (setf (cl-mpm::sim-multigrid-refinement current) refinement)
  (call-next-method)
  (pprint "Rebuild meshes")
  (cl-mpm:iterate-over-mps
   (cl-mpm:sim-mps current)
   (lambda (mp)
     (setf (cl-mpm/particle::mp-damage-position mp) nil)))
  (cl-mpm/setup::%post-make-simple-sim
   current
   (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh prev))
   (subseq
    (cl-mpm/mesh::mesh-count (cl-mpm:sim-mesh prev))
    0
    (cl-mpm/mesh::mesh-nd (cl-mpm:sim-mesh prev)))
   initargs)
  (apply #'cl-mpm/setup:setup-bcs
         current
         (cl-mpm/mesh::mesh-boundary-bcs (cl-mpm:sim-mesh prev))))

(defmethod update-instance-for-different-class ((prev mpm-sim-octree)
                                                (current mpm-sim-octree)
                                                &rest initargs)
  (pprint "We don't have do do anything")
  (call-next-method)
  (cl-mpm:iterate-over-mps
   (cl-mpm:sim-mps current)
   (lambda (mp)
     (setf (cl-mpm/particle::mp-damage-position mp) nil)))
  ;; (setf (cl-mpm::sim-multigrid-refinement current) refinement)
  ;; (cl-mpm/setup::%post-make-simple-sim
  ;;  current
  ;;  (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh prev))
  ;;  (subseq
  ;;   (cl-mpm/mesh::mesh-count (cl-mpm:sim-mesh prev))
  ;;   0
  ;;   (cl-mpm/mesh::mesh-nd (cl-mpm:sim-mesh prev))
  ;;   )
  ;;  initargs)
  ;; (apply #'cl-mpm/setup:setup-bcs
  ;;        current
  ;;        (cl-mpm/mesh::mesh-boundary-bcs (cl-mpm:sim-mesh prev)))
  )


(defmethod cl-mpm/setup::%post-make-simple-sim ((sim mpm-sim-octree) resolution element-count args-list)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (refinement cl-mpm::sim-multigrid-refinement)
                   (bcs-list cl-mpm::sim-bcs-list))
      sim
    (let* ((size (mapcar (lambda (x) (* x resolution)) element-count)))
      (let ((mesh-count (+ refinement 1)))
        (setf (cl-mpm::sim-mesh-list sim) (make-array mesh-count)
              (cl-mpm::sim-bcs-list sim) (make-array mesh-count))
        (dotimes (refine (+ (cl-mpm::sim-multigrid-refinement sim) 1))
          (let* ((m (cl-mpm::make-mesh size (/ resolution (expt 2 (+ refine 0))) nil)))
            (setf (cl-mpm::sim-mesh sim) m)
            (setf (aref mesh-list refine) m)
            (setf (aref bcs-list refine) (cl-mpm/bc:make-outside-bc m))
            ;; (setf (cl-mpm:sim-bcs sim) (aref bcs-list refine))
            (cl-mpm/setup::resolve-bc-nodes
             sim
             (cl-mpm::sim-mesh sim)
             (aref bcs-list refine))
            ))
        ;; (cl-mpm::apply-essential-bcs sim)

        (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                         (bcs-list cl-mpm::sim-bcs-list))
            sim
          (let ((nc 0))
            (dotimes (i (1+ refinement))
              (let ((mesh (aref mesh-list i)))
                (incf nc (length (cl-mpm/mesh::mesh-active-nodes mesh)))
                (cl-mpm::iterate-over-cells
                 mesh
                 (lambda (c)
                   (change-class c 'cell-octree)
                   (setf (cell-mesh-index c) i)))
                (cl-mpm::iterate-over-nodes
                 mesh
                 (lambda (n)
                   (change-class n 'node-octree)
                   (setf (node-mesh-index n) i)))))

            (set-mesh-default sim)

            (let ((new-node-array (make-array nc)))
              (let ((i 0))
                (loop for mesh across mesh-list
                      do (loop for n across (cl-mpm/mesh::mesh-active-nodes mesh)
                               do (progn
                                    (setf (aref new-node-array i) n)
                                    (incf i)))))
              (set-mesh-default sim)
              (setf (cl-mpm/mesh::mesh-active-nodes (cl-mpm:sim-mesh sim))
                    new-node-array))))

        ;; (set-mesh-default sim)
        sim)))
  ;; (call-next-method)
  )


(defmethod cl-mpm::iterate-over-neighbours-point ((sim cl-mpm/dynamic-relaxation::mpm-sim-octree) pos func)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (refinement cl-mpm::sim-multigrid-refinement))
      sim
    (let* ((resolved-mesh 0)
           (mesh (aref mesh-list resolved-mesh)))
      (loop for i from 0 to refinement
            do
               (let ()
                 (setf mesh (aref mesh-list resolved-mesh))
                 (labels ((get-cell-at-pos (pos)
                            (let ((p pos))
                              (cl-mpm/mesh::get-cell
                               mesh
                               (cl-mpm/mesh::position-to-index-cell
                                mesh
                                p
                                )))))
                   (let ((max-index (cell-octree-refine (get-cell-at-pos pos))))
                     (if (or
                          (> max-index 1))
                         (progn
                           (incf resolved-mesh))
                         (loop-finish))))))
      (cl-mpm::iterate-over-neighbours-point-linear mesh pos func))))


(defun setup-mp-iteration-cache (sim)
  (set-mesh-default sim)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (refinement cl-mpm::sim-multigrid-refinement))
      sim
      (cl-mpm::iterate-over-mps
       (cl-mpm:sim-mps sim)
       (lambda (mp)
         (if (> refinement 0)
           (let ((resolved-mesh 0)
                 (running t))
             (let ((mesh (aref mesh-list resolved-mesh)))
               (loop for i from 0 to refinement
                     do
                        (let ()
                          (setf mesh (aref mesh-list resolved-mesh))
                          (labels ((get-cell-at-pos (pos)
                                     (let ((p (cl-mpm/utils::vector-copy pos)))
                                       ;; (cl-mpm/mesh::clamp-point-to-bounds mesh p)
                                       (cl-mpm/mesh::get-cell
                                        mesh
                                        (mapcar (lambda (x) (min (max 0 x)))
                                                (cl-mpm/mesh::position-to-index-cell
                                                 mesh
                                                 p
                                                 ))))))

                            ;; (if (or
                            ;;      (= 1 (cell-octree-refine (get-cell-at-pos (cl-mpm/particle::mp-position mp))))
                            ;;      (= 2 (cell-octree-refine (get-cell-at-pos (cl-mpm/particle::mp-position mp)))))
                            ;;     (progn
                            ;;       (incf resolved-mesh)))
                            (let ((max-index -1)
                                  (min-index 3))
                              (cl-mpm::iterate-over-corners
                               mesh
                               mp
                               (lambda (corner)
                                 (let ((cell (get-cell-at-pos corner)))
                                   (setf max-index (max max-index (cell-octree-refine cell)))
                                   (setf min-index (min min-index (cell-octree-refine cell))))))
                              (if (or
                                   (> max-index 1)
                                   (and
                                    (= max-index 1)
                                    (= min-index 1))
                                   ;; (> min-index 0)
                                   )
                                  (progn
                                    (incf resolved-mesh))
                                  (loop-finish)))
                            )))

               ;; (setf resolved-mesh 0)
               ;; (pprint resolved-mesh)

               (cl-mpm/particle::reset-mp-node-cache mp)
               (setf (cl-mpm/particle::mp-index mp) resolved-mesh)
               (setf (cl-mpm/particle::mp-mesh-index mp) resolved-mesh)
               (cl-mpm::iterate-over-neighbours
                mesh
                mp
                (lambda (&rest args)))))
           (progn
             (setf (cl-mpm/particle::mp-index mp) 0)
             (setf (cl-mpm/particle::mp-mesh-index mp) 0)))))
    ;; (iterate-over-meshes
    ;;  sim
    ;;  (lambda (mesh mesh-index)
    ;;    (cl-mpm::iterate-over-cells
    ;;     mesh
    ;;     (lambda (cell)
    ;;       (when (= (cell-octree-refine cell) 1)
    ;;         (let ((any-active-cells nil))
    ;;           (with-accessors ((partial cl-mpm/mesh::cell-partial)
    ;;                            (nodes cl-mpm/mesh::cell-nodes)
    ;;                            (agg cl-mpm/mesh::cell-agg)
    ;;                            (neighbours cl-mpm/mesh::cell-cartesian-neighbours)
    ;;                            (active cl-mpm/mesh::cell-active))
    ;;               cell
    ;;             (loop for n in neighbours
    ;;                   do
    ;;                      (when (and (cl-mpm/mesh::cell-active n)
    ;;                                 (= (cell-octree-refine n) 0)
    ;;                                 (= (cell-octree-refine n) -1))
    ;;                        (setf any-active-cells t))))
    ;;           (unless any-active-cells
    ;;             (setf (cell-octree-refine cell) 2))))))))
    (set-mesh-default sim)))

(defun intra-mesh-agg (sim mesh)
  (setf (cl-mpm::sim-mesh sim) mesh)
  (let ((mesh-nodes
          (make-array (array-total-size (cl-mpm/mesh::mesh-nodes mesh))
                      :displaced-to
                      (cl-mpm/mesh::mesh-nodes mesh))))

    (cl-mpm::iterate-over-nodes-array
     mesh-nodes
     (lambda (n)
       (setf (cl-mpm/mesh::node-agg-building-flag n) nil
             (node-transition-zone n) nil)
       ;; (when  (and (cl-mpm/mesh::node-active n)
       ;;             (< (/ (cl-mpm/mesh::node-volume n)
       ;;                   (cl-mpm/mesh::node-volume-true n))
       ;;                0.1d0))
       ;;   (setf (cl-mpm/mesh::node-agg n) t))
       ))
    (cl-mpm::iterate-over-cells
     mesh
     (lambda (c)
       (setf (cl-mpm/mesh::cell-agg c) nil
             (cl-mpm/mesh::cell-interior c) nil)
       ;; (loop for n across (cl-mpm/mesh::cell-nodes c)
       ;;       do (when (and (cl-mpm/mesh::node-active n)
       ;;                     (< (/ (cl-mpm/mesh::node-volume n)
       ;;                           (cl-mpm/mesh::node-volume-true n))
       ;;                        0.1d0))
       ;;            (setf (cl-mpm/mesh::cell-partial c) t)))
       ))

    ;;First set all outside cells as aggregate
    (cl-mpm::iterate-over-cells
     mesh
     (lambda (cell)
       (with-accessors ((partial cl-mpm/mesh::cell-partial)
                        (nodes cl-mpm/mesh::cell-nodes)
                        (agg cl-mpm/mesh::cell-agg)
                        (neighbours cl-mpm/mesh::cell-cartesian-neighbours)
                        (active cl-mpm/mesh::cell-active))
           cell
         (when (and active partial)
           (setf agg t)
           (loop for n in neighbours
                 do
                    (when (cl-mpm/mesh::cell-active n)
                      (setf (cl-mpm/mesh::cell-agg n) t)))))))

    ;;Next set any nodes on agg elements to be agg
    (cl-mpm::iterate-over-cells
     mesh
     (lambda (cell)
       (with-accessors ((nodes cl-mpm/mesh::cell-nodes)
                        (agg cl-mpm/mesh::cell-agg)
                        (coupling cell-coupling)
                        (active cl-mpm/mesh::cell-active))
           cell
         (when (and active agg)
           (unless coupling
             (loop for n across nodes
                   do (when (and (cl-mpm/mesh::node-active n))
                        (setf (cl-mpm/mesh::node-agg n) t))))
           (when coupling
             (loop for n across nodes
                   do (when (cl-mpm/mesh::node-active n)
                        (setf (cl-mpm/mesh::node-interior n) nil
                              (cl-mpm/mesh::node-agg n) t
                              (node-transition-zone n) t)))))

         (when (and active coupling)
           (setf (cl-mpm/mesh::cell-interior cell) t)
           (loop for n across nodes
                 do (when (cl-mpm/mesh::node-active n)
                      (setf
                       (cl-mpm/mesh::node-agg n) t
                       (node-transition-zone n) t)))))))

    (cl-mpm::iterate-over-cells
     mesh
     (lambda (c)
       (loop for n across (cl-mpm/mesh::cell-nodes c)
             do (when (and (cl-mpm/mesh::node-active n)
                           (< (/ (cl-mpm/mesh::node-volume n)
                                 (cl-mpm/mesh::node-volume-true n))
                              0.1d0))
                  (setf (cl-mpm/mesh::cell-partial c) t)))))

    ;;For each aggregate node, locate the closest support cell
    (cl-mpm::iterate-over-nodes-array
     mesh-nodes
     (lambda (node)
       (with-accessors ((active cl-mpm/mesh::node-active)
                        (agg cl-mpm/mesh::node-agg)
                        (int cl-mpm/mesh::node-interior)
                        (building cl-mpm/mesh::node-agg-building-flag)
                        (transition node-transition-zone)
                        (int-cell cl-mpm/mesh::node-agg-interior-cell))
           node
         (when (and active agg (not transition))
           (let ((closest-elem (cl-mpm/aggregate::get-closest-cell
                                mesh
                                (cl-mpm/mesh::node-position node)
                                ;; :exclude int-cell
                                :filter (lambda (c)
                                          (and
                                           ;; (not (cell-coupling c))
                                           (or
                                            (= (cell-octree-refine c) 0)
                                            (= (cell-octree-refine c) 1))
                                           ;; (not (= (cell-octree-refine c) 2))
                                           ))
                                )))
             (if closest-elem
                 (progn
                   (setf int-cell closest-elem)
                   (setf (cl-mpm/mesh::cell-interior closest-elem) t))
                 (progn
                   (setf (cl-mpm/mesh::node-agg node) nil)
                   (format t "No closest elem? ~D - ~A~%" (node-mesh-index node) (cl-mpm/mesh::node-index node))
                   ;; (error "No closest elem? ~A~%")
                   ))))))))
  (cl-mpm::iterate-over-cells
   mesh
   (lambda (cell)
     (with-accessors ((nodes cl-mpm/mesh::cell-nodes)
                      (agg cl-mpm/mesh::cell-agg)
                      (int cl-mpm/mesh::cell-interior)
                      (active cl-mpm/mesh::cell-active)
                      (reagg cell-reaggregation-element)
                      (coupling cell-coupling)
                      )
         cell
       (when (and active agg coupling)
         (setf int nil)
         (let ((closest-elem (cl-mpm/aggregate::get-closest-cell
                              mesh
                              (cl-mpm/mesh::cell-centroid cell)
                              :exclude cell
                              :filter (lambda (c)
                                        (and
                                         ;; (not (cell-coupling c))
                                         ;; (not (= (cell-octree-refine c) 2))
                                         (not (= (cell-octree-refine c) 2))
                                         ;; (or
                                         ;;  (= (cell-octree-refine c) 0)
                                         ;;  (= (cell-octree-refine c) 1)
                                         ;;  )
                                         )))))
           (if closest-elem
               (progn
                 (reaggregate-cell cell closest-elem)
                 ;; (setf
                 ;;  (cell-reaggregation-element cell)
                 ;;  closest-elem)
                 ;; (setf (cl-mpm/mesh::cell-interior closest-elem) t)
                 )
               (progn
                 (setf (cl-mpm/mesh::cell-agg cell) nil)
                 (setf (cl-mpm/mesh::cell-interior cell) nil)
                 (format t "No closest reaggregation cell? ~D - ~A~%" (cell-mesh-index cell) (cl-mpm/mesh::cell-index cell)))))))))
  (cl-mpm::iterate-over-cells
   mesh
   (lambda (cell)
     (with-accessors ((nodes cl-mpm/mesh::cell-nodes)
                      (agg cl-mpm/mesh::cell-agg)
                      (int cl-mpm/mesh::cell-interior)
                      (active cl-mpm/mesh::cell-active)
                      (coupling cell-coupling)
                      )
         cell
       (when (and active int ;; (not coupling)
                  )
         (loop for n across nodes
               do (when (and (cl-mpm/mesh::node-active n)
                             ;; (not (node-transition-zone n))
                             ;; (cl-mpm/mesh::node-agg-building-flag n)
                             ;; (not (cl-mpm/mesh::node-agg-interior-cell n))
                             )
                    (setf (cl-mpm/mesh::node-agg n) t
                          (cl-mpm/mesh::node-agg-interior-cell n) cell
                          (cl-mpm/mesh::node-interior n) t))))))))

(defun reaggregate-cell (cell new-cell)
  (with-accessors ((nodes cl-mpm/mesh::cell-nodes)
                   (agg cl-mpm/mesh::cell-agg)
                   (int cl-mpm/mesh::cell-interior)
                   (active cl-mpm/mesh::cell-active)
                   (reagg cell-reaggregation-element)
                   (coupling cell-coupling)
                   )
      cell
    (setf int nil)
    (progn
      ;; (setf
      ;;  (cell-reaggregation-element cell)
      ;;  new-cell)
      (loop for n across nodes
            do (when (cl-mpm/mesh::node-active n)
                 (setf (cl-mpm/mesh::node-interior n) nil
                       (cl-mpm/mesh::node-agg n) t
                       (cl-mpm/mesh::node-agg-interior-cell n) new-cell
                       (node-transition-zone n) t)))
      (setf (cl-mpm/mesh::cell-interior new-cell) t)
      (setf (cl-mpm/mesh::cell-agg new-cell) t)
      )))

(defun inter-mesh-agg (sim mesh-index)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list))
      sim
    (let ((mesh (aref mesh-list mesh-index))
          (mesh-below (aref mesh-list (1+ mesh-index))))
        (cl-mpm::iterate-over-cells
         mesh
         (lambda (cell)
           (with-accessors ((nodes cl-mpm/mesh::cell-nodes)
                            (active cl-mpm/mesh::cell-active)
                            (index cl-mpm/mesh::cell-index))
               cell
             (when (= (cell-octree-refine cell) 1)
               (when active
                 (setf (cl-mpm/mesh::cell-interior cell) t)
                 (setf (cl-mpm/mesh::cell-agg cell) t)
                 (setf (cell-coupling cell) t)

                 (iterate-over-sub-cells
                  sim
                  mesh-index
                  cell
                  (lambda (c)
                    (reaggregate-cell c cell)
                    ;; (setf (cl-mpm/mesh::cell-interior c) nil)
                    ;; (when (cl-mpm/mesh::cell-agg-int c)
                    ;;   (setf (cell-coupling c) t)
                    ;;   (setf (cell-reaggregation-element c) cell))
                    ))

                 (let ((sub-node-active nil))
                   ;; (setf sub-node-active t)
                   (iterate-over-sub-nodes
                    sim
                    mesh-index
                    cell
                    (lambda (n)
                      (when (cl-mpm/mesh::node-active n)
                        (setf sub-node-active t)
                        ;;When we are "free" we can just grab them by pointing
                        (when t
                          (setf
                           (cl-mpm/mesh::node-agg n) t
                           (cl-mpm/mesh::node-agg-interior-cell n) cell
                           (cl-mpm/mesh::node-interior n) nil))
                        ;; Interior
                        (cl-mpm::iterate-over-cell-shape-local
                         mesh
                         cell
                         (cl-mpm/mesh::node-position n)
                         (lambda (n w grads)
                           (when (> (abs w) 1d-9)
                             (setf (cl-mpm/mesh::node-active n) t)))))))

                   (when t;sub-node-active
                     (loop for n across nodes
                           do (progn
                                (when (cl-mpm/mesh::node-active n)
                                  (setf (cl-mpm/mesh::node-agg n) t
                                        (node-transition-zone n) t
                                        (cl-mpm/mesh::node-interior n) t)
                                  ;; (unless (cl-mpm/mesh::node-agg-interior-cell n)
                                  ;;   (setf (cl-mpm/mesh::node-agg-interior-cell n) cell))
                                  (setf (cl-mpm/mesh::node-agg-interior-cell n) cell)
                                  )
                                ;; (setf (cl-mpm/mesh::node-active n) t)
                                ))))
                 )))))

      ;; (cl-mpm::iterate-over-cells
      ;;  mesh-below
      ;;  (lambda (cell)
      ;;    (when (cl-mpm/mesh::cell-interior cell)
      ;;      (loop for n across (cl-mpm/mesh::cell-nodes cell)
      ;;            do (progn
      ;;                 (setf (cl-mpm/mesh::node-agg n) t
      ;;                       (cl-mpm/mesh::node-interior n) nil))))))
      )))


(defmethod cl-mpm::update-cells ((sim mpm-sim-octree))
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list))
      sim
    (loop for mesh across mesh-list
          do
             (with-accessors ((dt cl-mpm::sim-dt))
                 sim
               (cl-mpm::iterate-over-cells
                mesh
                (lambda (cell)
                  (cl-mpm::filter-cell mesh cell dt)
                  (cl-mpm::update-cell mesh cell dt)))))

    (propogate-filter-cells-octree sim)))

(defun iterate-over-sub-cells-1d (sim mesh-index cell func)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (dt cl-mpm::sim-dt)
                   (agg cl-mpm/aggregate::sim-enable-aggregate))
      sim
    (let ((mesh (aref mesh-list (1+ mesh-index))))
      (destructuring-bind (x y z) (cl-mpm/mesh::cell-index cell)
        (dotimes (dx 2)
            (let ((nb (cl-mpm/mesh::get-cell
                       mesh
                       (list (+ dx (* 2 x))
                             0
                             0))))
              (funcall func nb)))))))

(defun iterate-over-sub-cells-2d (sim mesh-index cell func)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list))
      sim
    (let ((mesh (aref mesh-list (1+ mesh-index))))
      (destructuring-bind (x y z) (cl-mpm/mesh::cell-index cell)
        (dotimes (dx 2)
          (dotimes (dy 2)
            (let ((nb (cl-mpm/mesh::get-cell
                       mesh
                       (list (+ dx (* 2 x))
                             (+ dy (* 2 y))
                             0))))
              (funcall func nb))))))))
(defun iterate-over-sub-cells-3d (sim mesh-index cell func)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list))
      sim
    (let ((mesh (aref mesh-list (1+ mesh-index))))
      (destructuring-bind (x y z) (cl-mpm/mesh::cell-index cell)
        (dotimes (dx 2)
          (dotimes (dy 2)
            (dotimes (dz 2)
              (let ((nb (cl-mpm/mesh::get-cell
                         mesh
                         (list (+ dx (* 2 x))
                               (+ dy (* 2 y))
                               (+ dz (* 2 z))))))
                (funcall func nb)))))))))

(defun iterate-over-sub-cells (sim mesh-index cell func)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (dt cl-mpm::sim-dt)
                   (agg cl-mpm/aggregate::sim-enable-aggregate))
      sim
    (let ((mesh (aref mesh-list mesh-index)))
      (ecase (cl-mpm/mesh::mesh-nd mesh)
        (1 (iterate-over-sub-cells-1d sim mesh-index cell func))
        (2 (iterate-over-sub-cells-2d sim mesh-index cell func))
        (3 (iterate-over-sub-cells-3d sim mesh-index cell func))))))

(defun iterate-over-sub-nodes-1d (sim mesh-index cell func)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (dt cl-mpm::sim-dt)
                   (agg cl-mpm/aggregate::sim-enable-aggregate))
      sim
    (let ((mesh (aref mesh-list (1+ mesh-index))))
      (destructuring-bind (x y z) (cl-mpm/mesh::cell-index cell)
        (dotimes (dx 3)
          (let ((n (cl-mpm/mesh::get-node
                    mesh
                    (list (+ (* 2 x) dx)
                          0 0))))
            (funcall func n)))))))

(defun iterate-over-sub-nodes-2d (sim mesh-index cell func)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (dt cl-mpm::sim-dt)
                   (agg cl-mpm/aggregate::sim-enable-aggregate))
      sim
    (let ((mesh (aref mesh-list (1+ mesh-index))))
      (destructuring-bind (x y z) (cl-mpm/mesh::cell-index cell)
        (dotimes (dx 3)
          (dotimes (dy 3)
            (let ((n (cl-mpm/mesh::get-node
                      mesh
                      (list (+ (* 2 x) dx)
                            (+ (* 2 y) dy) 0))))
              (funcall func n))))))))
(defun iterate-over-sub-nodes-3d (sim mesh-index cell func)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (dt cl-mpm::sim-dt)
                   (agg cl-mpm/aggregate::sim-enable-aggregate))
      sim
    (let ((mesh (aref mesh-list (1+ mesh-index))))
      (destructuring-bind (x y z) (cl-mpm/mesh::cell-index cell)
        (dotimes (dx 3)
          (dotimes (dy 3)
            (dotimes (dz 3)
              (let ((n (cl-mpm/mesh::get-node
                        mesh
                        (list (+ (* 2 x) dx)
                              (+ (* 2 y) dy)
                              (+ (* 2 z) dz)))))
                (funcall func n)))))))))

(defun iterate-over-sub-nodes (sim mesh-index cell func)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (dt cl-mpm::sim-dt)
                   (agg cl-mpm/aggregate::sim-enable-aggregate))
      sim
    (let ((mesh (aref mesh-list mesh-index)))
      (ecase (cl-mpm/mesh::mesh-nd mesh)
        (1 (iterate-over-sub-nodes-1d sim mesh-index cell func))
        (2 (iterate-over-sub-nodes-2d sim mesh-index cell func))
        (3 (iterate-over-sub-nodes-3d sim mesh-index cell func))))))

(defun propogate-filter-cells-octree (sim)
  (with-accessors (;; (mesh cl-mpm::sim-mesh)
                   (mesh-list cl-mpm::sim-mesh-list)
                   (refinement cl-mpm::sim-multigrid-refinement)
                   (agg cl-mpm/aggregate::sim-enable-aggregate))
      sim
    ;;First when we propogate from the top down whether cells know they are fully populated
    (loop for i from 0 below refinement
          do
             (let ((mesh (aref mesh-list i)))
               (cl-mpm::iterate-over-cells
                mesh
                (lambda (cell)
                  (with-accessors ((active cl-mpm/mesh::cell-active)
                                   (coupling cell-coupling)
                                   (partial cl-mpm/mesh::cell-partial))
                      cell
                    (when (= (cell-octree-refine cell) 2)
                      (setf (cl-mpm/mesh::cell-active cell) nil))
                    (when (and active
                               ;; (not coupling)
                               )
                      (iterate-over-sub-cells
                       sim
                       i
                       cell
                       (lambda (c)
                         (setf (cl-mpm/mesh::cell-active c)
                               (or (cl-mpm/mesh::cell-active c)
                                   active))
                         (unless partial
                           (setf (cl-mpm/mesh::cell-partial c) nil))))))))))
    ;;Next we check cells that do not know if they are full by gathering that infomation recursivly
    (when (> refinement 0)
      (loop for i from (- refinement 1) downto 0
            do
               (let ((mesh (aref mesh-list i)))
                 (cl-mpm::iterate-over-cells
                  mesh
                  (lambda (cell)
                    (with-accessors ((active cl-mpm/mesh::cell-active)
                                     (partial cl-mpm/mesh::cell-partial))
                        cell
                      (when active
                        (let ((trial-partial nil))
                          (iterate-over-sub-cells
                           sim
                           i
                           cell
                           (lambda (c)
                             (setf trial-partial (or trial-partial (cl-mpm/mesh::cell-partial c)))))
                          (setf partial trial-partial)))))))))))

(defmethod cl-mpm::filter-cells ((sim mpm-sim-octree))
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (mesh-list cl-mpm::sim-mesh-list)
                   (dt cl-mpm::sim-dt)
                   (agg cl-mpm/aggregate::sim-enable-aggregate))
      sim
    (loop for mesh across mesh-list
          do (cl-mpm::iterate-over-cells
              mesh
              (lambda (cell)
                (cl-mpm::filter-cell mesh cell dt))))
    (propogate-filter-cells-octree sim)
    (when agg
      (cl-mpm/aggregate::update-aggregate-elements sim))))


(defun build-aggregate-bcs (sim)
  (when (cl-mpm/aggregate::sim-enable-aggregate sim)
    (let ((nd (cl-mpm/mesh::mesh-nd (cl-mpm:sim-mesh sim))))
      (setf (cl-mpm/aggregate::sim-global-bcs sim) (make-array nd :element-type t))
      (setf (cl-mpm/aggregate::sim-global-bcs-int sim) (make-array nd :element-type t))
      (cl-mpm/aggregate::iterate-over-dimensions
       nd
       (lambda (d)
         (setf (aref (cl-mpm/aggregate::sim-global-bcs sim) d) (cl-mpm/aggregate::assemble-global-bcs sim d))
         (setf (aref (cl-mpm/aggregate::sim-global-bcs-int sim) d) (cl-mpm/aggregate::assemble-internal-bcs sim d)))))))

(defmethod cl-mpm/aggregate::locate-aggregate-nodes ((sim mpm-sim-octree))
  (set-mesh-default sim)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)
                   (mesh-list cl-mpm::sim-mesh-list)
                   (intra-mesh-agg sim-intra-mesh-aggregation)
                   (refinement cl-mpm::sim-multigrid-refinement))
      sim
    ;;Reset aggregation structures
    (cl-mpm:iterate-over-nodes
     mesh
     (lambda (node)
       (setf (cl-mpm/mesh::node-agg-interior-cell node) nil
             (cl-mpm/mesh::node-interior node) nil
             (cl-mpm/mesh::node-agg node) nil
             (node-transition-zone node) nil)))
    (loop for mesh across mesh-list
          do (cl-mpm::iterate-over-cells
              mesh
              (lambda (node)
                (setf (cl-mpm/mesh::cell-agg node) nil
                      (cl-mpm/mesh::cell-aggregate-element node) nil
                      (cell-coupling node) nil
                      ;; (cell-coupled-from-above node) nil
                      (cell-reaggregation-element node) nil
                      (cl-mpm/mesh::cell-interior node) nil))))
    ;;First set all outside nodes as aggregate
    (when (> refinement 0)
      ;;Build the inter-mesh aggregation
      (loop for mesh-index from (- refinement 1) downto 0
            do (progn
                 (when intra-mesh-agg
                   (set-mesh sim (1+ mesh-index))
                   (intra-mesh-agg sim (aref mesh-list (1+ mesh-index))))
                 (set-mesh sim mesh-index)
                 (inter-mesh-agg sim mesh-index))))
    (when intra-mesh-agg
      (set-mesh sim 0)
      (intra-mesh-agg sim (aref mesh-list 0)))
    (set-mesh-default sim)

    (when (> refinement 0)
      (labels ((resolve-reaggregation (cell)
                 (let ((current-cell cell))
                   (loop while (cell-reaggregation-element current-cell)
                         do
                            (progn
                              (when (eq (cell-reaggregation-element current-cell) cell) nil
                                    (error "Recursive reaggregation setup? ~A at ~A" cell (cl-mpm/mesh::cell-index cell)))
                              (setf current-cell (cell-reaggregation-element current-cell))))
                   current-cell)))
        (set-mesh-default sim)
        (cl-mpm::iterate-over-nodes
         mesh
         (lambda (n)
           (when (cl-mpm/mesh::node-active n)
             (when (cl-mpm/mesh::node-agg n)
               (with-accessors ((int-cell cl-mpm/mesh::node-agg-interior-cell))
                   n
                 (setf int-cell (resolve-reaggregation int-cell)))))))
        ;; (loop for mesh across mesh-list
        ;;       do (progn
        ;;            (cl-mpm::iterate-over-cells
        ;;             mesh
        ;;             (lambda (c)
        ;;               (with-accessors ((reagg cell-reaggregation-element))
        ;;                   c
        ;;                 (when reagg
        ;;                   ;;We need to kick off the reaggregation process
        ;;                   (setf (cl-mpm/mesh::cell-interior c) nil)
        ;;                   (let ((resolved (resolve-reaggregation c)))
        ;;                     (loop for n across (cl-mpm/mesh::cell-nodes c)
        ;;                           do (progn
        ;;                                (setf (cl-mpm/mesh::node-agg-interior-cell n) resolved
        ;;                                      (cl-mpm/mesh::node-interior n) nil))))))))))
        ))
    (set-mesh-default sim)
    ;;For each aggregate node, locate the closest support cell
    (setf
     (cl-mpm/aggregate::sim-agg-nodes-fdc sim)
     (lparallel:premove-if-not
      (lambda (n) (and (cl-mpm/mesh::node-active n)
                       (cl-mpm/mesh::node-interior n)))
      (cl-mpm/mesh::mesh-active-nodes mesh))

     (cl-mpm/aggregate::sim-agg-nodes-fd sim)
     (lparallel:premove-if-not
      (lambda (n) (and (cl-mpm/mesh::node-active n)
                       (cl-mpm/mesh::node-agg n)))
      (cl-mpm/mesh::mesh-active-nodes mesh)))

    (let ((nodes (cl-mpm/aggregate::sim-agg-nodes-fdc sim)))
      (cl-mpm/utils::bpdotimes (fdc (length nodes))
        (let ((n (aref nodes fdc)))
             (setf (cl-mpm/mesh::node-agg-fdc n) fdc))))

    (let ((nodes (cl-mpm/aggregate::sim-agg-nodes-fd sim)))
      (cl-mpm/utils::bpdotimes (fd (length nodes))
        (let ((n (aref nodes fd)))
          (setf (cl-mpm/mesh::node-agg-fd n) fd))))
    ;; (format t "FDC list ~D~%" (length (cl-mpm/aggregate::sim-agg-nodes-fd sim)))

    (dotimes (i (1+ refinement))
      (let ((mesh (aref mesh-list i)))
        (cl-mpm::iterate-over-nodes-array
         (make-array (array-total-size (cl-mpm/mesh::mesh-nodes mesh))
                     :displaced-to
                     (cl-mpm/mesh::mesh-nodes mesh))
         (lambda (n)
           (when (cl-mpm/mesh::node-active n)
             (when (and (cl-mpm/mesh::node-agg n)
                        (cl-mpm/mesh::node-interior n))
               (unless (cl-mpm/mesh::node-agg-fdc n)
                 (format t "Node ~D ~A agg to ~A has no assigned FDC~%"
                         i
                         (cl-mpm/mesh::node-index n)
                         (cl-mpm/mesh::cell-index (cl-mpm/mesh::node-agg-interior-cell n))))
               ))
           )
         )))
    (set-mesh-default sim)
    ))

(defun resolve-mesh-coupling (sim)
  )

(defun set-mesh-refinement (sim filter &key (derefine nil) (project-displacement nil))
  (declare (function filter))
  (set-mesh-default sim)
  (let ((mesh-updated nil)
        (spaced-coupling nil))
    (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                     (active-refinement sim-octree-active-refinement)
                     (refinement cl-mpm::sim-multigrid-refinement))
        sim
      (iterate-over-meshes
       sim
       (lambda (mesh mi)
         (cl-mpm::iterate-over-cells
          mesh
          (lambda (c)
            (setf (cell-octree-refine-prev c)
                  (cell-octree-refine c))))))
      (setf active-refinement 0)
      (loop for mesh-index from 0 to refinement
            do
               (progn
                 ;; (setup-mp-iteration-cache sim)
                 (let* ((mesh (aref mesh-list mesh-index)))
                   ;; (cl-mpm::iterate-over-cells
                   ;;  mesh
                   ;;  (lambda (c)
                   ;;    (when (and (cl-mpm/mesh::cell-active c))
                   ;;      (loop for cn in (cl-mpm/mesh::cell-neighbours c)
                   ;;            do (when (and (= (cell-octree-refine cn) 1))
                   ;;                 (when spaced-coupling
                   ;;                   (setf (cell-coupled-from-above c) t)))))))
                   (when (< mesh-index refinement)
                     ;;When a mesh needs to refine, we set it to 2
                     (cl-mpm::iterate-over-cells
                      mesh
                      (lambda (c)
                        (when t;(cl-mpm/mesh::cell-active c)
                          (let ((coupling-nearby nil))
                            (loop for cn in (cl-mpm/mesh::cell-neighbours c)
                                  do
                                     (when (and
                                            (cl-mpm/mesh::cell-active cn)
                                            (cell-coupled-from-above cn))
                                       (setf coupling-nearby t)))
                            (let ((current-refine (cell-octree-refine c))
                                  (to-refine (funcall filter sim mesh c)))

                              (when (or (= (cell-octree-refine c) -1)
                                        (= (cell-octree-refine c) 1)
                                        (cell-coupled-from-above c)
                                        coupling-nearby)
                                (setf to-refine nil))

                              (unless derefine
                                (when (= current-refine 2)
                                  (setf to-refine t)))


                              (when (and
                                     project-displacement
                                     to-refine
                                     (not (= current-refine 2)))
                                (iterate-over-sub-nodes
                                 sim
                                 mesh-index
                                 c
                                 (lambda (n)
                                   (cl-mpm/fastmaths::fast-zero (cl-mpm/mesh::node-displacment n))
                                   (cl-mpm::iterate-over-cell-shape-local
                                    mesh
                                    c
                                    (cl-mpm/mesh::node-position n)
                                    (lambda (nw w grads)
                                      (cl-mpm/fastmaths:fast-fmacc
                                       (cl-mpm/mesh::node-displacment n)
                                       (cl-mpm/mesh::node-displacment nw)
                                       w)))))
                                )
                              (when to-refine
                                (setf active-refinement
                                      (max active-refinement
                                           refinement
                                           (1+ mesh-index)))
                                (setf (cell-octree-refine c) 2))

                              (unless to-refine
                                (setf (cell-octree-refine c) 0))


                              )))))

                     ;;Next we flood-fill transition cells that are next to 2's to 1's
                     (cl-mpm::iterate-over-cells
                      mesh
                      (lambda (c)
                        (let ((touching-refine nil)
                              (touching-derefine nil))
                          (when (and t;(cl-mpm/mesh::cell-active c)
                                   (= (cell-octree-refine c) 0))
                          (loop for cn in (cl-mpm/mesh::cell-neighbours c)
                                do
                                  (when (and
                                         (not (cell-coupled-from-above c))
                                         (= (cell-octree-refine cn) 2))
                                    (when spaced-coupling
                                      (setf (cell-coupled-from-above c) t))
                                    (setf (cell-octree-refine c) 1))))))))

                 (when (< mesh-index refinement)
                   (cl-mpm::iterate-over-cells
                    ;(aref mesh-list (- mesh-index 1))
                    mesh
                    (lambda (cell)
                      ;; (when (= (cell-octree-refine cell) 2)
                      ;;   (setf (cl-mpm/mesh::cell-active cell) nil))
                      (iterate-over-sub-cells
                       sim
                       ;; (- mesh-index 1)
                       mesh-index
                       cell
                       (lambda (c)
                         (sb-thread::with-mutex ((cl-mpm/mesh::cell-lock c))
                           (setf (cell-coupled-from-above c) nil)
                           (ecase (cell-octree-refine cell)
                             (-1 (setf (cell-octree-refine c) -1))
                             (0 (setf (cell-octree-refine c) -1))
                             (1
                              (setf (cell-octree-refine c) 0)
                              (when spaced-coupling
                                (setf (cell-coupled-from-above c) t))
                              )
                             (2 (setf (cell-octree-refine c) 0)))))))))
                   (set-mesh-default sim))))
      (iterate-over-meshes
       sim
       (lambda (mesh mi)
         (cl-mpm::iterate-over-cells
          mesh
          (lambda (c)
            (when (and
                   (= (cell-octree-refine c) 2)
                   (not (= (cell-octree-refine-prev c)
                           (cell-octree-refine c))))
              (setf mesh-updated t)))))))

    ;; (setf (cl-mpm/mesh::mesh-active-nodes mesh) (build-active-node-set sim))
    ;(build-active-node-set sim)
    (build-full-node-set sim)
    mesh-updated))
(defun build-full-node-set (sim)
  (let ((nc 0))
    (iterate-over-meshes
     sim
     (lambda (mesh mesh-index)
       (let ((mesh-nodes (make-array (array-total-size (cl-mpm/mesh::mesh-nodes mesh)) :displaced-to (cl-mpm/mesh::mesh-nodes mesh))))
         (loop for n across mesh-nodes
               do (when t;(cl-mpm/mesh::node-active n)
                    (incf nc))))))
    (set-mesh-default sim)
    (let ((new-node-array (make-array nc)))
      (let ((i 0))
        (iterate-over-meshes
         sim
         (lambda (mesh mesh-index)
           (loop for n across (make-array (array-total-size (cl-mpm/mesh::mesh-nodes mesh)) :displaced-to (cl-mpm/mesh::mesh-nodes mesh))
                 do (progn
                      (when t;(cl-mpm/mesh::node-active n)
                        (setf (aref new-node-array i) n)
                        (incf i)))))))
      (set-mesh-default sim)
      (setf (cl-mpm/mesh::mesh-active-nodes (cl-mpm:sim-mesh sim))
            new-node-array)))
  )
(defun build-active-node-set (sim)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (active-refinement sim-octree-active-refinement))
      sim
    (labels ((iterate-over-active-meshes
                 (sim func)
               (loop for i from 0 to active-refinement
                     do (progn
                          (set-mesh sim i)
                          (funcall func (aref mesh-list i) i)))))
      (let ((nc 0))
        (iterate-over-active-meshes
         sim
         (lambda (mesh mesh-index)
           (let ((mesh-nodes (make-array (array-total-size (cl-mpm/mesh::mesh-nodes mesh)) :displaced-to (cl-mpm/mesh::mesh-nodes mesh))))
             (loop for n across mesh-nodes
                   do (when (cl-mpm/mesh::node-active n)
                        (incf nc))))))
        (set-mesh-default sim)
        (let ((new-node-array (make-array nc)))
          (let ((i 0))
            (iterate-over-active-meshes
             sim
             (lambda (mesh mesh-index)
               (loop for n across (make-array (array-total-size (cl-mpm/mesh::mesh-nodes mesh)) :displaced-to (cl-mpm/mesh::mesh-nodes mesh))
                     do (progn
                          (when (cl-mpm/mesh::node-active n)
                            (setf (aref new-node-array i) n)
                            (incf i)))))))
          (set-mesh-default sim)
          (setf (cl-mpm/mesh::mesh-active-nodes (cl-mpm:sim-mesh sim))
                new-node-array)))))
  )

(defun resolve-h (sim mp)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list))
      sim
    (cl-mpm/mesh::mesh-resolution (aref mesh-list (cl-mpm/particle::mp-mesh-index mp)))))
(defun resolve-node-h (sim node)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list))
      sim
    (cl-mpm/mesh::mesh-resolution (aref mesh-list (node-mesh-index node)))))

(defmethod cl-mpm::split-mps-eigenvalue ((sim mpm-sim-octree))
  (declare (optimize (speed 0) (debug 3) (safety 3)))
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (split-factor (cl-mpm::sim-split-factor sim))
         )
    (cl-mpm::split-mps-vector
     sim
     (lambda (mp)
       (let ((split-dir nil))
         (let ((td (cl-mpm/particle::mp-true-domain mp)))
           (multiple-value-bind (l v) (cl-mpm/utils:eig
                                       (cl-mpm/utils::slice-matrix-nd td nd))
             (let* ((v (cl-mpm/utils::pad-matrix-nd v nd))
                    (abs-l (mapcar #'abs l))
                    (max-l (reduce #'max abs-l))
                    (min-l (reduce #'min (remove 0d0 abs-l)))
                    (crit (* (resolve-h sim mp) split-factor))
                    )
               (when (> max-l crit)
                 (let* ((pos (position max-l abs-l))
                        (vec (magicl::vector->column-matrix (magicl:column v pos))))
                   (setf split-dir (cl-mpm/fastmaths:norm
                                    (cl-mpm/utils:vector-from-list (list (cl-mpm/utils::varef vec 0) (cl-mpm/utils::varef vec 1) (cl-mpm/utils::varef vec 2))))))))))
         split-dir)))))

(defmethod cl-mpm::split-mps-cartesian ((sim mpm-sim-octree))
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)
                   (max-split-depth cl-mpm::sim-max-split-depth)
                   (split-factor cl-mpm::sim-split-factor))
      sim
    (declare (fixnum max-split-depth))
    (let* ((nd (cl-mpm/mesh::mesh-nd mesh))
           (mps-to-split (remove-if-not (lambda (mp)
                                          (and
                                           (< (cl-mpm/particle::mp-split-depth mp) max-split-depth)
                                           (cl-mpm::split-criteria-variable
                                            mp
                                            (resolve-h sim mp)
                                            split-factor
                                            nd))) mps))
           (split-direction (map 'list (lambda (mp) (cl-mpm::split-criteria-variable
                                                     mp
                                                     (resolve-h sim mp)
                                                     split-factor
                                                     nd
                                                     )) mps-to-split))
           )
      (cl-mpm::remove-mps-func sim (lambda (mp) (cl-mpm::split-criteria-variable mp (resolve-h sim mp) split-factor
                                                                                 nd)))
      (loop for mp across mps-to-split
            for direction in split-direction
            do (loop for new-mp in (cl-mpm::split-mp mp (resolve-h sim mp) direction)
                     do (cl-mpm::sim-add-mp sim new-mp))))))



(defmethod cl-mpm::update-sim :before ((sim mpm-sim-octree-usf))
  (set-mesh-default sim)
  (set-mesh-refinement sim (sim-octree-refinement-criteria sim)
                       :derefine t)
  (cl-mpm/dynamic-relaxation::setup-mp-iteration-cache sim)
  (split-until-stable sim)
  (set-mesh-default sim))

(defun split-until-stable (sim)
  (with-accessors ((mps cl-mpm::sim-mps)
                   (split-enable cl-mpm::sim-allow-mp-split))
      sim
    (when split-enable
      (let* ((mp-count (length (cl-mpm:sim-mps sim)))
             (mp-count-prev 0)
             (mp-count-0 mp-count))
        (loop while (not (= mp-count mp-count-prev))
              do (progn
                   (cl-mpm::split-mps sim)
                   ;; (cl-mpm::split-mps sim)
                   ;; (cl-mpm/dynamic-relaxation::setup-mp-iteration-cache sim)
                   (setf mp-count-prev mp-count)
                   (setf mp-count (length (cl-mpm:sim-mps sim)))))
        (unless (= mp-count mp-count-0)
          (cl-mpm/dynamic-relaxation::setup-mp-iteration-cache sim))))))

(defmethod cl-mpm/dynamic-relaxation::pre-step :before ((sim mpm-sim-octree-quasi-static))
  (set-mesh-refinement sim (sim-octree-refinement-criteria sim))
  (cl-mpm/dynamic-relaxation::setup-mp-iteration-cache sim)
  (split-until-stable sim)
  (set-mesh-default sim))

(defmethod cl-mpm/dynamic-relaxation::pre-step :after ((sim mpm-sim-octree))
  (build-active-node-set sim)
  )

(defmethod cl-mpm/dynamic-relaxation::map-stiffness ((sim cl-mpm/dynamic-relaxation::mpm-sim-octree-quasi-static))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps))
      sim
    (set-mesh-default sim)
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (node)
       (declare (cl-mpm/mesh::node node))
       (when (cl-mpm/mesh:node-active node)
         (setf (the double-float (cl-mpm/mesh::node-mass node)) 0d0))))

    (let* ((nd (cl-mpm/mesh::mesh-nd mesh))
           (mass-scale (the double-float (/ 1d0 (the double-float (cl-mpm::sim-dt-scale sim))))))
      (declare (double-float mass-scale)
               (fixnum nd))
      (cl-mpm::iterate-over-mps
       mps
       (lambda (mp)
         (let* ((h (resolve-h sim mp))
                (mp-volume (cl-mpm/particle::mp-volume mp))
                (mp-pmod (cl-mpm/particle::estimate-stiffness mp))
                (ul (estimate-ul-enhancement mp nd))
                (mp-factor (* mp-pmod mp-volume ul mass-scale
                              (/ 1d0 (* h h)))))
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
                               (* 0.25d0 mp-factor)))))))))))))
  (set-mesh-default sim)
  )

;; (defmethod cl-mpm/setup::%estimate-elastic-dt ((sim mpm-sim-octree))
;;   (* (call-next-method) (expt 2 (- (cl-mpm::sim-multigrid-refinement sim)))))
(defmethod cl-mpm/setup::%estimate-elastic-dt-mps ((sim mpm-sim-octree))
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mass-scale cl-mpm::sim-mass-scale))
      sim
    (* ;; (sqrt mass-scale)
       (cl-mpm::reduce-over-mps
        mps
        (lambda (mp)
          (cl-mpm/setup::estimate-elastic-dt-mp
           sim
           (cl-mpm/particle::mp-p-modulus mp)
           (/ (the double-float (cl-mpm/particle:mp-mass mp))
              (the double-float (cl-mpm/particle:mp-volume mp)))
           (resolve-h sim mp)))
        #'min))))

(in-package :cl-mpm/output)
(defmethod save-vtk-cells (filename (sim cl-mpm/dynamic-relaxation::mpm-sim-octree))
  (with-accessors ((mesh cl-mpm:sim-mesh)) sim
    (with-accessors ((cells cl-mpm/mesh::mesh-cells))
        mesh
        (with-open-file (fs filename :direction :output :if-exists :supersede)
          (format fs "# vtk DataFile Version 2.0~%")
          (format fs "Lisp generated vtk file, SJVS~%")
          (format fs "ASCII~%")
          (format fs "DATASET UNSTRUCTURED_GRID~%")

          (let (;; (node-count 0)
                (cells (remove-if-not (lambda (cell)
                                        t
                                        ;; (or
                                        ;;  (cl-mpm/mesh::cell-active cell)
                                        ;;  (and
                                        ;;   (cl-mpm/mesh::cell-mp-count cell)
                                        ;;   (= (cl-mpm/mesh::cell-mp-count cell) 1)))
                                        )
                                      (make-array (array-total-size cells) :displaced-to cells)))
                )
            ;; (cl-mpm::iterate-over-cells-serial
            ;;  mesh
            ;;  (lambda (n)
            ;;    (declare (ignore n))
            ;;    (incf node-count)))
            (format fs "POINTS ~d double~%" (array-total-size cells))
            (loop for i from 0 below (array-total-size cells)
                  do
                     (let ((cell (row-major-aref cells i)))
                       (when cell
                         ;; (incf node-count)
                         (let ((pos (cl-mpm/mesh::cell-centroid cell)))
                           (format fs "~E ~E ~E ~%"
                                   (coerce (cl-mpm/utils::varef pos 0) 'single-float)
                                   (coerce (cl-mpm/utils::varef pos 1) 'single-float)
                                   (coerce (cl-mpm/utils::varef pos 2) 'single-float))))))
            (format fs "~%")
            (let ((id 1))
              (declare (special id))
              (format fs "POINT_DATA ~d~%" (array-total-size cells))
              (save-parameter-cells "octree" (cl-mpm/dynamic-relaxation::cell-octree-refine cell))
              (save-parameter-cells "coupling" (if (cl-mpm/dynamic-relaxation::cell-coupling cell) 1 0))
              (save-parameter-cells "coupled-above" (if (cl-mpm/dynamic-relaxation::cell-coupled-from-above cell) 1 0))
              (save-parameter-cells "reaggregate" (if (cl-mpm/dynamic-relaxation::cell-reaggregation-element cell) 1 0))
              (save-parameter-cells "buoyancy" (if (cl-mpm/mesh::cell-boundary cell) 1 0))
              (save-parameter-cells "active" (if (cl-mpm/mesh::cell-active cell) 1 0))
              (save-parameter-cells "partial" (if (cl-mpm/mesh::cell-partial cell) 1 0))
              (save-parameter-cells "ghost" (if (cl-mpm/mesh::cell-ghost-element cell) 1 0))
              (save-parameter-cells "pressure" (cl-mpm/mesh::cell-pressure cell))
              (save-parameter-cells "cell-count" (cl-mpm/mesh::cell-mp-count cell))
              (save-parameter-cells "agg-int" (if (cl-mpm/mesh::cell-interior cell) 1d0 0d0))
              (save-parameter-cells "agg" (if (cl-mpm/mesh::cell-agg cell) 1d0 0d0))
              (save-parameter-cells "def-0" (if (cl-mpm/mesh::cell-def-list cell) (nth 0 (cl-mpm/mesh::cell-def-list cell)) 0))
              (save-parameter-cells "def-1" (if (cl-mpm/mesh::cell-def-list cell) (nth 1 (cl-mpm/mesh::cell-def-list cell)) 0))
              ;; (save-parameter-cells "def-2" (if (cl-mpm/mesh::cell-def-list cell) (nth 2 (cl-mpm/mesh::cell-def-list cell)) 0))
              (save-parameter-cells "j" (magicl:det (cl-mpm/mesh::cell-deformation-gradient cell)))
              ))))))

(defmethod save-vtk-nodes (filename (sim cl-mpm/dynamic-relaxation::mpm-sim-octree))
  (let ((nc 0))
    (loop for mesh across (cl-mpm::sim-mesh-list sim)
          do (incf nc (cl-mpm::reduce-over-nodes
                       mesh
                       (lambda (n) (if (cl-mpm/mesh::node-active n) 1 0))
                       #'+)))
    (let ((nodes (make-array nc))
          (iter 0))
      (loop for mesh across (cl-mpm::sim-mesh-list sim)
            do (cl-mpm::iterate-over-nodes-serial
                mesh
                (lambda (n)
                  (when (cl-mpm/mesh::node-active n)
                    (setf (aref nodes iter) n)
                    (incf iter)))))
      (with-open-file (fs filename :direction :output :if-exists :supersede)
        (format fs "# vtk DataFile Version 2.0~%")
        (format fs "Lisp generated vtk file, SJVS~%")
        (format fs "BINARY~%")
        (format fs "DATASET UNSTRUCTURED_GRID~%"))
      (with-open-file (fs filename :direction :output :if-exists :append)
        (with-open-file (fs-bin filename :direction :output :if-exists :append :element-type '(unsigned-byte 8))
          (let* ((node-count iter)
                 ;; (nodes (remove-if-not 'cl-mpm/mesh:node-active
                 ;;                       (make-array node-count :displaced-to full-nodes)))
                 )
            (when (= (length nodes) 0)
              (setf nodes (make-array node-count :displaced-to nodes)))
            (let ((node-count (array-total-size nodes)))
              (format fs "POINTS ~d double~%" (array-total-size nodes))
              (force-output fs)
              (loop for i from 0 below (array-total-size nodes)
                    do
                       (let ((node (row-major-aref nodes i)))
                         (when node
                           (let ((pos (cl-mpm/fastmaths:fast-.+
                                       (cl-mpm/mesh::node-position node)
                                       (cl-mpm/mesh::node-displacment node))))
                             (write-binary-float (cl-mpm/utils:varef pos 0) fs-bin)
                             (write-binary-float (cl-mpm/utils:varef pos 1) fs-bin)
                             (write-binary-float (cl-mpm/utils:varef pos 2) fs-bin)))))
              (force-output fs-bin)
              (format fs "~%")
              (let ((id 1)
                    (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim))))
                (declare (special id))
                (format fs "POINT_DATA ~d~%" (array-total-size nodes))

                (force-output fs)
                (dolist (f (cl-mpm::sim-output-list-nodes sim))
                  (destructuring-bind (type name accessor) f
                    (case type
                      (:BOOL
                       (cl-mpm/output::save-parameter-nodes name (if (funcall accessor node) 1d0 0d0)))
                      (:SCALAR
                       (cl-mpm/output::save-parameter-nodes name (funcall accessor node)))
                      (:VECTOR
                       (cl-mpm/output::save-parameter-nodes (format nil "~A_x" name) (varef (funcall accessor node) 0))
                       (cl-mpm/output::save-parameter-nodes (format nil "~A_y" name) (varef (funcall accessor node) 1))
                       (when (= nd 3)
                         (cl-mpm/output::save-parameter-nodes (format nil "~A_z" name) (varef (funcall accessor node) 2)))))))))))))))

(in-package :cl-mpm/dynamic-relaxation)



(defmethod cl-mpm/damage::calculate-damage :around ((sim mpm-sim-octree-damage) dt)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps)
                   (mesh-list cl-mpm::sim-mesh-list))
      sim
    (set-mesh-default sim)
    (call-next-method)
    (set-mesh-default sim)))

(defmethod cl-mpm::sim-add-mp :around ((sim mpm-sim-octree) mp)
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (mesh-list cl-mpm::sim-mesh-list)
                   )
      sim
    (set-mesh-default sim)
    (call-next-method)
    (set-mesh-default sim)))

(defmethod cl-mpm::calculate-min-dt-mps ((sim mpm-sim-octree))
  "Estimate minimum p-wave modulus"
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mass-scale cl-mpm::sim-mass-scale))
      sim
    (set-mesh-default sim)
    (let ((inner-factor most-positive-double-float))
      (declare (double-float inner-factor mass-scale))
      (setf inner-factor
            (cl-mpm::reduce-over-nodes
             mesh
             (lambda (node)
               (if (and (cl-mpm/mesh::node-active node)
                        (or (not (cl-mpm/mesh::node-agg node))
                            (cl-mpm/mesh::node-interior node)))
                   (with-accessors ((node-active  cl-mpm/mesh:node-active)
                                    (pmod cl-mpm/mesh::node-pwave)
                                    (mass cl-mpm/mesh::node-mass)
                                    (svp-sum cl-mpm/mesh::node-svp-sum)
                                    (vol cl-mpm/mesh::node-volume)
                                    (mesh-index node-mesh-index)
                                    (vel cl-mpm/mesh::node-velocity)
                                    ) node
                     (declare (double-float pmod mass svp-sum vol)
                              (fixnum mesh-index))
                     (if (and (> vol 0d0)
                              (> pmod 0d0)
                              (> svp-sum 0d0))
                         (let* ((h (resolve-node-h sim node))
                                (nf (* (/ mass pmod))))
                           (* h (sqrt nf)))
                         most-positive-double-float))
                   most-positive-double-float))
             #'min))
      (if (< inner-factor most-positive-double-float)
          (* (sqrt mass-scale) inner-factor)
          (cl-mpm:sim-dt sim)))))


(defun save-all-vtks (sim output-dir step)

  (cl-mpm:sim-format sim t "Save vtks ~D~%" step)
  (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" step)) sim)
  (dotimes (i (1+ (cl-mpm::sim-multigrid-refinement sim)))
    (setf (cl-mpm::sim-mesh sim) (aref (cl-mpm::sim-mesh-list sim) i))
    (cl-mpm/output::save-vtk-mesh-nodes (merge-pathnames output-dir (format nil "sim_nodes_~D_~5,'0d.vtk" i step)) (aref (cl-mpm::sim-mesh-list sim) i))
    (cl-mpm/output::save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~D_~5,'0d.vtk" i step)) sim))
  (setf (cl-mpm::sim-mesh sim) (aref (cl-mpm::sim-mesh-list sim) 0))
  (set-mesh-default sim)
  )

(defmethod cl-mpm/dynamic-relaxation::save-vtks ((sim cl-mpm/dynamic-relaxation::mpm-sim-octree) output-dir step)
  (cl-mpm:sim-format sim t "Save vtks ~D~%" step)
  (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" step)) sim)
  (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) sim)
  (cl-mpm/output::save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d.vtk" step)) sim)
  (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim_p_~5,'0d.vtk" step)) sim )
  (dotimes (i (1+ (cl-mpm::sim-multigrid-refinement sim)))
    (setf (cl-mpm::sim-mesh sim) (aref (cl-mpm::sim-mesh-list sim) i))
    (cl-mpm/output::save-vtk-mesh-nodes (merge-pathnames output-dir (format nil "sim_nodes_~D_~5,'0d.vtk" i step)) (aref (cl-mpm::sim-mesh-list sim) i) sim)
    (cl-mpm/output::save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~D_~5,'0d.vtk" i step)) sim))
  ;; (setf (cl-mpm::sim-mesh sim) (aref (cl-mpm::sim-mesh-list sim) 0))
  (set-mesh-default sim)
  )

(defmethod cl-mpm/dynamic-relaxation::save-vtks-dr-step ((sim cl-mpm/dynamic-relaxation::mpm-sim-octree) output-dir step iter)
  (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_step_~5,'0d_~5,'0d.vtk" step iter)) sim)
  (dotimes (i (1+ (cl-mpm::sim-multigrid-refinement sim)))
    (setf (cl-mpm::sim-mesh sim) (aref (cl-mpm::sim-mesh-list sim) i))
    (cl-mpm/output::save-vtk-mesh-nodes (merge-pathnames output-dir (format nil "sim_step_nodes_~D_~5,'0d_~5,'0d.vtk" i step iter)) (aref (cl-mpm::sim-mesh-list sim) i))
    (cl-mpm/output::save-vtk-cells (merge-pathnames output-dir (format nil "sim_step_cells_~D_~5,'0d_~5,'0d.vtk" i step iter)) sim))
  (set-mesh-default sim)
 )

;; (let* ((est-size 1000)
;;        ;; (v (make-array est-size :fill-pointer 0 :adjustable t :element-type 'double-float))
;;        ;; (r (make-array est-size :fill-pointer 0 :adjustable t :element-type 'fixnum))
;;        ;; (c (make-array est-size :fill-pointer 0 :adjustable t :element-type 'fixnum))
;;        (counter (make-array 1 :element-type '(UNSIGNED-BYTE 64)))
;;        )
;;   (defun assemble-sparse-multigrid-e (sim)
;;     (let* ((agg-nodes (sim-agg-nodes-fdc sim))
;;            (active-nodes (sim-agg-nodes-fd sim))
;;            (mesh (cl-mpm:sim-mesh sim))
;;            (nd (cl-mpm/mesh:mesh-nd mesh))
;;            (ndof (length active-nodes))
;;            (ndofC (length agg-nodes))
;;            (est-size (+ ndofC (* (expt 2 nd) (- ndof ndofC))))
;;            (lock (sb-thread:make-mutex))
;;            (v (make-array est-size :fill-pointer est-size :adjustable t :element-type 'double-float))
;;            (r (make-array est-size :fill-pointer est-size :adjustable t :element-type 'fixnum))
;;            (c (make-array est-size :fill-pointer est-size :adjustable t :element-type 'fixnum)))
;;       ;; (setf *assembly-counter* 0)
;;       (sb-ext:atomic-update (aref counter 0) (lambda (a) 0))
;;       ;; (setf v (adjust-array v est-size :fill-pointer 0))
;;       ;; (setf r (adjust-array r est-size :fill-pointer 0))
;;       ;; (setf c (adjust-array c est-size :fill-pointer 0))
;;       ;; (setf 
;;       ;;  (fill-pointer v) 0
;;       ;;  (fill-pointer r) 0
;;       ;;  (fill-pointer c) 0)
;;       (cl-mpm::iterate-over-nodes
;;        mesh
;;        (lambda (node)
;;          (when (cl-mpm/mesh:node-active node)
;;            (when (cl-mpm/mesh::node-agg node)
;;              ;; (sb-thread:with-mutex (lock))
;;              (if (cl-mpm/mesh::node-interior node)
;;                  (progn
;;                    ;;Interior -> populate 1s
;;                    (sb-thread:with-mutex (lock)
;;                      ;let ((place (sb-ext:atomic-incf (aref counter 0))))
;;                      ;; (setf (aref v place) 1d0)
;;                      ;; (setf (aref r place) (cl-mpm/mesh::node-agg-fd node))
;;                      ;; (setf (aref c place) (cl-mpm/mesh::node-agg-fdc node))
;;                      (vector-push-extend  1d0 v)
;;                      (vector-push-extend  (cl-mpm/mesh::node-agg-fd node) r)
;;                      (vector-push-extend  (cl-mpm/mesh::node-agg-fdc node) c)
;;                      ))
;;                  (progn
;;                    ;;Aggregate -> populate svp
;;                    (labels ((apply-svps (cell)
;;                               (cl-mpm::iterate-over-cell-shape-local
;;                                mesh
;;                                cell
;;                                (cl-mpm/mesh::node-position node)
;;                                (lambda (cn weight grads)
;;                                  (unless (cl-mpm/mesh::node-agg-fdc cn)
;;                                    (apply-svps (cl-mpm/mesh::node-agg-interior-cell cn)))
;;                                  (when (and (cl-mpm/mesh::node-active cn)
;;                                             (cl-mpm/mesh::node-agg-fdc cn)
;;                                             (> (abs weight) 1d-9))
;;                                    (sb-thread:with-mutex (lock)
;;                                         ;let ((place (sb-ext:atomic-incf (aref counter 0))))
;;                                      ;; (setf (aref v place) weight)
;;                                      ;; (setf (aref r place) (cl-mpm/mesh::node-agg-fd node))
;;                                      ;; (setf (aref c place) (cl-mpm/mesh::node-agg-fdc cn))
;;                                      (vector-push-extend  weight v)
;;                                      (vector-push-extend  (cl-mpm/mesh::node-agg-fd node) r)
;;                                      (vector-push-extend  (cl-mpm/mesh::node-agg-fdc cn) c)
;;                                      )
;;                                    )))))
;;                      (apply-svps 
;;                       (cl-mpm/mesh::node-agg-interior-cell node))
;;                      )))))))
;;       (values (cl-mpm/utils::build-sparse-matrix v r c ndof ndofC)
;;               (cl-mpm/utils::build-sparse-matrix v c r ndofC ndof)))))

(defclass mpm-sim-octree-implicit-dynamic (mpm-sim-implict-dynamic mpm-sim-octree-damage-quasi-static)
  ())

(defmethod cl-mpm/dynamic-relaxation::convergence-check ((sim cl-mpm/dynamic-relaxation::mpm-sim-octree))
  (let ((refinement-check
          (cl-mpm/dynamic-relaxation::set-mesh-refinement
           sim
           (cl-mpm/dynamic-relaxation::sim-octree-refinement-criteria sim)
           :derefine nil
           :project-displacement t)))
    (when refinement-check
      (with-accessors ((mesh cl-mpm:sim-mesh)
                       (mps cl-mpm:sim-mps)
                       (vel-algo cl-mpm::sim-velocity-algorithm)
                       (mass-filter cl-mpm::sim-mass-filter)
                       (ghost-factor cl-mpm::sim-ghost-factor))
          sim
        (cl-mpm/dynamic-relaxation::setup-mp-iteration-cache sim)
        (split-until-stable sim)
        (set-mesh-default sim)
        (cl-mpm::p2g mesh mps vel-algo)
        (setf (cl-mpm::sim-dt sim) 1d0)
        (when (> mass-filter 0d0)
          (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
        (cl-mpm::apply-essential-bcs sim)
        (cl-mpm::filter-cells sim)
        (cl-mpm::apply-essential-bcs sim)
        (cl-mpm::iterate-over-nodes
         mesh
         (lambda (n)
           (when (cl-mpm/mesh::node-active n)
             (setf
              (cl-mpm/mesh::node-true-mass n) (cl-mpm/mesh:node-mass n)))
           (cl-mpm/fastmaths:fast-zero (cl-mpm/mesh::node-true-velocity n)))))
      (cl-mpm::reset-nodes-force sim)
      ;; (midpoint-starter sim)

      (format t "Refining mesh inside convergence loop ~%"))
    (not refinement-check)))

(defmethod refine-mesh ((sim cl-mpm/dynamic-relaxation::mpm-sim-octree))
  (cl-mpm/dynamic-relaxation::set-mesh-refinement
   sim
   (cl-mpm/dynamic-relaxation::sim-octree-refinement-criteria sim)
   :derefine nil
   :project-displacement t))
