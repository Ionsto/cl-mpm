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
;;       (setf mesh m))))

;; (defmethod cl-mpm::sim-add-mp :around ((sim mpm-sim-dr-multigrid) mp)
;;   (with-accessors ((mesh cl-mpm::sim-mesh)
;;                    (mesh-list cl-mpm::sim-mesh-list)
;;                    )
;;       sim
;;     (let ((m (cl-mpm:sim-mesh sim)))
;;       (setf mesh (aref mesh-list (- (length mesh-list) 1)))
;;       (call-next-method)
;;       (setf mesh m))))

;; ;; (defmethod pre-step ((sim mpm-sim-dr-multigrid))
;; ;;   (with-accessors ((mesh cl-mpm:sim-mesh)
;; ;;                    (mps cl-mpm:sim-mps)
;; ;;                    (delocal-counter cl-mpm/damage::sim-damage-delocal-counter-max))
;; ;;       sim
;; ;;     (call-next-method)
;; ;;     (setf delocal-counter -1)
;; ;;     (progn
;; ;;       (cl-mpm/damage::update-delocalisation-list (first (last (cl-mpm::sim-mesh-list sim))) mps))))

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

(defun set-mesh-default (sim)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (mesh cl-mpm::sim-mesh)
                   )
      sim
      (setf (cl-mpm:sim-mesh sim) (aref mesh-list 0))))
(defun set-mesh (sim index)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (mesh cl-mpm::sim-mesh))
      sim
    (setf (cl-mpm:sim-mesh sim) (aref mesh-list index))))


(defclass cell-octree (cl-mpm/mesh::cell)
  ((octree-refine
    :initform 0
    :accessor cell-octree-refine)
   (mesh-index
    :initform 0
    :accessor cell-mesh-index)
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
   (intra-mesh-aggregation
    :initform t
    :initarg :intra-mesh-agg
    :accessor sim-intra-mesh-aggregation)))


(defclass mpm-sim-octree-usf (cl-mpm/aggregate::mpm-sim-agg-usf mpm-sim-octree )
  ())

(defclass mpm-sim-octree-damage-usf (cl-mpm/damage::mpm-sim-agg-damage mpm-sim-octree-usf)
  ())

(defclass mpm-sim-octree-quasi-static (cl-mpm/dynamic-relaxation::mpm-sim-dr-ul mpm-sim-octree)
  ())

(defclass mpm-sim-octree-damage-quasi-static (cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul mpm-sim-octree-quasi-static)
  ())


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
                              (break)
                              (setf node (cl-mpm/mesh:get-node mesh index))
                              (if node
                                  (cl-mpm/bc:apply-bc bc node mesh dt)
                                  (error "BC attempted to get a nil node ~A ~A" bc index)))))))))))))


(defmethod cl-mpm/setup::%setup-bcs ((sim cl-mpm/dynamic-relaxation::mpm-sim-octree)
                                     left
                                     right
                                     top
                                     bottom
                                     front
                                     back)
  (with-accessors ((bcs-list cl-mpm::sim-bcs-list )
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
                        mesh
                        left right top bottom front back))
                 ;;Avoid triggering the setf bcs active rebuild
                 (setf
                  (cl-mpm::sim-bcs sim)
                  (aref bcs-list i)
                  )
                 (cl-mpm/setup::resolve-bc-nodes sim mesh (aref bcs-list i))
                 (incf bc-count (length (aref bcs-list i))))))
    ;; (setf
    ;;  (cl-mpm::sim-bcs sim)
    ;;  (aref bcs-list 0))
    ))



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
          (let* ((m (cl-mpm::make-mesh size (/ resolution (expt 2 (+ refine 0))) nil))
                 (bcs (cl-mpm/bc:make-outside-bc m)))
            (setf (aref mesh-list refine) m)
            (setf (aref bcs-list refine) bcs)))

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

          sim))))



(defun setup-mp-iteration-cache (sim)
  (set-mesh-default sim)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (refinement cl-mpm::sim-multigrid-refinement))
      sim
      (cl-mpm::iterate-over-mps
       (cl-mpm:sim-mps sim)
       (lambda (mp)
         (when (> refinement 0)
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

                            ;; (if (= 2 (cell-octree-refine (get-cell-at-pos (cl-mpm/particle::mp-position mp))))
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
                              (if (and (>= max-index 1)
                                       (> min-index 0))
                                  (progn
                                    (incf resolved-mesh))
                                  (loop-finish))))))

               ;; (setf resolved-mesh 0)
               ;; (pprint resolved-mesh)

               (cl-mpm/particle::reset-mp-node-cache mp)
               (setf (cl-mpm/particle::mp-index mp) resolved-mesh)
               (setf (cl-mpm/particle::mp-mesh-index mp) resolved-mesh)
               (cl-mpm::iterate-over-neighbours
                mesh
                mp
                (lambda (&rest args))))))))
    (set-mesh-default sim)))

(defun intra-mesh-agg (sim mesh)
  (with-accessors ()
    sim
    (setf (cl-mpm::sim-mesh sim) mesh)
    (let ((mesh-nodes
            (make-array (array-total-size (cl-mpm/mesh::mesh-nodes mesh))
                        :displaced-to
                        (cl-mpm/mesh::mesh-nodes mesh))))

      (cl-mpm::iterate-over-nodes-array
       mesh-nodes
       (lambda (n)
         (setf (cl-mpm/mesh::node-agg-building-flag n) nil)))
      ;; (cl-mpm::iterate-over-cells
      ;;  mesh
      ;;  (lambda (c)
      ;;    (setf (cl-mpm/mesh::cell-agg c) nil
      ;;          (cl-mpm/mesh::cell-interior c) nil)))


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
                        (if (cl-mpm/mesh::cell-interior n)
                            (setf (cell-reaggregation-element n) t)
                            (setf (cl-mpm/mesh::cell-agg n) t))))))))

      ;;Next set any nodes on agg elements to be agg
      (cl-mpm::iterate-over-cells
       mesh
       (lambda (cell)
         (with-accessors ((nodes cl-mpm/mesh::cell-nodes)
                          (agg cl-mpm/mesh::cell-agg)
                          (active cl-mpm/mesh::cell-active))
             cell
           (when (and active agg (not (cell-reaggregation-element cell)))
             (loop for n across nodes
                   do (when (cl-mpm/mesh::node-active n)
                        (setf (cl-mpm/mesh::node-agg n) t
                              (cl-mpm/mesh::node-agg-building-flag n) t))))
           )))

      ;;For each aggregate node, locate the closest support cell
      (cl-mpm::iterate-over-nodes-array
       mesh-nodes
       (lambda (node)
         (with-accessors ((active cl-mpm/mesh::node-active)
                          (agg cl-mpm/mesh::node-agg)
                          (int cl-mpm/mesh::node-interior)
                          (building cl-mpm/mesh::node-agg-building-flag)
                          (int-cell cl-mpm/mesh::node-agg-interior-cell))
             node
           (when (and active agg ; (not int-cell)
                      )
             (if int-cell
                 (progn
                   ;;If we've already got aggregation happening, we need to reaggregate
                   ;; (setf (cell-reaggregation-element int-cell) closest-elem)
                   ;; (loop for n across (cl-mpm/mesh::cell-nodes int-cell)
                   ;;       do (setf (cl-mpm/mesh::node-interior n) nil))
                   ;; (setf (cl-mpm/mesh::cell-interior int-cell) nil)
                   ;; (setf (cl-mpm/mesh::cell-interior closest-elem) t)
                   (setf int nil))
                 (let ((closest-elem (cl-mpm/aggregate::get-closest-cell
                                      mesh
                                      (cl-mpm/mesh::node-position node)
                                      :exclude int-cell
                                      :filter (lambda (c)
                                                (not (= (cell-octree-refine c) -1))))))
                   (if closest-elem
                       (progn
                         (setf int-cell closest-elem)
                         (setf (cl-mpm/mesh::cell-interior closest-elem) t))
                       (error "No closest elem?")))))))))
  (cl-mpm::iterate-over-cells
   mesh
   (lambda (cell)
     (with-accessors ((nodes cl-mpm/mesh::cell-nodes)
                      (agg cl-mpm/mesh::cell-agg)
                      (int cl-mpm/mesh::cell-interior)
                      (active cl-mpm/mesh::cell-active)
                      (reagg cell-reaggregation-element)
                      )
         cell
       (when (and active int (not reagg))
         (loop for n across nodes
               do (when (and (cl-mpm/mesh::node-active n)
                             ;; (cl-mpm/mesh::node-agg-building-flag n)
                             )
                    (setf (cl-mpm/mesh::node-agg n) t
                          (cl-mpm/mesh::node-agg-interior-cell n) cell
                          (cl-mpm/mesh::node-interior n) t))))

       (when reagg
         (setf int nil)
         (setf
          (cell-reaggregation-element cell)
          (cl-mpm/aggregate::get-closest-cell
           mesh
           (cl-mpm/mesh::cell-centroid cell)
           :exclude cell
           :filter
           (lambda (c)
             (not (= (cell-octree-refine c) -1))))))
       )))))

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
                 (setf (cl-mpm/mesh::cell-agg cell) t))

               (iterate-over-sub-cells
                sim
                mesh-index
                cell
                (lambda (c)
                  (setf (cl-mpm/mesh::cell-interior c) nil)
                  (when (cl-mpm/mesh::cell-agg-int c)
                    (setf (cell-reaggregation-element c) cell))))

               (let ((sub-node-active nil))
                 (iterate-over-sub-nodes
                  sim
                  mesh-index
                  cell
                  (lambda (n)
                    (when (cl-mpm/mesh::node-active n)
                      (setf sub-node-active t)
                      (when t;(not (cl-mpm/mesh::node-agg n))
                        (setf
                         (cl-mpm/mesh::node-agg n) t
                         (cl-mpm/mesh::node-interior n) nil
                         (cl-mpm/mesh::node-agg-interior-cell n) cell
                         ))
                      (cl-mpm::iterate-over-cell-shape-local
                       mesh
                       cell
                       (cl-mpm/mesh::node-position n)
                       (lambda (n w grads)
                         (when (not (= w 0d0))
                           (setf (cl-mpm/mesh::node-active n) t)))))))
                 (when sub-node-active
                   (loop for n across nodes
                         do (progn
                              (when (cl-mpm/mesh::node-active n)
                                (setf (cl-mpm/mesh::node-agg n) t
                                      (cl-mpm/mesh::node-interior n) t)
                                (unless (cl-mpm/mesh::node-agg-interior-cell n)
                                  (setf (cl-mpm/mesh::node-agg-interior-cell n) cell)))
                              ;; (setf (cl-mpm/mesh::node-active n) t)
                              ))))))))

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

(defun iterate-over-sub-cells-2d (sim mesh-index cell func)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (dt cl-mpm::sim-dt)
                   (agg cl-mpm/aggregate::sim-enable-aggregate))
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

(defun iterate-over-sub-cells (sim mesh-index cell func)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (dt cl-mpm::sim-dt)
                   (agg cl-mpm/aggregate::sim-enable-aggregate))
      sim
    (let ((mesh (aref mesh-list mesh-index)))
      (ecase (cl-mpm/mesh::mesh-nd mesh)
        (1 (error "Not implemented"))
        (2 (iterate-over-sub-cells-2d sim mesh-index cell func))
        (3 (error "Not implemented"))))))

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

(defun iterate-over-sub-nodes (sim mesh-index cell func)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                   (dt cl-mpm::sim-dt)
                   (agg cl-mpm/aggregate::sim-enable-aggregate))
      sim
    (let ((mesh (aref mesh-list mesh-index)))
      (ecase (cl-mpm/mesh::mesh-nd mesh)
        (1 (error "Not implemented"))
        (2 (iterate-over-sub-nodes-2d sim mesh-index cell func))
        (3 (error "Not implemented"))))))

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
                                   (partial cl-mpm/mesh::cell-partial))
                      cell
                    (when (= (cell-octree-refine cell) 2)
                      (setf (cl-mpm/mesh::cell-active cell) nil))
                    (when t;active
                      (iterate-over-sub-cells
                       sim
                       i
                       cell
                       (lambda (c)
                         (setf (cl-mpm/mesh::cell-active c)
                               (or (cl-mpm/mesh::cell-active c)
                                   active))
                         (unless partial
                           (setf (cl-mpm/mesh::cell-partial c) nil))
                         ))))))))
    ;;Next we check cells that do not know if they are full by gathering that infomation recursivly
    (loop for i from 0 below refinement
          do
             (let ((mesh (aref mesh-list i)))
               (cl-mpm::iterate-over-cells
                mesh
                (lambda (cell)
                  (with-accessors ((active cl-mpm/mesh::cell-active)
                                   (partial cl-mpm/mesh::cell-partial))
                      cell
                    (when t;active
                      (let ((trial-partial nil))
                        (iterate-over-sub-cells
                         sim
                         i
                         cell
                         (lambda (c)
                           (setf trial-partial (or trial-partial (cl-mpm/mesh::cell-partial c)))))
                        (setf partial trial-partial))))))))
    ))

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
                      (cell-reaggregation-element node) nil
                      (cl-mpm/mesh::cell-interior node) nil))))


    ;;First set all outside nodes as aggregate

    (when (> refinement 0)
      ;;Build the inter-mesh aggregation
      (loop for mesh-index from (- refinement 1) downto 0
            do
               (progn
                 (when intra-mesh-agg
                   (intra-mesh-agg sim (aref mesh-list (1+ mesh-index))))
                 (set-mesh sim mesh-index)
                 (inter-mesh-agg sim mesh-index))))
    (when intra-mesh-agg
      (intra-mesh-agg sim (aref mesh-list 0)))

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
    (setf mesh (aref mesh-list 0))
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
    ))

(defun set-mesh-refinement (sim filter)
  (declare (function filter))
  (set-mesh-default sim)
  (let ((mesh-updated nil))
    (with-accessors ((mesh-list cl-mpm::sim-mesh-list)
                     (refinement cl-mpm::sim-multigrid-refinement))
        sim
      (loop for mesh-index from 0 below (- refinement 0)
            do
               (let* ((mesh (aref mesh-list mesh-index)))
                 (set-mesh sim mesh-index)
                 (cl-mpm::iterate-over-cells
                  mesh
                  (lambda (c)
                    (let ((current-refine (cell-octree-refine c))
                          (to-refine (funcall filter mesh c)))
                      (when (and to-refine
                                 (not (= (cell-octree-refine c) 2)))
                        (setf (cell-octree-refine c) 2)
                        (loop for cn in (cl-mpm/mesh::cell-neighbours c)
                              do (setf (cell-octree-refine cn)
                                       (max
                                        (cell-octree-refine cn)
                                        1))))
                      (when (not (= current-refine (cell-octree-refine c)))
                        (setf mesh-updated t)))))
                 (when (> (cl-mpm::sim-multigrid-refinement sim) 0)
                   (cl-mpm::iterate-over-cells
                    mesh
                    (lambda (cell)
                      (when (= (cell-octree-refine cell) 2)
                        (setf (cl-mpm/mesh::cell-active cell) nil))
                      (iterate-over-sub-cells
                       sim
                       mesh-index
                       cell
                       (lambda (c)
                         (ecase (cell-octree-refine cell)
                           (-1 nil)
                           (0 (setf (cell-octree-refine c) -1))
                           (1 (setf (cell-octree-refine c) 0))
                           (2 (setf (cell-octree-refine c) 0))))))))
                 (set-mesh-default sim))))
    mesh-updated))

(defun resolve-h (sim mp)
  (with-accessors ((mesh-list cl-mpm::sim-mesh-list))
      sim
    (cl-mpm/mesh::mesh-resolution (aref mesh-list (cl-mpm/particle::mp-mesh-index mp)))))

(defmethod cl-mpm::split-mps-eigenvalue ((sim mpm-sim-octree))
  (declare (optimize (speed 0) (debug 3) (safety 3)))
  (let* ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
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




(defmethod cl-mpm::update-sim :before ((sim mpm-sim-octree-usf))
  (set-mesh-default sim)
  (set-mesh-refinement sim (sim-octree-refinement-criteria sim))
  (cl-mpm/dynamic-relaxation::setup-mp-iteration-cache sim)
  (split-until-stable sim)
  (set-mesh-default sim))

(defun split-until-stable (sim)
  (with-accessors ((mps cl-mpm::sim-mps))
      sim
    (let* ((mp-count (length (cl-mpm:sim-mps sim)))
           (mp-count-prev 0)
           (mp-count-0 mp-count))
      (loop while (not (= mp-count mp-count-prev))
            do (progn
                 (cl-mpm::split-mps sim)
                 (setf mp-count-prev mp-count
                       mp-count (length (cl-mpm:sim-mps sim)))))
      (unless (= mp-count mp-count-0)
        (cl-mpm/dynamic-relaxation::setup-mp-iteration-cache sim)))))

(defmethod cl-mpm/dynamic-relaxation::pre-step :before ((sim mpm-sim-octree-quasi-static))
  (set-mesh-refinement sim (sim-octree-refinement-criteria sim))
  (cl-mpm/dynamic-relaxation::setup-mp-iteration-cache sim)
  (split-until-stable sim)
  (set-mesh-default sim)
  )

(defmethod cl-mpm/dynamic-relaxation::map-stiffness ((sim cl-mpm/dynamic-relaxation::mpm-sim-octree-quasi-static))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps))
      sim
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
                               (* 0.25d0 mp-factor))))))))))))))

(defmethod cl-mpm/setup::%estimate-elastic-dt ((sim mpm-sim-octree))
  (* (call-next-method) (expt 2 (- (cl-mpm::sim-multigrid-refinement sim)))))

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

(in-package :cl-mpm/dynamic-relaxation)

