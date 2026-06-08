(in-package :cl-mpm/buoyancy)
(format t "Loading octree-buoyancy link~%")
(defmethod populate-cells-volume ((sim cl-mpm/dynamic-relaxation::mpm-sim-octree) clip-function)
  (cl-mpm/dynamic-relaxation::iterate-over-meshes
   sim
   (lambda (mesh mesh-index)
     (cl-mpm::iterate-over-cells
      mesh
      (lambda (cell)
        (with-accessors ((neighbours cl-mpm/mesh::cell-neighbours)
                         (index cl-mpm/mesh::cell-index)
                         (nodes cl-mpm/mesh::cell-nodes)
                         (pruned cl-mpm/mesh::cell-pruned)
                         (boundary cl-mpm/mesh::cell-boundary)
                         (pos cl-mpm/mesh::cell-centroid)
                         (vt cl-mpm/mesh::cell-volume)
                         (active cl-mpm/mesh::cell-active))
            cell
          (setf boundary nil)
          (when (and
                 active
                 ;; (not (cl-mpm/dynamic-relaxation::cell-coupling cell))
                 (not (= (cl-mpm/dynamic-relaxation::cell-octree-refine cell) 2))
                 (not (= (cl-mpm/dynamic-relaxation::cell-octree-refine cell) 1))
                 ;; (not (= (cl-mpm/dynamic-relaxation::cell-octree-refine cell) -1))
                 )
            (flet ((check-cell (c)
                     (with-accessors ((pos cl-mpm/mesh::cell-centroid)
                                      (neighbours cl-mpm/mesh::cell-neighbours)
                                      (vt cl-mpm/mesh::cell-volume)
                                      (nns cl-mpm/mesh::cell-nodes)
                                      )
                         c
                       (when (and (funcall clip-function pos)
                                  (cl-mpm/mesh::cell-active c)
                                  ;; (not (cl-mpm/dynamic-relaxation::cell-coupling c))
                                  ;; (not (= (cl-mpm/dynamic-relaxation::cell-octree-refine c) 2))
                                  ;; (not (= (cl-mpm/dynamic-relaxation::cell-octree-refine c) 1))
                                  ;; (not (= (cl-mpm/dynamic-relaxation::cell-octree-refine c) -1))
                                  )
                         (let ((vest 0d0))
                           (loop for n across nns
                                 do
                                    (when (cl-mpm/mesh:node-active n)
                                      (incf vest
                                            (* 0.25d0 (/
                                                       (cl-mpm/mesh::node-volume n)
                                                       (cl-mpm/mesh::node-volume-true n))))))
                           (when (< vest 0.5d0)
                             (setf boundary t)
                             (loop for n across nodes
                                   do
                                      (when (cl-mpm/mesh:node-active n)
                                        (sb-thread:with-mutex ((cl-mpm/mesh:node-lock n))
                                          (setf (cl-mpm/mesh::node-boundary-node n) t))))
                             ))))))
              (check-cell cell)
              (loop for neighbour in neighbours
                    do (check-cell neighbour)))))))
     )))
(defmethod apply-force-cells-3d ((sim cl-mpm/dynamic-relaxation::mpm-sim-octree) func-stress func-div clip-func)
  "Update force on nodes, with virtual stress field from cells"
  (with-accessors ()
      sim
    (declare (function func-stress func-div))
    (cl-mpm/dynamic-relaxation::iterate-over-meshes
     sim
     (lambda (mesh mesh-index)
       (cl-mpm::iterate-over-cells
        mesh
        (lambda (cell)
          ;;Iterate over a cells nodes
          (with-accessors (;(pos cl-mpm/mesh::cell-centroid)
                           ;; (trial-pos cl-mpm/mesh::cell-trial-centroid)
                           ;; (pos cl-mpm/mesh::cell-trial-centroid)
                           (cell-active cl-mpm/mesh::cell-active)
                           (cell-partial cl-mpm/mesh::cell-partial)
                           (cell-buoyancy cl-mpm/mesh::cell-buoyancy)
                           (cell-pressure cl-mpm/mesh::cell-pressure)
                           (df cl-mpm/mesh::cell-deformation-gradient)
                           (disp cl-mpm/mesh::cell-displacement))
              cell
            (let* ((trial-pos (get-cell-position cell)))
              (cl-mpm::cell-iterate-over-neighbours
               mesh
               cell
               (lambda (p volume node svp grads)
                 (with-accessors ((node-pos cl-mpm/mesh::node-position)
                                  (node-lock  cl-mpm/mesh:node-lock)
                                  (node-active cl-mpm/mesh:node-active)
                                  (node-boundary cl-mpm/mesh::node-boundary-node)
                                  (node-volume cl-mpm/mesh::node-volume)
                                  (node-boundary-scalar cl-mpm/mesh::node-boundary-scalar))
                     node
                   (declare (double-float volume svp))
                   (let (;; (grads (cl-mpm::gradient-push-forwards grads df))
                         (volume (* volume (cl-mpm/fastmaths::det-3x3 df)))
                         )
                     (when (and node-active node-boundary (funcall clip-func node-pos))
                       (sb-thread:with-mutex (node-lock)
                         ;; (cl-mpm/fastmaths::fast-.+
                         ;;  (cl-mpm/fastmaths::fast-scale!
                         ;;   (cl-mpm/utils::vector-from-list
                         ;;    (list (cl-mpm/utils::gradients-dx grads)
                         ;;          (cl-mpm/utils::gradients-dy grads)
                         ;;          (cl-mpm/utils::gradients-dz grads)))
                         ;;   (* 1d0 volume))
                         ;;  (cl-mpm/mesh::node-boundary-vec node)
                         ;;  (cl-mpm/mesh::node-boundary-vec node))
                         ;; (incf node-boundary-scalar (* -1d0 volume svp (the double-float (calculate-val-cell cell #'melt-rate))))
                         ))))))))))))))
(defmethod apply-buoyancy :after ((sim cl-mpm/dynamic-relaxation::mpm-sim-octree) func-stress func-div clip-function datum)
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (let ((vol
            (cl-mpm/aggregate::assemble-global-scalar
             sim
             #'cl-mpm/mesh::node-volume)))
      ;; (cl-mpm:iterate-over-nodes
      ;;  mesh
      ;;  (lambda (node)
      ;;    (when t;(> (cl-mpm/mesh::node-volume node) 0d0)
      ;;      (setf (cl-mpm/mesh::node-boundary-scalar node)
      ;;            (max 0d0 (cl-mpm/fastmaths:mag (cl-mpm/mesh::node-boundary-vec node)))))))
      (let ((proj-vol
              (cl-mpm/aggregate::extend-vec-nobcs
               sim
               (cl-mpm/aggregate::aggregate-vec-nobcs
                sim
                vol 0) 0))
            (nd (cl-mpm/mesh:mesh-nd mesh)))
        ;; (cl-mpm/dynamic-relaxation::iterate-bottom-up-mesh
        ;;  sim
        ;;  (lambda (mesh mesh-index)
        ;;    (cl-mpm::iterate-over-cells
        ;;     mesh
        ;;     (lambda (cell)
        ;;       (when (and (= (cell-octree-refine cell) 1))
        ;;         (iterate-over-sub-nodes
        ;;          sim
        ;;          mesh-index
        ;;          cell
        ;;          (lambda (n)
        ;;            (cl-mpm::iterate-over-cell-shape-local
        ;;             mesh
        ;;             cell
        ;;             (cl-mpm/mesh::node-position n)
        ;;             (lambda (nw w grads)
        ;;               (declare (ignore grads))
        ;;               (dotimes (d nd)
        ;;                 (incf
        ;;                  (cl-mpm/utils:varef (cl-mpm/mesh::node-boundary-vec nw) d)
        ;;                  (*
        ;;                   w
        ;;                   (cl-mpm/utils:varef (cl-mpm/mesh::node-boundary-vec n) d)))))))))))))
        ;; (cl-mpm::iterate-over-nodes
        ;;  mesh
        ;;  (lambda (n)
        ;;    (when (> (cl-mpm/mesh::node-volume n) 0d0)
        ;;      (cl-mpm/fastmaths:fast-scale! (cl-mpm/mesh::node-boundary-vec n) (* 1d0 (cl-mpm/mesh::node-volume n))))))
        ;; (loop for d from 0 below (cl-mpm/mesh:mesh-nd mesh)
        ;;       do
        ;;          (progn
        ;;            (cl-mpm/dynamic-relaxation::aggregate-vector-up-grid
        ;;             sim
        ;;             #'cl-mpm/mesh::node-boundary-vec
        ;;             d)
        ;;            (cl-mpm/dynamic-relaxation::project-vector-down-grid
        ;;             sim
        ;;             #'cl-mpm/mesh::node-boundary-vec
        ;;             d)
        ;;            ))
        ;; (cl-mpm::iterate-over-nodes
        ;;  mesh
        ;;  (lambda (n)
        ;;    (when (> (cl-mpm/mesh::node-volume n) 0d0)
        ;;      (cl-mpm/fastmaths:fast-scale! (cl-mpm/mesh::node-boundary-vec n) (/ 1d0 (cl-mpm/mesh::node-volume n))))))
        ;; (cl-mpm/aggregate::iterate-over-dimensions
        ;;  (cl-mpm/mesh::mesh-nd mesh)
        ;;  (lambda (d)
        ;;    (cl-mpm/aggregate::project-global-vec
        ;;     sim
        ;;     (cl-mpm/aggregate::extend-vec-nobcs
        ;;      sim
        ;;      (cl-mpm/aggregate::aggregate-vec-nobcs
        ;;       sim
        ;;       (cl-mpm/aggregate::assemble-global-vec
        ;;        sim
        ;;        #'cl-mpm/mesh::node-boundary-vec
        ;;        d)
        ;;       d)
        ;;      d)
        ;;     #'cl-mpm/mesh::node-boundary-vec
        ;;     d)))
        ;; (cl-mpm/aggregate::iterate-over-dimensions
        ;;  (cl-mpm/mesh::mesh-nd mesh)
        ;;  (lambda (d)
        ;;    (cl-mpm/aggregate::project-global-vec
        ;;     sim
        ;;     (cl-mpm/fastmaths:fast-./
        ;;      (cl-mpm/aggregate::extend-vec
        ;;       sim
        ;;       (cl-mpm/aggregate::aggregate-vec
        ;;        sim
        ;;        (cl-mpm/fastmaths::fast-.*
        ;;         vol
        ;;         (cl-mpm/aggregate::assemble-global-vec
        ;;          sim
        ;;          #'cl-mpm/mesh::node-boundary-vec
        ;;          d))
        ;;        d)
        ;;       d)
        ;;      proj-vol
        ;;      )
        ;;     #'cl-mpm/mesh::node-boundary-vec
        ;;     d)))
        ;; (format t "Reprojection ~%")
        ;; (pprint vol)
        ;; (pprint proj-vol)
        ;; (let ((sim *sim*))
        ;;   (cl-mpm/aggregate::project-global-scalar
        ;;    sim
        ;;    (cl-mpm/aggregate::extend-vec-nobcs
        ;;     sim
        ;;     (cl-mpm/aggregate::aggregate-vec-nobcs
        ;;      sim
        ;;      (cl-mpm/aggregate::assemble-global-scalar
        ;;       sim
        ;;       #'cl-mpm/mesh::node-volume) 0) 0)
        ;;    cl-mpm/mesh::node-volume))
        ;; (cl-mpm/aggregate::project-global-scalar
        ;;  sim
        ;;   (cl-mpm/aggregate::extend-vec-nobcs
        ;;    sim
        ;;    (cl-mpm/aggregate::aggregate-vec-nobcs
        ;;     sim
        ;;     (cl-mpm/aggregate::assemble-global-scalar
        ;;      sim
        ;;      #'cl-mpm/mesh::node-boundary-scalar) 0) 0)
        ;;  cl-mpm/mesh::node-boundary-scalar)
        ;; (cl-mpm/aggregate::project-global-scalar
        ;;  sim
        ;;  (cl-mpm/fastmaths::fast-./
        ;;   (cl-mpm/aggregate::extend-vec-nobcs
        ;;    sim
        ;;    (cl-mpm/aggregate::aggregate-vec-nobcs
        ;;     sim
        ;;     (cl-mpm/fastmaths::fast-.*
        ;;      vol
        ;;      (cl-mpm/aggregate::assemble-global-scalar
        ;;       sim
        ;;       #'cl-mpm/mesh::node-boundary-scalar)) 0) 0)
        ;;   proj-vol)
        ;;  cl-mpm/mesh::node-boundary-scalar)
        )
      (cl-mpm:iterate-over-nodes
       mesh
       (lambda (node)
         (when t
                                        ;(> (cl-mpm/mesh::node-volume node) 0d0)
           (setf (cl-mpm/mesh::node-boundary-scalar node)
                 ;; (max 0d0 (abs (cl-mpm/mesh::node-boundary-scalar node)))
                 (max 0d0 (cl-mpm/fastmaths:mag (cl-mpm/mesh::node-boundary-vec node)))
                 ;; (* (cl-mpm/mesh::node-volume node)
                 ;;    (max 0d0 (cl-mpm/fastmaths:mag (cl-mpm/mesh::node-boundary-vec node)))
                 ;;    (max 0d0 (abs (cl-mpm/mesh::node-boundary-scalar node))))
                 )
           )))))
  )

(defmethod locate-mps-cells ((sim cl-mpm/dynamic-relaxation::mpm-sim-octree) clip-function)
  "mark boundary nodes based on neighbour mp inclusion"
  (with-accessors ((refinement cl-mpm::sim-multigrid-refinement)
                   (mesh-list cl-mpm::sim-mesh-list)
                   )
      sim
    (cl-mpm/dynamic-relaxation::iterate-over-meshes
      sim
      (lambda (mesh mesh-index)
        (cl-mpm::iterate-over-cells
         mesh
         (lambda (cell)
           (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                            (neighbours cl-mpm/mesh::cell-neighbours)
                            (index cl-mpm/mesh::cell-index)
                            ;; (nodes cl-mpm/mesh::cell-nodes)
                            (pruned cl-mpm/mesh::cell-pruned)
                            (active cl-mpm/mesh::cell-active)
                            (partial cl-mpm/mesh::cell-partial)
                            (boundary cl-mpm/mesh::cell-boundary)
                            (pos cl-mpm/mesh::cell-centroid)
                            (nodes cl-mpm/mesh::cell-nodes))
               cell
             ;; (let ((vest 0d0))
             ;;   (loop for n across nodes
             ;;         do (when t;(cl-mpm/mesh:node-active n)
             ;;              (incf vest
             ;;                    (* 0.25d0 (/
             ;;                               (cl-mpm/mesh::node-volume n)
             ;;                               (cl-mpm/mesh::node-volume-true n))))))
             ;;   (when (< vest 0.9d0)
             ;;     (setf boundary t)
             ;;     (loop for n across nodes
             ;;           do
             ;;              (when (cl-mpm/mesh:node-active n)
             ;;                (sb-thread:with-mutex ((cl-mpm/mesh:node-lock n))
             ;;                  (setf (cl-mpm/mesh::node-boundary-node n) t))))))

             (let ((any-boundary nil))
               (loop for n across nodes
                     do (when (and (cl-mpm/mesh:node-active n)
                                   (funcall clip-function (cl-mpm/fastmaths::fast-.+
                                                           (cl-mpm/mesh::node-position n)
                                                           (cl-mpm/mesh::node-displacment n))))
                          (setf any-boundary t)
                          (sb-thread:with-mutex ((cl-mpm/mesh:node-lock n))
                            (setf (cl-mpm/mesh::node-boundary-node n) t))))
               (setf boundary any-boundary))
             ;; (when (and
             ;;        active
             ;;        partial
             ;;        ;; (funcall clip-function trial-pos)
             ;;        )
             ;;   ;; (set-boundary cell)
             ;;   (loop for neighbour in neighbours
             ;;         do
             ;;            (when (cl-mpm/mesh::cell-active neighbour)
             ;;              (let ((any-boundary nil))
             ;;                (with-accessors ((nodes cl-mpm/mesh::cell-nodes)
             ;;                                 (boundary cl-mpm/mesh::cell-boundary))
             ;;                    neighbour
             ;;                  (loop for n across nodes
             ;;                        do (when (and (cl-mpm/mesh:node-active n)
             ;;                                      (funcall clip-function (cl-mpm/fastmaths::fast-.+
             ;;                                                              (cl-mpm/mesh::node-position n)
             ;;                                                              (cl-mpm/mesh::node-displacment n))))
             ;;                             (setf any-boundary t)
             ;;                             (sb-thread:with-mutex ((cl-mpm/mesh:node-lock n))
             ;;                               (setf (cl-mpm/mesh::node-boundary-node n) t)))))
             ;;                (setf (cl-mpm/mesh::cell-boundary neighbour) any-boundary)
             ;;                ;; (set-boundary neighbour)
             ;;                ))))
             )))))
    ;; (loop for i from 0 to refinement
    ;;       do (let ((mesh (aref mesh-list i)))
    ;;            (cl-mpm::iterate-over-cells
    ;;             mesh
    ;;             (lambda (cell)
    ;;               (with-accessors ((active cl-mpm/mesh::cell-active)
    ;;                                (coupling cell-coupling)
    ;;                                (boundary cl-mpm/mesh::cell-boundary))
    ;;                   cell
    ;;                 (when (and boundary)
    ;;                   (when (< i refinement)
    ;;                     (cl-mpm/dynamic-relaxation::iterate-over-sub-cells
    ;;                      sim
    ;;                      i
    ;;                      cell
    ;;                      (lambda (c)
    ;;                        (when boundary
    ;;                          (setf (cl-mpm/mesh::cell-boundary c) t)))))

    ;;                   (loop for n across (cl-mpm/mesh::cell-nodes cell)
    ;;                         do (when (and (cl-mpm/mesh:node-active n)
    ;;                                       (funcall clip-function (cl-mpm/fastmaths::fast-.+
    ;;                                                               (cl-mpm/mesh::node-position n)
    ;;                                                               (cl-mpm/mesh::node-displacment n))))
    ;;                              (setf (cl-mpm/mesh::node-boundary-node n) t)))))))))
    ))
