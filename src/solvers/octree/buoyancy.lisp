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
