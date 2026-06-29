(in-package :cl-mpm/mpi)


(defmethod cl-mpm/aggregate::locate-aggregate-nodes ((sim cl-mpm/mpi::mpm-sim-mpi))
  ;; (pprint "hello mpi aggregation")
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (cl-mpm:iterate-over-nodes
     mesh
     (lambda (node)
       (setf (cl-mpm/mesh::node-agg-interior-cell node) nil)
       (setf (cl-mpm/mesh::node-interior node) nil)
       (setf (cl-mpm/mesh::node-agg node) nil)))

    (cl-mpm::iterate-over-cells
     mesh
     (lambda (node)
       (setf (cl-mpm/mesh::cell-agg node) nil
             (cl-mpm/mesh::cell-aggregate-element node) nil
             (cl-mpm/mesh::cell-interior node) nil)))

    ;;First set all outside nodes as aggregate
    (cl-mpm::iterate-over-cells
     mesh
     (lambda (cell)
       (with-accessors ((partial cl-mpm/mesh::cell-partial)
                        (nodes cl-mpm/mesh::cell-nodes)
                        (agg cl-mpm/mesh::cell-agg)
                        (neighbours cl-mpm/mesh::cell-cartesian-neighbours)
                        (active cl-mpm/mesh::cell-active))
           cell
         (when (and active partial
                    ;; (in-computational-domain sim (cl-mpm/mesh::cell-centroid cell))
                    )
           (setf agg t)
           ;;Set all our neighbours to also be aggregate
           (loop for n in neighbours
                 do
                    (when (in-computational-domain sim (cl-mpm/mesh::cell-centroid n))
                      (setf (cl-mpm/mesh::cell-agg n) t)))))))
    ;;Next set any nodes on agg elements to be agg
    (cl-mpm::iterate-over-cells
     mesh
     (lambda (cell)
       (with-accessors ((nodes cl-mpm/mesh::cell-nodes)
                        (agg cl-mpm/mesh::cell-agg)
                        (active cl-mpm/mesh::cell-active))
           cell
         ;;Set all aggregated nodes
         (when (and active agg
                    ;; (in-computational-domain sim (cl-mpm/mesh::node-position cell))
                    )
           (loop for n across nodes
                 do (when (and (cl-mpm/mesh::node-active n)
                               (in-computational-domain sim (cl-mpm/mesh::node-position n)))
                      (setf (cl-mpm/mesh::node-agg n) t)))))))

    ;;For each aggregate node, locate the closest support cell
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (node)
       (with-accessors ((active cl-mpm/mesh::node-active)
                        (agg cl-mpm/mesh::node-agg))
           node
         (when (and active agg
                    (in-computational-domain sim (cl-mpm/mesh::node-position node))
                    )
           (let ((closest-elem
                   (cl-mpm/aggregate::get-closest-cell
                    mesh
                    (cl-mpm/mesh::node-position node)
                    :exclude
                    (lambda (c)
                      (not
                       (cl-mpm/mpi::in-computational-domain-buffer
                        sim
                        (cl-mpm/mesh::cell-centroid c)
                        (- (cl-mpm/mesh::mesh-resolution (cl-mpm::sim-mesh sim)))))))))
             (if closest-elem
                 (progn
                   (setf (cl-mpm/mesh::node-agg-interior-cell node) closest-elem)
                   (setf (cl-mpm/mesh::cell-interior closest-elem) t))
                 (format t "No closest elem?~%")))))))

    (cl-mpm::iterate-over-cells
     mesh
     (lambda (cell)
       (with-accessors ((nodes cl-mpm/mesh::cell-nodes)
                        (agg cl-mpm/mesh::cell-agg)
                        (int cl-mpm/mesh::cell-interior)
                        (active cl-mpm/mesh::cell-active))
           cell
         ;;Set all aggregated nodes - also set them all to be interior
         (when (and active int)
           (loop for n across nodes
                 do (when (cl-mpm/mesh::node-active n)
                      (setf (cl-mpm/mesh::node-agg n) t
                            (cl-mpm/mesh::node-interior n) t
                            (cl-mpm/mesh::node-agg-interior-cell n) nil)))))))
    (setf
     (cl-mpm/aggregate::sim-agg-nodes-fdc sim)
     (cl-mpm::filter-nodes sim (lambda (n) (and (cl-mpm/mesh::node-active n)
                                                (cl-mpm/mesh::node-interior n))))
     (cl-mpm/aggregate::sim-agg-nodes-fd sim)
     (cl-mpm::filter-nodes sim (lambda (n) (and (cl-mpm/mesh::node-active n)
                                                (cl-mpm/mesh::node-agg n)))))

    (let ((nodes (cl-mpm/aggregate::sim-agg-nodes-fdc sim)))
      (cl-mpm/utils::bpdotimes (fdc (length nodes))
        (let ((n (aref nodes fdc)))
             (setf (cl-mpm/mesh::node-agg-fdc n) fdc))))

    (let ((nodes (cl-mpm/aggregate::sim-agg-nodes-fd sim)))
      (cl-mpm/utils::bpdotimes (fd (length nodes))
        (let ((n (aref nodes fd)))
             (setf (cl-mpm/mesh::node-agg-fd n) fd))))))
