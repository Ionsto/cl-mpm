(in-package :cl-mpm/dynamic-relaxation)

(defun combi-stats-mpi (sim)
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
                    (not (cl-mpm/mesh::node-agg node))
                    (cl-mpm/mpi::node-in-computational-domain sim node))
               (with-accessors ((active cl-mpm/mesh::node-active)
                                (f-ext cl-mpm/mesh::node-external-force)
                                (f-rct cl-mpm/mesh::node-reaction-force)
                                ;(res cl-mpm/mesh::node-residual)
                                (res cl-mpm/mesh::node-force)
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
                    (cl-mpm/fastmaths::mag-squared res)
                    ;; (cl-mpm/fastmaths::mag-squared (cl-mpm/fastmaths::fast-.+ f-ext f-rct))
                    (cl-mpm/fastmaths::mag-squared f-ext)
                    (cl-mpm/fastmaths:dot disp f-ext))))
               (list 0d0 0d0 0d0 0d0)))
         (lambda (&rest args)
           (if args
               (mapcar (lambda (x y) (declare (double-float x y)) (+ x y)) (first args) (second args))
               (list 0d0 0d0 0d0 0d0))
           ))
      (declare (double-float energy oobf-num oobf-denom power))
      (declare (double-float mass energy oobf-num oobf-denom power))
      (when sim-agg
        (cl-mpm/aggregate::iterate-over-dimensions-with-mutex
         (cl-mpm/mesh::mesh-nd mesh)
         (lambda (d mut)
           (let* ((res (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-force d))
                  (f-ext (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-external-force d))
                  (f-rct (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-reaction-force d))
                  (ma (cl-mpm/aggregate::sim-global-ma sim))
                  (mv (cl-mpm/aggregate::assemble-global-scalar sim #'cl-mpm/mesh::node-mass))
                  (vi (cl-mpm/aggregate::assemble-internal-vec sim #'cl-mpm/mesh::node-velocity d))
                  (disp (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-displacment d)))

             (let ((doobf-num (cl-mpm/fastmaths::mag-squared
                               ;; res
                               (cl-mpm/aggregate::aggregate-vec sim res d)
                               ))
                   (doobf-denom
                     (+
                      ;; (cl-mpm/fastmaths::mag-squared
                      ;;  f-rct)
                      (cl-mpm/fastmaths::mag-squared
                       (cl-mpm/aggregate::aggregate-vec sim f-ext d))))
                   (dpower (cl-mpm/fastmaths:dot
                            disp f-ext))
                   (denergy (* 0.5d0 (cl-mpm/fastmaths::dot vi (cl-mpm/aggregate::@-mass-matrix-vec sim vi d)))))
               (sb-thread:with-mutex (mut)
                 (incf oobf-num doobf-num)
                 (incf oobf-denom doobf-denom)
                 (incf power dpower)
                 (incf energy denergy)))))))
      (let ((oobf 0d0)
            (oobf-num (cl-mpm/mpi::mpi-sum oobf-num))
            (oobf-denom (cl-mpm/mpi::mpi-sum oobf-denom))
            (power (cl-mpm/mpi::mpi-sum power))
            (energy (cl-mpm/mpi::mpi-sum energy)))
        (if (> oobf-denom 0d0)
            (setf oobf (sqrt (/ oobf-num oobf-denom)))
            (setf oobf (if (> oobf-num 0d0) sb-ext:double-float-positive-infinity 0d0)))
        (values energy oobf power)))))

(defmethod cl-mpm::update-dynamic-stats ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-mpi))
  (with-accessors ((stats-energy cl-mpm::sim-stats-energy)
                   (stats-oobf cl-mpm::sim-stats-oobf)
                   (stats-power cl-mpm::sim-stats-power)
                   (stats-work cl-mpm::sim-stats-work))
      sim
    (multiple-value-bind (e o p) (combi-stats-mpi sim)
      (setf stats-energy e
            stats-oobf o
            stats-power p)
      (incf stats-work p))))

