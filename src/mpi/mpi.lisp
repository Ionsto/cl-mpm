(in-package :cl-mpm/mpi)
(declaim #.cl-mpm/settings:*optimise-setting*)

(defmethod sim-add-mp ((sim mpm-sim-mpi) mp)
  (with-accessors ((uid-counter cl-mpm::sim-unique-index-counter)
                   (lock cl-mpm::sim-unique-index-lock))
      sim
    (let* ((size (cl-mpi:mpi-comm-size))
           (rank (cl-mpi:mpi-comm-rank))
           (shift (ceiling (log size 2))))
      (sb-thread:with-mutex (lock)
        (setf (cl-mpm/particle::mp-unique-index mp) (+ (ash uid-counter shift) rank))
        (incf uid-counter)))
    (call-next-method)))


(defgeneric update-min-domain-size (sim))
(defmethod update-min-domain-size ((sim mpm-sim-mpi))
  (setf (mpm-sim-mpi-min-size sim)
        (* (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) (mpm-sim-mpi-halo-depth sim))))

(defmethod update-min-domain-size ((sim mpm-sim-mpi-damage))
  (setf (mpm-sim-mpi-min-size sim)
        (max
         (* (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) (mpm-sim-mpi-halo-depth sim))
         (mpm-sim-mpi-halo-damage-size sim))))

(defmethod (setf cl-mpm::sim-mesh) :after (value (sim mpm-sim-mpi))
  (update-min-domain-size sim))
(defmethod (setf mpm-sim-mpi-halo-depth) :after (value (sim mpm-sim-mpi))
  (update-min-domain-size sim))
(defmethod (setf mpm-sim-mpi-halo-damage-size) :after (value (sim mpm-sim-mpi-damage))
  (update-min-domain-size sim))







(defmethod cl-mpm::calculate-min-dt-mps ((sim cl-mpm/mpi::mpm-sim-mpi))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mass-scale cl-mpm::sim-mass-scale))
      sim
    (let ((inner-factor
            ;;MPI not a fan of most-positive-double-float
            1d50
            ;most-positive-double-float
                        ))
      (iterate-over-nodes-serial
       mesh
       (lambda (node)
         (with-accessors ((node-active  cl-mpm/mesh:node-active)
                          (node-pos cl-mpm/mesh::node-position)
                          (pmod cl-mpm/mesh::node-pwave)
                          (mass cl-mpm/mesh::node-mass)
                          (svp-sum cl-mpm/mesh::node-svp-sum)
                          (vol cl-mpm/mesh::node-volume)
                          ) node
           (when (and node-active
                      ;(in-computational-domain sim node-pos)
                      (> vol 0d0)
                      (> pmod 0d0)
                      (> svp-sum 0d0))
             (let ((nf (/ mass (* vol (/ pmod svp-sum)))))
                 (when (< nf inner-factor)
                   (setf inner-factor nf)))))))
      (let ((rank (cl-mpi:mpi-comm-rank))
            (size (cl-mpi:mpi-comm-size)))
        (static-vectors:with-static-vector (source 1 :element-type 'double-float :initial-element inner-factor)
          (static-vectors:with-static-vector (dest 1 :element-type 'double-float :initial-element 0d0)
            (mpi-allreduce source dest cl-mpi:+mpi-min+ :type cl-mpi:+mpi-double+)
            (setf inner-factor (aref dest 0))
            (if (< inner-factor most-positive-double-float)
                (progn
                  ;; (format t "Rank ~D: dt - ~F~%" rank (* (sqrt mass-scale) (sqrt inner-factor) (cl-mpm/mesh:mesh-resolution mesh)))
                  ;; (format t "global : dt - ~F~%" (* (sqrt mass-scale) (sqrt inner-factor) (cl-mpm/mesh:mesh-resolution mesh)))
                  (* (sqrt mass-scale) (sqrt inner-factor) (cl-mpm/mesh:mesh-resolution mesh)))
                (cl-mpm::sim-dt sim))))))))

