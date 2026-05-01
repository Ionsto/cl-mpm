(in-package :cl-mpm/mpi)

(defun exchange-domain-bounds (sim)
  (with-accessors ((mesh sim-mesh)
                   (domain-bounds mpm-sim-mpi-domain-bounds ))
      sim
    (with-accessors ((h mesh-resolution))
        mesh
      (setf domain-bounds
            (mapcar (lambda (v) (mapcar (lambda (x) (* (round x h) h)) v))
                    domain-bounds)))))

(defun in-computational-domain (sim pos)
  (with-accessors ((bounds mpm-sim-mpi-domain-bounds))
      sim
      (let ((in-bounds t)
            (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim))))
        (loop for i from 0 below nd
              for bound in bounds
              while in-bounds
              do
                 (destructuring-bind (bl bu) bound
                   (declare (double-float bl bu))
                   (when (not (= bu bl))
                     (setf in-bounds
                           (and
                            in-bounds
                            (and
                             (> bu (cl-mpm/utils:varef pos i))
                             (<= bl (cl-mpm/utils:varef pos i)))
                            )))))
        in-bounds
        )))
(defun node-in-computational-domain (sim node)
  (with-accessors ((pos cl-mpm/mesh::node-position))
      node
    (with-accessors ((bounds mpm-sim-mpi-domain-bounds))
        sim
      (let ((in-bounds t)
            (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim))))
        (loop for i from 0 below nd
              for bound in bounds
              while in-bounds
              do
                 (destructuring-bind (bl bu) bound
                   (declare (double-float bl bu))
                   (when (not (= bu bl))
                     (setf in-bounds
                           (and
                            in-bounds
                            (and
                             (> bu (cl-mpm/utils:varef pos i))
                             (<= bl (cl-mpm/utils:varef pos i)))
                            )))))
        in-bounds
        ))))

(defun in-computational-domain-buffer (sim pos node-buffer)
  (let ((in-bounds t)
        (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim))))
    (loop for i from 0 below 3;(cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)) 
          do
             (destructuring-bind (bl bu) (nth i (mpm-sim-mpi-domain-bounds sim))
               (when (not (= bu bl))
                 (setf in-bounds
                       (and
                        in-bounds
                        (and
                         (> (+ bu (* node-buffer h)) (cl-mpm/utils:varef pos i))
                         (<= (- bl (* node-buffer h)) (cl-mpm/utils:varef pos i)))
                        )))))
    in-bounds
    ))

(defun calculate-domain-sizes (sim &optional size)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (divs mpm-sim-mpi-domain-count))
      sim
    (when size
      (setf divs size))
    (let ((mc (cl-mpm/mesh:mesh-mesh-size mesh)))
      ;; (loop for i from 0 below (- rank-count 1)
      ;;       do (incf (nth (mod i 3) divs)))
      ;; divs
      (mapcar #'/ mc divs))))

(defun setup-domain-bounds (sim &key (domain-scaler (lambda (x) x)))

  (with-accessors ((mesh-size cl-mpm/mesh::mesh-mesh-size)
                   (mesh-count cl-mpm/mesh::mesh-count)
                   (h cl-mpm/mesh::mesh-resolution))
      (cl-mpm:sim-mesh sim)
    (let* ((rank (cl-mpi:mpi-comm-rank))
           (count (cl-mpi:mpi-comm-size))
           (index (mpi-rank-to-index sim rank))
           (comp-size (mpm-sim-mpi-domain-count sim)))
      (setf (mpm-sim-mpi-domain-bounds sim)
            (mapcar (lambda (domain size) (mapcar (lambda (x) (* x size)) domain))
                    (funcall domain-scaler
                             (loop for i from 0 to 2
                                   collect
                                   (if (> (nth i comp-size) 1)
                                       (list (/ (* (/ (nth i mesh-size) (nth i comp-size)) (nth i index)) (nth i mesh-size))
                                             (/ (* (/ (nth i mesh-size) (nth i comp-size)) (+ (nth i index) 1)) (nth i mesh-size)))
                                       (list 0d0 1d0))))
                    mesh-size))
      (format t "Rank ~D: Domain bounds ~A~%" rank (mpm-sim-mpi-domain-bounds sim))
      (exchange-domain-bounds sim)
      (format t "Rank ~D: Exchanged domain bounds ~A~%" rank (mpm-sim-mpi-domain-bounds sim))))
  )

(defun origin-domain-scaler (sim real-x real-y real-z)
  "Build a domain scalar that re-scales all mpi domains to be tight packed from (0 - real-x) (0 - real-y) (0 - real-z)
leaves a hanging mpi domain at the back"
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (with-accessors ((size cl-mpm/mesh:mesh-mesh-size))
        mesh
      (let ((target-size (list real-x real-y real-z)))
        (lambda (domain)
          (loop for i from 0
                for dim in domain
                collect
                (mapcar (lambda (x)
                          (if (and (> x 0d0)
                                   (< x 1d0)
                                   (nth i target-size)
                                   )
                              (* x (/ (nth i target-size) (nth i size)))
                              x))
                        dim)))))))


(defun domain-decompose (sim &key (domain-scaler (lambda (x) x)))
  (%domain-decompose sim domain-scaler)
  )
(defgeneric %domain-decompose (sim domain-scaler)
  )
(defmethod %domain-decompose (sim domain-scaler)
  "The aim of domain decomposition is to take a full simulation and cut it into subsections for MPI"

  (with-accessors ((mesh-size cl-mpm/mesh::mesh-mesh-size)
                   (mesh-count cl-mpm/mesh::mesh-count)
                   (h cl-mpm/mesh::mesh-resolution))
      (cl-mpm:sim-mesh sim)
    (let* ((rank (cl-mpi:mpi-comm-rank))
           (count (cl-mpi:mpi-comm-size))
           (index (mpi-rank-to-index sim rank))
           (x-length (first mesh-size))
           (slice-size (/ x-length count))
           (slice-count count)
           (bound-lower (* rank slice-size))
           (bound-upper (* (+ 1 rank) slice-size))
           (comp-size (mpm-sim-mpi-domain-count sim)))
      (unless (mpm-sim-mpi-domain-bounds sim)
        (setup-domain-bounds sim :domain-scaler domain-scaler))
      (setf (mpm-sim-mpi-domain-index sim) index)
      (set-mp-mpi-index sim)
      (clear-ghost-mps sim)
      ;; (format t "Sim MPs: ~a~%" (length (cl-mpm:sim-mps sim)))
      (let ((min-mps (mpi-min (length (cl-mpm:sim-mps sim))))
            (max-mps (mpi-max (length (cl-mpm:sim-mps sim)))))
        (when (and (= rank 0)
                   (> min-mps 0))
          (format t "Occupancy ratio : ~F%~%" (* 100d0 (/ max-mps min-mps))))))
    ))
(defparameter *prune-nodes* t)
(defmethod %domain-decompose :after ((sim cl-mpm/mpi::mpm-sim-mpi-nodes) domain-scaler)
  (let* ((rank (cl-mpi:mpi-comm-rank))
         (size (cl-mpi:mpi-comm-size)))
    (with-accessors ((mesh cl-mpm:sim-mesh))
        sim
      (let* ((nd-nodes (cl-mpm/mesh:mesh-nodes mesh))
             (all-nodes (make-array (array-total-size nd-nodes) :displaced-to nd-nodes :displaced-index-offset 0))
             (index (mpi-rank-to-index sim rank))
             (bounds-list (mpm-sim-mpi-domain-bounds sim))
             (h (cl-mpm/mesh:mesh-resolution mesh))
             (halo-depth (cl-mpm/mpi::mpm-sim-mpi-halo-depth sim)))
        (loop for i from 0 to 2
              do
                 (let ((id-delta (list 0 0 0)))
                   (setf (nth i id-delta) 1)
                   (let ((left-neighbor (mpi-index-to-rank sim (mapcar #'- index id-delta)))
                         (right-neighbor (mpi-index-to-rank sim (mapcar #'+ index id-delta))))
                     (destructuring-bind (bl bu) (nth i bounds-list)
                       (labels
                           ((halo-filter (test)
                              (let ((res
                                      (lparallel:premove-if-not
                                       (lambda (node)
                                         (if node
                                             (funcall test (nth i (cl-mpm/mesh:node-index node)))
                                             nil))
                                       all-nodes)))
                                (make-array (length res) :initial-contents res)))
                            (left-filter ()
                              (halo-filter (lambda (pos)
                                             (and
                                              (>= pos (- (/ bl h) halo-depth))
                                              (<= pos (+ (/ bl h) halo-depth))
                                              )
                                             )))
                            (right-filter ()
                              (halo-filter (lambda (pos)
                                             (and
                                              (>= pos (- (/ bu h) halo-depth))
                                              (<= pos (+ (/ bu h) halo-depth)))
                                             ))))
                         (when t;(not (= left-neighbor -1))
                           (setf (nth 0 (nth i (mpm-sim-mpi-halo-node-list sim)))
                                 (left-filter)))

                         (when t;(not (= right-neighbor -1))
                           (setf (nth 1 (nth i (mpm-sim-mpi-halo-node-list sim)))
                                 (right-filter))))))))

        (with-accessors ((nodes cl-mpm/mesh:mesh-nodes)
                         (cells cl-mpm/mesh::mesh-cells))
            mesh

          (let* ((domain-sizes (mapcar (lambda (x) (abs (reduce #'- x)))
                                       (cl-mpm/mpi::mpm-sim-mpi-domain-bounds sim)))
                 (h (cl-mpm/mesh:mesh-resolution mesh))
                 (min-domain-length (reduce #'max (remove 0d0 domain-sizes)))
                 (buffer-size (+ 1 (max halo-depth (ceiling min-domain-length h)))))
            (format t "Compacting active nodes ~D nodes away~%" buffer-size)
            (format t "Domain size ~A~%" domain-sizes)
            (setf (cl-mpm/mesh::mesh-active-nodes mesh)
                  (cl-mpm::filter-nodes
                   sim
                   (lambda (node)
                     (when node
                       (not (in-computational-domain-buffer
                                    sim
                                    (cl-mpm/mesh::node-position node)
                                    buffer-size))))))
            (setf (cl-mpm/mesh::mesh-active-cells mesh)
                  (cl-mpm::filter-array-cells
                   sim
                   (lambda (cell)
                     (when cell
                       (not (in-computational-domain-buffer
                             sim
                             (cl-mpm/mesh::cell-centroid cell)
                             buffer-size))))))
            (setf (cl-mpm::sim-active-bcs sim)
                  (cl-mpm::filter-bcs
                   sim
                   (lambda (cell)
                     (when cell
                       (not (in-computational-domain-buffer
                             sim
                             (cl-mpm/mesh::index-to-position mesh (cl-mpm/bc::bc-index cell))
                             buffer-size))))))
            )
          )
        (when *prune-nodes*
          (with-accessors ((nodes cl-mpm/mesh:mesh-nodes)
                           (cells cl-mpm/mesh::mesh-cells))
              mesh
            (with-accessors ((bcs cl-mpm:sim-bcs))
                sim
              (let* ((domain-sizes (mapcar (lambda (x) (abs (reduce #'- x)))
                                           (cl-mpm/mpi::mpm-sim-mpi-domain-bounds sim)))
                     (h (cl-mpm/mesh:mesh-resolution mesh))
                     (min-domain-length (reduce #'max (remove 0d0 domain-sizes)))
                     (buffer-size (+ 1 (max halo-depth (ceiling min-domain-length h)))))
                (format t "Pruning nodes and BCs up to ~D nodes away~%" buffer-size)
                (format t "Domain size ~A~%" domain-sizes)
                ;;Remove nil bcs

                (setf bcs (delete nil bcs))
                (let ((prune-count 0))
                  (dotimes (i (length bcs))
                    (let ((bc (aref bcs i)))
                      (when bc
                        (let ((index (cl-mpm/bc:bc-index bc)))
                          ;;nil indexes indiciate global bcs
                          (when index
                            (progn
                              ;;SBCL is very sure that the result cannot be nil!
                              ;; (when (not (equal (cl-mpm/mesh:get-node mesh index) nil))
                              ;;   (let ((node (cl-mpm/mesh:get-node mesh index)))
                                  ;; (format t "Pruning BC at index: ~A node: ~A~%" index node)
                                  (when (not (in-computational-domain-buffer
                                              sim
                                              ;; (cl-mpm/mesh::node-position node)
                                              (cl-mpm/utils:vector-from-list  (cl-mpm/mesh:index-to-position mesh index))
                                              buffer-size))
                                    (incf prune-count)
                                    (setf (aref bcs i) nil))))))))

                  (format t "Rank ~D - Pruned ~D bcs~%" rank prune-count))
                (setf (cl-mpm:sim-bcs sim) (delete nil bcs))

                ;;Trim out all nodes that we can get rid of
                (let ((prune-count 0))
                  (dotimes (i (array-total-size nodes))
                    (let ((node (row-major-aref nodes i)))
                      (when node
                        (when (not (in-computational-domain-buffer
                                    sim
                                    (cl-mpm/mesh::node-position node)
                                    buffer-size))
                          (incf prune-count)
                          (setf (row-major-aref nodes i) nil)))))
                  (format t "Rank ~D - Pruned ~D nodes~%" rank prune-count))

                ;;Prune orphan bcs
                (let ((prune-count 0))
                  (dotimes (i (length bcs))
                    (let ((bc (aref bcs i)))
                      (when bc
                        (let ((index (cl-mpm/bc:bc-index bc)))
                          ;;nil indexes indiciate global bcs
                          (when index
                            ;;SBCL is very sure that the result cannot be nil!
                            (when (equal (cl-mpm/mesh:get-node mesh index) nil)
                              (incf prune-count)
                              (setf (aref bcs i) nil)))))))
                  (setf (cl-mpm:sim-bcs sim) (delete nil bcs))
                  (format t "Rank ~D - Pruned ~D orphan bcs~%" rank prune-count))

                (loop for bc across (cl-mpm:sim-bcs sim)
                      do
                        (when bc (when (equal (cl-mpm/mesh:get-node mesh (cl-mpm/bc:bc-index bc)) nil)
                                   (error "How on earth has bc ~A got a nil node rank ~D" bc rank))))
                ;;Cells
                (let ((prune-count 0))
                  (dotimes (i (array-total-size cells))
                    (let ((cell (row-major-aref cells i)))
                      (when cell
                        (when (not (in-computational-domain-buffer sim (cl-mpm/mesh::cell-centroid cell)
                                                                   0
                                                                   ;; buffer-size
                                                                   ))
                          (incf prune-count)
                          (setf (row-major-aref cells i) nil)))))
                  (format t "Rank ~D - Pruned ~D cells~%" rank prune-count))))))

        ))))



(defgeneric load-balance-metric (sim))

(defmethod load-balance-metric (sim)
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (set-mp-mpi-index sim)
    (lparallel:pcount-if (lambda (mp) (= (cl-mpm/particle::mp-mpi-index mp) rank))
                         (cl-mpm:sim-mps sim))))


(defun load-balance-setup (sim)
  (setf (mpm-sim-mpi-load-metric sim) (float (load-balance-metric sim) 0d0)))

(defun load-balance-dimension (sim dim &key (step-size 1d-2)
                                         (exchange nil))
  (let* ((rank (cl-mpi::mpi-comm-rank))
         (index-bottom-rank (mpi-rank-to-index sim rank))
         (dim-length (nth dim (mpm-sim-mpi-domain-count sim)))
         (dim-index (nth dim (mpi-rank-to-index sim rank)))
         (dim-size (abs (- (apply #'- (nth dim (mpm-sim-mpi-domain-bounds sim))))))
         (metric-array (make-array dim-length :element-type '(signed-byte 32) :initial-element 0))
         (size-array (make-array dim-length :element-type 'double-float :initial-element 0d0))
         (increment-length (- dim-length 1))
         (increment-array (make-array increment-length :element-type 'double-float :initial-element 0d0)))
    ;; (format t "Dims ~D length ~D~%" dim dim-length)
    (setf (nth dim index-bottom-rank) 0)
    (let ((bottom-rank (every #'zerop index-bottom-rank))
          (min-size (* 2 (mpm-sim-mpi-min-size sim))))
      ;; (when (= rank 0)
      ;;   (format t "Min size ~F~%" min-size)
      ;;   )
      (when (> dim-length 1)
        (setf (aref metric-array dim-index)
              (coerce (round (mpm-sim-mpi-load-metric sim)) 'integer)
              )
        ;; (format t "Metric ~E~%" (float (mpm-sim-mpi-load-metric sim) 0d0))
        (mpi-vector-integer-sum metric-array)
        ;; (format t "Dim index ~D~%" dim-index)
        ;; (when (= rank 0)
        ;;   (format t "Metric ~A~%" metric-array))
        (setf (aref size-array dim-index)
              (float dim-size 0d0))
        (mpi-vector-max size-array)
        (when bottom-rank
          (when (< dim-index increment-length)
            (setf (aref increment-array dim-index)
                  (float
                   (* step-size
                      (min (aref size-array dim-index)
                           (aref size-array (+ dim-index 1)))
                      (float
                       (if (or (> (aref metric-array dim-index) 0d0)
                               (> (aref metric-array (1+ dim-index)) 0d0))
                           (/ (- (aref metric-array dim-index)
                                 (aref metric-array (+ dim-index 1)))
                              (max (aref metric-array dim-index)
                                   (aref metric-array (+ dim-index 1))))
                           (signum
                            (-
                             (- (aref metric-array dim-index)
                                (aref metric-array (+ dim-index 1))))))
                       0d0))
                   0d0))
            (when (and
                   (<= (aref size-array dim-index) min-size)
                   (> (aref increment-array dim-index) 0d0))
              (setf (aref increment-array dim-index) 0d0)
              (format t "Rank ~D: Dim ~D hitting min size limit ~%" rank dim))
            (when (and
                   (<= (aref size-array (1+ dim-index)) min-size)
                   (< (aref increment-array dim-index) 0d0))
              (setf (aref increment-array dim-index) 0d0)
              (format t "Rank ~D: Dim ~D hitting min size limit ~%" (1+ rank) dim))
            ))
        (mpi-vector-sum increment-array)
        (let ((left-index (- dim-index 1)))
          (when (>= left-index 0)
            (incf (first (nth dim (mpm-sim-mpi-domain-bounds sim)))
                  (- (aref increment-array left-index)))
            ;; (format t "Rank ~D: Dim ~D moved lower by ~F ~%" (1+ rank) dim (- (aref increment-array left-index)))
            ))
        (let ((right-index dim-index))
          (when (< right-index increment-length)
            (incf (second (nth dim (mpm-sim-mpi-domain-bounds sim)))
                  (- (aref increment-array right-index)))
            ;; (format t "Rank ~D: Dim ~D moved upper by ~F ~%" (1+ rank) dim (- (aref increment-array right-index)))
            ))))
    (every (lambda (x) (= x 0d0)) increment-array)))

(defun load-balance-value (sim)
  (let ((rank (cl-mpi::mpi-comm-rank))
        (size (cl-mpi::mpi-comm-size)))
    (let* ((metric (float (mpm-sim-mpi-load-metric sim) 0d0))
           (sum-mps (mpi-sum metric))
           (max-mps (mpi-max metric))
           (balance nil))
      (when (> sum-mps 0)
        (setf balance
              (/
               max-mps
               (/ sum-mps
                  (cl-mpi:mpi-comm-size)))))
      ;; (when t;(> min-mps 0)
      ;;   ;; (setf balance (mpi-max (/ max-mps min-mps)))
      ;;   (when (= rank 0)
      ;;     (format t "Occupancy ratio : ~F%~%" (* 100d0 balance))))
      balance)
    ;; (let ((min-mps (mpi-min (float (load-balance-metric sim) 0d0)))
    ;;       (max-mps (mpi-max (float (load-balance-metric sim) 0d0)))
    ;;       (sum-mps (mpi-sum (float (load-balance-metric sim) 0d0)))
    ;;       (balance nil))
    ;;   ;; (when (= rank 0)
    ;;   ;;   (format t "Min : ~F%~%" min-mps)
    ;;   ;;   (format t "Max : ~F%~%" max-mps)
    ;;   ;;   )
    ;;   (when (> min-mps 0)
    ;;     (setf balance (mpi-max (/ max-mps min-mps)))
    ;;     (when (= rank 0)
    ;;       (format t "Occupancy ratio : ~F%~%" (* 100d0 balance))))
    ;;   ;; (format t "Rank ~D: Domain bounds ~A~%" rank (mpm-sim-mpi-domain-bounds sim))
    ;;   balance)
    ))


(defun load-balance (sim &key (substeps 10)
                           (step-size 1d-1)
                           (exchange-mps t)
                           (dims (list :x :y :z))
                           )
  (let ((rank (cl-mpi::mpi-comm-rank))
        (stagnent nil))
    (load-balance-setup sim)
    (loop for j from 0 to substeps
          while (not stagnent)
          do
             (progn
               (setf stagnent t)
               (loop for dim in dims ;(i (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
                     do
                        (let ((i (position dim (list :x :y :z))))
                          (when t;(< i (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
                            (setf
                             stagnent
                             (and
                              stagnent
                              (load-balance-dimension sim i :step-size step-size)))
                            (when exchange-mps
                              (set-mp-mpi-index sim)
                              (exchange-mps sim 0d0)
                              (set-mp-mpi-index sim)
                              (clear-ghost-mps sim))
                            (load-balance-setup sim))))))
    (let* ((metric (float (mpm-sim-mpi-load-metric sim) 0d0))
           (min-mps (mpi-min metric))
           (max-mps (mpi-max metric))
           (balance nil))
      (setf balance
            (/
             (mpi-max metric)
             (/ (mpi-sum metric)
                (cl-mpi:mpi-comm-size))))
      (when t;(> min-mps 0)
        ;; (setf balance (mpi-max (/ max-mps min-mps)))
        (when (= rank 0)
          (format t "Occupancy ratio : ~F%~%" (* 100d0 balance))))
      (format t "Rank ~D: Domain bounds ~A~%" rank (mpm-sim-mpi-domain-bounds sim))
      (values balance stagnent))))

(defun load-balance-algo (sim &key (substeps 10)
                                (max-iters 50)
                                (min-bounds 1.1d0)
                                (max-bounds 1.5d0)
                                (step-size 1d-1)
                                (dims (list :x :y :z))
                                )
  (when (typep sim 'cl-mpm/mpi:mpm-sim-mpi)
    (load-balance-setup sim)
    (let ((balance (load-balance-value sim))
          (rank (cl-mpi::mpi-comm-rank)))
      (when (= rank 0)
        (format t "Check balance~%"))
      (when (and balance (> balance max-bounds))
        (let ((balance nil)
              (stag nil))
          (loop repeat max-iters
                while (and (if balance (> balance min-bounds) t)
                           (not stag))
                do
                   (multiple-value-bind (balance stagnent) (cl-mpm/mpi::load-balance sim
                                                                                     :exchange-mps t
                                                                                     :step-size step-size
                                                                                     :substeps substeps
                                                                                     :dims dims)
                     (setf balance balance
                           stag stagnent))))
        (domain-decompose sim)))))
