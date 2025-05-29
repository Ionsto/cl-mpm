(defpackage :cl-mpm/aggregate
  (:use :cl
   :cl-mpm
   :cl-mpm/particle
   :cl-mpm/mesh
   :cl-mpm/utils
   :cl-mpm/fastmaths
   )
  ;; (:export
  ;;  #:apply-ghost)
  )
;; (declaim (optimize (debug 3) (safety 3) (speed 0)))
(declaim #.cl-mpm/settings:*optimise-setting*)

(in-package :cl-mpm/aggregate)

(defclass aggregate-element ()
  ((boundary-cell
    :initform nil
    :initarg :boundary-cell
    :accessor agg-boundary-cell)
   (interior-cell
    :initform nil
    :initarg :interior-cell
    :accessor agg-interior-cell)
   (node-list
    :initform nil
    :initarg :node-list
    :accessor agg-node-list)
   (shape-functions
    :initform nil
    :accessor agg-shape-functions)
   (mass-matrix
    :initform nil
    :accessor agg-mass-matrix
    )))

;; (defclass aggregate-node ()
;;   ((interior-cell
;;     :initform nil
;;     :initarg :interior-cell
;;     :accessor agg-interior-cell)
;;    (shape-functions
;;     :initform nil
;;     :accessor agg-shape-functions)
;;    (mass-matrix
;;     :initform nil
;;     :accessor agg-mass-matrix
;;     )))

(defclass mpm-sim-aggregated (mpm-sim)
  ((agg-elems
    :initform (make-array 0)
    :accessor sim-agg-elems
    )
   (enable-aggregate
    :initform t
    :initarg :enable-aggregate
    :accessor sim-enable-aggregate)))



(defclass mpm-sim-aggregated (mpm-sim)
  ((agg-elems
    :initform (make-array 0)
    :accessor sim-agg-elems
    )
   (enable-aggregate
    :initform t
    :initarg :enable-aggregate
    :accessor sim-enable-aggregate)))

(defun get-closest-cell (sim node)
  (let ((closest-elem nil)
        (dist 0d0)
        (mutex (sb-thread:make-mutex))
        (pos (cl-mpm/mesh::node-position node)))
    (cl-mpm::iterate-over-cells
     (cl-mpm:sim-mesh sim)
     (lambda (cell)
       (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                        (index cl-mpm/mesh::cell-index)
                        (centroid cl-mpm/mesh::cell-centroid)
                        (active cl-mpm/mesh::cell-active)
                        (ghost cl-mpm/mesh::cell-ghost-element))
           cell
         (when (and
                (cl-mpm/mesh::cell-active cell)
                (= (cl-mpm/mesh::cell-mp-count cell) 1)
                (not (cl-mpm/mesh::cell-ghost-element cell)))
           (let ((dist-tr (cl-mpm/fastmaths::diff-norm
                           pos
                           centroid)))
             (when (or
                    (not closest-elem)
                    (> dist dist-tr))
               (sb-thread:with-mutex (mutex)
                 (when (or
                        (not closest-elem)
                        (> dist dist-tr))
                   (setf dist dist-tr
                         closest-elem cell))))))
         )))
    closest-elem))

(defun locate-aggregate-node-elements (sim)
  (let ((agg-elem (list))
        (mutex (sb-thread:make-mutex)))
    (with-accessors ((mps cl-mpm:sim-mps)
                     (mesh cl-mpm:sim-mesh))
        sim
      (cl-mpm:iterate-over-nodes
       mesh
       (lambda (node)
         (setf (cl-mpm/mesh::node-agg node) nil)))

      ;;First set all outside nodes as aggregate
      (cl-mpm::iterate-over-cells
       mesh
       (lambda (cell)
         (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                          (index cl-mpm/mesh::cell-index)
                          (nodes cl-mpm/mesh::cell-nodes)
                          (active cl-mpm/mesh::cell-active))
             cell
           (when (and active (= mp-count 0))
             (loop for n in nodes
                   do (when (cl-mpm/mesh::node-active n)
                        (setf (cl-mpm/mesh::node-agg n) t)))))))

      ;;Next markup any cells that are touching aggregates as invalid as support
      (cl-mpm::iterate-over-cells
       mesh
       (lambda (cell)
         (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                          (index cl-mpm/mesh::cell-index)
                          (nodes cl-mpm/mesh::cell-nodes)
                          (ghost cl-mpm/mesh::cell-ghost-element)
                          (active cl-mpm/mesh::cell-active))
             cell
           (setf ghost nil)
           (when (and active)
             (loop for n in nodes
                   while (not ghost)
                   do (when (and (cl-mpm/mesh::node-active n)
                                 (cl-mpm/mesh::node-agg n))
                        (setf ghost t)))))))

      ;;For each aggregate node, locate the closest support cell
      (cl-mpm::iterate-over-nodes
       mesh
       (lambda (node)
         (with-accessors ((active cl-mpm/mesh::node-active)
                          (agg cl-mpm/mesh::node-agg))
             node
           (when (and active agg)
             (let ((closest-elem (get-closest-cell sim node)))
               (when closest-elem
                 (sb-thread:with-mutex (mutex)
                   (push
                    (make-instance 'aggregate-element
                                   :interior-cell closest-elem
                                   :node-list (append
                                               (list node)
                                               (cl-mpm/mesh::cell-nodes closest-elem)))
                    agg-elem)))))))))
    (make-array (length agg-elem) :initial-contents agg-elem)))

(defun locate-aggregate-elements (sim)
  (let ((agg-elem (list)))
    (with-accessors ((mps cl-mpm:sim-mps)
                     (mesh cl-mpm:sim-mesh))
        sim
      (cl-mpm:iterate-over-nodes
       mesh
       (lambda (node)
         (setf (cl-mpm/mesh::node-agg node) nil)))

      (cl-mpm::iterate-over-cells
       mesh
       (lambda (cell)
         (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                          (index cl-mpm/mesh::cell-index)
                          (nodes cl-mpm/mesh::cell-nodes)
                          (active cl-mpm/mesh::cell-active)
                          )
             cell
           (when (and active (= mp-count 0))
             (loop for n in nodes
                   do (when (cl-mpm/mesh::node-active n)
                        (setf (cl-mpm/mesh::node-agg n) t)))))))

      (cl-mpm::iterate-over-cells
       mesh
       (lambda (cell)
         (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                          ;; (neighbours cl-mpm/mesh::cell-neighbours)
                          (neighbours cl-mpm/mesh::cell-cartesian-neighbours)
                          (index cl-mpm/mesh::cell-index)
                          (nodes cl-mpm/mesh::cell-nodes)
                          (pruned cl-mpm/mesh::cell-pruned)
                          (active cl-mpm/mesh::cell-active)
                          (ghost cl-mpm/mesh::cell-ghost-element)
                          (pos cl-mpm/mesh::cell-centroid))
             cell
           (setf (cl-mpm/mesh::cell-ghost-element cell) nil)
           (when (and active (= mp-count 1))
             (loop for neighbour in neighbours
                   do
                      (when (cl-mpm/mesh::cell-active neighbour)
                        (when (= (cl-mpm/mesh::cell-mp-count neighbour) 0)
                          (setf (cl-mpm/mesh::cell-ghost-element cell) t))))))))

      (cl-mpm::iterate-over-cells-serial
       mesh
       (lambda (cell)
         (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                          (neighbours cl-mpm/mesh::cell-neighbours)
                          (index cl-mpm/mesh::cell-index)
                          (centroid cl-mpm/mesh::cell-centroid)
                          (nodes cl-mpm/mesh::cell-nodes)
                          (pruned cl-mpm/mesh::cell-pruned)
                          (active cl-mpm/mesh::cell-active)
                          (ghost cl-mpm/mesh::cell-ghost-element)
                          (pos cl-mpm/mesh::cell-centroid))
             cell
           (when ghost
             (let ((closest-elem nil)
                   (dist 0d0))
               (loop for neighbour in neighbours
                     do
                        (when (and
                               (cl-mpm/mesh::cell-active neighbour)
                               (= (cl-mpm/mesh::cell-mp-count neighbour) 1)
                               (not (cl-mpm/mesh::cell-ghost-element neighbour)))
                          (let ((dist-tr (cl-mpm/fastmaths::diff-norm
                                          centroid
                                          (cl-mpm/mesh::cell-centroid neighbour))))
                            (when (or
                                   (not closest-elem)
                                   (> dist dist-tr))
                              (setf dist dist-tr
                                    closest-elem neighbour)))))
               (when closest-elem
                 (setf (cl-mpm/mesh::cell-agg-int cell)
                       ;; 1
                       (nth 0 (cl-mpm/mesh::cell-index closest-elem))
                       )
                 (push
                  (make-instance 'aggregate-element
                                 :boundary-cell cell
                                 :interior-cell closest-elem
                                 :node-list (make-elem-node-list closest-elem cell))
                  agg-elem))))))))
    (make-array (length agg-elem) :initial-contents agg-elem)))

(defun iterate-over-agg-elem (agg-elems func)
  (declare (function func))
  (lparallel:pdotimes (i (length agg-elems))
    (funcall func (aref agg-elems i))))

(defun set-aggregate-nodes (sim)
  (with-accessors ((agg-elems sim-agg-elems)
                   (mesh cl-mpm::sim-mesh))
      sim
    ;; (cl-mpm:iterate-over-nodes
    ;;  mesh
    ;;  (lambda (node)
    ;;    (setf (cl-mpm/mesh::node-agg node) nil)))
    ;; (iterate-over-agg-elem
    ;;  agg-elems
    ;;  (lambda (agg-elem)
    ;;    (dolist (cell (list (agg-boundary-cell agg-elem)))
    ;;      (loop for n in (cl-mpm/mesh::cell-nodes cell)
    ;;            do
    ;;                                     ;(position n (cl-mpm/mesh::cell-nodes (agg-interior-cell agg-elem)))
    ;;               (when (cl-mpm/mesh::node-agg n)
    ;;                 (setf (cl-mpm/mesh::node-agg n) t))))))
    ))

(defun accumulate-over-shape-function (mesh pos)
  )


(defun find-local-coords-agg (mesh cell pos)
  (let ((tol 1d-9)
        (iters 0)
        (gp-loc (cl-mpm/utils:vector-zeros)))
    (loop for i from 0 to 2
          do
             (let ((tr-loc (cl-mpm/utils:vector-zeros))
                   (grad-loc (cl-mpm/utils:matrix-zeros)))
               (cl-mpm::iterate-over-cell-linear-3d
                mesh
                cell
                pos
                (lambda (node weight grads)
                  (cl-mpm/fastmaths:fast-fmacc tr-loc (cl-mpm/mesh::node-position node) weight)
                  (destructuring-bind (dx dy dz) grads
                    (fast-.+
                     (magicl:@
                      (vector-from-list (list dx dy dz))
                      (magicl:transpose (cl-mpm/mesh::node-position node)))
                     grad-loc
                     grad-loc))))
               (when (= (cl-mpm/mesh:mesh-nd mesh) 2)
                 (setf (mtref grad-loc 2 2) 1d0))
               (format t "Tr ~A~%" tr-loc)
               (cl-mpm/fastmaths:fast-.+
                (magicl:linear-solve grad-loc (cl-mpm/fastmaths:fast-.- pos tr-loc))
                gp-loc
                gp-loc)))
    (values gp-loc)))

(defun weight-local-pos (x y z nd dx dy dz)
  (let ((dx (- (* dx 2) 1))
        (dy (- (* dy 2) 1))
        (dz (- (* dz 2) 1)))
    (case nd
      (2 (* 0.25d0  (+ 1d0 (* dx x)) (+ 1d0 (* dy y))))
      (3 (* 0.125d0 (+ 1d0 (* dx x)) (+ 1d0 (* dy y)) (+ 1d0 (* dz z)))))))

(defun iterate-over-cell-shape-local (mesh cell local-position func)
  "Iterating over a given cell's basis functions"
  (declare (cl-mpm/mesh::mesh mesh))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (nd (cl-mpm/mesh::mesh-nd mesh))
           (pos-vec local-position)
           (cell-pos (cl-mpm/mesh::cell-centroid cell))
           (cell-index (cl-mpm/mesh::cell-index cell)))
      (declare (dynamic-extent cell-index pos-vec))
      (loop for dx from 0 to 1
            do (loop for dy from 0 to 1
                     do (loop for dz from 0 to 1
                              do (let* ((id (mapcar #'+ cell-index (list dx dy dz))))
                                   (declare (dynamic-extent id))
                                   (when (cl-mpm/mesh:in-bounds mesh id)
                                     (let* ((node (cl-mpm/mesh:get-node mesh id))
                                            (grads (list 0d0 0d0 0d0))
                                            (dist-x (* 2d0 (/ (- (varef pos-vec 0) (varef cell-pos 0)) h)))
                                            (dist-y (* 2d0 (/ (- (varef pos-vec 1) (varef cell-pos 1)) h)))
                                            (dist-z (* 2d0 (/ (- (varef pos-vec 2) (varef cell-pos 2)) h)))
                                            )
                                       (let* ((weight (weight-local-pos
                                                       dist-x
                                                       dist-y
                                                       dist-z
                                                       nd
                                                       dx
                                                       dy
                                                       dz)))
                                         (declare
                                          (double-float weight))
                                         (when t
                                           (funcall func node weight grads))))))))))))


(defun make-elem-node-list (boundary-cell interior-cell)
  (let ((node-list nil))
    (dolist (cell (list boundary-cell interior-cell))
      (dolist (n (cl-mpm/mesh::cell-nodes cell))
        (push n node-list)))
    (remove-duplicates node-list)))
(defun iterate-over-agg-elem-nodes (sim elem func)
  (declare (function func))
  (loop for n in (agg-node-list elem)
        do (funcall func n)))

(defun compute-extension-matrix (sim elem)
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (let* ((nd (cl-mpm/mesh:mesh-nd mesh))
           (nc (length (agg-node-list elem)))
           (E (cl-mpm/utils::arb-matrix (* nd nc)
                                        (* nd (expt 2 nd)))))
      (let ((iter 0))
        (iterate-over-agg-elem-nodes
         sim
         elem
         (lambda (n)
           (let ((rank 0))
             (iterate-over-cell-shape-local
              mesh
              (agg-interior-cell elem)
              (cl-mpm/mesh::node-position n)
              (lambda (cn weight grads)
                (loop for i from 0 below nd
                      do (setf (mtref E (+ iter i) (+ rank i)) weight))
                (incf rank nd)))
             (incf iter nd)))))
      (setf (agg-shape-functions elem) E))))

(defun update-aggregate-elements (sim)
  (setf
   (sim-agg-elems sim)
   (locate-aggregate-node-elements sim))
  ;; (set-aggregate-nodes sim)
  (iterate-over-agg-elem
   (sim-agg-elems sim)
   (lambda (elem)
     (compute-extension-matrix sim elem)
     (setf (agg-mass-matrix elem) (assemble-mass sim elem)))))


(defun assemble-mass (sim elem)
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (nc (length (agg-node-list elem)))
         (m (cl-mpm/utils::arb-matrix
             (* nd nc)
             (* nd nc)))
         (E (agg-shape-functions elem)))
    (let ((iter 0))
      (iterate-over-agg-elem-nodes
       sim
       elem
       (lambda (n)
         (let ((mi (cl-mpm/mesh::node-mass n)))
           (loop for i from 0 below nd
                 do (progn
                      (setf (mtref M iter iter) mi)
                      (incf iter)))))))
    (magicl:@
     (magicl:transpose E)
     m
     E)
    ))
(defun assemble-full-mass (sim elem)
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (nc (length (agg-node-list elem)))
         (m (cl-mpm/utils::arb-matrix
             (* nd nc)
             (* nd nc)))
         (E (agg-shape-functions elem)))
    (let ((iter 0))
      (iterate-over-agg-elem-nodes
       sim
       elem
       (lambda (n)
         (let ((mi (cl-mpm/mesh::node-mass n)))
           (loop for i from 0 below nd
                 do (progn
                      (setf (mtref M iter iter) mi)
                      (incf iter)))))))
    m))
(defun assemble-vector (sim elem accessor)
  (declare (function accessor))
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (nc (length (agg-node-list elem)))
         (f (cl-mpm/utils::arb-matrix (* nd nc) 1)))
    (let ((iter 0))
      (iterate-over-agg-elem-nodes
       sim
       elem
       (lambda (n)
         (let ((force (funcall accessor n)))
           (loop for i from 0 below nd
                 do (setf (mtref f (+ iter i) 0) (varef force i))))
         (incf iter nd))))
    f))
(defun assemble-force (sim elem)
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (nc (length (agg-node-list elem)))
         (f (cl-mpm/utils::arb-matrix
             (* nd nc)
             1)))
    (let ((iter 0))
      (iterate-over-agg-elem-nodes
       sim
       elem
       (lambda (n)
         (with-accessors ((force cl-mpm/mesh::node-force))
             n
           (loop for i from 0 below nd
                 do (setf (mtref f (+ iter i) 0) (varef force i)))
           (incf iter nd)))))
    f))
(defun assemble-damping-force (sim elem damping-factor)
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (nc (length (agg-node-list elem)))
         (f (cl-mpm/utils::arb-matrix
             (* nd nc)
             1)))
    (let ((iter 0))
      (iterate-over-agg-elem-nodes
       sim
       elem
       (lambda (n)
         (with-accessors ((vel cl-mpm/mesh::node-velocity))
             n
           (loop for i from 0 below nd
                 do (setf (mtref f (+ iter i) 0)
                          (*
                           (varef vel i)
                           damping-factor
                           -1d0)))
           (incf iter nd)))))
    f))
(defun assemble-residual (sim elem)
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (nc (length (agg-node-list elem)))
         (f (cl-mpm/utils::arb-matrix
             (* nd nc)
             1)))
    (let ((iter 0))
      (iterate-over-agg-elem-nodes
       sim
       elem
       (lambda (n)
         (with-accessors ((force cl-mpm/mesh::node-residual))
             n
           (loop for i from 0 below nd
                 do (setf (mtref f (+ iter i) 0) (varef force i)))
           (incf iter nd)))))
    f))

(defun project-acc (sim elem acc)
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (E (agg-shape-functions elem))
         (acc (magicl:@ E acc)))
    (let ((iter 0))
      (iterate-over-agg-elem-nodes
       sim
       elem
       (lambda (n)
         (when (cl-mpm/mesh::node-agg n)
           (with-accessors ((n-acc cl-mpm/mesh::node-acceleration))
               n
             (loop for i from 0 below nd
                   do (setf (varef n-acc i) (mtref acc (+ iter i) 0)))
             (incf iter nd)))))
      )))
(defun project-vel (sim elem acc)
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (E (agg-shape-functions elem))
         (acc (magicl:@ E acc)))
    (let ((iter 0))
      (iterate-over-agg-elem-nodes
       sim
       elem
       (lambda (n)
         (when (cl-mpm/mesh::node-agg n)
           (with-accessors ((n-acc cl-mpm/mesh::node-velocity))
               n
             (loop for i from 0 below nd
                   do (setf (varef n-acc i) (mtref acc (+ iter i) 0)))
             (incf iter nd))))))))

(defun project-sub-forces (sim elem f-int f-ext f-total)
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (E (agg-shape-functions elem))
         (f-int (magicl:@ E (magicl:transpose E) f-int))
         (f-ext (magicl:@ E (magicl:transpose E) f-ext))
         ;; (f-damp (magicl:@ E (magicl:transpose E) f-damp))
         )
    (let ((iter 0))
      (iterate-over-agg-elem-nodes
       sim
       elem
       (lambda (n)
         (when (cl-mpm/mesh::node-agg n)
           (with-accessors ((n-force-int cl-mpm/mesh::node-internal-force)
                            (n-force-ext cl-mpm/mesh::node-external-force)
                            (n-force-damp cl-mpm/mesh::node-damping-force)
                            (n-force cl-mpm/mesh::node-force)
                            )
               n
             (loop for i from 0 below nd
                   do (setf (varef n-force-int i) (mtref f-int (+ iter i) 0)
                            (varef n-force-ext i) (mtref f-ext (+ iter i) 0)
                            (varef n-force i) (mtref f-total (+ iter i) 0)
                            )
                      ;; (setf (varef))
                      ;; (incf (varef n-force i) (mtref f-damp (+ iter i) 0))
                   )
             (incf iter nd))))))))


(defun project-residual (sim elem residual)
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (E (agg-shape-functions elem))
         (force (magicl:@ E (magicl:transpose E) residual))
         )
    (let ((iter 0))
      (iterate-over-agg-elem-nodes
       sim
       elem
       (lambda (n)
         (when (cl-mpm/mesh::node-agg n)
           (with-accessors ((n-force cl-mpm/mesh::node-residual))
               n
             (loop for i from 0 below nd
                   do (setf (varef n-force i) (mtref force (+ iter i) 0)))
             (incf iter nd)))))
      )))

(defun calculate-forces-agg-elem (sim elem damping)
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (E (agg-shape-functions elem))
         (m (agg-mass-matrix elem))
         (f (assemble-force sim elem))
         (f-int (assemble-vector sim elem #'cl-mpm/mesh::node-internal-force))
         (f-ext (assemble-vector sim elem #'cl-mpm/mesh::node-external-force))
         ;; (f-damp (magicl:@ M (magicl:transpose E) (assemble-damping-force sim elem damping)))
         ;; (res (assemble-residual sim elem))
         )
    ;; (break)
    ;; (cl-mpm/fastmaths:fast-fmacc force-damp vel (* damping -1d0 mass))
    (let ((acc (magicl:linear-solve m (magicl:@ (magicl:transpose E)
                                                f
                                                ))))
      (project-acc sim elem acc)
      ;; (project-damping sim elem f-int f-ext)
      (project-sub-forces sim elem f-int f-ext (magicl:@ E (magicl:transpose E) f))
      )))

(declaim (notinline update-node-kinematics))
(defmethod update-node-kinematics ((sim mpm-sim-aggregated))
  ;;For non-aggregate nodes, use simple mass matrix inversion
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (agg-elems sim-agg-elems))
      sim
    (cl-mpm:iterate-over-nodes
     mesh
     (lambda (node)
       (unless (cl-mpm/mesh::node-agg node)
         (cl-mpm::calculate-kinematics node))))
    ;;For each aggregated element set solve mass matrix and velocity
    (iterate-over-agg-elem
     agg-elems
     (lambda (elem)
       (calculate-kinematics-agg-elem sim elem)))))

(defmethod cl-mpm::update-node-forces ((sim mpm-sim-aggregated))
  ;;For non-aggregate nodes, use simple mass matrix inversion
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mass-scale sim-mass-scale)
                   (damping sim-damping-factor)
                   (damping-algo sim-damping-algorithm)
                   (agg-elems sim-agg-elems)
                   (enable-aggregate sim-enable-aggregate)
                   (dt sim-dt))
      sim
    (iterate-over-nodes
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (if (cl-mpm/mesh::node-agg node)
             (cl-mpm::calculate-forces node damping 0d0 mass-scale)
             (cl-mpm::calculate-forces node damping dt mass-scale)))))
    ;;For each aggregated element set solve mass matrix and velocity
    (when enable-aggregate
      (iterate-over-agg-elem
       agg-elems
       (lambda (elem)
         (calculate-forces-agg-elem sim elem damping))))
    ;; (break)
    (iterate-over-nodes
     mesh
     (lambda (node)
       (when (and (cl-mpm/mesh:node-active node)
                  (cl-mpm/mesh::node-agg node))
         (cl-mpm/fastmaths:fast-fmacc
          (cl-mpm/mesh::node-velocity node)
          (cl-mpm/mesh::node-acceleration node) dt))))
    ;; (when enable-aggregate
    ;;   (iterate-over-agg-elem
    ;;    agg-elems
    ;;    (lambda (elem)
    ;;      (let* ((E (agg-shape-functions elem))
    ;;             (v-i (assemble-vector sim elem #'cl-mpm/mesh::node-velocity))
    ;;             (m-i (assemble-full-mass sim elem))
    ;;             (m (agg-mass-matrix elem)))
    ;;        ;; (break)
    ;;        (project-vel sim elem (magicl:linear-solve m (magicl:@ (magicl:transpose E) (magicl:@ m-i v-i))))))))
    ))


