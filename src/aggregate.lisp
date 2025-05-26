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
(declaim (optimize (debug 3) (safety 3) (speed 0)))
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
   (shape-functions
    :initform nil
    :accessor agg-shape-functions)
   (mass-matrix
    :initform nil
    :accessor agg-mass-matrix
    ))
  )
(defclass mpm-sim-aggregated (mpm-sim)
  ((agg-elems
    :initform (make-array 0)
    :accessor sim-agg-elems
    ))
  )


(defun locate-aggregate-elements (sim)
  (let ((agg-elem (list)))
    (with-accessors ((mps cl-mpm:sim-mps)
                     (mesh cl-mpm:sim-mesh))
        sim
      (let ((cells (cl-mpm/mesh::mesh-cells mesh)))
        (cl-mpm::iterate-over-cells-serial
         mesh
         (lambda (cell)
           (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                            (neighbours cl-mpm/mesh::cell-neighbours)
                            (index cl-mpm/mesh::cell-index)
                            (nodes cl-mpm/mesh::cell-nodes)
                            (pruned cl-mpm/mesh::cell-pruned)
                            (active cl-mpm/mesh::cell-active)
                            (ghost cl-mpm/mesh::cell-ghost-element)
                            (pos cl-mpm/mesh::cell-centroid))
               cell
             ;; (setf ghost nil)
             (when (and (= mp-count 1))
               ;; (setf ghost t)
               (loop for neighbour in neighbours
                     do
                        (when (cl-mpm/mesh::cell-active neighbour)
                          (when (= (cl-mpm/mesh::cell-mp-count neighbour) 0)
                            ;;Make aggregate element
                            (push
                             (make-instance 'aggregate-element
                                             :boundary-cell cell
                                             :interior-cell neighbour)
                             agg-elem))))))))))
    (make-array (length agg-elem) :initial-contents agg-elem)))

(defun iterate-over-agg-elem (agg-elems func)
  (declare (function func))
  (lparallel:pdotimes (i (length agg-elems))
    (funcall func (aref agg-elems i))))

(defun set-aggregate-nodes (sim)
  (with-accessors ((agg-elems sim-agg-elems))
      sim
      (iterate-over-agg-elem
       agg-elems
       (lambda (agg-elem)
         (dolist (cell (list (agg-boundary-cell agg-elem) (agg-interior-cell agg-elem)))
           (loop for n in (cl-mpm/mesh::cell-nodes cell)
                 do (setf (cl-mpm/mesh::node-agg n) t)))))))

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

(defun iterate-over-cell-shape-local (mesh cell local-position func)
  "Iterating over a given cell's basis functions"
  (declare (cl-mpm/mesh::mesh mesh))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (pos-vec local-position)
           (pos (list (varef pos-vec 0) (varef pos-vec 1) (varef pos-vec 2)))
           (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec #'floor))
           (cell-index (cl-mpm/mesh::cell-index cell))
           )
      (declare (dynamic-extent pos pos-index pos-vec))
      (loop for dx from 0 to 1
            do (loop for dy from 0 to 1
                     do (loop for dz from 0 to 1
                              do (let* ((id (mapcar #'+ cell-index (list dx dy dz))))
                                   (declare (dynamic-extent id))
                                   (when (cl-mpm/mesh:in-bounds mesh id)
                                     (let* ((dist (mapcar #'- pos (list dx dy dz)))
                                            (node (cl-mpm/mesh:get-node mesh id))
                                            (weights (mapcar (lambda (x) (cl-mpm/shape-function::shape-linear x h)) dist))
                                            (weight (reduce #'* weights))
                                            (lin-grads (mapcar (lambda (d)
                                                                 (cl-mpm/shape-function::shape-linear-dsvp d h))
                                                               dist))
                                            (grads (cl-mpm/shape-function::grads-3d weights lin-grads)))
                                       (declare
                                        (double-float weight)
                                        (dynamic-extent dist weights))
                                       (when t
                                         (funcall func node weight grads)))))))))))

(defun compute-extension-matrix (sim elem)
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (let* ((nd (cl-mpm/mesh:mesh-nd mesh))
           (E (cl-mpm/utils::arb-matrix (* 3 2 (expt 2 nd) ) (expt 2 nd))))
      (let ((iter 0))
        (dolist (cell (list (agg-boundary-cell elem)
                            (agg-interior-cell elem)))
          (dolist (n (cl-mpm/mesh::cell-nodes cell))
            (let ((rank 0))
              (iterate-over-cell-shape-local
               mesh
               (agg-interior-cell elem)
               (cl-mpm/mesh::node-position n)
               (lambda (cn weight grads)
                 (setf (mtref E iter rank) weight)
                 (setf (mtref E (+ iter 1) rank) weight)
                 (setf (mtref E (+ iter 2) rank) weight)
                 (incf rank))))
            (incf iter 3))))
      (setf (agg-shape-functions elem) E))))
(defun update-aggregate-elements (sim)

  
  )

(defun assemble-mass (sim elem)
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (m (cl-mpm/utils::arb-matrix
             (* 3 2 (expt 2 nd))
             (* 3 2 (expt 2 nd))))
         (E (agg-shape-functions elem)))
    (let ((iter 0))
      (dolist (cell (list (agg-boundary-cell elem)
                          (agg-interior-cell elem)))
        (dolist (n (cl-mpm/mesh::cell-nodes cell))
          (setf (mtref m iter iter) (cl-mpm/mesh::node-mass n))
          (incf iter)
          (setf (mtref m iter iter) (cl-mpm/mesh::node-mass n))
          (incf iter)
          (setf (mtref m iter iter) (cl-mpm/mesh::node-mass n))
          (incf iter))))
    ;; (pprint m)
    ;; (pprint E)
    (magicl:@
     (magicl:transpose E)
     m
     E
     )))
(defun assemble-force (sim elem)
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (f (cl-mpm/utils::arb-matrix
             (* 3 2 (expt 2 nd))
             1))
         (E (agg-shape-functions elem)))
    (let ((iter 0))
      (dolist (cell (list (agg-boundary-cell elem)
                          (agg-interior-cell elem)))
        (dolist (n (cl-mpm/mesh::cell-nodes cell))
          (with-accessors ((vel cl-mpm/mesh::node-force))
              n
            (setf (mtref f (+ iter 0) 0) (varef vel 0)
                  (mtref f (+ iter 1) 0) (varef vel 1)
                  (mtref f (+ iter 2) 0) (varef vel 2))
            (incf iter 3)))))
    f
    ;; (magicl:@
    ;;  (magicl:transpose E)
    ;;  m
    ;;  E
    ;;  )
    ))
(defun assemble-vel (sim elem)
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (v (cl-mpm/utils::arb-matrix
             (* 3 2 (expt 2 nd))
             1))
         (E (agg-shape-functions elem)))
    (let ((iter 0))
      (dolist (cell (list (agg-boundary-cell elem)
                          (agg-interior-cell elem)))
        (dolist (n (cl-mpm/mesh::cell-nodes cell))
          (with-accessors ((vel cl-mpm/mesh::node-velocity))
              n
            (setf (mtref v (+ iter 0) 0) (varef vel 0)
                  (mtref v (+ iter 1) 0) (varef vel 1)
                  (mtref v (+ iter 2) 0) (varef vel 2))
            (incf iter 3)))))
    (magicl:@
     (magicl:transpose E)
     m
     E
     )))
(defun project-acc (sim elem acc)
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (E (agg-shape-functions elem)))
    (let ((iter 0))
      (dolist (cell (list (agg-boundary-cell elem)
                          (agg-interior-cell elem)))
        (dolist (n (cl-mpm/mesh::cell-nodes cell))
          (with-accessors ((n-acc cl-mpm/mesh::node-acceleration))
              n
            (setf (varef n-acc 0) (mtref acc (+ iter 0) 0)
                  (varef n-acc 1) (mtref acc (+ iter 1) 0)
                  (varef n-acc 2) (mtref acc (+ iter 2) 0))
            (incf iter 3)))))))

(defun calculate-kinematics-agg-elem (sim elem)
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (acc (cl-mpm/utils::arb-matrix
               (* 3 2 (expt 2 nd))
               1))
         (E (agg-shape-functions elem))
         (m (assemble-mass sim elem))
         (f (assemble-force sim elem))
         )
    (let ((acc (magicl:linear-solve m f)))
      (project-acc sim elem acc))))

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
       (calculate-kinematics-agg-elem sim elem)
       ))))


