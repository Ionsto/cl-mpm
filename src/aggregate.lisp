(defpackage :cl-mpm/aggregate
  (:use :cl
   :cl-mpm
   :cl-mpm/particle
   :cl-mpm/mesh
   :cl-mpm/utils
   :cl-mpm/fastmaths
   )
  (:export
   #:mpm-sim-agg-usf)
  )
;; (declaim (optimize (debug 3) (safety 3) (speed 0)))
(declaim #.cl-mpm/settings:*optimise-setting*)

(in-package :cl-mpm/aggregate)

(defmacro project-global-vec (sim force accessor)
  `(progn
     ;; (declare (function accessor))
     (let* ((agg-nodes (filter-nodes ,sim #'cl-mpm/mesh::node-interior))
            (active-nodes (filter-nodes ,sim (lambda (n) (or (cl-mpm/mesh::node-agg n)
                                                             (cl-mpm/mesh::node-interior n)))))
            (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh ,sim)))
            (proj-val ,force)
            )
       (iterate-over-agg-elem
        (sim-agg-elems ,sim)
        (lambda (elem)
          (let* ((nc (length (agg-node-list elem))))
            (iterate-over-agg-elem-nodes
             ,sim
             elem
             (lambda (n)
               (let ((index (* nd (position n active-nodes))))
                 (with-accessors ((vec ,accessor))
                     n
                   (loop for i from 0 below nd
                         do (setf
                             (varef vec i)
                             (mtref proj-val (+ index i) 0)))))))))))))

(defclass aggregate-element ()
  ((boundary-cell
    :initform nil
    :initarg :boundary-cell
    :accessor agg-boundary-cell)
   (interior-cell
    :initform nil
    :initarg :interior-cell
    :accessor agg-interior-cell)
   (agg-node
    :initform nil
    :initarg :agg-node
    :accessor agg-node)
   (node-list
    :initform nil
    :initarg :node-list
    :accessor agg-node-list)
   (shape-functions
    :initform nil
    :accessor agg-shape-functions)
   (incremental-shape-functions
    :initform nil
    :accessor agg-inc-shape-functions)
   (contribution-list
    :initform nil
    :accessor agg-contribution-list)
   (internal-selector
    :initform nil
    :accessor agg-internal-selector)
   (mass-matrix
    :initform nil
    :accessor agg-mass-matrix)
   (init-mass-matrix
    :initform nil
    :accessor agg-init-mass-matrix)
   (vel-projection-matrix
    :initform nil
    :accessor agg-vel-projection-matrix)))



(defclass mpm-sim-agg-usf (cl-mpm::mpm-sim-usf cl-mpm/aggregate::mpm-sim-aggregated)
  ()
  (:documentation "Explicit simulation with update stress first update"))



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
    :accessor sim-enable-aggregate)
   (global-fd
    :initform nil
    :accessor sim-global-fd
    )
   (global-e
    :initform nil
    :accessor sim-global-e
    )
   (global-ma
    :initform nil
    :accessor sim-global-ma
    )
   )
  )



(defclass mpm-sim-aggregated (mpm-sim)
  ((agg-elems
    :initform (make-array 0)
    :accessor sim-agg-elems)
   (enable-aggregate
    :initform t
    :initarg :enable-aggregate
    :accessor sim-enable-aggregate)))


(defmacro project-vector-internal (sim elem accessor vector)
  `(let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh ,sim)))
          (contr (agg-contribution-list ,elem))
          )
     (let ((computed-v ,vector))
       (let ((iter 0))
         (declare (fixnum iter))
         (iterate-over-agg-elem-nodes
          ,sim
          ,elem
          (lambda (n)
            (with-accessors ((n-acc ,accessor))
                n
              (when (not (cl-mpm/mesh::node-agg n))
                (loop for i from 0 below nd
                      when (aref contr (+ iter i))
                      do (setf (the double-float (varef n-acc i)) (the double-float (mtref computed-v (+ iter i) 0))))
                (incf iter nd)))))))))


(defmacro project-vector (sim elem accessor vector)
  `(let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh ,sim))))
     (let ((computed-v ,vector))
       (let ((iter 0))
         (declare (fixnum iter))
         (iterate-over-agg-elem-nodes
          ,sim
          ,elem
          (lambda (n)
            (with-accessors ((n-acc ,accessor))
                n
              (when t;(cl-mpm/mesh::node-agg n)
                (loop for i from 0 below nd
                      do (setf (the double-float (varef n-acc i)) (the double-float (varef computed-v (+ iter i))))))
              (incf iter nd))))))))

(defmacro inc-project-vector (sim elem accessor vector)
  `(let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh ,sim))))
     (let ((computed-v ,vector))
       (let ((iter 0))
         (declare (fixnum iter))
         (iterate-over-agg-elem-nodes
          ,sim
          ,elem
          (lambda (n)
            (with-accessors ((n-acc ,accessor)
                             (lock cl-mpm/mesh::node-lock))
                n
              (sb-thread:with-mutex (lock)
                (loop for i from 0 below nd
                      do (incf (the double-float (varef n-acc i)) (the double-float (varef computed-v (+ iter i))))))
              (incf iter nd))))))))

(defmacro reproject-vector (sim elem accessor)
  `(let* ((E (agg-shape-functions ,elem))
          ;; (proj (agg-internal-selector ,elem))
          (d-i (assemble-vector ,sim ,elem #',accessor)))
     (project-vector
      ,sim
      ,elem
      ,accessor
      (magicl:@
       E
       (magicl:transpose E)
       d-i))))
(defmacro reproject-inc-vector (sim elem accessor)
  `(let* ((Einc (agg-inc-shape-functions ,elem))
          (E (agg-shape-functions ,elem))
          (d-i (assemble-vector ,sim ,elem #',accessor)))
     (inc-project-vector
      ,sim
      ,elem
      ,accessor
      (magicl:@
       E
       (magicl:transpose Einc)
       d-i))))

(defun get-closest-cell (sim node)
  (let ((closest-elem nil)
        (dist 0d0)
        ;; (mutex (sb-thread:make-mutex))
        (pos (cl-mpm/mesh::node-position node)))
    (cl-mpm::iterate-over-cells-serial
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
               ;; (sb-thread:with-mutex (mutex))
               (when (or
                      (not closest-elem)
                      (> dist dist-tr))
                 (setf dist dist-tr
                       closest-elem cell)))))
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
         (setf (cl-mpm/mesh::node-interior node) nil)
         (setf (cl-mpm/mesh::node-agg node) nil)))
      (cl-mpm::iterate-over-cells
       mesh
       (lambda (node)
         (setf (cl-mpm/mesh::cell-aggregate-element node) nil)))

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
             ;;Remove boundary cells
             ;; (loop for d from 0 below 2
             ;;       do (let ((di (list 0 0 0)))
             ;;            (setf (nth d di) 1)
             ;;            (when (or (not (cl-mpm/mesh::in-bounds-cell mesh (mapcar #'+ index di)))
             ;;                      (not (cl-mpm/mesh::in-bounds-cell mesh (mapcar #'- index di))))
             ;;              (setf ghost t))))
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
               (if closest-elem
                   (if (cl-mpm/mesh::cell-aggregate-element closest-elem)
                       (progn
                         ;; (push node (agg-node-list (cl-mpm/mesh::cell-aggregate-element closest-elem)))
                         (setf (agg-node-list (cl-mpm/mesh::cell-aggregate-element closest-elem))
                               (append
                                (agg-node-list (cl-mpm/mesh::cell-aggregate-element closest-elem))
                                (list node))))
                       (progn

                         (setf
                          (cl-mpm/mesh::cell-aggregate-element closest-elem)
                          (make-instance 'aggregate-element
                                         :interior-cell closest-elem
                                         :node-list (make-elem-node-list mesh node closest-elem)))
                         (sb-thread:with-mutex (mutex)
                           (push
                            (cl-mpm/mesh::cell-aggregate-element closest-elem)
                            agg-elem))))
                 (error "No closest elem?")
                 )))))))
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

(defun iterate-over-internal-cell (mesh elem func)
  (declare (cl-mpm/mesh::mesh mesh)
           (function func))
  (progn
    (let* ((cell (agg-interior-cell elem))
           (cell-index (cl-mpm/mesh::cell-index cell)))
      (declare (dynamic-extent cell-index))
      (loop for dx from 0 to 1
            do (loop for dy from 0 to 1
                     do (loop for dz from 0 to 1
                              do (let* ((id (mapcar #'+ cell-index (list dx dy dz))))
                                   (declare (dynamic-extent id))
                                   (when (cl-mpm/mesh:in-bounds mesh id)
                                     (let* ((node (cl-mpm/mesh:get-node mesh id)))
                                       (when t
                                         (funcall func node))))))))))
  )

(defun iterate-over-cell-shape-local (mesh cell local-position func)
  "Iterating over a given cell's basis functions"
  (declare (cl-mpm/mesh::mesh mesh)
           (function func))
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

(defun make-elem-node-list (mesh agg-node interior-cell)
  (let ((node-list
          (list)
          ;; (list agg-node)
                   ))
    (iterate-over-cell-shape-local
     mesh
     interior-cell
     (cl-mpm/mesh::cell-centroid interior-cell)
     (lambda (cn weight grads)
       (push cn node-list)))
    (setf node-list (reverse node-list))
    ;; (push agg-node node-list)
    ;; (nreverse node-list)
    (append node-list (list agg-node))
    ))

(defun iterate-over-agg-elem-agg-nodes (sim elem func)
  (declare (function func))
  (funcall func (agg-node elem)))

(defun iterate-over-agg-elem-nodes (sim elem func)
  (declare (function func))
  (loop for n in (agg-node-list elem)
        do (funcall func n)))

(defun compute-extension-matrix (sim elem)
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (let* ((nd (cl-mpm/mesh:mesh-nd mesh))
           (nc (length (agg-node-list elem)))
           (contrib (make-array (expt 2 nd) :initial-element nil))
           (E (cl-mpm/utils::arb-matrix (* nd nc)
                                        (* nd (expt 2 nd)))))
      (let ((iter 0)
            (ic 0))
        (iterate-over-agg-elem-nodes
         sim
         elem
         (lambda (n)
           (when t;(not (cl-mpm/mesh::node-agg n))
             (let ((rank 0))
               (iterate-over-cell-shape-local
                mesh
                (agg-interior-cell elem)
                (cl-mpm/mesh::node-position n)
                (lambda (cn weight grads)
                  (loop for i from 0 below nd
                        do (progn
                             (setf (mtref E (+ iter i) (+ rank i)) weight)))
                  (incf rank nd)))
               (incf iter nd))))))
      (setf (agg-contribution-list elem) contrib)
      (setf (agg-shape-functions elem) E))))

(defun compute-inc-extension-matrix (sim elem)
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
                (when (cl-mpm/mesh::node-agg n)
                  (loop for i from 0 below nd
                        do (setf (mtref E (+ iter i) (+ rank i)) weight)))
                (incf rank nd)))
             (incf iter nd)))))
      (setf (agg-inc-shape-functions elem) E))))



(defun compute-identity (sim elem)
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (let* ((nd (cl-mpm/mesh:mesh-nd mesh))
           (nc (length (agg-node-list elem)))
           (E (cl-mpm/utils::arb-matrix (* nd (expt 2 nd))
                                        (* nd (expt 2 nd)))))
      (let ((iter 0))
        (iterate-over-agg-elem-nodes
         sim
         elem
         (lambda (n)
           (when (not (cl-mpm/mesh::node-agg n))
             (let ((rank 0))
               (iterate-over-cell-shape-local
                mesh
                (agg-interior-cell elem) (cl-mpm/mesh::node-position n)
                (lambda (cn weight grads)
                  (loop for i from 0 below nd
                        do (setf (mtref E (+ iter i) (+ rank i)) weight))
                  (incf rank nd)))
               (incf iter nd))))))
      (setf (agg-shape-functions elem) E))))

(defun compute-internal-expander (sim elem)
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
                      do (setf (mtref E (+ iter i) (+ rank i)) (if (eq cn n) 1d0 0d0)))
                (incf rank nd)))
             (incf iter nd)))))
      E)))


(defun compute-internal-matrix (sim elem)
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
                      do (setf (mtref E (+ iter i) (+ rank i)) (if (eq cn n) 1d0 0d0)))
                (incf rank nd))) (incf iter nd)))))
      (setf (agg-internal-selector elem) (magicl:transpose E)))))


(defun update-aggregate-elements (sim)
  (when (sim-enable-aggregate sim)
    (setf
     (sim-agg-elems sim)
     (locate-aggregate-node-elements sim))
    ;; (set-aggregate-nodes sim)
    (iterate-over-agg-elem
     (sim-agg-elems sim)
     (lambda (elem)
       (iterate-over-cell-shape-local
        (cl-mpm:sim-mesh sim)
        (agg-interior-cell elem)
        (cl-mpm/mesh::cell-centroid (agg-interior-cell elem))
        (lambda (n w g) (setf (cl-mpm/mesh::node-interior n) t)))
       (compute-extension-matrix sim elem)
       (compute-inc-extension-matrix sim elem)
       (compute-internal-matrix sim elem)
       (setf (agg-mass-matrix elem) (assemble-mass sim elem)
             (agg-init-mass-matrix elem) (assemble-full-mass sim elem))
       (let* ((E (agg-shape-functions elem))
              (m-i (agg-init-mass-matrix elem))
              (m (agg-mass-matrix elem)))
         (setf (agg-vel-projection-matrix elem)
               (magicl:@
                E
                (magicl:inv m)
                (magicl:@ (magicl:transpose E) m-i))))))
    (let ((E (cl-mpm/aggregate::assemble-global-e sim)))
      (setf (sim-global-e sim) E)
      (setf (sim-global-ma sim) (magicl:@ (magicl:transpose E) (cl-mpm/aggregate::assemble-global-mass sim) E)))
    ))


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
         ;; (contrib (cl-mpm/mesh::))
         )
    (let ((iter 0))
      (declare (fixnum iter nd))
      (iterate-over-agg-elem-nodes
       sim
       elem
       (lambda (n)
         (with-accessors ((n-acc cl-mpm/mesh::node-acceleration)
                          (lock cl-mpm/mesh::node-lock))
             n
           (when t;(cl-mpm/mesh::node-agg n)
             (sb-thread:with-mutex (lock)
               (loop for i from 0 below nd
                     do (incf (varef n-acc i) (mtref acc (+ iter i) 0)))))
           (incf iter nd)))))))

(defun project-vel (sim elem acc)
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim))))
    (let ((iter 0))
      (iterate-over-agg-elem-nodes
       sim
       elem
       (lambda (n)
         (with-accessors ((n-acc cl-mpm/mesh::node-velocity)
                          (lock cl-mpm/mesh::node-lock))
             n
           (when (cl-mpm/mesh::node-agg n)
             (sb-thread:with-mutex (lock)
               (loop for i from 0 below nd
                     do (setf (varef n-acc i) (mtref acc (+ iter i) 0)))))
           (incf iter nd)))))))

(defun partial-project-force (sim elem damping)
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (E (agg-shape-functions elem))
         (m (agg-mass-matrix elem))
         (expander )
         (f (assemble-vector sim elem #'cl-mpm/mesh:node-force)))
    (project-vector-internal
     sim
     elem
     cl-mpm/mesh::node-force
     (magicl:@ (magicl:transpose E) f)
     )

    ;; (project-vector-partial
    ;;  sim
    ;;  elem
    ;;  cl-mpm/mesh::node-force
    ;;  (magicl:@ (magicl:transpose E) f)
    ;;  nil)
    ))

(defun partial-project-acc (sim elem damping)
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (E (agg-shape-functions elem))
         (proj (agg-internal-selector elem))
         (m (agg-mass-matrix elem))
         (f (assemble-vector sim elem #'cl-mpm/mesh::node-force))
         ;; (f-int (assemble-vector sim elem #'cl-mpm/mesh::node-internal-force))
         ;; (f-ext (assemble-vector sim elem #'cl-mpm/mesh::node-external-force))
         )
    (let ((acc (magicl:@ (magicl:inv m) proj f)))
      (project-acc sim elem acc)

      )))

(defun calculate-kinematic-global-agg (sim)
  (with-accessors ((bcs cl-mpm:sim-bcs)
                   (dt cl-mpm:sim-dt)
                   (mesh cl-mpm:sim-mesh))
      sim
    (let* ((E (cl-mpm/aggregate::assemble-global-e sim))
           (mii (cl-mpm/aggregate::assemble-global-mass sim))
           (vi (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-velocity)))
      (let* ((va
               (magicl:@
                E
                (magicl:inv (magicl:@ (magicl:transpose E) mii E))
                (magicl:transpose E)
                vi)))
        (cl-mpm/aggregate::project-global-vec sim va cl-mpm/mesh::node-velocity)))))

(defun calculate-forces-global-agg (sim)
  (with-accessors ((bcs cl-mpm:sim-bcs)
                   (dt cl-mpm:sim-dt)
                   (mesh cl-mpm:sim-mesh))
      sim
    (let* ((E (cl-mpm/aggregate::assemble-global-e sim))
           (mii (cl-mpm/aggregate::assemble-global-mass sim))
           (f (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-force))
           (acc
             (magicl:@
              E
              (magicl:inv (magicl:@ (magicl:transpose E) mii E))
              (magicl:transpose E)
              f)))
      (cl-mpm/aggregate::project-global-vec sim acc cl-mpm/mesh::node-acceleration)

      (cl-mpm/aggregate::project-global-vec sim (magicl:@ E (magicl:transpose E) (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-internal-force)) cl-mpm/mesh::node-internal-force)
      (cl-mpm/aggregate::project-global-vec sim (magicl:@ E (magicl:transpose E) (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-external-force)) cl-mpm/mesh::node-external-force)
      (cl-mpm/aggregate::project-global-vec sim (magicl:@ E (magicl:transpose E) (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-residual)) cl-mpm/mesh::node-residual)
      )
    ;; (let* ((E (cl-mpm/aggregate::assemble-global-e sim))
    ;;        (mii (cl-mpm/aggregate::assemble-global-mass sim))
    ;;        (f (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-force))
    ;;        (fa (magicl:@ E (magicl:transpose E) f)))

    ;;   (cl-mpm/aggregate::project-global-vec sim fa cl-mpm/mesh::node-force)

    ;;   (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt)
    ;;   (iterate-over-nodes
    ;;    mesh
    ;;    (lambda (node)
    ;;      (with-accessors ((agg cl-mpm/mesh::node-agg)
    ;;                       (active node-active)
    ;;                       (force node-force))
    ;;          node
    ;;        (when (and active agg)
    ;;          (fast-zero force)))))

    ;;   (let* ((f (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-force))
    ;;          (acc
    ;;            (magicl:@
    ;;             E
    ;;             (magicl:inv (magicl:@ (magicl:transpose E) mii E))
    ;;             (magicl:transpose E)
    ;;             f)))
    ;;     (cl-mpm/aggregate::project-global-vec sim acc cl-mpm/mesh::node-acceleration))
    ;;   (cl-mpm/aggregate::project-global-vec sim (magicl:@ E (magicl:transpose E) (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-internal-force)) cl-mpm/mesh::node-internal-force)
    ;;   (cl-mpm/aggregate::project-global-vec sim (magicl:@ E (magicl:transpose E) (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-external-force)) cl-mpm/mesh::node-external-force)
    ;;   (cl-mpm/aggregate::project-global-vec sim (magicl:@ E (magicl:transpose E) (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-residual)) cl-mpm/mesh::node-residual)
    ;;   )
    )
  )

(defun calculate-forces-agg-elem (sim elem damping)
  (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (E (agg-shape-functions elem))
         (Einc (agg-inc-shape-functions elem))
         (m (agg-mass-matrix elem))
         (f (assemble-vector sim elem #'cl-mpm/mesh:node-force)))
    (let ((acc (magicl:@ E (magicl:inv m) (magicl:transpose E) f)))
      (project-acc sim elem acc)
      ;; (reproject-inc-vector sim elem cl-mpm/mesh::node-internal-force)
      ;; (reproject-inc-vector sim elem cl-mpm/mesh::node-external-force)
      ;; (reproject-inc-vector sim elem cl-mpm/mesh::node-residual)
      )))

(defun calculate-kinematics-agg-elem (sim elem)
  (let* ((vel-project (agg-vel-projection-matrix elem))
         (v-i (assemble-vector sim elem #'cl-mpm/mesh::node-velocity)))
    (project-vel sim elem
                 (magicl:@
                  vel-project
                  v-i))))

(declaim (notinline update-node-kinematics))
(defmethod update-node-kinematics ((sim mpm-sim-aggregated))
  ;;For non-aggregate nodes, use simple mass matrix inversion
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (agg-elems sim-agg-elems)
                   (enable-aggregate sim-enable-aggregate))
      sim
    (cl-mpm:iterate-over-nodes
     mesh
     (lambda (node)
       (unless (cl-mpm/mesh::node-agg node)
         (cl-mpm::calculate-kinematics node))))
    ;;For each aggregated element set solve mass matrix and velocity
    (when enable-aggregate
      ;; (cl-mpm/aggregate::update-aggregate-elements sim)
      (calculate-kinematic-global-agg sim)
      ;; (iterate-over-agg-elem
      ;;  agg-elems
      ;;  (lambda (elem)
      ;;    (calculate-kinematics-agg-elem sim elem)))
      )))

(defun reproject-velocity (sim elem)
  (let* ((vel-project (agg-vel-projection-matrix elem))
         (v-i (assemble-vector sim elem #'cl-mpm/mesh::node-velocity)))
    (project-vel sim elem
                 (magicl:@
                  vel-project
                  v-i))))

(defun update-node-forces-agg (sim dt)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mass-scale sim-mass-scale)
                   (damping sim-damping-factor)
                   (damping-algo sim-damping-algorithm)
                   (agg-elems sim-agg-elems)
                   (enable-aggregate sim-enable-aggregate))
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
      (let* ((E (cl-mpm/aggregate::assemble-global-e sim))
             (mii (cl-mpm/aggregate::assemble-global-mass sim))
             (f (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-force))
             )
        (let* ((acc
                 (magicl:@
                  E
                  (magicl:inv (magicl:@ (magicl:transpose E) mii E))
                  (magicl:transpose E)
                  f)))
          (project-global-vec sim acc cl-mpm/mesh::node-acceleration))

        (project-global-vec sim (magicl:@ E (magicl:transpose E) (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-internal-force)) cl-mpm/mesh::node-internal-force)
        (project-global-vec sim (magicl:@ E (magicl:transpose E) (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-external-force)) cl-mpm/mesh::node-external-force)
        (project-global-vec sim (magicl:@ E (magicl:transpose E) (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-residual)) cl-mpm/mesh::node-residual)
        )
      ;; (calculate-forces-global-agg sim)
      ;; (iterate-over-agg-elem
      ;;  agg-elems
      ;;  (lambda (elem)
      ;;    (calculate-forces-agg-elem sim elem damping)))
      )
    (iterate-over-nodes
     mesh
     (lambda (node)
       (when (and (cl-mpm/mesh:node-active node)
                  (cl-mpm/mesh::node-agg node))
         (with-accessors ((mass node-mass)
                          (vel node-velocity)
                          (force node-force)
                          (acc node-acceleration))
             node
           (cl-mpm::integrate-vel-euler vel acc mass mass-scale dt damping)))))))

(defmethod cl-mpm::update-node-forces ((sim mpm-sim-aggregated))
  ;;For non-aggregate nodes, use simple mass matrix inversion
  (update-node-forces-agg sim (cl-mpm:sim-dt sim))
  )



(defmethod cl-mpm::filter-cells ((sim mpm-sim-aggregated))
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (dt cl-mpm::sim-dt)
                   (agg sim-enable-aggregate))
      sim
    (cl-mpm::iterate-over-cells
     mesh
     (lambda (cell)
       (cl-mpm::filter-cell mesh cell dt)))
    (when agg
      (cl-mpm/aggregate::update-aggregate-elements sim))
    ))

(defmethod cl-mpm::update-cells ((sim mpm-sim-aggregated))
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (dt cl-mpm::sim-dt)
                   (agg sim-enable-aggregate))
      sim
    (cl-mpm::iterate-over-cells
     mesh
     (lambda (cell)
       ;; (cl-mpm::filter-cell mesh cell dt)
       (cl-mpm::update-cell mesh cell dt)))
    ;; (when agg
    ;;   (cl-mpm/aggregate::update-aggregate-elements sim))
    ))


(defun reproject-displacements (sim elem)
  (let* ((E (agg-shape-functions elem))
         (proj (agg-internal-selector elem))
         (d-i (assemble-vector sim elem #'cl-mpm/mesh::node-displacment)))
    (project-vector
     sim
     elem
     cl-mpm/mesh::node-displacment
     (magicl:@
      E
      proj
      d-i))))

(defmethod cl-mpm::update-nodes (sim)
  (with-accessors ((mesh sim-mesh)
                   (dt sim-dt)
                   (agg sim-enable-aggregate))
      sim
    (iterate-over-nodes
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (cl-mpm::update-node node dt))))
    ;; (when agg
    ;;   (iterate-over-agg-elem
    ;;    (sim-agg-elems sim)
    ;;    (lambda (elem)
    ;;      (reproject-displacements sim elem))))
    ))


(defun filter-nodes (sim filter)
  (let ((nodes (cl-mpm/mesh::mesh-nodes (cl-mpm:sim-mesh sim))))
    (remove-if-not filter (make-array (array-total-size nodes) :displaced-to nodes))))


(defun assemble-global-fdofs (sim)
  (let* ((active-nodes (filter-nodes sim #'cl-mpm/mesh::node-active))
         (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (ndofs (* (length active-nodes) nd))
         (fd (make-array ndofs))
         )
    (iterate-over-agg-elem
     (sim-agg-elems sim)
     (lambda (elem)
       (let ((ma (agg-mass-matrix elem))
             (index-map nil)
             (iter 0)
             (node-count (length (agg-node-list elem))))
         (iterate-over-agg-elem-nodes
          sim
          elem
          (lambda (n)
            (let ((index (* nd (position n active-nodes))))
              (push index index-map))))
         (setf index-map (make-array (length index-map) :initial-contents (nreverse index-map)))
         )))))

(defun assemble-global-e (sim)
  (let* ((agg-nodes (filter-nodes sim #'cl-mpm/mesh::node-interior))
         (active-nodes (filter-nodes sim (lambda (n) (or (cl-mpm/mesh::node-agg n)
                                                         (cl-mpm/mesh::node-interior n)))))
         (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (ndof (* (length active-nodes) nd))
         (ndofC (* (length agg-nodes) nd))
         (mesh (cl-mpm:sim-mesh sim))
         (e (cl-mpm/utils::arb-matrix ndof ndofC)))
    (iterate-over-agg-elem
     (sim-agg-elems sim)
     (lambda (elem)
       (let* ((eshape (agg-shape-functions elem))
              (index-map nil)
              (int-map nil)
              (iter 0)
              (nc (length (agg-node-list elem)))
              (ni (expt 2 nd)))
         (iterate-over-agg-elem-nodes
          sim
          elem
          (lambda (n)
            (let ((index (* nd (position n active-nodes))))
              (push index index-map))
            (when (cl-mpm/mesh::node-interior n)
              (let ((index (* nd (position n agg-nodes))))
                (push index int-map)))
            ))
         ;; (iterate-over-internal-cell
         ;;  mesh
         ;;  elem
         ;;  (lambda (n)
         ;;    (let ((index (* nd (position n agg-nodes))))
         ;;      (push index int-map))))

         (setf index-map (make-array (length index-map) :initial-contents (nreverse index-map)))
         (setf int-map (make-array (length int-map) :initial-contents (nreverse int-map)))
         ;; (setf index-map (make-array (length index-map) :initial-contents index-map))
         ;; (setf int-map (make-array (length int-map) :initial-contents int-map))
         ;; (break)

         ;; (break)
         ;; (loop for x from 0 below ni
         ;;       do
         ;;          (loop for y from 0 below ni
         ;;                do
         ;;                   (loop for i from 0 below nd
         ;;                         do (setf (mtref e (+ (aref index-map x) i) (+ (aref int-map y) i))
         ;;                                  (mtref eshape (+ x i) (+ y i))))))
         (loop for x from 0 below nc
              do
                 (loop for y from 0 below ni
                       do
                          (loop for i from 0 below nd
                                do (setf (mtref e (+ (aref index-map x) i) (+ (aref int-map y) i))
                                         (mtref eshape (+ (* nd x) i) (+ (* nd y) i))))))
         )))
    (values e)))

(defun assemble-global-vec (sim accessor)
  (declare (function accessor))
  (let* (;; (agg-nodes (filter-nodes sim #'cl-mpm/mesh::node-interior))
         (active-nodes (filter-nodes sim (lambda (n) (or (cl-mpm/mesh::node-agg n)
                                                         (cl-mpm/mesh::node-interior n)))))
         (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (ndof (* (length active-nodes) nd))
         (v (cl-mpm/utils::arb-matrix ndof 1)))
    (iterate-over-agg-elem
     (sim-agg-elems sim)
     (lambda (elem)
       (iterate-over-agg-elem-nodes
        sim
        elem
        (lambda (n)
          (let ((index (* nd (position n active-nodes))))
            (let ((vec (funcall accessor n)))
              (loop for i from 0 below nd
                    do (setf (mtref v (+ index i) 0)
                             (varef vec i)))))))))
    (values v)))

(defun assemble-global-internal-vec (sim accessor)
  (declare (function accessor))
  (let* (;; (agg-nodes (filter-nodes sim #'cl-mpm/mesh::node-interior))
         (active-nodes (filter-nodes sim (lambda (n) (or ;(cl-mpm/mesh::node-agg n)
                                                         (cl-mpm/mesh::node-interior n)))))
         (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (ndof (* (length active-nodes) nd))
         (v (cl-mpm/utils::arb-matrix ndof 1)))
    (iterate-over-agg-elem
     (sim-agg-elems sim)
     (lambda (elem)
       (iterate-over-agg-elem-nodes
        sim
        elem
        (lambda (n)
          (when (cl-mpm/mesh::node-interior n)
            (let ((index (* nd (position n active-nodes))))
              (let ((vec (funcall accessor n)))
                (loop for i from 0 below nd
                      do (setf (mtref v (+ index i) 0)
                               (varef vec i))))))))))
    (values v)))


(defun assemble-global-mass (sim)
  (let* ((agg-nodes (filter-nodes sim #'cl-mpm/mesh::node-interior))
         (active-nodes (filter-nodes sim (lambda (n) (or (cl-mpm/mesh::node-agg n)
                                                         (cl-mpm/mesh::node-interior n)))))
         (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (ndof (* (length active-nodes) nd))
         (m (cl-mpm/utils::arb-matrix ndof ndof)))
    (iterate-over-agg-elem
     (sim-agg-elems sim)
     (lambda (elem)
       (let* ((eshape (agg-shape-functions elem))
              (index-map nil)
              (int-map nil)
              (iter 0)
              (nc (length (agg-node-list elem)))
              (ni (expt 2 nd)))
         (iterate-over-agg-elem-nodes
          sim
          elem
          (lambda (n)
            (let ((index (* nd (position n active-nodes))))
              (let ((n-mass (cl-mpm/mesh::node-mass n)))
                (loop for x from 0 below nc
                      do
                         (loop for i from 0 below nd
                               do (setf (mtref m (+ index i) (+ index i))
                                        n-mass))))))))))
    (values m)))

(defun assemble-global-full-mass (sim)
  (let* (;; (active-nodes (filter-nodes sim #'cl-mpm/mesh::node-active))
         (active-nodes (filter-nodes sim (lambda (n) (or (cl-mpm/mesh::node-agg n)
                                                         (cl-mpm/mesh::node-interior n)))))
         (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (ndofs (* (length active-nodes) nd))
         (m (cl-mpm/utils::arb-matrix ndofs ndofs)))
    (iterate-over-agg-elem
     (sim-agg-elems sim)
     (lambda (elem)
       (let ((ma (agg-mass-matrix elem))
             (index-map nil)
             (node-count (expt 2 nd)))
         (iterate-over-agg-elem-nodes
          sim
          elem
          (lambda (n)
            (let ((index (* nd (position n active-nodes))))
              (push index index-map))))
         (setf index-map (make-array (length index-map) :initial-contents (nreverse index-map)))
         (loop for x from 0 below node-count
              do
                 (loop for y from 0 below node-count
                       do
                          (loop for i from 0 below nd
                                do (setf (mtref m (+ (aref index-map x) i) (+ (aref index-map y) i))
                                         (mtref ma (+ x i) (+ y i)))))))))
    (values m)))

;; (defun assemble-global-vector (sim accessor)
;;   (declare (function accessor))
;;   (let* ((active-nodes (filter-nodes sim #'cl-mpm/mesh::node-active))
;;          (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
;;          (ndofs (* (length active-nodes) nd))
;;          (v (cl-mpm/utils::arb-matrix ndofs 1))
;;          )
;;     (iterate-over-agg-elem
;;      (sim-agg-elems sim)
;;      (lambda (elem)
;;        (let ((ma (agg-mass-matrix elem))
;;              (index-map nil)
;;              (iter 0)
;;              (node-count (length (agg-node-list elem))))
;;          (iterate-over-agg-elem-nodes
;;           sim
;;           elem
;;           (lambda (n)
;;             (let ((index (* nd (position n active-nodes))))
;;               (push index index-map))))
;;          (setf index-map (make-array (length index-map) :initial-contents (nreverse index-map)))
;;          (loop for x from 0 below node-count
;;               do
;;                  (loop for y from 0 below node-count
;;                        do
;;                           (loop for i from 0 below nd
;;                                 do (let ((index (aref index-map x)))
;;                                      (setf (mtref m (+ (aref index-map x) i) (+ (aref index-map y) i))
;;                                            (mtref ma (+ x i) (+ y i))))))))))
;;     ))



(defun list< (a b)
  (cond ((null a) (not (null b)))
        ((null b) nil)
        ((= (first a) (first b)) (list< (rest a) (rest b)))
        (t (< (first a) (first b)))))
(defun print-agg (sim)
  (let ((elems (sort (copy-seq (cl-mpm/aggregate::sim-agg-elems sim))
                     #'list<
                     :key (lambda (e) (cl-mpm/mesh::cell-index (agg-interior-cell e))))))
    (loop for elem across elems;(cl-mpm/aggregate::sim-agg-elems sim)
          do
             (let* (;(elem (aref (cl-mpm/aggregate::sim-agg-elems sim) 0))
                    (map-size 4)
                    (map (make-array (list map-size map-size) :initial-element "-"))
                    (center-index (cl-mpm/mesh::cell-index (cl-mpm/aggregate::agg-interior-cell elem))))
               (cl-mpm/aggregate::iterate-over-agg-elem-nodes
                sim
                elem
                (lambda (n)
                  (let ((ni (cl-mpm/mesh::node-index n)))
                    (let ((index-diff (mapcar #'+ (list 1 1) (mapcar #'- ni center-index))))
                      (setf (aref map (first index-diff) (second index-diff)) (if
                                                                               (cl-mpm/mesh::node-agg n)
                                                                               "a"
                                                                               "x")))
                    (format t "~a ~A~%" ni (cl-mpm/mesh::node-agg n))
                    )))
               (format t "Cell ~A~%" (cl-mpm/mesh::cell-index (agg-interior-cell elem)))
               (loop for x from (- map-size 1) downto 0
                     do (progn
                          (loop for y from 0 below map-size
                                do (format t "~A" (aref map y x)))
                          (format t "~%")))
               (let ((E (cl-mpm/aggregate::agg-shape-functions elem))
                     ;; (m (cl-mpm/aggregate::agg-mass-matrix elem))
                     ;; (proj (cl-mpm/aggregate::agg-internal-selector elem))
                     ;; (mii (cl-mpm/aggregate::assemble-full-mass sim elem))
                     ;; (f (cl-mpm/aggregate::assemble-vector sim elem #'cl-mpm/mesh::node-force))
                     ;; (id (cl-mpm/aggregate::compute-identity sim (aref (cl-mpm/aggregate::sim-agg-elems sim) 0)))
                     ;; (expander (cl-mpm/aggregate::compute-internal-expander sim elem))
                     )
                 ;; (pprint (magicl:transpose (cl-mpm/aggregate::agg-inc-shape-functions elem)))
                 (format t "~A~%" (cl-mpm/aggregate::agg-contribution-list elem))
                 )))))
