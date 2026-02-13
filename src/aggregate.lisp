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

(defclass mpm-sim-agg-usf (cl-mpm::mpm-sim-usf cl-mpm/aggregate::mpm-sim-aggregated)
  ()
  (:documentation "Explicit simulation with update stress first update"))



(defclass mpm-sim-aggregated (mpm-sim)
  ((agg-elems
    :initform (make-array 0)
    :accessor sim-agg-elems)

   (enable-aggregate
    :initform t
    :initarg :enable-aggregate
    :accessor sim-enable-aggregate)
   (agg-nodes-fd
    :initform nil
    :accessor sim-agg-nodes-fd)
   (agg-nodes-fdc
    :initform nil
    :accessor sim-agg-nodes-fdc)
   (global-e
    :initform nil
    :accessor sim-global-e)
   (global-ma
    :initform nil
    :accessor sim-global-ma)
   (internal-v
    :initform nil
    :accessor sim-internal-v)
   ))

(defun iterate-over-cell-patch-2d (sim node ring-size func)
  (declare (function func))
  (with-accessors ((mesh cl-mpm::sim-mesh))
      sim
    (with-accessors ((index cl-mpm/mesh::node-index))
        node
      (loop for dx from (- (+ ring-size 1)) to ring-size
            do
               (loop for dy from (- (+ ring-size 1)) to ring-size
                     do
                        (let ((cell-index (mapcar #'+ index (list dx dy 0))))
                          (when (cl-mpm/mesh::in-bounds-cell mesh cell-index)
                            (funcall func (cl-mpm/mesh::get-cell mesh cell-index)))))))))

(defun iterate-over-cell-patch-3d (sim node ring-size func)
  (declare (function func))
  (with-accessors ((mesh cl-mpm::sim-mesh))
      sim
    (with-accessors ((index cl-mpm/mesh::node-index))
        node
      (loop for dx from (- ring-size 1) to ring-size
            do
               (loop for dy from (- ring-size 1) to ring-size
                     do
                        (loop for dz from (- ring-size 1) to ring-size
                              do
                                 (let ((cell-index (mapcar #'+ index (list dx dy dz))))
                                   (when (cl-mpm/mesh::in-bounds-cell mesh cell-index)
                                     (funcall func (cl-mpm/mesh::get-cell mesh cell-index))))))))))

(defun iterate-over-cell-patch (sim node ring-size func)
  (declare (function func))
  (if (= (cl-mpm/mesh::mesh-nd (cl-mpm:sim-mesh sim)) 2)
      (iterate-over-cell-patch-2d sim node ring-size func)
      (iterate-over-cell-patch-3d sim node ring-size func)))

(defun get-closest-cell (sim node)
  (let ((closest-elem nil)
        (dist 0d0)
        ;; (mutex (sb-thread:make-mutex))
        (pos (cl-mpm/mesh::node-position node)))
    (flet ((check-cell (cell)
             (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                              (index cl-mpm/mesh::cell-index)
                              (centroid cl-mpm/mesh::cell-centroid)
                              (active cl-mpm/mesh::cell-active)
                              (agg cl-mpm/mesh::cell-agg))
                 cell
               (when (and
                      (cl-mpm/mesh::cell-active cell)
                      (not (cl-mpm/mesh::cell-partial cell))
                      (not (cl-mpm/mesh::cell-agg cell)))
                 (let ((dist-tr (cl-mpm/fastmaths::diff-norm
                                 pos
                                 centroid)))
                   (when (or
                          (not closest-elem)
                          (> dist dist-tr))
                     (when (or
                            (not closest-elem)
                            (> dist dist-tr))
                       (setf dist dist-tr
                             closest-elem cell))))))))
      (iterate-over-cell-patch
       sim
       node
       1
       #'check-cell)
      (unless closest-elem
        (iterate-over-cell-patch
         sim
         node
         2
         #'check-cell)
        (cl-mpm::iterate-over-cells-serial
         (cl-mpm:sim-mesh sim)
         (lambda (cell)
           (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                            (index cl-mpm/mesh::cell-index)
                            (centroid cl-mpm/mesh::cell-centroid)
                            (active cl-mpm/mesh::cell-active)
                            (agg cl-mpm/mesh::cell-agg))
               cell
             (when (and
                    (cl-mpm/mesh::cell-active cell)
                    (not (cl-mpm/mesh::cell-partial cell))
                    (not (cl-mpm/mesh::cell-agg cell)))
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
             )))))
    closest-elem))


(defun weight-local-pos (x y z nd dx dy dz)
  (let ((dx (- (* dx 2) 1))
        (dy (- (* dy 2) 1))
        (dz (- (* dz 2) 1)))
    (case nd
      (2 (* 0.25d0  (+ 1d0 (* dx x)) (+ 1d0 (* dy y))))
      (3 (* 0.125d0 (+ 1d0 (* dx x)) (+ 1d0 (* dy y)) (+ 1d0 (* dz z)))))))



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

(defun filter-nodes (sim filter)
  (let ((nodes (cl-mpm/mesh::mesh-nodes (cl-mpm:sim-mesh sim))))
    (remove-if-not filter (make-array (array-total-size nodes) :displaced-to nodes))))



(defun locate-aggregate-nodes (sim)
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
                        (neighbours cl-mpm/mesh::cell-neighbours)
                        (active cl-mpm/mesh::cell-active))
           cell
         (when (and active partial)
           (setf agg t)
           ;;Set all our neighbours to also be aggregate
           (loop for n in neighbours
                 do (setf (cl-mpm/mesh::cell-agg n) t))))))

    ;;Next set any nodes on agg elements to be agg
    (cl-mpm::iterate-over-cells
     mesh
     (lambda (cell)
       (with-accessors ((nodes cl-mpm/mesh::cell-nodes)
                        (agg cl-mpm/mesh::cell-agg)
                        (active cl-mpm/mesh::cell-active))
           cell
         ;;Set all aggregated nodes
         (when (and active agg)
           (loop for n in nodes
                 do (when (cl-mpm/mesh::node-active n)
                      (setf (cl-mpm/mesh::node-agg n) t)))))))

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
                 (progn
                   (setf (cl-mpm/mesh::node-agg-interior-cell node) closest-elem)
                   (setf (cl-mpm/mesh::cell-interior closest-elem) t))
                 (error "No closest elem?")))))))

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
           (loop for n in nodes
                 do (when (cl-mpm/mesh::node-active n)
                      (setf (cl-mpm/mesh::node-agg n) t
                            (cl-mpm/mesh::node-interior n) t)))))))
    (setf
     (cl-mpm/aggregate::sim-agg-nodes-fdc sim)
     (filter-nodes sim #'cl-mpm/mesh::node-interior)

     (cl-mpm/aggregate::sim-agg-nodes-fd sim)
     (filter-nodes sim (lambda (n) (or (cl-mpm/mesh::node-agg n)
                                       (cl-mpm/mesh::node-interior n)))))
    (let ((fdc 0))
      (loop for n across (cl-mpm/aggregate::sim-agg-nodes-fdc sim)
            do (progn
                 (setf (cl-mpm/mesh::node-agg-fdc n) fdc)
                 (incf fdc))))
    (let ((fd 0))
      (loop for n across (cl-mpm/aggregate::sim-agg-nodes-fd sim)
            do (progn
                 (setf (cl-mpm/mesh::node-agg-fd n) fd)
                 (incf fd))))
    )
  )

(defun assemble-global-scalar (sim accessor)
  (declare (function accessor))
  (let* ((active-nodes (sim-agg-nodes-fd sim))
         (ndof (length active-nodes))
         (v (cl-mpm/utils::arb-matrix ndof 1)))
    (cl-mpm::iterate-over-nodes-array
     active-nodes
     (lambda (n)
       (setf (varef v (cl-mpm/mesh::node-agg-fd n)) (funcall accessor n))))
    (values v)))

(defun assemble-global-vec (sim accessor dim)
  (declare (function accessor))
  (let* ((active-nodes (sim-agg-nodes-fd sim))
         (ndof (length active-nodes))
         (v (cl-mpm/utils::arb-matrix ndof 1)))
    (cl-mpm::iterate-over-nodes-array
     active-nodes
     (lambda (n)
       (setf (varef v (cl-mpm/mesh::node-agg-fd n)) (varef (funcall accessor n) dim))))
    (values v)))

(defun assemble-internal-vec (sim accessor dim)
  (declare (function accessor))
  (let* ((int-nodes (sim-agg-nodes-fdc sim))
         (ndof (length int-nodes))
         (v (cl-mpm/utils::arb-matrix ndof 1)))
    (cl-mpm::iterate-over-nodes-array
     int-nodes
     (lambda (n)
       (setf (varef v (cl-mpm/mesh::node-agg-fdc n)) (varef (funcall accessor n) dim))))
    (values v)))

(defmacro zero-global (sim accessor dim)
  `(progn
     (let* ((active-nodes (sim-agg-nodes-fd ,sim)))
       (cl-mpm::iterate-over-nodes-array
        active-nodes
        (lambda (n)
          (setf (varef (funcall ,accessor n) ,dim) 0d0))))))

(defmacro zero-global-scalar (sim accessor)
  `(progn
     (let* ((active-nodes (sim-agg-nodes-fd ,sim)))
       (cl-mpm::iterate-over-nodes-array
        active-nodes
        (lambda (n)
          (setf (,accessor n) 0d0))))))

(defmacro increment-global-vec (sim vector accessor dim)
  `(progn
     (let* ((active-nodes (sim-agg-nodes-fd ,sim))
            (proj-val ,vector))
       (cl-mpm::iterate-over-nodes-array
        active-nodes
        (lambda (n)
          (incf (varef (funcall ,accessor n) ,dim)
                (varef proj-val (cl-mpm/mesh::node-agg-fd n))))))))

(defmacro increment-global-scalar (sim vector accessor)
  `(progn
     (let* ((active-nodes (sim-agg-nodes-fd ,sim))
            (proj-val ,vector))
       (cl-mpm::iterate-over-nodes-array
        active-nodes
        (lambda (n)
          (incf (,accessor n)
                (varef proj-val (cl-mpm/mesh::node-agg-fd n))))))))

(defmacro increment-internal-scalar (sim vector accessor)
  `(progn
     (let* ((active-nodes (sim-agg-nodes-fdc ,sim))
            (proj-val ,vector))
       (cl-mpm::iterate-over-nodes-array
        active-nodes
        (lambda (n)
          (incf (,accessor n)
                (varef proj-val (cl-mpm/mesh::node-agg-fdc n))))))))

(defmacro project-global-vec (sim vector accessor dim)
  `(progn
     (let* ((active-nodes (sim-agg-nodes-fd ,sim))
            (proj-val ,vector))
       (cl-mpm::iterate-over-nodes-array
        active-nodes
        (lambda (n)
          (setf (varef (funcall ,accessor n) ,dim)
                (varef proj-val (cl-mpm/mesh::node-agg-fd n))))))))

(defmacro project-int-vec (sim vector accessor dim)
  `(progn
     (let* ((int-nodes (sim-agg-nodes-fdc ,sim))
            (proj-val ,vector))
       (cl-mpm::iterate-over-nodes-array
        int-nodes
        (lambda (n)
          (setf (varef (funcall ,accessor n) ,dim)
                (varef proj-val (cl-mpm/mesh::node-agg-fdc n))))))))
;;New version
(defun assemble-e (sim)
  (let* (;; (agg-nodes (filter-nodes sim #'cl-mpm/mesh::node-interior))
         ;; (active-nodes (filter-nodes sim (lambda (n) (or (cl-mpm/mesh::node-agg n)
         ;;                                                 (cl-mpm/mesh::node-interior n)))))
         (agg-nodes (sim-agg-nodes-fdc sim))
         (active-nodes (sim-agg-nodes-fd sim))
         (mesh (cl-mpm:sim-mesh sim))
         (ndof (length active-nodes))
         (ndofC (length agg-nodes))
         (e (cl-mpm/utils::arb-matrix ndof ndofC)))
    (cl-mpm::iterate-over-nodes-serial
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (when (cl-mpm/mesh::node-agg node)
           (if (cl-mpm/mesh::node-interior node)
               (progn
                 ;;Interior -> populate 1s
                 (setf (mtref e (cl-mpm/mesh::node-agg-fd node) (cl-mpm/mesh::node-agg-fdc node)) 1d0))
               (progn
                 ;;Aggregate -> populate svp
                 (iterate-over-cell-shape-local
                  mesh
                  (cl-mpm/mesh::node-agg-interior-cell node)
                  (cl-mpm/mesh::node-position node)
                  (lambda (cn weight grads)
                    (incf (mtref e (cl-mpm/mesh::node-agg-fd node) (cl-mpm/mesh::node-agg-fdc cn))
                          weight)))))))))
    (values e)))

(defun assemble-global-mass-matrix (sim)
  (let* ((active-nodes (sim-agg-nodes-fd sim))
         (ndofs (length active-nodes))
         (mesh (cl-mpm:sim-mesh sim))
         (m (cl-mpm/utils::arb-matrix ndofs ndofs))
         )
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (when (cl-mpm/mesh::node-agg node)
           (setf (mtref m
                        (cl-mpm/mesh::node-agg-fd node)
                        (cl-mpm/mesh::node-agg-fd node))
                 (cl-mpm/mesh::node-mass node))))))
    (values m)))

(defun update-aggregate-elements (sim)
  (when (sim-enable-aggregate sim)
    (locate-aggregate-nodes sim)
    (let ((E (cl-mpm/aggregate::assemble-e sim)))
      (setf (sim-global-e sim) E)
      (setf (sim-global-ma sim) (magicl:@ (magicl:transpose E) (cl-mpm/aggregate::assemble-global-mass-matrix sim) E)))))

(defun update-mass-matrix (sim)
  "Assuming that we already have an extension matrix, attempt to update the mass matrix directly"
  (when (sim-enable-aggregate sim)
    (let ((E (cl-mpm/aggregate::sim-global-e sim)))
      (setf (sim-global-ma sim) (magicl:@
                                 (magicl:transpose E)
                                 (cl-mpm/aggregate::assemble-global-mass-matrix sim)
                                 E)))))


(defun assemble-internal-bcs (sim d)
  (let* ((int-nodes (cl-mpm/aggregate::sim-agg-nodes-fdc sim))
         (ndof (length int-nodes))
         (v (cl-mpm/utils::arb-matrix ndof 1)))
    (cl-mpm::iterate-over-nodes-array
     int-nodes
     (lambda (n)
       (let ((index (cl-mpm/mesh::node-agg-fdc n)))
         (let ((bcs (cl-mpm/mesh::node-bcs n)))
           (if bcs
             (setf (varef v index) (varef bcs d))
             (setf (varef v index) 1d0))))))
    v))

(declaim (notinline linear-solve-with-bcs))
(defun linear-solve-with-bcs (ma v bcs &optional (target-vi nil))
  (let ((target-vi (if target-vi
                       target-vi
                       (cl-mpm/utils::arb-matrix (magicl:nrows v) 1)
                       ;; (cl-mpm/utils::deep-copy v)
                       ))
        (bc-map (make-array (magicl:nrows bcs) :fill-pointer 0)))
    (cl-mpm/fastmaths:fast-zero target-vi)
    (loop for i from 0
          for bc across (magicl::storage bcs)
          do (when (> bc 0d0)
               (vector-push-extend i bc-map)))
    ;; (print (magicl:nrows bcs))
    ;;TODO have a fallback if no bcs are applied (i.e. no resizing required)
    (let ((reduced-size (length bc-map)))
      ;;When we are solving a fully fixed system - i.e. out of plane dimensions
      (when (> reduced-size 0)
        (let* ((A-r (cl-mpm/utils::arb-matrix reduced-size reduced-size))
               (v-r (cl-mpm/utils::arb-matrix reduced-size 1)))
          (lparallel:pdotimes (i reduced-size)
            (dotimes (j reduced-size)
              (setf (mtref a-r i j) (mtref ma (aref bc-map i) (aref bc-map j))))
            (setf (varef v-r i) (varef v (aref bc-map i))))
          (let ((vs (magicl::linear-solve A-r v-r)))
            (lparallel:pdotimes (i reduced-size)
              (setf (varef target-vi (aref bc-map i)) (varef vs i))))
          )))
    target-vi
    ))


(declaim (notinline linear-solve-with-bcs))
(defun @-with-bcs (ma v bcs &optional (target-vi nil))
  (let ((target-vi (if target-vi
                       target-vi
                       (cl-mpm/utils::arb-matrix (magicl:nrows v) 1)
                       ;; (cl-mpm/utils::deep-copy v)
                       ))
        (bc-map (make-array (magicl:nrows bcs) :fill-pointer 0)))
    (cl-mpm/fastmaths:fast-zero target-vi)
    (loop for i from 0
          for bc across (magicl::storage bcs)
          do (when (> bc 0d0)
               (vector-push-extend i bc-map)))
    ;; (print (magicl:nrows bcs))
    ;;TODO have a fallback if no bcs are applied (i.e. no resizing required)
    (let ((reduced-size (length bc-map)))
      ;;When we are solving a fully fixed system - i.e. out of plane dimensions
      (when (> reduced-size 0)
        (let* ((A-r (cl-mpm/utils::arb-matrix reduced-size reduced-size))
               (v-r (cl-mpm/utils::arb-matrix reduced-size 1)))
          (lparallel:pdotimes (i reduced-size)
            (dotimes (j reduced-size)
              (setf (mtref a-r i j) (mtref ma (aref bc-map i) (aref bc-map j))))
            (setf (varef v-r i) (varef v (aref bc-map i))))
          (let ((vs (magicl::@ A-r v-r)))
            (lparallel:pdotimes (i reduced-size)
              (setf (varef target-vi (aref bc-map i)) (varef vs i))))
          )))
    target-vi
    ))


(defun calculate-kinematics-agg (sim)
  (with-accessors ((mesh cl-mpm::sim-mesh))
      sim
    (let ((ma (cl-mpm/aggregate::sim-global-ma sim))
          (E (cl-mpm/aggregate::sim-global-e sim)))
      (iterate-over-dimensions
       (cl-mpm/mesh:mesh-nd mesh)
       (lambda (d)
         (let* ((vi (magicl:@
                     (magicl:transpose E)
                     (assemble-global-vec sim #'cl-mpm/mesh::node-velocity d))))
           (project-global-vec
            sim
            (magicl:@ E (linear-solve-with-bcs ma vi (assemble-internal-bcs sim d)))
            #'cl-mpm/mesh::node-velocity
            d)))))))

(declaim (notinline update-node-kinematics))
(defmethod cl-mpm::update-node-kinematics ((sim mpm-sim-aggregated))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (agg-elems sim-agg-elems)
                   (enable-aggregate sim-enable-aggregate))
      sim
    ;;For non-aggregate nodes, use simple mass matrix inversion
    (cl-mpm:iterate-over-nodes
     mesh
     (lambda (node)
       (unless (cl-mpm/mesh::node-agg node)
         (cl-mpm::calculate-kinematics node))))
    ;; When we aggregate, do the global aggregation task
    (when enable-aggregate
      (calculate-kinematics-agg sim)
      (cl-mpm::fast-scale! (cl-mpm/aggregate::sim-global-ma sim) (cl-mpm::sim-mass-scale sim)))))

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
      (cl-mpm/aggregate::update-aggregate-elements sim))))



(defmethod cl-mpm::update-node-forces ((sim mpm-sim-aggregated))
  (update-node-forces-agg sim (cl-mpm:sim-dt sim)))

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
       (when (and (cl-mpm/mesh:node-active node))
         (cl-mpm::calculate-forces node damping 0d0 mass-scale)
         ;; (if (cl-mpm/mesh::node-agg node)
         ;;     (cl-mpm::calculate-forces node damping 0d0 mass-scale)
         ;;     (cl-mpm::calculate-forces node damping dt mass-scale)
         ;;     )
         )))

    (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt)
    ;;For each aggregated element set solve mass matrix and velocity
    (when enable-aggregate
      (let* ((E (cl-mpm/aggregate::sim-global-e sim))
             (ma (cl-mpm/aggregate::sim-global-ma sim)))
        (iterate-over-dimensions
         (mesh-nd mesh)
         (lambda (d)
           (let ((fa
                   (magicl:@
                    (magicl:transpose E)
                    (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-force d))))
             (apply-internal-bcs sim fa d)
             (let* ((acc
                      (linear-solve-with-bcs ma fa (assemble-internal-bcs sim d))))
               (cl-mpm/aggregate::apply-internal-bcs sim acc d)
               (project-global-vec sim (magicl:@ E acc) #'cl-mpm/mesh::node-acceleration d))))))
      ;; (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt)
      )
    (iterate-over-nodes
     mesh
     (lambda (node)
       (when (and (cl-mpm/mesh:node-active node))
         (with-accessors ((mass node-mass)
                          (vel node-velocity)
                          (force node-force)
                          (internal cl-mpm/mesh::node-interior)
                          (acc node-acceleration))
             node
           (when t
             (cl-mpm::integrate-vel-euler vel acc mass mass-scale dt damping))))))
    ;; (when enable-aggregate
    ;;   ;; (pprint "Hello")
    ;;   (project-velocity sim))
    (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt)
    ))


(defun project-displacement (sim)
  (with-accessors ((sim-mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)
                   (dt cl-mpm:sim-dt))
      sim
    (iterate-over-dimensions
     (cl-mpm/mesh::mesh-nd mesh)
     (lambda (d)
       (let ((E (sim-global-e sim))
             (vel-proj (cl-mpm/aggregate::assemble-internal-vec sim #'cl-mpm/mesh::node-displacment d)))
         (apply-internal-bcs sim vel-proj d)
         (cl-mpm/aggregate::project-global-vec
          sim
          (magicl:@ E vel-proj)
          #'cl-mpm/mesh::node-displacment
          d))))
    (cl-mpm::apply-bcs mesh (cl-mpm:sim-bcs sim) dt)))

(defun project-velocity (sim)
  (with-accessors ((sim-mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)
                   (dt cl-mpm:sim-dt))
      sim
    (iterate-over-dimensions
     (cl-mpm/mesh::mesh-nd mesh)
     (lambda (d)
       (let ((E (sim-global-e sim))
             (vel-proj (cl-mpm/aggregate::assemble-internal-vec sim #'cl-mpm/mesh::node-velocity d)))
         (apply-internal-bcs sim vel-proj d)
         (cl-mpm/aggregate::project-global-vec
          sim
          (magicl:@ E vel-proj)
          #'cl-mpm/mesh::node-velocity
          d))))
    (cl-mpm::apply-bcs mesh (cl-mpm:sim-bcs sim) dt)))

(defun project-acceleration (sim)
  (with-accessors ((sim-mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)
                   (dt cl-mpm:sim-dt))
      sim
    (iterate-over-dimensions
     (cl-mpm/mesh::mesh-nd mesh)
     (lambda (d)
       (let ((E (sim-global-e sim))
             (vel-proj (cl-mpm/aggregate::assemble-internal-vec sim #'cl-mpm/mesh::node-acceleration d)))
         (apply-internal-bcs sim vel-proj d)
         (cl-mpm/aggregate::project-global-vec
          sim
          (magicl:@ E vel-proj)
          #'cl-mpm/mesh::node-acceleration
          d))))
    (cl-mpm::apply-bcs mesh (cl-mpm:sim-bcs sim) dt)))

(defmethod cl-mpm::update-nodes ((sim cl-mpm/aggregate::mpm-sim-aggregated))
  (with-accessors ((mesh sim-mesh)
                   (dt sim-dt)
                   (agg sim-enable-aggregate))
      sim
    (iterate-over-nodes
     mesh
     (lambda (node)
       (when (and (cl-mpm/mesh:node-active node)
                  ;; (or (not (cl-mpm/mesh::node-agg node))
                  ;;     (cl-mpm/mesh::node-interior node))
                  )
         (cl-mpm::update-node node dt))))
    (when agg
      (project-displacement sim))
    ))

(defun apply-internal-bcs (sim vec d)
  (let* ((int-nodes (cl-mpm/aggregate::sim-agg-nodes-fdc sim)))
    (when (not (= (magicl:nrows vec) (length (cl-mpm/aggregate::sim-agg-nodes-fdc sim))))
      (error "Incorrect BCs applied"))
    (cl-mpm::iterate-over-nodes-array
     int-nodes
     (lambda (n)
       (let ((index (cl-mpm/mesh::node-agg-fdc n)))
         (let ((bcs (cl-mpm/mesh::node-bcs n)))
           (when bcs
             (setf (varef vec index) (* (varef bcs d) (varef vec index)))))))))
  vec)


(defun iterate-over-dimensions (nd func)
  (declare (fixnum nd)
           (function func))
  (lparallel:pdotimes (d nd)
    (funcall func d)))
(defun iterate-over-dimensions-serial (nd func)
  (declare (fixnum nd)
           (function func))
  (dotimes (d nd)
    (funcall func d)))


(defun test-proj (sim)
  (let* ((d 0)
         (E (cl-mpm/aggregate::sim-global-e sim))
         (vi (assemble-global-vec sim #'cl-mpm/mesh::node-internal-force d))
         (va (magicl:@
              (magicl:transpose E)
              vi))
         )
    (pprint vi)
    (pprint va)
    (project-int-vec
     sim
     va
     #'cl-mpm/mesh::node-internal-force
     d)
    (pprint (magicl:@ E (assemble-internal-vec sim #'cl-mpm/mesh::node-internal-force d)))
    ;; (pprint (magicl:@ E (assemble-internal-vec sim #'cl-mpm/mesh::node-internal-force d)))
    ))
