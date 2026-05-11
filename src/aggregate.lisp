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

(defstruct restrictor
  (mat))
(defstruct prologator
  (mat))

;; (defgeneric restrict (sim mat)
;;   (case (magicl:ncols vec 1)
;;     (cl-mpm/fastmaths::fast-@-sparse-mat-dense-vec (cl-mpm/aggregate::sim-global-e-sparse sim))))
(defgeneric prologate (sim mat)
  )

(defun resolve (sim op vec)
  )

(defun @-mass-matrix-vec (sim vec d)
  (let ((et (cl-mpm/aggregate::sim-global-sparse-et sim))
        (e (cl-mpm/aggregate::sim-global-sparse-e sim))
        (sma (cl-mpm/aggregate::sim-global-sparse-ma sim)))
    (cl-mpm/fastmaths::fast-@-sparse-mat-dense-vec-masked
     et
     (cl-mpm/fastmaths::fast-.*
      sma
      (cl-mpm/fastmaths::fast-@-sparse-mat-dense-vec-masked
       e
       vec
       (aref (sim-global-bcs sim) d)
       (aref (sim-global-bcs-int sim) d)
       ))
     (aref (sim-global-bcs-int sim) d)
     (aref (sim-global-bcs sim) d)
     )))

(defun aggregate-vec (sim vec d)
  (if (sim-global-bcs sim)
      (cl-mpm/fastmaths::fast-@-sparse-mat-dense-vec-masked
       (cl-mpm/aggregate::sim-global-sparse-et sim)
       vec
       (aref (sim-global-bcs-int sim) d)
       (aref (sim-global-bcs sim) d)
       )
      (cl-mpm/fastmaths::fast-@-sparse-mat-dense-vec
       (cl-mpm/aggregate::sim-global-sparse-et sim)
       vec))
  ;; (cl-mpm/fastmaths::fast-@-sparse-mat-dense-vec
  ;;  (cl-mpm/aggregate::sim-global-sparse-et sim)
  ;;  vec)
  ;; (magicl:@
  ;;  (magicl:transpose (cl-mpm/aggregate::sim-global-e sim))
  ;;  vec)
  )
(defun extend-vec (sim vec d)
  (if (sim-global-bcs sim)
      (cl-mpm/fastmaths::fast-@-sparse-mat-dense-vec-masked
       (cl-mpm/aggregate::sim-global-sparse-e sim)
       vec
       (aref (sim-global-bcs sim) d)
       (aref (sim-global-bcs-int sim) d)
       )
      (cl-mpm/fastmaths::fast-@-sparse-mat-dense-vec
       (cl-mpm/aggregate::sim-global-sparse-e sim)
       vec))
  ;; (cl-mpm/fastmaths::fast-@-sparse-mat-dense-vec
  ;;  (cl-mpm/aggregate::sim-global-sparse-e sim)
  ;;  vec)
  ;; (magicl:@
  ;;  (cl-mpm/aggregate::sim-global-e sim)
  ;;  vec)
  )


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
   (global-sparse-e
    :initform nil
    :accessor sim-global-sparse-e)
   (global-sparse-et
    :initform nil
    :accessor sim-global-sparse-et)
   (global-sparse-ma
    :initform nil
    :accessor sim-global-sparse-ma)
   (global-bcs
    :initform nil
    :accessor sim-global-bcs)
   (global-bcs-int
    :initform nil
    :accessor sim-global-bcs-int)
   (global-ma
    :initform nil
    :accessor sim-global-ma)
   (internal-v
    :initform nil
    :accessor sim-internal-v)))

(defun iterate-over-cell-patch-2d (sim node ring-size func)
  (declare (function func))
  (with-accessors ((mesh cl-mpm::sim-mesh))
      sim
    (with-accessors ((index cl-mpm/mesh::node-index))
        node
      (destructuring-bind (ix iy iz) index
        (declare (fixnum ix iy)
                 (ignore iz))
        (loop for dx from (- (+ ring-size 1)) to ring-size
              do
                 (when (cl-mpm/mesh::in-bounds-cell-1d mesh (+ ix dx) 0)
                   (loop for dy from (- (+ ring-size 1)) to ring-size
                         do
                            (when (cl-mpm/mesh::in-bounds-cell-1d mesh (+ iy dy) 1)
                              (let ((cell-index (mapcar #'+ index (list dx dy 0))))
                                ;; (when (cl-mpm/mesh::in-bounds-cell mesh cell-index))
                                (funcall func (cl-mpm/mesh::get-cell mesh cell-index)))))))))))

(defun iterate-over-cell-patch-3d (sim node ring-size func)
  (declare (function func))
  (with-accessors ((mesh cl-mpm::sim-mesh))
      sim
    (with-accessors ((index cl-mpm/mesh::node-index))
        node
      (destructuring-bind (ix iy iz) index
        (loop for dx from (- ring-size 1) to ring-size
              do
                 (when (cl-mpm/mesh::in-bounds-cell-1d mesh (+ ix dx) 0)
                   (loop for dy from (- ring-size 1) to ring-size
                         do
                            (when (cl-mpm/mesh::in-bounds-cell-1d mesh (+ iy dy) 1)
                              (loop for dz from (- ring-size 1) to ring-size
                                    do
                                       (when (cl-mpm/mesh::in-bounds-cell-1d mesh (+ iz dz) 2)
                                         (let ((cell-index (mapcar #'+ index (list dx dy dz))))
                                           ;; (when (cl-mpm/mesh::in-bounds-cell mesh cell-index))
                                           (funcall func (cl-mpm/mesh::get-cell mesh cell-index)))))))))))))

(defun iterate-over-cell-patch (sim node ring-size func)
  (declare (function func))
  (if (= (cl-mpm/mesh::mesh-nd (cl-mpm:sim-mesh sim)) 2)
      (iterate-over-cell-patch-2d sim node ring-size func)
      (iterate-over-cell-patch-3d sim node ring-size func)))

(defun get-closest-cell (sim node)
  (let* ((closest-elem nil)
         (dist 0d0)
         (volume-ratio-min 0.25d0)
         (mesh (cl-mpm:sim-mesh sim))
         (volume-t (expt (cl-mpm/mesh::mesh-resolution mesh) (cl-mpm/mesh:mesh-nd mesh)))
         (pos (cl-mpm/mesh::node-position node)))
    (flet ((check-cell (cell)
             (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                              (index cl-mpm/mesh::cell-index)
                              (centroid cl-mpm/mesh::cell-centroid)
                              (active cl-mpm/mesh::cell-active)
                              (volume cl-mpm/mesh::cell-volume)
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
                          ;; (> volume (* volume-ratio-min volume-t))
                          (> dist dist-tr))
                     (setf dist dist-tr
                           closest-elem cell)))))))
      ;; (iterate-over-cell-patch
      ;;  sim
      ;;  node
      ;;  1
      ;;  #'check-cell)
      (unless closest-elem
        (iterate-over-cell-patch
         sim
         node
         2
         #'check-cell)
        (unless closest-elem
          (let ((mutex (sb-thread:make-mutex)))
            (cl-mpm::iterate-over-cells
             (cl-mpm:sim-mesh sim)
             (lambda (cell)
               (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                                (index cl-mpm/mesh::cell-index)
                                (centroid cl-mpm/mesh::cell-centroid)
                                (active cl-mpm/mesh::cell-active)
                                (volume cl-mpm/mesh::cell-volume)
                                (agg cl-mpm/mesh::cell-agg))
                   cell
                 (when (and
                        (cl-mpm/mesh::cell-active cell)
                        (not (cl-mpm/mesh::cell-partial cell))
                        ;; (> volume (* volume-ratio-min volume-t))
                        (not (cl-mpm/mesh::cell-agg cell)))
                   (let ((dist-tr (cl-mpm/fastmaths::diff-norm
                                   pos
                                   centroid)))
                     ;;Double checked lock
                     (when (or
                            (not closest-elem)
                            (> dist dist-tr))
                       (sb-thread:with-mutex (mutex)
                         (when (or
                                (not closest-elem)
                                (> dist dist-tr))
                           (setf dist dist-tr
                                 closest-elem cell)))))))))))))
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





(defun locate-aggregate-nodes (sim)
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
      (lparallel:pdotimes (fdc (length nodes))
        do (let ((n (aref nodes fdc)))
             (setf (cl-mpm/mesh::node-agg-fdc n) fdc))))

    (let ((nodes (cl-mpm/aggregate::sim-agg-nodes-fd sim)))
      (lparallel:pdotimes (fd (length nodes))
        do (let ((n (aref nodes fd)))
             (setf (cl-mpm/mesh::node-agg-fd n) fd))))))

(defun assemble-global-scalar (sim accessor &optional (result nil))
  (declare (function accessor))
  (let* ((active-nodes (sim-agg-nodes-fd sim))
         (ndof (length active-nodes))
         (v (if result
                result
                (cl-mpm/utils::arb-matrix ndof 1))))
    (cl-mpm::iterate-over-nodes-array
     active-nodes
     (lambda (n)
       (setf (varef v (cl-mpm/mesh::node-agg-fd n)) (funcall accessor n))))
    (values v)))

(defun assemble-global-vec (sim accessor dim &optional (res nil))
  (declare (function accessor))
  (let* ((active-nodes (sim-agg-nodes-fd sim))
         (ndof (length active-nodes))
         (v (if res res (cl-mpm/utils::arb-matrix ndof 1))))
    (cl-mpm::iterate-over-nodes-array
     active-nodes
     (lambda (n)
       (setf (varef v (cl-mpm/mesh::node-agg-fd n)) (varef (funcall accessor n) dim))))
    (values v)))

(defun assemble-internal-vec (sim accessor dim &optional (res nil))
  (declare (function accessor))
  (let* ((int-nodes (sim-agg-nodes-fdc sim))
         (ndof (length int-nodes))
         (v (if res res (cl-mpm/utils::arb-matrix ndof 1))))
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
  (let* ((agg-nodes (sim-agg-nodes-fdc sim))
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
                 (cl-mpm::iterate-over-cell-shape-local
                  mesh
                  (cl-mpm/mesh::node-agg-interior-cell node)
                  (cl-mpm/mesh::node-position node)
                  (lambda (cn weight grads)
                    (incf (mtref e (cl-mpm/mesh::node-agg-fd node) (cl-mpm/mesh::node-agg-fdc cn))
                          weight)))))))))
    (values e)))


;; (eval-when (:compile-toplevel :load-toplevel :execute)
;;   (sb-ext::defglobal *assembly-counter* 0)
;;   (declaim (fixnum *assembly-counter*)))
(let* ((est-size 1000)
       ;; (v (make-array est-size :fill-pointer 0 :adjustable t :element-type 'double-float))
       ;; (r (make-array est-size :fill-pointer 0 :adjustable t :element-type 'fixnum))
       ;; (c (make-array est-size :fill-pointer 0 :adjustable t :element-type 'fixnum))
       (counter (make-array 1 :element-type '(UNSIGNED-BYTE 64)))
       )
  (defun assemble-sparse-e (sim)
    (let* ((agg-nodes (sim-agg-nodes-fdc sim))
           (active-nodes (sim-agg-nodes-fd sim))
           (mesh (cl-mpm:sim-mesh sim))
           (nd (cl-mpm/mesh:mesh-nd mesh))
           (ndof (length active-nodes))
           (ndofC (length agg-nodes))
           (est-size (+ ndofC (* (expt 2 nd) (- ndof ndofC))))
           (v (make-array est-size :fill-pointer est-size :adjustable t :element-type 'double-float))
           (r (make-array est-size :fill-pointer est-size :adjustable t :element-type 'fixnum))
           (c (make-array est-size :fill-pointer est-size :adjustable t :element-type 'fixnum))
           ;; (lock (sb-thread:make-mutex))
           )
      ;; (setf *assembly-counter* 0)
      (sb-ext:atomic-update (aref counter 0) (lambda (a) 0))
      ;; (setf v (adjust-array v est-size :fill-pointer 0))
      ;; (setf r (adjust-array r est-size :fill-pointer 0))
      ;; (setf c (adjust-array c est-size :fill-pointer 0))
      ;; (setf 
      ;;  (fill-pointer v) 0
      ;;  (fill-pointer r) 0
      ;;  (fill-pointer c) 0)
      (cl-mpm::iterate-over-nodes
       mesh
       (lambda (node)
         (when (cl-mpm/mesh:node-active node)
           (when (cl-mpm/mesh::node-agg node)
             ;; (sb-thread:with-mutex (lock))
             (if (cl-mpm/mesh::node-interior node)
                 (progn
                   ;;Interior -> populate 1s
                   (let ((place (sb-ext:atomic-incf (aref counter 0))))
                     (setf (aref v place) 1d0)
                     (setf (aref r place) (cl-mpm/mesh::node-agg-fd node))
                     (setf (aref c place) (cl-mpm/mesh::node-agg-fdc node))
                     ))
                 (progn
                   ;;Aggregate -> populate svp
                   (cl-mpm::iterate-over-cell-shape-local
                    mesh
                    (cl-mpm/mesh::node-agg-interior-cell node)
                    (cl-mpm/mesh::node-position node)
                    (lambda (cn weight grads)
                      (when (> (abs weight) 1d-9)
                        (let ((place (sb-ext:atomic-incf (aref counter 0))))
                          (setf (aref v place) weight)
                          (setf (aref r place) (cl-mpm/mesh::node-agg-fd node))
                          (setf (aref c place) (cl-mpm/mesh::node-agg-fdc cn))))))))))))
      (values (cl-mpm/utils::build-sparse-matrix v r c ndof ndofC)
              (cl-mpm/utils::build-sparse-matrix v c r ndofC ndof)))))

;; (defun assemble-implicit-e (sim)
;;   (let ((fdc 0))
;;     (loop for n across (cl-mpm/aggregate::sim-agg-nodes-fdc sim)
;;           do (progn
;;                (setf (cl-mpm/mesh::node-agg-fdc n) fdc)
;;                (incf fdc))))
;;   (let ((fd 0))
;;     (loop for n across (cl-mpm/aggregate::sim-agg-nodes-fd sim)
;;           do (progn
;;                (setf (cl-mpm/mesh::node-agg-fd n) fd)
;;                (incf fd))))
;;   )

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

(defgeneric update-aggregate-elements (sim))
(defmethod update-aggregate-elements ((sim mpm-sim-aggregated))
  (when (sim-enable-aggregate sim)
    (locate-aggregate-nodes sim)
    (let ()
      (multiple-value-bind (e et) (assemble-sparse-e sim)
        (setf (sim-global-sparse-e sim) e
              (sim-global-sparse-et sim) et))
      (let ((nd (cl-mpm/mesh::mesh-nd (cl-mpm:sim-mesh sim))))
        (setf (sim-global-bcs sim) (make-array nd :element-type t))
        (setf (sim-global-bcs-int sim) (make-array nd :element-type t))
        (iterate-over-dimensions
         nd
         (lambda (d)
           (setf (aref (sim-global-bcs sim) d) (assemble-global-bcs sim d))
           (setf (aref (sim-global-bcs-int sim) d) (assemble-internal-bcs sim d)))))
      (update-mass-matrix sim))))

(defun update-mass-matrix (sim)
  "Assuming that we already have an extension matrix, attempt to update the mass matrix directly"
  (when (sim-enable-aggregate sim)
    (let (;; (E (cl-mpm/aggregate::sim-global-e sim))
          )
      ;; (setf (sim-global-ma sim) (magicl:@
      ;;                            (magicl:transpose E)
      ;;                            (cl-mpm/aggregate::assemble-global-mass-matrix sim)
      ;;                            E))
      (when (sim-global-sparse-ma sim)
          (unless (= (magicl:nrows (sim-global-sparse-ma sim))
                     (length (cl-mpm/aggregate::sim-agg-nodes-fd sim)))
            (setf (sim-global-sparse-ma sim) nil)))
      (setf
       (sim-global-sparse-ma sim)
       (assemble-global-scalar sim
                               #'cl-mpm/mesh::node-mass
                               (sim-global-sparse-ma sim))))))

(defun assemble-global-bcs (sim d)
  (let* ((int-nodes (cl-mpm/aggregate::sim-agg-nodes-fd sim))
         (ndof (length int-nodes))
         (v (cl-mpm/utils::arb-matrix ndof 1)))
    (cl-mpm::iterate-over-nodes-array
     int-nodes
     (lambda (n)
       (let ((index (cl-mpm/mesh::node-agg-fd n)))
         (let ((bcs (cl-mpm/mesh::node-bcs n)))
           (if bcs
               (setf (varef v index) (varef bcs d))
               (setf (varef v index) 1d0))))))
    v))

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

(declaim (notinline linear-solve-with-bcs)
         (ftype (function (cl-mpm::mpm-sim
                           (or magicl::matrix/double-float null)
                           magicl::matrix/double-float
                           fixnum
                           &optional
                           (or magicl::matrix/double-float null))
                          magicl::matrix/double-float)
                linear-solve-with-bcs))
(defun linear-solve-with-bcs (sim ma v d &optional (target-vi nil))
  (let ((target-vi (if target-vi
                       target-vi
                       (cl-mpm/utils::arb-matrix (magicl:nrows v) 1))))
    (cl-mpm/fastmaths:fast-zero target-vi)
    (let* ((et (cl-mpm/aggregate::sim-global-sparse-et sim))
           (e (cl-mpm/aggregate::sim-global-sparse-e sim))
           (sma (cl-mpm/aggregate::sim-global-sparse-ma sim))
           (bcs (aref (sim-global-bcs-int sim) d))
           (gbcs (assemble-global-bcs sim d))
           (work-vec (cl-mpm/utils::arb-matrix (magicl:nrows sma) 1))
           (work-vec-agg (cl-mpm/utils::arb-matrix (magicl:nrows v) 1)))
      ;; (loop for i from 0 below (magicl:nrows bcs)
      ;;       do (setf (cl-mpm/utils:varef bcs i) 1d0))
      ;; (pprint bcs)
      ;; (pprint v)
      ;; (pprint sma)
      (cl-mpm/linear-solver::solve-conjugant-gradients
       (lambda (x)
         (@-mass-matrix-vec sim x d)
         ;; (cl-mpm/fastmaths::fast-@-sparse-mat-dense-vec
         ;;  e
         ;;  x
         ;;  work-vec)
         ;; (cl-mpm/fastmaths::fast-.*
         ;;  gbcs
         ;;  work-vec
         ;;  work-vec)
         ;; (cl-mpm/fastmaths::fast-.*
         ;;  sma
         ;;  work-vec
         ;;  work-vec)
         ;; (cl-mpm/fastmaths::fast-@-sparse-mat-dense-vec
         ;;  et
         ;;  work-vec
         ;;  work-vec-agg)
         ;; (cl-mpm/fastmaths::fast-@-sparse-mat-dense-vec
         ;;  et
         ;;  (cl-mpm/fastmaths::fast-.*
         ;;   gbcs
         ;;   (cl-mpm/fastmaths::fast-.*
         ;;    sma
         ;;    (cl-mpm/fastmaths::fast-.*
         ;;     gbcs
         ;;     (cl-mpm/fastmaths::fast-@-sparse-mat-dense-vec
         ;;      e
         ;;      x)))))
         )
       v
       :tol 1d-20
       :max-iters 10000
       :mask bcs
       ))))


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
    (let ((ma (cl-mpm/aggregate::sim-global-ma sim)))
      (iterate-over-dimensions
       (cl-mpm/mesh:mesh-nd mesh)
       (lambda (d)
         (let* ((vi (aggregate-vec
                     sim
                     (assemble-global-vec sim #'cl-mpm/mesh::node-velocity d)
                     d)))
           (project-global-vec
            sim
            (extend-vec sim (linear-solve-with-bcs sim ma vi d) d)
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
      ;; (cl-mpm::fast-scale! (cl-mpm/aggregate::sim-global-ma sim) (cl-mpm::sim-mass-scale sim))
      )))

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
      (update-aggregate-elements sim))))



(defmethod cl-mpm::update-node-forces ((sim mpm-sim-aggregated))
  (if (sim-enable-aggregate sim)
      (cl-mpm::update-node-forces-standard sim)
      (update-node-forces-agg sim (cl-mpm:sim-dt sim))))


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
         (cl-mpm::calculate-forces node damping 0d0 mass-scale))))

    ;; (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt)
    ;;For each aggregated element set solve mass matrix and velocity
    (when enable-aggregate
      (let* ((ma (cl-mpm/aggregate::sim-global-ma sim)))
        (iterate-over-dimensions
         (mesh-nd mesh)
         (lambda (d)
           (let ((fa
                   (aggregate-vec
                    sim
                    (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-force d)
                    d)))
             (let* ((acc (linear-solve-with-bcs sim ma fa d)))
               (cl-mpm/fastmaths::fast-scale! acc (/ 1d0 (sqrt mass-scale)))
               (project-global-vec
                sim
                (extend-vec sim acc d)
                #'cl-mpm/mesh::node-acceleration d)))))))

    ;; (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt)
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
    ;;   (project-velocity sim))
    ;; (cl-mpm::apply-bcs (cl-mpm:sim-mesh sim) (cl-mpm:sim-bcs sim) dt)
    ))


(defun project-displacement (sim)
  (with-accessors ((sim-mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)
                   (dt cl-mpm:sim-dt))
      sim
    (iterate-over-dimensions
     (cl-mpm/mesh::mesh-nd mesh)
     (lambda (d)
       (let ((vel-proj (cl-mpm/aggregate::assemble-internal-vec sim #'cl-mpm/mesh::node-displacment d)))
         (cl-mpm/aggregate::project-global-vec
          sim
          (extend-vec sim vel-proj d)
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
       (let (
             (vel-proj (cl-mpm/aggregate::assemble-internal-vec sim #'cl-mpm/mesh::node-velocity d)))
         (cl-mpm/aggregate::project-global-vec
          sim
          (apply-global-bcs
           sim
           (extend-vec sim vel-proj d)
           d)
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
       (let ((vel-proj (cl-mpm/aggregate::assemble-internal-vec sim #'cl-mpm/mesh::node-acceleration d)))
         (cl-mpm/aggregate::project-global-vec
          sim
          (apply-global-bcs
           sim
           (extend-vec sim vel-proj d)
           d)
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
    ;; (when agg
    ;;   (project-displacement sim))
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

(defun apply-global-bcs (sim vec d)
  (declare (magicl::matrix/double-float vec)
           (fixnum d))
  (let* ((agg-nodes (cl-mpm/aggregate::sim-agg-nodes-fd sim))
         (vs (cl-mpm/utils:fast-storage vec)))
    (declare ((simple-array double-float (*)) vs))
    ;; (when (not (= (magicl:nrows vec) (length (cl-mpm/aggregate::sim-agg-nodes-fd sim))))
    ;;   (error "Incorrect BCs applied"))
    (cl-mpm::iterate-over-nodes-array
     agg-nodes
     (lambda (n)
       (let ((index (cl-mpm/mesh::node-agg-fd n)))
         (let ((bcs (cl-mpm/mesh::node-bcs n)))
           (when bcs
             (setf (aref vs index)
                   (the double-float
                        (*
                         (varef bcs d)
                         (varef vec index)))))))
       (values))))
  vec)


(defun iterate-over-dimensions (nd func)
  (declare (fixnum nd)
           (function func))
  (lparallel:pdotimes (d nd)
    (funcall func d)))

(defun iterate-over-dimensions-with-mutex (nd func)
  (declare (fixnum nd)
           (function func))
  (let ((mutex (sb-thread:make-mutex)))
    (lparallel:pdotimes (d nd)
      (funcall func d mutex))))

(defun iterate-over-dimensions-serial (nd func)
  (declare (fixnum nd)
           (function func))
  (dotimes (d nd)
    (funcall func d)))


(defun test-proj (sim)
  (let* ((d 0)
         (E (cl-mpm/aggregate::sim-global-e sim))
         (vi (assemble-global-vec sim #'cl-mpm/mesh::node-internal-force d))
         (va (aggregate-vec
              sim
              vi
              d))
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

(defmethod cl-mpm::calculate-min-dt-mps ((sim cl-mpm/aggregate::mpm-sim-aggregated))
  "Estimate minimum p-wave modulus"
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mass-scale cl-mpm::sim-mass-scale)
                   (enable-agg cl-mpm/aggregate::sim-enable-aggregate))
      sim
    ;; (pprint "Hello")
    (let* ((inner-factor most-positive-double-float)
           (h (cl-mpm/mesh:mesh-resolution mesh))
           (h2 (* h h)))
      (declare (double-float inner-factor mass-scale))
      (setf inner-factor
            (cl-mpm::reduce-over-nodes
             mesh
             (lambda (node)
               (if (and (cl-mpm/mesh::node-active node)
                        (or
                         (not
                          (cl-mpm/mesh::node-agg node))
                         (cl-mpm/mesh::node-interior node)))
                   (with-accessors ((node-active  cl-mpm/mesh:node-active)
                                    (pmod cl-mpm/mesh::node-pwave)
                                    (mass cl-mpm/mesh::node-mass)
                                    (svp-sum cl-mpm/mesh::node-svp-sum)
                                    (vol cl-mpm/mesh::node-volume)
                                    (vel cl-mpm/mesh::node-velocity)
                                    ) node
                     (if (and (> vol 0d0)
                              (> pmod 0d0))
                         (let ((nf (/ mass pmod)))
                           nf)
                         most-positive-double-float))
                   most-positive-double-float))
             #'min))
      ;; (pprint "Hello")
      ;; (when enable-agg
      ;;   (when (cl-mpm/aggregate::sim-global-ma sim)
      ;;     (let* ((E (magicl:transpose (cl-mpm/aggregate::sim-global-e sim)))
      ;;            (ma (cl-mpm/aggregate::sim-global-ma sim))
      ;;            (agg-inner inner-factor))
      ;;       (let* ((p-mod (cl-mpm/aggregate::assemble-global-scalar sim #'cl-mpm/mesh::node-pwave))
      ;;              (p-ratio (magicl:@ (assemble-internal-identity sim (magicl:@ E p-mod)) ma)))
      ;;         (pprint p-mod)
      ;;         (pprint p-ratio)
      ;;         (magicl:map! #'abs p-ratio)
      ;;         (loop for i from 0 below (magicl:nrows p-ratio)
      ;;               do (setf agg-inner
      ;;                        (min
      ;;                         agg-inner
      ;;                         (magicl:sum (magicl:row p-ratio i)))
      ;;                        ;; (min agg-inner (max 0d0 (/ 1d0 v)))
      ;;                        ))
      ;;         ;; (loop for v across (cl-mpm/utils:fast-storage p-ratio)
      ;;         ;;       do (setf agg-inner (min agg-inner (max 0d0 (/ 1d0 v)))))
      ;;         (pprint agg-inner)
      ;;         )
      ;;       (setf inner-factor (min inner-factor agg-inner)))))
      (if (< inner-factor most-positive-double-float)
          (* (sqrt mass-scale) (sqrt inner-factor) h)
          (cl-mpm:sim-dt sim)))))
