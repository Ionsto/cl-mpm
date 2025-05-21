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
  ((cell-list
    :initform nil
    :initarg :cell-list
    :accessor agg-cell-list)
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
                             (make-instance 'aggregate-element :cell-list (list cell neighbour))
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
         (dolist (cell (agg-cell-list agg-elem))
           (loop for n in (cl-mpm/mesh::cell-nodes cell)
                 do (setf (cl-mpm/mesh::node-agg n) t)))))))

(defun update-aggregate-elements (sim)
  )

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
       ))))


