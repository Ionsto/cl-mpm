(defpackage :cl-mpm
  (:use :cl))

(in-package :cl-mpm)

(defclass node ()
  ((mass 
     :initform 0)
  (index :initarg :index)
  (acceleration
    :initarg :acceleration
    :initform (magicl:zeros '(1 1)))
  (force 
    :initarg :force
    :initform (magicl:zeros '(1 1)))
  (velocity
    :initarg :velocity
    :initform (magicl:zeros '(1 1))))
  (:documentation "A node on the computational mesh"))


(defclass mesh ()
  ( (nD :initarg :nD)
    (mesh-count :initarg :mesh-count)
    (mesh-size :initarg :mesh-size)
    (mesh-res :initarg :mesh-res)
    (nodes :initarg :nodes :accessor get-nodes)
    (shape-func :initarg :shape-func))
    (:documentation "MPM computational mesh"))

(defun make-node (nD pos)
  (make-instance 'node 
                 :force (magicl:zeros (list nD 1))
                 :velocity (magicl:zeros (list nD 1))
                 :acceleration (magicl:zeros (list nD 1))
                 :index (magicl:from-list pos (list nD 1))))



(defun make-nodes (nD size)
	(let* ((pos (loop for d from 0 to (- nD 1) 
					collect	(loop for x from 0 to (- (nth d size) 1) 
                        collect x)))
           (total-size (reduce #'* size)))
      (aops:reshape 
        (make-array total-size 
                    :initial-contents 
                    (loop for id in (apply #'alexandria:map-product #'list pos)
                        collect (apply #'make-node (cons nD (list id))))) size)))

(defun make-mesh (nD size resolution shape-function)
  (let* ((meshcount (loop for d in size collect (+ (floor d resolution) 1)))
         (nodes (make-nodes nd meshcount)))
    
    (make-instance 'mesh
      :nD nD
      :mesh-size size
      :mesh-count meshcount
      :mesh-res resolution
      :nodes nodes
      :shape-func shape-function
      )))

(defun in-bounds-value (mesh value dim)
  "Check a single dimension is inside a mesh"
  (and (>= value 0) (< value (nth dim (slot-value mesh 'mesh-count)))))
(defun in-bounds (mesh pos)
  "Check a single dimension is inside a mesh"
  (apply (lambda (x) x) (loop for d from 0 to (- (slot-value mesh 'nD) 1)
    collect (and (>= (nth d pos) 0) 
                 (< (nth d pos) (nth d (slot-value mesh 'mesh-count)))))))

(defun get-node (mesh pos)
  "Check bounds and get node"
  (if (in-bounds mesh pos)
      (apply #'aref (cons (get-nodes mesh) pos))
    (error "Access grid out of bounds")))

(defun reset-node (node)
  (with-slots ( (mass mass)
                 (vel velocity)
                 (force force))
                node
               (setf mass 0)
               (setf vel (magicl:scale vel 0))
               (setf force (magicl:scale force 0))))

(defun reset-grid (mesh)
  (let ((nodes (slot-value mesh 'nodes)))
  (dotimes (i (array-total-size nodes))
    (let ((node (row-major-aref nodes i)))
    (when (> (slot-value node 'mass) 0)
      (reset-node node))))))

(defun filter-grid (mesh mass-thresh)
  (let ((nodes (slot-value mesh 'nodes)))
  (dotimes (i (array-total-size nodes))
    (let ((node (row-major-aref nodes i)))
      (when (and (> (slot-value node 'mass) 0) 
                 (< (slot-value node 'mass) mass-thresh))
        (reset-node node))))))

