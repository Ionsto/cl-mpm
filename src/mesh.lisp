(defpackage :cl-mpm
  (:use :cl))

(in-package :cl-mpm)

(defclass node ()
  ((mass :initform 0)
  (index :initarg :index)
  (acceleration :initform (magicl:zeros '(2 1)))
  (force :initform (magicl:zeros '(2 1)))
  (velocity :initform (magicl:zeros '(2 1))))
  (:documentation "A node on the computational mesh"))


(defclass mesh ()
  ( (nD :initarg :nD)
    (mesh-count :initarg :mesh-count)
    (mesh-size :initarg :mesh-size)
    (mesh-res :initarg :mesh-res)
    (nodes :initarg :nodes :accessor get-nodes)
    (shape-func :initarg :shape-func))
    (:documentation "MPM computational mesh"))

(defun make-node (x y)
  (make-instance 'node 
                 :index (magicl:from-list (list x y) '(2 1))))

(defun make-mesh (nD size resolution shape-function)
  (let* ((meshcount (loop for d in size collect (+ (floor d resolution) 1)))
         (nodes (loop for x from 0 to (- (nth 0 meshcount) 1)
                      collect (loop for y from 0 to (- (nth 1 meshcount) 1)
                      collect (make-node x y)))))
    (progn
    (make-instance 'mesh
      :nD nD
      :mesh-size size
      :mesh-count meshcount
      :mesh-res resolution
      :nodes (make-array meshcount :initial-contents nodes)
      :shape-func shape-function
      ))))

(defun in-bounds (mesh value dim)
  "Check a single dimension is inside a mesh"
  (and (>= value 0) (< value (nth dim (slot-value mesh 'mesh-count)))))

(defun get-node (mesh x y)
  "Check bounds and get node"
  (if (and (in-bounds mesh x 0) (in-bounds mesh y 1))
    (aref (get-nodes mesh) x y)
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

