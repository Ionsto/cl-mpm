(defpackage :cl-mpm/mesh
  (:use :cl)
  (:export
    #:make-mesh
    #:mesh-mesh-size
    #:mesh-count
    #:mesh-nodes
    #:mesh-shape-func
    #:mesh-resolution
    #:mesh-nd
    #:node-velocity
    #:node-acceleration
    #:node-mass
    #:node-force
    #:node-lock
    #:node-index
    #:get-node
    #:reset-node
    #:in-bounds
    #:position-to-index
  ))

(in-package :cl-mpm/mesh)
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(defclass node ()
  ((mass 
     :accessor node-mass
     :initform 0
     )
  (index 
    :accessor node-index
    :initarg :index)
  (acceleration
    :accessor node-acceleration
    :initarg :acceleration
    :initform (magicl:zeros '(1 1))
    )
  (force 
    :accessor node-force
    :initarg :force
    :initform (magicl:zeros '(1 1)))
  (velocity
    :accessor node-velocity
    :initarg :velocity
    :initform (magicl:zeros '(1 1)))
  (lock 
    :accessor node-lock
    :initform (sb-thread:make-mutex)))
  (:documentation "A node on the computational mesh"))


(defclass mesh ()
  ( (nD
      :accessor mesh-nd
      :initarg :nD)
    (mesh-count 
      :accessor mesh-count
      :initarg :mesh-count)
    (mesh-size 
      :accessor mesh-mesh-size
      :initarg :mesh-size)
    (mesh-res 
      :accessor mesh-resolution
      :initarg :mesh-res)
    (nodes 
      :accessor mesh-nodes
      :initarg :nodes)
    (shape-func
      :accessor mesh-shape-func
      :initarg :shape-func))
    (:documentation "MPM computational mesh"))

(defun make-node (pos)
  "Default initialise a 2d node at pos"
  (make-instance 'node 
                 :force (magicl:zeros (list 2 1))
                 :velocity (magicl:zeros (list 2 1))
                 :acceleration (magicl:zeros (list 2 1))
                 :index (magicl:from-list pos (list 2 1))))

(defun make-nodes (size)
  "Make a 2d mesh of specific size"
  (make-array size :initial-contents 
      (loop for x from 0 to (- (nth 0 size) 1)
            collect (loop for y from 0 to (- (nth 1 size) 1)
                     collect (make-node (list x y))))))

(defun make-mesh (size resolution shape-function)
  "Create a 2D mesh and fill it with nodes"
  (let* ((size (mapcar (lambda (x) (coerce x 'double-float)) size))
         (resolution (coerce resolution 'double-float))
         (meshcount (loop for d in size collect (+ (floor d resolution) 1)))
         (nodes (make-nodes meshcount)))
    (make-instance 'mesh
      :nD 2 
      :mesh-size size
      :mesh-count meshcount
      :mesh-res resolution
      :nodes nodes
      :shape-func shape-function
      )))

(defun in-bounds-1d (mesh value dim)
  "Check a single dimension is inside a mesh"
  (and (>= value 0) (< value (nth dim (slot-value mesh 'mesh-count)))))
(defun in-bounds (mesh pos)
  "Check a position (list) is inside a mesh"
  (every (lambda (x) x) (loop for d from 0 to (- (mesh-nd mesh) 1)
      collect (in-bounds-1d mesh (nth d pos) d))))


(defun position-to-index (mesh pos &optional (round-operator 'round))
  "Turn a vector position into a list of indexes with rounding"
  (mapcar (lambda (x) (funcall round-operator (/ (magicl:tref pos x 0) (mesh-resolution mesh)))) '(0 1)))

(defun get-node (mesh pos)
  "Check bounds and get node"
  (if (in-bounds mesh pos)
      (apply #'aref (cons (mesh-nodes mesh) pos))
    (error (format nil "Access grid out of bounds at: ~a" pos))))

(defgeneric reset-node (node)
  (:documentation "Reset grid to default state"))

(defmethod reset-node (node)
  (with-slots ( (mass mass)
                (vel velocity)
                (force force))
                node
               (setf mass 0)
               (setf vel (magicl:scale vel 0))
               (setf force (magicl:scale force 0))))

