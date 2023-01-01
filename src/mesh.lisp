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
    #:node-temperature
    #:node-dtemp
    #:get-node
    #:reset-node
    #:in-bounds
    #:position-to-index
    #:index-to-position
  ))

(in-package :cl-mpm/mesh)
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(defclass node ()
  ((mass
     :accessor node-mass
     :type double-float
     :initform 0d0
     )
   (volume
    :accessor node-volume
    :type double-float
    :initform 0d0
    )
   (index
    :accessor node-index
    :initarg :index)
   (position
    :accessor node-position
    :initarg :position
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (magicl:zeros '(1 1) :type 'double-float))
  (acceleration
    :accessor node-acceleration
    :initarg :acceleration
     :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (magicl:zeros '(1 1) :type 'double-float)
    )
  (force
    :accessor node-force
    :initarg :force
     :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (magicl:zeros '(1 1) :type 'double-float))
  (velocity
    :accessor node-velocity
    :initarg :velocity
     :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (magicl:zeros '(1 1) :type 'double-float))
  (lock 
    :accessor node-lock
    :initform (sb-thread:make-mutex)))
  (:documentation "A node on the computational mesh"))

(defclass node-fracture (node)
  ((strain-energy-density
     :accessor node-strain-energy-density
     :type double-float
     :initform 0d0
     )))

(defclass node-thermal (node)
  ((temperature
     :accessor node-temperature
     :type double-float
     :initform 0d0
     )
   (temperature-gradient
     :accessor node-dtemp
     :type double-float
     :initform 0d0
     )))

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
      :initarg :shape-func)
   (boundary-shapes
    :accessor mesh-boundary-shapes
    :initform nil)
   (boundary-order
    :accessor mesh-boundary-order
    :initarg :boundary-order))
    (:documentation "MPM computational mesh"))

(defun make-node (index pos)
  "Default initialise a 2d node at pos"
  (make-instance 'node-thermal
                 :force (magicl:zeros (list 2 1) :type 'double-float)
                 :velocity (magicl:zeros (list 2 1) :type 'double-float)
                 :acceleration (magicl:zeros (list 2 1) :type 'double-float)
                 :index (magicl:from-list (mapcar (lambda (x) (coerce x 'double-float))
                                                  index)
                                          (list 2 1) :type 'double-float)
                 :position pos
                 ))

(defun make-nodes (mesh size)
  "Make a 2d mesh of specific size"
  (make-array size :initial-contents 
      (loop for x from 0 to (- (nth 0 size) 1)
            collect (loop for y from 0 to (- (nth 1 size) 1)
                          collect (make-node (list x y)
                                             (magicl:from-list (index-to-position mesh (list x y)) '(2 1)))))))

(defun make-mesh (size resolution shape-function)
  "Create a 2D mesh and fill it with nodes"
  (let* ((size (mapcar (lambda (x) (coerce x 'double-float)) size))
         (resolution (coerce resolution 'double-float))
         (boundary-order 0;(* 2 (- (cl-mpm/shape-function::order shape-function) 1))
           )
         (meshcount (loop for d in size collect (+ (floor d resolution) 1 (* boundary-order 2))))
         (nodes '()))
    (let ((mesh (make-instance 'mesh
                              :nD 2 
                              :mesh-size size
                              :mesh-count meshcount
                              :mesh-res resolution
                              :nodes nodes
                              :shape-func shape-function
                              :boundary-order boundary-order
                              )))
      (setf (mesh-nodes mesh) (make-nodes mesh meshcount))
      mesh)))

(defun query-boundary-shapes (mesh index)
  (not (member index (mesh-boundary-shapes mesh) :test #'equal)))

(declaim (inline in-bounds-1d))
(defun in-bounds-1d (mesh value dim)
  "Check a single dimension is inside a mesh"
  (and (>= value 0) (< value (nth dim (slot-value mesh 'mesh-count)))))
(defun in-bounds (mesh pos)
  "Check a position (list) is inside a mesh"
  (and (in-bounds-1d mesh (first pos) 0)
       (in-bounds-1d mesh (second pos) 1)
       (if (mesh-boundary-shapes mesh)
           (query-boundary-shapes mesh pos)
           t))
  ;; (every (lambda (x) x) (loop for d from 0 to (- (mesh-nd mesh) 1)
  ;;                             collect (in-bounds-1d mesh (nth d pos) d)))
  )


(defun position-to-index (mesh pos &optional (round-operator 'round))
  "Turn a vector position into a list of indexes with rounding"
  (mapcar (lambda (x) (+ (funcall round-operator (/ (magicl:tref pos x 0) (mesh-resolution mesh)))
                         (mesh-boundary-order mesh))) '(0 1)))

(defun index-to-position (mesh index)
  "Turn a vector position into a list of indexes with rounding"
  (mapcar (lambda (i) (* (- i (mesh-boundary-order mesh)) (mesh-resolution mesh))) index))

(defun get-node (mesh pos)
  "Check bounds and get node"
  (if (in-bounds mesh pos)
      (apply #'aref (cons (mesh-nodes mesh) pos))
    (error (format nil "Access grid out of bounds at: ~a" pos))))

(defgeneric node-g2p (mp node svp dsvp grads)
  (:documentation "G2P behaviour for specific nodes"))
(defmethod node-g2p (mp node svp dsvp grads))
;; (defmethod node-g2p (mp node svp dsvp grads)
;;   (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
;;                    (node-temp cl-mpm/mesh:node-temperature)) node
;;     (with-accessors ((vel mp-velocity)
;;                      (mass mp-mass)
;;                      (temp cl-mpm/particle::mp-temperature)
;;                      (strain-rate cl-mpm/particle:mp-strain-rate)) mp
;;       (progn 
;;         (setf vel (magicl:.+ vel (magicl:scale node-vel svp)))))))
(defmethod node-g2p (mp (node node-thermal) svp dsvp grads)
  (with-accessors ((node-temp cl-mpm/mesh:node-temperature)) node
    (with-accessors ((temp cl-mpm/particle::mp-temperature)) mp
      (progn 
        (setf temp (+ temp (* node-temp svp))))))
  (call-next-method))

(defgeneric reset-node (node)
  (:documentation "Reset grid to default state"))

(defmethod reset-node ((node node))
  (with-slots ( (mass mass)
                (vel velocity)
               (volume volume)
                (force force))
                node
    (setf mass 0d0)
    (setf volume 0d0)
    (setf vel (magicl:scale vel 0))
    (setf force (magicl:scale force 0))))

(defmethod reset-node ((node node-thermal))
  (with-accessors ((temperature node-temperature)
                   (dtemp node-dtemp))
                node
    (setf temperature 0d0)
    (setf dtemp 0d0))
  (call-next-method))


