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
    #:node-active
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
  ((active
    :accessor node-active
    :type boolean
    :initform nil
    )
   (mass
     :accessor node-mass
     :type double-float
     :initform 0d0
     )
   (volume
    :accessor node-volume
    :type double-float
    :initform 0d0
    )
   (volume-true
    :accessor node-volume-true
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
    :initform (magicl:zeros '(2 1) :type 'double-float))
  (acceleration
    :accessor node-acceleration
    :initarg :acceleration
     :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (magicl:zeros '(2 1) :type 'double-float)
    )
  (force
    :accessor node-force
    :initarg :force
     :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (magicl:zeros '(2 1) :type 'double-float))
  (velocity
    :accessor node-velocity
    :initarg :velocity
     :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (magicl:zeros '(2 1) :type 'double-float))
  (local-list
   :accessor node-local-list
   :initform (make-array 0 :fill-pointer 0 :adjustable t))
  (jacobian-inc
   :accessor node-jacobian-inc
   :initform 0d0)
   (jacobian
    :accessor node-jacobian
    :initform 0d0)
   (boundary-node
    :accessor node-boundary-node
    :type boolean
    :initform nil
    )
   (boundary-scalar
    :accessor node-boundary-scalar
    :initform 0d0)
   (p-wave
    :accessor node-pwave
    :type double-float
    :initform 0d0
    )
   (svp-sum
    :accessor node-svp-sum
    :type double-float
    :initform 0d0
    )

   (sdf
    :accessor node-sdf
    :initform 0d0)
   (damage
    :accessor node-damage
    :initform 0d0
    :type double-float)
  (lock
    :accessor node-lock
    :initform (sb-thread:make-mutex)))
  (:documentation "A node on the computational mesh"))

(defclass cell ()
  ((index
    :accessor cell-index
    :initarg :index)
   (nodes
    :accessor cell-nodes
    :initarg :nodes
    :type list
    :initform '())
   (neighbours
    :accessor cell-neighbours
    :type list
    :initform '())
   (mp-count
    :accessor cell-mp-count
    :initform 0)
   (volume
    :accessor cell-volume
    :initarg :volume
    :initform 0d0
    :type double-float)
   (centroid
    :accessor cell-centroid
    :initarg :centroid
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (magicl:zeros '(1 1) :type 'double-float))
   (deformation-gradient
       :accessor cell-deformation-gradient
       :type MAGICL:MATRIX/DOUBLE-FLOAT
       :initform (magicl:zeros '(2 2) :type 'double-float))
   (boundary
    :accessor cell-boundary
    :type boolean
    :initform nil)
   (pruned
    :accessor cell-pruned
    :type boolean
    :initform nil)
   )
   )

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
    (cells
      :accessor mesh-cells
      :initarg :cells)
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

(defun make-node (index pos h)
  "Default initialise a 2d node at pos"
  (make-instance 'node
                 ;; :force (magicl:zeros (list 2 1) :type 'double-float)
                 ;; :velocity (magicl:zeros (list 2 1) :type 'double-float)
                 ;; :acceleration (magicl:zeros (list 2 1) :type 'double-float)
                 :index (mapcar (lambda (x) (coerce x 'double-float))
                                 index)
                 :position pos
                 ))

(defun make-nodes (mesh size h)
  "Make a 2d mesh of specific size"
  (make-array size :initial-contents 
      (loop for x from 0 to (- (nth 0 size) 1)
            collect (loop for y from 0 to (- (nth 1 size) 1)
                          collect (make-node (list x y)
                                             (magicl:from-list (index-to-position mesh (list x y)) '(2 1))
                                             h
                                             )))))

(defun cell-calculate-centroid (nodes)
  (let ((centroid (magicl:zeros '(2 1))))
    (loop for node in nodes
          do
             (magicl:.+ centroid (node-position node) centroid))
    (magicl:scale centroid (/ 1d0 (length nodes)))))
(defun make-cell (mesh index h)
  ;;Get local nodes
  (let* ((nodes (loop for x from 0 to 1
                      append
                      (loop for y from 0 to 1
                            collect
                            (get-node mesh (mapcar #'+
                                                   index
                                                   (list x y))))))
         (volume (expt h 2))
         (centroid (cell-calculate-centroid nodes))
         (centroid (magicl:.+ (magicl:from-list (index-to-position mesh index) '(2 1) :type 'double-float)
                              (magicl:scale! (magicl:from-list (list h h) '(2 1) :type 'double-float) 0.5d0)))
         )
    (loop for n in nodes
          do (incf (node-volume-true n) (/ volume (length nodes))))
    (make-instance 'cell
                   :index index
                   :nodes nodes
                   :centroid centroid
                   :volume volume
                   )))
(defun make-cells (mesh size h)
  "Make a 2d mesh of specific size"
  (let ((cells (make-array (mapcar #'- size '(1 1)) :initial-contents
                          (loop for x from 0 to (- (nth 0 size) 2)
                                collect (loop for y from 0 to (- (nth 1 size) 2)
                                              collect (make-cell
                                                       mesh
                                                       (list x y)
                                                       h))))))
    (array-operations/utilities:nested-loop (i j) (array-dimensions cells)
      (let* ((cell (aref cells i j)))
        (with-accessors ((neighbours cell-neighbours))
            cell
          (setf neighbours (list))
          (loop for dx from -1 to 1
                do
                   (loop for dy from -1 to 1
                         do
                            (unless (and (= dx 0) (= dy 0))
                              (let ((di (mapcar #'+ (list i j) (list dx dy))))
                                (when (in-bounds-cell mesh di)
                                  (push (apply #'aref cells di) neighbours)))))))))
    cells))

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
      (setf (mesh-nodes mesh) (make-nodes mesh meshcount resolution))
      (setf (mesh-cells mesh) (make-cells mesh meshcount resolution))
      mesh)))

(defun in-bounds-cell (mesh index)
  (destructuring-bind (x y) index
    (and
     (and (>= x 0) (< x (- (nth 0 (mesh-count mesh)) 1)))
     (and (>= y 0) (< y (- (nth 1 (mesh-count mesh)) 1))))))

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

(declaim
 (inline position-to-index-array)
 (ftype (function (cl-mpm/mesh::mesh
                   magicl:matrix/double-float &optional function)
                  (simple-array double-float)) position-to-index-array))
(defun position-to-index-array (mesh pos &optional (round-operator #'round))
  "Turn a vector position into a list of indexes with rounding"
  (let ((h (mesh-resolution mesh))
        (p-a (magicl::matrix/double-float-storage pos)))
    (aops:each (lambda (x) (/ (the double-float x) h)) p-a))
  ;; (let ((res (make-array 2
  ;;                        :initial-contents (the (simple-array double-float) (magicl::matrix/double-float-storage pos))
  ;;                        :element-type 'double-float))
  ;;       (h (mesh-resolution mesh)))
  ;;   (declare (type double-float h))
  ;;   ;; (declare (type (simple-array double-float 2) res))
  ;;   (loop for i from 0 to 1
  ;;         do (setf (aref res i) (/ (the double-float (aref res i))
  ;;                                  h)))
  ;;   res)
  ;; (mapcar (lambda (x) (+ (funcall round-operator (/ (magicl:tref pos x 0) (mesh-resolution mesh)))
  ;;                        (mesh-boundary-order mesh))) '(0 1))
  )

(declaim
 (inline position-to-index)
 (ftype (function (cl-mpm/mesh::mesh
                   magicl:matrix/double-float &optional function)
                  list) position-to-index))
(defun position-to-index (mesh pos &optional (round-operator #'round))
  (declare (type function round-operator)
           (type magicl:matrix/double-float pos))
  "Turn a vector position into a list of indexes with rounding"
  (mapcar (lambda (x) (funcall round-operator (/
                                                  (the double-float (magicl:tref pos x 0))
                                                  (the double-float (mesh-resolution mesh))))
                         ) '(0 1)))

(defun index-to-position-array (mesh index)
  "Turn a vector position into a list of indexes with rounding"
  (let ((h (mesh-resolution mesh)))
    (declare (double-float h))
    (aops:each (lambda (x) (* (the double-float x) h)) index)))
(declaim (ftype (function (mesh list) list) index-to-position))
(defun index-to-position (mesh index)
  "Turn a vector position into a list of indexes with rounding"
  (let ((h (mesh-resolution mesh)))
    (mapcar (lambda (i) (* (the double-float h) (coerce i 'double-float))) index)))

(declaim (inline get-node)
         (ftype (function (mesh list) node) get-node))
(defun get-node (mesh pos)
  "Check bounds and get node"
  (policy-cond:policy-if (> safety speed)
                         (if (in-bounds mesh pos)
                             (apply #'aref (mesh-nodes mesh) pos)
                             (error (format nil "Access grid out of bounds at: ~a" pos)))
                         (apply #'aref (mesh-nodes mesh) pos)))

(declaim (inline get-cell)
         (ftype (function (mesh list) cell) get-cell))
(defun get-cell (mesh pos)
  "Check bounds and get cell"
  (policy-cond:policy-if (> safety speed)
                         (if (in-bounds-cell mesh pos)
                             (apply #'aref (mesh-cells mesh) pos)
                             (error (format nil "Access grid out of bounds at: ~a" pos)))
                         (apply #'aref (mesh-cells mesh) pos)))

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
               (acc acceleration)
               (j-inc jacobian-inc)
               (j jacobian)
               (p-wave p-wave)
               (svp-sum svp-sum)
               (volume volume)
               (damage damage)
               (active active)
               (boundary boundary-node)
               (boundary-scalar boundary-scalar)
                (force force))
                node
    (setf active nil)
    (setf boundary nil)
    (setf mass 0d0)
    (setf volume 0d0)
    (setf p-wave 0d0)
    (setf damage 0d0)
    (setf svp-sum 0d0)
    (setf j 0d0)
    (setf boundary-scalar 0d0)
    (setf j-inc 0d0)
    (magicl:scale! vel 0d0)
    (magicl:scale! acc 0d0)
    (magicl:scale! force 0d0)))

(defmethod reset-node ((node node-thermal))
  (with-accessors ((temperature node-temperature)
                   (dtemp node-dtemp))
                node
    (setf temperature 0d0)
    (setf dtemp 0d0))
  (call-next-method))

(defun gauss-points (n h)
  (cond
    ((= n 1)
     (list (magicl:zeros '(2 1))))
    ((= n 2)
     (loop for x from -1 to 1 by 2
           append
           (loop for y from -1 to 1 by 2
                 collect
                 (magicl:scale! (magicl:from-list (list x y) '(2 1) :type 'double-float)
                                (/ h (* 2 (sqrt 3)))))))
    )
  )

(defun cell-quadrature-iterate-over-neighbours (mesh cell gp func)
  (declare (function func))
  (let ((h (cl-mpm/mesh:mesh-resolution mesh)))
    (with-accessors ((nodes cell-nodes)
                     (centroid cell-centroid)
                     (volume cell-volume))
        cell
      (dolist (point (gauss-points gp h))
        (let ((quad (magicl:.+ centroid point))
              (volume-ratio (/ 1 (expt gp 2))))
          (loop for node in nodes
                do
                   (with-accessors ((n-pos node-position))
                       node
                     (let* ((dist-vec (magicl:.- quad n-pos))
                            (dist (list (magicl:tref dist-vec 0 0) (magicl:tref dist-vec 1 0)))
                            (weights (mapcar (lambda (x) (cl-mpm/shape-function::shape-linear x h)) dist))
                            (weight (reduce #'* weights))
                            (grads (mapcar (lambda (d w) (* (cl-mpm/shape-function::shape-linear-dsvp d h) w))
                                           dist (nreverse weights)))
                            )
                       (funcall func
                                mesh
                                cell
                                quad
                                volume
                                ;(* volume volume-ratio)
                                node
                                weight
                                grads)))))))))

(defun cell-iterate-over-neighbours (mesh cell func)
  (declare (function func))
  (let ((h (cl-mpm/mesh:mesh-resolution mesh)))
    (with-accessors ((nodes cell-nodes)
                     (centroid cell-centroid)
                     (volume cell-volume))
        cell
      (loop for node in nodes
            do
               (with-accessors ((n-pos node-position))
                   node
                 (let* ((dist-vec (magicl:.- centroid n-pos))
                        (dist (list (magicl:tref dist-vec 0 0) (magicl:tref dist-vec 1 0)))
                        (weights (mapcar (lambda (x) (cl-mpm/shape-function::shape-linear x h)) dist))
                        (weight (reduce #'* weights))
                        (grads (mapcar (lambda (d w) (* (cl-mpm/shape-function::shape-linear-dsvp d h) w))
                                       dist (nreverse weights)))
                        )
                   (when (< 0d0 weight)
                     (funcall func
                              mesh
                              cell
                              centroid
                              volume
                              node
                              weight
                              grads)))))))
  )


;;; Printing methods

(defmethod print-object ((obj node) stream)
  (print-unreadable-object (obj stream :type t)
    (with-accessors ((index node-index))
        obj
      (format stream "index: ~a" index))))
