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
;(declaim (optimize (debug 0) (safety 0) (speed 3)))
(declaim (optimize (debug 3) (safety 3) (speed 0)))

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
    :initform (cl-mpm/utils:vector-zeros))
  (acceleration
    :accessor node-acceleration
    :initarg :acceleration
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils:vector-zeros)
    )
  (force
    :accessor node-force
    :initarg :force
     :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils:vector-zeros))
   (internal-force
    :accessor node-internal-force
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils:vector-zeros))
   (external-force
    :accessor node-external-force
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils:vector-zeros))

   (buoyancy-force
    :accessor node-buoyancy-force
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils:vector-zeros))
  (velocity
    :accessor node-velocity
    :initarg :velocity
     :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils:vector-zeros))
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
   (pressure
    :accessor node-pressure
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
    :initform (cl-mpm/utils:vector-zeros))
   (deformation-gradient
       :accessor cell-deformation-gradient
       :type MAGICL:MATRIX/DOUBLE-FLOAT
       :initform (cl-mpm/utils:matrix-zeros))
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
                 :index (mapcar (lambda (x) (coerce x 'fixnum)) index)
                 :position pos
                 ))

(defun make-nodes (mesh size h)
  "Make a 2d mesh of specific size"
  (make-array size :initial-contents
      (loop for x from 0 to (- (nth 0 size) 1)
            collect
            (loop for y from 0 to (- (nth 1 size) 1)
                          collect
                          (loop for z from 0 to (- (nth 2 size) 1)
                                collect
                                (make-node (list x y z)
                                           (cl-mpm/utils:vector-from-list (index-to-position mesh (list x y z)))
                                           h))))))

(defun cell-calculate-centroid (nodes)
  (let ((centroid (cl-mpm/utils:vector-zeros)))
    (loop for node in nodes
          do
             (cl-mpm/fastmath::fast-.+ centroid (node-position node) centroid))
    (magicl:scale centroid (/ 1d0 (length nodes)))))
(defun make-cell (mesh index h)
  ;;Get local nodes
  (let* ((nodes
           (let ((res nil))
             (array-operations/utilities:nested-loop (x y z) '(2 2 2)
               (push (get-node mesh (mapcar #'+ index (list x y z))) res))
             res))
         (volume (expt h (mesh-nd mesh)))
         ;; (centroid (cell-calculate-centroid nodes))
         (centroid (cl-mpm/fastmath::fast-.+ (cl-mpm/utils:vector-from-list (index-to-position mesh index))
                                         (magicl:scale! (cl-mpm/utils:vector-from-list (list h h h)) 0.5d0)))
         )
    (loop for n in nodes
          do (incf (node-volume-true n) (/ volume (length nodes))))
    (make-instance 'cell
                   :index index
                   :nodes nodes
                   :centroid centroid
                   :volume volume
                   )))
(defun make-cell-2d (mesh index h)
  ;;Get local nodes
  (let* ((nodes
           (let ((res nil))
             (array-operations/utilities:nested-loop (x y z) '(2 2 1)
               (push (get-node mesh (mapcar #'+ index (list x y z))) res))
             res))
         (volume (expt h 2))
         ;; (centroid (cell-calculate-centroid nodes))
         (centroid (cl-mpm/fastmath::fast-.+ (cl-mpm/utils:vector-from-list (index-to-position mesh index))
                                         (magicl:scale! (cl-mpm/utils:vector-from-list (list h h 0d0)) 0.5d0)))
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
  "Make a 3d mesh of specific size"
  (let ((cells
          (make-array
           (mapcar (lambda (x) (- x 1)) size) :initial-contents
           (loop
             for x from 0 below (- (nth 0 size) 1)
             collect
             (loop
               for y from 0 below (- (nth 1 size) 1)
               collect
               (loop
                 for z from 0 below (- (nth 2 size) 1)
                 collect (make-cell mesh (list x y z) h)))))))
    (array-operations/utilities:nested-loop (i j k) (array-dimensions cells)
      (let* ((cell (aref cells i j k)))
        (with-accessors ((neighbours cell-neighbours))
            cell
          (setf neighbours (list))

          (array-operations/utilities:nested-loop (dxx dyy dzz) (list 3 3 3)
            (let ((dx (- dxx 1))
                  (dy (- dyy 1))
                  (dz (- dzz 1)))
              (unless (and (= dx 0)
                           (= dy 0)
                           (= dz 0))
                (let ((di (mapcar #'+ (list i j k) (list dx dy dz))))
                  (when (in-bounds-cell mesh di)
                    (push (apply #'aref cells di) neighbours)))))))))
    cells))
(defun make-cells-2d (mesh size h)
  "Make a 3d mesh of specific size"
  (let ((cells
          (make-array
           (list (- (nth 0 size) 1)
                 (- (nth 1 size) 1)
                 1
                 ) :initial-contents
           (loop
             for x from 0 below (- (nth 0 size) 1)
             collect
             (loop
               for y from 0 below (- (nth 1 size) 1)
               collect
               (loop
                 for z from 0 below 1
                 collect (make-cell-2d mesh (list x y 0) h)))))))
    (array-operations/utilities:nested-loop (i j k) (array-dimensions cells)
      (let* ((cell (aref cells i j k)))
        (with-accessors ((neighbours cell-neighbours))
            cell
          (setf neighbours (list))
          (array-operations/utilities:nested-loop (dxx dyy) (list 3 3)
            (let ((dx (- dxx 1))
                  (dy (- dyy 1))
                  (dz 0))
              (unless (and (= dx 0)
                           (= dy 0)
                           (= dz 0))
                (let ((di (mapcar #'+ (list i j k) (list dx dy dz))))
                  (when (in-bounds-cell mesh di)
                    (push (apply #'aref cells di) neighbours)))))))))
    cells))

(defun make-mesh (size resolution shape-function)
  "Create a 2D mesh and fill it with nodes"
  (let ((nD (length size)))
    (if (= nD 2)
        (setf size (append  size '(0))))
    (let* ((size (mapcar (lambda (x) (coerce x 'double-float)) size))
           (resolution (coerce resolution 'double-float))
           (boundary-order 0)
           (meshcount (loop for d in size collect (+ (floor d resolution) 1 (* boundary-order 2))))
           (nodes '()))
      (let ((mesh (make-instance 'mesh
                                  :nD nD
                                  :mesh-size size
                                  :mesh-count meshcount
                                  :mesh-res resolution
                                  :nodes nodes
                                  :shape-func shape-function
                                  :boundary-order boundary-order
                                  )))
        (setf (mesh-nodes mesh) (make-nodes mesh meshcount resolution))
        (if (= nD 2)
            (setf (mesh-cells mesh) (make-cells-2d mesh meshcount resolution))
            (setf (mesh-cells mesh) (make-cells mesh meshcount resolution)))
        mesh))))

(defun in-bounds-cell (mesh index)
  (destructuring-bind (x y z) index
    (and
     (and (>= x 0) (< x (max 1 (- (nth 0 (mesh-count mesh)) 1))))
     (and (>= y 0) (< y (max 1 (- (nth 1 (mesh-count mesh)) 1))))
     (and (>= z 0) (< z (max 1 (- (nth 2 (mesh-count mesh)) 1))))
     )))

(defun query-boundary-shapes (mesh index)
  (not (member index (mesh-boundary-shapes mesh) :test #'equal)))

(declaim (inline in-bounds-1d))
(defun in-bounds-1d (mesh value dim)
  (declare (fixnum value))
  "Check a single dimension is inside a mesh"
  (and (>= value 0) (< value (the fixnum (nth dim (mesh-count mesh))))))
(declaim
 (inline in-bounds)
 (ftype (function(mesh list) boolean) in-bounds))
(defun in-bounds (mesh pos)
  (declare (optimize (speed 3))
           (list pos))
  "Check a position (list) is inside a mesh"
  (destructuring-bind (x y z) pos
    (declare (fixnum x y z))
    (let ((mc (mesh-count mesh)))
      (declare (fixnum x y)
               (list mc))
      (and (>= x 0) (< x (the fixnum (nth 0 mc)))
           (>= y 0) (< y (the fixnum (nth 1 mc)))
           (>= z 0) (< z (the fixnum (nth 2 mc)))
           )))
  ;; (and (in-bounds-1d mesh (first pos) 0)
  ;;      (in-bounds-1d mesh (second pos) 1)
  ;;      )
  ;; (if (mesh-boundary-shapes mesh)
  ;;     (query-boundary-shapes mesh pos)
  ;;     t)
  ;; (every (lambda (x) x) (loop for d from 0 to (- (mesh-nd mesh) 1)
  ;;                             collect (in-bounds-1d mesh (nth d pos) d)))
  )
(declaim (ftype (function (mesh (simple-array fixnume))) in-bounds-array))
(defun in-bounds-array (mesh pos)
  (declare (type (simple-array fixnum) pos))
  "Check a position (list) is inside a mesh"
  (in-bounds mesh (list (aref pos 0)
                        (aref pos 1)
                        (aref pos 2)))
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
  (list
   (coerce (funcall round-operator (/ (the double-float (magicl:tref pos 0 0)) (the double-float (mesh-resolution mesh)))) 'fixnum)
   (coerce (funcall round-operator (/ (the double-float (magicl:tref pos 1 0)) (the double-float (mesh-resolution mesh)))) 'fixnum)
   (coerce (funcall round-operator (/ (the double-float (magicl:tref pos 2 0)) (the double-float (mesh-resolution mesh)))) 'fixnum)
   )
  ;; (mapcar (lambda (x) (funcall round-operator (/
  ;;                                                 (the double-float (magicl:tref pos x 0))
  ;;                                                 (the double-float (mesh-resolution mesh))))
  ;;                        ) '(0 1 2))
  )

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

(defun get-node-array (mesh pos)
  "Check bounds and get node"
  (policy-cond:policy-if (> safety speed)
                         (if (in-bounds-array mesh pos)
                             (aref (mesh-nodes mesh) (aref pos 0) (aref pos 1))
                             (error (format nil "Access grid out of bounds at: ~a" pos)))
                         (aref (mesh-nodes mesh) (aref pos 0) (aref pos 1))))

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
;;         (setf vel (cl-mpm/fastmath::fast-.+ vel (magicl:scale node-vel svp)))))))
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
               (pressure pressure)
               (force force)
               (int-force internal-force)
               (ext-force external-force)
               (buoyancy-force buoyancy-force)
               )
                node
    (setf active nil)
    (setf boundary nil)
    (setf mass 0d0)
    (setf volume 0d0)
    (setf p-wave 0d0)
    (setf damage 0d0)
    (setf svp-sum 0d0)
    (setf j 0d0)
    (setf j-inc 0d0)
    (setf boundary-scalar 0d0)
    (setf pressure 0d0)
    (magicl:scale! vel 0d0)
    (magicl:scale! acc 0d0)
    (magicl:scale! force 0d0)
    (magicl:scale! int-force 0d0)
    (magicl:scale! ext-force 0d0)
    (magicl:scale! buoyancy-force 0d0)
    ))

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
                                (/ h (* 2 (sqrt 3)))))))))

(defun gauss-weights (n h)
  (cond
    ((= n 1)
     (list 1d0)
     )
    ((= n 2)
     (loop for x from -1 to 1 by 2
           append
           (loop for y from -1 to 1 by 2
                 collect
                 0.25d0)))))

(defun newton-points (n h)
  (declare (fixnum n))
  (let ((spaceing (/ h (+ (coerce n 'double-float) 0d0))))
    (declare (double-float h spaceing))
    (cond
      ((= n 1)
       (list (magicl:zeros '(2 1))))
      ((> n 1)
       (let ((start-point
               (magicl:.+
                (magicl:scale! (cl-mpm/utils:vector-from-list (list -0.5d0 -0.5d0)) h)
                (magicl:scale! (cl-mpm/utils:vector-from-list (list 0.5d0 0.5d0)) spaceing)
                )
               ))
         ;; (break)
         (loop for x from 0 below n
               append
               (loop for y from 0 below n
                     collect
                     (magicl:.+
                      (magicl:scale!
                       (cl-mpm/utils:vector-from-list (mapcar
                                                       (lambda (x) (coerce x 'double-float))
                                                       (list x y))) spaceing
                       )
                      start-point))))))))

(defun newton-weights (n h)
  (cond
    ((= n 1)
     (list 1d0)
     )
    ((> n 1)
     (loop for x from 1 to n
           append
           (loop for y from 1 to n
                 collect
                 (/ 1d0 (* n n)))))))

(defun cell-quadrature-iterate-over-neighbours (mesh cell gp func)
  (declare (function func))
  (let ((h (cl-mpm/mesh:mesh-resolution mesh)))
    (with-accessors ((nodes cell-nodes)
                     (centroid cell-centroid)
                     (volume cell-volume))
        cell
      (loop for point in (newton-points gp h)
            for gweight in (newton-weights gp h)
            do
        (let ((quad (cl-mpm/fastmath::fast-.+ centroid point))
              (volume-ratio (/ 1 (expt gp 2))))
          (loop for node in nodes
                do
                   (with-accessors ((n-pos node-position))
                       node
                     (let* ((dist-vec (magicl:.- quad n-pos))
                            (dist (list (magicl:tref dist-vec 0 0) (magicl:tref dist-vec 1 0)))
                            (weights (mapcar (lambda (x) (cl-mpm/shape-function::shape-linear x h)) dist))
                            (weight (reduce #'* weights))
                            ;; (weight 1d0)
                            ;; (weights ())
                            (grads (mapcar (lambda (d w) (* (cl-mpm/shape-function::shape-linear-dsvp d h) w))
                                           dist (nreverse weights)))
                            )
                       (funcall func
                                mesh
                                cell
                                quad
                                (* volume gweight)
                                node
                                weight
                                grads)))))))))


(defun cell-iterate-over-neighbours (mesh cell func)
  (if (= (mesh-nd mesh) 2)
      (cell-iterate-over-neighbours-2d mesh cell func)
      (cell-iterate-over-neighbours-3d mesh cell func)
      )
  )
(defun cell-iterate-over-neighbours-2d (mesh cell func)
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
                        ;; (grads (mapcar (lambda (d w) (* (cl-mpm/shape-function::shape-linear-dsvp d h) w))
                        ;;                dist (nreverse weights)))
                        (lin-grads
                          (mapcar (lambda (d) (cl-mpm/shape-function::shape-linear-dsvp d h))
                                           dist))
                        (grads (cl-mpm/shape-function::grads-2d weights lin-grads))
                        )
                   (when (< 0d0 weight)
                     (funcall func
                              mesh
                              cell
                              centroid
                              volume
                              node
                              weight
                              (list (nth 0 grads)
                                    (nth 1 grads)
                                    0d0
                                    )))))))))
(defun cell-iterate-over-neighbours-3d (mesh cell func)
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
                        (dist (list (magicl:tref dist-vec 0 0) (magicl:tref dist-vec 1 0) (magicl:tref dist-vec 2 0)))
                        (weights (mapcar (lambda (x) (cl-mpm/shape-function::shape-linear x h)) dist))
                        (weight (reduce #'* weights))
                        ;; (grads (mapcar (lambda (d w) (* (cl-mpm/shape-function::shape-linear-dsvp d h) w))
                        ;;                dist (nreverse weights)))
                        (lin-grads
                          (mapcar (lambda (d) (cl-mpm/shape-function::shape-linear-dsvp d h))
                                           dist))
                        (grads (cl-mpm/shape-function::grads-3d weights lin-grads))
                        )
                   (when (< 0d0 weight)
                     (funcall func
                              mesh
                              cell
                              centroid
                              volume
                              node
                              weight
                              grads))))))))

;; (defun cell-iterate-over-neighbours-point (mesh cell pos func)
;;   (declare (function func))
;;     (with-accessors ((nodes cell-nodes)
;;                      (centroid cell-centroid)
;;                      (volume cell-volume))
;;         cell
;;       (cl-mpm::iterate-over-neighbours-point-linear
;;        mesh pos
;;        (lambda (mesh node weight grads)
;;          (funcall func mesh cell pos volume node weight grads)))))


;;; Printing methods

(defmethod print-object ((obj node) stream)
  (print-unreadable-object (obj stream :type t)
    (with-accessors ((index node-index)
                     (mass node-mass)
                     (vel node-velocity)
                     )
        obj
      (format stream "index: ~a, mass: ~a, vel: ~a" index mass vel))))
