(defpackage :cl-mpm/output
  (:use :cl
   :cl-mpm/utils)
  (:export
   #:save-vtk-mesh
   #:save-vtk
   #:save-vtk-nodes
   #:save-vtk-cells
   #:save-csv
   #:save-simulation-parameters
   ))

(in-package :cl-mpm/output)
(declaim (optimize (debug 0) (safety 0) (speed 3)))
;; (declaim (optimize (debug 3) (safety 3) (speed 0)))

(defun save-sim (sim)
  "Save a simulation such that it can be restarted from that point")
(defun load-sim (sim))

(defun sample-line (sim start end steps get-value)
	"Sample over a line from start to end with fixed steps count, get-value function gets given a sampling mp and a mesh and should return a scalar"
    (let (  (sample-mp (cl-mpm/particle:make-particle 2 :pos start))
            (direction (magicl:from-list (mapcar (lambda (s e) (/ (- e s) steps)) start end) '(2 1))))
      (loop :repeat (+ steps 1) collect 
              (prog1
                (funcall get-value (cl-mpm:sim-mesh sim) sample-mp)
                (setf (cl-mpm/particle:mp-position sample-mp) (cl-mpm/fastmaths::fast-.+ (cl-mpm/particle:mp-position sample-mp) direction))))))

(defun sample-point (sim point get-value)
    (let ((sample-mp (cl-mpm/particle:make-particle 2 :pos point)))
        (funcall get-value (cl-mpm:sim-mesh sim) sample-mp)))

(defmacro scalar-average (accessor)
    `(lambda (mesh mp)
        (let ((av 0))
          (progn
            (cl-mpm::iterate-over-neighbours mesh mp
                 (lambda (node svp dsvp &rest rest)
                   (setf av (+ av (* svp (,accessor node))))))
            av))))

(defmacro matrix-average (accessor shape)
    `(lambda (mesh mp)
        (let ((av (magicl:zeros ,shape)))
          (progn
            (cl-mpm::iterate-over-neighbours
             mesh mp
             (lambda (node svp dsvp  &rest rest)
               (setf av (cl-mpm/fastmaths::fast-.+ av (magicl:scale (,accessor node) svp)))))
            av))))

(defun sample-line-mass (sim start end steps)
    (sample-line sim start end steps
        (scalar-average cl-mpm/mesh:node-mass)))

(defun sample-line-velocity (sim start end steps)
    (sample-line sim start end steps
        (matrix-average cl-mpm/mesh:node-velocity '(2 1))))

(defun sample-point-mass (sim point)
    (sample-point sim point (scalar-average cl-mpm/mesh:node-mass)))

(defun sample-point-velocity (sim point)
    (sample-point sim point
        (matrix-average cl-mpm/mesh:node-velocity '(2 1))))

(defun format-scalar (stream stream-bin name id mps accessor)
  (format stream "SCALARS ~a double ~d~%" name 1)
  (format stream "LOOKUP_TABLE default~%")
  (force-output stream)
  (loop for mp across mps
        do (write-binary-float (coerce (funcall accessor mp) 'double-float) stream-bin)
           ;; (format stream "~E ~%" (coerce (min 1d30 (funcall accessor mp)) 'single-float))
        )
  (force-output stream-bin)
  (format stream "~%")
  )

(defmacro save-parameter (name accessor)
  (let ((mps (intern (symbol-name 'mps)))
        (fs (intern (symbol-name 'fs)))
        (fs-bin (intern (symbol-name 'fs-bin)))
        (mp (intern (symbol-name 'mp)))
        (id (intern (symbol-name 'id)))
        )
    `(progn
       (format-scalar ,fs ,fs-bin ,name ,id ,mps (lambda (,mp) ,accessor))
     (incf ,id))))


(defun format-vector (stream name id mps accessor)
  (format stream "VECTORS ~a FLOAT~%" name)
  (format stream "LOOKUP_TABLE default~%")
    (cl-mpm::iterate-over-mps-serial
     mps
     (lambda (mp)
       (let ((v (funcall accessor mp)))
         (format stream "~E ~E ~E~%"
                 (coerce (varef v 0) 'single-float)
                 (coerce (varef v 1) 'single-float)
                 (coerce (varef v 2) 'single-float)))))
  (format stream "~%"))

(defmacro save-parameter-vector (name accessor)
  (let ((mps (intern (symbol-name 'mps)))
        (fs (intern (symbol-name 'fs)))
        (mp (intern (symbol-name 'mp)))
        (id (intern (symbol-name 'id)))
        )
    `(progn
       (format-vector ,fs ,name ,id ,mps (lambda (,mp) ,accessor))
       (incf ,id))))

;; (defmacro save-parameter (name accessor)
;;   (let ((mps (intern (symbol-name 'mps)))
;;         (fs (intern (symbol-name 'fs)))
;;         (mp (intern (symbol-name 'mp)))
;;         (id (intern (symbol-name 'id)))
;;         )
;;   `(progn
;;      (format-scalar ,fs ,name ,id ,mps (lambda (,mp) ,accessor))
;;      (incf ,id))))


(defun save-simulation-parameters (filename sim &rest args)
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (str:to-file
     filename
     (jonathan:to-json
      (append
       (list
        :resolution (cl-mpm/mesh::mesh-resolution mesh)
        :domain-size (cl-mpm/mesh::mesh-mesh-size mesh))
       (cl-mpm::sim-settings sim)
       (apply #'append args))))))
(defun save-vtk-mesh (filename sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)) sim
    (with-open-file (fs filename :direction :output :if-exists :supersede)
      (format fs "# vtk DataFile Version 2.0~%")
      (format fs "Lisp generated vtk file, SJVS~%")
      (format fs "ASCII~%")
      (format fs "DATASET UNSTRUCTURED_GRID~%")
      (with-accessors ((nodes cl-mpm/mesh::mesh-nodes)
                       (size cl-mpm/mesh::mesh-count)
                       (h cl-mpm/mesh::mesh-resolution)
                       (border cl-mpm/mesh::mesh-boundary-order)) mesh
        ;; (format fs "POINTS ~d double~%" (floor (apply #'* size)))
        (format fs "POINTS ~d double~%" (array-total-size nodes))
        (array-operations/utilities:nested-loop
         (i j k) size
         (format fs "~E ~E ~E ~%"
                 (coerce (magicl:tref (cl-mpm/mesh::node-position (aref nodes i j k)) 0 0) 'single-float)
                 (coerce (magicl:tref (cl-mpm/mesh::node-position (aref nodes i j k)) 1 0) 'single-float)
                 (coerce (magicl:tref (cl-mpm/mesh::node-position (aref nodes i j k)) 2 0) 'single-float)))
        (let ((nels (floor (reduce #'* (mapcar (lambda (x) (max 1 (- x 1))) size))))
              (nen 4)
              (node-order (list 1 2 4 3 5 6 8 7))
              )
          (format fs "CELLS ~D ~D~%"
                  nels
                  (floor (* (+ nen 1) nels)))
          (flet ((id (x y z) (floor (+ z (* (+ y (* x (second size))) (third size)))))
                 ;; (local-id (x y z) (floor (+ z (* (+ y (* x 2)) 2))))
                 )
            (array-operations/utilities:nested-loop
             (i j k) (mapcar (lambda (x) (- x 1)) size)
             (let ((node-map '()))
               (array-operations/utilities:nested-loop
                (dx dy dz) '(2 2 2)
                (push (id (+ i dx)
                          (+ j dy)
                          (+ k dz)) node-map)
                )
               (setf node-map (reverse node-map))
               (format fs "~D " nen)
               (loop for lid in node-order
                     do (format fs "~D " (nth (- lid 1) node-map)))
               (format fs "~%")))
            ;; (array-operations/utilities:nested-loop
            ;;  (i j k) (mapcar (lambda (x) (- x 1)) size)
            ;;  (format fs "~D " nen)
            ;;  (array-operations/utilities:nested-loop
            ;;   (dx dy dz) '(2 2 2)
            ;;   (format fs "~D " (id (+ i dx)
            ;;                        (+ j dy)
            ;;                        (+ k dz)))
            ;;   )
            ;;  (format fs "~%")
            ;;  )
            (when (= 2 (cl-mpm/mesh:mesh-nd mesh))
              (loop for x from 0 to (- (first size) 2)
                    do
                       (loop for y from 0 to (- (second size) 2)
                             do
                                (format fs "~D ~D ~D ~D ~D~%"
                                        4
                                        (id x y 0)
                                        (id x (+ y 1) 0)
                                        (id (+ x 1) y 0)
                                        (id (+ x 1) (+ y 1) 0)
                                        )))
              (format fs "CELL_TYPES ~d~%" nels)
              (loop repeat nels
                    do (format fs "~D~%" 8)))
            (when (= 3 (cl-mpm/mesh:mesh-nd mesh))
              ;; (loop for x from 0 to (- (first size) 2)
              ;;       do
              ;;          (loop for y from 0 to (- (second size) 2)
              ;;                do
              ;;                   (loop for z from 0 to (- (third size) 2)
              ;;                         do
              ;;                            (format fs "~D ~D ~D ~D ~D~%"
              ;;                                    4
              ;;                                    (id x y)
              ;;                                    (id x (+ y 1))
              ;;                                    (id (+ x 1) y)
              ;;                                    (id (+ x 1) (+ y 1))
              ;;                                    ))))
              ;; (format fs "CELL_TYPES ~d~%" nels)
              ;; (loop repeat nels
              ;;       do (format fs "~D~%" 12))
              ))
          )))))

(defun save-csv (filename sim)
  (with-accessors ((mps cl-mpm:sim-mps)) sim
    (with-open-file (fs filename :direction :output :if-exists :supersede)
      (format fs "coord_x,coord_y,stress_xx,stress_yy,tau_xy,velocity_x,velocity_y,stress_1,eps,damage,lx,ly~%")
      (loop for mp across mps
            do (format fs "~E, ~E, ~F, ~F, ~F, ~F, ~F, ~F, ~F, ~F, ~F, ~F ~%"
                       (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 1 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle:mp-stress mp) 2 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle:mp-velocity mp) 0 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle:mp-velocity mp) 1 0) 'single-float)
                       (coerce (multiple-value-bind (l v)
                                   (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))
                                 (loop for sii in l maximize sii)) 'single-float)
                       (coerce
                        (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))
                          (- (apply #'max l) (cl-mpm/particle::mp-pressure mp))) 'single-float)
                       (coerce (if (slot-exists-p mp 'cl-mpm/particle::damage)
                                   (cl-mpm/particle:mp-damage mp)
                                   0d0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle::mp-domain-size mp) 0 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0) 'single-float)
                       ))
      )))

(defun format-scalar-node (stream stream-bin name id nodes accessor)
  (force-output stream-bin)
  (format stream "SCALARS ~a double ~d~%" name 1)
  (format stream "LOOKUP_TABLE default~%")
  (force-output stream)
  (loop for i from 0 below (array-total-size nodes)
        do
           (let ((node (row-major-aref nodes i)))
             (when node
               (write-binary-float (coerce (funcall accessor node) 'double-float) stream-bin)
               ;; (format stream "~E ~%" (coerce (min 1d30 (funcall accessor node)) 'single-float))
               )))
  (force-output stream-bin)
  (format stream "~%"))

(defmacro save-parameter-nodes (name accessor)
  `(progn
     (format-scalar-node fs fs-bin ,name id nodes (lambda (node) ,accessor))
     (incf id)))

(defun format-vector-node (stream name id nodes accessor)
  (format stream "VECTORS ~a FLOAT~%" name)
  ;; (format stream "LOOKUP_TABLE default~%")
  (loop for i from 0 below (array-total-size nodes)
        do
           (let ((node (row-major-aref nodes i)))
             (when node
               (let ((v (funcall accessor node)))
                 (format stream "~E ~E ~E~%"
                         (coerce (varef v 0) 'single-float)
                         (coerce (varef v 1) 'single-float)
                         (coerce (varef v 2) 'single-float))))))
  (format stream "~%"))

(defmacro save-vector-nodes (name accessor)
  `(progn
     (format-vector-node fs ,name id nodes (lambda (node) ,accessor))
     (incf id)))

(defgeneric save-vtk-mesh-nodes (filename mesh sim))
(defmethod save-vtk-mesh-nodes (filename mesh sim)
  (with-accessors ((full-nodes cl-mpm/mesh::mesh-nodes))
      mesh
    (with-open-file (fs filename :direction :output :if-exists :supersede)
      (format fs "# vtk DataFile Version 2.0~%")
      (format fs "Lisp generated vtk file, SJVS~%")
      (format fs "BINARY~%")
      (format fs "DATASET UNSTRUCTURED_GRID~%"))
    (with-open-file (fs filename :direction :output :if-exists :append)
      (with-open-file (fs-bin filename :direction :output :if-exists :append :element-type '(unsigned-byte 8))
        (let* ((node-count (array-total-size full-nodes))
               (nodes (remove-if-not 'cl-mpm/mesh:node-active
                                     (make-array node-count :displaced-to full-nodes))))
          (when (= (length nodes) 0)
            (setf nodes (make-array node-count :displaced-to full-nodes)))
          (let ((node-count (array-total-size nodes)))
            (format fs "POINTS ~d double~%" (array-total-size nodes))
            (force-output fs)
            (loop for i from 0 below (array-total-size nodes)
                  do
                     (let ((node (row-major-aref nodes i)))
                       (when node
                         (let ((pos (cl-mpm/fastmaths:fast-.+
                                     (cl-mpm/mesh::node-position node)
                                     (cl-mpm/mesh::node-displacment node))))
                           (write-binary-float (cl-mpm/utils:varef pos 0) fs-bin)
                           (write-binary-float (cl-mpm/utils:varef pos 1) fs-bin)
                           (write-binary-float (cl-mpm/utils:varef pos 2) fs-bin)))))
            (force-output fs-bin)
            (format fs "~%")
            (let ((id 1)
                  (nd (cl-mpm/mesh:mesh-nd mesh)))
              (declare (special id))
              (format fs "POINT_DATA ~d~%" (array-total-size nodes))

              (force-output fs)
              (dolist (f (cl-mpm::sim-output-list-nodes sim))
                (destructuring-bind (type name accessor) f
                  (case type
                    (:BOOL
                     (cl-mpm/output::save-parameter-nodes name (if (funcall accessor node) 1d0 0d0)))
                    (:SCALAR
                     (cl-mpm/output::save-parameter-nodes name (funcall accessor node)))
                    (:VECTOR
                     (cl-mpm/output::save-parameter-nodes (format nil "~A_mag" name)
                                                          (cl-mpm/fastmaths::mag (funcall accessor node)))
                     (cl-mpm/output::save-parameter-nodes (format nil "~A_x" name) (varef (funcall accessor node) 0))
                     (cl-mpm/output::save-parameter-nodes (format nil "~A_y" name) (varef (funcall accessor node) 1))
                     (when (= nd 3)
                       (cl-mpm/output::save-parameter-nodes (format nil "~A_z" name) (varef (funcall accessor node) 2)))))))
              #|
              (save-parameter-nodes "active" (if (cl-mpm/mesh::node-active node) 1 0))

              (save-parameter-nodes "agg" (if (cl-mpm/mesh::node-agg node) 1 0))
              (save-parameter-nodes "agg-int" (if (cl-mpm/mesh::node-interior node) 1 0))
              (save-parameter-nodes "agg-ext" (if (and (cl-mpm/mesh::node-agg node) (not (cl-mpm/mesh::node-interior node))) 1 0))
              (save-parameter-nodes "agg-element-x" (if (cl-mpm/mesh::node-agg-interior-cell node)
              (nth 0 (cl-mpm/mesh::cell-index (cl-mpm/mesh::node-agg-interior-cell node))) -1))
              (save-parameter-nodes "agg-element-y" (if (cl-mpm/mesh::node-agg-interior-cell node)
              (nth 1 (cl-mpm/mesh::cell-index (cl-mpm/mesh::node-agg-interior-cell node))) -1))


              (save-parameter-nodes "mass" (cl-mpm/mesh:node-mass node))
              (save-parameter-nodes "density" (if (> (cl-mpm/mesh::node-volume node) 0d0)
              (/ (cl-mpm/mesh:node-mass node) (cl-mpm/mesh::node-volume node))
              0d0))
              (save-parameter-nodes "mass-inv" (if (> (cl-mpm/mesh:node-mass node) 0d0)
              (/ 1d0 (cl-mpm/mesh:node-mass node))
              0d0))

              (save-parameter-nodes "vel_norm" (cl-mpm/fastmaths::mag (cl-mpm/mesh:node-velocity node)))
              (save-parameter-nodes "vel_x" (magicl:tref (cl-mpm/mesh:node-velocity node) 0 0))
              (save-parameter-nodes "vel_y" (magicl:tref (cl-mpm/mesh:node-velocity node) 1 0))
              (save-parameter-nodes "vel_z" (magicl:tref (cl-mpm/mesh:node-velocity node) 2 0))

              (save-parameter-nodes "disp_x" (magicl:tref (cl-mpm/mesh::node-displacment node) 0 0))
              (save-parameter-nodes "disp_y" (magicl:tref (cl-mpm/mesh::node-displacment node) 1 0))
              (save-parameter-nodes "disp_z" (magicl:tref (cl-mpm/mesh::node-displacment node) 2 0))

              (save-parameter-nodes "bc_x" (if (cl-mpm/mesh::node-bcs node) (varef (cl-mpm/mesh::node-bcs node) 0) 0))
              (save-parameter-nodes "bc_y" (if (cl-mpm/mesh::node-bcs node) (varef (cl-mpm/mesh::node-bcs node) 1) 0))
              (save-parameter-nodes "bc_z" (if (cl-mpm/mesh::node-bcs node) (varef (cl-mpm/mesh::node-bcs node) 2) 0))


              (save-vector-nodes "force" (cl-mpm/mesh:node-force node))
              (save-vector-nodes "force-ext" (cl-mpm/mesh::node-external-force node))
              (save-vector-nodes "force-int" (cl-mpm/mesh::node-internal-force node))
              (save-vector-nodes "force-rct" (cl-mpm/mesh::node-reaction-force node))
              (save-vector-nodes "force-buoyancy" (cl-mpm/mesh::node-buoyancy-force node))
              (save-vector-nodes "force-ghost" (cl-mpm/mesh::node-ghost-force node))
              (save-vector-nodes "force-damping" (cl-mpm/mesh::node-damping-force node))
              (save-vector-nodes "acc" (cl-mpm/mesh:node-acceleration node))

              (save-parameter-nodes "j-inc" (if (= (cl-mpm/mesh::node-jacobian-inc node) 0d0)
              1d0
              (cl-mpm/mesh::node-jacobian-inc node)
              ))
              (save-parameter-nodes "buoyancy_node" (if
              (cl-mpm/mesh::node-boundary-node node) 1 0))
              (save-parameter-nodes "buoyancy-scalar" (cl-mpm/mesh::node-boundary-scalar node))
              (save-vector-nodes "boundary-vec" (cl-mpm/mesh::node-boundary-vec node))
              (save-parameter-nodes "pressure" (cl-mpm/mesh::node-pressure node))
              (save-parameter-nodes "local-list-size" (length (cl-mpm/mesh::node-local-list node)))
              (save-parameter-nodes "energy"
              (*
              (cl-mpm/mesh::node-mass node)
              (cl-mpm/mesh::node-mass node)
              (cl-mpm/fastmaths::mag-squared (cl-mpm/mesh::node-velocity node))))

              (save-parameter-nodes "oobf" (cl-mpm/mesh::node-oobf node))
              (save-parameter-nodes "volume" (cl-mpm/mesh::node-volume node))
              (save-parameter-nodes "volume-true" (cl-mpm/mesh::node-volume-true node))
              (save-parameter-nodes "damage"
              (if (slot-exists-p node 'cl-mpm/mesh::damage)
              (cl-mpm/mesh::node-damage node)
              0d0))
              |#
              )))))))

(defgeneric save-vtk-nodes (filename sim))
(defmethod save-vtk-nodes (filename sim)
  (save-vtk-mesh-nodes filename (cl-mpm:sim-mesh sim) sim))


(defgeneric save-vtk (filename sim)
  (:documentation "Save vtk depending on type")
   (:method (f s)))

(defmacro optional-slot-access (slot mp)
  `(if (slot-exists-p ,mp ,slot)
       (slot-value ,mp ,slot)
       0d0))

(defun float-to-binary (data)
  (reverse
   (cl-intbytes:int64->octets (ieee-floats:encode-float64 (coerce data 'double-float))))
  ;; (cl-intbytes:int64->octets (ieee-floats:encode-float64 data))
  )

(defun write-binary-float (data stream)
  (loop for v across (float-to-binary data)
        do (write-byte v stream)))

(defmethod initialize-instance :after ((sim cl-mpm::mpm-sim) &key)
  (setf
   (cl-mpm::sim-output-list sim)
   (append
    (cl-mpm::sim-output-list sim)
    (list
     (list :BOOL "fric-contact" #'cl-mpm/particle::mp-penalty-contact-step)
     (list :SCALAR "i1" (lambda (mp) (cl-mpm/utils::trace-voigt (cl-mpm/particle:mp-stress mp))))
     (list :SCALAR "unique-id" #'cl-mpm/particle::mp-unique-index)
     (list :SCALAR "mass" #'cl-mpm/particle::mp-mass)
     (list :SCALAR "volume" #'cl-mpm/particle::mp-volume)
     (list :SCALAR "index" #'cl-mpm/particle::mp-index)
     (list :SCALAR "mpi-index" #'cl-mpm/particle::mp-mpi-index)
     (list :SCALAR "split-depth" #'cl-mpm/particle::mp-split-depth)
     (list :SCALAR "fric-normal" #'cl-mpm/particle::mp-penalty-normal-force)
     (list :SCALAR "p-wave-modulus" #'cl-mpm/particle::mp-p-modulus)
     (list :SCALAR "yield-function" (lambda (mp) (if (slot-exists-p mp 'cl-mpm/particle::yield-func) (cl-mpm/particle::mp-yield-func mp) 0d0)))
     (list :SCALAR "plastic-strain" (lambda (mp) (if (slot-exists-p mp 'cl-mpm/particle::strain-plastic-vm) (cl-mpm/particle::mp-strain-plastic-vm mp) 0d0)))
     (list :VECTOR "vel" #'cl-mpm/particle::mp-velocity)
     (list :VECTOR "acc" #'cl-mpm/particle::mp-acceleration)
     (list :VECTOR "disp" #'cl-mpm/particle::mp-displacement)
     (list :VECTOR "size" #'cl-mpm/particle::mp-domain-size)
     (list :VECTOR "fric-force" #'cl-mpm/particle::mp-penalty-frictional-force)
     (list :VOIGT "sig" #'cl-mpm/particle::mp-stress)
     (list :VOIGT "eps" #'cl-mpm/particle::mp-strain)
     ;; (list :MATRIX "size" (lambda (mp) (cl-mpm/particle::mp-true-domain mp)))
     )))
  (setf
   (cl-mpm::sim-output-list-nodes sim)
   (append
    (cl-mpm::sim-output-list-nodes sim)
    (list
     (list :BOOL "active" #'cl-mpm/mesh::node-active)
     (list :BOOL "agg" #'cl-mpm/mesh::node-agg)
     (list :BOOL "agg-int" #'cl-mpm/mesh::node-interior)
     (list :BOOL "agg-ext" (lambda (n) (and
                                        (cl-mpm/mesh::node-agg n)
                                        (not (cl-mpm/mesh::node-interior n)))))
     (list :BOOL "boundary" #'cl-mpm/mesh::node-boundary-node)
     (list :SCALAR "boundary-scalar" #'cl-mpm/mesh::node-boundary-scalar)
     (list :SCALAR "mass" #'cl-mpm/mesh::node-mass)
     (list :VECTOR "vel" #'cl-mpm/mesh::node-velocity)
     (list :VECTOR "disp" #'cl-mpm/mesh::node-displacment)
     (list :VECTOR "force" #'cl-mpm/mesh::node-force)
     (list :VECTOR "force-int" #'cl-mpm/mesh::node-internal-force)
     (list :VECTOR "force-ext" #'cl-mpm/mesh::node-external-force)
     (list :VECTOR "force-damping" #'cl-mpm/mesh::node-damping-force)
     (list :VECTOR "force-buoyancy" #'cl-mpm/mesh::node-buoyancy-force)
     ))))


(defmethod save-vtk (filename (sim cl-mpm::mpm-sim))
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)) sim
    (with-open-file (fs filename :direction :output :if-exists :supersede)
      (format fs "# vtk DataFile Version 2.0~%")
      (format fs "Lisp generated vtk file, SJVS~%")
      ;; (format fs "ASCII~%")
      (format fs "BINARY~%")
      (format fs "DATASET UNSTRUCTURED_GRID~%")
      (format fs "POINTS ~d double~%" (length mps)))
    (with-open-file (fs filename :direction :output :if-exists :append)
      (with-open-file (fs-bin filename :direction :output :if-exists :append :element-type '(unsigned-byte 8))
        (force-output fs)
        (loop for mp across mps
              do
                 (let ((pos (cl-mpm/particle::mp-position-trial mp)))
                   (write-binary-float (cl-mpm/utils:varef pos 0) fs-bin)
                   (write-binary-float (cl-mpm/utils:varef pos 1) fs-bin)
                   (write-binary-float (cl-mpm/utils:varef pos 2) fs-bin)
                   ;; (format fs "~E ~E ~E ~%"
                   ;;         (coerce (cl-mpm/utils:varef pos 0) 'single-float)
                   ;;         (coerce (cl-mpm/utils:varef pos 1) 'single-float)
                   ;;         (coerce (cl-mpm/utils:varef pos 2) 'single-float))
                   ))
        (force-output fs-bin)
        (format fs "~%")
        (let ((id 1)
              (nd (cl-mpm/mesh:mesh-nd mesh)))
          (declare (special id))
          (format fs "POINT_DATA ~d~%" (length mps))
          (dolist (f (cl-mpm::sim-output-list sim))
            (destructuring-bind (type name accessor) f
              (case type
                (:BOOL
                 (cl-mpm/output::save-parameter name (if (funcall accessor mp) 1d0 0d0)))
                (:SCALAR
                 (cl-mpm/output::save-parameter name (funcall accessor mp)))
                (:VECTOR
                 (cl-mpm/output::save-parameter (format nil "~A_mag" name) (cl-mpm/fastmaths::mag (funcall accessor mp)))
                 (cl-mpm/output::save-parameter (format nil "~A_x" name) (varef (funcall accessor mp) 0))
                 (cl-mpm/output::save-parameter (format nil "~A_y" name) (varef (funcall accessor mp) 1))
                 (when (= nd 3)
                   (cl-mpm/output::save-parameter (format nil "~A_z" name) (varef (funcall accessor mp) 2))))
                (:VOIGT
                 (cl-mpm/output::save-parameter (format nil "~A_xx" name) (varef (funcall accessor mp) 0))
                 (cl-mpm/output::save-parameter (format nil "~A_yy" name) (varef (funcall accessor mp) 1))
                 (when (= nd 3)
                   (cl-mpm/output::save-parameter (format nil "~A_zz" name) (varef (funcall accessor mp) 2))
                   (cl-mpm/output::save-parameter (format nil "~A_yz" name) (varef (funcall accessor mp) 3))
                   (cl-mpm/output::save-parameter (format nil "~A_xz" name) (varef (funcall accessor mp) 4)))
                 (cl-mpm/output::save-parameter (format nil "~A_xy" name) (varef (funcall accessor mp) 5)))
                (:MATRIX
                 (cl-mpm/output::save-parameter (format nil "~A_xx" name) (mtref (funcall accessor mp) 0 0))
                 (cl-mpm/output::save-parameter (format nil "~A_yy" name) (mtref (funcall accessor mp) 1 1))
                 (cl-mpm/output::save-parameter (format nil "~A_xy" name) (mtref (funcall accessor mp) 0 1))
                 (cl-mpm/output::save-parameter (format nil "~A_yx" name) (mtref (funcall accessor mp) 1 0))))))
          ;; (save-parameter "unique-id" (cl-mpm/particle::mp-unique-index mp))
          ;; (save-parameter "mass" (cl-mpm/particle:mp-mass mp))
          ;; (save-parameter "density" (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)))
          ;; (save-parameter "volume" (cl-mpm/particle:mp-volume mp))

          ;; (save-parameter "index" (cl-mpm/particle::mp-index mp))
          ;; (save-parameter "mpi-domain" (cl-mpm/particle::mp-mpi-index mp))
          ;; (save-parameter "split-depth" (cl-mpm/particle::mp-split-depth mp))
          ;; (cl-mpm/output::save-parameter "j2"
          ;;                                (sqrt (cl-mpm/fastmaths::voigt-j2 (cl-mpm/utils:deviatoric-voigt (cl-mpm/particle::mp-stress mp)))))

          ;; ;; (save-parameter-vector "vel" cl-mpm/particle::mp-velocity)
          ;; (save-parameter "vel_x" (magicl:tref (cl-mpm/particle:mp-velocity mp) 0 0))
          ;; (save-parameter "vel_y" (magicl:tref (cl-mpm/particle:mp-velocity mp) 1 0))
          ;; (save-parameter "p-wave-modulus" (cl-mpm/particle::mp-p-modulus mp))
          ;; ;; (save-parameter "vel_z" (magicl:tref (cl-mpm/particle:mp-velocity mp) 2 0))

          ;; ;; (save-parameter "acc_x" (magicl:tref (cl-mpm/particle::mp-acceleration mp) 0 0))
          ;; ;; (save-parameter "acc_y" (magicl:tref (cl-mpm/particle::mp-acceleration mp) 1 0))

          ;; (save-parameter "disp_x" (magicl:tref (cl-mpm/fastmaths:fast-.+
          ;;                                        (cl-mpm/particle::mp-displacement-increment mp)
          ;;                                        (cl-mpm/particle::mp-displacement mp)) 0 0))
          ;; (save-parameter "disp_y" (magicl:tref (cl-mpm/fastmaths:fast-.+
          ;;                                        (cl-mpm/particle::mp-displacement-increment mp)
          ;;                                        (cl-mpm/particle::mp-displacement mp)) 1 0))
          ;; (save-parameter "disp_z" (magicl:tref (cl-mpm/fastmaths:fast-.+
          ;;                                        (cl-mpm/particle::mp-displacement-increment mp)
          ;;                                        (cl-mpm/particle::mp-displacement mp)) 2 0))

          ;; (cl-mpm/output::save-parameter "J" (cl-mpm/fastmaths:det-3x3 (cl-mpm/particle::mp-deformation-gradient mp)))
          ;; (cl-mpm/output::save-parameter "sig_xx" (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0))
          ;; (cl-mpm/output::save-parameter "sig_yy" (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0))
          ;; (cl-mpm/output::save-parameter "sig_xy" (magicl:tref (cl-mpm/particle:mp-stress mp) 5 0))
          ;; (cl-mpm/output::save-parameter "i1" (cl-mpm/utils::trace-voigt (cl-mpm/particle:mp-stress mp)))
          ;; (when (= nd 3)
          ;;   (cl-mpm/output::save-parameter "sig_zz" (magicl:tref (cl-mpm/particle:mp-stress mp) 2 0))
          ;;   (cl-mpm/output::save-parameter "sig_yz" (magicl:tref (cl-mpm/particle:mp-stress mp) 3 0))
          ;;   (cl-mpm/output::save-parameter "sig_zx" (magicl:tref (cl-mpm/particle:mp-stress mp) 4 0)))

          ;; (cl-mpm/output::save-parameter "eps_xx" (magicl:tref (cl-mpm/particle:mp-strain mp) 0 0))
          ;; (cl-mpm/output::save-parameter "eps_yy" (magicl:tref (cl-mpm/particle:mp-strain mp) 1 0))
          ;; (cl-mpm/output::save-parameter "eps_xy" (magicl:tref (cl-mpm/particle:mp-strain mp) 5 0))
          ;; (save-parameter "eps_1"
          ;;                 (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))
          ;;                   (loop for sii in l maximize sii)))

          ;; (cl-mpm/output::save-parameter "boundary" (cl-mpm/particle::mp-boundary mp))
          ;; (cl-mpm/output::save-parameter "viscosity" (optional-slot-access 'cl-mpm/particle::viscosity mp))

          ;; ;; (cl-mpm/output::save-parameter "i-1" (cl-mpm/utils::trace-voigt (cl-mpm/particle::mp-stress mp)))

          ;; (when (= nd 3)
          ;;   (cl-mpm/output::save-parameter "eps_zz" (magicl:tref (cl-mpm/particle:mp-strain mp) 2 0))
          ;;   (cl-mpm/output::save-parameter "eps_yz" (magicl:tref (cl-mpm/particle:mp-strain mp) 3 0))
          ;;   (cl-mpm/output::save-parameter "eps_zx" (magicl:tref (cl-mpm/particle:mp-strain mp) 4 0)))


          ;; (save-parameter "size_x" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 0 0))
          ;; (save-parameter "size_y" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0))
          ;; (when (= nd 3)
          ;;   (save-parameter "size_z" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 2 0)))

          ;; (save-parameter "size_xx" (magicl:tref (cl-mpm/particle::mp-true-domain mp) 0 0))
          ;; (save-parameter "size_yy" (magicl:tref (cl-mpm/particle::mp-true-domain mp) 1 1))
          ;; (save-parameter "size_xy" (magicl:tref (cl-mpm/particle::mp-true-domain mp) 0 1))
          ;; (save-parameter "size_yx" (magicl:tref (cl-mpm/particle::mp-true-domain mp) 1 0))
          ;; (cl-mpm/output::save-parameter "fric-contact" (if (cl-mpm/particle::mp-penalty-contact-step mp) 1 0))
          ;; (cl-mpm/output::save-parameter "fric-contact-trial" (if (cl-mpm/particle::mp-penalty-contact mp) 1 0))
          ;; (cl-mpm/output::save-parameter "fric-contact-stick" (if (cl-mpm/particle::mp-penalty-friction-stick mp) 1 0))
          ;; (cl-mpm/output::save-parameter "fric-normal" (cl-mpm/particle::mp-penalty-normal-force mp))
          ;; (cl-mpm/output::save-parameter "fric-k" (cl-mpm/particle::mp-penalty-stiffness mp))
          ;; (cl-mpm/output::save-parameter "fric-x" (magicl:tref (cl-mpm/particle::mp-penalty-frictional-force mp) 0 0))
          ;; (cl-mpm/output::save-parameter "fric-y" (magicl:tref (cl-mpm/particle::mp-penalty-frictional-force mp) 1 0))
          ;; (cl-mpm/output::save-parameter "pressure" (cl-mpm/particle::mp-pressure mp))
          ;; ;; (save-parameter "pressure" (cl-mpm/particle::mp-pressure mp))
          ;; ;; (when (= (cl-mpm/mesh:mesh-nd mesh) 3)
          ;; ;;   (save-parameter "sig_zz" (magicl:tref (cl-mpm/particle:mp-stress mp) 2 0))
          ;; ;;   (save-parameter "size_z" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 2 0))
          ;; ;;   (save-parameter "sig_zx" (magicl:tref (cl-mpm/particle:mp-stress mp) 4 0))
          ;; ;;   (save-parameter "sig_xy" (magicl:tref (cl-mpm/particle:mp-stress mp) 5 0))
          ;; ;;   )
          ;; (save-parameter
          ;;  "plastic_strain"
          ;;  (if (slot-exists-p mp 'cl-mpm/particle::yield-func)
          ;;      (cl-mpm/particle::mp-strain-plastic-vm mp)
          ;;      0d0
          ;;      )
          ;;  )

          ;; (cl-mpm/output::save-parameter "plastic-iterations"
          ;;                                (if (slot-exists-p mp 'cl-mpm/particle::plastic-iterations)
          ;;                                    (cl-mpm/particle::mp-plastic-iterations mp)
          ;;                                    0d0))

          ;; (cl-mpm/output::save-parameter
          ;;  "rho"
          ;;  (if (slot-exists-p mp 'cl-mpm/particle::rho)
          ;;      (cl-mpm/particle::mp-rho mp)
          ;;      0d0))
          ;; (cl-mpm/output::save-parameter
          ;;  "c"
          ;;  (if (slot-exists-p mp 'cl-mpm/particle::c)
          ;;      (cl-mpm/particle::mp-c mp)
          ;;      0d0))
          ;; (cl-mpm/output::save-parameter
          ;;  "phi"
          ;;  (if (slot-exists-p mp 'cl-mpm/particle::phi)
          ;;      (* (cl-mpm/particle::mp-phi mp) (/ 180 pi))
          ;;      0d0))
          ;; (save-parameter
          ;;  "f"
          ;;  (if (slot-exists-p mp 'cl-mpm/particle::yield-func)
          ;;      (cl-mpm/particle::mp-yield-func mp)
          ;;      0d0)
          ;;  )
          )))
    ))

(defmacro with-parameter-list (file-stream mps &rest body)
  `(let ((id 1))
     (declare (special id))
     (format ,file-stream "POINT_DATA ~d~%" (length ,mps))
     ,(append
       (list 'progn)
       (loop for params in body
                   collect 
                   (destructuring-bind (name lamb) params
                     `(progn
                        (format-scalar ,file-stream ,name id ,mps
                                       ,lamb)
                        (incf id))
                     ))))

  )


;; (defmacro save-point-data (mp mps)
;;   (loop for )

;;   )

(defun format-scalar-cell (stream name id cells accessor)
  (format stream "SCALARS ~a FLOAT ~d~%" name 1)
  (format stream "LOOKUP_TABLE default~%")
  (loop for i from 0 below (array-total-size cells)
        do
           (let ((cell (row-major-aref cells i)))
             (when cell
               (format stream "~E ~%"
                       (coerce (funcall accessor cell) 'single-float)))))
  (format stream "~%"))
(defmacro save-parameter-cells (name accessor)
  `(progn
     (format-scalar-cell fs ,name id cells (lambda (cell) ,accessor))
     (incf id)))

(defgeneric save-vtk-cells (filename sim))

(defmethod save-vtk-cells (filename (sim cl-mpm::mpm-sim))
  (with-accessors ((mesh cl-mpm:sim-mesh)) sim
    (with-accessors ((cells cl-mpm/mesh::mesh-cells))
        mesh
        (with-open-file (fs filename :direction :output :if-exists :supersede)
          (format fs "# vtk DataFile Version 2.0~%")
          (format fs "Lisp generated vtk file, SJVS~%")
          (format fs "ASCII~%")
          (format fs "DATASET UNSTRUCTURED_GRID~%")

          (let (;; (node-count 0)
                (cells (remove-if-not (lambda (cell)
                                        (or
                                         (cl-mpm/mesh::cell-active cell)
                                         (and
                                          (cl-mpm/mesh::cell-mp-count cell)
                                          (= (cl-mpm/mesh::cell-mp-count cell) 1))))
                                      (make-array (array-total-size cells) :displaced-to cells)))
                )
            ;; (cl-mpm::iterate-over-cells-serial
            ;;  mesh
            ;;  (lambda (n)
            ;;    (declare (ignore n))
            ;;    (incf node-count)))
            (format fs "POINTS ~d double~%" (array-total-size cells))
            (loop for i from 0 below (array-total-size cells)
                  do
                     (let ((cell (row-major-aref cells i)))
                       (when cell
                         ;; (incf node-count)
                         (let ((pos (cl-mpm/mesh::cell-centroid cell)))
                           (format fs "~E ~E ~E ~%"
                                   (coerce (cl-mpm/utils::varef pos 0) 'single-float)
                                   (coerce (cl-mpm/utils::varef pos 1) 'single-float)
                                   (coerce (cl-mpm/utils::varef pos 2) 'single-float))))))
            (format fs "~%")
            (let ((id 1))
              (declare (special id))
              (format fs "POINT_DATA ~d~%" (array-total-size cells))
              ;; (save-parameter-cells "octree" (cl-mpm/dynamic-relaxation::cell-octree-refine cell))
              (save-parameter-cells "buoyancy" (if (cl-mpm/mesh::cell-boundary cell) 1 0))
              (save-parameter-cells "active" (if (cl-mpm/mesh::cell-active cell) 1 0))
              (save-parameter-cells "partial" (if (cl-mpm/mesh::cell-partial cell) 1 0))
              (save-parameter-cells "ghost" (if (cl-mpm/mesh::cell-ghost-element cell) 1 0))
              (save-parameter-cells "pressure" (cl-mpm/mesh::cell-pressure cell))
              (save-parameter-cells "cell-count" (cl-mpm/mesh::cell-mp-count cell))
              (save-parameter-cells "agg-int" (if (cl-mpm/mesh::cell-interior cell) 1d0 0d0))
              (save-parameter-cells "agg" (if (cl-mpm/mesh::cell-agg cell) 1d0 0d0))
              (save-parameter-cells "def-0" (if (cl-mpm/mesh::cell-def-list cell) (nth 0 (cl-mpm/mesh::cell-def-list cell)) 0))
              (save-parameter-cells "def-1" (if (cl-mpm/mesh::cell-def-list cell) (nth 1 (cl-mpm/mesh::cell-def-list cell)) 0))
              ;; (save-parameter-cells "def-2" (if (cl-mpm/mesh::cell-def-list cell) (nth 2 (cl-mpm/mesh::cell-def-list cell)) 0))
              (save-parameter-cells "j" (magicl:det (cl-mpm/mesh::cell-deformation-gradient cell)))
              ))))))
 
(defun save-vtk-bcs (filename sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)) sim
    (with-accessors ((bcs cl-mpm:sim-bcs))
        sim
        (with-open-file (fs filename :direction :output :if-exists :supersede)
          (format fs "# vtk DataFile Version 2.0~%")
          (format fs "Lisp generated vtk file, SJVS~%")
          (format fs "ASCII~%")
          (format fs "DATASET UNSTRUCTURED_GRID~%")

          (let ((node-count (length bcs))
                (h (cl-mpm/mesh:mesh-resolution mesh)))
            ;; (cl-mpm::iterate-over-bc-serial
            ;;  mesh
            ;;  (lambda (n)
            ;;    (declare (ignore n))
            ;;    (incf node-count)))
            (format fs "POINTS ~d double~%" node-count)
            (loop for bc across bcs
                  do
                     (let ((index (cl-mpm/bc:bc-index bc)))
                       (when index 
                         ;; (incf node-count)
                         (format fs "~E ~E ~E ~%"
                                 (coerce (* h (nth 0 index)) 'single-float)
                                 (coerce (* h (nth 1 index)) 'single-float)
                                 (coerce (* h (nth 2 index)) 'single-float)))))
            (format fs "~%")
            ;; (let ((id 1))
            ;;   (declare (special id))
            ;;   (format fs "POINT_DATA ~d~%" node-count)
            ;;   (save-parameter-cells "buoyancy" (if (cl-mpm/mesh::cell-boundary cell) 1 0))
            ;;   (save-parameter-cells "cell-count" (cl-mpm/mesh::cell-mp-count cell))
            ;;   )
            )
          ))))

(defun save-vtk-line (filename start end)

  )
(defun swap-endian (array)
  (reverse array)
  ;; (loop for i from 0 below (length array) by 2
  ;;       do (let ((msb (aref array i))
  ;;                (lsb (aref array (1+ i))))
  ;;            (setf (aref array i) lsb
  ;;                  (aref array (1+ i)) msb)))
  array)

(defmethod save-vtk-binary (filename (sim cl-mpm::mpm-sim))
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)) sim
    (with-open-file (fs filename :direction :output :if-exists :supersede)
      (with-open-file (fs-binary filename :direction :output :if-exists :append
                                          :element-type '(unsigned-byte 8)
                                          )
        
        (format fs "# vtk DataFile Version 2.0~%")
        (format fs "Lisp generated vtk file, SJVS~%")
        (format fs "BINARY~%")
        (format fs "DATASET UNSTRUCTURED_GRID~%")
        (format fs "POINTS ~d double~%" (length mps))
        (loop for mp across mps
              do
                 (let ((pos (cl-mpm/particle::mp-position mp)))
                   (dotimes (i 3)
                     (write-sequence (swap-endian (cl-intbytes:int64->octets (ieee-floats:encode-float64 (coerce (cl-mpm/utils:varef pos i) 'double-float)))) fs-binary)
                     )
                   ;; (format fs "~E ~E ~E ~%"
                   ;;         (coerce (cl-mpm/utils:varef pos 0) 'single-float)
                   ;;         (coerce (cl-mpm/utils:varef pos 1) 'single-float)
                   ;;         (coerce (cl-mpm/utils:varef pos 2) 'single-float))
                   )))
      (format fs "~%"))))
