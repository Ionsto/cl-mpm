(defpackage :cl-mpm/output
  (:use :cl
   :cl-mpm/utils)
  (:export
   #:save-vtk-mesh
   #:save-vtk
   #:save-vtk-nodes
   #:save-vtk-cells
   #:save-csv
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
                 (lambda (mesh mp node svp dsvp &rest rest) 
                   (setf av (+ av (* svp (,accessor node))))))
            av))))

(defmacro matrix-average (accessor shape)
    `(lambda (mesh mp)
        (let ((av (magicl:zeros ,shape)))
          (progn
            (cl-mpm::iterate-over-neighbours
             mesh mp
             (lambda (mesh mp node svp dsvp  &rest rest) 
               (setf av (cl-mpm/fastmaths::fast-.+ av (magicl:scale (,accessor node) svp)))))
            av))))

(defun sample-line-mass (sim start end steps)
    (sample-line sim start end steps 
        (scalar-average cl-mpm/mesh:node-mass)))

(defun sample-line-velocity (sim start end steps)
    (sample-line sim start end steps 
        (matrix-average cl-mpm/mesh:node-velocity '(2 1))))

(defun sample-point-mass (sim point)
    (sample-point sim point 
        (scalar-average cl-mpm/mesh:node-mass)))

(defun sample-point-velocity (sim point)
    (sample-point sim point
        (matrix-average cl-mpm/mesh:node-velocity '(2 1))))

(defun format-scalar (stream name id mps accessor)
  (format stream "SCALARS ~a FLOAT ~d~%" name 1)
  (format stream "LOOKUP_TABLE default~%")
    (loop for mp across mps
          do (format stream "~E ~%"
                     (coerce (funcall accessor mp) 'single-float)))
  (format stream "~%")
  )

(defmacro save-parameter (name accessor)
  (let ((mps (intern (symbol-name 'mps)))
        (fs (intern (symbol-name 'fs)))
        (mp (intern (symbol-name 'mp)))
        (id (intern (symbol-name 'id)))
        )
  `(progn
     (format-scalar ,fs ,name ,id ,mps (lambda (,mp) ,accessor))
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
;; (defun save-vtk (filename sim)
;;   (with-accessors ((mps cl-mpm:sim-mps)) sim
;;     (with-open-file (fs filename :direction :output :if-exists :supersede)
;;       (format fs "# vtk DataFile Version 2.0~%")
;;       (format fs "Lisp generated vtk file, WMC~%")
;;       (format fs "ASCII~%")
;;       (format fs "DATASET UNSTRUCTURED_GRID~%")
;;       (format fs "POINTS ~d double~%" (length mps))
;;       (loop for mp across mps
;;             do (format fs "~E ~E ~E ~%"
;;                        (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) 'single-float)
;;                        (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 1 0) 'single-float)
;;                        0e0))
;;       (format fs "~%")
;;       (let ((id 1))
;;         (declare (special id))
;;         (format fs "POINT_DATA ~d~%" (length mps))

;;         (save-parameter "mass" (cl-mpm/particle:mp-mass mp))
;;         (save-parameter "density" (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)))
;;         (save-parameter "index" (cl-mpm/particle::mp-index mp))
;;         (save-parameter "vel_x" (magicl:tref (cl-mpm/particle:mp-velocity mp) 0 0))
;;         (save-parameter "vel_y" (magicl:tref (cl-mpm/particle:mp-velocity mp) 1 0))
;;         (save-parameter "acc_x" (magicl:tref (cl-mpm/particle::mp-acceleration mp) 0 0))
;;         (save-parameter "acc_y" (magicl:tref (cl-mpm/particle::mp-acceleration mp) 1 0))
;;         (save-parameter "disp_x" (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0))
;;         (save-parameter "disp_y" (magicl:tref (cl-mpm/particle::mp-displacement mp) 1 0))
;;         (save-parameter "sig_xx" (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0))
;;         (save-parameter "sig_yy" (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0))
;;         (save-parameter "sig_xy" (magicl:tref (cl-mpm/particle:mp-stress mp) 2 0))

;;         (save-parameter "e_xx" (magicl:tref (cl-mpm/particle::mp-strain mp) 0 0))
;;         (save-parameter "e_yy" (magicl:tref (cl-mpm/particle::mp-strain mp) 1 0))
;;         (save-parameter "e_xy" (magicl:tref (cl-mpm/particle::mp-strain mp) 2 0))
;;         (save-parameter "temp" (magicl:tref (cl-mpm/particle::mp-velocity-rate mp) 2 0))

;;         (save-parameter "damage-inc-average"
;;                         (let ((v (/ (cl-mpm/particle::mp-time-averaged-damage-inc mp)
;;                                     (max 1d0
;;                                          (cl-mpm/particle::mp-time-averaged-counter mp)))))
;;                           (setf (cl-mpm/particle::mp-time-averaged-damage-inc mp) 0d0)
;;                           v))
;;         (save-parameter "ybar-average"
;;                         (let ((v (/ (cl-mpm/particle::mp-time-averaged-ybar mp)
;;                                     (max 1d0
;;                                          (cl-mpm/particle::mp-time-averaged-counter mp)))))
;;                                          (setf (cl-mpm/particle::mp-time-averaged-counter mp) 0d0
;;                                                (cl-mpm/particle::mp-time-averaged-ybar mp) 0d0)
;;                                          v))
;;         ;; (save-parameter "viscosity" (cl-mpm/particle::mp-time-averaged-visc mp))
;;         ;; (save-parameter "visc-plastic" (cl-mpm/particle::mp-visc-plastic mp))
;;         ;; (save-parameter "visc-glen" (cl-mpm/particle::mp-visc-glen mp))

;;         (save-parameter "strain_rate"
;;                         (cl-mpm/constitutive::effective-strain-rate (cl-mpm/particle::mp-eng-strain-rate mp))
;;                         ;; (multiple-value-bind (l v)
;;                         ;;     (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle::mp-velocity-rate mp)))
;;                         ;;   (reduce #'+ (mapcar #'* l l)))
;;                         )
;;         (save-parameter "pressure" (cl-mpm/particle::mp-pressure mp))
;;         ;; (save-parameter "pressure" (/ (+ (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0)
;;         ;;                                 (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0)) 2d0))
;;         (labels ((dot (a b) (magicl::sum (cl-mpm/fastmaths::fast-.* a b)))
;;                  (norm (a) (magicl:scale a (/ 1d0 (sqrt (dot a a)))))
;;                  (radial-stress (mp)
;;                    (with-accessors ((stress cl-mpm/particle:mp-stress)
;;                                     (pos cl-mpm/particle:mp-position))
;;                        mp
;;                      (let ((normal (norm (magicl:.- pos (magicl:from-list '(250d0 250d0) '(2 1))))))
;;                        (dot normal
;;                             (magicl:@ (cl-mpm/utils:voight-to-matrix stress)
;;                                       normal))))))
;;           (save-parameter "s_rr" (radial-stress mp)))

;;         (save-parameter "s_1"
;;                         (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))
;;                           (loop for sii in l maximize sii)))
;;         (save-parameter "s_vm"
;;                         (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))

;;                           (let* ((l (sort l #'>))
;;                                  (s_1 (max 0 (- (first l) (cl-mpm/particle::mp-pressure mp))))
;;                                  (s_2 (max 0 (- (second l) (cl-mpm/particle::mp-pressure mp)))))

;;                             (* (sqrt (/ 3 4)) (- s_1 s_2))
;;                             )
;;                           ))
;;         (save-parameter "s_vm_t"
;;                         (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))

;;                           (let* ((l (sort l #'>))
;;                                  (s_1 (first l))
;;                                  (s_2 (second l)))
;;                             (* ;(sqrt (/ 3 4))
;;                               1d0
;;                                (- s_1 s_2))
;;                             )
;;                           ))


;;         (save-parameter "EPS"
;;                         (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))
;;                           (- (loop for sii in l maximize sii) (cl-mpm/particle::mp-pressure mp))))
;;         (save-parameter "EPS-pd"
;;                         (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))
;;                           (- (loop for sii in l maximize sii) (* (cl-mpm/particle::mp-damage mp)
;;                                                                  (cl-mpm/particle::mp-pressure mp)))))
;;         (save-parameter "size_x" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 0 0))
;;         (save-parameter "size_y" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0))
;;         (save-parameter "damage"
;;                         (if (slot-exists-p mp 'cl-mpm/particle::damage)
;;                             (cl-mpm/particle:mp-damage mp)
;;                             0d0))
;;         (save-parameter "damage_inc"
;;                         (if (slot-exists-p mp 'cl-mpm/particle::damage-increment)
;;                             (cl-mpm/particle::mp-damage-increment mp)
;;                             0d0))
;;         (save-parameter "damage_ybar"
;;                         (if (slot-exists-p mp 'cl-mpm/particle::damage-ybar)
;;                             (cl-mpm/particle::mp-damage-ybar mp)
;;                             0d0))
;;         (save-parameter "local_length"
;;                         (if (slot-exists-p mp 'cl-mpm/particle::true-local-length)
;;                             (cl-mpm/particle::mp-true-local-length mp)
;;                             0d0))
;;         )
;;       )))

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

(defun format-scalar-node (stream name id nodes accessor)
  (format stream "SCALARS ~a FLOAT ~d~%" name 1)
  (format stream "LOOKUP_TABLE default~%")
  (loop for i from 0 below (array-total-size nodes)
        do
           (let ((node (row-major-aref nodes i)))
             (when node
               (format stream "~E ~%" (coerce (funcall accessor node) 'single-float)))))
  (format stream "~%"))

(defmacro save-parameter-nodes (name accessor)
  `(progn
     (format-scalar-node fs ,name id nodes (lambda (node) ,accessor))
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

(defun save-vtk-mesh-nodes (filename mesh)
  (with-accessors ((nodes cl-mpm/mesh::mesh-nodes))
      mesh
    (with-open-file (fs filename :direction :output :if-exists :supersede)
      (format fs "# vtk DataFile Version 2.0~%")
      (format fs "Lisp generated vtk file, SJVS~%")
      (format fs "ASCII~%")
      (format fs "DATASET UNSTRUCTURED_GRID~%")

      (let* ((node-count (array-total-size nodes))
             (nodes (remove-if-not #'cl-mpm/mesh:node-active
                                   (make-array node-count :displaced-to nodes)))
             (node-count (array-total-size nodes)))
        (format fs "POINTS ~d double~%" (array-total-size nodes))
        (loop for i from 0 below (array-total-size nodes)
              do
                 (let ((node (row-major-aref nodes i)))
                   (when node
                     ;; (incf node-count)
                     (let ((pos (cl-mpm/fastmaths:fast-.+
                                 (cl-mpm/mesh::node-position node)
                                 (cl-mpm/mesh::node-displacment node))))
                       (format fs "~E ~E ~E ~%"
                               (coerce (cl-mpm/utils::varef pos 0) 'single-float)
                               (coerce (cl-mpm/utils::varef pos 1) 'single-float)
                               (coerce (cl-mpm/utils::varef pos 2) 'single-float))))))
        (format fs "~%")
        (let ((id 1))
          (declare (special id))
          (format fs "POINT_DATA ~d~%" (array-total-size nodes))


          (save-parameter-nodes "active" (if (cl-mpm/mesh::node-active node) 1 0))

          (save-parameter-nodes "agg" (if (cl-mpm/mesh::node-agg node) 1 0))
          (save-parameter-nodes "agg-int" (if (cl-mpm/mesh::node-interior node) 1 0))

          (save-parameter-nodes "mass" (cl-mpm/mesh:node-mass node))
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

          ;; (save-parameter-nodes "disp_x" (magicl:tref (cl-mpm/mesh::node-displacment node) 0 0))
          ;; (save-parameter-nodes "disp_y" (magicl:tref (cl-mpm/mesh::node-displacment node) 1 0))
          ;; (save-parameter-nodes "disp_z" (magicl:tref (cl-mpm/mesh::node-displacment node) 2 0))

          (save-vector-nodes "force" (cl-mpm/mesh:node-force node))
          (save-vector-nodes "force-buoyancy" (cl-mpm/mesh::node-buoyancy-force node))
          (save-vector-nodes "force-ghost" (cl-mpm/mesh::node-ghost-force node))
          (save-vector-nodes "acc" (cl-mpm/mesh:node-acceleration node))
          ;; (save-vector-nodes "force-rct" #'cl-mpm/mesh::node-force-)

          ;; (save-parameter-nodes "force_x" (magicl:tref (cl-mpm/mesh:node-force node) 0 0))
          ;; (save-parameter-nodes "force_y" (magicl:tref (cl-mpm/mesh:node-force node) 1 0))
          ;; (save-parameter-nodes "force_z" (magicl:tref (cl-mpm/mesh:node-force node) 2 0))
          ;; (save-parameter-nodes "force_buoy_x" (magicl:tref (cl-mpm/mesh::node-buoyancy-force node) 0 0))
          ;; (save-parameter-nodes "force_buoy_y" (magicl:tref (cl-mpm/mesh::node-buoyancy-force node) 1 0))
          ;; (save-parameter-nodes "force_buoy_z" (magicl:tref (cl-mpm/mesh::node-buoyancy-force node) 2 0))

          ;; (save-parameter-nodes "force_damp_x" (magicl:tref (cl-mpm/mesh::node-damping-force node) 0 0))
          ;; (save-parameter-nodes "force_damp_y" (magicl:tref (cl-mpm/mesh::node-damping-force node) 1 0))
          ;; (save-parameter-nodes "force_damp_z" (magicl:tref (cl-mpm/mesh::node-damping-force node) 2 0))

          ;; (save-parameter-nodes "force_ghost_x" (varef (cl-mpm/mesh::node-ghost-force node) 0))
          ;; (save-parameter-nodes "force_ghost_y" (varef (cl-mpm/mesh::node-ghost-force node) 1))
          (save-parameter-nodes "j-inc" (if (= (cl-mpm/mesh::node-jacobian-inc node) 0d0)
                                            1d0
                                            (cl-mpm/mesh::node-jacobian-inc node)
                                            ))
          (save-parameter-nodes "buoyancy_node" (if
                                                 (cl-mpm/mesh::node-boundary-node node) 1 0))
          (save-parameter-nodes "buoyancy-scalar" (cl-mpm/mesh::node-boundary-scalar node))
          (save-parameter-nodes "pressure" (cl-mpm/mesh::node-pressure node))
          (save-parameter-nodes "local-list-size" (length (cl-mpm/mesh::node-local-list node)))
          (save-parameter-nodes "energy"
                                (*
                                 ;; (/ (cl-mpm/mesh::node-volume node) (cl-mpm/mesh::node-volume-true node))
                                 (cl-mpm/mesh::node-mass node)
                                 (cl-mpm/mesh::node-mass node)
                                 (cl-mpm/fastmaths::mag-squared (cl-mpm/mesh::node-velocity node))))

          (save-parameter-nodes "oobf" (cl-mpm/mesh::node-oobf node))
          (save-parameter-nodes "volume" (cl-mpm/mesh::node-volume node))
          (save-parameter-nodes "volume-true" (cl-mpm/mesh::node-volume-true node))
          ;; (save-parameter-nodes "lift" (- (magicl:tref (cl-mpm/mesh:node-force node) 1 0)
          ;;                                 (cl-mpm/mesh::node-pressure node)))
          (save-parameter-nodes "damage"
                                (if (slot-exists-p node 'cl-mpm/mesh::damage)
                                    (if (= 0d0 (cl-mpm/mesh::node-svp-sum node))
                                        0d0
                                        (/ (cl-mpm/mesh::node-damage node)
                                           (cl-mpm/mesh::node-svp-sum node))) 0d0))
          )))))

(defun save-vtk-nodes (filename sim)
  (save-vtk-mesh-nodes filename (cl-mpm:sim-mesh sim)))


(defgeneric save-vtk (filename sim)
  (:documentation "Save vtk depending on type")
   (:method (f s)))

(defmacro optional-slot-access (slot mp)
  `(if (slot-exists-p ,mp ,slot)
       (slot-value ,mp ,slot)
       0d0))

(defmethod save-vtk (filename (sim cl-mpm::mpm-sim))
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)) sim
    (with-open-file (fs filename :direction :output :if-exists :supersede)
      (format fs "# vtk DataFile Version 2.0~%")
      (format fs "Lisp generated vtk file, SJVS~%")
      (format fs "ASCII~%")
      (format fs "DATASET UNSTRUCTURED_GRID~%")
      (format fs "POINTS ~d double~%" (length mps))
      (loop for mp across mps
            do
               (let ((pos (cl-mpm/particle::mp-position-trial mp)))
                 (format fs "~E ~E ~E ~%"
                         (coerce (cl-mpm/utils:varef pos 0) 'single-float)
                         (coerce (cl-mpm/utils:varef pos 1) 'single-float)
                         (coerce (cl-mpm/utils:varef pos 2) 'single-float))))
      (format fs "~%")

      (let ((id 1)
            (nd (cl-mpm/mesh:mesh-nd mesh)))
        (declare (special id))
        (format fs "POINT_DATA ~d~%" (length mps))

        (save-parameter "mass" (cl-mpm/particle:mp-mass mp))
        (save-parameter "density" (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)))
        (save-parameter "volume" (cl-mpm/particle:mp-volume mp))
        (save-parameter "index" (cl-mpm/particle::mp-index mp))
        (save-parameter "mpi-domain" (cl-mpm/particle::mp-mpi-index mp))
        (save-parameter "split-depth" (cl-mpm/particle::mp-split-depth mp))

        ;; (save-parameter-vector "vel" cl-mpm/particle::mp-velocity)
        (save-parameter "vel_x" (magicl:tref (cl-mpm/particle:mp-velocity mp) 0 0))
        (save-parameter "vel_y" (magicl:tref (cl-mpm/particle:mp-velocity mp) 1 0))
        ;; (save-parameter "vel_z" (magicl:tref (cl-mpm/particle:mp-velocity mp) 2 0))

        ;; (save-parameter "acc_x" (magicl:tref (cl-mpm/particle::mp-acceleration mp) 0 0))
        ;; (save-parameter "acc_y" (magicl:tref (cl-mpm/particle::mp-acceleration mp) 1 0))

        (save-parameter "disp_x" (magicl:tref (cl-mpm/fastmaths:fast-.+
                                               (cl-mpm/particle::mp-displacement-increment mp)
                                               (cl-mpm/particle::mp-displacement mp)) 0 0))
        (save-parameter "disp_y" (magicl:tref (cl-mpm/fastmaths:fast-.+
                                               (cl-mpm/particle::mp-displacement-increment mp)
                                               (cl-mpm/particle::mp-displacement mp)) 1 0))
        (save-parameter "disp_z" (magicl:tref (cl-mpm/fastmaths:fast-.+
                                               (cl-mpm/particle::mp-displacement-increment mp)
                                               (cl-mpm/particle::mp-displacement mp)) 2 0))

        (cl-mpm/output::save-parameter "sig_xx" (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0))
        (cl-mpm/output::save-parameter "sig_yy" (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0))
        (cl-mpm/output::save-parameter "sig_xy" (magicl:tref (cl-mpm/particle:mp-stress mp) 5 0))
        (cl-mpm/output::save-parameter "i1" (cl-mpm/utils::trace-voigt (cl-mpm/particle:mp-stress mp)))
        (when (= nd 3)
          (cl-mpm/output::save-parameter "sig_zz" (magicl:tref (cl-mpm/particle:mp-stress mp) 2 0))
          (cl-mpm/output::save-parameter "sig_yz" (magicl:tref (cl-mpm/particle:mp-stress mp) 3 0))
          (cl-mpm/output::save-parameter "sig_zx" (magicl:tref (cl-mpm/particle:mp-stress mp) 4 0)))

        (cl-mpm/output::save-parameter "eps_xx" (magicl:tref (cl-mpm/particle:mp-strain mp) 0 0))
        (cl-mpm/output::save-parameter "eps_yy" (magicl:tref (cl-mpm/particle:mp-strain mp) 1 0))
        (cl-mpm/output::save-parameter "eps_xy" (magicl:tref (cl-mpm/particle:mp-strain mp) 5 0))
        (save-parameter "eps_1"
                        (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))
                          (loop for sii in l maximize sii)))

        (cl-mpm/output::save-parameter "boundary" (cl-mpm/particle::mp-boundary mp))
        (cl-mpm/output::save-parameter "viscosity" (optional-slot-access 'cl-mpm/particle::viscosity mp))

        ;; (cl-mpm/output::save-parameter "i-1" (cl-mpm/utils::trace-voigt (cl-mpm/particle::mp-stress mp)))

        (when (= nd 3)
          (cl-mpm/output::save-parameter "eps_zz" (magicl:tref (cl-mpm/particle:mp-strain mp) 2 0))
          (cl-mpm/output::save-parameter "eps_yz" (magicl:tref (cl-mpm/particle:mp-strain mp) 3 0))
          (cl-mpm/output::save-parameter "eps_zx" (magicl:tref (cl-mpm/particle:mp-strain mp) 4 0)))


        (save-parameter "size_x" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 0 0))
        (save-parameter "size_y" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0))
        (when (= nd 3)
          (save-parameter "size_z" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 2 0)))
        (cl-mpm/output::save-parameter "fric-contact" (if (cl-mpm/particle::mp-penalty-contact-step mp) 1 0))
        (cl-mpm/output::save-parameter "fric-contact-trial" (if (cl-mpm/particle::mp-penalty-contact mp) 1 0))
        (cl-mpm/output::save-parameter "fric-contact-stick" (if (cl-mpm/particle::mp-penalty-friction-stick mp) 1 0))
        (cl-mpm/output::save-parameter "fric-normal" (cl-mpm/particle::mp-penalty-normal-force mp))
        (cl-mpm/output::save-parameter "fric-k" (cl-mpm/particle::mp-penalty-stiffness mp))
        (cl-mpm/output::save-parameter "fric-x" (magicl:tref (cl-mpm/particle::mp-penalty-frictional-force mp) 0 0))
        (cl-mpm/output::save-parameter "fric-y" (magicl:tref (cl-mpm/particle::mp-penalty-frictional-force mp) 1 0))
        (cl-mpm/output::save-parameter "pressure" (cl-mpm/particle::mp-pressure mp))
        ;; (save-parameter "pressure" (cl-mpm/particle::mp-pressure mp))
        ;; (when (= (cl-mpm/mesh:mesh-nd mesh) 3)
        ;;   (save-parameter "sig_zz" (magicl:tref (cl-mpm/particle:mp-stress mp) 2 0))
        ;;   (save-parameter "size_z" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 2 0))
        ;;   (save-parameter "sig_zx" (magicl:tref (cl-mpm/particle:mp-stress mp) 4 0))
        ;;   (save-parameter "sig_xy" (magicl:tref (cl-mpm/particle:mp-stress mp) 5 0))
        ;;   )
        (save-parameter
         "plastic_strain"
         (if (slot-exists-p mp 'cl-mpm/particle::yield-func)
             (cl-mpm/particle::mp-strain-plastic-vm mp)
             0d0
             )
         )

        (cl-mpm/output::save-parameter "plastic-iterations"
                                       (if (slot-exists-p mp 'cl-mpm/particle::plastic-iterations)
                                           (cl-mpm/particle::mp-plastic-iterations mp)
                                           0d0))

        (cl-mpm/output::save-parameter
         "rho"
         (if (slot-exists-p mp 'cl-mpm/particle::rho)
             (cl-mpm/particle::mp-rho mp)
             0d0))
        (cl-mpm/output::save-parameter
         "c"
         (if (slot-exists-p mp 'cl-mpm/particle::c)
             (cl-mpm/particle::mp-c mp)
             0d0))
        (cl-mpm/output::save-parameter
         "phi"
         (if (slot-exists-p mp 'cl-mpm/particle::phi)
             (* (cl-mpm/particle::mp-phi mp) (/ 180 pi))
             0d0))
        (save-parameter
         "f"
         (if (slot-exists-p mp 'cl-mpm/particle::yield-func)
             (cl-mpm/particle::mp-yield-func mp)
             0d0)
         )
        )
      )))

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

(defun save-vtk-cells (filename sim)
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
                         (let ((pos (cl-mpm/mesh::cell-trial-centroid cell)))
                           (format fs "~E ~E ~E ~%"
                                   (coerce (cl-mpm/utils::varef pos 0) 'single-float)
                                   (coerce (cl-mpm/utils::varef pos 1) 'single-float)
                                   (coerce (cl-mpm/utils::varef pos 2) 'single-float))))))
            (format fs "~%")
            (let ((id 1))
              (declare (special id))
              (format fs "POINT_DATA ~d~%" (array-total-size cells))
              (save-parameter-cells "buoyancy" (if (cl-mpm/mesh::cell-boundary cell) 1 0))
              (save-parameter-cells "active" (if (cl-mpm/mesh::cell-active cell) 1 0))
              (save-parameter-cells "pressure" (cl-mpm/mesh::cell-pressure cell))
              (save-parameter-cells "cell-count" (cl-mpm/mesh::cell-mp-count cell))
              (save-parameter-cells "agg-int" (if (cl-mpm/mesh::cell-interior cell) 1d0 0d0))
              (save-parameter-cells "agg" (if (cl-mpm/mesh::cell-agg cell) 1d0 0d0))
              (save-parameter-cells "j" (magicl:det (cl-mpm/mesh::cell-deformation-gradient cell)))
              (save-parameter-cells "ghost" (if (cl-mpm/mesh::cell-ghost-element cell) 1 0))))))))
 
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
  (loop for i from 0 below (length array) by 2
        do (let ((msb (aref array i))
                 (lsb (aref array (1+ i))))
             (setf (aref array i) lsb
                   (aref array (1+ i)) msb)))
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
                     (write-sequence (swap-endian (cl-intbytes:int32->octets (ieee-floats:encode-float32 (coerce (cl-mpm/utils:varef pos i) 'single-float)))) fs-binary)
                     )
                   ;; (format fs "~E ~E ~E ~%"
                   ;;         (coerce (cl-mpm/utils:varef pos 0) 'single-float)
                   ;;         (coerce (cl-mpm/utils:varef pos 1) 'single-float)
                   ;;         (coerce (cl-mpm/utils:varef pos 2) 'single-float))
                   )))
      (format fs "~%"))))
