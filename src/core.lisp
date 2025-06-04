;    #:make-shape-function
(in-package :cl-mpm)

;(declaim (optimize (debug 3) (safety 3) (speed 0)))
;; (declaim (optimize (debug 3) (safety 3) (speed 0)))
;; (declaim (optimize (debug 0) (safety 0) (speed 3)))
(declaim #.cl-mpm/settings:*optimise-setting*)

(eval-when (:compile-toplevel :load-toplevel :execute)
  #+cl-mpm-special (print "Compiled with special hooks")
  #-cl-mpm-special (print "Compiled without special hooks")
  (format t "Compiled with optimize ~A~%" (sb-ext:restrict-compiler-policy)))

(defun FLIP-status ()
  #+cl-mpm-pic (print "Compiled with PIC")
  #-cl-mpm-pic (print "Compiled with FLIP")
  )

(declaim (notinline check-mps))
(defun check-mps (sim)
  "Function to check that stresses and positions are sane, deleting mps moving very fast"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)
                   (removal-factor sim-mp-removal-size))
      sim
  (with-accessors ((h cl-mpm/mesh::mesh-resolution))
      mesh
    (when removal-factor
      (let ((h (* h removal-factor)))
        (remove-mps-func
         sim
         (lambda (mp)
           (with-accessors ((damage cl-mpm/particle:mp-damage)
                            (def cl-mpm/particle::mp-deformation-gradient))
               mp
             (or
              (gimp-removal-criteria mp h))))))))))

(defun check-single-mps (sim)
  "Function to check and remove single material points"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
  (with-accessors ((h cl-mpm/mesh::mesh-resolution))
      mesh
    (iterate-over-mps
     mps
     (lambda (mp)
       (setf (cl-mpm/particle::mp-single-particle mp)
             (single-particle-criteria mesh mp))))
    (remove-mps-func
     sim
     #'cl-mpm/particle::mp-single-particle))))


(defgeneric update-sim (sim)
  (:documentation "Update an mpm simulation by one timestep"))



(defun update-node (node dt)
  "Calculate velocity from momentum on a single node"
  (when (cl-mpm/mesh:node-active node)
    (with-accessors ((vel   node-velocity)
                     (disp   cl-mpm/mesh::node-displacment))
        node
      (cl-mpm/fastmaths:fast-fmacc disp vel dt))))

(defun get-cell-df (mesh point)
  (let ((dF (cl-mpm/utils:matrix-eye 1d0)))
    (cl-mpm::iterate-over-neighbours-point-linear
     mesh
     point
     (lambda (mesh node weight grads)
       (when (cl-mpm/mesh::node-active node)
         (cl-mpm/shape-function::@-combi-assemble-dstretch-3d grads (cl-mpm/mesh::node-displacment node) df))))
    dF))

(defun filter-cell (mesh cell dt)
  (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                   (active cl-mpm/mesh::cell-active)
                   (neighbours cl-mpm/mesh::cell-neighbours)
                   (index cl-mpm/mesh::cell-index)
                   (nodes cl-mpm/mesh::cell-nodes)
                   (df cl-mpm/mesh::cell-deformation-gradient)
                   (disp cl-mpm/mesh::cell-displacement)
                   (centroid cl-mpm/mesh::cell-centroid)
                   (trial-pos cl-mpm/mesh::cell-trial-centroid)
                   )
      cell
    (setf mp-count t)
    (setf active nil)
    (loop for n in nodes
          do
             (let ((oc (cl-mpm/mesh:node-active n)))
               (setf mp-count (and mp-count oc))
               (setf active (or active oc))))
    (setf mp-count (if mp-count 1 0))))

(defgeneric filter-cells (sim))
(defmethod filter-cells (sim)
  (with-accessors ((mesh sim-mesh)
                   (dt sim-dt))
      sim
    (iterate-over-cells
     mesh
     (lambda (cell)
       (filter-cell mesh cell dt)))))

(defun update-cell (mesh cell dt)
  "Update cell data, useful for ghost penalty"
  (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                   (active cl-mpm/mesh::cell-active)
                   (neighbours cl-mpm/mesh::cell-neighbours)
                   (index cl-mpm/mesh::cell-index)
                   (nodes cl-mpm/mesh::cell-nodes)
                   (df cl-mpm/mesh::cell-deformation-gradient)
                   (disp cl-mpm/mesh::cell-displacement)
                   (centroid cl-mpm/mesh::cell-centroid)
                   (trial-pos cl-mpm/mesh::cell-trial-centroid)
                   )
      cell
    (when active
      (cl-mpm/fastmaths:fast-zero disp)
      (cl-mpm/fastmaths:fast-zero df)
      (setf (varef df 0) 1d0
            (varef df 4) 1d0
            (varef df 8) 1d0)
      (let ((w 0d0))
        (declare (double-float w))
        (cl-mpm::iterate-over-neighbours-point-linear
         mesh
         centroid
         (lambda (mesh node weight grads)
           (declare (double-float weight))
           (when (cl-mpm::node-active node)
             (incf w weight)
             (let ((ndisp (cl-mpm/mesh::node-displacment node)))
               (cl-mpm/shape-function::@-combi-assemble-dstretch-3d grads ndisp df)
               (cl-mpm/fastmaths:fast-fmacc
                disp
                ndisp weight)))))
        (when (> w 0d0)
          (cl-mpm::fast-scale! disp (/ 1d0 w)))
        )
      (cl-mpm/fastmaths:fast-.+ centroid disp trial-pos)))
  )

(defgeneric update-cells (sim))
(defmethod update-cells (sim)
  (with-accessors ((mesh sim-mesh)
                   (dt sim-dt))
      sim
    (iterate-over-cells
     mesh
     (lambda (cell)
       (filter-cell mesh cell dt)
       (update-cell mesh cell dt)))))

(defgeneric update-nodes (sim))
(defmethod update-nodes (sim)
  (with-accessors ((mesh sim-mesh)
                   (dt sim-dt))
      sim
    (iterate-over-nodes
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (update-node node dt))))))

(defun reset-nodes-force (sim)
  (cl-mpm::iterate-over-nodes
   (cl-mpm:sim-mesh sim)
   (lambda (n)
     (when (cl-mpm/mesh::node-active n)
       (cl-mpm/mesh::reset-node-force n)))))


(defgeneric pre-particle-update-hook (particle dt))
(defmethod pre-particle-update-hook (particle dt))

(defgeneric special-update-node (mesh dt node damping)
  (:documentation "Update node method")
  (:method (mesh dt node damping)))

(defmethod special-update-node (mesh dt (node cl-mpm/mesh::node-thermal) damping)
  (with-accessors ((mass  node-mass)
                   (temp node-temperature)
                   (dtemp node-dtemp))
      node
    (progn
      (setf temp (+ (/ temp mass) (* dtemp dt)))
      )))

(defun integrate-vel-euler (vel acc force mass mass-scale dt damping)
  (declare (double-float mass mass-scale dt damping))
  (unless (= dt 0d0)
    (if (= damping 0d0)
        (cl-mpm/fastmaths:fast-fmacc vel acc dt)
        (let ((exp-fac (exp (- (* damping dt)))))
          ;; (cl-mpm/fastmaths:fast-fmacc vel acc dt)
          (cl-mpm/fastmaths:fast-scale! vel exp-fac)
          (cl-mpm/fastmaths:fast-fmacc vel acc (* (/ 1d0 damping) (- 1d0 exp-fac)))
          ))))

(declaim (notinline calculate-forces)
         (ftype (function (cl-mpm/mesh::node double-float double-float double-float) (vaules)) calculate-forces))
(defun calculate-forces (node damping dt mass-scale)
  "Update forces and nodal velocities with viscous damping"
  (when (cl-mpm/mesh:node-active node)
    (with-accessors ((mass node-mass)
                     (vel node-velocity)
                     (force node-force)
                     (force-ext cl-mpm/mesh::node-external-force)
                     (force-int cl-mpm/mesh::node-internal-force)
                     (force-damp cl-mpm/mesh::node-damping-force)
                     (force-ghost cl-mpm/mesh::node-ghost-force)
                     (residual cl-mpm/mesh::node-residual)
                     (residual-prev cl-mpm/mesh::node-residual-prev)
                     (acc node-acceleration))
        node
      (declare (double-float mass dt damping mass-scale))
      (progn
        (cl-mpm/fastmaths:fast-zero acc)
        ;;Set acc to f/m
        (cl-mpm/fastmaths::fast-.+-vector force-int force force)
        (cl-mpm/fastmaths::fast-.+-vector force-ext force force)
        ;;Include velocity prop damping
        (cl-mpm/fastmaths::fast-.+-vector force-damp force force)
        (cl-mpm/fastmaths::fast-.+-vector force-ghost force force)
        (cl-mpm/fastmaths:fast-fmacc acc force (/ 1d0 (* mass mass-scale)))
        (integrate-vel-euler vel acc force mass mass-scale dt damping)
        (cl-mpm/utils::vector-copy-into residual residual-prev)
        (cl-mpm/utils::vector-copy-into force-int residual))))
  (values))


(defun integrate-vel-midpoint (vel acc force mass mass-scale dt damping)
  (declare (double-float mass mass-scale dt damping))
  (unless (= dt 0d0)
    (let ((damp-dt (/(* dt damping) mass-scale)))
      (cl-mpm/fastmaths:fast-scale! vel (/ (- 2d0 damp-dt)
                                           (+ 2d0 damp-dt)))
      (cl-mpm/fastmaths:fast-fmacc vel acc (/ (* dt 2d0) (+ 2d0 damp-dt))))))

(declaim (notinline calculate-forces-midpoint)
         (ftype (function (cl-mpm/mesh::node double-float double-float double-float) (vaules)) calculate-forces-midpoint))
(defun calculate-forces-midpoint (node damping dt mass-scale)
  "Update forces and nodal velocities with viscous damping"
  (when (cl-mpm/mesh:node-active node)
    (with-accessors ((mass node-mass)
                     (vel node-velocity)
                     (force node-force)
                     (force-ext cl-mpm/mesh::node-external-force)
                     (force-int cl-mpm/mesh::node-internal-force)
                     (force-damp cl-mpm/mesh::node-damping-force)
                     (force-ghost cl-mpm/mesh::node-ghost-force)
                     (residual cl-mpm/mesh::node-residual)
                     (residual-prev cl-mpm/mesh::node-residual-prev)
                     (acc node-acceleration))
        node
      (declare (double-float mass dt damping mass-scale))
      (progn
        (cl-mpm/fastmaths:fast-zero acc)
        ;;Set acc to f/m
        (cl-mpm/fastmaths::fast-.+-vector force-int force force)
        (cl-mpm/fastmaths::fast-.+-vector force-ext force force)
        ;; (cl-mpm/fastmaths:fast-fmacc force-damp vel (* damping -1d0 mass))
        (cl-mpm/fastmaths::fast-.+-vector force-damp force force)
        (cl-mpm/fastmaths::fast-.+-vector force-ghost force force)
        (cl-mpm/fastmaths:fast-fmacc acc force (/ 1d0 (* mass mass-scale)))

        (integrate-vel-midpoint vel acc force mass mass-scale dt damping)

        (cl-mpm/utils::vector-copy-into residual residual-prev)
        (cl-mpm/utils::vector-copy-into force-int residual))))
  (values))

(defun calculate-forces-psudo-viscous (node damping dt mass-scale)
  "Update forces and nodal velocities with viscous damping - except without scaling by mass
This allows for a non-physical but viscous damping scheme that is robust to GIMP domains "
  (when (cl-mpm/mesh:node-active node)
    (with-accessors ( (mass  node-mass)
                     (vel   node-velocity)
                     (force node-force)
                     (force-ext cl-mpm/mesh::node-external-force)
                     (force-int cl-mpm/mesh::node-internal-force)
                     (acc   node-acceleration))
        node
      (declare (double-float mass dt damping mass-scale))
      (progn
        (magicl:scale! acc 0d0)
        ;;Set acc to f/m
        (cl-mpm/fastmaths::fast-.+-vector force-int force-ext force)
        (cl-mpm/fastmaths:fast-fmacc acc force (/ 1d0 (* mass mass-scale)))
        (cl-mpm/fastmaths:fast-fmacc acc vel (* damping -1d0))
        (cl-mpm/fastmaths:fast-fmacc vel acc dt)
        )))
  (values))

(defun calculate-forces-cundall (node damping dt mass-scale)
  "Apply cundall damping to the system"
  (when (cl-mpm/mesh:node-active node)
    (with-accessors ((mass  node-mass)
                     (vel   node-velocity)
                     (force node-force)
                     (force-ext cl-mpm/mesh::node-external-force)
                     (force-int cl-mpm/mesh::node-internal-force)
                     (acc   node-acceleration)
                     )
        node
      (declare (double-float mass dt damping mass-scale))
      (progn
        (magicl:scale! acc 0d0)
        ;;Set acc to f/m
        (cl-mpm/fastmaths::fast-.+-vector force-int force-ext force)
        (let* ((vel-sign (cl-mpm/utils:vector-zeros))
               (f-s (magicl::matrix/double-float-storage force))
               (v-s (magicl::matrix/double-float-storage vel))
               (fnorm (cl-mpm/fastmaths::mag force))
               (vs-s (magicl::matrix/double-float-storage vel-sign)))
          (loop for i from 0 to 2
                do
                   ;;Possible form of damping
                   (when t;(> (* (aref f-s i) (aref v-s i)) 0)
                     (incf (aref f-s i)
                           (*
                            (signum (aref v-s i))
                            (abs (aref f-s i))
                            damping
                            -1d0))))
          (cl-mpm/fastmaths:fast-fmacc acc force (/ 1d0 (* mass mass-scale))))
        (cl-mpm/fastmaths:fast-fmacc vel acc dt))))
  (values))
(defun calculate-forces-cundall-conservative (node damping dt mass-scale)
  "Apply cundall damping to the system only when doing negative work"
  (when (cl-mpm/mesh:node-active node)
    (with-accessors ((mass  node-mass)
                     (vel   node-velocity)
                     (force node-force)
                     (acc   node-acceleration)
                     )
        node
      (declare (double-float mass dt damping))
      (progn
        (magicl:scale! acc 0d0)
        ;;Set acc to f/m

        (let* ((vel-sign (cl-mpm/utils:vector-zeros))
               (f-s (magicl::matrix/double-float-storage force))
               (v-s (magicl::matrix/double-float-storage vel))
               (fnorm (cl-mpm/fastmaths::mag force))
               (vs-s (magicl::matrix/double-float-storage vel-sign)))
          (loop for i from 0 to 2
                do
                   ;;Possible form of damping
                   (when (> (* (aref f-s i) (aref v-s i)) 0)
                     (incf (aref f-s i)
                           (*
                            (signum (aref v-s i))
                            (abs (aref f-s i))
                            damping
                            -1d0)))
                )
          (cl-mpm/fastmaths:fast-fmacc acc force (/ 1d0 (* mass mass-scale))))
        (cl-mpm/fastmaths:fast-fmacc vel acc dt))))
  (values))

(declaim (notinline calculate-forces-water-viscous)
         (ftype (function (cl-mpm/mesh::node double-float double-float double-float double-float) (vaules)) calculate-forces-water-viscous))
(defun calculate-forces-water-viscous (node damping damping-water dt mass-scale)
  "Update forces and nodal velocities with viscous damping"
  (when (cl-mpm/mesh:node-active node)
    (with-accessors ((mass  node-mass)
                     (vel   node-velocity)
                     (boundary   cl-mpm/mesh::node-boundary-node)
                     (force node-force)
                     (force-ext cl-mpm/mesh::node-external-force)
                     (force-int cl-mpm/mesh::node-internal-force)
                     (acc   node-acceleration))
        node
      (declare (double-float mass dt damping mass-scale))
      (progn
        (cl-mpm/fastmaths:fast-zero acc)
        ;;Set acc to f/m
        ;; (cl-mpm/fastmaths:fast-fmacc force-ext vel (* -2d0 damping mass))
        (cl-mpm/fastmaths::fast-.+-vector force-int force-ext force)
        ;; (cl-mpm/fastmaths:fast-fmacc acc vel (/ (* damping -1d0) 1d0))
        (cl-mpm/fastmaths:fast-fmacc acc force (/ 1d0 (* mass mass-scale)))
        (cl-mpm/fastmaths:fast-fmacc acc vel (/ (* damping -1d0) mass-scale))
        (when boundary
          (cl-mpm/fastmaths:fast-fmacc acc vel (/ (* damping-water -1d0) mass-scale)))
        ;; (cl-mpm/fastmaths:fast-fmacc acc vel (/ (* damping -1d0) (sqrt mass-scale)))
        ;; (cl-mpm/fastmaths:fast-fmacc acc vel (* damping -1d0))
        (cl-mpm/fastmaths:fast-fmacc vel acc dt)
        )))
  (values))

(defun calculate-kinematics (node)
  "Calculate velocity from momentum on a single node"
  (when (cl-mpm/mesh:node-active node)
    (with-accessors ((mass  node-mass)
                     (vel   node-velocity)
                     (disp   cl-mpm/mesh::node-displacment)
                     )
        node
      (declare (double-float mass))
      (progn
        (cl-mpm/fastmaths::fast-scale! vel (/ 1.0d0 mass))))))

(defgeneric update-node-kinematics (sim))
(defmethod update-node-kinematics ((sim mpm-sim))
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (iterate-over-nodes
     mesh
     (lambda (node)
       (calculate-kinematics node)))))

(defgeneric update-node-forces (sim)
  (:documentation "Update the acceleration from forces and apply any damping"))
(defmethod update-node-forces ((sim mpm-sim))
  (with-accessors ((damping sim-damping-factor)
                   (mass-scale sim-mass-scale)
                   (mesh sim-mesh)
                   (damping-algo sim-damping-algorithm)
                   (dt sim-dt))
      sim
    (iterate-over-nodes
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (calculate-forces node damping dt mass-scale))))))



(defmethod update-node-forces ((sim mpm-sim-quasi-static))
  (with-accessors ((damping sim-damping-factor)
                   (mass-scale sim-mass-scale)
                   (mesh sim-mesh)
                   (dt sim-dt))
      sim
    (iterate-over-nodes
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (calculate-forces-cundall node damping dt mass-scale))))))



(defun apply-bcs (mesh bcs dt)
  "Apply all normal bcs onto the mesh"
  (declare (cl-mpm/mesh::mesh mesh))
  (with-accessors ((nodes  mesh-nodes)
                   (nD     mesh-nD)
                   (mc     mesh-count)) mesh
                                        ;each bc is a list (pos value)
    (lparallel:pdotimes (i (length bcs))
      (let ((bc (aref bcs i)))
        (when bc
          (with-accessors ((node cl-mpm/bc::bc-node)
                           (index cl-mpm/bc::bc-index))
              bc
            (if (or node
                    (not index))
                (cl-mpm/bc:apply-bc bc node mesh dt)
                (progn
                  (setf node (cl-mpm/mesh:get-node mesh index))
                  (if node
                      (cl-mpm/bc:apply-bc bc node mesh dt)
                      (error "BC attempted to get a nil node ~A ~A" bc index))))))))))




(defun rotation-matrix (degrees)
  (let ((angle (/ (* pi degrees) 180)))
    (magicl:from-list (list (cos angle) (- (sin angle))
                            (sin angle) (cos angle))
                      '(2 2)
                      :type 'double-float
                      )))

(defun map-jacobian (mesh mp dt)
  (with-accessors ((stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (volume cl-mpm/particle:mp-volume)
                   (def cl-mpm/particle:mp-deformation-gradient))
      mp
    (let* ((df (cl-mpm/utils::matrix-from-list '(1d0 0d0 0d0
                                                 0d0 1d0 0d0
                                                 0d0 0d0 1d0))))
      (cl-mpm/fastmaths::fast-.+ df stretch-tensor df)
      (let ((j-inc (cl-mpm/fastmaths:det-3x3 df))
            (j-n
              1d0
              ;; (cl-mpm/fastmaths:det-3x3 def)
              ))
        (iterate-over-neighbours
         mesh mp
         (lambda (mesh mp node svp grads fsvp fgrads)
           (with-accessors ((node-active cl-mpm/mesh:node-active)
                            (node-j-inc cl-mpm/mesh::node-jacobian-inc)
                            (node-volume cl-mpm/mesh::node-volume)
                            (node-svp cl-mpm/mesh::node-svp-sum)
                            (node-lock cl-mpm/mesh::node-lock))
               node
             (declare (double-float node-j-inc svp volume j-inc j-n node-volume))
             (when node-active
               (sb-thread:with-mutex (node-lock)
                 ;(incf node-j-inc (* svp j-inc j-n (/ volume node-volume)))
                 (incf node-j-inc (* j-inc j-n (/ svp node-svp)))
                 )))))))))


(defun update-particle-kirchoff (mesh mp dt)
  (with-accessors ((disp cl-mpm/particle::mp-displacement)
                   (disp-inc cl-mpm/particle::mp-displacement-increment)
                   (pos cl-mpm/particle::mp-position)
                   (nc cl-mpm/particle::mp-cached-nodes)
                   (contact-step cl-mpm/particle::mp-penalty-contact-step)
                   (friction-force cl-mpm/particle:mp-penalty-frictional-force)
                   (normal-force cl-mpm/particle::mp-penalty-normal-force)
                   (contact cl-mpm/particle::mp-penalty-contact)
                   )
      mp
    (setf (fill-pointer nc) 0)
    (cl-mpm/fastmaths::fast-.+-vector pos disp-inc pos)
    (cl-mpm/fastmaths::fast-.+-vector disp disp-inc disp)
    ;; (break)
    (setf contact-step contact)
    (unless contact
      (cl-mpm/fastmaths:fast-zero friction-force)
      (setf normal-force 0d0))
    (setf contact nil)))

(defgeneric update-particle (mesh mp dt)
  (:documentation "End of step update"))
(defmethod update-particle (mesh (mp cl-mpm/particle:particle) dt)
  (update-particle-kirchoff mesh mp dt)
  (update-domain-corner mesh mp dt)
  ;; (update-domain-stretch mesh mp dt)
  (cl-mpm::scale-domain-size mesh mp)
  )

(declaim (notinline p2g))
(defun update-particles (sim)
  "Map particle momentum to the grid"
  (with-accessors ((mps sim-mps)
                   (mesh sim-mesh)
                   (dt sim-dt))
      sim
    (iterate-over-mps
     mps
     (lambda (mp)
       (update-particle mesh mp dt)))))



(declaim (ftype (function (cl-mpm/mesh::mesh
                           (array cl-mpm/particle:particle)
                           double-float
                           &optional boolean) (values)) update-stress))
(defun update-stress (mesh mps dt &optional (fbar nil))
  "Update all stresses, with optional f-bar"
  (declare ((array cl-mpm/particle:particle) mps) (cl-mpm/mesh::mesh mesh))
  ;; (iterate-over-mps
  ;;  mps
  ;;  (lambda (mp)
  ;;    (calculate-strain-rate mesh mp dt)
  ;;    (map-jacobian mesh mp dt)))
  (iterate-over-mps
   mps
   (lambda (mp)
     (update-stress-mp mesh mp dt fbar)
     ;; (post-stress-step mesh mp dt)
     ))
  (values))


(defgeneric update-stress-mp (mesh mp dt fbar)
  (:documentation "A mp dependent stress update scheme"))
(defmethod update-stress-mp (mesh (mp cl-mpm/particle::particle) dt fbar)
  (update-stress-kirchoff mesh mp dt fbar))

(defgeneric post-stress-step (mesh mp dt)
  (:documentation "Mp dependent step called after full stress update"))
(defmethod post-stress-step (mesh (mp cl-mpm/particle::particle) dt)
  (declare (ignore mesh mp dt)))

(declaim (notinline reset-grid))
(defun reset-grid (mesh)
  "Reset all nodes on the grid"
  (declare (cl-mpm/mesh::mesh mesh))
  (iterate-over-nodes
   mesh
   (lambda (node)
     (when (cl-mpm/mesh:node-active node)
       (cl-mpm/mesh:reset-node node)))))
(declaim (notinline reset-grid-velocity))
(defun reset-grid-velocity (mesh)
  "Reset all velocity map on grid for MUSL"
  (declare (cl-mpm/mesh::mesh mesh))
  (iterate-over-nodes
   mesh
   (lambda (node)
     (when (cl-mpm/mesh:node-active node)
       (cl-mpm/mesh::reset-node-velocity node)))))

(defun zero-grid-velocity (mesh)
  "Reset all velocity map on grid for MUSL"
  (declare (cl-mpm/mesh::mesh mesh))
  (iterate-over-nodes
   mesh
   (lambda (node)
     (when (cl-mpm/mesh:node-active node)
       (cl-mpm/fastmaths::fast-zero (cl-mpm/mesh::node-velocity node))))))


(defun filter-grid (mesh mass-thresh)
  "Filter out nodes with very small massess"
  (declare (cl-mpm/mesh::mesh mesh) (double-float mass-thresh))
  (iterate-over-nodes
   mesh
   (lambda (node)
     (when (and (cl-mpm/mesh:node-active node)
                (<= (cl-mpm/mesh:node-mass node) mass-thresh))
       (cl-mpm/mesh:reset-node node)))))

(defun filter-grid-velocity (mesh mass-thresh)
  "Filter out nodes with very small massess"
  (declare (cl-mpm/mesh::mesh mesh) (double-float mass-thresh))
  (iterate-over-nodes
   mesh
   (lambda (node)
     (when (and (cl-mpm/mesh:node-active node)
                (<= (cl-mpm/mesh:node-mass node) mass-thresh))
       (cl-mpm/mesh::reset-node-velocity node)))))

(defun filter-grid-volume (mesh volume-ratio)
  (declare (cl-mpm/mesh::mesh mesh) (double-float volume-ratio))
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
    (dotimes (i (array-total-size nodes))
      (let ((node (row-major-aref nodes i)))
        (when (and (cl-mpm/mesh:node-active node)
                   (< (the double-float (cl-mpm/mesh:node-mass node))
                      (* volume-ratio (the double-float (cl-mpm/mesh::node-volume-true node)))))
          (setf (cl-mpm/mesh::node-active node) nil)
          (cl-mpm/mesh:reset-node node)
          )))))

(defgeneric sim-add-mp (sim mp)
  (:documentation "A function to append an mp to a simulation"))

(defmethod sim-add-mp (sim mp)
  (with-accessors ((uid-counter sim-unique-index-counter)
                   (lock sim-unique-index-lock))
    sim
      (sb-thread:with-mutex (lock)
        (setf (cl-mpm/particle::mp-unique-index mp) uid-counter)
        (incf uid-counter))
    (vector-push-extend mp (cl-mpm:sim-mps sim))))

(defgeneric remove-mps-func (sim func)
  (:documentation "A function for removing mps from a sim"))

(defmethod remove-mps-func (sim func)
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (declare ((vector cl-mpm/particle::particle) mps))
    (when (> (length mps) 0)
      ;; (delete-mps-func mps func)
      (setf mps
            (delete-if func mps)))
    ;;Sometimes when compacting the array; sbcl will just discard make and unadjustable array in place which is a bit wild
    (when (and (not (adjustable-array-p mps))
               (= (length mps) 0))
      (setf mps (make-array 0 :adjustable t :fill-pointer 0)))
    )
  (values))

(defun remove-material-damaged (sim)
  "Remove material points that have strain energy density exceeding fracture toughness"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)
                   (instant-damage-removal cl-mpm::sim-mp-damage-removal-instant)
                   )
      sim
    (let ((h (cl-mpm/mesh:mesh-resolution mesh)))
      (remove-mps-func
       sim
       (lambda (mp)
         (with-accessors ((damage cl-mpm/particle:mp-damage)
                          (def cl-mpm/particle::mp-deformation-gradient))
             mp
           (and (>= damage 0.9d0)
                (or instant-damage-removal
                    (damage-removal-criteria mp h) )
                )))))
    ))

(defun damage-removal-criteria (mp h)
  "Some criteria for when we should remove fully damaged MPs"
  (with-accessors ((def cl-mpm/particle:mp-deformation-gradient)
                   (lens cl-mpm/particle::mp-domain-size)
                   (strain cl-mpm/particle::mp-strain)
                   (lens-0 cl-mpm/particle::mp-domain-size-0)
                   ;; (ft cl-mpm/particle::mp-initiation-stress)
                   ;; (ductility cl-mpm/particle::mp-ductility)
                   ;; (E cl-mpm/particle::mp-e)
                   )
      mp
    ;; (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix strain)))
    (let* ((l-factor 2.00d0)
           ;; (e0 (/ ft E))
           ;; (ef (/ (* ft (+ ductility 1d0)) (* 2d0 E)))
           ;; (ef (* ef 50d0))
           ;; (e_max (reduce #' max l))
           )
      (cond
        ;; ((< ef e_max) :x)
        ;; ((< ef (tref strain 0 0)) :x)
        ;; ((< ef (tref strain 1 0)) :y)
        ;; ((< ef (tref strain 2 0)) :z)
        ((< (* l-factor (varef lens-0 0)) (varef lens 0)) :x)
        ((< (* l-factor (varef lens-0 1)) (varef lens 1)) :y)
        ;; ((and (< (* l-factor (tref lens-0 2 0)) (tref lens 2 0))
        ;;       (> (tref lens-0 2 0) 0d0)
        ;;       ) :z)
        (t nil)
        ))))

(defun split-criteria (mp h)
  "Some numerical splitting estimates"
  (with-accessors ((def cl-mpm/particle:mp-deformation-gradient)
                   (lens cl-mpm/particle::mp-domain-size)
                   (lens-0 cl-mpm/particle::mp-domain-size-0)
                   (split-depth cl-mpm/particle::mp-split-depth)
                   )
      mp
    (when t;(< split-depth *max-split-depth*)
      (let ((l-factor 1.00d0)
            (h-factor (* 0.55d0 h))
            (s-factor 1.5d0))
        (cond
          ((< h-factor (varef lens 0)) :x)
          ((< h-factor (varef lens 1)) :y)
          ((< h-factor (varef lens 2)) :z)
          ;; ((< (* l-factor (tref lens-0 0 0)) (tref lens 0 0)) :x)
          ;; ((< (* l-factor (tref lens-0 1 0)) (tref lens 1 0)) :y)
          ;; ((< (* l-factor (tref lens-0 2 0)) (tref lens 2 0)) :z)
          ;; ((< l-factor (/ (tref lens 0 0) (tref lens-0 0 0))) t)
          ;; ((< l-factor (/ (tref lens 1 0) (tref lens-0 1 0))) t)
                                        ;((< 2.0d0 (tref def 1 1)) t)
                                        ;((> 1.5d0 (tref def 0 0)) t)
          (t nil)
          )))))

(defun single-particle-criteria (mesh mp)
  "Criteria for checking if material point is unconnected to another MP"
  (let ((svp-sum 0d0)
        (alone t))
    (iterate-over-neighbours
      mesh mp
      (lambda (mesh mp node svp grads fsvp fgrads)
        (declare
          (cl-mpm/mesh::node node)
          (cl-mpm/particle:particle mp)
          (type double-float svp))
        (with-accessors ((node-svp cl-mpm/mesh::node-svp-sum)
                         (node-active cl-mpm/mesh:node-active)
                         ) node
          (when node-active
            (incf svp-sum node-svp)))))
    (and (< svp-sum 2d0) (not (= svp-sum 0d0)))
    ;; (setf alone t)
    ;; alone
    ))
(defun gimp-removal-criteria (mp h)
  "Criteria for removal of gimp mps based on domain length"
  (with-accessors ((lens cl-mpm/particle::mp-domain-size))
      mp
    (declare (double-float h))
    (let ((h-factor
            h
            ;; (* 0.7d0 h)
                    )
          (aspect 0.01d0))
      (cond
        ((< h-factor (the double-float (varef lens 0))) :x)
        ((< h-factor (the double-float (varef lens 1))) :y)
        ((< h-factor (the double-float (varef lens 2))) :z)
        ;; ((> aspect (/ (the double-float (varef lens 0))
        ;;               (the double-float (varef lens 1))
        ;;               )) :x)
        ;; ((> aspect (/ (the double-float (varef lens 1))
        ;;               (the double-float (varef lens 0))
        ;;               )) :x)
        (t nil)
        ))))

(defgeneric calculate-min-dt-mps (sim)
  (:documentation "A function for calculating an approximate stable timestep based on MPs"))
(defmethod calculate-min-dt-mps ((sim mpm-sim))
  "Estimate minimum p-wave modulus"
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mass-scale cl-mpm::sim-mass-scale))
      sim
    (let ((inner-factor most-positive-double-float))
      (declare (double-float inner-factor mass-scale))
      (setf inner-factor
            (reduce-over-nodes
             mesh
             (lambda (node)
               (if (and (cl-mpm/mesh::node-active node)
                        (not (cl-mpm/mesh::node-agg node)))
                   (with-accessors ((node-active  cl-mpm/mesh:node-active)
                                    (pmod cl-mpm/mesh::node-pwave)
                                    (mass cl-mpm/mesh::node-mass)
                                    (svp-sum cl-mpm/mesh::node-svp-sum)
                                    (vol cl-mpm/mesh::node-volume)
                                    (vel cl-mpm/mesh::node-velocity)
                                    ) node
                     (if (and (> vol 0d0)
                              (> pmod 0d0)
                              (> svp-sum 0d0))
                         (let ((nf (+ (/ mass (* vol (+ (/ pmod svp-sum)))))))
                           nf)
                         most-positive-double-float))
                   most-positive-double-float))
             #'min))
      (if (< inner-factor most-positive-double-float)
          (* (sqrt mass-scale) (sqrt inner-factor) (cl-mpm/mesh:mesh-resolution mesh))
          (cl-mpm:sim-dt sim)))))

(defgeneric calculate-min-dt-bcs (sim))
(defmethod calculate-min-dt-bcs (sim)
  (let ((min-dt nil))
    (cl-mpm::iterate-over-bcs-force-serial
     sim
     (lambda (bc)
       (let ((dt-req (cl-mpm/bc::calculate-min-dt-bc sim bc)))
         (pprint dt-req)
         (when (and
                dt-req
                (if min-dt
                 (< dt-req min-dt)
                 t))
           (setf min-dt dt-req)))))
    min-dt))


(defgeneric calculate-min-dt (sim)
  (:documentation "A function for calculating an approximate stable timestep"))
(defmethod calculate-min-dt (sim)
  ;; (let ((dt-req (calculate-min-dt-mps sim)))
  ;;   (let ((bc-min-dt (calculate-min-dt-bcs sim)))
  ;;     (when bc-min-dt
  ;;       (setf dt-req (min dt-req bc-min-dt))))
  ;;   dt-req)
  ;; (calculate-min-dt-bcs sim)
  (calculate-min-dt-mps sim)
  )

(defun calculate-adaptive-time (sim target-time &key (dt-scale 1d0))
  "Given that the system has previously been mapped to, caluclate an estimated dt and substep for a given target time
The value dt-scale allows for the estimated dt to be scaled up or down
This modifies the dt of the simulation in the process
"
  (let* ((dt-e (* dt-scale (calculate-min-dt sim)))
         (substeps-e (floor target-time dt-e)))
    ;; (format t "CFL dt estimate: ~f~%" dt-e)
    ;; (format t "CFL step count estimate: ~D~%" substeps-e)
    (setf (cl-mpm:sim-dt sim) dt-e)
    ;; (setf substeps substeps-e)
    (values dt-e substeps-e)
    )
  )


(defgeneric add-mps (sim mps-array))
(defmethod add-mps (sim mps-array)
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
      (if (> (length mps) 0)
          (progn
            (loop for mp across mps-array
                  do (sim-add-mp sim mp)))
          (progn
            (setf (cl-mpm:sim-mps sim) (make-array (length mps-array) :adjustable t :fill-pointer 0))
            (loop for mp across mps-array
                  do (sim-add-mp sim mp))))))

(defun add-bcs (sim bcs-array)
  "Add nodal essential bcs"
  (with-accessors ((bcs cl-mpm:sim-bcs))
      sim
    (if (> (length bcs) 0)
        (progn
          (loop for bc across bcs-array
                do (vector-push-extend bc bcs)))
        (setf bcs bcs-array))))
(defun add-bcs-force-list (sim new-bcs)
  "Add bcs that apply forces, ordered in a FILO"
  (with-accessors ((bcs-force-list cl-mpm:sim-bcs-force-list))
      sim
    ;; (when (= (length bcs-force-list 0)))
    (typecase new-bcs
      (list
       ;;Common case, we have a list, or arrays of BCs, running bcs in serial down the list, but in parrallel over the array
       (setf bcs-force-list (nconc new-bcs bcs-force-list)))
      (array
       ;;We have an array of bcs, that we want to add to the list -> add it to the back in its own list
       (setf bcs-force-list (nconc (list new-bcs) bcs-force-list)))
      (cl-mpm/bc::bc
       ;;Uncommon case we have a singluar bc we wish to add, so double wrap it
       (setf bcs-force-list (nconc (list (cl-mpm/bc:make-bcs-from-list (list new-bcs))) bcs-force-list)))
      (t (error "BCs must be a list, array or singluar bc")))))

(defun get-mp (sim index)
  (aref (sim-mps sim) index))

(defmethod update-sim ((sim mpm-sim-sd))
  "Update stress first algorithm"
  (with-slots (
               (mesh mesh)
               (mesh-p mesh-p)
               (mps mps)
               (bcs bcs)
               (bcs-p bcs-p)
               (bcs-force bcs-force)
               (bcs-force-list bcs-force-list)
               (dt dt)
               (mass-filter mass-filter)
               (split allow-mp-split)
               (enable-damage enable-damage)
               (nonlocal-damage nonlocal-damage)
               (remove-damage allow-mp-damage-removal)
               (fbar enable-fbar)
               (time time)
               (vel-algo velocity-algorithm)
               )
                sim
    (declare (double-float mass-filter dt time))
    (progn
      (reset-grid mesh)
      (reset-grid mesh-p)
      (when (> (length mps) 0)
        (p2g mesh mps)
        (p2g mesh-p mps)
        (when (> mass-filter 0d0)
          (filter-grid mesh (sim-mass-filter sim))
          (filter-grid mesh-p (sim-mass-filter sim)))
        (update-node-kinematics mesh dt)
        (update-node-kinematics mesh-p dt)
        (apply-bcs mesh bcs dt)
        (apply-bcs mesh-p bcs-p dt)
        (update-stress-p mesh-p mps dt)
        (update-stress mesh mps dt fbar)
        (p2g-force mesh mps)
        (loop for bcs-f in bcs-force-list
              do (apply-bcs mesh bcs-f dt))

        (update-node-forces sim)
        ;; ;; Reapply velocity BCs
        (apply-bcs mesh bcs dt)
        ;; (apply-bcs mesh-p bcs-p dt)
        ;; Also updates mps inline
        (g2p mesh mps dt vel-algo)
        (when split
          (split-mps sim))
        (check-mps sim)
        (check-single-mps sim)
        )
      (incf time dt))))


(defgeneric update-dynamic-stats (sim))
(defmethod update-dynamic-stats (sim))

(defmethod update-node-forces ((sim mpm-sim-sd))
  (with-accessors ((damping sim-damping-factor)
                   (mass-scale sim-mass-scale)
                   (mesh sim-mesh)
                   (mesh-p sim-mesh-p)
                   (dt sim-dt))
      sim
    (iterate-over-nodes
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (calculate-forces node damping dt mass-scale))))))


(defgeneric finalise-loadstep (sim)
  )
(defmethod finalise-loadstep ((sim mpm-sim))
  ;;Standard MPM do not need to do anything fancy for finalisation as we set each step as a new loadstep anyway
  ;; (incf (sim-time sim) (cl-mpm:sim-dt sim))
  )

(defgeneric reset-loadstep (sim)
  )
(defmethod reset-loadstep ((sim mpm-sim))
  (cl-mpm:iterate-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (cl-mpm/particle::reset-loadstep-mp mp)))
  ;; (when (cl-mpm::sim-allow-mp-split sim)
  ;;   (split-mps sim))
  ;; (check-mps sim)
  (reset-node-displacement sim)
  )

(defgeneric new-loadstep (sim)
  )

(defun reset-node-displacement (sim)
  (cl-mpm::iterate-over-cells
   (cl-mpm:sim-mesh sim)
   (lambda (cell)
     (cl-mpm/fastmaths:fast-zero (cl-mpm/mesh::cell-displacement cell))))
  (cl-mpm:iterate-over-nodes
   (cl-mpm:sim-mesh sim)
   (lambda (node)
     (cl-mpm/fastmaths:fast-zero (cl-mpm/mesh::node-displacment node)))))

(defmethod new-loadstep ((sim mpm-sim))
  (update-particles sim)
  (cl-mpm:iterate-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (cl-mpm/particle::new-loadstep-mp mp)))
  (when (cl-mpm::sim-allow-mp-split sim)
    (split-mps sim))
  (check-mps sim)
  (reset-node-displacement sim)
  )

(defun gradient-push-forwards (grads df)
  (let* ((grads-vec (magicl:linear-solve df (cl-mpm/utils:vector-from-list grads)))
         (grads (list (varef grads-vec 0)
                      (varef grads-vec 1)
                      (varef grads-vec 2))))
    grads))
(defun gradient-push-forwards-cached (grads df-inv)
  (let* ((grads-vec
           (magicl:@ (magicl:transpose! (cl-mpm/utils:vector-from-list grads)) df-inv)
           ;; (cl-mpm/fastmaths::fast-@-matrix-vector df-inv (cl-mpm/utils:vector-from-list grads))
                    )
         (grads (list (varef grads-vec 0)
                      (varef grads-vec 1)
                      (varef grads-vec 2))))
    grads))


