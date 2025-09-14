(defpackage :cl-mpm/implicit
  (:use :cl
   :cl-mpm/utils)
  (:export
   #:mpm-sim-implicit))

(in-package :cl-mpm/implicit)
(declaim (optimize (debug 0) (safety 0) (speed 3)))


(in-package :cl-mpm/mesh)
(defclass cl-mpm/mesh::node-implicit (node)
  ((displacement
    :accessor node-displacement
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils:vector-zeros))))

(in-package :cl-mpm/implicit)
(defclass mpm-sim-implicit (cl-mpm::mpm-sim)
  ()
  (:documentation "Explicit simulation with update stress first update"))



(defmethod update-sim ((sim mpm-sim-implicit))
  "Update stress first algorithm"
  (declare (cl-mpm::mpm-sim-usf sim))
  (with-slots ((mesh mesh)
               (mps mps)
               (bcs bcs)
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
               (gravity gravity)
               )
                sim
    (declare (double-float mass-filter dt time))
    (progn
      (reset-grid mesh)
      (when (> (length mps) 0)


        (p2g mesh mps)
        (when (> mass-filter 0d0)
          (filter-grid mesh (sim-mass-filter sim)))
        (update-node-kinematics mesh dt )
        (apply-bcs mesh bcs dt)
        (update-stress mesh mps dt fbar)
        ;; Map forces onto nodes
        (p2g-force sim)
        (loop for bcs-f in bcs-force-list
              do (apply-bcs mesh bcs-f dt))
        (update-node-forces sim)
        ;; Reapply velocity BCs
        (apply-bcs mesh bcs dt)
        ;; Also updates mps inline
        (g2p mesh mps dt vel-algo)
        (when split
          (split-mps sim)))
      ;; (check-mps sim)
      (incf time dt))))


(defun det-external-force (sim))

(defun det-mp-info (sim)
  )
(defun lin-solve (sim))
(defun update-kinematics (sim))
(defun det-mps (sim))
(defun compute-residuals (sim))
(defun update-mps (sim)
  )
(defun update-mps-g2p (mesh mp dt)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp)
           (double-float dt))
  "Map one MP from the grid"
  (with-accessors ((vel mp-velocity)
                   (pos mp-position)
                   (disp cl-mpm/particle::mp-displacement)
                   (acc cl-mpm/particle::mp-acceleration)
                   (nc cl-mpm/particle::mp-cached-nodes))
      mp
    (let* ((mapped-disp (cl-mpm/utils:vector-zeros)))
      (cl-mpm/fastmaths::fast-zero acc)
      (iterate-over-neighbours
       mesh mp
       (lambda (mesh mp node svp grads fsvp fgrads)
         (declare
          (ignore mp mesh fsvp fgrads)
          (cl-mpm/mesh::node node)
          (cl-mpm/particle:particle mp)
          (double-float svp))
         (with-accessors ((node-disp cl-mpm/mesh::node-displacement)
                          (node-active cl-mpm/mesh:node-active))
             node
           (declare (boolean node-active))
           (when node-active
             (cl-mpm/fastmaths::fast-fmacc mapped-disp node-disp svp)))))
      (progn
        ;;Invalidate shapefunction/gradient cache
        (update-particle mesh mp dt)
        (setf (fill-pointer nc) 0)
        (cl-mpm/fastmaths::fast-.+-vector pos mapped-disp pos)
        (cl-mpm/fastmaths::fast-.+-vector disp mapped-disp disp)))))
