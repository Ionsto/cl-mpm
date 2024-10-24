(defpackage :cl-mpm/penalty
  (:use :cl
        :cl-mpm
        :cl-mpm/particle
        :cl-mpm/mesh
        :cl-mpm/utils
        :cl-mpm/fastmaths
        )
  (:import-from
    :magicl tref .+ .-)
  (:export
   #:make-bc-penalty
   #:bc-penalty-friction))
(declaim (optimize (debug 0) (safety 0) (speed 3)))
;; (declaim (optimize (debug 3) (safety 3) (speed 0)))
(in-package :cl-mpm/penalty)

(defun ssqrt (a)
  (* (signum a) (sqrt (abs a))))

(declaim
 (ftype (function
         (magicl::matrix/double-float
          double-float
          magicl::matrix/double-float)
         double-float)
        penetration-distance-point))
(defun penetration-distance-point (point datum normal)
  "Get linear penetration distance"
  (let* ((ypos (cl-mpm/fastmaths::dot point normal))
         (dist (- datum ypos)))
    (the double-float dist)))

(declaim
 (ftype (function
         (cl-mpm/particle::particle
          double-float
          magicl::matrix/double-float)
         double-float)
        penetration-distance))
(defun penetration-distance (mp datum normal)
  "Get linear penetration distance"
  (let* ((ypos (cl-mpm/fastmaths::dot (cl-mpm/particle:mp-position mp) normal))
         (yheight (cl-mpm/fastmaths::dot
                   (magicl:scale
                    (cl-mpm/particle::mp-domain-size mp) 0.5d0)
                   (cl-mpm/fastmaths::norm
                    (magicl:.* normal normal))))
         (dist (- datum (- ypos yheight))))
    (the double-float dist)))

(defun penetration-point (mp pen datum normal)
  "Get linear penetration distance"
  (let* ((pos (cl-mpm/particle:mp-position mp))
         (domain (cl-mpm/particle::mp-domain-size mp)))
    (magicl:.- pos
               (magicl:.* normal (magicl:scale domain 0.5d0))
               )))

(defun trial-corner (mp normal)
  "Get linear penetration distance"
  (let* ((corner (cl-mpm/fastmaths:fast-.+
                  (cl-mpm/particle:mp-position mp)
                  (cl-mpm/fastmaths:fast-.*
                   (cl-mpm/fastmaths::fast-scale-vector (cl-mpm/particle::mp-domain-size mp) -0.5d0)
                   normal))))
    corner))

(defparameter *debug-mutex* (sb-thread:make-mutex))
(defparameter *debug-force* 0d0)
(defparameter *debug-force-count* 0)
(defparameter *debug-force-mp-count* 0)
;;Only vertical condition
(declaim (notinline apply-force-mps))
(defun apply-force-mps (mesh mps dt normal datum epsilon friction &optional (func-clip (lambda (mp) t))
                                                                    &key (damping 0.1d0))
  "Update force on nodes, with virtual stress field from mps"
  ;;If we lose contact we need to zero out our friction force
  (with-accessors ((nd cl-mpm/mesh::mesh-nd))
      mesh
    (cl-mpm:iterate-over-mps
     mps
     (lambda (mp)
       (let* ((penetration-dist (penetration-distance mp datum normal)))
         (declare (double-float penetration-dist))
         (when (> penetration-dist 0d0)
             (progn
               ;;Contact
               (with-accessors ((volume cl-mpm/particle:mp-volume)
                                (pressure cl-mpm/particle::mp-pressure)
                                (mp-vel cl-mpm/particle::mp-velocity)
                                (mp-mass cl-mpm/particle::mp-mass)
                                (mp-contact cl-mpm/particle::mp-penalty-contact)
                                (mp-friction cl-mpm/particle::mp-penalty-frictional-force)
                                (mp-normal-force cl-mpm/particle::mp-penalty-normal-force)
                                )
                   mp
                 (let* ((pen-point (penetration-point mp penetration-dist datum normal))
                        (normal-force (* (expt penetration-dist 1d0)
                                         epsilon
                                         ;; (expt volume (/ (- nd 1) nd))
                                         (expt volume (/ (- nd 1) nd))
                                         )))
                   (when (funcall func-clip
                                  (trial-corner mp normal)
                                  )
                     (sb-thread:with-mutex (*debug-mutex*)
                       (incf *debug-force* (* normal-force 1d0))
                       (incf *debug-force-count* 1))
                     (setf mp-contact t)
                     ;;Iterate over neighbour nodes
                     (let* ((force (cl-mpm/utils:vector-zeros))
                            (rel-vel (cl-mpm/fastmaths:dot normal mp-vel))
                            (tang-vel (cl-mpm/fastmaths:fast-.- mp-vel (magicl:scale normal rel-vel)))
                            (tang-vel-norm-squared (cl-mpm/fastmaths::mag-squared tang-vel))
                            ;; (normal-damping (* damping (sqrt (/ epsilon (/ mp-mass volume)))))
                            (normal-damping (* (/ pi 2) damping (sqrt (* epsilon mp-mass))))
                            ;; (normal-damping (* (/ pi 2) damping (sqrt epsilon)))
                            (damping-force (* normal-damping rel-vel))
                            (force-friction mp-friction)
                            (stick-friction (* friction (abs normal-force))))
                       ;; update trial frictional force
                       (when (> friction 0d0)
                         (when (> tang-vel-norm-squared 0d0)
                           ;; We have sliding behaviour
                           (let* (;(tang-vel-norm (sqrt tang-vel-norm-squared))
                                  ;; (tang-normal (cl-mpm/fastmaths:norm tang-vel))
                                  ;;Trial friction
                                        ;(trial-friction-force (* (/ epsilon 2d0) tang-vel-norm dt))
                                  )
                             (cl-mpm/fastmaths::fast-fmacc
                              force-friction
                              tang-vel
                              (* -1d0 (/ epsilon 2d0) dt))))
                         (when (> (cl-mpm/fastmaths::mag-squared force-friction) 0d0)
                           (if (> (cl-mpm/fastmaths::mag force-friction) stick-friction)
                               (progn
                                 ;; (cl-mpm/fastmaths::fast-scale! force-friction
                                 ;;                              (/ (cl-mpm/fastmaths::mag force-friction)
                                 ;;                                 stick-friction))
                                 (setf force-friction
                                       (magicl:scale
                                        (cl-mpm/fastmaths:norm force-friction)
                                        stick-friction))
                                 (setf (cl-mpm/particle::mp-penalty-friction-stick mp) t)
                                 )
                               (progn
                                 (setf (cl-mpm/particle::mp-penalty-friction-stick mp) nil)
                                 ))
                           (cl-mpm/fastmaths::fast-.+ force force-friction force))
                         (setf mp-friction force-friction))
                       (setf mp-normal-force (- normal-force damping-force))
                       (cl-mpm/fastmaths::fast-fmacc force
                                                    normal
                                                    (- normal-force damping-force))
                       (cl-mpm::iterate-over-neighbours-point-linear-3d
                        mesh
                        pen-point
                        (lambda (mesh node svp grads)
                          (with-accessors ((node-force cl-mpm/mesh:node-force)
                                           (node-ext-force cl-mpm/mesh::node-external-force)
                                           (node-lock  cl-mpm/mesh:node-lock)
                                           (node-vel  cl-mpm/mesh:node-velocity)
                                           (node-mass  cl-mpm/mesh:node-mass)
                                           (node-boundary cl-mpm/mesh::node-boundary-node)
                                           (node-boundary-scalar cl-mpm/mesh::node-boundary-scalar)
                                           (node-active  cl-mpm/mesh:node-active))
                              node
                            (declare (double-float volume svp))
                            ;;Lock node for multithreading
                            (when node-active
                              ;;Lock node for multithreading
                              (sb-thread:with-mutex (node-lock)
                                ;; (let* ((svp (* svp (expt volume (/ (- nd 1) nd))))))

                                ;; (cl-mpm/fastmaths::fast-fmacc node-force
                                ;;                              force
                                ;;                              svp)
                                (cl-mpm/fastmaths::fast-fmacc node-ext-force
                                                             force
                                                             svp)
                                ))))))))))))))))




(defun collect-contact-points (mesh mps normal datum)
  (loop for mp across mps
        when (> (penetration-distance mp datum normal) 0d0)
          collect
          (let* ((penetration-dist (penetration-distance mp datum normal)))
            (declare (double-float penetration-dist))
            (when (> penetration-dist 0d0)
              (with-accessors ((volume cl-mpm/particle:mp-volume)
                               (pressure cl-mpm/particle::mp-pressure)
                               (mp-vel cl-mpm/particle::mp-velocity)
                               )
             mp
           (penetration-point mp penetration-dist datum normal))))))
(defun collect-contact-points-bc (mesh mps bc)
  (with-accessors ((normal bc-penalty-normal)
                   (datum bc-penalty-datum)
                   )
      bc
    (collect-contact-points mesh mps normal datum)))


(defclass bc-penalty (cl-mpm/bc::bc)
  ((normal
    :accessor bc-penalty-normal
    :initarg :normal
    :initform (cl-mpm/utils:vector-zeros))
   (datum
    :accessor bc-penalty-datum
    :initarg :datum
    :initform 0d0)
   (sim
    :accessor bc-penalty-sim
    :initarg :sim)
   (friction
    :accessor bc-penalty-friction
    :initform 0d0
    :initarg :friction)
   (epsilon
    :accessor bc-penalty-epsilon
    :initform 0d0
    :initarg :epsilon)
   (damping
    :accessor bc-penalty-damping
    :initform 0d0
    :initarg :damping)
   (contact-points
    :accessor bc-penalty-contact-points
    :initform (list))
   (load
    :accessor bc-penalty-load
    :initform 0d0)
   (load-lock
    :accessor bc-penalty-load-lock
    :initform (sb-thread:make-mutex)))
  (:documentation "A nonconforming neumann bc"))

(defclass bc-penalty-distance (bc-penalty)
  ((center-point
    :accessor bc-penalty-distance-center-point
    :initarg :center-point
    :initform (cl-mpm/utils:vector-zeros))
   (radius
    :accessor bc-penalty-distance-radius
    :initarg :radius
    :initform 0d0)))

(defun make-bc-penalty-distance-point (sim normal point radius epsilon friction damping)
  (let* ((normal (cl-mpm/fastmaths::norm normal))
         (datum (- (penetration-distance-point point 0d0 normal))))
    (make-instance 'bc-penalty-distance
                   :index nil
                   :sim sim
                   :datum datum
                   :normal normal
                   :epsilon epsilon
                   :friction friction
                   :center-point point
                   :radius radius)))

(defun 2d-orthog (vec)
  (cl-mpm/utils:vector-from-list (list (- (varef vec 1)) (varef vec 0) 0d0)))

(defun make-bc-penalty-line-segment (sim point-a point-b epsilon friction damping)
  (let* ((diff-a (cl-mpm/fastmaths:fast-.- point-b point-a))
         (point (cl-mpm/fastmaths:fast-scale! (cl-mpm/fastmaths:fast-.+ point-a point-b) 0.5d0))
         (radius (* 0.5d0 (cl-mpm/fastmaths:mag diff-a)))
         (normal (cl-mpm/fastmaths:norm (2d-orthog diff-a)))
         (datum (- (penetration-distance-point point 0d0 normal))))
    (make-instance 'bc-penalty-distance
                   :index nil
                   :sim sim
                   :datum datum
                   :normal normal
                   :epsilon epsilon
                   :friction friction
                   :center-point point
                   :radius radius)))


(defun bc-set-center (bc new-center)
  (with-accessors ((normal bc-penalty-normal)
                   (datum bc-penalty-datum)
                   (center-point bc-penalty-distance-center-point))
      bc
    (let ((new-datum (- (penetration-distance-point new-center 0d0 normal))))
      (setf center-point new-center
            datum new-datum))))
(defun bc-increment-center (bc delta-center)
  (with-accessors ((normal bc-penalty-normal)
                   (datum bc-penalty-datum)
                   (center-point bc-penalty-distance-center-point))
      bc
    (let* ((new-center (cl-mpm/fastmaths:fast-.+ center-point delta-center))
           (new-datum (- (penetration-distance-point new-center 0d0 normal))))
      (setf center-point new-center
            datum new-datum))))


(defgeneric penalty-contact-valid (bc point))

(defmethod penalty-contact-valid ((bc bc-penalty) point)
  t)
(defmethod penalty-contact-valid ((bc bc-penalty-distance) point)
  (let* ((normal (bc-penalty-normal bc))
         (center-point (bc-penalty-distance-center-point bc))
         (p-d (cl-mpm/fastmaths:dot point normal))
         (c-d (cl-mpm/fastmaths:dot center-point normal))
         (diff
          (cl-mpm/fastmaths::fast-.-
           (cl-mpm/fastmaths:fast-.- point (cl-mpm/fastmaths:fast-scale-vector normal p-d))
           (cl-mpm/fastmaths:fast-.- center-point (cl-mpm/fastmaths:fast-scale-vector normal c-d)))))
    (<=
     (cl-mpm/fastmaths::mag
      diff)
     (bc-penalty-distance-radius bc))))



(defun make-bc-penalty (sim datum epsilon friction)
  (let ((normal (cl-mpm/fastmaths::norm (cl-mpm/utils:vector-zeros))))
    (setf (magicl:tref normal (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)) 0) 1d0)
    (setf normal (cl-mpm/fastmaths::norm normal))
    (make-instance 'bc-penalty
                   :index nil
                   :sim sim
                   :datum datum
                   :normal normal
                   :epsilon epsilon
                   :friction friction)))

(defun make-bc-penalty-point-normal (sim normal point epsilon friction &optional (damping 0d0))
  (let* ((normal (cl-mpm/fastmaths::norm normal))
         ;(datum (cl-mpm/fastmaths::dot normal point))
         (datum (- (penetration-distance-point point 0d0 normal)))
         )
    ;; (format t "Normal ~F ~F ~%" (magicl:tref normal 0 0) (magicl:tref normal 1 0))
    (make-instance 'bc-penalty
                   :index nil
                   :sim sim
                   :datum datum
                   :normal normal
                   :epsilon epsilon
                   :friction friction
                   :damping damping
                   )))

(defun disp-distance (mp datum normal)
  "Get linear penetration distance"
  (let* ((ypos (cl-mpm/fastmaths::dot (cl-mpm/particle:mp-position mp) normal))
         (yheight (cl-mpm/fastmaths::dot (magicl:scale (cl-mpm/particle::mp-domain-size mp) 0.5d0)
                                        (cl-mpm/fastmaths::norm (magicl:.* normal normal))))
         (dist (- datum (- ypos yheight))))
    (the double-float dist)))

(defun disp-penetration-point (mp pen datum normal)
  "Get linear penetration distance"
  (let* ((pos (cl-mpm/particle:mp-position mp))
         (domain (cl-mpm/particle::mp-domain-size mp)))
    (magicl:.- pos
               (magicl:.* normal (magicl:scale domain 0.5d0)))))

(defun apply-displacement-control-mps (mesh mps dt normal datum epsilon friction)
  "Update force on nodes, with virtual stress field from mps"
  (declare (double-float datum dt epsilon friction))
  (cl-mpm:iterate-over-mps
   mps
   (lambda (mp)
     (let* ((penetration-dist (disp-distance mp datum normal)))
       (declare (double-float penetration-dist))
       (when t;(> (abs penetration-dist) 0d0)
         (with-accessors ((volume cl-mpm/particle:mp-volume)
                          (pressure cl-mpm/particle::mp-pressure)
                          (mp-vel cl-mpm/particle::mp-velocity)
                          (mp-mass cl-mpm/particle::mp-mass)
                          )
             mp
           (let* ((pen-point (disp-penetration-point mp penetration-dist datum normal))
                  (h (cl-mpm/mesh:mesh-resolution mesh))
                  (nd (cl-mpm/mesh:mesh-nd mesh))
                  (normal-force (* (signum penetration-dist)
                                   (expt (abs penetration-dist) 1d0)
                                   epsilon
                                   ;; h
                                   (expt volume (/ (- nd 1) nd))
                                   ;(expt h (/ 1 nd))
                                   )))
             (sb-thread:with-mutex (*debug-mutex*)
               (incf *debug-force* (* normal-force 1d0))
               (incf *debug-force-count* 1))
             ;;Iterate over neighbour nodes
             (cl-mpm::iterate-over-neighbours-point-linear
              mesh pen-point
              (lambda (mesh node svp grads)
                (when (cl-mpm/mesh:node-active node)
                  (with-accessors ((node-force cl-mpm/mesh:node-force)
                                   (node-ext-force cl-mpm/mesh::node-external-force)
                                   (node-lock  cl-mpm/mesh:node-lock)
                                   (node-vel  cl-mpm/mesh:node-velocity)
                                   (node-mass  cl-mpm/mesh:node-mass)
                                   (node-boundary cl-mpm/mesh::node-boundary-node)
                                   (node-boundary-scalar cl-mpm/mesh::node-boundary-scalar))
                      node
                    (declare (double-float volume svp))
                    ;;Lock node for multithreading
                    (let* ((force (cl-mpm/utils:vector-zeros))
                           (rel-vel (cl-mpm/fastmaths::dot normal mp-vel))
                           ;; (svp (* svp 1))
                           ;; (normal-force (* dt normal-force))
                           (normal-damping 0d0)
                           (damping-force 0d0
                                          ;(* normal-damping rel-vel)
                                          ))
                      (declare (double-float rel-vel normal-damping damping-force))

                      (cl-mpm/fastmaths::fast-fmacc force
                                                   normal
                                                   (- normal-force
                                                      damping-force))

                      (sb-thread:with-mutex (node-lock)
                        ;; (cl-mpm/fastmaths::fast-fmacc
                        ;;  node-force
                        ;;  force
                        ;;  svp)
                        (cl-mpm/fastmaths::fast-fmacc
                         node-ext-force
                         force
                         svp))))))))))))))


(defstruct line-segment
  normal
  datum
  point-start
  point-end)

(defclass bc-penalty-structure (cl-mpm/bc::bc)
  ((sim
    :accessor bc-penalty-sim
    :initarg :sim)
   (friction
    :accessor bc-penalty-friction
    :initarg :friction)
   (epsilon
    :accessor bc-penalty-epsilon
    :initarg :epsilon)
   (damping
    :accessor bc-penalty-damping
    :initarg :damping)
   (load
    :accessor bc-penalty-load
    :initform 0d0)
   (load-lock
    :accessor bc-penalty-load-lock
    :initform (sb-thread:make-mutex))
   (sub-bcs
    :accessor bc-penalty-structure-sub-bcs
    :initarg :sub-bcs
    :initform (list))
   (contact-points
    :accessor bc-penalty-structure-contact-points
    :initform (list)))
  (:documentation "A single multi-surface structure that should resolve contact through a closest-point algorithm"))

(defun make-bc-penalty-structure (sim epsilon friction damping sub-bcs)
  (make-instance 'bc-penalty-structure
                  :index nil
                  :sim sim
                  :epsilon epsilon
                  :friction friction
                  :damping damping
                  :sub-bcs sub-bcs))

(defstruct contact
  point
  datum
  penetration
  normal
  sub-bc
  )


;;
(defun apply-penalty-point (mesh bc mp point dt)
  (with-accessors ((nd cl-mpm/mesh:mesh-nd))
      mesh
    (with-accessors ((datum bc-penalty-datum)
                     (damping bc-penalty-damping)
                     (epsilon bc-penalty-epsilon)
                     (friction bc-penalty-friction)
                     (normal bc-penalty-normal)
                     (debug-mutex bc-penalty-load-lock)
                     (debug-load bc-penalty-load))
        bc
      (with-accessors ((volume cl-mpm/particle:mp-volume)
                       (pressure cl-mpm/particle::mp-pressure)
                       (mp-vel cl-mpm/particle::mp-velocity)
                       (mp-mass cl-mpm/particle::mp-mass)
                       (mp-contact cl-mpm/particle::mp-penalty-contact)
                       (mp-friction cl-mpm/particle::mp-penalty-frictional-force)
                       (mp-normal-force cl-mpm/particle::mp-penalty-normal-force))
          mp
        (let* ((penetration (penetration-distance-point point datum normal))
               (normal-force (* (expt penetration 1d0)
                                epsilon
                                (expt volume (/ (- nd 1) nd)))))
          (sb-thread:with-mutex (debug-mutex)
            (incf debug-load normal-force))
          (setf mp-contact t)
          (let* ((force (cl-mpm/utils:vector-zeros))
                 (rel-vel (cl-mpm/fastmaths:dot normal mp-vel))
                 (tang-vel (cl-mpm/fastmaths:fast-.- mp-vel (magicl:scale normal rel-vel)))
                 (tang-vel-norm-squared (cl-mpm/fastmaths::mag-squared tang-vel))
                 ;(normal-damping (* (/ pi 2) damping (sqrt (* epsilon mp-mass))))
                 (normal-damping (* 2d0 damping (sqrt (* epsilon mp-mass))))
                 (damping-force (* normal-damping rel-vel))
                 (force-friction mp-friction)
                 (stick-friction (* friction (abs normal-force))))
            ;; update trial frictional force

            (when (> friction 0d0)
              (when (> tang-vel-norm-squared 0d0)
                (cl-mpm/fastmaths::fast-fmacc
                 force-friction
                 tang-vel
                 (* -1d0 (/ epsilon 2d0) dt)
                 ;; (* -1d0 epsilon dt)
                 ))

              ;; (break)
              (when (> (cl-mpm/fastmaths::mag-squared force-friction) 0d0)
                ;; (break)
                (if (> (cl-mpm/fastmaths::mag force-friction) stick-friction)
                    (progn
                      (setf force-friction
                            (magicl:scale
                             (cl-mpm/fastmaths:norm force-friction)
                             stick-friction))
                      (setf (cl-mpm/particle::mp-penalty-friction-stick mp) t))
                    (progn
                      (setf (cl-mpm/particle::mp-penalty-friction-stick mp) nil)))
                (cl-mpm/fastmaths::fast-.+ force force-friction force))
              (setf mp-friction force-friction))
            (setf mp-normal-force (- normal-force damping-force))
            (cl-mpm/fastmaths::fast-fmacc force
                                         normal
                                         (- normal-force damping-force))
            (cl-mpm::iterate-over-neighbours-point-linear
             mesh
             point
             (lambda (mesh node svp grads)
               (with-accessors ((node-force cl-mpm/mesh:node-force)
                                (node-ext-force cl-mpm/mesh::node-external-force)
                                (node-lock  cl-mpm/mesh:node-lock)
                                (node-vel  cl-mpm/mesh:node-velocity)
                                (node-mass  cl-mpm/mesh:node-mass)
                                (node-boundary cl-mpm/mesh::node-boundary-node)
                                (node-boundary-scalar cl-mpm/mesh::node-boundary-scalar)
                                (node-active  cl-mpm/mesh:node-active))
                   node
                 (declare (double-float volume svp))
                 ;;Lock node for multithreading
                 (when node-active
                   ;;Lock node for multithreading
                   (sb-thread:with-mutex (node-lock)
                     ;; (cl-mpm/fastmaths::fast-fmacc node-force
                     ;;                              force
                     ;;                              svp)
                     (cl-mpm/fastmaths::fast-fmacc node-ext-force
                                                  force
                                                  svp))))))
            (values normal-force)))))))

(defmethod cl-mpm/bc::apply-bc ((bc bc-penalty) node mesh dt)
  (with-accessors ((epsilon bc-penalty-epsilon)
                   (friction bc-penalty-friction)
                   (normal bc-penalty-normal)
                   (datum bc-penalty-datum)
                   (damping bc-penalty-damping)
                   (sub-bcs bc-penalty-structure-sub-bcs)
                   (debug-mutex bc-penalty-load-lock)
                   (debug-force bc-penalty-load)
                   (sim bc-penalty-sim))
      bc
    (setf (bc-penalty-contact-points bc) nil)
    (with-accessors ((mps cl-mpm:sim-mps)
                     (mesh cl-mpm:sim-mesh))
        sim
      (cl-mpm:iterate-over-mps
       mps
       (lambda (mp)
         (let ((in-contact nil)
               (closest-point (make-contact :penetration 0d0)))
           (cl-mpm/penalty::iterate-over-corners
            mesh
            mp
            (lambda (corner)
              (let* ((penetration-dist (penetration-distance-point corner datum normal)))
                (declare (double-float penetration-dist))
                (when (and
                       (>= penetration-dist 0d0)
                       (penalty-contact-valid bc corner))
                  (if in-contact
                      (cond
                        ((and (< (abs (- penetration-dist (contact-penetration closest-point))) 1d-3))
                         (setf (contact-point closest-point)
                               (cl-mpm/fastmaths:fast-.+
                                corner
                                (contact-point closest-point)))
                         (cl-mpm/fastmaths:fast-scale! (contact-point closest-point) 0.5d0))
                        ((< (abs penetration-dist) (abs (contact-penetration closest-point)))
                         (setf closest-point (make-contact
                                              :point corner
                                              :penetration penetration-dist))))
                      (progn
                        (setf in-contact t)
                        (setf closest-point (make-contact
                                             :point corner
                                             :penetration penetration-dist))))))))
           (when in-contact
             ;; (push (contact-point closest-point) (bc-penalty-contact-points bc))
             (apply-penalty-point mesh bc mp (contact-point closest-point) dt))))))))

(defmethod cl-mpm/bc::apply-bc ((bc bc-penalty-structure) node mesh dt)
  (with-accessors ((epsilon bc-penalty-epsilon)
                   (friction bc-penalty-friction)
                   (damping bc-penalty-damping)
                   (sub-bcs bc-penalty-structure-sub-bcs)
                   (debug-mutex bc-penalty-load-lock)
                   (debug-force bc-penalty-load)
                   (sim bc-penalty-sim))
      bc
    ;; (setf (bc-penalty-structure-contact-points bc) nil)
    (with-accessors ((mps cl-mpm:sim-mps)
                     (mesh cl-mpm:sim-mesh))
        sim
      (cl-mpm:iterate-over-mps
       mps
       (lambda (mp)
         (cl-mpm/penalty::iterate-over-corners
          mesh
          mp
          (lambda (corner)
            (cl-mpm/mesh::clamp-point-to-bounds mesh corner)
            (let ((in-contact nil)
                  (closest-point (make-contact :penetration 0d0)))
              (loop for sub-bc in sub-bcs
                    do (with-accessors ((datum bc-penalty-datum)
                                        (normal bc-penalty-normal))
                           sub-bc
                         (let* ((penetration-dist (penetration-distance-point corner datum normal)))
                           (declare (double-float penetration-dist))
                           (when (and
                                  (>= penetration-dist 0d0)
                                  (penalty-contact-valid sub-bc corner))
                             (if in-contact
                                 (cond
                                   ((< (abs penetration-dist) (abs (contact-penetration closest-point)))
                                    (setf closest-point (make-contact
                                                         :point corner
                                                         :datum (bc-penalty-datum sub-bc)
                                                         :penetration penetration-dist
                                                         :normal (bc-penalty-normal sub-bc)
                                                         :sub-bc sub-bc))))
                                 (progn
                                   (setf in-contact t)
                                   (setf closest-point (make-contact
                                                        :point corner
                                                        :datum (bc-penalty-datum sub-bc)
                                                        :penetration penetration-dist
                                                        :normal (bc-penalty-normal sub-bc)
                                                        :sub-bc sub-bc))))))))
              (when in-contact
                (let ((load (apply-penalty-point mesh
                                                 (contact-sub-bc closest-point)
                                                 mp
                                                 (contact-point closest-point) dt)))
                  ;; (break
                  ;;  (contact-sub-bc closest-point))
                  ;; (when (< (varef (contact-point closest-point) 1) 0d0)
                  ;;   (break "Corner out of bounds"))
                  (sb-thread:with-mutex (debug-mutex)
                    ;; (push (contact-point closest-point) (bc-penalty-structure-contact-points bc))
                    (incf debug-force load)))))))
         (setf (cl-mpm/particle::mp-penalty-energy mp) (calculate-bc-energy-mp sim bc mp)))))))

(defun calculate-bc-energy-mp (sim bc mp)
  (let ((total-energy 0d0))
    (with-accessors ((epsilon bc-penalty-epsilon)
                     (friction bc-penalty-friction)
                     (damping bc-penalty-damping)
                     (sub-bcs bc-penalty-structure-sub-bcs)
                     (debug-mutex bc-penalty-load-lock)
                     (debug-force bc-penalty-load))
        bc
      (with-accessors ((mps cl-mpm:sim-mps)
                       (mesh cl-mpm:sim-mesh)
                       (dt cl-mpm:sim-dt)
                       )
          sim
        (cl-mpm/penalty::iterate-over-corners
         mesh
         mp
         (lambda (corner)
           (cl-mpm/mesh::clamp-point-to-bounds mesh corner)
           (let ((in-contact nil)
                 (closest-point (make-contact :penetration 0d0)))
             (loop for sub-bc in sub-bcs
                   do (with-accessors ((datum bc-penalty-datum)
                                       (normal bc-penalty-normal))
                          sub-bc
                        (let* ((penetration-dist (penetration-distance-point corner datum normal)))
                          (declare (double-float penetration-dist))
                          (when (and
                                 (>= penetration-dist 0d0)
                                 (penalty-contact-valid sub-bc corner))
                            (if in-contact
                                (cond
                                  ((< (abs penetration-dist) (abs (contact-penetration closest-point)))
                                   (setf closest-point (make-contact
                                                        :point corner
                                                        :datum (bc-penalty-datum sub-bc)
                                                        :penetration penetration-dist
                                                        :normal (bc-penalty-normal sub-bc)
                                                        :sub-bc sub-bc))))
                                (progn
                                  (setf in-contact t)
                                  (setf closest-point (make-contact
                                                       :point corner
                                                       :datum (bc-penalty-datum sub-bc)
                                                       :penetration penetration-dist
                                                       :normal (bc-penalty-normal sub-bc)
                                                       :sub-bc sub-bc))))))))
             (when in-contact
               (let* ((contact-bc (contact-sub-bc closest-point))
                      (pen (contact-penetration closest-point))
                      (epsilon (bc-penalty-epsilon contact-bc))
                      (nd (cl-mpm/mesh:mesh-nd mesh))
                      (volume (cl-mpm/particle:mp-volume mp))
                      (load (* pen epsilon (expt volume (/ (- nd 1) nd)))))
                 (incf total-energy (expt load 2)))))))))
    total-energy))

(defun calculate-bc-energy (sim bc)
  (let ((total-energy 0d0)
        (mutex (sb-thread:make-mutex)))
    (with-accessors ((epsilon bc-penalty-epsilon)
                     (friction bc-penalty-friction)
                     (damping bc-penalty-damping)
                     (sub-bcs bc-penalty-structure-sub-bcs)
                     (debug-mutex bc-penalty-load-lock)
                     (debug-force bc-penalty-load))
        bc
      (with-accessors ((mps cl-mpm:sim-mps)
                       (mesh cl-mpm:sim-mesh)
                       (dt cl-mpm:sim-dt)
                       )
          sim
        (cl-mpm:iterate-over-mps
         mps
         (lambda (mp)
           (incf total-energy (calculate-bc-energy-mp sim bc mp))))))
    total-energy))

(defun calculate-mp-energy-gradient-2d (sim bc mp)
  (with-accessors ((pos cl-mpm/particle:mp-position))
      mp
    (let ((dx 1d-7)
          (energy-grads (cl-mpm/utils:vector-zeros))
          (e-0 (calculate-bc-energy-mp sim bc mp))
          (initial-pos (cl-mpm/utils::vector-copy (cl-mpm/particle:mp-position mp))))
      (declare (double-float e-0 dx))
      (loop for d from 0 to 1
            do (progn
                 (incf (cl-mpm/utils:varef pos d) dx)
                 (setf (varef energy-grads d) (/ (- (the double-float (calculate-bc-energy-mp sim bc mp)) e-0) dx))
                 (cl-mpm/utils:vector-copy-into pos initial-pos)))
      energy-grads)))


(defun apply-penalty-nudge (sim bc)
  (let ((total-energy 0d0)
      (mutex (sb-thread:make-mutex)))
    (with-accessors ((epsilon bc-penalty-epsilon)
                     (friction bc-penalty-friction)
                     (damping bc-penalty-damping)
                     (sub-bcs bc-penalty-structure-sub-bcs)
                     (debug-mutex bc-penalty-load-lock)
                     (debug-force bc-penalty-load))
        bc
      (with-accessors ((mps cl-mpm:sim-mps)
                       (mesh cl-mpm:sim-mesh)
                       (dt cl-mpm:sim-dt)
                       )
          sim
        (cl-mpm:iterate-over-mps
         mps
         (lambda (mp)
           (with-accessors ((pos cl-mpm/particle:mp-position)
                            (energy-old cl-mpm/particle::mp-penalty-energy))
               mp
             (when (> energy-old 0d0)
               (let ((normal (calculate-mp-energy-gradient-2d sim bc mp))
                     (energy-new energy-old)
                     (dstep (/ (- (cl-mpm/mesh:mesh-resolution mesh)) 1000)))
                 (loop for i from 0 to 1
                       while (and (> energy-new (* 0.1d0 energy-old))
                                  (> (cl-mpm/fastmaths:mag normal) 1d-15)
                                  )
                       do (progn
                            (setf normal (calculate-mp-energy-gradient-2d sim bc mp))
                            (setf normal (cl-mpm/fastmaths:fast-scale! (cl-mpm/fastmaths:norm normal) dstep))
                            (setf energy-new (calculate-bc-energy-mp sim bc mp))
                            (cl-mpm/fastmaths:fast-.+ pos normal pos))))))))))))
