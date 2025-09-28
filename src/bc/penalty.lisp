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
  (:import-from
   :cl-mpm/fastmaths  fast-.+ fast-.- fast-.*)
  (:export
   #:make-bc-penalty
   #:make-bc-penalty-distance-point
   #:save-vtk-penalties
   #:bc-penalty-friction))
;; (declaim (optimize (debug 0) (safety 0) (speed 3)))
;; (declaim (optimize (debug 3) (safety 3) (speed 0)))
(declaim #.cl-mpm/settings:*optimise-setting*)
(in-package :cl-mpm/penalty)

(defclass bc-penalty-structure (bc-penalty)
  ((sim
    :accessor bc-penalty-sim
    :initarg :sim)
   (sub-bcs
    :accessor bc-penalty-structure-sub-bcs
    :initarg :sub-bcs
    :initform (list)))
  (:documentation "A single multi-surface structure that should resolve contact through a closest-point algorithm"))

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
   (mp-stiffness
    :initform nil
    :accessor bc-mp-stiffness)
   (load
    :accessor bc-penalty-load
    :initform 0d0)
   (velocity
    :accessor bc-penalty-velocity
    :initform (cl-mpm/utils:vector-zeros))
   (load-lock
    :accessor bc-penalty-load-lock
    :initform (sb-thread:make-mutex)))
  (:documentation "A nonconforming neumann bc"))
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
                   (cl-mpm/fastmaths:fast-scale-vector
                    (cl-mpm/particle::mp-domain-size mp) 0.5d0)
                   (cl-mpm/fastmaths::norm
                    (cl-mpm/fastmaths:fast-.* normal normal))))
         (dist (- datum (- ypos yheight))))
    (the double-float dist)))

(defun penetration-point (mp pen datum normal)
  "Get linear penetration distance"
  (let* ((pos (cl-mpm/particle:mp-position mp))
         (domain (cl-mpm/particle::mp-domain-size mp)))
    (cl-mpm/fastmaths:fast-.- pos
               (fast-.* normal (cl-mpm/fastmaths:fast-scale-vector domain 0.5d0))
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
;; (defun apply-force-mps (mesh mps dt normal datum epsilon friction &optional (func-clip (lambda (mp) t))
;;                                                                     &key (damping 0.1d0))
;;   "Update force on nodes, with virtual stress field from mps"
;;   ;;If we lose contact we need to zero out our friction force
;;   (with-accessors ((nd cl-mpm/mesh::mesh-nd))
;;       mesh
;;     (cl-mpm:iterate-over-mps
;;      mps
;;      (lambda (mp)
;;        (let* ((penetration-dist (penetration-distance mp datum normal)))
;;          (declare (double-float penetration-dist))
;;          (when (> penetration-dist 0d0)
;;              (progn
;;                ;;Contact
;;                (with-accessors ((volume cl-mpm/particle:mp-volume)
;;                                 (pressure cl-mpm/particle::mp-pressure)
;;                                 (mp-vel cl-mpm/particle::mp-velocity)
;;                                 (mp-disp-inc cl-mpm/particle::mp-displacement-increment)
;;                                 (mp-mass cl-mpm/particle::mp-mass)
;;                                 (mp-contact cl-mpm/particle::mp-penalty-contact)
;;                                 (mp-friction cl-mpm/particle::mp-penalty-frictional-force)
;;                                 (mp-normal-force cl-mpm/particle::mp-penalty-normal-force))
;;                    mp
;;                  (let* ((pen-point (penetration-point mp penetration-dist datum normal))
;;                         (normal-force (* (expt (* epsilon penetration-dist) 1d0)
;;                                          (expt volume (/ (- nd 1) nd))
;;                                          )))
;;                    (when (funcall func-clip
;;                                   (trial-corner mp normal)
;;                                   )
;;                      (sb-thread:with-mutex (*debug-mutex*)
;;                        (incf *debug-force* (* normal-force 1d0))
;;                        (incf *debug-force-count* 1))
;;                      (setf mp-contact t)
;;                      ;;Iterate over neighbour nodes
;;                      (let* ((force (cl-mpm/utils:vector-zeros))
;;                             (rel-vel (cl-mpm/fastmaths:dot normal mp-vel))
;;                             (rel-disp (cl-mpm/fastmaths:dot normal mp-disp-inc))
;;                             (tang-disp (cl-mpm/fastmaths:fast-.- mp-disp-inc (cl-mpm/fastmaths:fast-scale-vector normal rel-disp)))
;;                             (tang-disp-norm-squared (cl-mpm/fastmaths::mag-squared tang-disp))
;;                             ;; (normal-damping (* damping (sqrt (/ epsilon (/ mp-mass volume)))))
;;                             (normal-damping (* (/ pi 2) damping (sqrt (* epsilon mp-mass))))
;;                             ;; (normal-damping (* (/ pi 2) damping (sqrt epsilon)))
;;                             (damping-force (* normal-damping rel-vel))
;;                             (force-friction mp-friction)
;;                             (stick-friction (* friction (abs normal-force))))
;;                        ;; update trial frictional force
;;                        (when (> friction 0d0)
;;                          (pprint mp-disp-inc)
;;                          (break)
;;                          (when (> tang-disp-norm-squared 0d0)
;;                            ;; We have sliding behaviour
;;                            (let* (;(tang-vel-norm (sqrt tang-vel-norm-squared))
;;                                   ;; (tang-normal (cl-mpm/fastmaths:norm tang-vel))
;;                                   ;;Trial friction
;;                                         ;(trial-friction-force (* (/ epsilon 2d0) tang-vel-norm dt))
;;                                   )
;;                              (cl-mpm/fastmaths::fast-fmacc
;;                               force-friction
;;                               tang-disp
;;                               (* -1d0 (/ epsilon 2d0)))))
;;                          (when (> (cl-mpm/fastmaths::mag-squared force-friction) 0d0)
;;                            (if (> (cl-mpm/fastmaths::mag force-friction) stick-friction)
;;                                (progn
;;                                  ;; (cl-mpm/fastmaths::fast-scale! force-friction
;;                                  ;;                              (/ (cl-mpm/fastmaths::mag force-friction)
;;                                  ;;                                 stick-friction))
;;                                  (setf force-friction
;;                                        (cl-mpm/fastmaths:fast-scale-vector
;;                                         (cl-mpm/fastmaths:norm force-friction)
;;                                         stick-friction))
;;                                  (setf (cl-mpm/particle::mp-penalty-friction-stick mp) t)
;;                                  )
;;                                (progn
;;                                  (setf (cl-mpm/particle::mp-penalty-friction-stick mp) nil)
;;                                  ))
;;                            (cl-mpm/fastmaths::fast-.+ force force-friction force))
;;                          (setf mp-friction force-friction))
;;                        (setf mp-normal-force (- normal-force damping-force))
;;                        (cl-mpm/fastmaths::fast-fmacc force
;;                                                     normal
;;                                                     (- normal-force damping-force))
;;                        (cl-mpm::iterate-over-neighbours-point-linear-3d
;;                         mesh
;;                         pen-point
;;                         (lambda (mesh node svp grads)
;;                           (with-accessors ((node-ext-force cl-mpm/mesh::node-external-force)
;;                                            (node-lock  cl-mpm/mesh:node-lock)
;;                                            (node-vel  cl-mpm/mesh:node-velocity)
;;                                            (node-mass  cl-mpm/mesh:node-mass)
;;                                            (node-boundary cl-mpm/mesh::node-boundary-node)
;;                                            (node-boundary-scalar cl-mpm/mesh::node-boundary-scalar)
;;                                            (node-active  cl-mpm/mesh:node-active))
;;                               node
;;                             (declare (double-float volume svp))
;;                             ;;Lock node for multithreading
;;                             (when node-active
;;                               ;;Lock node for multithreading
;;                               (sb-thread:with-mutex (node-lock)
;;                                 (cl-mpm/fastmaths::fast-fmacc node-ext-force
;;                                                              force
;;                                                              svp)))))))))))))))))

(defun collect-contact-points (mesh mps normal datum)
  (loop for mp across mps
        when (> (penetration-distance mp datum normal) 0d0)
          collect
          (let* ((penetration-dist (penetration-distance mp datum normal)))
            (declare (double-float penetration-dist))
            (when (> penetration-dist 0d0)
              (with-accessors ((volume cl-mpm/particle:mp-volume)
                               (pressure cl-mpm/particle::mp-pressure)
                               (mp-vel cl-mpm/particle::mp-velocity))
             mp
           (penetration-point mp penetration-dist datum normal))))))
(defun collect-contact-points-bc (mesh mps bc)
  (with-accessors ((normal bc-penalty-normal)
                   (datum bc-penalty-datum)
                   )
      bc
    (collect-contact-points mesh mps normal datum)))



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
                   :radius radius
                   :damping damping)))

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
                   :radius radius
                   :damping damping
                   )))

(defun make-bc-line-segments (sim positions epsilon friction damping)
  (loop for pa in (butlast positions)
        for pb in (rest positions)
        collect (make-bc-penalty-line-segment sim pa pb epsilon friction damping)))


(defun chaikin-smooth (point-list)
  (let ((new-list (list (first point-list))))
    (loop for p in (butlast point-list)
          for j in (rest point-list)
          do (let ((dir (cl-mpm/fastmaths:fast-.- j p)))
               (push (cl-mpm/fastmaths:fast-.+ p (cl-mpm/fastmaths:fast-scale-vector dir 0.25d0)) new-list)
               (push (cl-mpm/fastmaths:fast-.+ p (cl-mpm/fastmaths:fast-scale-vector dir 0.75d0)) new-list)
               ))
    (push (first (last point-list)) new-list)
    (reverse new-list)))

;; (pprint
;;  (cl-mpm/penalty::chaikin-smooth (list (cl-mpm/utils:vector-from-list (list 0d0 0d0 0d0))
;;                                        (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0))
;;                                        (cl-mpm/utils:vector-from-list (list 1d0 1d0 0d0))
;;                                        )))

(defun make-bc-penalty-smooth-corner (sim p-a corner p-b smooth-steps epsilon friction damping)
  (let ((points (list p-a corner p-b)))
    (dotimes (i smooth-steps)
      (setf points (chaikin-smooth points)))
    (let ((penalties (make-bc-line-segments sim points epsilon friction damping)))
      (cl-mpm/penalty::make-bc-penalty-structure
       sim
       epsilon
       friction
       damping
       penalties))))


(defgeneric bc-set-center (bc new-center))
(defmethod bc-set-center ((bc bc-penalty-distance) new-center)
  (with-accessors ((normal bc-penalty-normal)
                   (datum bc-penalty-datum)
                   (center-point bc-penalty-distance-center-point))
      bc
    (let ((new-datum (- (penetration-distance-point new-center 0d0 normal))))
      (setf center-point new-center
            datum new-datum))))

(defgeneric bc-increment-center (bc delta-center))

(defmethod bc-increment-center ((bc bc-penalty) delta-center)
  (with-accessors ((normal bc-penalty-normal)
                   (datum bc-penalty-datum))
      bc
    (let* ((center-point (cl-mpm/fastmaths:fast-scale-vector normal datum))
           (new-center (cl-mpm/fastmaths:fast-.+ center-point delta-center))
           (new-datum (- (penetration-distance-point new-center 0d0 normal))))
      (format t "Datum ~E - ~E~%" datum new-datum)
      (setf datum new-datum))))

(defmethod bc-increment-center ((bc bc-penalty-distance) delta-center)
  (with-accessors ((normal bc-penalty-normal)
                   (datum bc-penalty-datum)
                   (center-point bc-penalty-distance-center-point))
      bc
    (let* ((new-center (cl-mpm/fastmaths:fast-.+ center-point delta-center))
           (new-datum (- (penetration-distance-point new-center 0d0 normal))))
      (setf center-point new-center
            datum new-datum))))

(defmethod bc-increment-center ((bc bc-penalty-structure) delta-center)
  (loop for sub-bc in (bc-penalty-structure-sub-bcs bc)
        do (bc-increment-center sub-bc delta-center)))


(defmethod (setf bc-penalty-velocity) :after (value (bc bc-penalty-structure))
  (loop for sub-bc in (bc-penalty-structure-sub-bcs bc)
        do (setf (bc-penalty-velocity sub-bc) value)))

(defmethod (setf bc-penalty-friction) :after (value (bc bc-penalty-structure))
  (loop for sub-bc in (bc-penalty-structure-sub-bcs bc)
        do (setf (bc-penalty-friction sub-bc) value)))

(defgeneric early-sweep-intersection (bc mp))
(defmethod early-sweep-intersection ((bc cl-mpm/bc::bc) mp)
  t)

(defmethod early-sweep-intersection ((bc bc-penalty) mp)
  (with-accessors ((datum bc-penalty-datum)
                   (normal bc-penalty-normal))
      bc
      (let ((penetration-dist (penetration-distance-point (cl-mpm/particle::mp-position mp) datum normal))
            (domain (cl-mpm/particle::mp-domain-size mp)))
        (< penetration-dist
           (sqrt
            (+
             (expt (varef domain 0) 2)
             (expt (varef domain 1) 2)
             (expt (varef domain 2) 2)))))))

(defmethod early-sweep-intersection ((bc bc-penalty-distance) mp)
  (with-accessors ((datum bc-penalty-datum)
                   (normal bc-penalty-normal))
      bc
    (let ((penetration-dist (penetration-distance-point (cl-mpm/particle::mp-position mp) datum normal))
          (domain (cl-mpm/particle::mp-domain-size mp)))
      (< penetration-dist
          (sqrt
           (+
            (expt (varef domain 0) 2)
            (expt (varef domain 1) 2)
            (expt (varef domain 2) 2)))))))

(defmethod early-sweep-intersection ((bc bc-penalty-structure) mp)
  (let ((contact nil))
    (loop for sub-bc in (bc-penalty-structure-sub-bcs bc)
          while (not contact)
          do (with-accessors ((datum bc-penalty-datum)
                              (normal bc-penalty-normal))
                 sub-bc
               (when (early-sweep-intersection sub-bc mp)
                 (setf contact t))))
    contact))



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
         (yheight (cl-mpm/fastmaths::dot (cl-mpm/fastmaths:fast-scale-vector (cl-mpm/particle::mp-domain-size mp) 0.5d0)
                                        (cl-mpm/fastmaths::norm (fast-.* normal normal))))
         (dist (- datum (- ypos yheight))))
    (the double-float dist)))

(defun disp-penetration-point (mp pen datum normal)
  "Get linear penetration distance"
  (let* ((pos (cl-mpm/particle:mp-position mp))
         (domain (cl-mpm/particle::mp-domain-size mp)))
    (fast-.- pos
               (fast-.* normal (cl-mpm/fastmaths:fast-scale-vector domain 0.5d0)))))


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
                                   (expt volume (/ (- nd 1) nd)))))
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



(defun make-bc-penalty-structure (sim epsilon friction damping sub-bcs)
  (make-instance 'bc-penalty-structure
                 :index nil
                 :sim sim
                 :epsilon epsilon
                 :friction friction
                 :damping damping
                 :sub-bcs sub-bcs))

(defstruct (contact
            (:constructor make-contact
                (point datum penetration normal sub-bc)))
  (point nil :type (or magicl:matrix/double-float null))
  (trial-point nil :type (or magicl:matrix/double-float null))
  (disp nil :type (or magicl:matrix/double-float null))
  (datum 0d0 :type double-float)
  (penetration 0d0 :type double-float)
  (normal nil :type (or magicl:matrix/double-float null))
  (sub-bc nil :type (or cl-mpm/bc::bc null)))
(declaim (inline make-contact))


(defun compute-mp-stiffness (mp)
  (sqrt
   (/
    (/ (the double-float (cl-mpm/particle:mp-mass mp))
       (the double-float (cl-mpm/particle:mp-volume mp)))
    1d0
    ;; (cl-mpm/particle::mp-e mp)
    ;; (cl-mpm/particle::mp-p-modulus mp)
    )))

(defun apply-penalty-point (mesh bc mp
                            point
                            trial-point
                            disp-inc
                            dt)
  (with-accessors ((nd cl-mpm/mesh:mesh-nd))
      mesh
    (with-accessors ((datum bc-penalty-datum)
                     (damping bc-penalty-damping)
                     (epsilon bc-penalty-epsilon)
                     (friction bc-penalty-friction)
                     (normal bc-penalty-normal)
                     (pen-vel bc-penalty-velocity)
                     (mp-stiffness bc-mp-stiffness)
                     (debug-mutex bc-penalty-load-lock)
                     (debug-load bc-penalty-load))
        bc
      (with-accessors ((volume cl-mpm/particle:mp-volume)
                       (pressure cl-mpm/particle::mp-pressure)
                       (mp-vel cl-mpm/particle::mp-velocity)
                       ;; (mp-disp-inc cl-mpm/particle::mp-displacement-increment)
                       (mp-mass cl-mpm/particle::mp-mass)
                       (mp-contact cl-mpm/particle::mp-penalty-contact)
                       (mp-friction cl-mpm/particle::mp-penalty-frictional-force)
                       (mp-friction-n cl-mpm/particle::mp-penalty-frictional-force-prev)
                       (mp-normal-force cl-mpm/particle::mp-penalty-normal-force))
          mp
        (let* ((penetration (penetration-distance-point point datum normal))
               (contact-area (expt volume (/ (- nd 1) nd)))
               (normal-force (* ;(signum penetration) (expt (* (abs penetration)) 1.5d0)
                              (expt penetration 1)
                              epsilon
                              contact-area)))
          (when (> penetration 0d0)
              (sb-thread:with-mutex (debug-mutex)
                (incf debug-load normal-force))
            (setf mp-contact t)
            ;; (break)
            (let ((new-stiff (compute-mp-stiffness mp)))
              (setf mp-stiffness (min new-stiff (if mp-stiffness mp-stiffness new-stiff))))

            (setf (cl-mpm/particle::mp-penalty-stiffness mp) (* 2d0 epsilon contact-area (* 2d0 (+ 1d0 friction))))
            (setf (cl-mpm/particle::mp-penalty-contact-point mp) trial-point)

            (let* (
                   ;; (mp-disp-inc (cl-mpm/fastmaths:fast-scale-vector mp-vel dt))
                   (mp-disp-inc disp-inc)
                   (force (cl-mpm/utils:vector-zeros))
                   (resultant-vel (cl-mpm/fastmaths:fast-.- mp-vel pen-vel))
                   (rel-vel (cl-mpm/fastmaths:dot normal resultant-vel))
                   (rel-disp (cl-mpm/fastmaths:dot normal mp-disp-inc))
                   (tang-disp (cl-mpm/fastmaths:fast-.-
                               mp-disp-inc
                               (cl-mpm/fastmaths:fast-scale-vector normal rel-disp)))
                   (tang-disp-norm-squared (cl-mpm/fastmaths::mag-squared tang-disp))
                   ;; (tang-vel (cl-mpm/fastmaths:fast-.- resultant-vel
                   ;;                                     (cl-mpm/fastmaths:fast-scale-vector normal rel-vel)))
                   ;; (tang-vel-norm-squared (cl-mpm/fastmaths::mag-squared tang-vel))

                   (normal-damping (* 2d0 damping (sqrt (* epsilon mp-mass)) contact-area))
                   (damping-force (* normal-damping rel-vel))
                   (force-friction mp-friction)
                   (stick-friction (* friction (abs normal-force))))
              ;; update trial frictional force

              (when (> friction 0d0)
                ;; (pprint tang-disp-norm-squared)
                (when (> tang-disp-norm-squared 0d0)
                  ;; (pprint tang-disp)
                  (cl-mpm/utils:vector-copy-into mp-friction-n force-friction)
                  (cl-mpm/fastmaths::fast-fmacc
                   force-friction
                   tang-disp
                   (* -1d0 (/ epsilon 2d0)))

                  (setf (cl-mpm/particle::mp-penalty-stiffness mp) (+ (* 2d0 epsilon contact-area) (* contact-area (cl-mpm/fastmaths::mag force-friction))))
                  )

                (when (> (cl-mpm/fastmaths::mag-squared force-friction) 0d0)
                  (if (> (cl-mpm/fastmaths::mag force-friction) stick-friction)
                      (progn
                        (setf force-friction
                              (cl-mpm/fastmaths:fast-scale-vector
                               (cl-mpm/fastmaths:norm force-friction)
                               stick-friction))
                        (setf (cl-mpm/particle::mp-penalty-friction-stick mp) t))
                      (progn
                        (setf (cl-mpm/particle::mp-penalty-friction-stick mp) nil)))
                  (cl-mpm/fastmaths::fast-.+ force force-friction force))
                (setf mp-friction force-friction))
              (setf mp-normal-force (- normal-force damping-force))
              ;; (cl-mpm/fastmaths::fast-fmacc force
              ;;                              normal
              ;;                              (- normal-force damping-force))
              (cl-mpm/fastmaths::fast-fmacc force
                                            normal
                                            normal-force)
              (cl-mpm::iterate-over-neighbours-point-linear
               mesh
               trial-point
               (lambda (mesh node svp grads)
                 (with-accessors ((node-ext-force cl-mpm/mesh::node-external-force)
                                  (node-damp-force cl-mpm/mesh::node-damping-force)
                                  (node-lock  cl-mpm/mesh:node-lock)
                                  (node-mass  cl-mpm/mesh:node-mass)
                                  (node-volume  cl-mpm/mesh::node-volume)
                                  (node-volume-t  cl-mpm/mesh::node-volume-true)
                                  (node-active  cl-mpm/mesh:node-active))
                     node
                   (declare (double-float volume svp))
                   ;;Lock node for multithreading
                   (when node-active
                     ;; (break)
                     ;;Lock node for multithreading
                     ;; (let ((svp (* svp (/ node-volume node-volume-t)))))
                     (sb-thread:with-mutex (node-lock)
                       (cl-mpm/fastmaths::fast-fmacc node-ext-force
                                                     force
                                                     svp)
                       (cl-mpm/fastmaths::fast-fmacc node-damp-force
                                                     normal
                                                     (* -1d0 svp damping-force)))))))
              (values normal-force))))))))

(defgeneric resolve-load-direction (bc direction))

(defmethod resolve-load-direction ((bc bc-penalty) direction)
  (* (cl-mpm/fastmaths:dot
      direction
      (bc-penalty-normal bc))
     (bc-penalty-load bc)))

(defmethod resolve-load-direction ((bc bc-penalty-structure) direction)
  (loop for sub-bc in (bc-penalty-structure-sub-bcs bc)
        sum (resolve-load-direction sub-bc direction)))

(defgeneric resolve-load (bc))

(defmethod resolve-load ((bc bc-penalty))
  (bc-penalty-load bc))

(defmethod resolve-load ((bc bc-penalty-structure))
  (loop for sub-bc in (bc-penalty-structure-sub-bcs bc)
        sum (resolve-load sub-bc)))

(defgeneric reset-load (bc))
(defmethod reset-load ((bc bc-penalty))
  (setf (bc-penalty-load bc) 0d0))
(defmethod reset-load ((bc bc-penalty-structure))
  (setf (bc-penalty-load bc) 0d0)
  (loop for sub-bc in (bc-penalty-structure-sub-bcs bc)
        sum (reset-load sub-bc)))

(defun compute-corner-displacement (mesh mp corner)
  ;; (cl-mpm/fastmaths:fast-.+ corner (cl-mpm/particle::mp-displacement-increment mp))
  (let ((corner-disp (cl-mpm/utils::vector-zeros)))
    (cl-mpm::iterate-over-neighbours-point-linear
     mesh
     corner
     (lambda (mesh node svp grads)
       (with-accessors ((node-disp cl-mpm/mesh::node-displacment)
                        (node-active  cl-mpm/mesh:node-active))
           node
         (when node-active
           (cl-mpm/fastmaths:fast-fmacc corner-disp node-disp svp)))))
    corner-disp))

(defmethod cl-mpm/bc::apply-bc ((bc bc-penalty) node mesh dt)
  (with-accessors ((epsilon bc-penalty-epsilon)
                   (friction bc-penalty-friction)
                   (normal bc-penalty-normal)
                   (datum bc-penalty-datum)
                   (damping bc-penalty-damping)
                   (sub-bcs bc-penalty-structure-sub-bcs)
                   (debug-mutex bc-penalty-load-lock)
                   (debug-force bc-penalty-load)
                   (mp-stiffness bc-mp-stiffness)
                   (sim bc-penalty-sim))
      bc

    (reset-load bc)
    (setf (bc-penalty-contact-points bc) nil)
    (setf mp-stiffness nil)
    (with-accessors ((mps cl-mpm:sim-mps)
                     (mesh cl-mpm:sim-mesh))
        sim
      (cl-mpm:iterate-over-mps
       mps
       (lambda (mp)
         (when t;(early-sweep-intersection bc mp)
           (let ((in-contact nil)
                 (closest-point (make-contact
                                 nil
                                 0d0
                                 0d0
                                 nil
                                 nil)))
             (;;cl-mpm::iterate-over-corners
              cl-mpm::iterate-over-midpoints
              mesh
              mp
              (lambda (corner-trial)
                ;; (pprint corner-trial)
                (let* ((disp (compute-corner-displacement mesh mp corner-trial))
                       (corner (cl-mpm/fastmaths:fast-.+ corner-trial disp)))
                  ;; (break)
                  (let* ((penetration-dist (penetration-distance-point corner datum normal)))
                    (declare (double-float penetration-dist))
                    (when (and
                           (> penetration-dist 0d0)
                           (penalty-contact-valid bc corner))
                      (if in-contact
                          (cond
                            ;; ((and (< (abs (- penetration-dist (contact-penetration closest-point))) 1d-3))
                            ;;  (setf (contact-point closest-point)
                            ;;        (cl-mpm/fastmaths:fast-.+
                            ;;         corner
                            ;;         (contact-point closest-point)))
                            ;;  (cl-mpm/fastmaths:fast-scale! (contact-point closest-point) 0.5d0))
                            ((< (abs penetration-dist) (abs (contact-penetration closest-point)))
                             (setf
                              (contact-point closest-point) corner
                              (contact-penetration closest-point) penetration-dist
                              )))
                          (progn
                            (setf in-contact t)
                            (setf
                             (contact-point closest-point) corner
                             (contact-trial-point closest-point) corner-trial
                             (contact-disp closest-point) disp
                             (contact-penetration closest-point) penetration-dist))))))))
             (when in-contact
               (push (contact-point closest-point) (bc-penalty-contact-points bc))
               ;; (break)
               (apply-penalty-point mesh bc mp
                                    (contact-point closest-point)
                                    (contact-trial-point closest-point)
                                    (contact-disp closest-point)
                                    dt)))))))))


(defgeneric resolve-closest-contact (bc corner))
(defmethod resolve-closest-contact ((bc cl-mpm/penalty::bc-penalty) corner)
  (with-accessors ((datum bc-penalty-datum)
                   (normal bc-penalty-normal))
      bc
    (let* ((penetration-dist (penetration-distance-point corner datum normal)))
      (declare (double-float penetration-dist))
      (if (and
           (>= penetration-dist 0d0)
           (penalty-contact-valid bc corner))
          (progn
            (values t
                    penetration-dist
                    (make-contact
                     corner
                     (bc-penalty-datum bc)
                     penetration-dist
                     (bc-penalty-normal bc)
                     bc)
                    ))
          (values nil penetration-dist nil)))))

(defmethod resolve-closest-contact ((bc cl-mpm/penalty::bc-penalty-structure) corner)
  (with-accessors ((sub-bcs bc-penalty-structure-sub-bcs))
      bc
    (let ((in-contact nil)
          (closest-point (make-contact
                          nil
                          0d0
                          0d0
                          nil
                          nil)))
      (loop for sub-bc in sub-bcs
            do (with-accessors ((datum bc-penalty-datum)
                                (normal bc-penalty-normal))
                   sub-bc
                 (multiple-value-bind (inc pen point) (resolve-closest-contact sub-bc corner)
                   (when inc
                     (if in-contact
                         (cond
                           ((< (abs pen) (abs (contact-penetration closest-point)))
                            (setf closest-point point)))
                         (progn
                           (setf in-contact t)
                           (setf closest-point point)))))))
      (values in-contact (contact-penetration closest-point) closest-point))))

(defmethod cl-mpm/bc::apply-bc ((bc bc-penalty-structure) node mesh dt)
  (with-accessors ((epsilon bc-penalty-epsilon)
                   (friction bc-penalty-friction)
                   (damping bc-penalty-damping)
                   (sub-bcs bc-penalty-structure-sub-bcs)
                   (debug-mutex bc-penalty-load-lock)
                   (debug-force bc-penalty-load)
                   (sim bc-penalty-sim))
      bc
    (reset-load bc)
    (with-accessors ((mps cl-mpm:sim-mps)
                     (mesh cl-mpm:sim-mesh))
        sim
      (cl-mpm:iterate-over-mps
       mps
       (lambda (mp)
         (let (;; (disp-inc (cl-mpm/particle::mp-displacement-increment mp))
               )
           (when t;(early-sweep-intersection bc mp)
             (cl-mpm::iterate-over-midpoints
              ;; cl-mpm::iterate-over-corners
              mesh
              mp
              (lambda (corner-trial)
                (let* ((disp (compute-corner-displacement mesh mp corner-trial))
                       (corner (cl-mpm/fastmaths:fast-.+ corner-trial disp)))
                  (cl-mpm/mesh::clamp-point-to-bounds mesh corner)
                  (multiple-value-bind (in-contact pen closest-point) (resolve-closest-contact bc corner)
                    ;; (cl-mpm/fastmaths:fast-.- (contact-point closest-point) disp-inc (contact-point closest-point))
                    (when in-contact
                      (setf (contact-trial-point closest-point) corner-trial
                             (contact-disp closest-point) disp)
                      (let ((load (apply-penalty-point mesh
                                                       (contact-sub-bc closest-point)
                                                       mp
                                                       (contact-point closest-point)
                                                       (contact-trial-point closest-point)
                                                       (contact-disp closest-point)
                                                       dt)))
                        (push (contact-point closest-point) (bc-penalty-contact-points bc))
                        (sb-thread:with-mutex (debug-mutex)
                          (incf debug-force load)))))))))))))))

;; (defun calculate-bc-energy-mp (sim bc mp)
;;   (let ((total-energy 0d0))
;;     (with-accessors ((epsilon bc-penalty-epsilon)
;;                      (friction bc-penalty-friction)
;;                      (damping bc-penalty-damping)
;;                      (sub-bcs bc-penalty-structure-sub-bcs)
;;                      (debug-mutex bc-penalty-load-lock)
;;                      (debug-force bc-penalty-load))
;;         bc
;;       (with-accessors ((mps cl-mpm:sim-mps)
;;                        (mesh cl-mpm:sim-mesh)
;;                        (dt cl-mpm:sim-dt)
;;                        )
;;           sim
;;         (cl-mpm::iterate-over-corners
;;          mesh
;;          mp
;;          (lambda (corner)
;;            (cl-mpm/mesh::clamp-point-to-bounds mesh corner)
;;            (let ((in-contact nil)
;;                  (closest-point (make-contact :penetration 0d0)))
;;              (loop for sub-bc in sub-bcs
;;                    do (with-accessors ((datum bc-penalty-datum)
;;                                        (normal bc-penalty-normal))
;;                           sub-bc
;;                         (let* ((penetration-dist (penetration-distance-point corner datum normal)))
;;                           (declare (double-float penetration-dist))
;;                           (when (and
;;                                  (>= penetration-dist 0d0)
;;                                  (penalty-contact-valid sub-bc corner))
;;                             (if in-contact
;;                                 (cond
;;                                   ((< (abs penetration-dist) (abs (contact-penetration closest-point)))
;;                                    (setf closest-point (make-contact
;;                                                         :point corner
;;                                                         :datum (bc-penalty-datum sub-bc)
;;                                                         :penetration penetration-dist
;;                                                         :normal (bc-penalty-normal sub-bc)
;;                                                         :sub-bc sub-bc))))
;;                                 (progn
;;                                   (setf in-contact t)
;;                                   (setf closest-point (make-contact
;;                                                        :point corner
;;                                                        :datum (bc-penalty-datum sub-bc)
;;                                                        :penetration penetration-dist
;;                                                        :normal (bc-penalty-normal sub-bc)
;;                                                        :sub-bc sub-bc))))))))
;;              (when in-contact
;;                (let* ((contact-bc (contact-sub-bc closest-point))
;;                       (pen (contact-penetration closest-point))
;;                       (epsilon (bc-penalty-epsilon contact-bc))
;;                       (nd (cl-mpm/mesh:mesh-nd mesh))
;;                       (volume (cl-mpm/particle:mp-volume mp))
;;                       (load (* pen epsilon (expt volume (/ (- nd 1) nd)))))
;;                  (incf total-energy (expt load 2)))))))))
;;     total-energy))

;; (defun calculate-bc-energy (sim bc)
;;   (let ((total-energy 0d0)
;;         (mutex (sb-thread:make-mutex)))
;;     (with-accessors ((epsilon bc-penalty-epsilon)
;;                      (friction bc-penalty-friction)
;;                      (damping bc-penalty-damping)
;;                      (sub-bcs bc-penalty-structure-sub-bcs)
;;                      (debug-mutex bc-penalty-load-lock)
;;                      (debug-force bc-penalty-load))
;;         bc
;;       (with-accessors ((mps cl-mpm:sim-mps)
;;                        (mesh cl-mpm:sim-mesh)
;;                        (dt cl-mpm:sim-dt)
;;                        )
;;           sim
;;         (cl-mpm:iterate-over-mps
;;          mps
;;          (lambda (mp)
;;            (incf total-energy (calculate-bc-energy-mp sim bc mp))))))
;;     total-energy))

;; (defun calculate-mp-energy-gradient-2d (sim bc mp)
;;   (with-accessors ((pos cl-mpm/particle:mp-position))
;;       mp
;;     (let ((dx 1d-7)
;;           (energy-grads (cl-mpm/utils:vector-zeros))
;;           (e-0 (calculate-bc-energy-mp sim bc mp))
;;           (initial-pos (cl-mpm/utils::vector-copy (cl-mpm/particle:mp-position mp))))
;;       (declare (double-float e-0 dx))
;;       (loop for d from 0 to 1
;;             do (progn
;;                  (incf (cl-mpm/utils:varef pos d) dx)
;;                  (setf (varef energy-grads d) (/ (- (the double-float (calculate-bc-energy-mp sim bc mp)) e-0) dx))
;;                  (cl-mpm/utils:vector-copy-into pos initial-pos)))
;;       energy-grads)))


;; (defun apply-penalty-nudge (sim bc)
;;   (let ((total-energy 0d0)
;;       (mutex (sb-thread:make-mutex)))
;;     (with-accessors ((epsilon bc-penalty-epsilon)
;;                      (friction bc-penalty-friction)
;;                      (damping bc-penalty-damping)
;;                      (sub-bcs bc-penalty-structure-sub-bcs)
;;                      (debug-mutex bc-penalty-load-lock)
;;                      (debug-force bc-penalty-load))
;;         bc
;;       (with-accessors ((mps cl-mpm:sim-mps)
;;                        (mesh cl-mpm:sim-mesh)
;;                        (dt cl-mpm:sim-dt)
;;                        )
;;           sim
;;         (cl-mpm:iterate-over-mps
;;          mps
;;          (lambda (mp)
;;            (with-accessors ((pos cl-mpm/particle:mp-position)
;;                             (energy-old cl-mpm/particle::mp-penalty-energy))
;;                mp
;;              (when (> energy-old 0d0)
;;                (let ((normal (calculate-mp-energy-gradient-2d sim bc mp))
;;                      (energy-new energy-old)
;;                      (dstep (/ (- (cl-mpm/mesh:mesh-resolution mesh)) 1000)))
;;                  (loop for i from 0 to 1
;;                        while (and (> energy-new (* 0.1d0 energy-old))
;;                                   (> (cl-mpm/fastmaths:mag normal) 1d-15)
;;                                   )
;;                        do (progn
;;                             (setf normal (calculate-mp-energy-gradient-2d sim bc mp))
;;                             (setf normal (cl-mpm/fastmaths:fast-scale! (cl-mpm/fastmaths:norm normal) dstep))
;;                             (setf energy-new (calculate-bc-energy-mp sim bc mp))
;;                             (cl-mpm/fastmaths:fast-.+ pos normal pos))))))))))))

(defun save-vtk-penalties (filename sim)
  (with-accessors ((bcs cl-mpm:sim-bcs-force-list))
      sim
    (let ((bc-save-list (list)))
      (loop for bc-list in bcs
            do
               (loop for bc across bc-list
                     do
                        (typecase bc
                          (cl-mpm/penalty::bc-penalty-distance
                           (push bc bc-save-list)
                           )
                          (cl-mpm/penalty::bc-penalty-structure
                           (setf bc-save-list (append bc-save-list (get-all-bcs bc))))
                          (t ))))
      (let ((json-object (list)))
        (loop for bc in bc-save-list
              do
                 (with-accessors ((center cl-mpm/penalty::bc-penalty-distance-center-point )
                                  (normal cl-mpm/penalty::bc-penalty-normal)
                                  (radius cl-mpm/penalty::bc-penalty-distance-radius))
                     bc
                   (let* ((line-dir (cl-mpm/penalty::2d-orthog normal))
                          (p1 (cl-mpm/fastmaths:fast-.+ center (cl-mpm/fastmaths:fast-scale-vector line-dir radius)))
                          (p2 (cl-mpm/fastmaths:fast-.- center (cl-mpm/fastmaths:fast-scale-vector line-dir radius)))
                          (load (cl-mpm/penalty::bc-penalty-load bc))
                          (size (cl-mpm/penalty::bc-penalty-distance-radius bc))
                          )
                     (push (list (list (list (cl-mpm/utils:varef p1 0)
                                             (cl-mpm/utils:varef p1 1)
                                             (cl-mpm/utils:varef p1 2))
                                       (list (cl-mpm/utils:varef p2 0)
                                             (cl-mpm/utils:varef p2 1)
                                             (cl-mpm/utils:varef p2 2)))
                                 ;(/ load size)
                                 load
                                 (list (cl-mpm/utils:varef normal 0)
                                       (cl-mpm/utils:varef normal 1)
                                       (cl-mpm/utils:varef normal 2))
                                 ) json-object))))
        (with-open-file (fs filename :direction :output :if-exists :supersede)
          (format fs "# vtk DataFile Version 2.0~%")
          (format fs "Lisp generated vtk file, SJVS~%")
          (format fs "ASCII~%")
          (format fs "DATASET UNSTRUCTURED_GRID~%")

          (let* ((line-count (length json-object))
                 (point-count (* 2 line-count)))
            (format fs "POINTS ~d double~%" point-count)
            (loop for bc in json-object
                  do (loop for p in (first bc)
                           do (format fs "~E ~E ~E ~%"
                                      (coerce (first p) 'single-float)
                                      (coerce (second p) 'single-float)
                                      0e0)))

            (format fs "~%")
            (format fs "CELLS ~D ~D~%" line-count (* 3 line-count))
            (loop for bc in json-object
                  for i from 0
                  do (format fs "~D ~D ~D~%" 2 (* i 2) (+ 1 (* i 2))))

            (format fs "CELL_TYPES ~D~%" line-count)
            (loop for bc in json-object
                  do (format fs "3~%"))

            (format fs "CELL_DATA ~d~%" line-count)
            (format fs "SCALARS ~a FLOAT ~d~%" "LOAD" 1)
            (format fs "LOOKUP_TABLE default~%")
            (loop for bc in json-object
                  do (progn
                       (format fs "~E ~%"
                               (coerce (second bc) 'single-float))
                       ))
            (format fs "NORMALS ~a FLOAT~%" "NORMAL")
            (loop for bc in json-object
                  do (progn
                       (destructuring-bind (x y z) (third bc)
                         (format fs "~f ~f ~f~%"
                                 (coerce x 'single-float)
                                 (coerce y 'single-float)
                                 (coerce z 'single-float)
                                 ))))))))))

(defmethod cl-mpm/bc::calculate-min-dt-bc (sim (bc bc-penalty))
  (when (bc-mp-stiffness bc)
    (*
     (the double-float (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
     (the double-float (sqrt (the double-float (cl-mpm::sim-mass-scale sim))))
     (/ 1d0 (sqrt (bc-penalty-epsilon bc)))
     (the double-float (the double-float (bc-mp-stiffness bc))))))


(defgeneric get-contact-points (bc))
(defmethod get-contact-points ((bc bc-penalty))
  (bc-penalty-contact-points bc))

(defmethod get-contact-points ((bc bc-penalty-structure))
  (append 
   (bc-penalty-contact-points bc)
   (loop for sub-bc in (bc-penalty-structure-sub-bcs bc)
         append (bc-penalty-contact-points bc))))



(defclass bc-penalty-square (bc-penalty)
  ((center-point
    :accessor bc-penalty-square-center-point
    :initarg :center-point
    :initform (cl-mpm/utils:vector-zeros))
   (diag-a
    :accessor bc-penalty-square-diag-a
    :initarg :diag-a
    :initform (cl-mpm/utils:vector-zeros))
   (diag-b
    :accessor bc-penalty-square-diag-b
    :initarg :diag-b
    :initform (cl-mpm/utils:vector-zeros))
   ))

(defun make-bc-penalty-square (sim normal point diag-a diag-b epsilon friction damping)
  "Make a square penalty condition, where its is required that point-a.normal == point-b.normal"
  (let* ((normal (cl-mpm/fastmaths::norm normal))
         ;; (diag-a (cl-mpm/fastmaths:fast-.- point-a point))
         ;; (diag-b (cl-mpm/fastmaths:fast-.- point-b point))
         (datum (- (penetration-distance-point point 0d0 normal))))
    (make-instance 'bc-penalty-square
                   :index nil
                   :sim sim
                   :datum datum
                   :normal normal
                   :epsilon epsilon
                   :friction friction
                   :center-point point
                   :diag-a diag-a
                   :diag-b diag-b)))

(defmethod bc-set-center ((bc bc-penalty-square) new-center)
  (with-accessors ((normal bc-penalty-normal)
                   (datum bc-penalty-datum)
                   (center-point bc-penalty-square-center-point))
      bc
    (let ((new-datum (- (penetration-distance-point new-center 0d0 normal))))
      (setf center-point new-center
            datum new-datum))))
(defmethod bc-increment-center ((bc bc-penalty-square) delta-center)
  (with-accessors ((normal bc-penalty-normal)
                   (datum bc-penalty-datum)
                   (center-point bc-penalty-square-center-point))
      bc
    (let* ((new-center (cl-mpm/fastmaths:fast-.+ center-point delta-center))
           (new-datum (- (penetration-distance-point new-center 0d0 normal))))
      (setf center-point new-center
            datum new-datum))))

(defmethod penalty-contact-valid ((bc bc-penalty-square) point)
  (let* ((normal (bc-penalty-normal bc))
         (center-point (bc-penalty-square-center-point bc))
         (c-d (cl-mpm/fastmaths:dot center-point normal))
         (p-d (cl-mpm/fastmaths:dot point normal))
         (trial-center-point (cl-mpm/fastmaths:fast-.- center-point (cl-mpm/fastmaths:fast-scale-vector normal c-d)))
         (trial-plane-point (cl-mpm/fastmaths:fast-.- point (cl-mpm/fastmaths:fast-scale-vector normal p-d)))
         (a-b (bc-penalty-square-diag-a bc))
         (b-c (bc-penalty-square-diag-b bc))
         (a-trial (cl-mpm/fastmaths:fast-.- (cl-mpm/fastmaths:fast-.+ trial-center-point a-b) trial-plane-point ))
         (b-trial (cl-mpm/fastmaths:fast-.- (cl-mpm/fastmaths:fast-.+ trial-center-point b-c) trial-plane-point ))
         )
    ;; (break)
    (and
     (<= 0d0 (cl-mpm/fastmaths:dot a-b a-trial))
     (<= (cl-mpm/fastmaths:dot a-b a-trial) (cl-mpm/fastmaths:dot a-b a-b))
     (<= 0d0 (cl-mpm/fastmaths:dot b-c b-trial))
     (<= (cl-mpm/fastmaths:dot b-c b-trial) (cl-mpm/fastmaths:dot b-c b-c))
     )))

(defun iterate-over-penalty-point-normal ())

(defun assemble-penalty-stiffness-matrix (sim)
  "For dynamic relaxation, it is helpful to assemble the estimated maximum stiffness onto the mass matrix"
  (with-accessors ((mps sim-mps)
                   (mesh sim-mesh)
                   (dt-scale cl-mpm::sim-dt-scale))
      sim
    ;; (cl-mpm::iterate-over-bcs-force
    ;;  sim
    ;;  (lambda (bc)
    ;;    ()

    ;;    ))
    (cl-mpm::iterate-over-mps
       mps
       (lambda (mp)
         (when (cl-mpm/particle::mp-penalty-contact mp)
           (with-accessors ((mp-stiffness cl-mpm/particle::mp-penalty-stiffness))
               mp
             ;; (cl-mpm::iterate-over-midpoints
             ;;  mesh mp
             ;;  (lambda (point)
             ;;    (cl-mpm::iterate-over-neighbours-point-linear
             ;;     mesh
             ;;     point
             ;;     (lambda (mesh node svp grads)
             ;;       (with-accessors ((node-active cl-mpm/mesh:node-active)
             ;;                        (node-volume cl-mpm/mesh::node-volume)
             ;;                        (node-mass cl-mpm/mesh::node-mass)
             ;;                        (node-lock cl-mpm/mesh::node-lock))
             ;;           node
             ;;         (declare (double-float node-mass node-volume mp-stiffness))
             ;;         (when node-active
             ;;           (sb-thread:with-mutex (node-lock)
             ;;             (setf
             ;;              node-mass
             ;;              (+
             ;;               node-mass
             ;;               (* svp mp-stiffness))
             ;;              ;; (max
             ;;              ;;  node-mass
             ;;              ;;  (* 1d0
             ;;              ;;     (/ (* ;; (expt (cl-mpm/mesh::mesh-resolution mesh) 2)
             ;;              ;;         node-volume
             ;;              ;;         mp-stiffness) dt-scale)))
             ;;              ))))))))
             ;; (cl-mpm::iterate-over-neighbours-point-linear
             ;;  mesh
             ;;  (cl-mpm/particle::mp-penalty-contact-point mp)
             ;;  (lambda (mesh node svp grads)
             ;;    (with-accessors ((node-active cl-mpm/mesh:node-active)
             ;;                     (node-volume cl-mpm/mesh::node-volume)
             ;;                     (node-mass cl-mpm/mesh::node-mass)
             ;;                     (node-lock cl-mpm/mesh::node-lock))
             ;;        node
             ;;      (declare (double-float node-mass node-volume mp-stiffness))
             ;;      (when node-active
             ;;        (sb-thread:with-mutex (node-lock)
             ;;          (setf
             ;;           node-mass
             ;;           ;; mp-stiffness
             ;;           (max
             ;;            node-mass
             ;;            (/ (* 10 svp mp-stiffness) dt-scale)
             ;;            ;; (/ mp-stiffness dt-scale)
             ;;            )
             ;;           ))))))
             (iterate-over-neighbours
              mesh mp
              (lambda (mesh mp node svp grads fsvp fgrads)
                (with-accessors ((node-active cl-mpm/mesh:node-active)
                                 (node-volume cl-mpm/mesh::node-volume)
                                 (node-mass cl-mpm/mesh::node-mass)
                                 (node-lock cl-mpm/mesh::node-lock))
                    node
                  (declare (double-float node-mass node-volume mp-stiffness))
                  (when node-active
                    (sb-thread:with-mutex (node-lock)
                      (setf
                       node-mass
                       ;; (+
                       ;;  node-mass
                       ;;  (* svp mp-stiffness))
                       (+
                        node-mass
                        (* 1d0 (/ (* node-volume mp-stiffness) dt-scale))
                        ;; (* 1.5d0 (/ (* node-volume mp-stiffness) dt-scale))
                        )
                       ))))))
             ))))
    ))


(defun reset-penalty (sim)
  (cl-mpm:iterate-over-mps
   (sim-mps sim)
   (lambda (mp)
     (cl-mpm/fastmaths:fast-zero (cl-mpm/particle::mp-penalty-frictional-force mp))
     ;; (setf (cl-mpm/particle::mp-penalty-contact mp) nil)
     (setf (cl-mpm/particle::mp-penalty-stiffness mp) (* 0.99d0 (cl-mpm/particle::mp-penalty-stiffness mp)))
     ))
  )
