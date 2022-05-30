(defpackage :cl-mpm
  (:use :cl)
  (:import-from 
    :magicl tref .+ .-)
  (:export
    #:mpm-sim
    #:make-mpm-sim
    #:make-shape-function-linear
    #:sim-mesh
    #:sim-mps
    #:sim-bcs
    #:sim-dt
    ))
;    #:make-shape-function
(in-package :cl-mpm)

(defclass mpm-sim ()
  ((dt 
     :accessor sim-dt
     :initarg :dt)
   (mesh 
     :accessor sim-mesh
     :initarg :mesh)
   (mps 
     :accessor sim-mps
     :initarg :mps)
   (bcs
     :accessor sim-bcs
     :initarg :bcs
     :initform '())
   ))

(defun make-mpm-sim (size resolution dt shape-function)
  (make-instance 'mpm-sim
                 :dt dt
                 :mesh (make-mesh size resolution shape-function)
                 :mps '()))

(defun voight-to-matrix (vec)
  (let* ( (exx (tref vec 0 0))
          (eyy (tref vec 1 0))
          (exy (tref vec 2 0)))
    (magicl:from-list (list exx exy 
                            exy eyy)
                      '(2 2) :type 'double-float)))

(defun update-sim (sim)
  (with-slots ((mesh mesh)
               (mps mps)
               (bcs bcs)
               (dt dt))
                sim
                (progn
                    (reset-grid mesh)
                    (p2g mesh mps)
                    (filter-grid mesh 1e-3)
                    (update-nodes mesh dt)
                    (apply-bcs mesh bcs)
                    (g2p mesh mps)
                    (update-particle mps dt)
                    (update-stress mesh mps dt) 
                  )))

(defgeneric iterate-over-neighbours-shape (mesh shape-func mp func)
  (:documentation "For a given shape function iterate over relevant nodes and call func with mesh, mp, node, weight and gradients"))

(defun iterate-over-neighbours (mesh mp func)
  (iterate-over-neighbours-shape mesh (slot-value mesh 'shape-func) mp func))

(defmethod iterate-over-neighbours-shape (mesh shape-func mp func)
  (let* ((pos (slot-value mp 'position))
         (pos-x (tref pos 0 0))
         (pos-y (tref pos 1 0))
         (nD (slot-value mesh 'nD))
         (h (slot-value mesh 'mesh-res))
         (order (slot-value shape-func 'order))
         (pos-i (magicl:map! (lambda (x) (floor x)) (magicl:scale pos (/ 1 h))))
         (pos-ix (tref pos-i 0 0))
         (pos-iy (tref pos-i 1 0)))
    (loop for dx from 0 to order
          do (loop for dy from 0 to order
                   do
                   (let* (
                          (id-x (floor (+ pos-ix dx)))
                          (id-y (floor (+ pos-iy dy)))
                          (dist-x (- pos-x (* id-x h)))
                          (dist-y (- pos-y (* id-y h)))
                          (node (get-node mesh (list id-x id-y)))
                          (weight (apply (svp shape-func) (list dist-x dist-y)))
                          (grads (apply (dsvp shape-func) (list dist-x dist-y))) 
                          )
                          (funcall func mesh mp node weight (assemble-dsvp nD grads)))))))

;(defmethod iterate-over-neighbours-shape (mesh (shape-func shape-function-linear) mp func)
;  (let* ((pos (slot-value mp 'position))
;          (h (slot-value mesh 'mesh-res))
;          (order (slot-value shape-func 'order))
;          (x-mp (tref pos 0 0))
;          (y-mp (tref pos 1 0))
;          (xi (floor (/ x-mp h)))
;          (yi (floor (/ y-mp h))))
;          (loop for dx from 0 to 1
;          do (loop for dy from 0 to 1
;            do (let ( (vx (- x-mp xi dx))
;                      (vy (- y-mp yi dy))
;                      (node (get-node mesh (+ xi dx) (+ yi dy))))
;                        (multiple-value-bind (svp grads) 
;                        (funcall (slot-value shape-func 'svp) vx vy)
;                          (funcall func mesh mp node svp (assemble-dsvp grads))))))))
;


(defun p2g (mesh mps)
  (loop for mp in mps
    do (iterate-over-neighbours mesh mp
      (lambda (mesh mp node svp dsvp) 
      (progn
        (with-slots ( (node-vel velocity)
                      (node-mass mass)
                      (node-force force)) node
          (with-slots ( (mp-vel velocity)
                        (mp-mass mass)) mp 
            (progn 
                (setf node-mass 
                  (+ node-mass (* mp-mass svp)))
                (setf node-vel 
                  (magicl:.+ node-vel (magicl:scale mp-vel (* mp-mass svp))))
                (setf node-force 
                  (magicl:.- node-force  (det-int-force mp node dsvp)))
                (setf node-force 
                  (magicl:.+ node-force  (det-ext-force mp node svp)))
                         ))))))))

(defun g2p-mp-node (mesh mp node svp dsvp) 
  "Map grid to particle for one mp-node pair"
  (with-slots ((node-vel velocity)) node
    (with-slots ((vel velocity)
                 ;(strain-rate strain-rate)
                 (mass mass)) mp
      (progn 
        (setf vel (magicl:.+ vel (magicl:scale node-vel svp)))
        ;(setf strain-rate (magicl:.+ strain-rate (magicl:@ dsvp node-vel)))
        ))))
(defun g2p-mp (mesh mp)
  "Map one MP onto the grid"
  (with-slots ((vel velocity))
    mp
    (progn
      (setf vel (magicl:scale vel 0d0)))
    (iterate-over-neighbours mesh mp #'g2p-mp-node)))

(defun g2p (mesh mps)
  "Map grid values to all particles"
  (lparallel:pdotimes (i (length mps)) 
    (g2p-mp mesh (nth i mps))))

(defun g2p-serial (mesh mps)
  "Map grid values to all particles"
  (dolist (mp mps)
    (g2p-mp mesh mp)))


(defun update-node (mesh dt node)
    (when (> (slot-value node 'mass) 0)
      (with-slots ( (mass mass)
                    (vel velocity)
                    (force force)
                    (acc acceleration))
        node
        (progn 
          (setf vel (magicl:scale vel (/ 1 mass)))
          (setf acc (magicl:scale force (/ 1.0 mass)))
          (setf acc (magicl:.- acc (magicl:scale vel 0.5)))
          (setf vel (magicl:.+ vel (magicl:scale acc dt)))
          ))))

(defun update-nodes (mesh dt)
  "Explicitly integrate nodes"
  (let ((nodes (slot-value mesh 'nodes)))
  (lparallel:pdotimes (i (array-total-size nodes))
    (let ((node (row-major-aref nodes i)))
      (update-node mesh dt node)))))

(defun update-nodes-serial (mesh dt)
  "Explicitly integrate nodes"
  (let ((nodes (slot-value mesh 'nodes)))
  (dotimes (i (array-total-size nodes))
    (let ((node (row-major-aref nodes i)))
      (update-node mesh dt node)))))

(defun apply-bcs (mesh bcs)
  (with-slots ((nodes nodes)
               (nD nD)
               (mc mesh-count)) mesh
    ;each bc is a list (pos value)
    (dolist (bc bcs)
      (let ((index (cl-mpm/bc:bc-index bc)))
        (when (in-bounds mesh index)
          (cl-mpm/bc:apply-bc bc (get-node mesh index)))))))

(defun update-particle (mps dt)
  (loop for mp in mps
    do (with-slots ((vel velocity)
                    (pos position)) mp
      (progn 
      (setf pos (magicl:.+ pos (magicl:scale vel dt)))))))

;Could include this in p2g but idk
(defun calculate-strain-rate (mesh mp dt)
  (let ((dstrain (magicl:zeros '(3 1))))
    (progn
      (iterate-over-neighbours mesh mp 
       (lambda (mesh mp node svp dsvp)
         (setf dstrain (.+ dstrain (magicl:@ dsvp (node-velocity node))))))
      (magicl:scale  dstrain dt))))

(defun update-stress-mp (mesh mp dt)
    (with-slots ((stress stress)
                    (volume volume)
                    (strain strain)
                    (def deformation-matrix)) mp
         
           (let ((dstrain (calculate-strain-rate mesh mp dt)))
             (progn
               (setf strain (magicl:.+ strain dstrain))
               (setf stress (constitutive-model mp strain))
                (let ((df (.+ (magicl:eye 2) (voight-to-matrix dstrain))))
                   (progn
                     (setf def (magicl:@ df def))
                     (setf volume (* volume (magicl:det df)))
                     ))))))

(defun update-stress (mesh mps dt)
  (lparallel:pdotimes (i (length mps))
     (update-stress-mp mesh (nth i mps) dt)))
(defun update-stress-serial (mesh mps dt)
  (loop for mp in mps
    do (update-stress-mp mesh mp dt)))
