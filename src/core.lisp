(defpackage :cl-mpm
  (:use :cl )
  (:export
    #:mpm-sim
    #:make-mpm-sim
    #:make-shape-function-linear
    ))
;    #:make-shape-function
(in-package :cl-mpm)

(defclass mpm-sim ()
  ((dt 
     :accessor :get-dt
     :initarg :dt)
   (mesh 
     :accessor :get-mesh
     :initarg :mesh)
   (mps 
     :accessor :get-mps
     :initarg :mps)
   (nD 
     :accessor :get-nd
     :initarg :nD)))

(defun make-mpm-sim (nD size resolution shape-function)
  (make-instance 'mpm-sim
                 :dt 1
                 :mesh (make-mesh nD size resolution shape-function)
                 :mps '()))



(defun update-sim (sim)
  (with-slots ((mesh mesh)
               (mps mps)
               (dt dt))
                sim
                (progn
                    (reset-grid mesh)
                    (p2g mesh mps)
                    (filter-grid mesh 1e-3)
                    (update-nodes mesh dt)
                    (apply-bc mesh)
                    (g2p mesh mps)
                    (update-particle mps dt)
                    (update-stress mps dt) 
                  )))

(defgeneric iterate-over-neighbours-shape (mesh shape-func mp func)
  (:documentation "For a given shape function iterate over relevant nodes and call func with mesh, mp, node, weight and gradients"))

(defun iterate-over-neighbours (mesh mp func)
  (iterate-over-neighbours-shape mesh (slot-value mesh 'shape-func) mp func))

(defmethod iterate-over-neighbours-shape (mesh shape-func mp func)
  (let* ((pos (slot-value mp 'position))
          (h (slot-value mesh 'mesh-res))
          (nD (slot-value mesh 'nD))
          (order (slot-value shape-func 'order))
          (pos-i (magicl:map! (lambda (x) (floor x)) (magicl:scale pos (/ 1 h)))))
    (dolist (dindex (apply #'alexandria:map-product #'list
                         (loop for d from 1 to nD
                               collect (loop for v from 0 to order collect v))))
      (progn
        ;(print dindex)
        (let* ((idpos (loop for i from 0 to (- nD 1) 
                            collect (floor (+ (magicl:tref pos-i i 0) 
                                              (nth i dindex)))))
               (node (apply #'get-node (cons mesh (list idpos))))
               (vpos (magicl:.- pos 
                                (magicl:.+ pos-i (magicl:from-list dindex (list nD 1)))))
               (vpos-list (loop for i from 0 to (- nD 1) collect (magicl:tref vpos i 0)))

               (weight (apply (svp shape-func) vpos-list))
               (grads (apply (dsvp shape-func) vpos-list)))
          (funcall func mesh mp node weight (assemble-dsvp nD grads))
        )))))
;(defmethod iterate-over-neighbours-shape (mesh (shape-func shape-function-linear) mp func)
;  (let* ((pos (slot-value mp 'position))
;          (h (slot-value mesh 'mesh-res))
;          (order (slot-value shape-func 'order))
;          (x-mp (magicl:tref pos 0 0))
;          (y-mp (magicl:tref pos 1 0))
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
                ;(setf node-force 
                ;  (magicl:.- node-force  (det-int-force mp node dsvp)))
                (setf node-force 
                  (magicl:.+ node-force  (det-ext-force mp node svp)))
                         ))))))))

(defun g2p-mp (mesh mp node svp dsvp) 
  "Map grid to particle for one mp-node pair"
  (with-slots ((node-vel velocity)) node
    (with-slots ((vel velocity)
                 (strain-rate strain-rate)
                 (mass mass)) mp
      (progn 
        (setf vel (magicl:.+ vel (magicl:scale node-vel svp)))
        (setf strain-rate (magicl:.+ strain-rate (magicl:@ dsvp node-vel)))
        ))))

(defun g2p (mesh mps)
  "Map grid values to all particles"
  (dolist (mp mps)
    (progn
      (setf (slot-value mp 'velocity) (magicl:zeros '(2 1)))
      (setf (slot-value mp 'strain-rate) (magicl:zeros '(3 1)))
      (iterate-over-neighbours mesh mp #'g2p-mp))))

(defun update-nodes (mesh dt)
  "Explicitly integrate nodes"
  (let ((nodes (slot-value mesh 'nodes)))
  (dotimes (i (array-total-size nodes))
    (let ((node (row-major-aref nodes i)))
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
          )))))))

(defun apply-bc (mesh)
  (with-slots ((nodes nodes)
                (mc mesh-count)) mesh
    (dotimes (x   (nth 0 mc))
      (dotimes (y (nth 1 mc))
      (progn
        (with-slots ((vel velocity)) (get-node mesh x y)
          (when (or (<= y 0) (>= y (- (nth 1 mc) 1)))
            (setf (magicl:tref vel 1 0)
                  0d0))
          (when (or (<= x 0) (>= x (- (nth 0 mc) 1)))
            (setf (magicl:tref vel 0 0)
                  0d0))))))))

(defun update-particle (mps dt)
  (loop for mp in mps
    do (with-slots ((vel velocity)
                    (pos position)) mp
      (progn 
      (setf pos (magicl:.+ pos (magicl:scale vel dt)))))))
(defun update-stress (mps dt)
  (loop for mp in mps
    do (with-slots ((stress stress)
                    (volume volume)
                    (strain-rate strain-rate)
                    (def deformation-matrix)) mp
         
           (let ((eng-strain (magicl:scale strain-rate dt)))
             (progn
               ;Engineering strain increment
               (setf eng-strain (magicl:.* eng-strain (magicl:from-list '(1d0 1d0 1d0) '(3 1))))
               ;Update stress incrementally
               (setf stress (magicl:.+ stress (constitutive-model mp eng-strain)))
               (let* (
                      (exx (magicl:tref eng-strain 0 0))
                      (eyy (magicl:tref eng-strain 1 0))
                      (exy (magicl:tref eng-strain 2 0))
                      )
                 (let ((df (magicl:from-list (list (+ 1 exx) exy 
                                                    exy (+ 1 eyy)) '(2 2) :type 'double-float)))
                   (progn
                     (setf def (magicl:@ df def))
                     (setf volume (* volume (magicl:det df)))
                     ))))
             ))))
