(ql:quickload "magicl")
(ql:quickload "cl-autodiff")
(ql:quickload "vgplot")
(ql:quickload "array-operations")
(ql:quickload "swank.live")

;(ql:quickload "ltk")
;(ltk:with-ltk ()  (let ((b (make-instance 'ltk:button
;                                        :master nil
;                                        :text "Hello")))
;                                        (ltk:pack b)))

(defpackage :cl-mpm
  (:use :cl))
(in-package :cl-mpm)

(defclass particle ()
  ((mass :initarg :mass 
         :initform 1)
   (volume :initarg :volume 
           :initform 1)
   (position :initarg :position)
   (velocity :initarg velocity
             :initform (magicl:zeros '(2 1)))
   (stress :initform (magicl:zeros '(3 1)))
   (strain-rate :initform (magicl:zeros '(3 1)))
   (deformation-matrix :initform (magicl:eye 2)))
  (:documentation "A single material point"))

(defclass node ()
  ((mass :initform 0)
  (index :initarg :index)
  (acceleration :initform (magicl:zeros '(2 1)))
  (force :initform (magicl:zeros '(2 1)))
  (velocity :initform (magicl:zeros '(2 1))))
  (:documentation "A node on the computational mesh"))

(defclass mesh ()
  (
    (nD :initarg :nD)
    (mesh-count :initarg :mesh-count)
    (mesh-size :initarg :mesh-size)
    (mesh-res :initarg :mesh-res)
    (nodes :initarg :nodes :accessor get-nodes)
    (svp :initarg :svp))
    (:documentation "MPM computational mesh"))

(defmacro shape-linear (x)
  `(- 1d0 (abs ,x)))

(defclass shape-function ()
  ((nD :initarg :nD)
   (order :initarg :order)
   (svp :initarg :svp)))
(defclass shape-function-linear (shape-function)
  ( (order :initform 2)
   (svp :initform (autodiff:lambda-ad (x y) (* (shape-linear x) (shape-linear y))))))

(defmacro make-shape-function (shape order)
  `(make-instance 'shape-function
                 :nD 2
                 :order ,order 
                 :svp (autodiff:lambda-ad (x y) (* (,shape x) (,shape y)))))

(defun make-node (x y)
  (make-instance 'node 
                 :index (magicl:from-list (list x y) '(2 1))))

(defun make-particle (&key (x 0) (y 0) (volume 1))
  (make-instance 'particle
                 :volume volume
                 :position (magicl:from-list (list x y) '(2 1))))

(defun make-mesh (nD size resolution shape-function)
  (let* ((meshcount (loop for d in size collect (+ (floor d resolution) 1)))
         (nodes (loop for x from 0 to (- (nth 0 meshcount) 1)
                      collect (loop for y from 0 to (- (nth 1 meshcount) 1)
                      collect (make-node x y)))))
    (progn
    (make-instance 'mesh
      :nD nD
      :mesh-size size
      :mesh-count meshcount
      :mesh-res resolution
      :nodes (make-array meshcount :initial-contents nodes)
      :svp shape-function
      ))))

(defun in-bounds (mesh value dim)
  "Check a single dimension is inside a mesh"
  (and (>= value 0) (< value (nth dim (slot-value mesh 'mesh-count)))))

(defun get-node (mesh x y)
  "Check bounds and get node"
  (if (and (in-bounds mesh x 0) (in-bounds mesh y 1))
    (aref (get-nodes mesh) x y)
    (error "Access grid out of bounds")))

(defun assemble-dsvp (dsvp)
  "Assemble d/di to the strain-displacement matrix"
  (let ((nD 2)
        (dx (aref dsvp 0))
        (dy (aref dsvp 1)))
    (magicl:from-list (list dx 0d0 0d0 dy dy 0d0) '(3 2) :type 'double-float)))
(defun det-int-force (mp node dsvp)
  "Calculate internal force contribution from mp at node"
  (with-slots ( (stress stress)
                (volume volume)) mp

    (magicl:@  (magicl:transpose dsvp) (magicl:scale stress volume))
    ))
(defun det-ext-force (mp node svp)
  "Calculate external force contribution from mp at node"
  (with-slots ( (mass mass)) mp
    (magicl:scale (magicl:from-list '(0d0 -9.8d0) '(2 1)) (* mass svp))))

(defgeneric iterate-over-neighbours-shape (mesh shape-func mp func)
  (:documentation "For a given shape function iterate over relevant nodes and call func with mesh, mp, node, weight and gradients"))

(defun iterate-over-neighbours (mesh mp func)
  (iterate-over-neighbours-shape mesh (slot-value mesh 'svp) mp func))

(defmethod iterate-over-neighbours-shape (mesh (shape-func shape-function-linear) mp func)
  (let* ((pos (slot-value mp 'position))
          (h (slot-value mesh 'mesh-res))
          (order (slot-value shape-func 'order))
          (x-mp (magicl:tref pos 0 0))
          (y-mp (magicl:tref pos 1 0))
          (xi (floor (/ x-mp h)))
          (yi (floor (/ y-mp h))))
          (loop for dx from 0 to 1
          do (loop for dy from 0 to 1
            do (let ( (vx (- x-mp xi dx))
                      (vy (- y-mp yi dy))
                      (node (get-node mesh (+ xi dx) (+ yi dy))))
                        (multiple-value-bind (svp grads) 
                        (funcall (slot-value shape-func 'svp) vx vy)
                          (funcall func mesh mp node svp (assemble-dsvp grads))))))))

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
(defun reset-node (node)
  (with-slots ( (mass mass)
                 (vel velocity)
                 (force force))
                node
               (setf mass 0)
               (setf vel (magicl:scale vel 0))
               (setf force (magicl:scale force 0))))
(defun reset-grid (mesh)
  (let ((nodes (slot-value mesh 'nodes)))
  (dotimes (i (array-total-size nodes))
    (let ((node (row-major-aref nodes i)))
    (when (> (slot-value node 'mass) 0)
      (reset-node node))))))
(defun filter-grid (mesh mass-thresh)
  (let ((nodes (slot-value mesh 'nodes)))
  (dotimes (i (array-total-size nodes))
    (let ((node (row-major-aref nodes i)))
      (when (and (> (slot-value node 'mass) 0) (< (slot-value node 'mass) mass-thresh))
        (reset-node node))))))

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

(defun linear-elastic (E nu)
  (magicl:scale  (magicl:from-list (list (- 1d0 nu) nu 0d0 
                                         nu (- nu) 0d0 
                                         0d0 0d0 (- 1d0 (* 2 nu))) '(3 3) :type 'double-float) (/ E (* (+ 1 nu) (- 1 nu nu)))))

(defun constitutive-model (strain)
  (magicl:@ (linear-elastic 1e4 0.1) strain))

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
               (setf stress (magicl:.+ stress (constitutive-model eng-strain)))
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
  "Map "
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


;(defparameter *mps* (list (make-particle :x 1 :y 1)))
;(setf (slot-value (first *mps*) 'velocity) (cl-math:matrix-from-data '((2 4))))

(vgplot:figure)
(vgplot:axis (list 0 (nth 0 (slot-value *mesh* 'mesh-size)) 
                   0 (nth 1 (slot-value *mesh* 'mesh-size))))
(defun mesh-step (k substeps) 
  (dotimes (n k)
                (progn
                  (format t "Step ~d ~%" n)
                  (dotimes (i substeps)
                    (reset-grid *mesh*)
                    (p2g *mesh* *mps*)
                    ;(filter-grid *mesh* 1e-3)
                    (update-nodes *mesh* *dt*)
                    (apply-bc *mesh*)
                    (g2p *mesh* *mps*)
                    (update-particle *mps* *dt*)
                    (update-stress *mps* *dt*))
                    (print (loop for mp in *mps* collect (magicl:tref (slot-value (first *mps*) 'stress) 1 0)))
                    (format t "~%")
                  ;(print (loop for mp in *mps* collect (cl-math:[][] 0 1 (slot-value mp 'velocity)))))
                  (let* ( (pos (loop for mp in *mps* collect (slot-value mp 'position)))
                         (x (loop for p in pos collect (magicl:tref p 0 0)))
                         (y (loop for p in pos collect (magicl:tref p 1 0))))
                    (vgplot:plot x y ";;with points pt 7")
                    )
                  (vgplot:replot)
                  (swank.live:update-swank)
                  (sleep .1)
                  )))

(defparameter *mesh* (make-mesh 2 '(1 10) 1 (make-instance 'shape-function-linear)))
(defparameter *mps* (loop for id from 0 to 1 collect (make-particle :x 0.5 :y (+ (* 0.5 id) 4) :volume 0.5)))
(defparameter *mps* (list (make-particle :x 0.5 :y 5)))
(defparameter *dt* 1e-3)
(mesh-step 100 10)
(setf (slot-value (nth 2 *mps*) 'velocity) (magicl:from-list '(100d0 10d0) '(2 1)))
(print (slot-value  (first *mps*) 'deformation-matrix))
(print (slot-value  (first *mps*) 'strain-rate))

(dotimes (n 11)
  (print (magicl:tref (slot-value (get-node *mesh* 0 n) 'force) 1 0)))
