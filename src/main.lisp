;(ql:quickload "l-math")
(ql:quickload "cl-autodiff")
(ql:quickload "cl-math")
(ql:quickload "magicl")
(ql:quickload "vgplot")
(ql:quickload "array-operations")
(ql:quickload "ltk")

(ltk:with-ltk ()  (let ((b (make-instance 'button
                                        :master nil
                                        :text "Hello")))
                                        (ltk:pack b)))
(defpackage :cl-mpm
  (:use :cl))
(in-package :cl-mpm)

(defun 2d-array-to-list (array)
  (loop for i below (array-dimension array 0)
        collect (loop for j below (array-dimension array 1)
                      collect (aref array i j))))

(defclass particle ()
  ((mass :initarg :mass :initform 1)
   (volume :initarg :volume :initform 1)
   (position :initarg :position)
   (velocity :initform (cl-math:zeros 2))
   (stress :initform (cl-math:zeros 3))
   (strain-rate :initform (cl-math:zeros 3))
   (deformation-matrix :initform :deformation-matrix)))

(defclass node ()
  ((mass :initform 0)
  (index :initarg :index)
  (acceleration :initform (cl-math:zeros 2))
  (force :initform (cl-math:zeros 2))
  (velocity :initform (cl-math:zeros 2))))

(defclass mesh ()
  (
    (nD :initarg :nD)
    (mesh-count :initarg :mesh-count)
    (mesh-size :initarg :mesh-size)
    (mesh-res :initarg :mesh-res)
    (nodes :initarg :nodes :accessor get-nodes)
    (svp :initarg :svp)))

(defmacro shape-linear (x)
  `(- 1 (abs ,x)))
(defun make-node (x y)
  (make-instance 'node 
  :index (cl-math:matrix-from-data '((x y)))))
(defun make-particle (&key (x 0) (y 0))
  (make-instance 'particle
    :position (cl-math:matrix-from-data (list (list x y)))))
          ;(nodes (loop for x from 0 to (- mesh-count 1)
          ;          collect (loop for y from 0 to (- mesh-count 1)
          ;            collect (make-node x y)))))
(defun make-mesh (nD size resolution)
  (let* ( (meshcount (floor (/ size resolution)))
          (nodes (loop for x from 0 to (- meshcount 1)
                    collect (loop for y from 0 to (- meshcount 1)
                      collect (make-node x y)))))
    (progn
    (make-instance 'mesh
      :nD nD
      :mesh-size size
      :mesh-count meshcount
      :mesh-res resolution
      :nodes (make-array (list meshcount meshcount) :initial-contents nodes)
      :svp (autodiff:lambda-ad (x y) (* (shape-linear x) (shape-linear y)))
      ))))
(defun in-bounds (mesh dim)
  (and (>= dim 0) (< dim (slot-value mesh 'mesh-count))))
(defun get-node (mesh x y)
  (if (and (in-bounds mesh x) (in-bounds mesh y))
    (get-node mesh x y)
    (error "Access grid out of bounds")))
  
(defparameter *test-svp* (autodiff:lambda-ad (x y) (* (shape-linear x) (shape-linear y))))
(defun assemble-dsvp (dsvp)
  (let ((nD 2)
        (dx (aref dsvp 0))
        (dy (aref dsvp 1)))
    (cl-math:matrix-from-data (list (list dx 0) (list 0 dy) (list dy 0)))))
(defparameter *dv* (cl-math:matrix-from-data '((0 1))))
(defparameter *dsvp* (assemble-dsvp (make-array 2 :initial-contents '(1 1))))
*dsvp*
(cl-math:transpose *dv*)
(cl-math:dot *dsvp* (cl-math:transpose *dv*))

(defun det-int-force (mp node dsvp)
  (with-slots ( (stress stress)
                (volume volume)) mp
    (cl-math:dot (cl-math:*mv stress volume) dsvp)))
(defun det-ext-force (mp node svp)
  (with-slots ( (mass mass)) mp
    (cl-math:*mv (cl-math:matrix-from-data '((0 -9.8))) (* mass svp))))
(defun iterate-over-neighbours (mesh mp func)
  (let* ((pos (slot-value mp 'position))
          (h (slot-value mesh 'mesh-res))
          (xi (floor (/ (cl-math:[][] 0 0 pos) h)))
          (yi (floor (/ (cl-math:[][] 0 1 pos) h))))
          (loop for dx from 0 to 1
          do (loop for dy from 0 to 1
            do (let ( (vx (- (cl-math:[][] 0 0 pos) xi dx))
                      (vy (- (cl-math:[][] 0 1 pos) yi dy))
                      (node (get-node mesh (+ xi dx) (+ yi dy))))
                        (multiple-value-bind (svp grads) 
                        (funcall (slot-value mesh 'svp) vx vy)
                          (funcall func mesh mp node svp (assemble-dsvp grads))))))))
(defun apply-bc (mesh)
  (with-slots ((nodes nodes)
                (mc mesh-count)) mesh
  (dotimes (x mc)
    (dotimes (y mc)
      (progn
      (when (or (= 0 y) (= (- mc 1) y))
        (setf (cl-math:[][] 0 1 (slot-value (get-node mesh x y) 'velocity))
          0))
      (when (or (= 0 x) (= (- mc 1) x))
        (setf (cl-math:[][] 0 0 (slot-value (get-node mesh x y) 'velocity))
          0)))))))
(defun reset-grid (mesh)
  (let ((nodes (slot-value mesh 'nodes)))
  (dotimes (i (array-total-size nodes))
    (let ((node (row-major-aref nodes i)))
    (when (> (slot-value node 'mass) 0)
      (with-slots ( (mass mass)
                    (vel velocity)
                    (force force))
        node
        (progn 
        (setf mass 0)
        (setf vel (cl-math:*mv vel 0))
        (setf force (cl-math:*mv force 0)))))))))

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
            (format t "svp ~d" svp)
            (setf node-mass 
              (+ node-mass (* mp-mass svp)))
            (setf node-vel 
              (cl-math:+mm node-vel (cl-math:*mv mp-vel (* mp-mass svp))))
            ;(setf node-force 
            ;  (cl-math:+mm node-force  (det-int-force mp node dsvp)))
            (setf node-force 
              (cl-math:+mm node-force  (det-ext-force mp node svp)))
                         ))))))))
(defun constitutive-model (stress strain)
  (cl-math:matrix-from-data '((0 0 0))))
(defun update-particle (mps dt)
  (loop for mp in mps
    do (with-slots ((vel velocity)
                    (pos position)) mp
      (progn 
      (setf pos (cl-math:+mm pos (cl-math:*mv vel dt)))))))
                        
(defun update-stress (mps dt)
  (loop for mp in mps
    do (with-slots ((stress stress)
                    (strain strain)
                    (def deformation-matrix)) mp
      (progn 
      (setf stress (constitutive-model stress strain)))
      (let ((ddf (cl-math:matrix-from-data '((1 0) (0 1)))))
      (setf def (cl-math:*mm ddf def)))
    )))
(defun integrate-node (node dt)
  (with-slots (acceleration velocity) node
    (setf velocity (cl-math:+mm velocity (cl-math:*mv acceleration dt)))))

;(defparameter *node* (make-instance 'node))
;(describe *node*)
;(setf (lm:x (slot-value *node* 'acceleration)) 1)




(defun g2p (mesh mps)
  (dolist (mp mps)
      (let ((strain (cl-math:zeros 2))
            (vel (slot-value mp 'velocity)))
            (progn
              (setf vel (cl-math:*mv vel 0))
              (iterate-over-neighbours mesh mp
              (lambda (mesh mp node svp dsvp) 
                (with-slots ((node-vel velocity)
                              (strain-rate strain-rate)) node
                  (with-slots ((vel velocity)
                               (mass mass)) mp
                  (setf vel (cl-math:+mm vel (cl-math:*mv node-vel svp)))
                  (let ((partial-strain (cl-math:dot dsvp (cl-math:transpose node-vel))))
                  (progn
                  (setf strain (cl-math:+mm strain partial-strain))
                  ))))))))))
(g2p *mesh* *mps*)
(defun update-nodes (mesh dt)
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
          (setf vel (cl-math:/mv vel mass))
          (setf acc (cl-math:*mv force (/ 1.0 mass)))
          (setf vel (cl-math:+mm vel (cl-math:*mv acc dt)))
          )))))))

(defparameter *mesh* (make-mesh 2 3 1))
(defparameter *mps* (loop for x from 0 to 1 collect (make-particle)))
(defparameter *mps* '())
(push (make-particle :x 0.5 :y 1) *mps*)
(defparameter *dt* 0.01)

(describe *mesh*)
(describe (first *mps*))
(vgplot:figure)
(vgplot:axis (list 0 (slot-value *mesh* 'mesh-size) 
                   0 (slot-value *mesh* 'mesh-size)))
(reset-grid *mesh*)
(aops:flatten (get-nodes *mesh*))
(aops:each (lambda (x) (cl-math:[][] 0 1 (slot-value x 'velocity)))
  (aops:flatten (get-nodes *mesh*)))
(aops:each (lambda (x) (cl-math:[][] 0 1 (slot-value x 'force)))
  (aops:flatten (get-nodes *mesh*)))

(print (loop for mp in *mps* collect (cl-math:[][] 0 1 (slot-value mp 'velocity))))
(print (loop for mp in *mps* collect (cl-math:[][] 0 1 (slot-value mp 'position))))
(reset-grid *mesh*)
(p2g *mesh* *mps*)
(update-nodes *mesh* *dt*)
(apply-bc *mesh*)
(g2p *mesh* *mps*)
(update-particle *mps* *dt*)
(let* ( (pos (loop for mp in *mps* collect (slot-value mp 'position)))
        (x (loop for p in pos collect (cl-math:[][] 0 0 p)))
        (y (loop for p in pos collect (cl-math:[][] 0 1 p))))
      (vgplot:plot x y ";;with points pt 7")
      )


*mps*
(progn 
(dotimes (n 1)
  (reset-grid *mesh*)
  (p2g *mesh* *mps*)
  (update-nodes *mesh* *dt*)
  (apply-bc *mesh*)
  (g2p *mesh* *mps*)
  (update-particle *mps* *dt*))
(print (loop for mp in *mps* collect (cl-math:[][] 0 1 (slot-value mp 'velocity))))
(print (loop for mp in *mps* collect (cl-math:[][] 0 1 (slot-value mp 'position))))
(let* ( (pos (loop for mp in *mps* collect (slot-value mp 'position)))
        (x (loop for p in pos collect (cl-math:[][] 0 0 p)))
        (y (loop for p in pos collect (cl-math:[][] 0 1 p))))
      (vgplot:plot x y ";;with points pt 7")
      (vgplot:replot)
      ))
(dotimes (n 5)
  (progn
    (format t "Step ~d" n)
    (reset-grid *mesh*)
    (p2g *mesh* *mps*)
    (apply-bc *mesh*)
    (update-nodes *mesh* *dt*)
    (g2p *mesh* *mps*)
    ;(setf (cl-math:[][] 0 0 (slot-value (first *mps*) 'position)) n)
    (let* ( (pos (loop for mp in *mps* collect (slot-value mp 'position)))
            (x (loop for p in pos collect (cl-math:[][] 0 0 p)))
            (y (loop for p in pos collect (cl-math:[][] 0 1 p))))
          (vgplot:plot x y ";;with points pt 7")
          )
    (vgplot:replot)
    (vgplot:format-plot nil "e\n")
    (sleep 0.1)
    ))