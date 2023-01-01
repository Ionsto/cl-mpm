(defpackage :cl-mpm
  (:use :cl
        :cl-mpm/particle
        :cl-mpm/mesh
        :cl-mpm/utils
;        :cl-mpm/shape-function
        )
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
    #:sim-damping-factor
    #:sim-mass-filter
    #:post-stress-step
    ))
(declaim (optimize (debug 0) (safety 0) (speed 3)))
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
   (bcs-force
    :accessor sim-bcs-force
    :initarg :bcs-force
    :initform '())
   (damping-factor
     :type double-float
     :accessor sim-damping-factor
     :initarg :damping-factor
     :initform 1d0)
   (mass-filter
     :type double-float
     :accessor sim-mass-filter
     :initarg :mass-filter
     :initform 0.01d0))
  (:documentation "A self contained mpm simulation"))

(defun make-mpm-sim (size resolution dt shape-function)
  (make-instance 'mpm-sim
                 :dt (coerce dt 'double-float)
                 :mesh (make-mesh size resolution shape-function)
                 :mps '()))
(defun check-mps (mps)
  "Function to check that stresses and positions are sane"
  (loop for mp across mps
        do
        (progn
          (with-accessors ((vel cl-mpm/particle::mp-velocity)) mp
            (loop for i from 0 to 1
                 do (progn
                      (when (equal (abs (magicl:tref vel i 0)) #.sb-ext:double-float-positive-infinity)
                        (print mp)
                        (error "Infinite velocity found"))
                      (when (> (abs (magicl:tref vel i 0)) 1e10)
                        (print mp)
                        (error "High velocity found"))))))))

(defgeneric update-sim (sim)
  (:documentation "Update an mpm simulation by one timestep"))

(defmethod update-sim (sim)
  (declare (cl-mpm::mpm-sim sim))
  (with-slots ((mesh mesh)
               (mps mps)
               (bcs bcs)
               (bcs-force bcs-force)
               (dt dt))
                sim
                (progn
                    (reset-grid mesh)
                    (p2g mesh mps)
                    (when (> (sim-mass-filter sim) 0d0)
                      (filter-grid mesh (sim-mass-filter sim)))
                    (apply-bcs mesh bcs-force dt)
                    (update-nodes mesh dt (sim-damping-factor sim))
                    (apply-bcs mesh bcs dt)
                    (g2p mesh mps)
                    (update-particle mps dt)
                    (update-stress mesh mps dt) 
                    (check-mps mps)
                    ;(remove-material-damaged sim)
                    )))

(defgeneric iterate-over-neighbours-shape (mesh shape-func mp func)
  (:documentation "For a given shape function iterate over relevant nodes and call func with mesh, mp, node, weight and gradients"))

(declaim
; (inline iterate-over-neighbours)
         (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle function) (values)) iterate-over-neighbours))
(defun iterate-over-neighbours (mesh mp func)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  ;(iterate-over-neighbours-shape-bspline mesh mp func)
  ;(iterate-over-neighbours-shape-linear mesh mp func)
  (iterate-over-neighbours-shape-gimp mesh mp func)
  ;; (iterate-over-neighbours-shape mesh (cl-mpm/mesh:mesh-shape-func mesh) mp func)
  (values)
  )

(defmethod iterate-over-neighbours-shape (mesh (shape-func cl-mpm/shape-function:shape-function-linear) mp func)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/shape-function::shape-function-linear shape-func)
           (cl-mpm/particle:particle mp))
  (let* ((nd (cl-mpm/mesh:mesh-nd mesh))
         (h (cl-mpm/mesh:mesh-resolution mesh))
         ;(order (slot-value shape-func 'order))
         (pos-vec (cl-mpm/particle:mp-position mp));Material point position 
         (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
         (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec 'floor)))
    (loop for dx from 0 to 1
          do (loop for dy from 0 to 1
                   do
                   (let* ((id (mapcar #'+ pos-index (list dx dy))))
                     (when (cl-mpm/mesh:in-bounds mesh id)
                       (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                              (node (cl-mpm/mesh:get-node mesh id))
                              (weight (apply (cl-mpm/shape-function::svp shape-func) dist))
                              (grads (apply (cl-mpm/shape-function::dsvp shape-func) dist))) 
                         (funcall func mesh mp node weight (cl-mpm/shape-function::assemble-dsvp nd grads) grads))))))))

(defmethod iterate-over-neighbours-shape (mesh (shape-func cl-mpm/shape-function:shape-function-bspline) mp func)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/shape-function::shape-function-bspline shape-func)
           (cl-mpm/particle:particle mp))
  (let* ((nd (cl-mpm/mesh:mesh-nd mesh))
         (h (cl-mpm/mesh:mesh-resolution mesh))
         (pos-vec (cl-mpm/particle:mp-position mp));Material point position 
         (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
         (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec 'round))
         (border (not (and (cl-mpm/mesh:in-bounds (mapcar #'+ pos-index (list 1 1)))
                          (cl-mpm/mesh:in-bounds (mapcar #'- pos-index (list 1 1)))
                          ))))
    (if border
        ;; Do something janky if we are in a border element
        ()
        (loop for dx from -1 to 1
              do (loop for dy from -1 to 1
                       do
                          (let* (
                                 (id (mapcar #'+ pos-index (list dx dy)))
                                 (dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                 (weight (apply (cl-mpm/shape-function:svp shape-func) dist))
                                 (grads (apply (cl-mpm/shape-function:dsvp shape-func) dist))
                                 )
                            (when (cl-mpm/mesh:in-bounds mesh id)
                              (funcall func mesh mp (cl-mpm/mesh:get-node mesh id) weight
                                       (cl-mpm/shape-function:assemble-dsvp nd grads)))))))))
;; (declaim (inline dispatch)
;;          (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle function) (values)) iterate-over-neighbours-shape-linear))
(defmacro iterate-over-neighbours-general (mesh mp order body)
  `(progn
    (let* ((h (cl-mpm/mesh:mesh-resolution ,mesh))
           (pos-vec (cl-mpm/particle:mp-position ,mp))
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
           (pos-index (magicl:scale pos-vec (/ 1 h)))
           )
      (loop for dx from (- ,order) to ,order
            do (loop for dy from (- ,order) to ,order
                     do (let* ((id (list (+ (round (/ (tref pos-index 0 0) h)) dx)
                                         (+ (round (/ (tref pos-index 1 0) h)) dy))))
                          (when (cl-mpm/mesh:in-bounds ,mesh id)
                            (funcall ,body
                                     ,mesh
                                     ,mp
                                     id
                                     (mapcar (lambda (p i) (- p (* i h))) pos id)))))))))
(defmacro iterate-over-neighbours-general-pure (mesh mp order args body)
  (let ((h (gensym))
        (pos-vec (gensym))
        (pos (gensym))
        (pos-index (gensym))
        (id (gensym))
        (dx (gensym))
        (dy (gensym))
        (dist (gensym)))
    `(progn
       (let* ((,h (cl-mpm/mesh:mesh-resolution ,mesh))
              (,pos-vec (cl-mpm/particle:mp-position ,mp))
              (,pos (list (tref ,pos-vec 0 0) (tref ,pos-vec 1 0)))
              (,pos-index (magicl:scale ,pos-vec (/ 1 ,h)))
              )
         (loop for ,dx from (- ,order) to ,order
               do (loop for ,dy from (- ,order) to ,order
                        do (let* ((,id (list (+ (round (tref ,pos-index 0 0)) ,dx)
                                             (+ (round (tref ,pos-index 1 0)) ,dy))))
                             (when (cl-mpm/mesh:in-bounds ,mesh ,id)
                               (let ((,dist (mapcar (lambda (p i) (- p (* i ,h))) ,pos ,id)))
                                 (symbol-macrolet ((,(first args) ,mesh)
                                                   (,(second args) ,mp)
                                                   (,(third args) ,id)
                                                   (,(fourth args) ,dist))
                                   ,body))))))))))
(defmacro ion-linear (number body)
  (let ((test (expt number 2)))
    `(+ ,test ,body))
  ;; `(with-gensyms (exp-number-name)
  ;;    (let (exp-number-name ,)))
  )
(declaim (inline iterate-over-neighbours-shape-linear))
(defun iterate-over-neighbours-shape-linear (mesh mp func)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (pos-vec (cl-mpm/particle:mp-position mp))
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
           (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec 'floor)))
      (loop for dx from 0 to 1
            do (loop for dy from 0 to 1
                     do (let* ((id (mapcar #'+ pos-index (list dx dy))))
                          (when (cl-mpm/mesh:in-bounds mesh id)
                            (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                   (node (cl-mpm/mesh:get-node mesh id))
                                   (weights (mapcar (lambda (x) (cl-mpm/shape-function::shape-linear x h)) dist))
                                   (weight (reduce #'* weights))
                                   (grads (mapcar (lambda (d w) (* (cl-mpm/shape-function::shape-linear-dsvp d h) w))
                                                  dist (nreverse weights)))
                                   (dsvp (cl-mpm/shape-function::assemble-dsvp-2d grads)))
                              (funcall func mesh mp node weight dsvp)))))))))

;(declaim (inline iterate-over-neighbours-shape-gimp))
(defun iterate-over-neighbours-shape-gimp (mesh mp func)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (pos-vec (cl-mpm/particle:mp-position mp))
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
           (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec 'round))
           (domain (list (* 0.5d0 (magicl:tref (cl-mpm/particle::mp-domain-size mp) 0 0))
                         (* 0.5d0 (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0)) )))
      (loop for dx from -1 to 1
            do (loop for dy from -1 to 1
                     do (let* ((id (mapcar #'+ pos-index (list dx dy))))
                          (when (cl-mpm/mesh:in-bounds mesh id)
                            (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                   (node (cl-mpm/mesh:get-node mesh id))
                                   (weights (mapcar (lambda (x l)
                                                      (cl-mpm/shape-function::shape-gimp x l h)) dist domain))
                                   (weight (reduce #'* weights))
                                   (grads (mapcar (lambda (d l w) (* (cl-mpm/shape-function::shape-gimp-dsvp d l h)
                                                                     w))
                                                  dist domain (nreverse weights)))
                                   (dsvp (cl-mpm/shape-function::assemble-dsvp-2d grads)))
                              ;; (format t "~a ~%" weights)
                              ;; (format t "~a ~%" grads)
                              (funcall func mesh mp node weight dsvp grads)))))))))

(declaim (inline iterate-over-neighbours-shape-bspline))
(defun iterate-over-neighbours-shape-bspline (mesh mp func)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (pos-vec (cl-mpm/particle:mp-position mp))
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
           ;; (pos-index (magicl:scale pos-vec (/ 1 h)))
           (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec 'round))
           (border (not (and (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list 2 2)))
                             (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list -2 -2)))
                             (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list 2 -2)))
                             (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list -2 2)))
                             ))))
      (if nil ;border
          (let (
                (knots-x
                  (loop for dy from -1 to 1
                        collect
                        (loop for dx from -4 to 4
                              collect (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list dx dy))))))
                (knots-y
                  (loop for dx from -1 to 1
                        collect
                        (loop for dy from -4 to 4
                              collect (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list dx dy))))))
                )
            (loop for dx from -1 to 1
                  do (loop for dy from -1 to 1
                           do (let* ((id (mapcar #'+ pos-index (list dx dy))))
                                (when (cl-mpm/mesh:in-bounds mesh id)
                                  ;; Ill defined bspline
                                  (when (and (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list 0 dy)))
                                             (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list dx 0)))
                                             )
                                    (let* ((knot-list (list (nth (+ dy 1) knots-x)
                                                            (nth (+ dx 1) knots-y)))
                                           (local-id-list (list (+ dx 0) (+ dy 0)))
                                           (dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh pos-index)))
                                           (node (cl-mpm/mesh:get-node mesh id))
                                           (weights (mapcar (lambda (k x n)
                                                              (cl-mpm/shape-function::nodal-bspline k x n h))
                                                            knot-list dist local-id-list))
                                           (weight (reduce #'* weights))
                                           (grads (mapcar (lambda (k d n w)
                                                            (* (cl-mpm/shape-function::nodal-bspline-dsvp k d n h) w))
                                                          knot-list
                                                          dist
                                                          local-id-list
                                                          (nreverse weights)))
                                           (dsvp (cl-mpm/shape-function:assemble-dsvp-2d grads)))
                                      (funcall func mesh mp node weight dsvp))))))))
          (loop for dx from -1 to 1
                do (loop for dy from -1 to 1
                         do (let* ((id (mapcar #'+ pos-index (list dx dy))))
                              (when (cl-mpm/mesh:in-bounds mesh id)
                                (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                       (node (cl-mpm/mesh:get-node mesh id))
                                       (weights (mapcar (lambda (x) (cl-mpm/shape-function:shape-bspline x h)) dist))
                                       (weight (reduce #'* weights))
                                       (grads (mapcar (lambda (d w)
                                                        (* (cl-mpm/shape-function:shape-bspline-dsvp d h) w))
                                                      dist (nreverse weights)))
                                       (dsvp (cl-mpm/shape-function:assemble-dsvp-2d grads)))
                                  (funcall func mesh mp node weight dsvp))))))))))
(defmacro iterate-over-neighbours-shape-linear-macro (mesh mp func-body)
  `(progn
    (let* ((h (cl-mpm/mesh:mesh-resolution ,mesh))
           (pos-vec (cl-mpm/particle:mp-position ,mp))
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
           (pos-index (magicl:scale pos-vec (/ 1 h)))
           ;; (pos-index (cl-mpm/mesh:position-to-index mesh pos-vec 'floor)))
           )
      (loop for dx from 0 to 1
            do (loop for dy from 0 to 1
                     do (let* (
                               (id (list (+ (floor (tref pos-index 0 0)) dx)
                                         (+ (floor (tref pos-index 1 0)) dy)
                                         ))
                               (dist (mapcar (lambda (p i) (- p (* i h))) pos id))
                               (node (cl-mpm/mesh:get-node ,mesh id))
                               (weights (mapcar #'shape-linear dist))
                               (weight (reduce #'* weights))
                               (dsvp (mapcar (lambda (d w) (* (shape-linear-dsvp d) w)) dist (nreverse weights)))
                               ;; (dsvp (assemble-dsvp-2d (list (* (shape-linear-dsvp (first dist)) wy)
                               ;;                               (* (shape-linear-dsvp (second dist)) wx))
                               ;;        ))
                               )
                          ,func-body
                          )
                     )))))


;(defparameter *mesh* (sim-mesh *sim*))
;(defparameter *sf* (cl-mpm/mesh:mesh-shape-func *mesh*))
;(defparameter *mp* (cl-mpm/particle:make-particle 2 :pos '(0.5 0.1)))
;(iterate-over-neighbours-shape *mesh* *sf* *mp* 
;                               (lambda (mesh mp node w dsvp) (print (cl-mpm/mesh:node-index  node))))
(declaim
; (inline p2g-mp-node)
 (ftype (function (cl-mpm/particle:particle
                           cl-mpm/mesh::node
                           double-float
                           magicl:matrix/double-float) (values))
                p2g-mp-node
                          ))
(defun p2g-mp-node (mp node svp dsvp)
        (declare
         (cl-mpm/particle:particle mp)
         (cl-mpm/mesh::node node)
         (double-float svp)
         )
  (with-accessors ((node-vel   cl-mpm/mesh:node-velocity)
                   (node-mass  cl-mpm/mesh:node-mass)
                   (node-temp  cl-mpm/mesh:node-temperature)
                   (node-dtemp  cl-mpm/mesh:node-dtemp)
                   (node-volume  cl-mpm/mesh::node-volume)
                   (node-force cl-mpm/mesh:node-force)
                   (node-lock  cl-mpm/mesh:node-lock)) node
    (with-accessors ((mp-vel  cl-mpm/particle:mp-velocity)
                     (mp-mass cl-mpm/particle:mp-mass)
                     (mp-temp cl-mpm/particle::mp-temperature)
                     (mp-heat-capaciy cl-mpm/particle::mp-heat-capacity)
                     (mp-thermal-conductivity cl-mpm/particle::mp-thermal-conductivity)
                     (mp-volume cl-mpm/particle:mp-volume)
                     ;; (strain-rate cl-mpm/particle:mp-strain-rate)
                     ) mp 
            ;; (declare (double-float mp-mass
            ;;                        node-mass))
            (progn
              (let* ((weighted-mass (* mp-mass svp))
                     (weighted-volume (* mp-volume svp))
                     (weighted-temp (* mp-temp mp-mass svp))
                     (weighted-dtemp (* (/ mp-volume
                                           (* mp-mass
                                              mp-heat-capaciy))
                                        mp-thermal-conductivity
                                        mp-temp
                                        (magicl::sum
                                         (magicl:.* dsvp dsvp))))
                       (weighted-vel (magicl:scale mp-vel weighted-mass))
                       (force (magicl:.- (det-ext-force mp node svp) (det-int-force mp node dsvp)))
                       )
                  (sb-thread:with-mutex (node-lock)
                    (setf node-mass
                          (+ node-mass weighted-mass))
                    (setf node-temp
                          (+ node-temp weighted-temp))
                    (setf node-dtemp
                          (+ node-dtemp weighted-dtemp))
                    (setf node-volume
                          (+ node-volume weighted-volume))
                    (setf node-vel
                          (magicl:.+ node-vel weighted-vel))
                    (setf node-force
                          (magicl:.+ node-force force))
                    ))
              )))
  (values))
(declaim (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle) (values))
                p2g-mp))
(defun p2g-mp (mesh mp)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (iterate-over-neighbours
   mesh mp
      (lambda (mesh mp node svp dsvp grads) 
        (p2g-mp-node mp node svp dsvp)))
  (values))
;; (declaim (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle) (values))
;;                 p2g-mp-new))
;; (defun p2g-mp-new (mesh mp)
;;   (declare (cl-mpm/mesh::mesh mesh)
;;            (cl-mpm/particle:particle mp))
;;   (iterate-over-neighbours-shape-linear
;;    mesh mp
;;       (lambda (mesh mp node svp dsvp) 
;;         ;; (p2g-mp-node mp node svp dsvp)
;;         ))
;;   (values))

(defun p2g (mesh mps)
  (declare (type (array cl-mpm/particle:particle) mps) (cl-mpm/mesh::mesh mesh))
  (lparallel:pdotimes (i (length mps)) 
    (p2g-mp mesh (aref mps i))))

(declaim
; (inline g2p-mp-node)
 (ftype (function (cl-mpm/mesh::mesh
                           cl-mpm/particle:particle
                           cl-mpm/mesh::node double-float magicl:matrix/double-float) (values))
                g2p-mp-node
                          ))
(defun g2p-mp-node (mesh mp node svp dsvp grads) 
  (declare (ignore mesh))
  "Map grid to particle for one mp-node pair"
  (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
                   (node-temp cl-mpm/mesh:node-temperature)) node
    (with-accessors ((vel mp-velocity)
                     (mass mp-mass)
                     (temp cl-mpm/particle::mp-temperature)
                     (strain-rate cl-mpm/particle:mp-strain-rate)) mp
      (progn 
        (setf vel (magicl:.+ vel (magicl:scale node-vel svp)))
        (setf temp (+ temp (* node-temp svp)))
        ;; (setf strain-rate (magicl:.+ strain-rate (magicl:@ dsvp node-vel)))
        ))))

(declaim (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle) (values))
                g2p-mp))
(defun g2p-mp (mesh mp)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  "Map one MP onto the grid"
  (with-accessors ((vel mp-velocity)
                   (temp cl-mpm/particle::mp-temperature)
                   (strain-rate mp-strain-rate))
    mp
    (progn
      (setf vel (magicl:scale vel 0d0))
      (setf temp 0d0)
      ;; (setf strain-rate (magicl:scale strain-rate 0d0))
      )
    ;; (iterate-over-neighbours mesh mp)
    (iterate-over-neighbours mesh mp #'g2p-mp-node)))

(defun g2p (mesh mps)
  (declare (cl-mpm/mesh::mesh mesh) (array mps))
  "Map grid values to all particles"
  (lparallel:pdotimes (i (length mps)) 
    (g2p-mp mesh (aref mps i))))

(defun update-node (mesh dt node damping)
    (when (> (cl-mpm/mesh:node-mass node) 0d0)
      (with-accessors ( (mass  node-mass)
                        (vel   node-velocity)
                        (force node-force)
                        (temp node-temperature)
                        (dtemp node-dtemp)
                        (acc   node-acceleration))
        node
        (progn 
          (setf vel (magicl:scale vel (/ 1.0 mass)))
          (setf temp (+ (/ temp mass) (* dtemp dt)))
          ;(setf temp (+ (/ temp mass)))
          (setf acc (magicl:scale
                     (magicl:.- force
                                (magicl:scale vel damping))
                     (/ 1.0 mass)))
          (setf vel (magicl:.+ vel (magicl:scale acc dt)))
          ))))

(defun update-nodes (mesh dt damping)
  "Explicitly integrate nodes"
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
  (lparallel:pdotimes (i (array-total-size nodes))
    (let ((node (row-major-aref nodes i)))
      (update-node mesh dt node damping)))))

(defun apply-bcs (mesh bcs dt)
  (declare (cl-mpm/mesh::mesh mesh))
  (with-accessors ( (nodes  mesh-nodes)
                    (nD     mesh-nD)
                    (mc     mesh-count)) mesh
    ;each bc is a list (pos value)
    (dolist (bc bcs)
      (when bc
        (let ((index (cl-mpm/bc:bc-index bc)))
          (when (cl-mpm/mesh:in-bounds mesh index);Potentially throw here
            (cl-mpm/bc:apply-bc bc (cl-mpm/mesh:get-node mesh index) mesh dt)))))))


;; (defun update-particle (mps dt)
;;   (declare (array mps))
;;   (loop for mp across mps
;;     do (with-accessors ((vel cl-mpm/particle:mp-velocity)
;;                         (pos cl-mpm/particle:mp-position)) mp
;;       (progn 
;;       (setf pos (magicl:.+ pos (magicl:scale vel dt)))))))
(declaim (inline update-particle-mp))
(defun update-particle-mp (mp dt)
  (declare (cl-mpm/particle:particle mp) (double-float dt))
  (with-accessors ((vel cl-mpm/particle:mp-velocity)
                   (pos cl-mpm/particle:mp-position)) mp
      (progn 
        (setf pos (magicl:.+ pos (magicl:scale vel dt))))))

(defun update-particle (mps dt)
  (declare (array mps) (double-float dt))
  (lparallel:pdotimes (i (length mps))
     (update-particle-mp (aref mps i) dt)))

;Could include this in p2g but idk
(defun calculate-strain-rate (mesh mp dt)
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt))
  (with-accessors ((strain-rate cl-mpm/particle:mp-strain-rate)
                   (vorticity cl-mpm/particle:mp-vorticity)
                   ) mp
        (progn
          (setf strain-rate (magicl:scale strain-rate 0d0))
          (setf vorticity (magicl:scale vorticity 0d0))
            (iterate-over-neighbours mesh mp 
                (lambda (mesh mp node svp dsvp grads)
                  ;;(setf strain-rate (.+ strain-rate (magicl:@ dsvp (cl-mpm/mesh:node-velocity node))))
                  (let* ((curl (cl-mpm/shape-function::assemble-vorticity-2d grads)))
                      (setf strain-rate (.+ strain-rate (magicl:@ dsvp (cl-mpm/mesh:node-velocity node))))
                      (setf vorticity (.+ vorticity (magicl:@ curl (cl-mpm/mesh:node-velocity node))))
                    )
                  ;; (let* ((dv-voigt (magicl:@ dsvp (cl-mpm/mesh:node-velocity node)))
                  ;;       (dv (voight-to-matrix dv-voigt))
                  ;;       )
                  ;;   (setf strain-rate (.+ strain-rate (matrix-to-voight (.+ dv (magicl:transpose dv)))))
                  ;;   (setf vorticity (.+ vorticity (matrix-to-voight (.- dv (magicl:transpose dv))))))
                  ))
            (setf strain-rate (magicl:scale strain-rate dt))
            (setf vorticity (magicl:scale vorticity dt)))))


(defun rotation-matrix (degrees)
  (let ((angle (/ (* pi degrees) 180)))
    (magicl:from-list (list (cos angle) (- (sin angle))
                            (sin angle) (cos angle))
                      '(2 2)
                      :type 'double-float
                      )))
(defun update-strain-linear (mp dstrain)
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   ) mp
    (let ((df (.+ (magicl:eye 2) (voight-to-matrix dstrain))))
                   (progn
                     (setf def (magicl:@ df def))
                     (setf strain (magicl:.+ strain dstrain))
                     (setf volume (* volume (magicl:det df)))
                     ))))

(defun update-strain-kirchoff (mp dstrain)
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   (domain cl-mpm/particle::mp-domain-size)
                   ) mp
    (progn
      (let ((df (.+ (magicl:eye 2) (voight-to-matrix dstrain))))
        (progn
          (setf def (magicl:@ df def))
          (multiple-value-bind (l v) (magicl:eig (voight-to-matrix strain))
            (let ((trial-lgs (magicl:@ df
                                       (magicl:@
                                        v
                                        (magicl:from-diag (mapcar (lambda (x) (exp (* x 2d0))) l) :type 'double-float)
                                        (magicl:transpose v))
                                       (magicl:transpose df))))
              (multiple-value-bind (lf vf) (magicl:eig trial-lgs)
                (setf strain (magicl:scale
                              (matrix-to-voight
                               (magicl:@ vf
                                         (magicl:from-diag (mapcar (lambda (x) (log x)) lf) :type 'double-float)
                                         (magicl:transpose vf)))
                              0.5d0)))))
          (setf volume (* volume (magicl:det df)))
          (multiple-value-bind (l v) (magicl:eig (magicl:@ df (magicl:transpose df)))
            (let ((stretch (magicl:@ v
                                     (magicl:from-diag (mapcar (lambda (x) (sqrt x)) l) :type 'double-float)
                                     (magicl:transpose v))))
              (setf (tref domain 0 0) (* (tref domain 0 0) (tref stretch 0 0)))
              (setf (tref domain 1 0) (* (tref domain 1 0) (tref stretch 1 1)))
              ;(setf domain (magicl:.* domain (magicl:diag stretch)))
              )
          ))))))

(defun update-stress-mp (mesh mp dt)
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt))
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (volume cl-mpm/particle:mp-volume)
                     (strain cl-mpm/particle:mp-strain)
                     (def    cl-mpm/particle:mp-deformation-gradient)
                     (strain-rate cl-mpm/particle:mp-strain-rate)
                     ;; (damage cl-mpm/particle:mp-damage)
                        ) mp
      (progn
        (calculate-strain-rate mesh mp dt)
        (let ((dstrain strain-rate))
          (progn
            ;; (setf strain-rate dstrain)
            ;;(setf strain (magicl:scale (magicl:.+ strain dstrain) (expt (- 1 damage) 2)))
            (progn
              ;; (update-strain-linear mp dstrain)
              (update-strain-kirchoff mp dstrain)
              (setf stress (cl-mpm/particle:constitutive-model mp strain dt))
              (when (<= volume 0d0)
                (error "Negative volume"))
              (setf stress (magicl:scale stress (/ 1.0 (magicl:det def))))
              (cl-mpm/particle:post-stress-step mesh mp dt)
              ))))))

(defun update-stress (mesh mps dt)
  (declare (array mps) (cl-mpm/mesh::mesh mesh))
  (lparallel:pdotimes (i (length mps))
     (update-stress-mp mesh (aref mps i) dt)))

(defun reset-grid (mesh)
  (declare (cl-mpm/mesh::mesh mesh))
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
  (dotimes (i (array-total-size nodes))
    (let ((node (row-major-aref nodes i)))
    (when (> (cl-mpm/mesh:node-mass node) 0)
      (cl-mpm/mesh:reset-node node))))))

(defun filter-grid (mesh mass-thresh)
  (declare (cl-mpm/mesh::mesh mesh) (double-float mass-thresh))
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
  (dotimes (i (array-total-size nodes))
    (let ((node (row-major-aref nodes i)))
      (when (and (> (cl-mpm/mesh:node-mass node) 0d0) 
                 (< (cl-mpm/mesh:node-mass node) mass-thresh))
        (cl-mpm/mesh:reset-node node))))))

(defun remove-material-damaged (sim)
  "Remove material points that have strain energy density exceeding fracture toughness"
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (delete-if (lambda (mp)
                           (with-accessors ((damage cl-mpm/particle:mp-damage))
                               mp
                             (>= damage 1d0))) mps)))
#||
(progn
(ql:quickload :cl-mpm/examples/fracture)
(in-package :cl-mpm/examples/fracture))
||#

