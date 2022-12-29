(defpackage :cl-mpm/bc
  (:use :cl)
  (:export
    #:apply-bc
    #:bc-index
    #:bc-value
    #:make-outside-bc
    #:make-outside-bc-nostick
    #:make-outside-bc-roller
    #:make-bc-fixed
    #:make-bc-surface
    #:make-bc-friction
    #:make-bc-buoyancy
    #:make-bc-fixed-temp
    )
  )
(declaim (optimize (debug 0) (safety 0) (speed 3)))
(in-package :cl-mpm/bc)
(defclass bc ()
  ( (index 
      :accessor bc-index
      :initarg :index))
  (:documentation "A boundary condition that applies some operation at an index"))

(defclass bc-fixed (bc)
  ((value
     :accessor bc-value
     :initarg :value
     :initform '(nil nil)))
  (:documentation "A fixed velocity BC can be 1/0 dof"))

(defclass bc-fixed-temp (bc-fixed)
  ()
  (:documentation "A fixed velocity BC can be 1/0 dof"))

(defclass bc-ambient-temp (bc-fixed)
  ((threshold
    :accessor bc-threshold
    :initarg :threshold
    :initform 0d0
    :type DOUBLE-FLOAT)
   )
  (:documentation "BC that applies the temperature if theres no mass"))

(defun make-bc-fixed (index value)
  (make-instance 'bc-fixed
                 :index index 
                 :value value))

(defclass bc-surface (bc)
  ((normal
     :accessor bc-normal
     :initarg :normal
     :initform (magicl:zeros '(2 1))
     :type MAGICL:MATRIX/DOUBLE-FLOAT))
  (:documentation "A BC applied on a normal"))

(defclass bc-friction (bc)
  (
   (normal
    :accessor bc-normal
    :initarg :normal
    :initform (magicl:zeros '(2 1))
    :type MAGICL:MATRIX/DOUBLE-FLOAT)
   (friction-coefficent
    :accessor bc-friction-coefficent
    :initarg :friction-coefficent
    :initform 0.25d0
    :type DOUBLE-FLOAT)
   )
  (:documentation "A BC applied on a normal"))
(defclass bc-body-force (bc)
  ((force
    :accessor bc-force
    :initarg :force
    :initform (magicl:zeros '(2 1))
    :type MAGICL:MATRIX/DOUBLE-FLOAT))
  (:documentation "A BC applied like a body force"))

(defun make-bc-surface (index normal)
  (make-instance 'bc-surface
                 :index index 
                 :normal normal))

(defun make-bc-friction (index normal friction-coefficent)
  (make-instance 'bc-friction
                 :index index
                 :normal normal
                 :friction-coefficent friction-coefficent))
(defun make-bc-buoyancy (index force)
  (make-instance 'bc-body-force
                 :index index
                 :force force))

(defgeneric apply-bc (bc node mesh)
  (:documentation "Apply a boundary condition onto a node"))

(defmethod apply-bc (bc node mesh)
  "Natural BC condition is nothing")

(defmethod apply-bc ((bc bc-fixed) node mesh)
  "Fixed velocity BC over some dimensions"
  (with-slots ((value value))
    bc
    (loop for d from 0 to (length value)
            do (when (nth d value)
                 (setf (magicl:tref (cl-mpm/mesh:node-velocity node) d 0) (nth d value))))))


(defmethod apply-bc ((bc bc-fixed-temp) node mesh)
  "Fixed velocity BC over some dimensions"
  (with-slots ((value value))
    bc
    (when value
      (setf (cl-mpm/mesh::node-temperature node)
            value))))
(defmethod apply-bc ((bc bc-ambient-temp) node mesh)
  "Fixed velocity BC over some dimensions"
  (with-slots ((value value)
               (index index)
               (threshold threshold))
    bc
    (when (<= (cl-mpm/mesh::node-mass node) threshold)
      (loop for dx from -1 to 1
            do (loop for dy from -1 to 1
                     do
                        (let* ((ix (mapcar #'+ index (list dx dy))))
                          (when (cl-mpm/mesh::in-bounds mesh ix)
                            (setf (cl-mpm/mesh::node-temperature
                                   (cl-mpm/mesh:get-node mesh ix))
                                  value)
                            ;; (setf (cl-mpm/mesh::node-temperature
                            ;;        (cl-mpm/mesh::get-node mesh ix))
                            ;;       value)
                            )))))))

(defmethod apply-bc ((bc bc-surface) node mesh)
  "Fixed velocity BC over some non-stick surface"
  (with-slots ((normal normal))
    bc
    (with-accessors ((node-vel cl-mpm/mesh:node-velocity)) node
      (let ((rel-vel (magicl:tref (magicl:@ (magicl:transpose normal) node-vel) 0 0)))
        (when (< rel-vel 0)
          (setf node-vel (magicl:.- (magicl:scale normal rel-vel) node-vel)))))))

(defmethod apply-bc ((bc bc-friction) node mesh)
  "Frictional velocity BC over some non-stick surface"
  (with-slots ((normal normal)
               (mu friction-coefficent))
      bc
    (with-accessors ((node-vel cl-mpm/mesh:node-velocity)) node
      (let ((rel-vel (magicl::sum (magicl:.* node-vel normal))))
        (when (< rel-vel 0)
          (let* ((tang-vel (magicl:.- node-vel (magicl:scale normal rel-vel)))
                 (friction-impulse (magicl:scale tang-vel mu))
                 )
            (setf node-vel (magicl:.- node-vel friction-impulse))
            ))
        (setf node-vel (magicl:.- node-vel (magicl:scale normal rel-vel)))
        ))))


(defun make-bc-fixed-temp (index value)
  (make-instance 'bc-fixed-temp
                 :index index 
                 :value value))

(defun make-bc-ambient-temp (index value threshold)
  (make-instance 'bc-ambient-temp
                 :index index 
                 :threshold threshold
                 :value value))

(defmethod apply-bc ((bc bc-body-force) node mesh)
  "Fixed velocity BC over some dimensions"
  (with-slots ((value value))
      bc
    (setf (cl-mpm/mesh:node-force node) (magicl:.+ (cl-mpm/mesh:node-force node)
                                                 (magicl:scale (bc-force bc) (cl-mpm/mesh::node-volume node))
                                                 ))))
(defclass bc-volume (bc)
  ((index
     :accessor bc-index
     :initarg :index))
  (:documentation "A boundary condition that applies some operation over a set of nodes"))
(defclass bc-inflow (bc-fixed)
  ((vel-bcs
    :accessor bc-inflow-vel
    :initarg :vel-bcs)
   (volume-threshold
    :accessor bc-inflow-volume-threshold
    :initarg :volume-threshold)
   (particle-constructor
    :accessor bc-inflow-particle-constructor
    :initarg :particle-constructor
    ))
  (:documentation "A boundary condition that applies some operation over a set of nodes"))

(defun make-bc-inflow (index vel volume-threshold p-cons)
    (make-instance 'bc-inflow
                   :index index
                   :particle-constructor p-cons
                   :volume-threshold volume-threshold
                   :value vel)
  )
(defmethod apply-bc ((bc bc-inflow) node mesh)
  "Nodal inflow"
  (with-slots ((value value)
               (volume-threshold volume-threshold)
               (pcons particle-constructor)
               )
      bc
    (when (<= (cl-mpm/mesh::node-volume node) volume-threshold)
      ;; Add particle?
      (funcall pcons (list (magicl:tref (cl-mpm/mesh::node-position node) 0 0)
                           (magicl:tref (cl-mpm/mesh::node-position node) 1 0)
                           ))
      )
    )
  (call-next-method))

;; (defun make-wall-bc (mesh &rest args)
;;   (loop for v from 0 to ()))
(defun make-outside-bc (mesh-count)
  "Construct fixed bcs over the outside of a mesh"
  (destructuring-bind (xsize ysize) (mapcar (lambda (x) (- x 1)) mesh-count)
    (append 
      (loop for x from 0 to xsize 
            append 
            (list (make-bc-fixed (list x 0)     '(0d0 0d0))
                  (make-bc-fixed (list x ysize) '(0d0 0d0))))
       (loop for y from 0 to ysize 
            append 
            (list (make-bc-fixed (list 0     y) '(0d0 0d0))
                  (make-bc-fixed (list xsize y) '(0d0 0d0)))))))

(defun make-outside-bc-nostick (mesh-count)
    "Construct nostick bcs over the outside of a mesh"
    (destructuring-bind (xsize ysize) (mapcar (lambda (x) (- x 1)) mesh-count)
        (append 
            (loop for x from 0 to xsize 
                append 
                (list (make-bc-surface (list x 0)     (magicl:from-list '(0d0  1d0) '(2 1)))
                      (make-bc-surface (list x ysize) (magicl:from-list '(0d0 -1d0) '(2 1)))))
            (loop for y from 0 to ysize 
                append 
                (list (make-bc-surface (list 0     y) (magicl:from-list '( 1d0 0d0) '(2 1)))
                      (make-bc-surface (list xsize y) (magicl:from-list '(-1d0 0d0) '(2 1))))))))

(defun make-outside-bc-roller (mesh-count order)
  "Construct fixed bcs over the outside of a mesh"
  (destructuring-bind (xsize ysize) (mapcar (lambda (x) (- x 1)) mesh-count)
    (loop for o from 0 to order
          append
          (append
           (loop for x from 0 to xsize 
                 append 
                 (list (make-bc-fixed (list x order)     '(nil 0d0))
                       (make-bc-fixed (list x (- ysize order)) '(nil 0d0))))
           (loop for y from 0 to ysize 
                 append 
                 (list (make-bc-fixed (list order     y) '(0d0 nil))
                       (make-bc-fixed (list (- xsize order) y) '(0d0 nil))))))))

(defun make-outside-bc-var (mesh left right top bottom)
  "Construct fixed bcs over the outside of a mesh"
  (with-accessors ((mesh-count cl-mpm/mesh:mesh-count)
                   (order cl-mpm/mesh::mesh-boundary-order))
      mesh
    (destructuring-bind (xsize ysize) (mapcar (lambda (x) (- x 1)) mesh-count)
      (loop for o from 0 to order
            append
            (append
             (loop for x from 0 to xsize 
                   append 
                   (list (funcall bottom (list x order))
                         (funcall top (list x (- ysize order)))))
             (loop for y from 0 to ysize 
                   append 
                   (list (funcall left (list order y))
                         (funcall right (list (- xsize order) y)))))))))
(defun make-sub-domain-bcs (mesh start end make-bc)
  "Construct  bcs over the outside of a mesh"
  (with-accessors ((mesh-count cl-mpm/mesh:mesh-count)
                   (order cl-mpm/mesh::mesh-boundary-order))
      mesh
    (destructuring-bind (xstart ystart) start
      (destructuring-bind (xend yend) end
        (loop for x from
              xstart to xend
              append
              (loop for y from
                    ystart to yend
                    collect
                    (funall make-bc (list x y))))))))
(defun make-domain-bcs (mesh make-bc)
  "Construct  bcs over the outside of a mesh"
  (with-accessors ((mesh-count cl-mpm/mesh:mesh-count)
                   (order cl-mpm/mesh::mesh-boundary-order))
      mesh
    (destructuring-bind (xsize ysize) (mapcar (lambda (x) (- x 1)) mesh-count)
      (loop for x from
            0 to xsize
            append
            (loop for y from
                  0 to ysize
                  collect
                       (funcall make-bc (list x y)))))))
