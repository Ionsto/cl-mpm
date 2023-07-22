(defpackage :cl-mpm/bc
  (:use :cl)
  (:export
    #:apply-bc
    #:bc-index
    #:bc-value
    #:make-outside-bc
    #:make-outside-bc-nostick
    #:make-outside-bc-roller
    #:make-bcs-from-list
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
      :initarg :index)
   (node
    :accessor bc-node
    :initarg :node
    :initform nil)
   )
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

(defclass bc-closure (bc)
  ((func
    :accessor bc-func
    :initarg :func))
  (:documentation "A bc with a function closure"))

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

(defgeneric apply-bc (bc node mesh dt)
  (:documentation "Apply a boundary condition onto a node"))

(defmethod apply-bc (bc node mesh dt)
  "Natural BC condition is nothing")

(defmethod apply-bc ((bc bc-fixed) node mesh dt)
  "Fixed velocity BC over some dimensions"
  (with-slots ((value value))
    bc
    (loop for d from 0 to (length value)
            do (when (nth d value)
                 (setf (magicl:tref (cl-mpm/mesh:node-velocity node) d 0) (nth d value))
                 (setf (magicl:tref (cl-mpm/mesh:node-acceleration node) d 0) (nth d value))
                 ))))


(defmethod apply-bc ((bc bc-fixed-temp) node mesh dt)
  "Fixed velocity BC over some dimensions"
  (with-slots ((value value))
    bc
    (when value
      (setf (cl-mpm/mesh::node-temperature node)
            value))))
(defmethod apply-bc ((bc bc-ambient-temp) node mesh dt)
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
                            (let ((n (cl-mpm/mesh:get-node mesh ix)))
                              (with-accessors ((mass cl-mpm/mesh::node-mass)
                                               (temp cl-mpm/mesh::node-temperature)
                                               (dtemp cl-mpm/mesh::node-dtemp)
                                               )
                                  n
                                ;(setf temp (* value mass))
                                ;(setf dtemp 0)
                                (when (> mass 0d0)
                                  (incf dtemp
                                        (* 5d1
                                           (- value (/ temp mass)))
                                        ))
                                  )
                              )
                            ;; (setf (cl-mpm/mesh::node-temperature
                            ;;        (cl-mpm/mesh::get-node mesh ix))
                            ;;       value)
                            )))))))

(defmethod apply-bc ((bc bc-surface) node mesh dt)
  "Fixed velocity BC over some non-stick surface"
  (with-slots ((normal normal))
    bc
    (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
                     (node-acc cl-mpm/mesh:node-acceleration)
                     (node-force cl-mpm/mesh:node-force)
                     (node-lock cl-mpm/mesh:node-lock)
                     (node-pressure cl-mpm/mesh::node-pressure)
                     ) node
      (let ((rel-vel (cl-mpm/fastmath::mag (magicl:.* normal node-vel)))
            (rel-force (cl-mpm/fastmath::mag (magicl:.* normal node-force)))
            (rel-acc (cl-mpm/fastmath::mag (magicl:.* normal node-acc)))
            )
        ;; (when (< rel-vel 0)
          ;; (magicl:.+ node-vel (magicl:scale normal (sqrt (abs rel-vel))) node-vel)
          ;; (magicl:.+ node-force (magicl:scale normal (sqrt (abs rel-force))) node-force)
          ;; )
        ;; (when (< rel-force 0)
        ;;   (magicl:.+ node-force (magicl:scale normal (sqrt (abs rel-force))) node-force))
          ;; (magicl:.+ node-acc (magicl:scale normal (sqrt (abs rel-acc))) node-acc)
        (when (> rel-force 0)
          ;; (format t "~f ~%" (+ (sqrt rel-force) node-pressure))
          (when (> (+ (sqrt rel-force) node-pressure) 0)
                                        ;(magicl:.+ node-vel (magicl:scale normal (sqrt rel-vel)) node-vel)
            (setf (magicl:tref node-vel 1 0) 0d0)
            ;; (setf (magicl:tref node-acc 1 0) (* -1d0 (magicl:tref node-acc 1 0)))
            ;; (setf (magicl:tref node-force 1 0) (* -1d0 (magicl:tref node-force 1 0)))
            ))
        ;; (when (> rel-acc 0)
        ;;   (setf (magicl:tref node-acc 1 0) 0d0)
        ;;   )
        ;; (when (> rel-force 0)
        ;;   (setf (magicl:tref node-force 1 0) 0d0)
        ;;   )
        ))))

(defmethod apply-bc ((bc bc-friction) node mesh dt)
  "Frictional velocity BC over some non-stick surface"
  (with-slots ((normal normal)
               (mu friction-coefficent))
      bc
    (with-accessors (
                     (node-vel cl-mpm/mesh:node-velocity)
                     (node-force cl-mpm/mesh:node-force)
                     ) node
      (let ((rel-vel (magicl::sum (magicl:.* node-vel normal))))
        (when t;(< rel-vel 0)
          ;Note this is momentum
          (let* ((tang-p (magicl:.- node-vel (magicl:scale normal rel-vel)))
                 (h (cl-mpm/mesh:mesh-resolution mesh))
                 ;; (normal-force (magicl::sum (magicl:.* node-force normal)))
                 (friction-force (magicl:scale tang-p (* mu h)))
                 ;(friction-force (magicl:scale tang-p (* mu (cl-mpm/mesh:mesh-resolution mesh))))
                 ;; (friction-force (magicl:scale (cl-mpm/fastmath::norm tang-p) (* mu h)))
                 )
            ;(setf node-vel (magicl:.- node-vel friction-impulse))
            ;; (setf node-force (magicl:.- node-force friction-force))
            (setf node-force (magicl:.- node-force friction-force))
            ;; (setf node-vel (magicl:.+ node-vel (magicl:scale tang-p -0.5d0)))
            ))
        ;(setf node-vel (magicl:.- node-vel (magicl:scale normal rel-vel)))
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

(defmethod apply-bc ((bc bc-body-force) node mesh dt)
  "Fixed velocity BC over some dimensions"
  (with-slots ((value value))
      bc
    (when (cl-mpm/mesh:node-active node)
      (cl-mpm/fastmath:fast-fmacc (cl-mpm/mesh:node-force node) (bc-force bc) (cl-mpm/mesh::node-volume node)))
    ;; (setf (cl-mpm/mesh:node-force node)
    ;;       (magicl:.+ (cl-mpm/mesh:node-force node)
    ;;                  (magicl:scale (bc-force bc) (cl-mpm/mesh::node-volume node))))
    ))
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
(defmethod apply-bc ((bc bc-inflow) node mesh dt)
  "Nodal inflow"
  (with-slots ((value value)
               (volume-threshold volume-threshold)
               (pcons particle-constructor)
               )
      bc
    (when (<= (cl-mpm/mesh::node-volume node) volume-threshold)
      (print "Created particle")
      ;; Add particle?
      (funcall pcons (list (magicl:tref (cl-mpm/mesh::node-position node) 0 0)
                           (magicl:tref (cl-mpm/mesh::node-position node) 1 0)
                           ))
      )
    )
  (call-next-method))

;; (defun make-wall-bc (mesh &rest args)
;;   (loop for v from 0 to ()))

(defun make-bcs-from-list (bc-list)
  (make-array (length bc-list) :initial-contents bc-list :adjustable t :fill-pointer (length bc-list)))

(defun make-outside-bc (mesh-count)
  "Construct fixed bcs over the outside of a mesh"
  (destructuring-bind (xsize ysize) (mapcar (lambda (x) (- x 1)) mesh-count)
    (make-bcs-from-list
     (append
      (loop for x from 0 to xsize
            append
            (list (make-bc-fixed (list x 0)     '(0d0 0d0))
                  (make-bc-fixed (list x ysize) '(0d0 0d0))))
      (loop for y from 0 to ysize
            append
            (list (make-bc-fixed (list 0     y) '(0d0 0d0))
                  (make-bc-fixed (list xsize y) '(0d0 0d0))))))))

(defun make-outside-bc-nostick (mesh-count)
    "Construct nostick bcs over the outside of a mesh"
    (destructuring-bind (xsize ysize) (mapcar (lambda (x) (- x 1)) mesh-count)
      (make-bcs-from-list
       (append 
        (loop for x from 0 to xsize 
              append 
              (list (make-bc-surface (list x 0)     (magicl:from-list '(0d0  1d0) '(2 1)))
                    (make-bc-surface (list x ysize) (magicl:from-list '(0d0 -1d0) '(2 1)))))
        (loop for y from 0 to ysize 
              append 
              (list (make-bc-surface (list 0     y) (magicl:from-list '( 1d0 0d0) '(2 1)))
                    (make-bc-surface (list xsize y) (magicl:from-list '(-1d0 0d0) '(2 1)))))))))

(defun make-outside-bc-roller (mesh-count order)
  "Construct fixed bcs over the outside of a mesh"
  (destructuring-bind (xsize ysize) (mapcar (lambda (x) (- x 1)) mesh-count)
    (make-bcs-from-list
     (remove nil
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
                                (make-bc-fixed (list (- xsize order) y) '(0d0 nil))))))))))

(defun make-outside-bc-var-list (mesh left right top bottom)
  "Construct fixed bcs over the outside of a mesh"
  (with-accessors ((mesh-count cl-mpm/mesh:mesh-count)
                   (order cl-mpm/mesh::mesh-boundary-order))
      mesh
    (destructuring-bind (xsize ysize) (mapcar (lambda (x) (- x 1)) mesh-count)
      (remove nil
              (loop for o from 0 to order
                    append
                    (append
                     (loop for x from 0 to xsize
                           append
                           (list (when bottom (funcall bottom (list x o)))
                                 (when top (funcall top (list x (- ysize o))))))
                     (loop for y from 0 to ysize
                           append
                           (list (when left (funcall left (list o y)))
                                 (when right (funcall right (list (- xsize o) y)))))))))))
(defun make-outside-bc-var (mesh left right top bottom)
  "Construct fixed bcs over the outside of a mesh"
      (make-bcs-from-list
       (make-outside-bc-var-list mesh left right top bottom)))
(defun make-sub-domain-bcs (mesh start end make-bc)
  "Construct  bcs over the outside of a mesh"
  (with-accessors ((mesh-count cl-mpm/mesh:mesh-count)
                   (order cl-mpm/mesh::mesh-boundary-order))
      mesh
    (destructuring-bind (xstart ystart) start
      (destructuring-bind (xend yend) end
        (make-bcs-from-list
         (loop for x from
               xstart to xend
               append
               (loop for y from
                     ystart to yend
                     collect
                     (funall make-bc (list x y)))))))))
(defun make-domain-bcs (mesh make-bc)
  "Construct  bcs over the outside of a mesh"
  (with-accessors ((mesh-count cl-mpm/mesh:mesh-count)
                   (order cl-mpm/mesh::mesh-boundary-order))
      mesh
    (destructuring-bind (xsize ysize) (mapcar (lambda (x) (- x 1)) mesh-count)
      (make-bcs-from-list
       (loop for x from
             0 to xsize
             append
             (loop for y from
                   0 to ysize
                   collect
                   (funcall make-bc (list x y))))))))

(defun make-bc-closure (index func)
  (make-instance 'bc-closure
                 :index index
                 :func func))

(defmethod apply-bc ((bc bc-closure) node mesh dt)
  "Arbitrary closure BC"
  (with-slots ((func func))
    bc
    (declare (function func))
    (funcall func)))

;;; Printing methods

(defmethod print-object ((obj bc) stream)
  (print-unreadable-object (obj stream :type t)
    (with-accessors ((index bc-index))
        obj
      (format stream "index: ~a" index))))
