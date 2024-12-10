(defpackage :cl-mpm
  (:use :cl
   :cl-mpm/particle
   :cl-mpm/mesh
   :cl-mpm/utils
   :cl-mpm/fastmaths)
  (:import-from
   :magicl tref .+ .-
   )
  (:import-from
   :cl-mpm/forces
   det-int-force
   det-ext-force
   det-ext-force-2d
   det-int-force-unrolled
   det-int-force-unrolled-2d
   )
  (:export
   #:mpm-sim
   #:make-mpm-sim
   #:update-sim
   #:make-shape-function-linear
   #:sim-mesh
   #:sim-mps
   #:sim-bcs
   #:sim-bcs-force-list
   #:sim-dt
   #:sim-damping-factor
   #:sim-mass-filter
   #:sim-mass-scale
   #:sim-allow-mp-split
   #:post-stress-step
   #:iterate-over-nodes
   #:iterate-over-nodes-serial
   #:iterate-over-neighbours
   #:iterate-over-mps
   #:iterate-over-corners
   #:iterate-over-neighbours-point-linear
   #:calculate-adaptive-time
   #:add-mps
   #:add-bcs
   #:add-bcs-force-list))
;; System definition for cl-mpm
(in-package :cl-mpm)
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(defclass mpm-sim ()
  ((dt
     :accessor sim-dt
     :initarg :dt)
   (time
    :accessor sim-time
    :initform 0d0
    :initarg :time)
   (stats-mps-removed
    :accessor sim-stats-mps-removed
    :initform 0d0
    :initarg :stats-mps-removed)
   (stats-mps-added
    :accessor sim-stats-mps-added
    :initform 0d0
    :initarg :stats-mps-added)
   (mesh
     :accessor sim-mesh
     :initarg :mesh)
   (mps
     :accessor sim-mps
     :initarg :mps)
   (bcs
     :accessor sim-bcs
     :initarg :bcs
     :initform (make-array 0))
   (bcs-force
    :accessor sim-bcs-force
    :initarg :bcs-force
    :initform (make-array 0))
   (bcs-force-list
    :accessor sim-bcs-force-list
    :initarg :bcs-force-list
    :initform nil)
   (damping-factor
     :type double-float
     :accessor sim-damping-factor
     :initarg :damping-factor
     :initform 0d0)
   (mass-scale
    :type double-float
    :accessor sim-mass-scale
    :initarg :mass-scale
    :initform 1d0)
   (allow-mp-split
    :type boolean
    :accessor sim-allow-mp-split
    :initarg :splitting
    :initform nil)
   (enable-damage
    :type boolean
    :accessor sim-enable-damage
    :initarg :enable-damage
    :initform nil)
   (enable-fbar
    :type boolean
    :accessor sim-enable-fbar
    :initarg :enable-fbar
    :initform nil)
   (nonlocal-damage
    :type boolean
    :accessor sim-nonlocal-damage
    :initform t)
   (allow-mp-damage-removal
    :type boolean
    :accessor sim-allow-mp-damage-removal
    :initarg :damage-removal
    :initform nil)
   (mp-damage-removal-instant
    :type boolean
    :accessor sim-mp-damage-removal-instant
    :initform nil)
   (velocity-algorithm
    :type symbol
    :accessor sim-velocity-algorithm
    :initform :BLEND)
   (mass-filter
     :type double-float
     :accessor sim-mass-filter
     :initarg :mass-filter
     :initform 1d-15))
  (:documentation "A self contained mpm simulation"))
(defclass mpm-nd-2d ()())
(defclass mpm-nd-3d ()())
(defclass mpm-sf-mpm ()())
(defclass mpm-sf-gimp ()())
(defclass mpm-sim-quasi-static (mpm-sim) ())

(defclass mpm-sim-usf (mpm-sim)
  ()
  (:documentation "Explicit simulation with update stress first update"))
(defclass mpm-sim-usl (mpm-sim)
  ()
  (:documentation "Explicit simulation with update stress last update"))


(defun make-mpm-sim (size resolution dt shape-function &key (sim-type 'mpm-sim-usf))
  "Constructs an mp with critical infomation like mesh and number of dimentions"
  (make-instance sim-type
                 :dt (coerce dt 'double-float)
                 :mesh (make-mesh size resolution shape-function)
                 :mps '()))
