(defpackage :cl-mpm/dynamic-relaxation
  (:use :cl)
  (:export
   #:estimate-energy-norm
   #:estimate-power-norm
   #:estimate-oobf
   #:converge-quasi-static))
(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)

(defclass mpm-sim-dr (cl-mpm/aggregate::mpm-sim-aggregated)
  ((dt-loadstep
    :initform 1d0
    :initarg :dt-loadstep
    :accessor sim-dt-loadstep)
   (damping-scale
    :initform 1d0
    :initarg :damping-scale
    :accessor sim-damping-scale)
   (kinetic-damping
    :initform nil
    :accessor sim-kinetic-damping)
   (ke
    :initform 0d0
    :accessor sim-ke)
   (ke-prev
    :initform 0d0
    :accessor sim-ke-prev)
   )
  (:default-initargs
   :vel-algo :QUASI-STATIC)
  (:documentation "DR psudo-linear step with update stress last update"))

(defclass mpm-sim-dr-ul (mpm-sim-dr cl-mpm/aggregate::mpm-sim-aggregated)
  ((initial-setup
    :initform nil
    :accessor sim-initial-setup))
  (:documentation "DR solved quasi-static implicit mpm"))

(defclass mpm-sim-dr-dynamic (mpm-sim-dr-ul cl-mpm/damage:mpm-sim-damage)
  ((enable-dynamics
    :initform t
    :accessor sim-enable-dynamics))
  (:documentation "DR implicit-dynamic mpm base class"))


(defclass mpm-sim-quasi-static (mpm-sim-dr-ul)
  ()
  (:documentation "DR implicit quasi-static wrapper class"))

(defclass mpm-sim-implict-dynamic (mpm-sim-dr-dynamic)
  ()
  (:documentation "DR implicit-dynamic wrapper class - acts just like an explicit solver"))

(defclass mpm-sim-dr-usf (mpm-sim-dr)
  ()
  (:documentation "DR explicit-dynamics with damage"))

(defclass mpm-sim-dr-damage-usf (mpm-sim-dr cl-mpm/damage:mpm-sim-damage)
  ()
  (:documentation "DR explicit-dynamics with damage"))

(defclass mpm-sim-dr-damage-ul (mpm-sim-dr-ul cl-mpm/damage:mpm-sim-damage)
  ()
  (:documentation "DR implicit damage"))

(defclass mpm-sim-dr-multigrid (mpm-sim-dr-damage-ul cl-mpm::mpm-sim-multigrid)
  ()
  (:documentation "Dynamic relaxation with multigrid"))

;; (defclass mpm-sim-dr-multigrid (mpm-sim-dr-ul cl-mpm::mpm-sim-multigrid)
;;   ()
;;   (:documentation "Dynamic relaxation with multigrid"))
