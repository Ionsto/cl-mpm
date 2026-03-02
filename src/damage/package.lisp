(defpackage :cl-mpm/damage
  (:use :cl
   :cl-mpm/utils)
  (:export
   #:mpm-sim-damage
   #:calculate-damage
   #:update-damage
   #:post-stress-step
   #:update-damage
   #:damage-model-calculate-y))

(in-package :cl-mpm/damage)

(defclass mpm-sim-damage (cl-mpm::mpm-sim-usf)
  ((delocal-counter
    :accessor sim-damage-delocal-counter
    :type fixnum
    :initarg :delocal-counter
    :initform 0)
   (delocal-counter-max
    :accessor sim-damage-delocal-counter-max
    :type fixnum
    :initarg :delocal-counter-max
    :initform 10)
   (enable-length-localisation
    :type boolean
    :accessor sim-enable-length-localisation
    :initarg :enable-length-localisation
    :initform nil)
   (enable-ekl
    :type boolean
    :accessor sim-enable-ekl
    :initarg :enable-ekl
    :initform nil)
   )
  (:documentation "Explicit simulation with update stress first update"))

(defclass mpm-sim-usl-damage (mpm-sim-damage cl-mpm::mpm-sim-usl) ())

(defclass mpm-sim-musl-damage (mpm-sim-damage cl-mpm::mpm-sim-musl) ())

(defclass mpm-sim-damage-nd-2 (mpm-sim-damage cl-mpm::mpm-nd-2d) ())

(defclass mpm-sim-agg-damage (mpm-sim-damage cl-mpm/aggregate::mpm-sim-aggregated)
  ()
  (:documentation "Explicit damage simulation"))
