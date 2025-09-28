(defpackage :cl-mpm/dynamic-relaxation-mpi
  (:use :cl)
  (:export))
(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)

(defclass mpm-sim-quasi-static-mpi (cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
                                    cl-mpm/mpi::mpm-sim-mpi-nodes)
  ()
  (:default-initargs
   :vel-algo :QUASI-STATIC)
  (:documentation "DR psudo-linear step with update stress last update"))
