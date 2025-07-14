(defpackage :cl-mpm/ssa
  (:use :cl :cl-mpm/particle :cl-mpm/mesh :cl-mpm/utils :cl-mpm/fastmaths))

(in-package :cl-mpm/ssa)
(declaim #.cl-mpm/settings:*optimise-setting*)

(defclass mpm-sim-ssa (cl-mpm::mpm-sim)
  ())
