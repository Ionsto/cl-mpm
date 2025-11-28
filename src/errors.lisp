(defpackage :cl-mpm/errors
  (:use :cl)
  (:export
   #:error-simulation
   #:error-lagrangian))
(in-package :cl-mpm/errors)


(define-condition error-simulation (error)
  ())

(define-condition error-lagrangian (error-simulation)
  ())


