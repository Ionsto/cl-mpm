(defpackage :cl-mpm/errors
  (:use :cl)
  (:export
   #:error-simulation
   #:error-lagrangian
   #:error-dF-negative
   #:error-volume-negative
   ))
(in-package :cl-mpm/errors)


(define-condition error-simulation (error)
  ())

(define-condition error-lagrangian (error-simulation)
  ())

(define-condition error-dF-negative (error-simulation)
  ())

(define-condition error-volume-negative (error-simulation)
  ())

