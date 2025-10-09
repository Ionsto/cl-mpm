(defpackage :cl-mpm/errors
  (:use :cl)
  (:export
   #:error-lagrangian))
(in-package :cl-mpm/errors)
(define-condition error-lagrangian (error)
  ())

