(defpackage :cl-mpm/settings
  (:use :cl)
  (:export #:*optimise-setting*)
  )
(in-package :cl-mpm/settings)
(defvar *optimise-debug*
  '(optimize
    (speed 0)
    (safety 3)
    (debug 3)
    (compilation-speed 0)))
(defvar *optimise-speed*
  '(optimize
    (speed 3)
    (safety 0)
    (debug 0)
    (compilation-speed 0)))

;; (defvar *optimise-setting* *optimise-debug*)
(defvar *optimise-setting* *optimise-speed*)
