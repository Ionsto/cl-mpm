(defpackage :cl-mpm/example
  (:use :cl)
  (:import-from
   :cl-mpm/utils varef)
  (:export
   #:plot-domain
   #:stop
   #:run
   #:setup
   #:*sim*
   #:*run-sim*
   )
  )
(in-package :cl-mpm/example)
(defparameter *run-sim* nil)
(defparameter *sim* nil)
(defun plot-domain ()
  (when *sim*
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     ;; :colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :xy))
     )))
(defun stop ()
  (setf *run-sim* nil))
(defun run ())
(defun setup (&key (refine 1) (mps 2))
  )
