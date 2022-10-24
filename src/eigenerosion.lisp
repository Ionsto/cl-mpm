(defpackage :cl-mpm/eigenerosion
  (:use :cl
        :cl-mpm
        :cl-mpm/utils)
  (:export
    #:update-fracture
    )
  )
(in-package :cl-mpm/eigenerosion)
(declaim (optimize (debug 0) (safety 0) (speed 3)))
(defun update-fracture (sim)
  (with-accessors ((mps sim-mps))
      sim
      (loop for mp across mps
            do
               (when (> (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0) 1e4)
                 (progn
                   (setf (cl-mpm/particle:mp-damage mp) 1)
                   )))))
