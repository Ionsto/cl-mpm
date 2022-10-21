(defpackage :cl-mpm/damage
  (:use :cl
        :cl-mpm/utils)
  (:export
    #:update-damage
    #:post-stress-step
    )
  )
(in-package :cl-mpm/damage)
(declaim (optimize (debug 0) (safety 0) (speed 3)))
(defun damage-rate-profile (stress damage)
  (let ((critical-stress 1e3))
    (if (> stress (* 1/2 critical-stress))
        (* (/ (max 0 stress) critical-stress) 0.1)
        0)))
(defun damage-profile (damage)
  (expt (- 1 damage) 2)
  )
(defun update-damage (mp dt)
  (let ((damage-increment 0))
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (damage cl-mpm/particle:mp-damage)) mp
        (progn
          (multiple-value-bind (l v) (magicl:eig (voight-to-matrix stress))
            (loop for i from 0 to 1
                  do (let ((sii (nth i l)))
                       (progn
                         (incf damage-increment (* (damage-rate-profile sii damage) dt)))))
            (incf damage damage-increment)
            (setf damage (min 1 (max 0 damage)))
            (setf stress (magicl:scale stress (damage-profile damage)))
            )))))
(defmethod cl-mpm/particle:post-stress-step (mesh (mp cl-mpm/particle:particle-damage) dt)
  (update-damage mp dt))
