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
(defun damage-rate-profile (critical-stress stress damage)
  "Function that controls how damage evolves with principal stresses"
  (if (> stress (* 1/2 critical-stress))
      (* (/ (max 0 stress) critical-stress) 0.1)
      0)
  0
  )

(defun damage-profile (damage)
  "Constitive law describing the scalar stress decrease as a function of damage"
  (expt (- 1 damage) 2)
  )
(defun update-damage (mp dt)
  (let ((damage-increment 0))
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (critical-stress cl-mpm/particle:mp-critical-stress)
                     ) mp
        (progn
          (multiple-value-bind (l v) (magicl:eig (voight-to-matrix stress))
            (loop for i from 0 to 1
                  do (let ((sii (nth i l)))
                       (progn
                         (incf damage-increment (* (damage-rate-profile critical-stress sii damage) dt)))))
            (incf damage damage-increment)
            (setf damage (min 1 (max 0 damage)))
            (loop for i from 0 to 1
                  do (let ((sii (nth i l)))
                       (progn
                         (when (> sii 0)
                           (setf (nth i l) (* (nth i l) (damage-profile damage)))))))

            ;; (loop for sii in l
            ;;       do (when (> sii 0)
            ;;            (setf sii (* sii (damage-profile damage)))))

            ;; (setf stress (matrix-to-voight (magicl:@ v
            ;;                                          (magicl:from-list (list (first l) 0 0 (second l)) '(2 2) :type 'double-float)
            ;;                                          (magicl:inv v))))
            (when (> damage 0)
              (setf stress (matrix-to-voight (magicl:@ v
                                                       (magicl:from-diag l :type 'double-float)
                                                       (magicl:inv v)))))
            ;; (setf stress (magicl:scale stress (damage-profile damage)))
            )))))
(defmethod cl-mpm/particle:post-stress-step (mesh (mp cl-mpm/particle:particle-damage) dt)
  (update-damage mp dt)
  )
