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
  (if (> stress (* 0.0d0 critical-stress))
      (/ (* (/ (max 0d0 stress) critical-stress) 1d-2) (max 1d-5 (expt (- 1d0 damage) 3)))
      ;; (* (/ (max 0d0 stress) critical-stress) 1d-2)
      0d0)
  )

(defun damage-profile (damage)
  "Constitive law describing the scalar stress decrease as a function of damage"
  (- 1 damage)
  ;(max (expt (- 1 damage) 2) 1d-6)
  )
(defun calculate-damage-increment (mp dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (critical-stress cl-mpm/particle:mp-critical-stress)
                     ) mp
        (progn
          (multiple-value-bind (l v) (magicl:eig (voight-to-matrix stress))
            (let ((s1 (loop for i from 0 to 1
                            maximize (nth i l))))
              (when (> s1 0d0)
                (incf damage-increment (* (damage-rate-profile critical-stress s1 damage) dt))))
            (setf damage-increment (max 0d0 (min damage-increment (- 1d0 damage))))
            (setf (cl-mpm/particle::mp-local-damage-increment mp)
                  damage-increment))))))
(defun apply-damage (mp dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (damage-inc cl-mpm/particle::mp-local-damage-increment)
                     ) mp
        (progn
          (incf damage damage-inc)
          (setf damage (max 0d0 (min 1d0 damage)))
          (setf undamaged-stress (magicl:scale stress 1d0))

          (when (> damage 0.0d0)
            (multiple-value-bind (l v) (magicl:eig (voight-to-matrix stress))
              (loop for i from 0 to 1
                    do (let ((sii (nth i l)))
                           (when (> sii 0d0)
                             (setf (nth i l)
                                   (* (nth i l) (damage-profile damage))
                                   ))))
              (setf stress (matrix-to-voight (magicl:@ v
                                                       (magicl:from-diag l :type 'double-float)
                                                       (magicl:transpose v))))
              )
            ))))
(defun update-damage (mp dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (damage-inc cl-mpm/particle:mp-damage)
                     (critical-stress cl-mpm/particle:mp-critical-stress)
                     ) mp
        (progn
          (setf undamaged-stress stress)
          (multiple-value-bind (l v) (magicl:eig (voight-to-matrix undamaged-stress))
            ;; (loop for i from 0 to 1
            ;;       do (let ((sii (nth i l)))
            ;;            (progn
            ;;              (incf damage-increment (* (damage-rate-profile critical-stress sii damage) dt)))))
            (let ((s1 (loop for i from 0 to 1
                            maximize (nth i l))))
              (when (> s1 0d0)
                (incf damage-increment (* (damage-rate-profile critical-stress s1 damage) dt))))
            (setf damage-increment (min damage-increment (- 1d0 damage)))
            (setf damage-inc damage-increment)
            (incf damage damage-increment)

            (setf damage (min 1d0 (max 0d0 damage)))

            ;; (loop for sii in l
            ;;       do (when (> sii 0)
            ;;            (setf sii (* sii (damage-profile damage)))))

            ;; (setf stress (matrix-to-voight (magicl:@ v
            ;;                                          (magicl:from-list (list (first l) 0 0 (second l)) '(2 2) :type 'double-float)
            ;;                                          (magicl:inv v))))
            (when (> damage 0d0)
              ;(setf stress (magicl:scale stress (- 1d0 damage)))
              (loop for i from 0 to 1
                    do (let ((sii (nth i l)))
                         (progn
                           (when (> sii 0)
                             (setf (nth i l) (* (nth i l) (damage-profile damage)))))))
              (setf stress (matrix-to-voight (magicl:@ v
                                                       (magicl:from-diag l :type 'double-float)
                                                       (magicl:inv v))))
              )
            )))))
(defmethod cl-mpm/particle:post-stress-step (mesh (mp cl-mpm/particle:particle-damage) dt)
  ;(update-damage mp dt)
  )

(defun calculate-damage (mesh mps dt len)
    (create-delocalisation-list mesh mps len)
    (lparallel:pdotimes (i (length mps))
      (calculate-damage-increment (aref mps i) dt))
    (delocalise-damage mesh mps dt len)
    (lparallel:pdotimes (i (length mps))
      (apply-damage (aref mps i) dt))
  )
(defun create-delocalisation-list (mesh mps length)
  (with-accessors ((nodes cl-mpm/mesh:mesh-nodes))
        mesh
      (lparallel:pdotimes (i (array-total-size nodes))
        (let ((node (row-major-aref nodes i)))
          (setf (fill-pointer (cl-mpm/mesh::node-local-list node)) 0)))
    (lparallel:pdotimes (i (length mps))
      (let ((mp (aref mps i)))
            (let ((node-id (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp))))
              (when (cl-mpm/mesh:in-bounds mesh node-id)
                (let ((node (cl-mpm/mesh:get-node mesh node-id)))
                (sb-thread:with-mutex ((cl-mpm/mesh:node-lock node))
                  (vector-push-extend mp (cl-mpm/mesh::node-local-list node))))))))
    (lparallel:pdotimes (i (array-total-size nodes))
      (let ((node (row-major-aref nodes i)))
        (when (= 0 (length (cl-mpm/mesh::node-local-list node)))
          (adjust-array (cl-mpm/mesh::node-local-list node) 0 :fill-pointer 0))))))

(declaim
 (inline diff-squared)
 (ftype (function (cl-mpm/particle:particle cl-mpm/particle:particle) double-float) diff-squared))
(defun diff-squared (mp-a mp-b)
  (let ((dist (magicl:.- (cl-mpm/particle:mp-position mp-a)
                         (cl-mpm/particle:mp-position mp-b)))
        )
    (values (the double-float (magicl::sum (magicl:.* dist dist))))))
(defun weight-func (dist-squared length)
  (exp (- (* (/ 4 (* length length)) dist-squared))))
(defun weight-func-mps (mp-a mp-b length)
  (weight-func (diff-squared mp-a mp-b) length))

(defun calculate-delocalised-damage (mesh mp length)
  (let ((node-id (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp)))
        (node-reach (+ 1 (ceiling length (cl-mpm/mesh:mesh-resolution mesh)))))
    (let ((damage (cl-mpm/particle:mp-damage mp))
          (damage-inc 0d0)
          (mass-total 0d0))
      (loop for dx from (- node-reach) to node-reach
            do (loop for dy from (- node-reach) to node-reach
                     do (when (cl-mpm/mesh:in-bounds mesh (mapcar #'+ node-id (list dx dy)))
                          (let ((node (cl-mpm/mesh:get-node mesh (mapcar #'+ node-id (list dx dy)))))
                            (loop for mp-other across (cl-mpm/mesh::node-local-list node)
                                  do
                                     (with-accessors ((l-d cl-mpm/particle::mp-local-damage)
                                                      (m cl-mpm/particle:mp-mass)
                                                      (p cl-mpm/particle:mp-position))
                                         mp-other
                                       (let ((weight (weight-func-mps mp mp-other length)))
                                         (incf mass-total weight)
                                         (incf damage-inc
                                               (* (cl-mpm/particle::mp-local-damage-increment mp-other)
                                                  weight)
                                               )))
                                  )))))
      (setf damage-inc (/ damage-inc mass-total))
      (setf damage (max 0d0 (min 1d0 (+ damage damage-inc)))))))
(defun delocalise-damage (mesh mps dt len)
  ;; Move local damage into own temporary
  ;; (lparallel:pdotimes (i (length mps)) 
  ;;   (let ((mp (aref mps i)))
  ;;     (with-accessors ((damage cl-mpm/particle::mp-damage)
  ;;                      (local-damage cl-mpm/particle::mp-local-damage))
  ;;         mp
  ;;       (setf local-damage damage))))
  (lparallel:pdotimes (i (length mps)) 
    (let ((mp (aref mps i)))
      (with-accessors ((damage cl-mpm/particle::mp-damage)
                       (local-damage cl-mpm/particle::mp-local-damage))
          mp
        (setf damage (calculate-delocalised-damage mesh mp len))))))
