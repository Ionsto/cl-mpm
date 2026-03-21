(in-package :cl-mpm/damage)
(declaim
 (notinline diff-damaged)
 (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle cl-mpm/particle:particle) double-float) diff-damaged))
(defun diff-damaged (mesh mp-a mp-b)
  (labels ((node-dam (node)
             (if (cl-mpm/mesh::node-active node)
                 (min 1d0 (max 0d0 (the double-float (cl-mpm/mesh::node-damage node))))
                 0d0))
           (get-damage (position)
             (let ((damage 0d0)
                   (weight 0d0))
               (cl-mpm::iterate-over-neighbours-point-linear
                mesh
                position
                (lambda (m node weight grads)
                  (declare (double-float damage weight))
                  (if (cl-mpm/mesh::node-active node)
                    (incf damage
                          (* (node-dam node)
                             weight))
                    (incf damage (* 1d0 weight)))))
               damage)
             )
           )
    (with-accessors
          ((h cl-mpm/mesh::mesh-resolution))
        mesh
      (let* ((pos-a (cl-mpm/particle:mp-position mp-a))
             (pos-b (cl-mpm/particle:mp-position mp-b))
             (diff (cl-mpm/utils:vector-zeros))
             (epsilon 1d-8)
             (final-distance 0d0)
             )
        (declare (double-float h final-distance epsilon))
        (cl-mpm/fastmaths:fast-.- pos-b pos-a diff)
        (let* ((length (sqrt (diff-squared mp-a mp-b))))

          ;; Sample damage at midpoint and integrated with constant damage assumption
          ;; Can utilise linear shape function damage assumption but produced a nasty looking intergral
          (declare (double-float h length final-distance))
          (when (> length 0d0)
            (let* ((step-norm (cl-mpm/utils::vector-copy diff))
                   ;; (step-norm (magicl:scale diff (/ 1d0 length)))
                   ;;Start at point a and step through to b
                   (step-point (cl-mpm/utils::vector-copy pos-a))
                   ;;Resolution of our midpoint integration
                   (step-size (/ h 2d0))
                   )
              (cl-mpm/fastmaths:fast-scale! step-norm (/ 1d0 length))
              (multiple-value-bind (steps remainder) (floor length step-size)
                (when (> steps 0)
                  (let* ((dhstep (cl-mpm/utils::vector-copy step-norm)))
                    (cl-mpm/fastmaths:fast-scale! dhstep (* 0.5d0 step-size))
                    (loop for i fixnum from 0 below steps
                          do
                             (progn
                               (cl-mpm/fastmaths::fast-.+ step-point dhstep step-point)
                               (let ((damage (get-damage step-point)))
                                 (declare (double-float damage))
                                 (incf final-distance
                                       (* step-size
                                          (/ 1d0 (max epsilon (sqrt (max 0d0 (min 1d0 (- 1d0 damage)))))))))
                               (cl-mpm/fastmaths::fast-.+ step-point dhstep step-point)
                               )
                          )))
                (let* ((step-size remainder)
                       (dhstep (cl-mpm/utils::vector-copy step-norm)))
                  (magicl:scale! dhstep (* 0.5d0 step-size))
                  (progn
                    (cl-mpm/fastmaths::fast-.+ step-point dhstep step-point)
                    (let ((damage 0d0))
                      (cl-mpm::iterate-over-neighbours-point-linear
                       mesh
                       step-point
                       (lambda (m node weight grads)
                         (declare (double-float damage weight))
                         (incf damage
                               (* (node-dam node) weight))))
                      (incf final-distance (* step-size
                                              (/ 1d0
                                                 (max epsilon
                                                      (the double-float (sqrt (max 0d0 (min 1d0 (- 1d0 damage)))))))))))))))
          (max 0d0 (* final-distance final-distance)))))))


(declaim
 (notinline weight-func-mps-damaged)
 (ftype (function (cl-mpm/mesh::mesh
                   cl-mpm/particle:particle
                   cl-mpm/particle:particle
                   double-float
                   ) double-float)
        weight-func-mps-damaged
        ))
(defun weight-func-mps-damaged (mesh mp-a mp-b length)
  (weight-func (diff-damaged mesh mp-a mp-b) length))

(defun delocalise-damage-ekl (sim)
  ;; (declare (optimize (speed 3)))
  (with-accessors ((mps cl-mpm::sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (g2p-damage sim)
    ;; Calculate the delocalised damage for each damage particle
    (cl-mpm:iterate-over-mps
     mps
     (lambda (mp)
       (when (typep mp 'cl-mpm/particle:particle-damage)
         (with-accessors ((damage-ybar cl-mpm/particle::mp-damage-ybar)
                          (damage cl-mpm/particle::mp-damage)
                          (local-length-t cl-mpm/particle::mp-local-length))
             mp
           (setf damage-ybar (calculate-delocalised-damage-ekl mesh mp local-length-t)))))))
  (values))

(defun calculate-delocalised-damage-ekl (mesh mp length)
  (let ((damage-inc 0d0)
        (mass-total 0d0))
    (declare (double-float damage-inc mass-total))
    (iterate-over-neighour-mps
     mesh mp
     (cl-mpm/particle::mp-local-length mp)
     (lambda (mesh mp mp-other dist)
       (with-accessors ((d cl-mpm/particle::mp-damage)
                        (m cl-mpm/particle:mp-volume)
                        (ll cl-mpm/particle::mp-local-length)
                        (p cl-mpm/particle:mp-position))
           mp-other
         (let ((weight (weight-func-mps-damaged mesh mp mp-other (cl-mpm/particle::mp-local-length mp))))
           (declare (double-float weight m d mass-total damage-inc))
           (incf mass-total (* weight m))
           (incf damage-inc
                 (* (the double-float (cl-mpm/particle::mp-damage-y-local mp-other))
                    weight m))))))
    (when (> mass-total 0d0)
      (setf damage-inc (/ damage-inc mass-total)))
    damage-inc))


(defun calculate-debug-ekl (mesh mp length)
  )
