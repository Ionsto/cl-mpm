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
(defun damage-rate-profile (critical-stress stress damage rate init-stress)
  "Function that controls how damage evolves with principal stresses"
  (if (> stress init-stress)
      (* (expt (max 0d0 (- stress init-stress)) 2d0) rate)
      0d0))

(defun damage-profile (damage damage-crit)
  "Constitive law describing the scalar stress decrease as a function of damage"
  (if (< damage damage-crit)
    (expt (- 1d0 damage) 1d0)
    0d0))

(defun calculate-y (de strain)
  (magicl:tref
   (magicl:@ (magicl:transpose strain) de strain)
   0 0))

(defun calculate-e-mazar (strain)
  (magicl:tref
   (magicl:@ (magicl:transpose strain) de strain)
   0 0))
(defun calculate-e-rankine (stress E)
  (multiple-value-bind (l v) (magicl:eig (voight-to-matrix stress))
    (/ (apply #'max l) E)))

(defun damage-profile-linear (k e-0 e-f)
  (* (/ e-f (- e-f e-0)) (- 1d0 (/ e-0 (max k e-0)))))

(defun damage-profile-exp (k e-0 e-f)
  (- 1d0 (* (/ e-0 k) (exp (- (/ (- k e-0) (- e-f e-0)))))))

(defun calculate-damage-increment (mp dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (E-modulus cl-mpm/particle::mp-E)
                     (def cl-mpm/particle:mp-deformation-gradient)
                     (e cl-mpm/particle::mp-effective-strain)
                     ) mp
        (progn
          (setf e (calculate-e-rankine (magicl:scale stress (/ 1 (magicl:det def))) E-modulus))
          ))))
(declaim
 (inline apply-damage)
 (ftype (function (cl-mpm/particle:particle double-float) (values)) apply-damage))
(defun apply-damage (mp dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (damage-inc cl-mpm/particle::mp-damage-increment)
                     (e-bar cl-mpm/particle::mp-e-bar)
                     (k cl-mpm/particle::mp-k)
                     (e-0 cl-mpm/particle::mp-e-0)
                     (e-f cl-mpm/particle::mp-e-f)
                     (pressure cl-mpm/particle::mp-pressure)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     ) mp
        (progn
          (let ((d-0 damage))
            ;;Update our internal strain parameter
            (when (> e-bar k)
              (setf k e-bar)
              (let ((k-linear (- (exp k) 1)))
                (setf damage (damage-profile-exp k e-0 e-f))))
            ;;Apply critical damage parameter
            (when (> damage critical-damage)
              (setf damage 1d0))
            ;; Clamp damage
            (setf damage (max 0d0 (min 1d0 damage)))
            (setf damage-inc (- damage d-0)))

          ;;Only do stress decomposition when damage is nonzero
          (when (> damage 0.0d0)
            ;;Apply damage onto stress
            (multiple-value-bind (l v) (magicl:eig
                                        (voight-to-matrix stress))
              (loop for i from 0 to 1
                    do (let* ((sii (nth i l)))
                         (when (> sii 0d0)
                           (setf (nth i l) (* sii (- 1d0 damage))))))
              (setf stress (matrix-to-voight (magicl:@
                                              v
                                              (magicl:from-diag l :type 'double-float)
                                              (magicl:transpose v))))))))
  (values)
  )

(defmethod cl-mpm/particle:post-stress-step (mesh (mp cl-mpm/particle:particle-damage) dt)
  ;(update-damage mp dt)
  )

(defun calculate-damage (mesh mps dt len)
    (create-delocalisation-list mesh mps len)
    (lparallel:pdotimes (i (length mps))
      (when (typep (aref mps i) 'cl-mpm/particle:particle-damage)
        (calculate-damage-increment (aref mps i) dt)))
    (delocalise-damage mesh mps dt len)
    (lparallel:pdotimes (i (length mps))
      (when (typep (aref mps i) 'cl-mpm/particle:particle-damage)
        (apply-damage (aref mps i) dt))))

(defun create-delocalisation-list (mesh mps length)
  (with-accessors ((nodes cl-mpm/mesh:mesh-nodes))
        mesh
      (lparallel:pdotimes (i (array-total-size nodes))
        (let ((node (row-major-aref nodes i)))
          (setf (fill-pointer (cl-mpm/mesh::node-local-list node)) 0)))
    (lparallel:pdotimes (i (length mps))
      (let ((mp (aref mps i)))
        (when (typep (aref mps i) 'cl-mpm/particle:particle-damage)
            (let ((node-id (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp))))
              (when (cl-mpm/mesh:in-bounds mesh node-id)
                (let ((node (cl-mpm/mesh:get-node mesh node-id)))
                (sb-thread:with-mutex ((cl-mpm/mesh:node-lock node))
                  (vector-push-extend mp (cl-mpm/mesh::node-local-list node)))))))))
    (lparallel:pdotimes (i (array-total-size nodes))
      (let ((node (row-major-aref nodes i)))
        (when (= 0 (length (cl-mpm/mesh::node-local-list node)))
          (adjust-array (cl-mpm/mesh::node-local-list node) 0 :fill-pointer 0))))))

#- :sb-simd
(progn
  (declaim
   (inline diff-squared)
   (ftype (function (cl-mpm/particle:particle cl-mpm/particle:particle) double-float) diff-squared))
  (defun diff-squared (mp-a mp-b)
    (let ((dist (magicl:.- (cl-mpm/particle:mp-position mp-a)
                           (cl-mpm/particle:mp-position mp-b)))
          )
      (values (the double-float (magicl::sum (magicl:.* dist dist)))))))

;;This is a simd dot product
#+ :sb-simd
(progn
  (declaim
   (inline simd-accumulate)
   (ftype (function ((simple-array double-float) (simple-array double-float)) (values double-float)) simd-accumulate))
  (defun simd-accumulate (a b)
    (declare (type sb-simd:f64vec a b))
    (let ((diff
            (sb-simd-avx:f64.2-
             (sb-simd-avx:f64.2-aref a 0)
             (sb-simd-avx:f64.2-aref b 0)
             )))
      (values (sb-simd-avx::f64.2-horizontal+
               (sb-simd-avx:f64.2* diff diff))))
    )
  (declaim
   (inline diff-squared)
   (ftype (function (cl-mpm/particle:particle cl-mpm/particle:particle) double-float) diff-squared))
  (defun diff-squared (mp-a mp-b)
    (let ((pos-a (magicl::storage (cl-mpm/particle:mp-position mp-a)))
          (pos-b (magicl::storage (cl-mpm/particle:mp-position mp-b)))
          )
      (values (the double-float (simd-accumulate pos-a pos-b))))))

(declaim
 (inline weight-func)
 (ftype (function (double-float
                   double-float
                   ) (values double-float))
        weight-func
        ))
(defun weight-func (dist-squared length)
  (values (the double-float (exp (- (* (/ 4d0 (* length length)) dist-squared))))))
(declaim
 (inline weight-func-mps)
 (ftype (function (cl-mpm/particle:particle
                   cl-mpm/particle:particle
                   double-float
                   ) double-float)
        weight-func-mps
        ))
(defun weight-func-mps (mp-a mp-b length)
  (weight-func (diff-squared mp-a mp-b) length))

(declaim
 (inline calculate-delocalised-damage)
 (ftype (function (cl-mpm/mesh::mesh
                   cl-mpm/particle:particle-damage
                   double-float
                   ) double-float)
        calculate-delocalised-damage
        ))
(defun calculate-delocalised-damage (mesh mp length)
  (let ((node-id (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp)))
        (node-reach (+ 0 (ceiling (* length 2d0) (cl-mpm/mesh:mesh-resolution mesh)))))
    (let ((damage-inc 0d0)
          (mass-total 0d0))
      (loop for dx from (- node-reach) to node-reach
            do (loop for dy from (- node-reach) to node-reach
                     do
                        (let ((idx (mapcar #'+ node-id (list dx dy))))
                          (when (cl-mpm/mesh:in-bounds mesh idx)
                            (let ((node (cl-mpm/mesh:get-node mesh idx)))
                              (loop for mp-other across (cl-mpm/mesh::node-local-list node)
                                    do
                                       (with-accessors ((l-d cl-mpm/particle::mp-local-damage)
                                                        (v cl-mpm/particle:mp-volume)
                                                        (p cl-mpm/particle:mp-position))
                                           mp-other
                                         (let ((weight (* v (weight-func-mps mp mp-other length))))
                                           (incf mass-total weight)
                                           (incf damage-inc
                                                 (* (cl-mpm/particle::mp-effective-strain mp-other)
                                                    weight)
                                                 )))))))))
      (setf damage-inc (/ damage-inc mass-total)))))

(defun delocalise-damage (mesh mps dt len)
  "Calculate the delocalised effective strain of each MP"
  (lparallel:pdotimes (i (length mps))
    (let ((mp (aref mps i)))
      (when (typep mp 'cl-mpm/particle:particle-damage)
        (with-accessors ((e-bar cl-mpm/particle::mp-e-bar)
                         (local-length cl-mpm/particle::mp-local-length)
                         )
            mp
          (setf e-bar (calculate-delocalised-damage mesh mp local-length)))))))
