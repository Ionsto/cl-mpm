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
      ;; (/ (* (/ (max 0d0 stress) critical-stress) 1d-1) (max 1d-5 (expt (- 1d0 damage) 3)))
      ;(* (expt (/ (max 0d0 stress) critical-stress) 2d0) 1d-1 (/ 1 (max (/ 1 1) (expt (- 1d0 damage) 3))))
      ;; (* (expt (/ (max 0d0 stress) critical-stress) 2d0) 1d-1 (/ 1 (max (/ 1 100) (expt (- 1d0 damage) 3))))
      ;; (* (expt (/ (max 0d0 (- stress init-stress)) (- critical-stress init-stress)) 1d0) (/ 1 (max (/ 1 100) (expt (- 1d0 damage) 1.5))) rate)
      ;; (* (expt (/ (max 0d0 (- stress init-stress)) (- critical-stress init-stress)) 3d0)
      ;;    rate)
      (* (expt (max 0d0 (- stress init-stress)) 2d0) rate
         ;; (/ 1d0 (- 1d0 (min damage 0.4d0)))
         )
      ;; (* (expt (/ (max 0d0 stress) critical-stress) 2d0) 1d1)
      0d0)
  )

(defun damage-profile (damage damage-crit)
  "Constitive law describing the scalar stress decrease as a function of damage"
  (if (< damage damage-crit)
    (expt (- 1d0 damage) 1d0)
    0d0))

(defun calculate-damage-increment (mp dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (strain-rate cl-mpm/particle::mp-velocity-rate)
                     (critical-stress cl-mpm/particle:mp-critical-stress)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (pressure cl-mpm/particle::mp-pressure)
                     ) mp
        (progn
          (multiple-value-bind (l v) (magicl:eig (voight-to-matrix stress))
            ;; (let ((s_1 (apply #'max l)))
            ;;   (when (> s_1 0d0)
            ;;     (setf damage-increment (* (damage-rate-profile 0d0 s_1 damage 1d-8 init-stress) dt))))
            (let* ((l (sort l #'>))
                   ;;Effective principal stress
                   ;; (s_1 (- (nth 0 l) pressure))
                   (s_1 (- (nth 0 l) pressure))
                   ;; (s_1 (apply #' max l))
                   (s_v (sqrt (apply #'+ (mapcar (lambda (a b)
                                             (expt (- a b) 2))
                                           l
                                           (list (nth 1 l)
                                                 (nth 0 l))))))
                   (s_kk (magicl:trace (voight-to-matrix stress)))
                   (alpha 0.21d0)
                   (beta 0.63d0)
                   (hayhurst (+ (* alpha s_1)
                                (* beta s_v)
                                (* (- 1d0 alpha beta) s_kk)))
                   (B 5.23d-7)
                   (r 0.43d0)
                   )
              (when (> s_1 0d0)
                (multiple-value-bind (l v) (magicl:eig (voight-to-matrix strain-rate))
                  (let ((strain-rate (reduce #'+ (mapcar #'* l l))))
                    (incf damage-increment (* ;(expt (max 1d0 (* 1d4 strain-rate)) 1d0)
                                              (damage-rate-profile critical-stress
                                                                   ;;I think we calculate damage evolution on this?
                                                                   ;; (* s_1 (damage-profile damage critical-damage))
                                                                   s_1
                                                                   damage damage-rate init-stress) dt))))
                )
              ;; (when (> hayhurst 0d0)
              ;;   (incf damage-increment (* B (expt hayhurst r))))
              )
            ;(setf damage-increment (max 0d0 (min damage-increment (- critical-damage damage))))

            ;;Does transitioning from damage to fracture cause damage?
            ;; (when (> damage critical-damage)
            ;;   (setf damage 1d0))

            (when (>= damage 1d0)
              (setf damage-increment 0d0))
            (setf (cl-mpm/particle::mp-local-damage-increment mp)
                  damage-increment))))))
(declaim
 (inline apply-damage)
 (ftype (function (cl-mpm/particle:particle double-float) (values)) apply-damage))
(defun apply-damage (mp dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (damage-inc cl-mpm/particle::mp-damage-increment)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (pressure cl-mpm/particle::mp-pressure)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     ) mp
        (progn

          (incf damage damage-inc)
          (setf damage (max 0d0 (min 1d0 damage)))
          (when (> damage critical-damage)
            (setf damage 1d0))
          ;; (setf undamaged-stress (magicl:scale stress 1d0))
          ;; (setf undamaged-stress (magicl:scale stress (magicl:det def)))
          (when (> damage 0.0d0)
            ;; (setf stress
            ;;       (magicl:scale stress (- 1 damage)))
            ;;       (magicl:.- (magicl:scale stress (- 1 damage))
            ;;                  (magicl:from-list (list (* damage pressure)
            ;;                                          (* damage pressure)
            ;;                                          0d0) '(3 1))
            ;;                  )
            ;;       )
            (multiple-value-bind (l v) (magicl:eig
                                        (voight-to-matrix stress))
              (loop for i from 0 to 1
                    do (let* ((sii (nth i l))
                              (esii (- sii (* pressure damage)))
                              )
                         (when (> esii 0d0)
                           ;; (> sii 0d0)
                           (setf (nth i l)
                                 (+
                                  ;; (* sii (damage-profile damage critical-damage))
                                  ;; (* esii (damage-profile damage critical-damage))
                                  (* esii (- 1d0 damage))
                                  ;; pressure
                                  ;; 0d0
                                  (* damage pressure)
                                  )
                                 ))))
              (setf stress (matrix-to-voight (magicl:@ v
                                                       (magicl:from-diag l :type 'double-float)
                                                       (magicl:transpose v))))
              ))))
  (values)
  )
(defun update-damage (mp dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (damage-inc cl-mpm/particle:mp-damage)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (critical-stress cl-mpm/particle:mp-critical-stress)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     ) mp
        (progn
          ;; (setf undamaged-stress (magicl:scale stress 1d0))
          (multiple-value-bind (l v) (magicl:eig (voight-to-matrix undamaged-stress))
            ;; (loop for i from 0 to 1
            ;;       do (let ((sii (nth i l)))
            ;;            (progn
            ;;              (incf damage-increment (* (damage-rate-profile critical-stress sii damage) dt)))))
            (let ((s1 (loop for i from 0 to 1
                            maximize (nth i l))))
              (when (> s1 0d0)
                (incf damage-increment (* (damage-rate-profile critical-stress s1 damage damage-rate) dt))))
            (setf damage-increment (min damage-increment (- 1d0 damage)))
            (setf damage-inc damage-increment)
            (incf damage damage-increment)
            (setf damage (min 1d0 (max 0d0 damage)))
            ;;Jumping damage when hit a specific crack density
            (when (> damage critical-damage)
              (setf damage 1d0))

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
                           (when (> sii 0d0)
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
        (node-reach (+ 0 (ceiling length (cl-mpm/mesh:mesh-resolution mesh)))))
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
                                                        (m cl-mpm/particle:mp-volume)
                                                        (p cl-mpm/particle:mp-position))
                                           mp-other
                                         (let (
                                               (weight (weight-func-mps mp mp-other length))
                                               )
                                           (incf mass-total (* weight m
                                                               ))
                                           (incf damage-inc
                                                 (* (cl-mpm/particle::mp-local-damage-increment mp-other)
                                                    weight m)
                                                 )))))))))
      (setf damage-inc (* (/ damage-inc
                            mass-total)
                          )))))

(defun delocalise-damage (mesh mps dt len)
  ;; Calculate the delocalised damage for each damage particle
  (lparallel:pdotimes (i (length mps))
    (let ((mp (aref mps i)))
      (when (typep mp 'cl-mpm/particle:particle-damage)
        (with-accessors ((damage-inc cl-mpm/particle::mp-damage-increment)
                         (damage-inc-local cl-mpm/particle::mp-local-damage-increment)
                         (local-length cl-mpm/particle::mp-local-length)
                         )
            mp
          (setf damage-inc (calculate-delocalised-damage mesh mp local-length))
          )))))
