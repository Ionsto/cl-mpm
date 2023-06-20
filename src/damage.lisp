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
(declaim
 ;(inline damage-rate-profile)
 (ftype (function (double-float double-float double-float double-float ) (double-float)) damage-rate-profile))
(defun damage-rate-profile (stress damage rate init-stress)
  (declare (double-float stress damage rate init-stress))
  "Function that controls how damage evolves with principal stresses"
  (if (> stress init-stress)
      ;(* (expt (max 0d0 (- stress init-stress)) 0.43d0) rate)
      ;(* (expt (max 0d0 (- stress init-stress)) 0.50d0) rate)
      (* (expt (max 0d0 (- stress init-stress)) 1d0) rate)
      0d0))

(defun damage-profile (damage damage-crit)
  "Constitive law describing the scalar stress decrease as a function of damage"
  (if (< damage damage-crit)
    (expt (- 1d0 damage) 1d0)
    0d0))

(declaim
 ;(inline calculate-damage-increment)
 (ftype (function (cl-mpm/particle:particle double-float) (values)) calculate-damage-increment))
(defun calculate-damage-increment (mp dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
                     ;(stress cl-mpm/particle::mp-stress)
                     ;; (strain cl-mpm/particle::mp-strain)
                     (damage cl-mpm/particle:mp-damage)
                     (strain-rate cl-mpm/particle::mp-velocity-rate)
                     (critical-stress cl-mpm/particle:mp-critical-stress)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (pressure cl-mpm/particle::mp-pressure)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     ;(damage-driving-factor cl-mpm/particle::mp-damage-driving-factor)
                     ) mp
      (declare (double-float pressure damage))
        (progn
          (progn
                                        ;multiple-value-bind (l v) (magicl:eig (magicl:scale (voight-to-matrix stress) (/ 1d0 (magicl:det def))))
            (let* (;(l (sort l #'>))
                   ;(s_1 (nth 0 l))
                   ;; (s_2 (nth 1 l))
                   ;(cauchy-stress (magicl:scale! (voight-to-matrix stress) (/ 1d0 (magicl:det def))))
                   ;; (cauchy-stress (magicl:scale! (voight-to-matrix stress) ))

                   ;(j (/ 1d0 (the double-float (magicl:det def))))
                   (j 1d0)
                   (av (/ (* j (+ (the double-float (magicl:tref stress 0 0))
                                  (the double-float (magicl:tref stress 1 0)))) 2))
                   (diff (the double-float
                              (sqrt (the double-float
                                         (+ (the double-float
                                                 (expt
                                                  (/
                                                   (* j (the double-float
                                                             (- (the double-float (magicl:tref stress 0 0))
                                                                (the double-float (magicl:tref stress 1 0)))))
                                                   2d0)
                                                  2d0))
                                            (the double-float
                                                 (expt (* j (the double-float (magicl:tref stress 2 0))) 2d0)))))))
                   (s_1 (+ av diff))
                   (s_2 (- av diff))
                   (pressure-effective (* pressure 1d0))
                   (s_1 (- s_1 pressure-effective))
                   (s_2 (- s_2 pressure-effective))
                   (s_1 (max 0d0 s_1))
                   (s_2 (max 0d0 s_2))
                   ;; (vm (* (sqrt (/ 3 4)) (- s_1 s_2)))
                   (vm (- s_1 s_2))
                   (s_1 vm)
                   ;(damage-inv (- 1d0 damage))
                   )
              (when (> s_1 0d0)
                ;;(setf damage-increment (* s_1 (- 1d0 damage)))
                (setf damage-increment s_1)
                ;; (if (< damage 1d0)
                ;;   (setf damage-increment (/ s_1 (- 1d0 damage)))
                ;;   (setf damage-increment s_1)
                ;;   )
                ))
            (when (>= damage 1d0)
              (setf damage-increment 0d0))

            ;; (when (= damage 1d0)
            ;;   (setf damage-increment 0d0))
            ;;(setf ybar damage-increment)

            ;; (setf damage-increment (* dt (damage-rate-profile damage-increment damage damage-rate init-stress)))
            (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
            ))))
  (values))

(declaim
 (inline apply-damage)
 (ftype (function (cl-mpm/particle:particle double-float) (values)) apply-damage))
(defun apply-damage (mp dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (damage-inc cl-mpm/particle::mp-damage-increment)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (pressure cl-mpm/particle::mp-pressure)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     ) mp
      (declare (double-float damage damage-inc critical-damage))
        (progn
          ;;Damage increment holds the delocalised driving factor
          (setf ybar damage-inc)
          ;; (when (< damage 1d0)
          (setf damage-inc (* damage-inc (/ 1d0 (expt (- 1d0 damage) 1d0))));3
          (setf damage-inc (* dt (damage-rate-profile damage-inc damage damage-rate init-stress)))
          (when (>= damage 1d0)
            (setf damage-inc 0d0)
            (setf ybar 0d0))
          (incf damage damage-inc)
          (setf damage (max 0d0 (min 1d0 damage)))
          (when (> damage critical-damage)
            (setf damage 1d0)
            (setf damage-inc 0d0))
          ;; (when (> damage 0.0d0)
          ;;   (multiple-value-bind (l v) (magicl:eig
          ;;                               (voight-to-matrix stress))
          ;;     (loop for i from 0 to 1
          ;;           do (let* ((sii (nth i l))
          ;;                     (esii (- sii
          ;;                              (* pressure 1)))
          ;;                     )
          ;;                (when (> esii 0d0)
          ;;                  (setf (nth i l)
          ;;                        (+
          ;;                         (* esii (- 1d0 damage))
          ;;                         (* 1 pressure)
          ;;                         ))
          ;;                  ))
          ;;     (setf stress (matrix-to-voight (magicl:@ v
          ;;                                              (magicl:from-diag l :type 'double-float)
          ;;                                              (magicl:transpose v))))
          ;;     )))
          )
  (values)
  ))
(defmethod cl-mpm/particle:post-stress-step (mesh (mp cl-mpm/particle:particle-damage) dt)
  ;(update-damage mp dt)
  )

(defparameter *delocal-counter* 0)
(defparameter *delocal-counter-max* 1)
(defun calculate-damage (mesh mps dt len non-local-damage)
  (when non-local-damage
    (when (<= *delocal-counter* 0)
      (create-delocalisation-list mesh mps len)
      (setf *delocal-counter* *delocal-counter-max*))
    (decf *delocal-counter*))
  (lparallel:pdotimes (i (length mps))
                      (when (typep (aref mps i) 'cl-mpm/particle:particle-damage)
                        (calculate-damage-increment (aref mps i) dt)))
  (if non-local-damage
    (delocalise-damage mesh mps dt len)
    (localise-damage mesh mps dt))
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
        (when (typep mp 'cl-mpm/particle:particle-damage)
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
    (let ((pos-a (magicl::matrix/double-float-storage (cl-mpm/particle:mp-position mp-a)))
          (pos-b (magicl::matrix/double-float-storage (cl-mpm/particle:mp-position mp-b)))
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
  (values (the double-float (exp (the double-float (- (* (/ 4d0 (* length length)) dist-squared)))))))
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
        (node-reach (the fixnum (+ 0 (truncate (ceiling (* length 2d0)
                                                        (the double-float (cl-mpm/mesh:mesh-resolution mesh))))))))
    (declare (dynamic-extent node-id))
    (let ((damage-inc 0d0)
          (mass-total 0d0))
      (declare (double-float damage-inc mass-total))
      (loop for dx from (- node-reach) to node-reach
            do (loop for dy from (- node-reach) to node-reach
                     do
                        (let ((idx (mapcar #'+ node-id (list dx dy))))
                          (declare (dynamic-extent idx))
                          (when (cl-mpm/mesh:in-bounds mesh idx)
                            (let ((node (cl-mpm/mesh:get-node mesh idx)))
                              (loop for mp-other across (the (vector T *) (cl-mpm/mesh::node-local-list node))
                                    do
                                       (with-accessors ((d cl-mpm/particle::mp-damage)
                                                        (m cl-mpm/particle:mp-volume)
                                                        (p cl-mpm/particle:mp-position))
                                           mp-other
                                         (when (< (the double-float d) 1d0)
                                           (let ((weight (weight-func-mps mp mp-other length)))
                                             (declare (double-float weight m d mass-total damage-inc))
                                             (incf mass-total (* weight m))
                                             (incf damage-inc
                                                   (* (the double-float (cl-mpm/particle::mp-local-damage-increment mp-other))
                                                      weight m)))))))))))
      (when (> mass-total 0d0)
        (setf damage-inc (/ damage-inc mass-total)))
      damage-inc
      ;; (setf damage-inc (cl-mpm/particle::mp-local-damage-increment mp))
      )))

(defun length-localisation (local-length local-length-damaged damage)
  (+ (* local-length (- 1d0 damage)) (* local-length-damaged damage)))
(declaim
 (ftype
  (function (cl-mpm/mesh::mesh
             (array cl-mpm/particle::particle)
             double-float
             double-float
             )
            (values))
  delocalise-damage))
(defun delocalise-damage (mesh mps dt len)
  ;; Calculate the delocalised damage for each damage particle
  (lparallel:pdotimes (i (length mps))
    (let ((mp (aref mps i)))
      (when (typep mp 'cl-mpm/particle:particle-damage)
        (with-accessors ((damage-inc cl-mpm/particle::mp-damage-increment)
                         (damage-inc-local cl-mpm/particle::mp-local-damage-increment)
                         (damage cl-mpm/particle::mp-damage)
                         (local-length cl-mpm/particle::mp-local-length)
                         (local-length-damaged cl-mpm/particle::mp-local-length-damaged)
                         (local-length-t cl-mpm/particle::mp-true-local-length)
                         )
            mp
          (setf local-length-t (length-localisation local-length local-length-damaged damage))
          (setf damage-inc (calculate-delocalised-damage mesh mp local-length-t))
          ))))
  (values))

(declaim
 (ftype
  (function (cl-mpm/mesh::mesh
             (array cl-mpm/particle::particle)
             double-float
             )
            (values))
  localise-damage))
(defun localise-damage (mesh mps dt)
  ;(declare ((array cl-mpm/particle::particle *) mps))
  "Apply local damage model"
  (lparallel:pdotimes (i (length mps))
                      (let ((mp (aref mps i)))
                        (when (typep mp 'cl-mpm/particle:particle-damage)
                          (with-accessors ((damage-inc cl-mpm/particle::mp-damage-increment)
                                           (damage-inc-local cl-mpm/particle::mp-local-damage-increment)
                                           (local-length cl-mpm/particle::mp-local-length)
                                           )
                              mp
                            (setf damage-inc damage-inc-local)))))
  (values))
