(defpackage :cl-mpm/damage
  (:use :cl
        :cl-mpm/utils)
  (:export
    #:update-damage
    #:post-stress-step
    )
  )
;; (defconstant +damage-exp+ 0.5d0)
;; (defconstant +stress-exp+ 2d0)

(in-package :cl-mpm/damage)
(declaim (optimize (debug 3) (safety 3) (speed 0)))
(declaim
 ;(inline damage-rate-profile)
 (ftype (function (double-float double-float double-float double-float ) (double-float)) damage-rate-profile))
(defun damage-rate-profile (stress damage rate init-stress)
  (declare (double-float stress damage rate init-stress))
  "Function that controls how damage evolves with principal stresses"
  (if (> stress init-stress)
      ;(* (expt (max 0d0 (- stress init-stress)) 0.43d0) rate)
      ;(* (expt (max 0d0 (- stress init-stress)) 0.50d0) rate)
      (* (expt (max 0d0 (/ (- stress init-stress) init-stress)) 2d0) rate)
      ;; (* (expt (max 0d0 (- stress init-stress)) 3d0) rate)
      0d0))

(defun principal-stresses (stress)
  (multiple-value-bind (l v) (cl-mpm/utils::eig (voight-to-matrix stress))
    (declare (ignore v))
    (values (apply #'max l) (apply #'min l))))

(defun principal-stresses-3d (stress)
  (multiple-value-bind (l v) (cl-mpm/utils::eig (voight-to-matrix stress))
    (declare (ignore v))
    (setf l (sort l #'>))
    (values (nth 0 l) (nth 1 l) (nth 2 l))))

(defun damage-profile (damage damage-crit)
  "Constitive law describing the scalar stress decrease as a function of damage"
  (if (< damage damage-crit)
    (expt (- 1d0 damage) 2d0)
    0d0))

(declaim
 ;(inline calculate-damage-increment)
 (ftype (function (cl-mpm/particle:particle double-float) (values)) calculate-damage-increment))
(defun calculate-damage-increment (mp dt)
  (let ((damage-increment 0d0))
    (with-accessors (;(stress cl-mpm/particle::mp-stress)
                     (stress cl-mpm/particle::mp-undamaged-stress)
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
                     (damage-tensor cl-mpm/particle::mp-damage-tensor)
                     (ybar-tensor cl-mpm/particle::mp-damage-ybar-tensor)
                     (local-length cl-mpm/particle::mp-local-length)
                     (local-length-damaged cl-mpm/particle::mp-local-length-damaged)
                     (local-length-t cl-mpm/particle::mp-true-local-length)
                     ) mp
      (declare (double-float pressure damage))
        (progn
          (progn
                                        ;multiple-value-bind (l v) (cl-mpm/utils::eig (magicl:scale (voight-to-matrix stress) (/ 1d0 (magicl:det def))))
            (when (< damage 1d0)
              (let ((cauchy-undamaged (magicl:scale stress (/ 1d0 (* 1d0;(- 1d0 damage)
                                                                     (magicl:det def))))))
                (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d cauchy-undamaged)
                  (let* (;;Only allow tensile damage
                         (pressure-effective (* 1d0 damage pressure))
                         ;; (pressure-effective (* 1d0 pressure))
                         (j2 (sqrt (cl-mpm/constitutive::voigt-j2 (cl-mpm/constitutive::deviatoric-voigt
                                                                   cauchy-undamaged))))
                         (p (* 0.3d0 (+ s_1 s_2 s_3)))
                         (s_1 (- s_1 pressure-effective))
                         (s_2 (- s_2 pressure-effective))
                         (s_3 (- s_3 pressure-effective))
                         (s_1 (max 0d0 s_1))
                         (s_2 (max 0d0 s_2))
                         (s_3 (max 0d0 s_3))
                         ;; (vm (* (sqrt (/ 3 4)) (- s_1 s_2)))
                         ;; (vm (sqrt (* 0.5d0
                         ;;              (+ (expt (- s_1 s_2) 2)
                         ;;                 (expt (- s_2 s_3) 2)
                         ;;                 (expt (- s_3 s_1) 2)))))
                         (s_1 (- j2 (* 0.2d0 p)))

                         ;; (s_1 (- s_1 s_3))
                         ;; (s_1 vm)
                                        ;(damage-inv (- 1d0 damage))
                         )
                    (when (> s_1 0d0)
                      (setf damage-increment (* (- 1d0 damage) s_1))
                      ;; (if (< damage 1d0)
                      ;;     (setf damage-increment (/ s_1 (expt (- 1d0 damage) 0.5d0)))
                      ;;     (setf damage-increment s_1)
                      ;;     )
                      )))))
            (when (>= damage 1d0)
              (setf damage-increment 0d0))

            ;;Delocalisation switch
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
                     (log-damage cl-mpm/particle::mp-log-damage)
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
          ;; (setf damage-inc (* damage-inc (/ 1d0 (expt (- 1d0 damage) 1d0))));3
          ;;Normal
          ;; (setf damage-inc (* dt (- critical-damage damage) (damage-rate-profile damage-inc damage damage-rate init-stress)))
          (setf damage-inc (* dt (damage-rate-profile-chalk damage-inc damage damage-rate init-stress)))
          ;;Mohr coloumb
          ;; (let ((angle 30d0))
          ;;   (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d (magicl:scale stress (/ 1d0 (magicl:det def))))
          ;;     (let ((init-stress (* (+ s_1 s_3) (sin angle))))
          ;;       (setf damage-inc (* dt (damage-rate-profile damage-inc damage damage-rate init-stress))))
          ;;     ))


          (when (>= damage 1d0)
            (setf damage-inc 0d0)
            (setf ybar 0d0))
          (incf (cl-mpm/particle::mp-time-averaged-damage-inc mp) damage-inc)
          (incf (cl-mpm/particle::mp-time-averaged-ybar mp) ybar)
          (incf (cl-mpm/particle::mp-time-averaged-counter mp))
          ;;Transform to log damage
          (incf damage damage-inc)
          ;;Transform to linear damage
          (setf damage (max 0d0 (min 1d0 damage)))
          (when (> damage critical-damage)
            (setf damage 1d0)
            (setf damage-inc 0d0)))
  (values)
  ))
(defmethod cl-mpm/particle:post-stress-step (mesh (mp cl-mpm/particle:particle-damage) dt)
  ;(update-damage mp dt)
  )
(defun find-nodal-local-length (mesh mp)
  (let ((nodal-damage 0d0))
    (cl-mpm::iterate-over-neighbours
     mesh mp
     (lambda (mesh mp node svp grads fsvp fgrads)
       (with-accessors ((node-damage cl-mpm/mesh::node-damage)
                        (node-svp cl-mpm/mesh::node-svp-sum))
           node
         (when (> node-svp 0d0)
           (incf nodal-damage (/ (* node-damage svp)
                                 node-svp))))
       ))
    (with-accessors ((local-length cl-mpm/particle::mp-local-length)
                     (local-length-damaged cl-mpm/particle::mp-local-length-damaged)
                     (local-length-t cl-mpm/particle::mp-true-local-length)
                     ) mp
      ;;Damage length control
      (setf local-length-t (length-localisation local-length local-length-damaged nodal-damage))
      ;; (setf local-length-t local-length)
      )))

(defparameter *delocal-counter* 0)
(defparameter *delocal-counter-max* 10)
(defun calculate-damage (mesh mps dt len non-local-damage)
  (when non-local-damage
    (when (<= *delocal-counter* 0)
      ;; (create-delocalisation-list mesh mps)
      ;;Ensure they have a home
      (update-delocalisation-list mesh mps)
      (setf *delocal-counter* *delocal-counter-max*))
    (decf *delocal-counter*))
  (lparallel:pdotimes (i (length mps))
    (let ((mp (aref mps i)))
      (when (typep mp 'cl-mpm/particle:particle-damage)
        ;(find-nodal-local-length mesh mp)
        ;(setf (cl-mpm/particle::mp-true-local-length mp) (cl-mpm/particle::mp-local-length mp))
        (setf (cl-mpm/particle::mp-true-local-length mp)
              (length-localisation (cl-mpm/particle::mp-local-length mp)
                                   (cl-mpm/particle::mp-local-length-damaged mp)
                                   (cl-mpm/particle::mp-damage mp)))
        ;; (calculate-damage-increment (aref mps i) dt)
        (damage-model-calculate-y mp dt)
        )))
  (if non-local-damage
    (delocalise-damage mesh mps dt len)
    (localise-damage mesh mps dt))
  (lparallel:pdotimes (i (length mps))
    (let ((mp (aref mps i)))
      (when (typep mp 'cl-mpm/particle:particle-damage)
        ;; (find-nodal-local-length mesh (aref mps i))
                                        ;(apply-damage (aref mps i) dt)
        (update-damage mp dt)
        ;; (damage-model-update-damage mp (cl-mpm/particle::mp-damage-model mp) dt)
        ))))
(defun create-delocalisation-list (mesh mps)
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

(defun local-list-add-particle (mesh mp)
  "A function for putting an MP into the nodal MP lookup table"
  (let ((node-id (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp))))
    (when (cl-mpm/mesh:in-bounds mesh node-id)
      (let ((node (cl-mpm/mesh:get-node mesh node-id)))
        (sb-thread:with-mutex ((cl-mpm/mesh:node-lock node))
          (setf (cl-mpm/particle::mp-damage-position mp) (magicl:copy-tensor (cl-mpm/particle:mp-position mp)))
          (vector-push-extend mp (cl-mpm/mesh::node-local-list node)))))))

(defun local-list-remove-particle (mesh mp)
  "A function for removing an MP into the nodal MP lookup table"
  (flet ((remove-mp-ll (node)
           (with-accessors ((ll cl-mpm/mesh::node-local-list)
                            (lock cl-mpm/mesh:node-lock)
                            ) node
             (if (position mp ll)
                 (progn
                   (sb-thread:with-mutex (lock)
                     (setf ll (delete mp ll)))
                   t)
                 nil))))
    ;;Check if the particle has been inserted by checking nil equality of damage position
    (when (cl-mpm/particle::mp-damage-position mp)
      (let ((node-id (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle::mp-damage-position mp))))
        (when (cl-mpm/mesh:in-bounds mesh node-id)
          (let ((node (cl-mpm/mesh:get-node mesh node-id)))
            (if (remove-mp-ll node)
                t
                (cl-mpm::iterate-over-nodes-serial
                 mesh
                 (lambda (node)
                   (remove-mp-ll node))))))))))

(defun update-delocalisation-list (mesh mps)
  (with-accessors ((nodes cl-mpm/mesh:mesh-nodes)
                   (h cl-mpm/mesh:mesh-resolution))
        mesh
      ;; (lparallel:pdotimes (i (array-total-size nodes))
      ;;   (let ((node (row-major-aref nodes i)))
      ;;     (setf (fill-pointer (cl-mpm/mesh::node-local-list node)) 0)))
    (lparallel:pdotimes (i (length mps))
      (let ((mp (aref mps i)))
        (when (typep mp 'cl-mpm/particle:particle-damage)
          (if (eq (cl-mpm/particle::mp-damage-position mp) nil)
              ;;New particle - needs to be added to the mesh
              (local-list-add-particle mesh mp)
              ;; Already inserted mesh - sanity check to see if it should be recalced
              (let* ((delta (diff-squared-mat (cl-mpm/particle:mp-position mp)
                                              (cl-mpm/particle::mp-damage-position mp)
                                             )))
                (declare (double-float delta h))
                (when (> delta (/ h 4d0))
                  (when (not (eq
                              (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp))
                              (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle::mp-damage-position mp))
                              ))
                  (local-list-remove-particle mesh mp)
                  (local-list-add-particle mesh mp))))))))
    (lparallel:pdotimes (i (array-total-size nodes))
      (let ((node (row-major-aref nodes i)))
        (when (= 0 (length (cl-mpm/mesh::node-local-list node)))
          (adjust-array (cl-mpm/mesh::node-local-list node) 0 :fill-pointer 0))))))

;; (defgeneric update-delocalisation-list (sim))
;; (defmethod update-delocalisation-list ((sim cl-mpm::mpm-sim))
;;   (update-delocalisation-list-mps (cl-mpm:sim-mesh *sim*)
;;                                   (cl-mpm:sim-mps *sim*)))


#- :sb-simd
(progn
  (declaim
   (inline diff-squared)
   (ftype (function (cl-mpm/particle:particle cl-mpm/particle:particle) double-float) diff-squared))
  (defun diff-squared (mp-a mp-b)
    (let ((dist (magicl:.- (cl-mpm/particle:mp-position mp-a)
                           (cl-mpm/particle:mp-position mp-b)))
          )
      (values (the double-float (magicl::sum (magicl.simd::.*-simd dist dist)))))))

;;This is a simd dot product
#+ :sb-simd
(progn
  (declaim
   (inline simd-accumulate)
   (ftype (function ((simple-array double-float) (simple-array double-float)) (values double-float)) simd-accumulate))
  ;; (defun simd-accumulate (a b)
  ;;   (declare (type sb-simd:f64vec a b))
  ;;   (let ((diff
  ;;           (sb-simd-avx:f64.2-
  ;;            (sb-simd-avx:f64.2-aref a 0)
  ;;            (sb-simd-avx:f64.2-aref b 0)
  ;;            )))
  ;;     (values (sb-simd-avx::f64.2-horizontal+
  ;;              (sb-simd-avx:f64.2* diff diff))))
  ;;   )
  (defun simd-accumulate (a b)
    (+
     (the double-float (expt (- (aref a 0) (aref b 0)) 2))
     (the double-float (expt (- (aref a 1) (aref b 1)) 2))
     (the double-float (expt (- (aref a 2) (aref b 2)) 2))
     )
    ;; (loop for i fixnum from 0 to 2
    ;;       sum
    ;;       (the double-float (expt (- (aref a i) (aref b i)) 2))
    ;;       double-float
    ;;       )
    )
  (defun diff-squared-mat (pos-a pos-b)
    (let ((pos-a (magicl::matrix/double-float-storage pos-a))
          (pos-b (magicl::matrix/double-float-storage pos-b))
          )
      (values (the double-float (simd-accumulate pos-a pos-b)))))
  (declaim
   (inline diff-squared)
   (ftype (function (cl-mpm/particle:particle cl-mpm/particle:particle) double-float) diff-squared))
  (defun diff-squared (mp-a mp-b)
    (let ((pos-a (magicl::matrix/double-float-storage (cl-mpm/particle:mp-position mp-a)))
          (pos-b (magicl::matrix/double-float-storage (cl-mpm/particle:mp-position mp-b)))
          )
      (values (the double-float (simd-accumulate pos-a pos-b)))))
  )

(declaim
 (inline diff-damaged)
 (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle cl-mpm/particle:particle) double-float) diff-damaged))
(defun diff-damaged (mesh mp-a mp-b)
  (declare (optimize (speed 3)))
  (flet ((node-dam (node)
           (if (cl-mpm/mesh::node-active node)
               (/ (the double-float (cl-mpm/mesh::node-damage node))
                  (the double-float (cl-mpm/mesh::node-svp-sum node))
                  )
               0d0)))
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
        (magicl:.- pos-b pos-a diff)
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
                   (step-size (/ h 1d0))
                   )
              (magicl:scale! step-norm (/ 1d0 length))
              (multiple-value-bind (steps remainder) (floor length step-size)
                (when (> steps 0)
                  (let* ((dhstep (cl-mpm/utils::vector-copy step-norm)))
                    (magicl:scale! dhstep (* 0.5d0 step-size))
                    (loop for i fixnum from 0 below steps
                          do
                             (progn
                               (magicl.simd::.+-simd step-point dhstep step-point)
                               (let ((damage 0d0))
                                 (declare (double-float damage))
                                 (cl-mpm::iterate-over-neighbours-point-linear-3d
                                  mesh
                                  step-point
                                  (lambda (m node weight grads)
                                    (declare (double-float damage weight))
                                    (incf damage
                                          (* (node-dam node)
                                             weight))))
                                 ;; (print damage)
                                 (incf final-distance (* step-size
                                                         (/ 1d0 (max epsilon (sqrt (- 1d0 damage)))))))
                               (magicl.simd::.+-simd step-point dhstep step-point)
                               )
                          )))
                (let* ((step-size remainder)
                       (dhstep (cl-mpm/utils::vector-copy step-norm)))
                  (magicl:scale! dhstep (* 0.5d0 step-size))
                  (progn
                    (magicl.simd::.+-simd step-point dhstep step-point)
                    (let ((damage 0d0))
                      (cl-mpm::iterate-over-neighbours-point-linear-3d
                       mesh
                       step-point
                       (lambda (m node weight grads)
                         (declare (double-float damage weight))
                         (incf damage
                               (* (node-dam node) weight))))
                                        ; (print damage)
                      (incf final-distance (* step-size (/ 1d0
                                                           (max epsilon
                                                                (the double-float (sqrt (- 1d0 damage))))))))
                    )
                  )
                )
              ))
                                         ;(setf final-distance length)
          (* final-distance final-distance))
        ;; (diff-squared mp-a mp-b)
        ))))



(declaim
 (inline weight-func)
 (ftype (function (double-float
                   double-float
                   ) (values double-float))
        weight-func
        ))
(defun weight-func (dist-squared length)
  ;(values (the double-float (exp (the double-float (- (* (/ 4d0 (* length length)) dist-squared))))))
  (if (< dist-squared (* 4d0 length length))
      (values (the double-float (exp (the double-float (* 1d0 (/ (- dist-squared) (* 2d0 length length)))))))
      0d0)
  ;; (values (the double-float (exp (the double-float (* 4d0 (/ (- dist-squared) (* 1.00d0 length length)))))))
  ;; (if (= dist-squared 0d0)
  ;;   1d0
  ;;   0d0)
  )
(declaim
 (inline weight-func-mps)
 (ftype (function (cl-mpm/mesh::mesh
                   cl-mpm/particle:particle
                   cl-mpm/particle:particle
                   double-float
                   ) double-float)
        weight-func-mps
        ))
(defun weight-func-mps (mesh mp-a mp-b length)
  (weight-func (diff-squared mp-a mp-b) length))

(defun weight-func-pos (mesh pos-a pos-b length)
  (let ((ps-a (magicl::matrix/double-float-storage pos-a))
        (ps-b (magicl::matrix/double-float-storage pos-b))
        )
    (weight-func (the double-float (simd-accumulate ps-a ps-b)) length))
  )
(declaim
 (inline weight-func-mps-damaged)
 (ftype (function (cl-mpm/mesh::mesh
                   cl-mpm/particle:particle
                   cl-mpm/particle:particle
                   double-float
                   ) double-float)
        weight-func-mps-damaged
        ))
(defun weight-func-mps-damaged (mesh mp-a mp-b length)
  (weight-func (diff-damaged mesh mp-a mp-b) length))

(defun iterate-over-damage-bounds (mesh mp length func)
  (if (= (cl-mpm/mesh:mesh-nd mesh) 2)
      (iterate-over-damage-bounds-2d mesh mp length func)
      (iterate-over-damage-bounds-3d mesh mp length func)))
(defun iterate-over-damage-bounds-2d (mesh mp length func)
  (let ((node-id (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp)))
        (node-reach (the fixnum (+ 0 (truncate (ceiling (* length 4d0)
                                                        (the double-float (cl-mpm/mesh:mesh-resolution mesh))))))))
    (declare (dynamic-extent node-id))
    (loop for dx fixnum from (- node-reach) to node-reach
          do (loop for dy fixnum from (- node-reach) to node-reach
                   do
                      (let ((idx (mapcar #'+ node-id (list dx dy 0))))
                        (declare (dynamic-extent idx))
                        (when (cl-mpm/mesh:in-bounds mesh idx)
                          (let ((node (cl-mpm/mesh:get-node mesh idx)))
                            (funcall func mesh mp node))))))))
(defun iterate-over-damage-bounds-3d (mesh mp length func)
  (let ((node-id (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp)))
        (node-reach (the fixnum (+ 0 (truncate (ceiling (* length 4d0)
                                                        (the double-float (cl-mpm/mesh:mesh-resolution mesh))))))))
    (declare (dynamic-extent node-id))
    (loop for dx fixnum from (- node-reach) to node-reach
          do (loop for dy fixnum from (- node-reach) to node-reach
                   do
                      (loop for dz fixnum from (- node-reach) to node-reach
                            do
                               (let ((idx (mapcar #'+ node-id (list dx dy dz))))
                                 (declare (dynamic-extent idx))
                                 (when (cl-mpm/mesh:in-bounds mesh idx)
                                   (let ((node (cl-mpm/mesh:get-node mesh idx)))
                                     (funcall func mesh mp node)))))))))

(defparameter *enable-reflect-x* nil)

(declaim
 (notinline calculate-delocalised-damage)
 (ftype (function (cl-mpm/mesh::mesh
                   cl-mpm/particle:particle-damage
                   double-float
                   ) double-float)
        calculate-delocalised-damage
        ))
(defun calculate-delocalised-damage (mesh mp length)
  (let ((damage-inc 0d0)
        (mass-total 0d0))
    (declare (double-float damage-inc mass-total))
    (iterate-over-damage-bounds mesh mp length
                                (lambda (mesh mp node)
                                  (loop for mp-other across (the (vector cl-mpm/particle::particle *) (cl-mpm/mesh::node-local-list node))
                                        do
                                           (with-accessors ((d cl-mpm/particle::mp-damage)
                                                            (m cl-mpm/particle:mp-volume)
                                                            (ll cl-mpm/particle::mp-true-local-length)
                                                            (p cl-mpm/particle:mp-position))
                                               mp-other
                                             (when (< (the double-float d) 1d0)
                                               (let (
                                                     ;;Nodally averaged local funcj
                                                     ;; (weight (weight-func-mps mesh mp mp-other (* 0.5d0 (+ length ll))))
                                                     ;; (weight (weight-func-mps mesh mp mp-other (sqrt (* length ll))))
                                                     ;;
                                                     (weight
                                                       (weight-func-mps mesh mp mp-other (sqrt (* length ll)))
                                                       ;; (if (and (<= d 0) (<= (cl-mpm/particle::mp-damage mp-other) 0))
                                                       ;;     (weight-func-mps mesh mp mp-other (* 0.5d0 (+ length ll)))
                                                       ;; (weight-func-mps-damaged mesh mp mp-other
                                                       ;;                          (cl-mpm/particle::mp-local-length mp)
                                                       ;;                          ;; (* 0.5d0 (+ length ll))
                                                       ;;                          )
                                                       )
                                                     )
                                                 (declare (double-float weight m d mass-total damage-inc))
                                                 (incf mass-total (* weight m))
                                                 (incf damage-inc
                                                       (* (the double-float (cl-mpm/particle::mp-local-damage-increment mp-other))
                                                          weight m)))
                                               (when (and *enable-reflect-x*
                                                      (< (magicl:tref (cl-mpm/particle::mp-position mp) 0 0) (* 4 length)))
                                                 (let ((weight (weight-func-pos mesh
                                                                                (cl-mpm/particle::mp-position mp)
                                                                                (magicl:.* (cl-mpm/particle:mp-position mp-other)
                                                                                           (cl-mpm/utils::vector-from-list (list -1d0 0d0 0d0)))
                                                                                (sqrt (* length ll)))))
                                                 (declare (double-float weight m d mass-total damage-inc))
                                                 (incf mass-total (* weight m))
                                                 (incf damage-inc
                                                       (* (the double-float (cl-mpm/particle::mp-local-damage-increment mp-other))
                                                          weight m)))

                                                 ))
                                             ))))
    (when (> mass-total 0d0)
      (setf damage-inc (/ damage-inc mass-total)))
    damage-inc
    ;; (setf damage-inc (cl-mpm/particle::mp-local-damage-increment mp))
    ))
;; (defun calculate-delocalised-damage (mesh mp length)
;;   (let ((node-id (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp)))
;;         (node-reach (the fixnum (+ 0 (truncate (ceiling (* length 2d0)
;;                                                         (the double-float (cl-mpm/mesh:mesh-resolution mesh))))))))
;;     (declare (dynamic-extent node-id))
;;     (let ((damage-inc 0d0)
;;           (mass-total 0d0))
;;       (declare (double-float damage-inc mass-total))
;;       (loop for dx fixnum from (- node-reach) to node-reach
;;             do (loop for dy fixnum from (- node-reach) to node-reach
;;                      do (loop for dz fixnum from (- node-reach) to node-reach
;;                               do
;;                                  (let ((idx (mapcar #'+ node-id (list dx dy dz))))
;;                                    (declare (dynamic-extent idx))
;;                                    (when (cl-mpm/mesh:in-bounds mesh idx)
;;                                      (let ((node (cl-mpm/mesh:get-node mesh idx)))
;;                                        (loop for mp-other across (the (vector cl-mpm/particle::particle *) (cl-mpm/mesh::node-local-list node))
;;                                              do
;;                                                 (with-accessors ((d cl-mpm/particle::mp-damage)
;;                                                                  (m cl-mpm/particle:mp-volume)
;;                                                                  (ll cl-mpm/particle::mp-true-local-length)
;;                                                                  (p cl-mpm/particle:mp-position))
;;                                                     mp-other
;;                                                   (when (< (the double-float d) 1d0)
;;                                                     (let (
;;                                                           ;;Nodally averaged local funcj
;;                                                           ;; (weight (weight-func-mps mesh mp mp-other (* 0.5d0 (+ length ll))))
;;                                                           (weight (weight-func-mps mesh mp mp-other (sqrt (* length ll))))
;;                                                           ;;
;;                                                           ;; (weight (weight-func-mps-damaged mesh mp mp-other
;;                                                           ;;                                  (cl-mpm/particle::mp-local-length mp)
;;                                                           ;;                                  ;(* 0.5d0 (+ length ll))
;;                                                           ;;                                  ))
;;                                                           )
;;                                                       (declare (double-float weight m d mass-total damage-inc))
;;                                                       (incf mass-total (* weight m))
;;                                                       (incf damage-inc
;;                                                             (* (the double-float (cl-mpm/particle::mp-local-damage-increment mp-other))
;;                                                                weight m))))))))))))
;;       (when (> mass-total 0d0)
;;         (setf damage-inc (/ damage-inc mass-total)))
;;       damage-inc
;;       ;; (setf damage-inc (cl-mpm/particle::mp-local-damage-increment mp))
;;       )))

(defun length-localisation (local-length local-length-damaged damage)
  (+ (* local-length (- 1d0 damage)) (* local-length-damaged damage))
  ;; (* local-length (max (sqrt (- 1d0 damage)) 1d-10))
  ;; local-length
  )
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
                         (local-length-t cl-mpm/particle::mp-true-local-length)
                         )
            mp
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

(defclass mpm-sim-damage (cl-mpm::mpm-sim-usf)
  ()
  (:documentation "Explicit simulation with update stress first update"))
(defclass mpm-sim-damage-nd-2 (mpm-sim-damage cl-mpm::mpm-nd-2d) ())

(defun remove-mp-damage (sim func)
  )

(defmethod cl-mpm::remove-mps-func ((sim mpm-sim-damage) func)
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (mps cl-mpm::sim-mps))
      sim
    (loop for mp across mps
          when (and
                (typep mp 'cl-mpm/particle:particle-damage)
                (funcall func mp))
            do (local-list-remove-particle mesh mp)))
  (call-next-method))

(defmethod (setf cl-mpm::sim-mps) (mps (sim mpm-sim-damage))
  ;; (format t "Resetting mps~%")
  ;; (with-accessors ((mesh cl-mpm::sim-mesh))
  ;;     sim
  ;;   (loop for mp across mps
  ;;         when (and
  ;;               (typep mp 'cl-mpm/particle:particle-damage)
  ;;               (eq (cl-mpm/particle::mp-damage-position mp) nil))
  ;;           do (local-list-add-particle mesh mp)))
  (call-next-method))

(defmethod cl-mpm::update-sim ((sim mpm-sim-damage))
  (declare (cl-mpm::mpm-sim-usf sim))
  (with-slots ((mesh cl-mpm::mesh)
               (mps  cl-mpm::mps)
               (bcs  cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (bcs-force-list cl-mpm::bcs-force-list)
               (dt cl-mpm::dt)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (nonlocal-damage cl-mpm::nonlocal-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (fbar cl-mpm::enable-fbar)
               )
                sim
    (declare (type double-float mass-filter))
                (progn
                    (cl-mpm::reset-grid mesh)
                    (cl-mpm::p2g mesh mps)
                    ;; (cl-mpm::check-single-mps sim)
                    (when (> mass-filter 0d0)
                      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
                      ;; (cl-mpm::filter-grid-volume mesh 1d-8)
                    (cl-mpm::update-node-kinematics mesh dt )
                    (cl-mpm::apply-bcs mesh bcs dt)
                    (cl-mpm::update-stress mesh mps dt fbar)
                    (when enable-damage
                     (cl-mpm/damage::calculate-damage mesh
                                                      mps
                                                      dt
                                                      50d0
                                                      nonlocal-damage
                                                      ))
                    ;Map forces onto nodes
                    ;;13.5
                    (cl-mpm::p2g-force mesh mps)
                    ;;15
                    (loop for bcs-f in bcs-force-list
                          do
                             (cl-mpm::apply-bcs mesh bcs-f dt))
                    ;;15
                    (cl-mpm::update-node-forces mesh (cl-mpm::sim-damping-factor sim) dt (cl-mpm::sim-mass-scale sim))
                    ;Reapply velocity BCs
                    ;;15
                    (cl-mpm::apply-bcs mesh bcs dt)
                    ;Also updates mps inline
                    ;;15
                    (cl-mpm::g2p mesh mps dt)

                    ;;18
                    (when remove-damage
                      (cl-mpm::remove-material-damaged sim))
                    (when split
                      (cl-mpm::split-mps sim))
                    (cl-mpm::check-mps sim)
                    )))


(defstruct damage-model
  name)

(defstruct (damage-model-creep (:include damage-model))
           (initiation-stress
            0d0
            :type DOUBLE-FLOAT)
           (damage-rate
            0d0
            :type DOUBLE-FLOAT)
           )

(defstruct (damage-model-rankine (:include damage-model))
  (initiation-stress
   0d0
   :type DOUBLE-FLOAT)
  (failure-stress
   0d0
   :type DOUBLE-FLOAT)
  (y-history
   0d0
   :type DOUBLE-FLOAT)
  )

(defstruct (damage-model-chalk (:include damage-model))
  (coheasion
   0d0
   :type DOUBLE-FLOAT)
  (friction-angle
   0d0
   :type DOUBLE-FLOAT)
  (damage-rate
   0d0
   :type DOUBLE-FLOAT)
  )


(defgeneric damage-model-calculate-y (mp dt))

(defmethod damage-model-calculate-y ((mp cl-mpm/particle::particle) dt)
  0d0)
(defmethod damage-model-calculate-y ((mp cl-mpm/particle::particle-damage) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (pressure cl-mpm/particle::mp-pressure)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     ) mp
      (declare (double-float pressure damage))
        (progn
          (when (< damage 1d0)
            (let ((cauchy-undamaged (magicl:scale stress (/ 1d0 (* 1d0;(- 1d0 damage)
                                                                   (magicl:det def))))))
              (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d cauchy-undamaged)
                (let* (;;Only allow tensile damage
                       (pressure-effective (* 1d0 damage pressure))
                       ;; (pressure-effective (* 1d0 pressure))
                       (j2 (sqrt (cl-mpm/constitutive::voigt-j2 (cl-mpm/constitutive::deviatoric-voigt
                                                                 cauchy-undamaged))))
                       (p (* 0.3d0 (+ s_1 s_2 s_3)))
                       (s_1 (- s_1 pressure-effective))
                       (s_2 (- s_2 pressure-effective))
                       (s_3 (- s_3 pressure-effective))
                       (s_1 (max 0d0 s_1))
                       (s_2 (max 0d0 s_2))
                       (s_3 (max 0d0 s_3))
                       ;; (s_1 (- j2 (* 0.2d0 p)))
                       )
                  (when (> s_1 0d0)
                    (setf damage-increment s_1)
                    )))))
          (when (>= damage 1d0)
            (setf damage-increment 0d0))
          ;;Delocalisation switch
          (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
          ))))


;; (defmethod damage-model-calculate-y (mp (dm damage-model-rankine) dt)
;;   (let ((damage-increment 0d0))
;;     (with-accessors ((stress cl-mpm/particle::mp-stress)
;;                      (damage cl-mpm/particle:mp-damage)
;;                      (pressure cl-mpm/particle::mp-pressure)
;;                      (ybar cl-mpm/particle::mp-local-damage-increment)
;;                      ) mp
;;       (declare (double-float pressure damage))
;;       (progn
;;         (progn
;;           (multiple-value-bind (s_1 s_2) (principal-stresses stress)
;;             (let* ((pressure-effective (* pressure 1d0))
;;                    (s_1 (- s_1 pressure-effective))
;;                    (s_2 (- s_2 pressure-effective))
;;                    (s_1 (max 0d0 s_1))
;;                    (s_2 (max 0d0 s_2))
;;                    ;; (vm (- s_1 s_2))
;;                                         ;(s_1 vm)
;;                    )
;;               (when (> s_1 0d0)
;;                 (if (< damage 1d0)
;;                     (setf damage-increment (/ s_1 (expt (- 1d0 damage) 2d0)))
;;                     (setf damage-increment s_1)))))
;;           (when (>= damage 1d0)
;;             (setf damage-increment 0d0))
;;           ;; damage-increment
;;           (setf ybar damage-increment)
;;           ;; (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
;;           )))))

;; (defmethod damage-model-calculate-y (mp (dm damage-model-creep) dt)
;;   (let ((damage-increment 0d0))
;;     (with-accessors ((stress cl-mpm/particle::mp-stress)
;;                      (damage cl-mpm/particle:mp-damage)
;;                      (pressure cl-mpm/particle::mp-pressure)
;;                      (ybar cl-mpm/particle::mp-local-damage-increment)
;;                      ) mp
;;       (with-accessors
;;             ((init-stress damage-model-creep-initiation-stress)
;;              (damage-rate damage-model-creep-damage-rate)
;;              )
;;           dm
;;         (declare (double-float pressure damage))
;;         (progn
;;           (progn
;;             (multiple-value-bind (s_1 s_2) (principal-stresses stress)
;;               (let* ((pressure-effective (* pressure 1d0))
;;                      ;; (s_1 (- s_1 pressure-effective))
;;                      ;; (s_2 (- s_2 pressure-effective))
;;                      (s_1 (max 0d0 s_1))
;;                      (s_2 (max 0d0 s_2))
;;                      ;; (vm (- s_1 s_2))
;;                                         ;(s_1 vm)
;;                      )
;;                 (when (> s_1 0d0)
;;                   (if (< damage 1d0)
;;                       (setf damage-increment (/ s_1 (expt (- 1d0 damage) 2d0)))
;;                       (setf damage-increment s_1)))))
;;             (when (>= damage 1d0)
;;               (setf damage-increment 0d0))
;;             ;; damage-increment
;;             (setf ybar damage-increment)
;;             ;; (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
;;             ))))))


(defgeneric damage-model-update-damage (mp damage-model dt))

(defmethod damage-model-update-damage (mp damage-model dt))

(defmethod damage-model-update-damage (mp (dm damage-model-creep) dt)
  (with-accessors ((stress cl-mpm/particle:mp-stress)
                   (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                   (damage cl-mpm/particle:mp-damage)
                   (log-damage cl-mpm/particle::mp-log-damage)
                   (damage-inc cl-mpm/particle::mp-damage-increment)
                   (ybar cl-mpm/particle::mp-damage-ybar)
                   (critical-damage cl-mpm/particle::mp-critical-damage)
                   (pressure cl-mpm/particle::mp-pressure)
                   ) mp
    (with-accessors ((damage-rate damage-model-creep-damage-rate)
                     (init-stress damage-model-creep-initiation-stress))
        dm
      (declare (double-float damage damage-inc dt critical-damage))
      (progn
        ;;Damage increment holds the delocalised driving factor
        (setf ybar damage-inc)
        (setf damage-inc (* dt (damage-rate-profile damage-inc damage damage-rate init-stress)))
        (when (>= damage 1d0)
          (setf damage-inc 0d0)
          (setf ybar 0d0))
        (incf damage damage-inc)
        (setf damage (max 0d0 (min 1d0 damage)))
        (when (> damage critical-damage)
          (setf damage 1d0)
          (setf damage-inc 0d0))))
  (values)))

(defmethod damage-model-update-damage (mp (dm damage-model-rankine) dt)
  (with-accessors ((stress cl-mpm/particle:mp-stress)
                   (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                   (damage cl-mpm/particle:mp-damage)
                   (log-damage cl-mpm/particle::mp-log-damage)
                   (damage-inc cl-mpm/particle::mp-damage-increment)
                   (ybar cl-mpm/particle::mp-damage-ybar)
                   (critical-damage cl-mpm/particle::mp-critical-damage)
                   (pressure cl-mpm/particle::mp-pressure)
                   ) mp
    (with-accessors (;(damage-rate damage-model-creep-damage-rate)
                     (init-stress damage-model-rankine-initiation-stress)
                     (crit-stress damage-model-rankine-failure-stress)
                     (y-history damage-model-rankine-y-history)
                     )
        dm
      (declare (double-float damage damage-inc dt critical-damage))
      (progn
        ;;Damage increment holds the delocalised driving factor
        (setf ybar damage-inc)
        (setf y-history (max y-history ybar))
        (setf damage (/ (- y-history init-stress) crit-stress))
        (when (>= damage 1d0)
          (setf damage-inc 0d0)
          (setf ybar 0d0))

        (incf damage damage-inc)
        (setf damage (max 0d0 (min 1d0 damage)))
        (when (> damage critical-damage)
          (setf damage 1d0)
          (setf damage-inc 0d0))))
  (values)))


(defmethod cl-mpm/output::save-vtk (filename (sim cl-mpm/damage::mpm-sim-damage))
  (with-accessors ((mps cl-mpm:sim-mps)) sim
    (with-open-file (fs filename :direction :output :if-exists :supersede)
      (format fs "# vtk DataFile Version 2.0~%")
      (format fs "Lisp generated vtk file, SJVS~%")
      (format fs "ASCII~%")
      (format fs "DATASET UNSTRUCTURED_GRID~%")
      (format fs "POINTS ~d double~%" (length mps))
      (loop for mp across mps
            do (format fs "~E ~E ~E ~%"
                       (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 1 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 2 0) 'single-float)
                       ))
      (format fs "~%")

      ;; (cl-mpm/output::with-parameter-list fs mps
      ;;   ("mass" 'cl-mpm/particle:mp-mass)
      ;;   ("density" (lambda (mp) (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp))))
      ;;   )
      (let ((id 1))
        (declare (special id))
        (format fs "POINT_DATA ~d~%" (length mps))

        (cl-mpm/output::save-parameter "mass" (cl-mpm/particle:mp-mass mp))
        (cl-mpm/output::save-parameter "density" (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)))
        (cl-mpm/output::save-parameter "index" (cl-mpm/particle::mp-index mp))
        (cl-mpm/output::save-parameter "mpi-index" (cl-mpm/particle::mp-mpi-index mp))
        (cl-mpm/output::save-parameter "vel_x" (magicl:tref (cl-mpm/particle:mp-velocity mp) 0 0))
        (cl-mpm/output::save-parameter "vel_y" (magicl:tref (cl-mpm/particle:mp-velocity mp) 1 0))
        ;; (cl-mpm/output::save-parameter "acc_x" (magicl:tref (cl-mpm/particle::mp-acceleration mp) 0 0))
        ;; (cl-mpm/output::save-parameter "acc_y" (magicl:tref (cl-mpm/particle::mp-acceleration mp) 1 0))

        (cl-mpm/output::save-parameter "disp_x" (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0))
        (cl-mpm/output::save-parameter "disp_y" (magicl:tref (cl-mpm/particle::mp-displacement mp) 1 0))

        (cl-mpm/output::save-parameter "sig_xx" (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0))
        (cl-mpm/output::save-parameter "sig_yy" (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0))
        (cl-mpm/output::save-parameter "sig_zz" (magicl:tref (cl-mpm/particle:mp-stress mp) 2 0))
        (cl-mpm/output::save-parameter "sig_yz" (magicl:tref (cl-mpm/particle:mp-stress mp) 3 0))
        (cl-mpm/output::save-parameter "sig_zx" (magicl:tref (cl-mpm/particle:mp-stress mp) 4 0))
        (cl-mpm/output::save-parameter "sig_xy" (magicl:tref (cl-mpm/particle:mp-stress mp) 5 0))

        (cl-mpm/output::save-parameter "e_xx" (magicl:tref (cl-mpm/particle::mp-strain mp) 0 0))
        (cl-mpm/output::save-parameter "e_yy" (magicl:tref (cl-mpm/particle::mp-strain mp) 1 0))
        (cl-mpm/output::save-parameter "e_zz" (magicl:tref (cl-mpm/particle::mp-strain mp) 2 0))

        ;; (cl-mpm/output::save-parameter "temp" (magicl:tref (cl-mpm/particle::mp-velocity-rate mp) 2 0))

        (cl-mpm/output::save-parameter "damage-inc-average"
                                       (if (slot-exists-p mp 'cl-mpm/particle::time-averaged-damage-inc)
                                           (let ((v (/ (cl-mpm/particle::mp-time-averaged-damage-inc mp)
                                                       (max 1d0
                                                            (cl-mpm/particle::mp-time-averaged-counter mp)))))
                                             (setf (cl-mpm/particle::mp-time-averaged-damage-inc mp) 0d0)
                                             v)
                                           0d0
                                           ))
        (cl-mpm/output::save-parameter "damage-ybar-average"
                                       (if (slot-exists-p mp 'cl-mpm/particle::time-averaged-damage-inc)
                                           (let ((v (/ (cl-mpm/particle::mp-time-averaged-ybar mp)
                                                       (max 1d0
                                                            (cl-mpm/particle::mp-time-averaged-counter mp)))))
                                             (setf (cl-mpm/particle::mp-time-averaged-counter mp) 0d0
                                                   (cl-mpm/particle::mp-time-averaged-ybar mp) 0d0)
                                             v)
                                           0))
        ;; (cl-mpm/output::save-parameter "viscosity" (cl-mpm/particle::mp-time-averaged-visc mp))
        ;; (save-parameter "visc-plastic" (cl-mpm/particle::mp-visc-plastic mp))
        ;; (save-parameter "visc-glen" (cl-mpm/particle::mp-visc-glen mp))

        (cl-mpm/output::save-parameter "strain-rate"
                                       (cl-mpm/constitutive::effective-strain-rate (cl-mpm/particle::mp-eng-strain-rate mp))
                                       ;; (multiple-value-bind (l v)
                                       ;;     (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle::mp-velocity-rate mp)))
                                       ;;   (reduce #'+ (mapcar #'* l l)))
                                       )
        (cl-mpm/output::save-parameter "pressure" (cl-mpm/particle::mp-pressure mp))

        (cl-mpm/output::save-parameter "s_1"
                                       (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))
                                         (loop for sii in l maximize sii)))
        (cl-mpm/output::save-parameter "su_1"
                                       (if (slot-exists-p mp 'cl-mpm/particle::undamaged-stress)
                                           (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle::mp-undamaged-stress mp)))
                                             (loop for sii in l maximize sii))
                                           0d0
                                           ))
        (cl-mpm/output::save-parameter "s_vm"
                                       (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))

                                         (let* ((l (sort l #'>))
                                                (s_1 (max 0 (- (first l) (cl-mpm/particle::mp-pressure mp))))
                                                (s_2 (max 0 (- (second l) (cl-mpm/particle::mp-pressure mp)))))

                                           (* (sqrt (/ 3 4)) (- s_1 s_2))
                                           )
                                         ))
        (cl-mpm/output::save-parameter "s_vm_t"
                                       (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))

                                         (let* ((l (sort l #'>))
                                                (s_1 (first l))
                                                (s_2 (second l)))
                                           (* ;(sqrt (/ 3 4))
                                            1d0
                                            (- s_1 s_2))
                                           )
                                         ))


        (cl-mpm/output::save-parameter "EPS"
                                       (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))
                                         (- (loop for sii in l maximize sii) (cl-mpm/particle::mp-pressure mp))))
        (cl-mpm/output::save-parameter "EPS-pd"
                                       (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))
                                         (- (loop for sii in l maximize sii) (* (cl-mpm/particle::mp-damage mp)
                                                                                (cl-mpm/particle::mp-pressure mp)))))
        (cl-mpm/output::save-parameter "size_x" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 0 0))
        (cl-mpm/output::save-parameter "size_y" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0))

        (cl-mpm/output::save-parameter "damage"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage)
                                           (cl-mpm/particle:mp-damage mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-inc"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-increment)
                                           (cl-mpm/particle::mp-damage-increment mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-ybar"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-ybar)
                                           (cl-mpm/particle::mp-damage-ybar mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-y"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-y-local)
                                           (cl-mpm/particle::mp-damage-y-local mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-xx"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-tensor)
                                           (magicl:tref (cl-mpm/particle::mp-damage-tensor mp) 0 0)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-yy"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-tensor)
                                           (magicl:tref (cl-mpm/particle::mp-damage-tensor mp) 1 1)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-y-xx"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-tensor)
                                           (magicl:tref (cl-mpm/particle::mp-damage-ybar-tensor mp) 0 0)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-y-yy"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-tensor)
                                           (magicl:tref (cl-mpm/particle::mp-damage-ybar-tensor mp) 1 1)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-k"
                                       (if (slot-exists-p mp 'cl-mpm/particle::history-stress)
                                           (cl-mpm/particle::mp-history-stress mp)
                                           0d0))
        (cl-mpm/output::save-parameter "local-length"
                                       (if (slot-exists-p mp 'cl-mpm/particle::true-local-length)
                                           (cl-mpm/particle::mp-true-local-length mp)
                                           0d0))
        (cl-mpm/output::save-parameter "j2"
                                       (with-accessors ((damage cl-mpm/particle::mp-damage)
                                                        (stress cl-mpm/particle::mp-undamaged-stress)
                                                        (def cl-mpm/particle::mp-deformation-gradient)
                                                        )
                                           mp
                                         (if (< damage 1d0)
                                             (sqrt (cl-mpm/constitutive::voigt-j2
                                                    (cl-mpm/utils::deviatoric-voigt stress)))
                                             0d0)))
        (cl-mpm/output::save-parameter "i1"
                                       (with-accessors ((damage cl-mpm/particle::mp-damage)
                                                        (stress cl-mpm/particle::mp-undamaged-stress)
                                                        (def cl-mpm/particle::mp-deformation-gradient)
                                                        )
                                           mp
                                         (if (< damage 1d0)
                                             (cl-mpm/utils::trace-voigt stress)
                                             0d0)))
        (cl-mpm/output::save-parameter "split-depth"
                                       (cl-mpm/particle::mp-split-depth mp))
        (cl-mpm/output::save-parameter
         "plastic_strain"
         (if (slot-exists-p mp 'cl-mpm/particle::strain-plastic-vm)
             (cl-mpm/particle::mp-strain-plastic-vm mp)
             0d0)
         )
        (cl-mpm/output::save-parameter
         "yield-func"
         (if (slot-exists-p mp 'cl-mpm/particle::yield-func)
             (cl-mpm/particle::mp-yield-func mp)
             0d0))
        (cl-mpm/output::save-parameter
         "rho"
         (if (slot-exists-p mp 'cl-mpm/particle::rho)
             (cl-mpm/particle::mp-rho mp)
             0d0))
        )
      )))

(defgeneric update-damage (mp dt))
(defmethod update-damage ((mp cl-mpm/particle::particle) dt))

(defmethod update-damage ((mp cl-mpm/particle::particle-damage) dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (log-damage cl-mpm/particle::mp-log-damage)
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
          ;; (setf damage-inc (* dt (- critical-damage damage) (damage-rate-profile damage-inc damage damage-rate init-stress)))
          (setf damage-inc (* dt (damage-rate-profile damage-inc damage damage-rate init-stress)))
          (when (>= damage 1d0)
            (setf damage-inc 0d0)
            (setf ybar 0d0))
          (incf (cl-mpm/particle::mp-time-averaged-damage-inc mp) damage-inc)
          (incf (cl-mpm/particle::mp-time-averaged-ybar mp) ybar)
          (incf (cl-mpm/particle::mp-time-averaged-counter mp))
          (incf damage damage-inc)
          (setf damage (max 0d0 (min 1d0 damage)))
          (when (> damage critical-damage)
            (setf damage 1d0)
            (setf damage-inc 0d0)))
  (values)
  ))


(defmethod update-damage ((mp cl-mpm/particle::particle-chalk) dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (log-damage cl-mpm/particle::mp-log-damage)
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
          (when (< damage 1d0)
            (setf damage-inc (* dt
                                ;; (/ 1d0 (- 1d0 damage))
                                (damage-rate-profile-chalk damage-inc damage damage-rate init-stress))))
          (when (>= damage 1d0)
            (setf damage-inc 0d0)
            (setf ybar 0d0))
          (incf (cl-mpm/particle::mp-time-averaged-damage-inc mp) damage-inc)
          (incf (cl-mpm/particle::mp-time-averaged-ybar mp) ybar)
          (incf (cl-mpm/particle::mp-time-averaged-counter mp))
          ;;Transform to log damage
          (incf damage damage-inc)
          ;;Transform to linear damage
          (setf damage (max 0d0 (min 1d0 damage)))
          (when (> damage critical-damage)
            (setf damage 1d0)
            (setf damage-inc 0d0)))
  (values)
  ))

(defmethod damage-model-calculate-y ((mp cl-mpm/particle::particle-chalk) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (pressure cl-mpm/particle::mp-pressure)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     (angle cl-mpm/particle::mp-friction-angle)
                     (c cl-mpm/particle::mp-coheasion)
                     ) mp
      (declare (double-float pressure damage))
        (progn
          (when (< damage 1d0)
            (let ((cauchy-undamaged (magicl:scale stress (/ 1d0 (magicl:det def)))))
              (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d cauchy-undamaged)
                (let* ((pressure-effective (* 1d0 damage pressure))
                       (j2 (sqrt (cl-mpm/constitutive::voigt-j2
                                  (cl-mpm/utils::deviatoric-voigt cauchy-undamaged))))
                       (p (+ s_1 s_2 s_3))
                       (s_1 (- s_1 pressure-effective))
                       (s_2 (- s_2 pressure-effective))
                       (s_3 (- s_3 pressure-effective))
                       (s_1 (max 0d0 s_1))
                       (s_2 (max 0d0 s_2))
                       (s_3 (max 0d0 s_3))
                       (angle (* angle (/ pi 180d0)))
                       ;; (c 1d4)
                       (A (/ (* 6 c (cos angle))
                             (* (sqrt 3) (- 3 (sin angle)))))
                       (B (/ (* 2d0 (sin angle)) (* (sqrt 3) (- 3 (sin angle)))))
                       ;; (s_1 (+ j2 (- (* B p) A)))
                       ;; (s_1 j2)
                       )
                  (when (> s_1 0d0)
                    ;; (setf damage-increment (* s_1 (expt (- 1d0 damage) -2d0)))
                    (setf damage-increment s_1)
                    )))))
          (when (>= damage 1d0)
            (setf damage-increment 0d0))
          ;;Delocalisation switch
          (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
          ))))

(defun damage-rate-profile-chalk (stress damage rate init-stress)
  (declare (double-float stress damage rate init-stress))
  "Function that controls how damage evolves with principal stresses"
  (if (> stress init-stress)
      (* (expt (max 0d0 (- stress init-stress)) 0.5d0) rate)
      0d0))


(defmethod damage-model-calculate-y ((mp cl-mpm/particle::particle-chalk-creep) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (pressure cl-mpm/particle::mp-pressure)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     (angle cl-mpm/particle::mp-friction-angle)
                     (c cl-mpm/particle::mp-coheasion)
                     ) mp
      (declare (double-float pressure damage))
        (progn
          (when (< damage 1d0)
            (let ((cauchy-undamaged (magicl:scale stress (/ 1d0 (magicl:det def)))))
              (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d cauchy-undamaged)
                (let* ((pressure-effective (* 1d0 damage pressure))
                       (j2 (sqrt (cl-mpm/constitutive::voigt-j2
                                  (cl-mpm/utils::deviatoric-voigt cauchy-undamaged))))
                       (p (+ s_1 s_2 s_3))
                       (s_1 (- s_1 pressure-effective))
                       (s_2 (- s_2 pressure-effective))
                       (s_3 (- s_3 pressure-effective))
                       (s_1 (max 0d0 s_1))
                       (s_2 (max 0d0 s_2))
                       (s_3 (max 0d0 s_3))
                       (angle (* angle (/ pi 180d0)))
                       ;; (c 1d4)
                       (A (/ (* 6 c (cos angle))
                             (* (sqrt 3) (- 3 (sin angle)))))
                       (B (/ (* 2d0 (sin angle)) (* (sqrt 3) (- 3 (sin angle)))))
                       ;; (s_1 (+ j2 (- (* B p) A)))
                       ;; (s_1 j2)
                       )
                  (when (> s_1 0d0)
                    ;; (setf damage-increment (* s_1 (expt (- 1d0 damage) -2d0)))
                    (setf damage-increment s_1)
                    )))))
          (when (>= damage 1d0)
            (setf damage-increment 0d0))
          ;;Delocalisation switch
          (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
          ))))

(defmethod update-damage ((mp cl-mpm/particle::particle-chalk-creep) dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (log-damage cl-mpm/particle::mp-log-damage)
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
          (when (< damage 1d0)
            (setf damage-inc (* dt
                                ;; (/ 1d0 (- 1d0 damage))
                                (damage-rate-profile-chalk damage-inc damage damage-rate init-stress))))
          (when (>= damage 1d0)
            (setf damage-inc 0d0)
            (setf ybar 0d0))
          (incf (cl-mpm/particle::mp-time-averaged-damage-inc mp) damage-inc)
          (incf (cl-mpm/particle::mp-time-averaged-ybar mp) ybar)
          (incf (cl-mpm/particle::mp-time-averaged-counter mp))
          ;;Transform to log damage
          (incf damage damage-inc)
          ;;Transform to linear damage
          (setf damage (max 0d0 (min 1d0 damage)))
          (when (> damage critical-damage)
            (setf damage 1d0)
            (setf damage-inc 0d0)))
  (values)
  ))

;;Chalk brittle

(defun drucker-prager-criterion (stress angle)
  (let ((p (cl-mpm/utils::trace-voigt stress))
        (j2 (cl-mpm/constitutive::voigt-j2
             (cl-mpm/utils::deviatoric-voigt stress))))
    (* (/ 3d0 (+ 3d0 (tan angle))) (+ (sqrt (* 3 j2)) (* 1/3 (tan angle) p))))
  )
(defun modified-vm-criterion (stress nu k)
  (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d stress)
    (let* ((j2 (cl-mpm/constitutive::voigt-j2
                (cl-mpm/utils::deviatoric-voigt stress)))
           (i1 (+ s_1 s_2 s_3))
           (k-factor (/ (- k 1d0)
                        (- 1d0 (* 2d0 nu))))
           (s_1 (+ (* i1 (/ k-factor (* 2d0 k)))
                   (* (/ 1d0 (* 2d0 k))
                      (sqrt (+ (expt (* k-factor i1) 2)
                               (* (/ (* 12 k) (expt (- 1d0 nu) 2))j2)
                               )))))
           )
      s_1
      )))

(defmethod damage-model-calculate-y ((mp cl-mpm/particle::particle-chalk-brittle) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (pressure cl-mpm/particle::mp-pressure)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     (angle cl-mpm/particle::mp-friction-angle)
                     (c cl-mpm/particle::mp-coheasion)
                     (nu cl-mpm/particle::mp-nu)
                     ) mp
      (declare (double-float pressure damage))
        (progn
          (when (< damage 1d0)
            (let ((cauchy-undamaged (magicl:scale stress (/ 1d0 (magicl:det def)))))
              (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d cauchy-undamaged)
                (let* ((pressure-effective (* 1d0 damage pressure))
                       (j2 (cl-mpm/constitutive::voigt-j2
                                   (cl-mpm/utils::deviatoric-voigt cauchy-undamaged)))
                       (p (+ s_1 s_2 s_3))
                       (s_1 (- s_1 pressure-effective))
                       (s_2 (- s_2 pressure-effective))
                       (s_3 (- s_3 pressure-effective))
                       ;; (s_1 (max 0d0 s_1))
                       ;; (s_2 (max 0d0 s_2))
                       ;; (s_3 (max 0d0 s_3))
                       ;; (s_1 (sqrt
                       ;;       (+ (expt s_1 2)
                       ;;          (expt s_2 2)
                       ;;          (expt s_3 2))))
                       ;; (angle (* angle (/ pi 180d0)))
                       (ft 300d3)
                       (fc 500d3)
                       (angle (atan (* 3 (/ (- fc ft) (+ fc ft)))))
                       (k (* (/ 2d0 (sqrt 3)) (/ (* ft fc) (+ ft fc))))
                       ;; ;; (c 1d4)
                       ;;Another dp
                       ;; (alpha (/ (- (/ fc ft) 1) (+ (/ fc ft) 1)))
                       ;; (s_1 (/ (+ (sqrt (* 3 j2)) (* alpha p)) (+ 1 alpha)))

                       ;; (A (/ (* 6 c (cos angle))
                       ;;       (* (sqrt 3) (- 3 (sin angle)))))
                       ;; (B (/ (* 2d0 (sin angle)) (* (sqrt 3) (- 3 (sin angle)))))
                       ;; (s_1 (sqrt
                       ;;       (+ (expt s_1 2)
                       ;;          (expt s_2 2)
                       ;;          (expt s_3 2))))
                       (s_1 (* (/ 3d0 (+ 3 (tan angle)))
                               (+ (sqrt (* 3 j2)) (* 1/3 (tan angle) p))))
                       ;; (s_1 (+ (sqrt j2) (- (* B p) A)))
                       ;; (s_1 j2)

                       ;; (k (/ fc ft))
                       ;; (i1 (+ s_1 s_2 s_3))
                       ;; (k-factor (/ (- k 1d0)
                       ;;              (- 1d0 (* 2d0 nu))))
                       ;; (s_1 (+ (* i1 (/ k-factor (* 2d0 k)))
                       ;;         (* (/ 1d0 (* 2d0 k))
                       ;;            (sqrt (+ (expt (* k-factor i1) 2)
                       ;;                     (* (/ (* 12 k) (expt (- 1d0 nu) 2))j2)
                       ;;                     )))))
                       )
                  ;; (setf damage-increment s_1)
                  ;; (when (> s_1 0d0)
                    ;; (setf damage-increment (* s_1 (expt (- 1d0 damage) -2d0)))
                    (setf damage-increment s_1)
                    ;; )
                  ))))
          (when (>= damage 1d0)
            (setf damage-increment 0d0))
          ;;Delocalisation switch
          (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
          (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
          ))))

(defun brittle-chalk-d (stress E Gf length init-stress)
  (declare (double-float stress E Gf length init-stress))
  "Function that controls how damage evolves with principal stresses"
  ;; (if (> stress init-stress)
  ;;     (* (expt (max 0d0 (- stress init-stress)) 0.5d0) rate) 0d0)
  ;; (when (> length (/ (* 2 Gf E) (expt init-stress 2)))
  ;;   (error "Length scale is too long"))
  (let* ((hs (/ (expt init-stress 2) (* 2d0 E Gf)))
         (hsl (/ (* hs length) (- 1d0 (* hs length)))))
    (if (> stress init-stress)
        (min 1d0 (max 0d0
                      ;; (* (+ 1d0 hsl)
                      ;;    (- 1d0 (/ init-stress stress)
                      ;;       ))
                      (- 1d0 (* (/ init-stress stress) (exp (* -2d0 hs (/ (- stress init-stress) stress)))))
                      ;; ()
                      ;; (exp (* -2d0 hs (/ (abs (- stress init-stress)) init-stress)))
                      ))
        0d0)
    ))

(defun delayed-damage-y (stress E Gf length init-stress)
  (declare (double-float stress E Gf length init-stress))
  "Function that controls how damage evolves with principal stresses"
  (let* ((ft init-stress)
         (e0 (/ ft E))
         (ef (+ (/ e0 2) (/ Gf (* length ft))))
         (k (/ stress E)))
    (if (> k e0)
        (/ (- stress init-stress) (- (* ef E) init-stress))
        0d0
        ))
  ;; (let* ((hs (/ (expt stress 2) (* 2 E Gf)))
  ;;        (hsl (/ (* hs length) (- 1d0 (* hs length)))))
  ;;   (min 1d0 (max 0d0 (- 1d0 (exp (* -1d0 hs (/ init-stress (max init-stress stress)))))))
  ;;   )
  )
(defmethod update-damage ((mp cl-mpm/particle::particle-chalk-brittle) dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (E cl-mpm/particle::mp-e)
                     (Gf cl-mpm/particle::mp-Gf)
                     (log-damage cl-mpm/particle::mp-log-damage)
                     (damage-inc cl-mpm/particle::mp-damage-increment)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (pressure cl-mpm/particle::mp-pressure)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     ;; (length-t cl-mpm/particle::mp-true-local-length)
                     (length cl-mpm/particle::mp-local-length)
                     (k cl-mpm/particle::mp-history-stress)
                     (tau cl-mpm/particle::mp-delay-time)
                     ) mp
      (declare (double-float damage damage-inc critical-damage))
        (progn
          ;;Damage increment holds the delocalised driving factor
          (setf ybar damage-inc)
          (setf k (max k ybar))
          (setf damage-inc 0d0)
          (let ((new-damage (max damage
                                 (damage-response-exponential k E Gf length init-stress)
                                 )))
            (setf damage-inc (- new-damage damage)))
          (when (>= damage 1d0)
            (setf damage-inc 0d0)
            (setf ybar 0d0))
          (incf (cl-mpm/particle::mp-time-averaged-damage-inc mp) damage-inc)
          (incf (cl-mpm/particle::mp-time-averaged-ybar mp) ybar)
          (incf (cl-mpm/particle::mp-time-averaged-counter mp))
          ;;Transform to log damage
          (incf damage damage-inc)
          ;;Transform to linear damage
          (setf damage (max 0d0 (min 1d0 damage)))
          (when (> damage critical-damage)
            (setf damage 1d0)
            (setf damage-inc 0d0)))
  (values)
  ))
(defmethod update-damage ((mp cl-mpm/particle::particle-chalk-delayed) dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (E cl-mpm/particle::mp-e)
                     (Gf cl-mpm/particle::mp-Gf)
                     (log-damage cl-mpm/particle::mp-log-damage)
                     (damage-inc cl-mpm/particle::mp-damage-increment)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (pressure cl-mpm/particle::mp-pressure)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     ;; (length-t cl-mpm/particle::mp-true-local-length)
                     (length cl-mpm/particle::mp-local-length)
                     (k cl-mpm/particle::mp-history-stress)
                     (tau cl-mpm/particle::mp-delay-time)
                     ) mp
      (declare (double-float damage damage-inc critical-damage))
        (progn
          ;;Damage increment holds the delocalised driving factor
          (setf ybar damage-inc)
          ;; (setf k (max k ybar))
          (setf damage-inc 0d0)

          ;; (incf k (* dt (/ (max 0d0 (- ybar k)) tau)))

          (let ((new-damage (max damage
                                 ;; (test-damage k 1d7 init-stress)
                                 ;; (delayed-damage-y ybar E Gf length init-stress)
                                 (damage-response-exponential ybar E Gf length init-stress)
                                 )))
            ;; (setf damage-inc (- new-damage damage))
            (setf damage-inc (* (/ dt tau) (- 1d0 (exp (- (* 1d0 (abs (- new-damage damage))))))))
            )
          ;; (setf damage (max damage (brittle-chalk-d k E Gf length init-stress))
          ;;       damage-inc 0d0)
          (when (>= damage 1d0)
            (setf damage-inc 0d0)
            (setf ybar 0d0))
          (incf (cl-mpm/particle::mp-time-averaged-damage-inc mp) damage-inc)
          (incf (cl-mpm/particle::mp-time-averaged-ybar mp) ybar)
          (incf (cl-mpm/particle::mp-time-averaged-counter mp))
          ;;Transform to log damage
          (incf damage damage-inc)
          ;;Transform to linear damage
          (setf damage (max 0d0 (min 1d0 damage)))
          (when (> damage critical-damage)
            (setf damage 1d0)
            (setf damage-inc 0d0)))
  (values)
  ))
(defmethod update-damage ((mp cl-mpm/particle::particle-chalk-anisotropic) dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (damage-tensor cl-mpm/particle::mp-damage-tensor)
                     (ybar-tensor cl-mpm/particle::mp-damage-ybar-tensor)
                     (E cl-mpm/particle::mp-e)
                     (Gf cl-mpm/particle::mp-Gf)
                     (log-damage cl-mpm/particle::mp-log-damage)
                     (damage-inc cl-mpm/particle::mp-damage-increment)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (pressure cl-mpm/particle::mp-pressure)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     ;; (length-t cl-mpm/particle::mp-true-local-length)
                     (length cl-mpm/particle::mp-local-length)
                     (k cl-mpm/particle::mp-history-stress)
                     (tau cl-mpm/particle::mp-delay-time)
                     ) mp
      (declare (double-float damage damage-inc critical-damage))
        (progn
          (setf ybar damage-inc)
          (magicl:scale! ybar-tensor 0d0)
          (when (< damage 1d0)
            (let ((damage-inc-mat (cl-mpm/utils:matrix-zeros))
                  (cauchy-undamaged (magicl:scale stress (/ 1d0 (magicl:det def)))))
              (multiple-value-bind (ls v) (cl-mpm/utils:eig (cl-mpm/utils:voigt-to-matrix cauchy-undamaged))
                (loop for i from 0 to 2
                      do
                         (let* ((sii (nth i ls))
                                (vii (magicl::column v i))
                                (vsi (magicl:@ vii (magicl:transpose vii)))
                                (dii (magicl::trace (magicl:@ damage-tensor vsi)))
                                )
                           (when (< sii 0d0)
                             (setf sii 0d0))
                           (when (> sii 0d0)
                             (setf sii damage-inc))
                           (let* ((new-damage (damage-response-exponential sii E Gf length init-stress))
                                  (damage-increment ;(- (max dii new-damage) dii)
                                                    (* (/ dt tau) (- 1d0 (exp (- (* 1d0 (abs (- new-damage dii))))))))
                                  )
                             (magicl:.+ ybar-tensor
                                        (magicl:scale vsi sii)
                                        ybar-tensor)
                             (magicl:.+ damage-inc-mat
                                        (magicl:scale vsi damage-increment)
                                        damage-inc-mat)))))
              ;; (break)
              (magicl:.+ damage-tensor
                         damage-inc-mat
                         damage-tensor)))
          ;;Damage increment holds the delocalised driving factor
          ;; (setf ybar damage-inc)
          ;; (setf k (max k ybar))
          ;; (setf damage-inc 0d0)

          ;; (incf k (* dt (/ (max 0d0 (- ybar k)) tau)))

          ;; (let ((new-damage (max damage
          ;;                        ;; (test-damage k 1d7 init-stress)
          ;;                        ;; (delayed-damage-y ybar E Gf length init-stress)
          ;;                        (damage-response-exponential k E Gf length init-stress)
          ;;                        )))
          ;;   (setf damage-inc (- new-damage damage))
          ;;   ;; (setf damage-inc (* (/ dt tau) (- 1d0 (exp (- (* 1d0 (abs (- new-damage damage))))))))
          ;;   )
          ;; (setf damage (max damage (brittle-chalk-d k E Gf length init-stress))
          ;;       damage-inc 0d0)
          ;; (when (>= damage 1d0)
          ;;   (setf damage-inc 0d0)
          ;;   (setf ybar 0d0))
          ;; (incf (cl-mpm/particle::mp-time-averaged-damage-inc mp) damage-inc)
          ;; (incf (cl-mpm/particle::mp-time-averaged-ybar mp) ybar)
          ;; (incf (cl-mpm/particle::mp-time-averaged-counter mp))
          ;;Transform to log damage
          ;; (incf damage damage-inc)
          ;;Transform to linear damage
          (multiple-value-bind (ls v) (cl-mpm/utils:eig damage-tensor)
            (loop for i from 0 to 2
                  do
                     (let* ()
                       (setf (nth i ls) (max 0d0 (min 1d0 (nth i ls))))
                       (when (> (nth i ls) critical-damage)
                         (setf (nth i ls) 1d0))
                       (setf damage (max damage (nth i ls)))
                       )
                  )
            (setf damage-tensor (magicl:@ v
                                          (magicl:from-diag ls :type 'double-float)
                                          (magicl:transpose v)))
            )
          )
  (values)
  ))


(defmethod damage-model-calculate-y ((mp cl-mpm/particle::particle-concrete) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (pressure cl-mpm/particle::mp-pressure)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     (angle cl-mpm/particle::mp-friction-angle)
                     (c cl-mpm/particle::mp-coheasion)
                     ) mp
      (declare (double-float pressure damage))
        (progn
          (when (< damage 1d0)
            (let ((cauchy-undamaged (magicl:scale stress (/ 1d0 (magicl:det def)))))
              (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d cauchy-undamaged)
                (let* (;(s_1 (max 0d0 s_1))
                       )
                  (when (> s_1 0d0)
                    ;; (setf damage-increment s_1)
                    (setf damage-increment (sqrt
                                            (+ (expt s_1 2)
                                               (expt s_2 2)
                                               (expt s_3 2))))
                    )))))
          (when (>= damage 1d0)
            (setf damage-increment 0d0))
          ;;Delocalisation switch
          (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
          (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
          ))))

(defun brittle-concrete-d (stress E Gf length init-stress)
  (declare (double-float stress E Gf length init-stress))
  "Function that controls how damage evolves with principal stresses"
  (let* (
         (ft init-stress)
         (e0 (/ ft E))
         (ef (+ (/ e0 2) (/ Gf (* length ft))))
         (k (/ stress E))
         )

    (when (> length (/ (* 2 Gf E) (expt ft 2)))
      (error "Length scale is too long"))
    (if (> k e0)
        (- 1d0 (* (/ e0 k) (exp (- (/ (- k e0) (- ef e0))))))
        0d0
        ))
  ;; (let* ((hs (/ (expt stress 2) (* 2 E Gf)))
  ;;        (hsl (/ (* hs length) (- 1d0 (* hs length)))))
  ;;   (min 1d0 (max 0d0 (- 1d0 (exp (* -1d0 hs (/ init-stress (max init-stress stress)))))))
  ;;   )
  )

(defun test-damage (stress ef init-stress)
  (if (> stress init-stress)
        (* (/ ef (- ef init-stress)) (- 1d0 (/ init-stress stress)))
        0d0
        ))

(defun brittle-concrete-linear-d (stress E Gf length init-stress)
  (let* ((hsl (/ (expt init-stress 2) (* 2 E Gf)))
         (hs (/ (* hsl length) (- 1 (* hsl length)))))
    (if (> stress init-stress)
        (* (+ 1d0 hs) (- 1d0 (/ init-stress stress)))
        0d0
        )))

(defmethod update-damage ((mp cl-mpm/particle::particle-concrete) dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (E cl-mpm/particle::mp-e)
                     (Gf cl-mpm/particle::mp-Gf)
                     (log-damage cl-mpm/particle::mp-log-damage)
                     (damage-inc cl-mpm/particle::mp-damage-increment)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (pressure cl-mpm/particle::mp-pressure)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     ;; (length cl-mpm/particle::mp-true-local-length)
                     (length cl-mpm/particle::mp-local-length)
                     (k cl-mpm/particle::mp-history-stress)
                     ) mp
      (declare (double-float damage damage-inc critical-damage))
        (progn
          ;;Damage increment holds the delocalised driving factor
          (setf ybar damage-inc)
          (setf k (max k ybar))
          (let ((new-damage (max damage
                                 (brittle-concrete-d k E Gf length init-stress)
                                 ;; (brittle-concrete-linear-d k E Gf length init-stress)
                                 )))
            (setf damage-inc (- new-damage damage)))
          ;; (setf damage (max damage (brittle-chalk-d k E Gf length init-stress))
          ;;       damage-inc 0d0)
          (when (>= damage 1d0)
            (setf damage-inc 0d0)
            (setf ybar 0d0))
          (incf (cl-mpm/particle::mp-time-averaged-damage-inc mp) damage-inc)
          (incf (cl-mpm/particle::mp-time-averaged-ybar mp) ybar)
          (incf (cl-mpm/particle::mp-time-averaged-counter mp))
          ;;Transform to log damage
          (incf damage damage-inc)
          ;;Transform to linear damage
          (setf damage (max 0d0 (min 1d0 damage)))
          (when (> damage critical-damage)
            (setf damage 1d0)
            (setf damage-inc 0d0)))
  (values)
  ))


(defun damage-response-exponential (stress E Gf length init-stress ductility)
  (declare (double-float stress E Gf length init-stress))
  "Function that controls how damage evolves with principal stresses"
  (let* ((ft init-stress)
         (e0 (/ ft E))
         ;; (ef (+ (/ e0 2) (/ Gf (* length ft))))
         (ef (/ (* ft (+ ductility 1d0)) (* 2d0 E)))
         (k (/ stress E))
         ;; (beta (/ (* E e0 length) Gf))
         (beta (/ 1d0 (- ef e0)))
         ;(beta (/ (* E e0 length) Gf))
         )
    (when (> length (/ (* 2 Gf E) (expt ft 2)))
      (error "Length scale is too long, e0 ~F, Ef ~F, beta: ~F" e0 ef beta))
    (if (> k e0)
        (- 1d0 (* (/ e0 k) (+ 0d0 (exp (- (* beta (- k e0)))))))
        0d0
        ))
  ;; (let* ((hs (/ (expt stress 2) (* 2 E Gf)))
  ;;        (hsl (/ (* hs length) (- 1d0 (* hs length)))))
  ;;   (min 1d0 (max 0d0 (- 1d0 (exp (* -1d0 hs (/ init-stress (max init-stress stress)))))))
  ;;   )
  )


(defmethod update-damage ((mp cl-mpm/particle::particle-limestone) dt)
    (with-accessors ((stress cl-mpm/particle:mp-stress)
                     (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                     (strain cl-mpm/particle::mp-strain)
                     (damage cl-mpm/particle:mp-damage)
                     (E cl-mpm/particle::mp-e)
                     (Gf cl-mpm/particle::mp-Gf)
                     (log-damage cl-mpm/particle::mp-log-damage)
                     (damage-inc cl-mpm/particle::mp-damage-increment)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (pressure cl-mpm/particle::mp-pressure)
                     (volume cl-mpm/particle::mp-volume)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     (length cl-mpm/particle::mp-internal-length)
                     (k cl-mpm/particle::mp-history-stress)
                     (eng-inc cl-mpm/particle::mp-dissipated-energy-inc)
                     (eng-int cl-mpm/particle::mp-dissipated-energy)
                     (p cl-mpm/particle::mp-p-modulus)
                     (nu cl-mpm/particle::mp-nu)
                     (ductility cl-mpm/particle::mp-ductility)
                     ) mp
      (declare (double-float damage damage-inc critical-damage))
        (progn
          ;;Damage increment holds the delocalised driving factor

          (setf ybar damage-inc)
          (setf k (max k ybar))
          (let ((new-damage (max damage (damage-response-exponential k E Gf length init-stress ductility))))
            (setf damage-inc (- new-damage damage)))
          (when (>= damage 1d0)
            (setf damage-inc 0d0)
            (setf ybar 0d0))
          (incf (cl-mpm/particle::mp-time-averaged-damage-inc mp) damage-inc)
          (incf (cl-mpm/particle::mp-time-averaged-ybar mp) ybar)
          (incf (cl-mpm/particle::mp-time-averaged-counter mp))
          ;;Transform to log damage
          (incf damage damage-inc)
          ;;Transform to linear damage
          (setf damage (max 0d0 (min 1d0 damage)))
          (when (> damage critical-damage)
            (setf damage 1d0)
            ;(setf damage-inc 0d0)
            )
          (setf p (/ (* (- 1d0 damage) E) (* (+ 1d0 nu) (- 1d0 nu))))
          (incf eng-inc
                (* damage-inc
                   ;volume
                   0.5d0 (magicl:tref
                          (magicl:@
                           (magicl:transpose strain)
                           undamaged-stress
                           ) 0 0)))
          (incf eng-int eng-inc)
          )
  (values)
  ))


(defmethod damage-model-calculate-y ((mp cl-mpm/particle::particle-limestone) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((strain cl-mpm/particle::mp-strain)
                     (stress cl-mpm/particle::mp-undamaged-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (pressure cl-mpm/particle::mp-pressure)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     (de cl-mpm/particle::mp-elastic-matrix)
                     (E cl-mpm/particle::mp-e)
                     (nu cl-mpm/particle::mp-nu)
                     (k cl-mpm/particle::mp-compression-ratio)
                     ) mp
      (declare (double-float pressure damage))
        (progn
          (when (< damage 1d0)
            (let ((cauchy-undamaged (magicl:scale stress (/ 1d0 (* 1d0;(- 1d0 damage)
                                                                   (magicl:det def))))))
              (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d cauchy-undamaged)
                (let* (
                       (j2 (cl-mpm/constitutive::voigt-j2 (cl-mpm/constitutive::deviatoric-voigt
                                                                  cauchy-undamaged)))
                       (i1 (+ s_1 s_2 s_3))
                       (k-factor (/ (- k 1d0)
                                    (- 1d0 (* 2d0 nu))))
                       (s_1 (+ (* i1 (/ k-factor (* 2d0 k)))
                               (* (/ 1d0 (* 2d0 k))
                                  (sqrt (+ (expt (* k-factor i1) 2)
                                           (* (/ (* 12 k) (expt (- 1d0 nu) 2))j2)
                                           )))))
                       )
                  (setf damage-increment s_1)
                  ;; (when (> s_1 0d0)
                  ;;   (setf damage-increment s_1)
                  ;;   )
                  ))))


            ;; (multiple-value-bind (ls v) (cl-mpm/utils:eig (cl-mpm/utils:voigt-to-matrix stress))
            ;;   (multiple-value-bind (le v) (cl-mpm/utils:eig (cl-mpm/utils:voigt-to-matrix strain))
            ;;     (let* ((eps+ (sqrt (/ (loop for ps in ls
            ;;                                  for pe in le
            ;;                                 sum (* (max 0d0 ps) (max 0d0 pe))) E)))
            ;;            (eps- (sqrt (/ (loop for ps in ls
            ;;                                for pe in le
            ;;                                 sum (* (min 0d0 ps) (min 0d0 pe))) E)))
            ;;           (eps-eq (/ (+ (* k eps+) eps-) (* (+ k 1))))
            ;;           )
            ;;       (setf damage-increment (* eps-eq E))

            ;;       )
            ;;     ))

            ;; (multiple-value-bind (le v) (cl-mpm/utils:eig (cl-mpm/utils:voigt-to-matrix strain))
            ;;   (let ((strain+
            ;;           (magicl:@ v
            ;;                     (magicl:from-diag (mapcar (lambda (x) (max 0d0 x)) le) :type 'double-float)
            ;;                     (magicl:transpose v)))
            ;;         (strain-
            ;;           (magicl:@ v
            ;;                     (magicl:from-diag (mapcar (lambda (x) (min 0d0 x)) le) :type 'double-float)
            ;;                     (magicl:transpose v)))
            ;;           );))
            ;;   (let* ((eps+ (sqrt (/ (magicl:tref (magicl:@ (magicl:transpose! (matrix-to-voigt strain+))
            ;;                                                de
            ;;                                                (matrix-to-voight strain+)
            ;;                                                )0 0) E)))
            ;;          (eps- (sqrt (/ (magicl:tref (magicl:@ (magicl:transpose! (matrix-to-voigt strain-))
            ;;                                                de
            ;;                                                (matrix-to-voight strain-)
            ;;                                                )0 0) E)))
            ;;          (eps-eq (/ (+ (* k eps+) eps-) (* (+ k 1))))
            ;;          )
            ;;     (setf damage-increment (* eps-eq E))))))

                ;;       )
                ;;     ))
                ;; (let ((cauchy-undamaged (magicl:scale stress (/ 1d0 (magicl:det def)))))
                ;;   (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d cauchy-undamaged)
                ;;     (let* (;(s_1 (max 0d0 s_1))
                ;;            )
                ;;       (when (> s_1 0d0)
                ;;         ;; (setf damage-increment s_1)
                ;;         (setf damage-increment (sqrt
                ;;                                 (+ (expt s_1 2)
                ;;                                    (expt s_2 2)
                ;;                                    (expt s_3 2))))
                ;;         ))))
                )
              (when (>= damage 1d0)
                (setf damage-increment 0d0))
              ;;Delocalisation switch
              (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
              (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
              )))
