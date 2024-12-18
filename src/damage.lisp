(defpackage :cl-mpm/damage
  (:use :cl
        :cl-mpm/utils)
  (:export
    #:update-damage
    #:post-stress-step
    #:update-damage
    #:damage-model-calculate-y))

(in-package :cl-mpm/damage)
(declaim (optimize (debug 0) (safety 0) (speed 3)))
;; (declaim (optimize (debug 3) (safety 3) (speed 0)))

(defclass mpm-sim-damage (cl-mpm::mpm-sim-usf)
  ((delocal-counter
    :accessor sim-damage-delocal-counter
    :type fixnum
    :initarg :delocal-counter
    :initform 0)
   (delocal-counter-max
    :accessor sim-damage-delocal-counter-max
    :type fixnum
    :initarg :delocal-counter-max
    :initform 10))
  (:documentation "Explicit simulation with update stress first update"))

(defclass mpm-sim-usl-damage (mpm-sim-damage cl-mpm::mpm-sim-usl) ())

(defclass mpm-sim-damage-nd-2 (mpm-sim-damage cl-mpm::mpm-nd-2d) ())

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

            (when t;(< damage 1d0)
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

(defun calculate-damage (sim)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)
                   (dt cl-mpm:sim-dt)
                   (delocal-counter sim-damage-delocal-counter)
                   (delocal-counter-max sim-damage-delocal-counter-max)
                   (non-local-damage cl-mpm::sim-nonlocal-damage))
      sim

    (when non-local-damage
      (when (<= delocal-counter 0)
        ;;Ensure they have a home
        (update-delocalisation-list mesh mps)
        (setf delocal-counter delocal-counter-max))
      (decf delocal-counter))

    (cl-mpm:iterate-over-mps
     mps
     (lambda (mp)
       (when (typep mp 'cl-mpm/particle:particle-damage)
         ;; (find-nodal-local-length mesh mp)
         (setf (cl-mpm/particle::mp-true-local-length mp)
               (length-localisation (cl-mpm/particle::mp-local-length mp)
                                    (cl-mpm/particle::mp-local-length-damaged mp)
                                    (cl-mpm/particle::mp-damage mp)))
         (damage-model-calculate-y mp dt)
         )))

    (if non-local-damage
        (progn
          (delocalise-damage sim)
          )
        (localise-damage mesh mps dt))
    (cl-mpm:iterate-over-mps
     mps
      (lambda (mp)
        (when (typep mp 'cl-mpm/particle:particle-damage)
          (update-damage mp dt))))))

(defun create-delocalisation-list (mesh mps)
  (with-accessors ((nodes cl-mpm/mesh:mesh-nodes))
        mesh
    (cl-mpm:iterate-over-nodes
     mesh
     (lambda (node)
       (setf (fill-pointer (cl-mpm/mesh::node-local-list node)) 0)))
    (cl-mpm:iterate-over-mps
     mps
      (lambda (mp)
        (when (typep mp 'cl-mpm/particle:particle-damage)
            (let ((node-id (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp))))
              (when (cl-mpm/mesh:in-bounds mesh node-id)
                (let ((node (cl-mpm/mesh:get-node mesh node-id)))
                (sb-thread:with-mutex ((cl-mpm/mesh:node-lock node))
                  (vector-push-extend mp (cl-mpm/mesh::node-local-list node)))))))))
    ;;Smart in thought but actually super dumb
    ;; (cl-mpm:iterate-over-nodes
    ;;  mesh
    ;;  (lambda (node)
    ;;    (when (= 0 (length (cl-mpm/mesh::node-local-list node)))
    ;;      (adjust-array (cl-mpm/mesh::node-local-list node) 0 :fill-pointer 0))))
    ))

(defun local-list-add-particle (mesh mp)
  "A function for putting an MP into the nodal MP lookup table"
  (let ((node-id (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp))))
    (when (cl-mpm/mesh:in-bounds mesh node-id)
      (let ((node (cl-mpm/mesh:get-node mesh node-id)))
        (setf (cl-mpm/particle::mp-damage-position mp) (cl-mpm/utils::vector-copy (cl-mpm/particle:mp-position mp)))
        (sb-thread:with-mutex ((cl-mpm/mesh:node-lock node))
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
    (if (cl-mpm/particle::mp-damage-position mp)
        (let ((node-id (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle::mp-damage-position mp))))
          (if (cl-mpm/mesh:in-bounds mesh node-id)
            (let ((node (cl-mpm/mesh:get-node mesh node-id)))
              (when (not (remove-mp-ll node))
                  (cl-mpm::iterate-over-nodes-serial
                   mesh
                   (lambda (node)
                     (remove-mp-ll node))))
              t)
            nil))
        nil)))

(defun update-delocalisation-list (mesh mps)
  (with-accessors ((nodes cl-mpm/mesh:mesh-nodes)
                   (h cl-mpm/mesh:mesh-resolution))
        mesh
    (cl-mpm:iterate-over-mps
     mps
     (lambda (mp)
       (when (typep mp 'cl-mpm/particle:particle-damage)
         (if (eq (cl-mpm/particle::mp-damage-position mp) nil)
             ;;New particle - needs to be added to the mesh
             (local-list-add-particle mesh mp)
             ;; Already inserted mesh - sanity check to see if it should be recalced
             (let* ((delta (cl-mpm/fastmaths::diff-norm (cl-mpm/particle:mp-position mp)
                                                        (cl-mpm/particle::mp-damage-position mp))))
               (declare (double-float delta h))
               (when (> delta (/ h 4d0))
                 (when (not (eq
                             (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp))
                             (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle::mp-damage-position mp))
                             ))
                   (local-list-remove-particle mesh mp)
                   (local-list-add-particle mesh mp))))))))
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (node)
       (when (= 0 (length (cl-mpm/mesh::node-local-list node)))
         (adjust-array (cl-mpm/mesh::node-local-list node) 0 :fill-pointer 0))))))

;; (defgeneric update-delocalisation-list (sim))
;; (defmethod update-delocalisation-list ((sim cl-mpm::mpm-sim))
;;   (update-delocalisation-list-mps (cl-mpm:sim-mesh *sim*)
;;                                   (cl-mpm:sim-mps *sim*)))


;;This is a simd dot product
(progn

  (defun diff-squared-mat (pos-a pos-b)
    (let ((pos-a (magicl::matrix/double-float-storage pos-a))
          (pos-b (magicl::matrix/double-float-storage pos-b)))
      (values (the double-float (cl-mpm/fastmaths::dot-vector pos-a pos-b)))))

  (defun test-simd-acc (a b)
    (let (
          ;; (a (cl-mpm/utils:vector-from-list (list 1d0 2d0 3d0)))
          ;; (b (cl-mpm/utils:vector-from-list (list -1d0 5d0 8d0)))
          (a (cl-mpm/utils:vector-from-list (list 0d0 0d0 0d0)))
          (b (cl-mpm/utils:vector-from-list (list 1d0 1d0 1d0)))
          )
      (pprint (diff-squared-mat a b))
      ))

  (declaim
   (inline diff-squared)
   (ftype (function (cl-mpm/particle:particle cl-mpm/particle:particle) double-float) diff-squared))
  (defun diff-squared (mp-a mp-b)
    (cl-mpm/fastmaths::diff-norm
     (cl-mpm/particle:mp-position mp-a)
     (cl-mpm/particle:mp-position mp-b))
    )
  )

(declaim
 (inline diff-damaged)
 (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle cl-mpm/particle:particle) double-float) diff-damaged))
(defun diff-damaged (mesh mp-a mp-b)
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
                               (cl-mpm/fastmaths::fast-.+ step-point dhstep step-point)
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
                               (cl-mpm/fastmaths::fast-.+ step-point dhstep step-point)
                               )
                          )))
                (let* ((step-size remainder)
                       (dhstep (cl-mpm/utils::vector-copy step-norm)))
                  (magicl:scale! dhstep (* 0.5d0 step-size))
                  (progn
                    (cl-mpm/fastmaths::fast-.+ step-point dhstep step-point)
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
                   ) double-float)
        weight-func
        ))
(defun weight-func (dist-squared length)
  (if (< dist-squared (* length length))
      (the double-float (expt (- 1d0 (expt (/ dist-squared (* length length)) 2)) 2))
      ;; (the double-float (exp (the double-float (* 1d0 (/ (- dist-squared) (* 2d0 length length))))))
      0d0)
  )
;; (defun weight-func (dist-squared length)
;;   ;(values (the double-float (exp (the double-float (- (* (/ 4d0 (* length length)) dist-squared))))))
;;   (if (< dist-squared (* 4d0 length length))
;;       (values (the double-float (exp (the double-float (* 1d0 (/ (- dist-squared) (* 2d0 length length)))))))
;;       0d0)
;;   ;; (values (the double-float (exp (the double-float (* 4d0 (/ (- dist-squared) (* 1.00d0 length length)))))))
;;   ;; (if (= dist-squared 0d0)
;;   ;;   1d0
;;   ;;   0d0)
;;   )
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
  (weight-func (diff-squared mp-a mp-b) length)
  )

(defun weight-func-pos (mesh pos-a pos-b length)
  (weight-func (the double-float (cl-mpm/fastmaths::dot-vector pos-a pos-b)) length))

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

(defun patch-in-bounds-2d (mesh pos bound)
  (destructuring-bind (x y z) pos
    (declare (fixnum x y z bound))
    (and
     (cl-mpm/mesh:in-bounds mesh (list (+ x bound) (+ y bound) z))
     (cl-mpm/mesh:in-bounds mesh (list (- x bound) (- y bound) z)))))
(defun patch-in-bounds-3d (mesh pos bound)
  (destructuring-bind (x y z) pos
    (declare (fixnum x y z bound))
    (and
     (cl-mpm/mesh:in-bounds mesh (list (+ x bound) (+ y bound) (+ z bound)))
     (cl-mpm/mesh:in-bounds mesh (list (- x bound) (- y bound) (- z bound))))))

(defun iterate-over-damage-bounds (mesh mp length func)
  "Function for calling a function over every node that could contain mps within 2*length
Calls the function with the mesh mp and node"
  (if (= (cl-mpm/mesh:mesh-nd mesh) 2)
      (iterate-over-damage-bounds-2d mesh mp length func)
      (iterate-over-damage-bounds-3d mesh mp length func)))

(defun iterate-over-damage-bounds-2d (mesh mp length func)
  (declare (optimize (speed 3)))
  (declare (function func)
           (double-float length))
  (let* ((node-id (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp)))
         (node-reach (the fixnum (ceiling (* length 1d0) (the double-float (cl-mpm/mesh:mesh-resolution mesh)))))
         (potentially-in-bounds (patch-in-bounds-2d mesh node-id node-reach)))
    (declare (dynamic-extent node-id))
    (loop for dx fixnum from (- node-reach) to node-reach
          do (loop for dy fixnum from (- node-reach) to node-reach
                   do
                      (progn
                        (destructuring-bind (x y z) node-id
                          (declare (fixnum dx dy x y)
                                   (ignore z))
                          (let ((idx
                                (list (+ x dx) (+ y dy) 0)))
                            (declare (dynamic-extent idx))
                            (when (or
                                   potentially-in-bounds
                                   (cl-mpm/mesh:in-bounds mesh idx))
                              (let ((node (cl-mpm/mesh:get-node mesh idx)))
                                (funcall func mesh mp node))))))))))

(defun iterate-over-damage-bounds-3d (mesh mp length func)
  (let ((node-id (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp)))
        (node-reach (the fixnum (+ 0 (truncate (ceiling (* length 1d0)
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
(defparameter *enable-reflect-y* nil)
(defparameter *enable-reflect-z* nil)

(defun iterate-over-neighour-mps (mesh mp length func)
  "Search for mps in a patch sized 2*length, then return the "
  (declare (double-float length)
           (function func))
  (let ((len-squared (* length length)))
    (declare (double-float length len-squared))
    (iterate-over-damage-bounds
     mesh mp length
     (lambda (mesh mp node)
       (loop for mp-other across (cl-mpm/mesh::node-local-list node)
                                        ;(the (vector cl-mpm/particle::particle *))
             do
                (with-accessors ((d cl-mpm/particle::mp-damage)
                                 (m cl-mpm/particle:mp-volume)
                                 (ll cl-mpm/particle::mp-true-local-length)
                                 (p cl-mpm/particle:mp-position))
                    mp-other
                  (let ((distance (diff-squared mp mp-other)))
                    (when (< distance len-squared)
                      (funcall func mesh mp mp-other (sqrt distance)))))))))
  (values))


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
    (iterate-over-neighour-mps
     mesh mp length
     (lambda (mesh mp mp-other dist)
       (with-accessors ((d cl-mpm/particle::mp-damage)
                        (m cl-mpm/particle:mp-volume)
                        (ll cl-mpm/particle::mp-true-local-length)
                        (p cl-mpm/particle:mp-position))
           mp-other
         (when t;(< (the double-float d) 1d0)
           (let (
                 ;;Nodally averaged local funcj
                 ;; (weight (weight-func-mps mesh mp mp-other (* 0.5d0 (+ length ll))))
                 ;; (weight (weight-func-mps mesh mp mp-other (sqrt (* length ll))))
                 ;;
                 (weight
                   (weight-func-mps mesh mp mp-other (sqrt (* length ll)))
                   ;; (weight-func-mps mesh mp mp-other (* 0.5d0 (+ length ll)))
                   ;; (weight-func-mps-damaged mesh mp mp-other
                   ;;                          (cl-mpm/particle::mp-local-length mp)
                   ;;                          )
                   ))
             (declare (double-float weight m d mass-total damage-inc))
             (incf mass-total (* weight m))
             (incf damage-inc
                   (* (the double-float (cl-mpm/particle::mp-local-damage-increment mp-other))
                      weight m)))
           (macrolet ((reflect-axis (axis enable)
                        (declare (fixnum axis))
                        `(when (and ,enable
                                    (< (magicl:tref (cl-mpm/particle::mp-position mp) ,axis 0) (* 0.5d0 length)))
                           (let ((weight (weight-func-pos
                                          mesh
                                          (cl-mpm/particle::mp-position mp)
                                          (cl-mpm/fastmaths::fast-.* (cl-mpm/particle:mp-position mp-other)
                                                                     (cl-mpm/utils::vector-from-list
                                                                      (list ,(if (= axis 0) -1d0 0d0)
                                                                            ,(if (= axis 1) -1d0 0d0)
                                                                            ,(if (= axis 2) -1d0 0d0)
                                                                            )))
                                          (sqrt (* length ll))
                                          )))
                             (declare (double-float weight m d mass-total damage-inc))
                             (incf mass-total (* weight m))
                             (incf damage-inc
                                   (* (the double-float (cl-mpm/particle::mp-local-damage-increment mp-other))
                                      weight m))))))
             (reflect-axis 0 *enable-reflect-x*)
             (reflect-axis 1 *enable-reflect-y*)
             (reflect-axis 2 *enable-reflect-z*)
             )
           ))
       ))
    (when (> mass-total 0d0)
      (setf damage-inc (/ damage-inc mass-total)))
    damage-inc
    ;; (setf damage-inc (cl-mpm/particle::mp-local-damage-increment mp))
    ))

(declaim (notinline length-localisation))
(defun length-localisation (local-length local-length-damaged damage)
  ;; (+ (* local-length (- 1d0 damage)) (* local-length-damaged damage))
  ;; (* local-length (max (sqrt (- 1d0 damage)) 1d-10))
  local-length
  )


(defgeneric delocalise-damage (sim))

(defmethod delocalise-damage ((sim cl-mpm/damage::mpm-sim-damage))
  (declare (optimize (speed 3)))
  (with-accessors ((mps cl-mpm::sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
        ;; Calculate the delocalised damage for each damage particle
        (cl-mpm:iterate-over-mps
         mps
          (lambda (mp)
            (when (typep mp 'cl-mpm/particle:particle-damage)
              (with-accessors ((damage-inc cl-mpm/particle::mp-damage-increment)
                               (damage-inc-local cl-mpm/particle::mp-local-damage-increment)
                               (damage cl-mpm/particle::mp-damage)
                               (local-length-t cl-mpm/particle::mp-true-local-length))
                  mp
                (setf damage-inc (calculate-delocalised-damage mesh mp local-length-t)))))))
  (values))


(defun localise-damage (mesh mps dt)
  "Apply local damage model"
  (cl-mpm:iterate-over-mps
   mps
    (lambda (mp)
      (when (typep mp 'cl-mpm/particle:particle-damage)
        (with-accessors ((damage-inc cl-mpm/particle::mp-damage-increment)
                         (damage-inc-local cl-mpm/particle::mp-local-damage-increment)
                         (local-length cl-mpm/particle::mp-local-length)
                         )
            mp
          (setf damage-inc damage-inc-local)))))
  (values))


(defun remove-mp-damage (sim func)
  )

(defmethod cl-mpm::remove-mps-func ((sim mpm-sim-damage) func)
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (mps cl-mpm::sim-mps))
      sim
    (when (> (length mps) 0)
      (loop for mp across mps
            when (and
                  (typep mp 'cl-mpm/particle:particle-damage)
                  (funcall func mp))
              do (local-list-remove-particle mesh mp))
      (call-next-method))))
(defmethod cl-mpm::sim-add-mp ((sim mpm-sim-damage) mp)
  (local-list-add-particle (cl-mpm:sim-mesh sim) mp)
  (call-next-method))

(defmethod (setf cl-mpm::sim-mps) (mps (sim cl-mpm/damage::mpm-sim-damage))
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
               (update-type cl-mpm::update-type)
               (vel-algo cl-mpm::velocity-algorithm)
               )
                sim
    (declare (type double-float mass-filter))
                (progn
                    (cl-mpm::reset-grid mesh)
                    (cl-mpm::p2g mesh mps)
                    (cl-mpm::check-single-mps sim)
                    (when (> mass-filter 0d0)
                      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
                    (cl-mpm::update-node-kinematics mesh dt )
                    (cl-mpm::apply-bcs mesh bcs dt)

                    (cl-mpm::update-stress mesh mps dt fbar)

                    (when enable-damage
                     (cl-mpm/damage::calculate-damage sim))
                    ;; ;Map forces onto nodes
                    (cl-mpm::p2g-force mesh mps)

                    (loop for bcs-f in bcs-force-list
                          do (cl-mpm::apply-bcs mesh bcs-f dt))

                    (cl-mpm::update-node-forces sim)
                    ;; ;Reapply velocity BCs
                    (cl-mpm::apply-bcs mesh bcs dt)
                    ;; ;Also updates mps inline
                    (cl-mpm::g2p mesh mps dt vel-algo)

                    (when remove-damage
                      (cl-mpm::remove-material-damaged sim))
                    (when split
                      (cl-mpm::split-mps sim))
                    ;; (cl-mpm::check-mps sim)
                    )))

(defmethod cl-mpm::update-sim ((sim mpm-sim-usl-damage))
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
               (vel-algo cl-mpm::velocity-algorithm)
               )
                sim
    (declare (type double-float mass-filter))
                (progn
                    (cl-mpm::reset-grid mesh)
                    (cl-mpm::p2g mesh mps)
                    (cl-mpm::check-single-mps sim)
                    ;;Do optional mass filter
                    (when (> mass-filter 0d0)
                      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
                    (cl-mpm::update-node-kinematics mesh dt )
                    (cl-mpm::apply-bcs mesh bcs dt)
                    ;Map forces onto nodes
                    (cl-mpm::p2g-force mesh mps)
                    (loop for bcs-f in bcs-force-list
                          do
                             (cl-mpm::apply-bcs mesh bcs-f dt))
                    (cl-mpm::update-node-forces sim)
                    ;Reapply velocity BCs
                    (cl-mpm::apply-bcs mesh bcs dt)
                    ;Also updates mps inline


                    (cl-mpm::g2p mesh mps dt vel-algo)

                    ;;Update stress last

                    (cl-mpm::reset-grid-velocity mesh)
                    (cl-mpm::p2g mesh mps)
                    (cl-mpm::check-single-mps sim)
                    ;;Do optional mass filter
                    (when (> mass-filter 0d0)
                      (cl-mpm::filter-grid-velocity mesh (cl-mpm::sim-mass-filter sim)))
                    (cl-mpm::update-node-kinematics mesh dt )
                    (cl-mpm::apply-bcs mesh bcs dt)

                    (cl-mpm::update-stress mesh mps dt fbar)
                    (when enable-damage
                      (cl-mpm/damage::calculate-damage sim))

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
            ;; (setf damage-increment 0d0)
            )
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
          ;; (setf damage-inc 0d0)
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
          ;; (setf damage-inc 0d0)
          (setf ybar 0d0))

        (incf damage damage-inc)
        (setf damage (max 0d0 (min 1d0 damage)))
        (when (> damage critical-damage)
          (setf damage 1d0)
          (setf damage-inc 0d0))))
  (values)))


(defmethod cl-mpm/output::save-vtk (filename (sim cl-mpm/damage::mpm-sim-damage))
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)) sim
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
      (let ((id 1)
            (nd (cl-mpm/mesh:mesh-nd mesh)))
        (declare (special id))
        (format fs "POINT_DATA ~d~%" (length mps))

        (cl-mpm/output::save-parameter "mass" (cl-mpm/particle:mp-mass mp))
        (cl-mpm/output::save-parameter "density" (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)))
        (cl-mpm/output::save-parameter "index" (cl-mpm/particle::mp-index mp))
        (cl-mpm/output::save-parameter "mpi-index" (cl-mpm/particle::mp-mpi-index mp))
        (cl-mpm/output::save-parameter "j" (magicl:det (cl-mpm/particle::mp-deformation-gradient mp)))
        (cl-mpm/output::save-parameter "volume" (cl-mpm/particle::mp-volume mp))
        (cl-mpm/output::save-parameter "vel_x" (magicl:tref (cl-mpm/particle:mp-velocity mp) 0 0))
        (cl-mpm/output::save-parameter "vel_y" (magicl:tref (cl-mpm/particle:mp-velocity mp) 1 0))
        ;; (cl-mpm/output::save-parameter "acc_x" (magicl:tref (cl-mpm/particle::mp-acceleration mp) 0 0))
        ;; (cl-mpm/output::save-parameter "acc_y" (magicl:tref (cl-mpm/particle::mp-acceleration mp) 1 0))

        (cl-mpm/output::save-parameter "disp_x" (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0))
        (cl-mpm/output::save-parameter "disp_y" (magicl:tref (cl-mpm/particle::mp-displacement mp) 1 0))
        (when (= 3 nd)
          (cl-mpm/output::save-parameter "disp_z" (magicl:tref (cl-mpm/particle::mp-displacement mp) 2 0)))

        (cl-mpm/output::save-parameter "sig_xx" (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0))
        (cl-mpm/output::save-parameter "sig_yy" (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0))
        (cl-mpm/output::save-parameter "sig_xy" (magicl:tref (cl-mpm/particle:mp-stress mp) 5 0))
        (when (= 3 nd)
          (cl-mpm/output::save-parameter "sig_zz" (magicl:tref (cl-mpm/particle:mp-stress mp) 2 0))
          (cl-mpm/output::save-parameter "sig_yz" (magicl:tref (cl-mpm/particle:mp-stress mp) 3 0))
          (cl-mpm/output::save-parameter "sig_zx" (magicl:tref (cl-mpm/particle:mp-stress mp) 4 0)))

        (cl-mpm/output::save-parameter "eps_xx" (magicl:tref (cl-mpm/particle:mp-strain mp) 0 0))
        (cl-mpm/output::save-parameter "eps_yy" (magicl:tref (cl-mpm/particle:mp-strain mp) 1 0))
        (cl-mpm/output::save-parameter "eps_xy" (magicl:tref (cl-mpm/particle:mp-strain mp) 5 0))
        (when (= 3 nd)
          (cl-mpm/output::save-parameter "eps_zz" (magicl:tref (cl-mpm/particle:mp-strain mp) 2 0))
          (cl-mpm/output::save-parameter "eps_yz" (magicl:tref (cl-mpm/particle:mp-strain mp) 3 0))
          (cl-mpm/output::save-parameter "eps_zx" (magicl:tref (cl-mpm/particle:mp-strain mp) 4 0)))

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
        (cl-mpm/output::save-parameter "pressure" (cl-mpm/particle::mp-pressure mp))
        (cl-mpm/output::save-parameter "boundary" (cl-mpm/particle::mp-boundary mp))

        (cl-mpm/output::save-parameter "s_1"
                                       (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))
                                         (loop for sii in l maximize sii)))
        (cl-mpm/output::save-parameter "e_1"
                                       (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-strain mp)))
                                         (loop for sii in l maximize sii)))


        (cl-mpm/output::save-parameter "EPS"
                                       (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))
                                         (- (loop for sii in l maximize sii) (cl-mpm/particle::mp-pressure mp))))
        (cl-mpm/output::save-parameter "size_x" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 0 0))
        (cl-mpm/output::save-parameter "size_y" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0))
        (when (= 3 nd)
          (cl-mpm/output::save-parameter "size_z" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 2 0)))

        (cl-mpm/output::save-parameter "damage"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage)
                                           (cl-mpm/particle:mp-damage mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-shear"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-shear)
                                           (cl-mpm/particle::mp-damage-shear mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-compression"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-compression)
                                           (cl-mpm/particle::mp-damage-compression mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-tension"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-tension)
                                           (cl-mpm/particle::mp-damage-tension mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-inc"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-increment)
                                           (cl-mpm/particle::mp-damage-increment mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-ybar"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-ybar)
                                           (cl-mpm/particle::mp-damage-ybar mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-ybar-scaled"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-ybar)
                                           (*
                                            (- 1d0 (cl-mpm/particle::mp-damage mp))
                                            (cl-mpm/particle::mp-damage-ybar mp)
                                            )
                                           0d0))
        (cl-mpm/output::save-parameter "damage-y"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-y-local)
                                           (cl-mpm/particle::mp-damage-y-local mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-k"
                                       (if (slot-exists-p mp 'cl-mpm/particle::history-stress)
                                           (cl-mpm/particle::mp-history-stress mp)
                                           0d0))
        (cl-mpm/output::save-parameter "viscosity"
                                       (if (slot-exists-p mp 'cl-mpm/particle::true-visc)
                                           (cl-mpm/particle::mp-true-visc mp)
                                           0d0))
        ;; (cl-mpm/output::save-parameter "damage-xx"
        ;;                                (if (slot-exists-p mp 'cl-mpm/particle::damage-tensor)
        ;;                                    (magicl:tref (cl-mpm/particle::mp-damage-tensor mp) 0 0)
        ;;                                    0d0))
        ;; (cl-mpm/output::save-parameter "damage-yy"
        ;;                                (if (slot-exists-p mp 'cl-mpm/particle::damage-tensor)
        ;;                                    (magicl:tref (cl-mpm/particle::mp-damage-tensor mp) 1 1)
        ;;                                    0d0))
        ;; (cl-mpm/output::save-parameter "damage-y-xx"
        ;;                                (if (slot-exists-p mp 'cl-mpm/particle::damage-tensor)
        ;;                                    (magicl:tref (cl-mpm/particle::mp-damage-ybar-tensor mp) 0 0)
        ;;                                    0d0))
        ;; (cl-mpm/output::save-parameter "damage-y-yy"
        ;;                                (if (slot-exists-p mp 'cl-mpm/particle::damage-tensor)
        ;;                                    (magicl:tref (cl-mpm/particle::mp-damage-ybar-tensor mp) 1 1)
        ;;                                    0d0))
        
        ;; (cl-mpm/output::save-parameter "damage-plastic-phi"
        ;;                                (if (and (slot-exists-p mp 'cl-mpm/particle::damage-ybar)
        ;;                                         (slot-exists-p mp 'cl-mpm/particle::phi))
        ;;                                    (with-accessors ((kt-r cl-mpm/particle::mp-k-tensile-residual-ratio)
        ;;                                                     (kc-r cl-mpm/particle::mp-k-compressive-residual-ratio)
        ;;                                                     (g-r cl-mpm/particle::mp-shear-residual-ratio)
        ;;                                                     (angle-plastic cl-mpm/particle::mp-phi)
        ;;                                                     (stress cl-mpm/particle::mp-undamaged-stress)
        ;;                                                     (damage cl-mpm/particle::mp-damage)
        ;;                                                     )
        ;;                                        mp
        ;;                                      (let ((p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
        ;;                                            (s (cl-mpm/constitutive::deviatoric-voigt stress))
        ;;                                            (p-r 0d0))
        ;;                                        (setf p-r
        ;;                                              (if (> p 0d0)
        ;;                                                  (- 1d0 (* (- 1d0 kt-r) damage))
        ;;                                                  (- 1d0 (* (- 1d0 kc-r) damage))))
        ;;                                        (* (/ 180d0 pi) (atan (* (/ (- 1d0 (* (- 1d0 g-r) damage)) p-r) (tan angle-plastic)))))
        ;;                                      )
        ;;                                    0d0))
        ;; (cl-mpm/output::save-parameter "damage-plastic-c"
        ;;                                (if (and (slot-exists-p mp 'cl-mpm/particle::damage-ybar)
        ;;                                         (slot-exists-p mp 'cl-mpm/particle::phi))
        ;;                                    (with-accessors ((kt-r cl-mpm/particle::mp-k-tensile-residual-ratio)
        ;;                                                     (kc-r cl-mpm/particle::mp-k-compressive-residual-ratio)
        ;;                                                     (g-r cl-mpm/particle::mp-shear-residual-ratio)
        ;;                                                     (c cl-mpm/particle::mp-c)
        ;;                                                     (stress cl-mpm/particle::mp-undamaged-stress)
        ;;                                                     (damage cl-mpm/particle::mp-damage)
        ;;                                                     )
        ;;                                        mp
        ;;                                      (let ((p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
        ;;                                            (s (cl-mpm/constitutive::deviatoric-voigt stress))
        ;;                                            (p-r 0d0))
        ;;                                        (setf p-r
        ;;                                              (if (> p 0d0)
        ;;                                                  (- 1d0 (* (- 1d0 kt-r) damage))
        ;;                                                  (- 1d0 (* (- 1d0 kc-r) damage))))
        ;;                                        (* c p-r)))
        ;;                                    0d0))
        
        (cl-mpm/output::save-parameter "local-length"
                                       (if (slot-exists-p mp 'cl-mpm/particle::true-local-length)
                                           (cl-mpm/particle::mp-true-local-length mp)
                                           0d0))
        ;; (cl-mpm/output::save-parameter "j2"
        ;;                                (with-accessors ((damage cl-mpm/particle::mp-damage)
        ;;                                                 (stress cl-mpm/particle::mp-undamaged-stress)
        ;;                                                 (def cl-mpm/particle::mp-deformation-gradient)
        ;;                                                 )
        ;;                                    mp
        ;;                                  (if (< damage 1d0)
        ;;                                      (sqrt (cl-mpm/constitutive::voigt-j2
        ;;                                             (cl-mpm/utils::deviatoric-voigt stress)))
        ;;                                      0d0)))
        (cl-mpm/output::save-parameter "i1"
                                       (if (slot-exists-p mp 'cl-mpm/particle::undamaged-stress)
                                           (cl-mpm/utils::trace-voigt (cl-mpm/particle::mp-undamaged-stress mp))
                                           (cl-mpm/utils::trace-voigt (cl-mpm/particle::mp-stress mp))))
        (cl-mpm/output::save-parameter "stress-pressure"
                                       (/ (if (slot-exists-p mp 'cl-mpm/particle::undamaged-stress)
                                              (cl-mpm/utils::trace-voigt (cl-mpm/particle::mp-undamaged-stress mp))
                                              (cl-mpm/utils::trace-voigt (cl-mpm/particle::mp-stress mp))) 3))
        (cl-mpm/output::save-parameter "e_i1"
                                       (with-accessors ((strain cl-mpm/particle::mp-strain))
                                           mp
                                         (multiple-value-bind (e_1 e_2 e_3) (principal-stresses-3d strain)
                                           (+ e_1 e_2 e_3))))
        (cl-mpm/output::save-parameter "e_j2"
                                       (with-accessors ((strain cl-mpm/particle::mp-strain))
                                           mp
                                         (multiple-value-bind (e_1 e_2 e_3) (principal-stresses-3d strain)
                                           (let* ((j2
                                                    (* 1/6 (+
                                                            (expt (- e_1 e_2) 2)
                                                            (expt (- e_2 e_3) 2)
                                                            (expt (- e_3 e_1) 2)))))
                                             (sqrt j2)))))


        (cl-mpm/output::save-parameter "split-depth"
                                       (cl-mpm/particle::mp-split-depth mp))
        (cl-mpm/output::save-parameter
         "plastic_strain"
         (if (slot-exists-p mp 'cl-mpm/particle::yield-func)
             (cl-mpm/particle::mp-strain-plastic-vm mp)
             0d0)
         )
        (cl-mpm/output::save-parameter
         "yield-func"
         (if (slot-exists-p mp 'cl-mpm/particle::yield-func)
             (cl-mpm/particle::mp-yield-func mp)
             0d0))

        (cl-mpm/output::save-parameter
         "plastic-c"
         (if (slot-exists-p mp 'cl-mpm/particle::c)
             (cl-mpm/particle::mp-c mp)
             0d0))
        (cl-mpm/output::save-parameter
         "plastic-phi"
         (if (slot-exists-p mp 'cl-mpm/particle::phi)
             (* (cl-mpm/particle::mp-phi mp) (/ 180 pi))
             0d0))

        (cl-mpm/output::save-parameter
         "energy"
         ;; (cl-mpm/particle::mp-penalty-energy mp)
         (* (cl-mpm/particle::mp-mass mp)
            (cl-mpm/fastmaths::mag-squared (cl-mpm/particle::mp-velocity mp)))
         )

        (cl-mpm/output::save-parameter
         "fbar-j"
         (cl-mpm/particle::mp-debug-j mp)
         )
        (cl-mpm/output::save-parameter
         "fbar-j-gather"
         (cl-mpm/particle::mp-debug-j-gather mp)
         )
        (cl-mpm/output::save-parameter
         "fbar-j-diff"
         (if (> (cl-mpm/particle::mp-debug-j mp) 0d0) 
             (/ (- (cl-mpm/particle::mp-debug-j-gather mp) (cl-mpm/particle::mp-debug-j mp)) 
                (cl-mpm/particle::mp-debug-j mp)
                )
             0d0)
         )
        (cl-mpm/output::save-parameter "fric-contact" (if (cl-mpm/particle::mp-penalty-contact-step mp) 1 0))
        (cl-mpm/output::save-parameter "fric-contact-stick" (if (cl-mpm/particle::mp-penalty-friction-stick mp) 1 0))
        (cl-mpm/output::save-parameter "fric-normal" (cl-mpm/particle::mp-penalty-normal-force mp))
        (cl-mpm/output::save-parameter "fric-x" (magicl:tref (cl-mpm/particle::mp-penalty-frictional-force mp) 0 0))
        (cl-mpm/output::save-parameter "fric-y" (magicl:tref (cl-mpm/particle::mp-penalty-frictional-force mp) 1 0))
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
            ;; (setf damage-inc 0d0)
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



(defun damage-rate-profile-chalk (stress damage rate init-stress)
  (declare (double-float stress damage rate init-stress))
  "Function that controls how damage evolves with principal stresses"
  (if (> stress init-stress)
      (* (expt (max 0d0 (- stress init-stress)) 0.5d0) rate)
      0d0))




;;Chalk brittle

(defun drucker-prager-criterion (stress angle)
  (let ((p (cl-mpm/utils::trace-voigt stress))
        (j2 (cl-mpm/constitutive::voigt-j2
             (cl-mpm/utils::deviatoric-voigt stress))))
    (* (/ 3d0 (+ 3 (tan angle)))
       (+ (sqrt (* 3 j2)) (* 1/3 (tan angle) p)))))


(defun criterion-dp (stress angle)
  (criterion-dp-coheasion stress angle)
  ;; (criterion-dp-tensile stress angle)
  )

(defun criterion-dp-coheasion (stress angle)
  "Return some drucker-prager damage criterion from a stress level and and angle (radians)"
  (let ((p (cl-mpm/utils::trace-voigt stress))
        (j2 (cl-mpm/constitutive::voigt-j2
             (cl-mpm/utils::deviatoric-voigt stress)))
        (A (/ (* 6 (cos angle))
              (* (sqrt 3) (- 3d0 (sin angle)))))
        (B (/ (* 2 (sin angle))
              (* (sqrt 3) (- 3d0 (sin angle))))))
    (* (/ 1d0 A)
       (+ (* B p) (sqrt j2)))))

(defun criterion-dp-middle-circumscribe-coheasion (stress angle)
  "Return some drucker-prager damage criterion from a stress level and and angle (radians)"
  (let ((p (cl-mpm/utils::trace-voigt stress))
        (j2 (cl-mpm/constitutive::voigt-j2
             (cl-mpm/utils::deviatoric-voigt stress)))
        (A (/ (* 6 (cos angle))
              (* (sqrt 3) (+ 3d0 (sin angle)))))
        (B (/ (* 2 (sin angle))
              (* (sqrt 3) (+ 3d0 (sin angle))))))
    (* (/ 1d0 A)
       (+ (* B p) (sqrt j2)))))


(defun criterion-dp-inscribe-coheasion (stress angle)
  "Return some drucker-prager damage criterion from a stress level and and angle (radians)"
  (let* ((p (/ (cl-mpm/utils::trace-voigt stress) 3))
         (j2 (cl-mpm/constitutive::voigt-j2
              (cl-mpm/utils::deviatoric-voigt stress)))
         (factor (sqrt (+ 9 (* 3 (expt (sin angle) 2)))))
         (A (/ (* 3 (cos angle))
               factor))
         (B (/ (* (sin angle))
               factor)))
    (* (/ 1d0 A)
       ;; (/ 1d0 (- (/ 1d0 (sqrt 3)) B))
       (+ (* B p) (sqrt j2)))))


(defun degredation-function (stress damage kc kt ks)
  (let ((p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
        (s (cl-mpm/constitutive::deviatoric-voigt stress))
        (damage-t damage)
        (damage-c (* kc damage))
        (damage-s (* ks damage)))
    (declare (double-float damage-t damage-c damage-s))
    (setf p
          (if (> p 0d0)
              (* (- 1d0 damage-t) p)
              (* (- 1d0 damage-c) p)))
    (cl-mpm/fastmaths:fast-.+
     (cl-mpm/constitutive::voight-eye p)
     (cl-mpm/fastmaths:fast-scale! s (- 1d0 damage-s)))))

(defun criterion-y (stress strain damage
                    E
                    kc
                    kt
                    ks)
  (let ((d_0 damage)
        (ddamage (- (min 1d0 (+ damage 1d-3))
                    damage)))
    (if (> ddamage 0d0)
        (let* ((d_1 (+ damage ddamage))
               (se_0 (* 0.5d0
                        (cl-mpm/fastmaths:dot
                         stress
                         strain)))
               (se_1 (* 0.5d0
                        (cl-mpm/fastmaths:dot
                         (degredation-function stress
                                               d_1
                                               kc
                                               kt
                                               ks)
                         strain))))
          (sqrt (* E (max 0d0 (- (/ (- se_1 se_0) ddamage))))))
        0d0)))

;; (let* ((strain (cl-mpm/utils:voigt-from-list (list
;;                                               1d0 1d0 1d0
;;                                               0d0 0d0 1d0
;;                                               )))
;;        (E 1d0)
;;        (nu 0d0)
;;        (de (cl-mpm/constitutive::linear-elastic-matrix E nu))
;;        (stress (magicl:@ de strain))
;;        (kt-r 1d0)
;;        (kc-r 1d0)
;;        (g-r 1d0))
;;   ;; (pprint (degredation-function stress
;;   ;;                               0.1d0
;;   ;;                               1d0
;;   ;;                               1d0
;;   ;;                               1d0
;;   ;;                               ))
;;   (pprint (* 0.5d0
;;              (cl-mpm/fastmaths:dot
;;               stress
;;               strain)))
;;   ;(pprint (criterion-y stress strain damage E kt-r kc-r g-r))
;;   (pprint (volumetric-deviatoric-norm strain E nu kt-r kc-r g-r))
;;   )


(defun criterion-dp-shear-coheasion (stress angle)
  "Return some drucker-prager damage criterion from a stress level and and angle (radians)"
  (let* ((p (cl-mpm/utils::trace-voigt stress))
         (j2 (cl-mpm/constitutive::voigt-j2
              (cl-mpm/utils::deviatoric-voigt stress)))
         (A (cos angle))
         (B (/ (sin angle) 3d0)))
    (* (/ 1d0 A)
       (+ (* B p) (sqrt j2)))))


(defun dp-tension-from-DP-circumscribed (angle coheasion)
  (declare (double-float angle))
  (let* ((factor (* (sqrt 3) (- 3d0 (sin angle))))
         (A (/ (* 6 coheasion (cos angle))
               factor))
         (B (/ (* 2 (sin angle))
               factor)))
    (declare (double-float A B factor))
    (* (/ A (- (/ 1d0 (sqrt 3)) B)))))

(defun dp-tension-from-DP-inscribed (angle coheasion)
  (declare (double-float angle))
  (let* ((factor (sqrt (+ 9 (* 3 (expt (sin angle) 2)))))
         (A (/ (* 3 coheasion (cos angle)) factor))
         (B (/ (sin angle) factor)))
    (declare (double-float A B factor))
    (* (/ A (- (/ 1d0 (sqrt 3)) B)))))


(defun criterion-dp-tensile (stress angle)
  (declare (double-float angle))
  (let* ((p (cl-mpm/utils::trace-voigt stress))
         (j2 (cl-mpm/constitutive::voigt-j2
              (cl-mpm/utils::deviatoric-voigt stress)))
         (factor (* (sqrt 3) (- 3d0 (sin angle))))
         (A (/ (* 6 (cos angle))
               factor))
         (B (/ (* 2 (sin angle))
               factor)))
    (declare (double-float A B factor j2 p))
    (* (/ 1d0 (- (/ 1d0 (sqrt 3)) B))
       (+ (* B p) (the double-float (sqrt j2))))))

(defun criterion-dp-tensile-circumscribed (stress angle)
  (declare (double-float angle))
  (let* ((p (cl-mpm/utils::trace-voigt stress))
         (j2 (cl-mpm/constitutive::voigt-j2
              (cl-mpm/utils::deviatoric-voigt stress)))
         (factor (* (sqrt 3) (- 3d0 (sin angle))))
         (A (/ (* 6 (cos angle))
               factor))
         (B (/ (* 2 (sin angle))
               factor)))
    (declare (double-float A B factor j2 p))
    (* (/ 1d0 (- (/ 1d0 (sqrt 3)) B))
       (+ (* B p) (the double-float (sqrt j2))))))

(defun criterion-dp-tensile-middle-circumscribed (stress angle)
  (declare (double-float angle))
  (let* ((p (cl-mpm/utils::trace-voigt stress))
         (j2 (cl-mpm/constitutive::voigt-j2
              (cl-mpm/utils::deviatoric-voigt stress)))
         (factor (* (sqrt 3) (+ 3d0 (sin angle))))
         (A (/ (* 6 (cos angle))
               factor))
         (B (/ (* 2 (sin angle))
               factor)))
    (declare (double-float A B factor j2 p))
    (* (/ 1d0 (- (/ 1d0 (sqrt 3)) B))
       (+ (* B p) (the double-float (sqrt j2))))))

(defun criterion-dp-tensile-inscribed (stress angle)
  (declare (double-float angle))
  (let* ((p (cl-mpm/utils::trace-voigt stress))
         (j2 (cl-mpm/constitutive::voigt-j2
              (cl-mpm/utils::deviatoric-voigt stress)))
         (factor (* (sqrt (+ 9 (* 3 (expt (sin angle) 2))))))
         (A (/ (* 3 (cos angle))
               factor))
         (B (/ (sin angle)
               factor)))
    (declare (double-float A B factor j2 p))
    (* (/ 1d0 (- (/ 1d0 (sqrt 3)) B))
       (+ (* B p) (the double-float (sqrt j2))))))


(defun criterion-dp-pressure (stress angle pressure)
  (let ((i1 (+ (cl-mpm/utils::trace-voigt stress) (* 3d0 pressure)))
        (j2 (cl-mpm/constitutive::voigt-j2
             (cl-mpm/utils::deviatoric-voigt stress)))
        (B (/ (* 2 (sin angle))
              (* (sqrt 3) (- 3d0 (sin angle))))))
    (* (+ (* B i1) (sqrt j2)))))

(defun criterion-dp-pressure-tensile (stress angle pressure tension)
  (let ((i1 (+ (cl-mpm/utils::trace-voigt stress) (* 3d0 pressure)))
        (j2 (cl-mpm/constitutive::voigt-j2
             (cl-mpm/utils::deviatoric-voigt stress)))
        (B (/ (* 2 (sin angle))
              (* (sqrt 3) (- 3d0 (sin angle)))))
        )
    (if nil;(> i1 (* tension 3))
        (/ i1 3)
        (*  (+ (* B i1) (sqrt j2)))
        )))



;; (defun criterion-dp-strain (stress ft fc)
;;   (let ((p (cl-mpm/utils::trace-voigt stress))
;;         (j2 (cl-mpm/constitutive::voigt-j2
;;              (cl-mpm/utils::deviatoric-voigt
;;               stress
;;               )))
;;         )
;;     (*
;;      (/ 1d0 (* 2 fc))
;;      (- (* (+ fc ft) (sqrt (* 3d0 j2)))
;;         (* (- ft fc) p)))))

(defun criterion-dp-strain (stress k)
  (let ((p (cl-mpm/utils::trace-voigt stress))
        (j2 (cl-mpm/constitutive::voigt-j2
             (cl-mpm/utils::deviatoric-voigt
              stress
              ))))
    (* (/ (sqrt 3) (* 2 (+ k 1)))
       (-
        (sqrt j2)
        (*
         (/ (- k 1)
            (+ k 1))
         p)))
    ))

(defun criterion-j2 (stress)
  (let ((j2 (cl-mpm/constitutive::voigt-j2
             (cl-mpm/utils::deviatoric-voigt stress))))
     (sqrt j2)))

(defun modified-vm-criterion (stress nu k)
  (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d stress)
    (declare (double-float nu k s_1 s_2 s_3))
    (let* ((j2 (cl-mpm/constitutive::voigt-j2
                (cl-mpm/utils::deviatoric-voigt stress)))
           (i1 (+ s_1 s_2 s_3))
           (k-factor (/ (- k 1d0)
                        (- 1d0 (* 2d0 nu))))
           (s_1 (+ (* i1 (/ k-factor (* 2d0 k)))
                   (* (/ 1d0 (* 2d0 k))
                      (sqrt (+ (expt (* k-factor i1) 2)
                               (* (/ (* 12 k) (expt (- 1d0 nu) 2))j2)
                               ))))))
      s_1
      )))
(defun smooth-rankine-criterion (stress)
  (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d stress)
     (sqrt
      (+ (expt (max 0d0 s_1) 2)
         (expt (max 0d0 s_2) 2)
         (expt (max 0d0 s_3) 2)))
     ))
(defun criterion-max-principal-stress (stress)
  (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d stress)
    (declare (double-float s_1 s_2 s_3))
    (max 0d0 s_1 s_2 s_3)))

;; (defun modified-vm-strain (strain nu k E)
;;   (multiple-value-bind (e1 e2 e3)
;;       (principal-stresses-3d
;;        (magicl:.*
;;         strain
;;         (cl-mpm/utils:voigt-from-list
;;          (list 1d0 1d0 1d0
;;                0.5d0
;;                0.5d0
;;                0.5d0))))
;;     (let* ((i1 (+ e1 e2 e3))
;;            (j2
;;              ;; (cl-mpm/constitutive::voigt-j2 (cl-mpm/utils::deviatoric-voigt strain))
;;              (* 1/6 (+
;;                      (expt (- e1 e2) 2)
;;                      (expt (- e2 e3) 2)
;;                      (expt (- e3 e1) 2)))
;;              )
;;            (k-factor (/ (- k 1d0)
;;                         (- 1d0 (* 2d0 nu))))
;;            (s_1 (* E
;;                    (+ (* i1 (/ k-factor (* 2d0 k)))
;;                       (* (/ 1d0 (* 2d0 k))
;;                          (sqrt (+ (expt (* k-factor i1) 2)
;;                                   (* (/ (* 12 k) (expt (- 1d0 nu) 2)) j2)
;;                                   )))))))
;;       s_)))

(defun principal-strains (strains)
  (multiple-value-bind (l v) (cl-mpm/utils::eig (voigt-to-matrix strains))
    (declare (ignore v))
    (setf l (sort l #'>))
    (values (nth 0 l) (nth 1 l) (nth 2 l))))

(defun criterion-mohr-coloumb (strain angle E nu)
  (multiple-value-bind (e1 e2 e3) (principal-strains strain)
    (let ((a (/ (sin angle)
                (- 1d0 (* 2 nu))))
          (G (/ E (* 2 (- 1d0 nu))))
          )
      (* G
         (- (* (- 1d0 a) e3)
            (* (+ 1d0 a) e1))
         (/ 1d0 (cos angle)))
      ;; (/
      ;;  (* E (- (* (- 1d0 a) e1)
      ;;          (* (+ 1d0 a) e3)
      ;;          )
      ;;     (/ 1d0 (cos angle)))
      ;;  (* (- 1d0 a) (/ 1d0 (cos angle))))
      )))

(defun criterion-mohr-coloumb-3d (stress angle)
  (multiple-value-bind (s1 s2 s3) (principal-stresses-3d stress)
    (declare (double-float angle s1 s2 s3))
    (let ((tang (tan angle))
          (cang (cos angle)))
      (* 0.5d0
         (max
          (- (* (+ s3 s1) tang) (/ (- s3 s1) cang))
          ;; (- (- (* (+ s2 s3) tang) (/ (abs (- s2 s3)) cang)))
          ;; (- (- (* (+ s1 s2) tang) (/ (abs (- s1 s2)) cang)))
          )))))

(defun criterion-mohr-coloumb-stress (stress angle)
  (multiple-value-bind (s1 s2 s3) (principal-stresses-3d stress)
    (declare (double-float angle s1 s2 s3))
    (let ()
      (* 0.5d0
         (-
          (* (+ s3 s1) (tan angle))
          (/ (- s3 s1) (cos angle)))))))

(defun criterion-mohr-coloumb-stress-tensile (stress angle)
  (declare (double-float angle))
  (multiple-value-bind (s1 s2 s3) (principal-stresses-3d stress)
    (declare (double-float angle s1 s2 s3))
    (let ((k (/ (+ 1d0 (sin angle))
                (- 1d0 (sin angle)))))
      (/ (- (* k s1) s3)
         k)))
  )
(defun criterion-mohr-coloumb-stress-will (stress angle)
  (multiple-value-bind (s1 s2 s3) (principal-stresses-3d stress)
    (declare (double-float angle s1 s2 s3))
    (let ((k (/ (+ 1d0 (sin angle))
                (- 1d0 (sin angle))
                )))
      (* 0.5d0
         (/ (- (* k s1) s3)
            (sqrt k))))))

(defun mohr-coloumb-coheasion-to-tensile (coheasion angle)
  (let (
        ;; (a (/ (sin angle)
        ;;       (- 1d0 (* 2 nu))))
        ;; (G (/ E (* 2 (- 1d0 nu))))
        (k (/ (+ 1d0 (sin angle))
              (- 1d0 (sin angle))
              )))
    ;; (/ (* 2 coheasion) (sqrt k))
    (the double-float
         (/
          (* 2 coheasion (cos angle))
          (+ 1d0 (sin angle))))
    ;; (* E
    ;;    0.5d0
    ;;    (/
    ;;     (/ coheasion G)
    ;;     (* (- 1d0 a) (/ 1d0 (cos angle)))))
    ))


;; (let* ((stress (cl-mpm/utils:voigt-from-list (list 1d0 0d0 0d0
;;                                                    0d0 0d0 1d0)))
;;        (angle 42d0))
;;   (print (cl-mpm/damage::criterion-mohr-coloumb-stress stress (* angle (/ pi 180))))
;;   (print (cl-mpm/damage::criterion-mohr-coloumb-3d stress (* angle (/ pi 180))))
;;   ;; (print (cl-mpm/damage::criterion-mohr-coloumb-stress-will stress (* angle (/ pi 180))))
;;   ;; (print (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile stress (* angle (/ pi 180))))
;;   )

;;   (let ((mc (cl-mpm/damage::criterion-mohr-coloumb-stress stress (* angle (/ pi 180)))))
;;     (format t "MC: ~%~E~%" mc)
;;     (format t "Insc: ~E~%"  (/ (cl-mpm/damage::criterion-dp-inscribe-coheasion stress (* angle (/ pi 180))) mc))
;;     (format t "Shear: ~E~%"  (/ (cl-mpm/damage::criterion-dp-shear-coheasion stress (* angle (/ pi 180))) mc))
;;     (format t "Mid: ~E~%"  (/ (cl-mpm/damage::criterion-dp-middle-circumscribe-coheasion stress (* angle (/ pi 180))) mc))
;;     (format t "Circ: ~E~%"  (/ (cl-mpm/damage::criterion-dp-coheasion stress (* angle (/ pi 180))) mc)))
;;   )

(defun volumetric-deviatoric-norm (strain E nu
                                   res-t res-c res-s)
  (let* ((etr (cl-mpm/utils:trace-voigt strain))
         (e-tension (max 0d0 etr))
         (e-compression (- (max -0d0 (- etr))))
         (ed (cl-mpm/constitutive::deviatoric-voigt strain))
         (K (/ E (* 3d0 (- 1d0 (* 2d0 nu)))))
         (G (/ E (* 2 (+ 1d0 nu))))
         (energy
           (* 0.5d0
              (+
               (* res-t K (expt e-tension 2))
               (* res-c K (expt e-compression 2))
               (*
                G
                (cl-mpm/fastmaths::dot ed ed))))))
    (sqrt (max 0d0 (* E energy)))
    ))

(defun tensile-energy-norm (strain E de)
  (let* ((strain+
           (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix strain))
             (loop for i from 0 to 2
                   do
                      (setf (nth i l) (max (nth i l) 0d0)))
             (cl-mpm/utils:matrix-to-voigt
              (magicl:@ v
                        (cl-mpm/utils::matrix-from-diag l)
                        (magicl:transpose v))))))
    (sqrt (max 0d0 (* E (cl-mpm/fastmaths::dot strain+ (magicl:@ de strain+)))))))

(defun tensile-energy-norm-pressure (strain E nu de pressure)
  (let* ((K (/ E (* 3d0 (- 1d0 (* 2d0 nu)))))
         (ep (* -1/3 (/ pressure K)))
         (strain+
           (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix strain))
             (loop for i from 0 to 2
                   do
                      (setf (nth i l) (max (+ (nth i l) ep) 0d0)))
             (cl-mpm/utils:matrix-to-voigt
              (magicl:@
               v
               (cl-mpm/utils::matrix-from-diag l)
               (magicl:transpose v))))))
    (sqrt (max 0d0 (* E (cl-mpm/fastmaths::dot strain+ (magicl:@ de strain+)))))))

(defun criterion-effective-principal-stress (stress pressure)
  (multiple-value-bind (s1 s2 s3) (principal-stresses-3d stress)
    (- (max s1 s2 s3) pressure))
  ;; (let* ((K (/ E (* 3d0 (- 1d0 (* 2d0 nu)))))
  ;;        (ep (* -1/3 (/ pressure K)))
  ;;        (strain+
  ;;          (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix strain))
  ;;            (loop for i from 0 to 2
  ;;                  do
  ;;                     (setf (nth i l) (max (+ (nth i l) ep) 0d0)))
  ;;            (cl-mpm/utils:matrix-to-voigt
  ;;             (magicl:@
  ;;              v
  ;;              (cl-mpm/utils::matrix-from-diag l)
  ;;              (magicl:transpose v))))))
  ;;   (sqrt (max 0d0 (* E (cl-mpm/fastmaths::dot strain+ (magicl:@ de strain+))))))
  )

(defun criterion-enhanced-bi-energy-norm (stress strain E nu k)
  (let* ((strain+
           (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix strain))
             (loop for i from 0 to 2
                   do
                      (setf (nth i l) (max (nth i l) 0d0)))
             (cl-mpm/utils:matrix-to-voigt
              (magicl:@ v
                        (cl-mpm/utils::matrix-from-diag l)
                        (magicl:transpose v)))))

         (i1 (cl-mpm/utils:trace-voigt strain))
         (j2
           (cl-mpm/constitutive::voigt-j2
            (cl-mpm/utils::deviatoric-voigt
             (cl-mpm/utils::voigt-contra->covar strain))))
         (j2-factor (sqrt (* 3d0 j2)))

         (et (+ (/ i1 (* 2d0 (- 1d0 (* 2d0 nu))))
                (/ j2-factor (* 2d0 (+ 1d0 nu)))
                ))
         (ec (+ (/ i1 (* 5d0 (- 1d0 (* 2d0 nu))))
                (/ (* 6d0 j2-factor) (* 5d0 (+ 1d0 nu)))
                ))
         ;; (r
         ;;   (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix stress))
         ;;     (loop for i from 0 to 2
         ;;           do
         ;;              (setf (nth i l) (max (nth i l) 0d0)))
         ;;     (/ )(reduce #'+ (mapcar (lambda (x) (max 0d0 x)) l))
         ;;     )
         ;;   )
         )
    ))

(defun criterion-mc (strain angle E nu)
  (multiple-value-bind (e_1 e_2 e_3) (principal-stresses-3d
                                      ;; strain
                                      (magicl:scale!
                                       ;(cl-mpm/utils::voigt-contra->covar strain)
                                       ;(cl-mpm/utils::voigt-contra->covar strain)
                                       strain
                                       -1d0
                                       )
                                      )
    (let ((a (/ (sin angle) (- 1d0 (* 2d0 nu)))))
      (* E
         (/ 1d0 (cos angle))
         (-
          (* (- 1d0 a) e_1)
          (* (+ 1d0 a) e_2)
          )))))




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


(declaim (ftype (function (double-float double-float double-float double-float )
                          double-float) damage-response-exponential))
(defun damage-response-exponential (stress E init-stress ductility)
  (declare (double-float stress E init-stress ductility))
  "Function that controls how damage evolves with principal stresses"
  (let* ((ft init-stress)
         (e0 (/ ft E))
         (ef (/ (* ft (+ ductility 1d0)) (* 2d0 E)))
         (k (/ stress E))
         (beta (/ 1d0 (- ef e0)))
         )
    (declare (double-float ft e0 ef k beta))
    (if (> k e0)
        (- 1d0 (* (/ e0 k) (exp (- (* beta (- k e0))))))
        0d0)))

(defun find-k-damage-mp (mp damage)

  (with-accessors
        ((E cl-mpm/particle::mp-e)
         (init-stress cl-mpm/particle::mp-initiation-stress)
         (ductility cl-mpm/particle::mp-ductility))
      mp
    (find-k-damage E init-stress ductility damage))
  )
(defun find-k-damage (E init-stress ductility damage)
  (let* ((e0 init-stress)
         (d-est 0d0)
         (k-est e0)
         (k-max nil)
         (k-prev e0)
         )
    (loop for i from 0 to 10000
          while (> (abs (- damage d-est)) 1d-8)
          do
             (progn
               (setf d-est (damage-response-exponential k-est E init-stress ductility))
               ;; (format t "~F ~F~%" d-est k-est)
               (when (> damage d-est)
                 (setf k-prev k-est)
                 (setf k-est
                       (if k-max
                           (/ (+ k-prev k-max) 2d0)
                           (* k-est 2))))
               (when (< damage d-est)
                 (setf k-max k-est)
                 (setf k-est (/ (+ k-prev k-max) 2d0)))
               ))
    k-est))
;; (let* ((E 1d9)
;;        (in 20d3)
;;        (duct 100d0)
;;        (k-est (cl-mpm/damage::find-k-damage E in duct 0.9999d0)))
;;   (pprint (damage-response-exponential k-est E in duct)))

;; (defun get-k-damage (mp damage)
;;   (with-accessors
;;         ((E cl-mpm/particle::mp-e)
;;          (init-stress cl-mpm/particle::mp-initiation-stress)
;;          (ductility cl-mpm/particle::mp-ductility))
;;       mp
;;     ;; (declare (double-float stress E init-stress ductility))
;;     (let* ((ft init-stress)
;;            (e0 (/ ft E))
;;            (ef (/ (* ft (+ ductility 1d0)) (* 2d0 E)))
;;            ;; (k (/ stress E))
;;            (beta (/ 1d0 (- ef e0)))
;;            )
;;       (if (> damage 0d0)

;;           0d0)
;;       (if (> k e0)
;;           (- 1d0 (* (/ e0 k) (exp (- (* beta (- k e0))))))
;;           0d0))))


(declaim (ftype (function (double-float double-float double-float double-float double-float)
                          double-float) damage-response-exponential-peerlings-residual))
(defun damage-response-exponential-peerlings-residual (stress E init-stress ductility residual)
  (declare (double-float stress E init-stress ductility residual))
  "Function that controls how damage evolves with principal stresses"
  (let* ((ft init-stress)
         (e0 (/ ft E))
         (ef (/ (* ft (+ ductility 1d0)) (* 2d0 E)))
         (k (/ stress E))
         (beta (/ 1d0 (- ef e0))))
    (declare (double-float ft e0 ef k beta))
    (if (> k e0)
        (- 1d0 (* (/ e0 k) (+ (- 1d0 residual) (* residual (exp (- (* beta (- k e0))))))))
        0d0)))
(declaim (ftype (function (double-float double-float double-float double-float double-float)
                          double-float) damage-response-exponential-residual))
(defun damage-response-exponential-residual (stress E init-stress ductility residual)
  (declare (double-float stress E init-stress ductility residual))
  "Function that controls how damage evolves with principal stresses"
  (let* ((ft init-stress)
         (e0 (/ ft E))
         (ef (/ (* ft (+ ductility 1d0)) (* 2d0 E)))
         (k (/ stress E))
         (beta (/ 1d0 (- ef e0))))
    (declare (double-float ft e0 ef k beta))
    (if (> k e0)
        (- 1d0 (+ (- 1d0 residual) (* residual (exp (- (* beta (- k e0)))))))
        0d0)))

(defun damage-response-exponential-beta (E init-stress ductility)
  (declare (double-float E init-stress ductility))
  "Function that controls how damage evolves with principal stresses"
  (let* ((ft init-stress)
         (e0 (/ ft E))
         (ef (/ (* ft (+ ductility 1d0)) (* 2d0 E)))
         (beta (/ 1d0 (- ef e0))))
    beta))

;; (defun damage-response-exponential-inverse (damage E length init-stress ductility)
;;   (declare (double-float stress E length init-stress ductility))
;;   "Function that controls how damage evolves with principal stresses"
;;   (let* ((ft init-stress)
;;          (e0 (/ ft E))
;;          (ef (/ (* ft (+ ductility 1d0)) (* 2d0 E)))
;;          (k (/ stress E))
;;          (beta (/ 1d0 (- ef e0)))
;;          )
;;     (if (> k e0)
;;         (- 1d0 (* (/ e0 k) (+ 0d0 (exp (- (* beta (- k e0)))))))
;;         0d0))
;;   )


(defun modified-vm-stress (stress k init-stress)
  (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d stress)
    (let ((i1 (+ s_1 s_2 s_3)))
      (* 0.5d0
         (+
          (* (/ 1d0 (* k init-stress))
             (+
              (expt (- s_1 s_2) 2)
              (expt (- s_2 s_3) 2)
              (expt (- s_3 s_1) 2))
             )
          (* 2d0 (- 1 (/ 1d0 k)) i1)))
      )))
;; (defun modified-vm-strain (stress k init-stress nu)
;;   (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d stress)
;;     (let (
;;           (j2 (cl-mpm/constitutive::voigt-j2 (cl-mpm/constitutive::deviatoric-voigt
;;                                               stress)))
;;           (i1 (+ s_1 s_2 s_3))
;;           (k-factor (/ (- k 1d0)
;;                        (- 1d0 (* 2d0 nu)))))
;;       (+ (* i1 (/ k-factor (* 2d0 k)))
;;          (* (/ 1d0 (* 2d0 k))
;;             (sqrt (+ (expt (* k-factor i1) 2)
;;                      (* (/ (* 12 k) (expt (- 1d0 nu) 2)) j2)
;;                      )))))))
(defun criterion-modified-vm (strain k nu)
  (multiple-value-bind (s_1 s_2 s_3) (principal-stresses-3d strain)
    (let ((i1 (+ s_1 s_2 s_3))
          (j2 (cl-mpm/constitutive::voigt-j2 (cl-mpm/constitutive::deviatoric-voigt strain)))
          (k-factor (/ (- k 1d0)
                       (- 1d0 (* 2d0 nu))))
          )
      (* (/ 1d0 (* 2d0 k))
         (+ (* i1 k-factor)
            (sqrt (+ (expt (* k-factor i1) 2)
                     (* (/ (* 12d0 k) (expt (- 1d0 nu) 2)) j2))))))))





(defun estimate-ductility-jirsek2004 (GF R ft E &optional (k 1d0))
  "Simple estimate of ductility from key elastic parameters"
  (declare (double-float GF R ft E k))
  (let* ((e0 (/ ft E))
         (ef (+ (/ GF (* k R E e0)) (/ e0 2d0))))
    (- (* 2d0 (/ ef e0)) 1d0)))

(defun gf-from-ductility (ductility R ft E &optional (k 1d0))
  (declare (double-float ductility R ft E k))
  "Sanity check for estimating fracture energy from ductility and elastic parameters"
  (let ((e0 (/ ft E)))
    (/ (* ductility R E (expt e0 2))
       2d0)))




