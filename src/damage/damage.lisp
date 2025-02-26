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
    :initform 10)
   (enable-length-localisation
    :type boolean
    :accessor sim-enable-length-localisation
    :initarg :enable-length-localisation
    :initform nil))
  (:documentation "Explicit simulation with update stress first update"))

(defclass mpm-sim-usl-damage (mpm-sim-damage cl-mpm::mpm-sim-usl) ())

(defclass mpm-sim-damage-nd-2 (mpm-sim-damage cl-mpm::mpm-nd-2d) ())

(defgeneric cl-mpm/particle::post-damage-step (mp dt))
(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle) dt))
(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-damage) dt))

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
(defun find-intergral-local-length (mesh mp)
  (with-accessors ((tll cl-mpm/particle::mp-true-local-length)
                   (ll cl-mpm/particle::mp-local-length)
                   (lld cl-mpm/particle::mp-local-length-damaged))
      mp
    (setf tll
          (length-localisation
           ll
           lld
           (calculate-average-damage mesh mp ll)))))



(defun calculate-damage (sim)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)
                   (dt cl-mpm:sim-dt)
                   (enable-damage cl-mpm::sim-enable-damage)
                   (delocal-counter sim-damage-delocal-counter)
                   (delocal-counter-max sim-damage-delocal-counter-max)
                   (non-local-damage cl-mpm::sim-nonlocal-damage))
      sim

    (when enable-damage
      (when non-local-damage
        (when (<= delocal-counter 0)
          ;;Ensure they have a home
          ;; (create-delocalisation-list mesh mps)
          (update-delocalisation-list mesh mps)
          (setf delocal-counter delocal-counter-max))
        ;;Worst case we want dont want our delocal counter to exceed the max incase we get adjusted
        (setf delocal-counter (min (- delocal-counter 1)
                                   delocal-counter-max)))
      (cl-mpm:iterate-over-mps
       mps
       (lambda (mp)
         (when (typep mp 'cl-mpm/particle:particle-damage)
           (damage-model-calculate-y mp dt)
           )))

      (if non-local-damage
          (progn
            (when (sim-enable-length-localisation sim)
              (update-localisation-lengths sim))
            (delocalise-damage sim))
          (localise-damage mesh mps dt))
      (cl-mpm:iterate-over-mps
       mps
       (lambda (mp)
         (when (typep mp 'cl-mpm/particle:particle-damage)
           (update-damage mp dt)))))
    (cl-mpm:iterate-over-mps
     mps
     (lambda (mp)
       (when (typep mp 'cl-mpm/particle:particle-damage)
         (cl-mpm/particle::post-damage-step mp dt))))))

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
    (cl-mpm:iterate-over-nodes
     mesh
     (lambda (node)
       (when (= 0 (length (cl-mpm/mesh::node-local-list node)))
         (adjust-array (cl-mpm/mesh::node-local-list node) 0 :fill-pointer 0))))
    ))

(defun local-list-add-particle (mesh mp)
  "A function for putting an MP into the nodal MP lookup table"
  (let ((node-id (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp))))
    (if (cl-mpm/mesh:in-bounds mesh node-id)
      (let ((node (cl-mpm/mesh:get-node mesh node-id)))
        (setf (cl-mpm/particle::mp-damage-position mp)
              (cl-mpm/utils::vector-copy (cl-mpm/particle:mp-position mp)))
        (sb-thread:with-mutex ((cl-mpm/mesh:node-lock node))
          (vector-push-extend mp (cl-mpm/mesh::node-local-list node))))
      (error "Could not add damage MP to local list as it is out of bounds"))))

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
                  (cl-mpm::iterate-over-nodes
                   mesh
                   (lambda (node)
                     (remove-mp-ll node))))
              t)
            nil))
        nil)))

(declaim (notinline update-delocalisation-list))
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
               ;; (format t "~A ~A~%" (cl-mpm/particle:mp-position mp) (cl-mpm/particle::mp-damage-position mp))
               (when (> delta (/ h 4d0))
                 (when (not (equal
                             (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp))
                             (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle::mp-damage-position mp))))
                   ;; (pprint "Moving MP")
                   (local-list-remove-particle mesh mp)
                   (local-list-add-particle mesh mp))))))))
    ;;Smart in thought, presumably creates lots of garbage
    ;; (cl-mpm::iterate-over-nodes
    ;;  mesh
    ;;  (lambda (node)
    ;;    (when (= 0 (length (cl-mpm/mesh::node-local-list node)))
    ;;      (adjust-array (cl-mpm/mesh::node-local-list node) 0 :fill-pointer 0))))
    ))

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

(defun calculate-average-damage (mesh mp length)
  (let ((damage-average 0d0)
        (volume-average 0d0))
    (declare (double-float damage-average))
    (iterate-over-neighour-mps
     mesh mp length
     (lambda (mesh mp mp-other dist)
       (with-accessors ((d cl-mpm/particle::mp-damage)
                        (m cl-mpm/particle:mp-volume)
                        (ll cl-mpm/particle::mp-true-local-length)
                        (p cl-mpm/particle:mp-position))
           mp-other
         (when t
           (let ((weight (weight-func-mps mesh mp mp-other length)))
             (declare (double-float weight m d damage-average volume-average))
             (incf volume-average (* weight m))
             (incf damage-average (* d weight m)))))))
    (when (> volume-average 0d0)
      (setf damage-average (/ damage-average volume-average)))
    damage-average))

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
                   ;; (weight-func-mps mesh mp mp-other length)
                   ;; (weight-func-mps mesh mp mp-other (* 0.5d0 (+ length ll)))
                   ;; (weight-func-mps mesh mp mp-other ll)
                   ;; (weight-func-mps mesh mp mp-other ll)
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
  ;; (* local-length (max (/ (- (exp (- damage)) (exp -1d0))
  ;;                         (- 1d0 (exp -1d0))) 1d-10))
  )
(defun find-mp-local-length (mesh mp)
  (with-accessors ((tll cl-mpm/particle::mp-true-local-length)
                   (ll cl-mpm/particle::mp-local-length)
                   (lld cl-mpm/particle::mp-local-length-damaged)
                   (damage cl-mpm/particle::mp-damage))
      mp
    (setf tll
          (length-localisation
           ll
           lld
           damage))))

(defgeneric update-localisation-lengths (sim))
(defmethod update-localisation-lengths ((sim cl-mpm/damage::mpm-sim-damage))
  (with-accessors ((mps cl-mpm::sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (cl-mpm:iterate-over-mps
     mps
     (lambda (mp)
       (when (typep mp 'cl-mpm/particle:particle-damage)
         (find-intergral-local-length mesh mp)
         ;; (find-mp-local-length mesh mp)
         )))))

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
         ;; (find-nodal-local-length mesh mp)
         ;; (find-intergral-local-length mesh mp)
         ;; (setf (cl-mpm/particle::mp-true-local-length mp)
         ;;       (length-localisation (cl-mpm/particle::mp-local-length mp)
         ;;                            (cl-mpm/particle::mp-local-length-damaged mp)
         ;;                            (cl-mpm/particle::mp-damage mp)))
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
  (when (typep mp 'cl-mpm/particle::particle-damage)
    (local-list-add-particle (cl-mpm:sim-mesh sim) mp))
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
               (time cl-mpm::time)
               )
                sim
    (declare (type double-float mass-filter))
                (progn
                    (cl-mpm::reset-grid mesh)
                    (cl-mpm::p2g mesh mps)
                    (when (> mass-filter 0d0)
                      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
                    (cl-mpm::update-node-kinematics mesh dt)
                    (cl-mpm::apply-bcs mesh bcs dt)

                    (cl-mpm::update-stress mesh mps dt fbar)

                    (cl-mpm/damage::calculate-damage sim)
                    ;; ;Map forces onto nodes
                    (cl-mpm::p2g-force mesh mps)

                    (loop for bcs-f in bcs-force-list
                          do (cl-mpm::apply-bcs mesh bcs-f dt))

                    (cl-mpm::update-node-forces sim)
                    ;; ;Reapply velocity BCs
                    (cl-mpm::apply-bcs mesh bcs dt)
                    ;; ;Also updates mps inline
                    (cl-mpm::g2p mesh mps dt vel-algo)
                    (cl-mpm::update-dynamic-stats sim)

                    (when remove-damage
                      (cl-mpm::remove-material-damaged sim))
                    (when split
                      (cl-mpm::split-mps sim))
                    (cl-mpm::check-mps sim)
                    (cl-mpm::check-single-mps sim)
                    (incf time dt)
                    )))

(defmethod cl-mpm::update-sim ((sim mpm-sim-usl-damage))
  (declare (cl-mpm::mpm-sim-usl sim))
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
               (time cl-mpm::time)
               )
      sim
    (declare (type double-float mass-filter))
    (progn
      (cl-mpm::check-single-mps sim)
      (cl-mpm::reset-grid mesh)
      (cl-mpm::p2g mesh mps)
      ;;Do optional mass filter
      (when (> mass-filter 0d0)
        (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
      (cl-mpm::update-node-kinematics mesh dt)
      (cl-mpm::apply-bcs mesh bcs dt)
                                        ;Map forces onto nodes
      (cl-mpm::p2g-force mesh mps)
      (loop for bcs-f in bcs-force-list
            do (cl-mpm::apply-bcs mesh bcs-f dt))
      (cl-mpm::update-node-forces sim)
                                        ;Reapply velocity BCs
      (cl-mpm::apply-bcs mesh bcs dt)
                                        ;Also updates mps inline
      (cl-mpm::update-dynamic-stats sim)
      (cl-mpm::g2p mesh mps dt vel-algo)

      ;;Update stress last
      (cl-mpm::reset-grid-velocity mesh)
      (cl-mpm::p2g mesh mps)
      ;; (cl-mpm::check-single-mps sim)
      ;;Do optional mass filter
      (when (> mass-filter 0d0)
        (cl-mpm::filter-grid-velocity mesh (cl-mpm::sim-mass-filter sim)))
      (cl-mpm::update-node-kinematics mesh dt)
      (cl-mpm::apply-bcs mesh bcs dt)

      (cl-mpm::update-stress mesh mps dt fbar)
      (cl-mpm/damage::calculate-damage sim)
      (when remove-damage
        (cl-mpm::remove-material-damaged sim))
      (when split
        (cl-mpm::split-mps sim))
      (cl-mpm::check-mps sim)
      (incf time dt)
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
        (cl-mpm/output::save-parameter "unique-id" (cl-mpm/particle::mp-unique-index mp))
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
        (cl-mpm/output::save-parameter "eps_1"
                        (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-strain mp)))
                          (loop for sii in l maximize sii)))
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
        (cl-mpm/output::save-parameter "erosion"
                                       (if (slot-exists-p mp 'cl-mpm/particle::eroded-volume)
                                           (cl-mpm/particle::mp-eroded-volume mp)
                                           0))
        (cl-mpm/output::save-parameter "pressure" (cl-mpm/particle::mp-pressure mp))
        (cl-mpm/output::save-parameter "boundary" (cl-mpm/particle::mp-boundary mp))

        (cl-mpm/output::save-parameter "s_1"
                                       (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))

                                         (nth 0 (sort l #'>))))
        (cl-mpm/output::save-parameter "s_3"
                                       (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))

                                         (nth 2 (sort l #'>))))

        (cl-mpm/output::save-parameter "su_1"
                                       (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle::mp-undamaged-stress mp)))

                                         (nth 0 (sort l #'>))))
        (cl-mpm/output::save-parameter "su_3"
                                       (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle::mp-undamaged-stress mp)))

                                         (nth 2 (sort l #'>))))
        ;; (cl-mpm/output::save-parameter "e_1"
        ;;                                (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-strain mp)))
        ;;                                  (loop for sii in l maximize sii)))


        ;; (cl-mpm/output::save-parameter "EPS"
        ;;                                (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))
        ;;                                  (- (loop for sii in l maximize sii) (cl-mpm/particle::mp-pressure mp))))
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
        (cl-mpm/output::save-parameter "local-length"
                                       (if (slot-exists-p mp 'cl-mpm/particle::true-local-length)
                                           (cl-mpm/particle::mp-true-local-length mp)
                                           0d0))
        (cl-mpm/output::save-parameter "p-undamaged"
                                       (if (slot-exists-p mp 'cl-mpm/particle::undamaged-stress)
                                           (cl-mpm/utils::trace-voigt (cl-mpm/particle::mp-undamaged-stress mp))
                                           (cl-mpm/utils::trace-voigt (cl-mpm/particle::mp-stress mp))))
        (cl-mpm/output::save-parameter "p-damaged"
                                       (if (slot-exists-p mp 'cl-mpm/particle::undamaged-stress)
                                           (cl-mpm/utils::trace-voigt (cl-mpm/particle::mp-stress mp))
                                           0d0))

        (cl-mpm/output::save-parameter "q-undamaged"
                                       (if (slot-exists-p mp 'cl-mpm/particle::undamaged-stress)
                                           (sqrt (cl-mpm/constitutive::voigt-j2 (cl-mpm/utils:deviatoric-voigt (cl-mpm/particle::mp-undamaged-stress mp))))
                                           (sqrt (cl-mpm/constitutive::voigt-j2 (cl-mpm/utils:deviatoric-voigt (cl-mpm/particle::mp-stress mp))))))
        (cl-mpm/output::save-parameter "q-damaged"
                                       (if (slot-exists-p mp 'cl-mpm/particle::undamaged-stress)
                                           (sqrt (cl-mpm/constitutive::voigt-j2 (cl-mpm/utils:deviatoric-voigt (cl-mpm/particle::mp-stress mp))))
                                           0d0))

        (cl-mpm/output::save-parameter "split-depth"
                                       (cl-mpm/particle::mp-split-depth mp))

        (cl-mpm/output::save-parameter "plastic-iterations" (cl-mpm/particle::mp-plastic-iterations mp))
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

        ;; (cl-mpm/output::save-parameter
        ;;  "plastic-c"
        ;;  (if (slot-exists-p mp 'cl-mpm/particle::c)
        ;;      (cl-mpm/particle::mp-c mp)
        ;;      0d0))
        ;; (cl-mpm/output::save-parameter
        ;;  "plastic-phi"
        ;;  (if (slot-exists-p mp 'cl-mpm/particle::phi)
        ;;      (* (cl-mpm/particle::mp-phi mp) (/ 180 pi))
        ;;      0d0))

        (cl-mpm/output::save-parameter
         "energy"
         ;; (cl-mpm/particle::mp-penalty-energy mp)
         (* (cl-mpm/particle::mp-mass mp)
            (cl-mpm/fastmaths::mag-squared (cl-mpm/particle::mp-velocity mp)))
         )

        ;; (cl-mpm/output::save-parameter
        ;;  "fbar-j"
        ;;  (cl-mpm/particle::mp-debug-j mp)
        ;;  )
        ;; (cl-mpm/output::save-parameter
        ;;  "fbar-j-gather"
        ;;  (cl-mpm/particle::mp-debug-j-gather mp)
        ;;  )
        ;; (cl-mpm/output::save-parameter
        ;;  "fbar-j-diff"
        ;;  (if (> (cl-mpm/particle::mp-debug-j mp) 0d0) 
        ;;      (/ (- (cl-mpm/particle::mp-debug-j-gather mp) (cl-mpm/particle::mp-debug-j mp)) 
        ;;         (cl-mpm/particle::mp-debug-j mp)
        ;;         )
        ;;      0d0)
        ;;  )
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














;; (defun solve-quasi-static ()
;;   (let ((damage-0 0d0)
;;         (damage-1 0d0))
;;     (loop for i from 0 to 100
;;           while (or (= i 0) (> damage-1 (* 1.1d0 damage-0)))
;;           do
;;              (progn
;;                (setf damage-0 damage-1)
;;                (format t "Staggered solve: ~D~%" i)
;;                (cl-mpm/dynamic-relaxation:converge-quasi-static
;;                 *sim*
;;                 :energy-crit 1d-2
;;                 :oobf-crit 1d-2
;;                 :substeps 50
;;                 :conv-steps 200
;;                 :dt-scale dt-scale
;;                 :post-iter-step
;;                 (lambda (i energy oobf)
;;                   (format t "Surcharge load ~E~%" (/ *piston-confinement* *piston-steps*))
;;                   (setf *piston-confinement* 0d0
;;                         *piston-steps* 0)))
;;                (cl-mpm/damage::calculate-damage *sim*)
;;                ;; (plot *sim*)
;;                (setf damage-1 (get-damage))
;;                ;; (plot-domain)
;;                (format t "Staggered damage diff: ~E~%" (- damage-1 damage-0))
;;                ))
;;     (when (> damage-1 (* damage-0 1.1d0))
;;       (break "Staggered error"))))
