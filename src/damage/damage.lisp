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
    :initform nil)
   (enable-ekl
    :type boolean
    :accessor sim-enable-ekl
    :initarg :enable-ekl
    :initform nil)
   )
  (:documentation "Explicit simulation with update stress first update"))

(defclass mpm-sim-usl-damage (mpm-sim-damage cl-mpm::mpm-sim-usl) ())

(defclass mpm-sim-damage-nd-2 (mpm-sim-damage cl-mpm::mpm-nd-2d) ())

(defclass mpm-sim-agg-damage (mpm-sim-damage cl-mpm/aggregate::mpm-sim-aggregated)
  ()
  (:documentation "Explicit damage simulation"))


(defgeneric set-mp-damage (mp d))
(defmethod set-mp-damage ((mp cl-mpm/particle::particle-damage) d)
  (setf
   (cl-mpm/particle::mp-damage mp) d
   (cl-mpm/particle::mp-damage-n mp) d))

(defmethod cl-mpm::reset-loadstep ((sim mpm-sim-damage))
  (call-next-method)
  ;; (calculate-damage sim 1d0)
  )

(defmethod cl-mpm/particle::reset-loadstep-mp ((mp cl-mpm/particle::particle-damage))
  (with-accessors ((ybar      cl-mpm/particle::mp-damage-ybar)
                   (ybar-prev cl-mpm/particle::mp-damage-ybar-prev)
                   (y         cl-mpm/particle::mp-damage-y-local)
                   ;; (y-prev    cl-mpm/particle::mp-damage-y-local-prev)
                   (damage    cl-mpm/particle::mp-damage)
                   (damage-n  cl-mpm/particle::mp-damage-n)
                   (stress cl-mpm/particle::mp-undamaged-stress))
      mp
    (setf ybar ybar-prev)
    (setf damage damage-n)
    (cl-mpm/fastmaths::fast-zero stress)
    (cl-mpm/damage::compute-damage mp)
    (call-next-method)))

(defgeneric cl-mpm/particle::post-damage-step (mp dt))
(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle) dt))

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

;; (declaim
;;  ;(inline calculate-damage-increment)
;;  (ftype (function (cl-mpm/particle:particle double-float) t) calculate-damage-increment))
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
  )
(defun find-nodal-local-length (mesh mp)
  (let ((nodal-damage 0d0))
    (cl-mpm::iterate-over-neighbours
     mesh mp
     (lambda (mesh mp node svp grads fsvp fgrads)
       (with-accessors ((node-damage cl-mpm/mesh::node-damage)
                        (node-svp cl-mpm/mesh::node-svp-sum))
           node
         (incf nodal-damage (* node-damage svp)))
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


(defgeneric calculate-damage (sim dt))

(defmethod calculate-damage ((sim mpm-sim-damage) dt)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)
                   (enable-damage cl-mpm::sim-enable-damage)
                   (delocal-counter sim-damage-delocal-counter)
                   (delocal-counter-max sim-damage-delocal-counter-max)
                   (non-local-damage cl-mpm::sim-nonlocal-damage))
      sim
    (when enable-damage
      (when non-local-damage
        ;;When delocal counter max is -1, we want to update never -> 0 is always
        (unless (= delocal-counter-max -1)
          (when (= delocal-counter 0)
            ;;Ensure they have a home
            ;; (create-delocalisation-list mesh mps)
            (update-delocalisation-list mesh mps)
            (setf delocal-counter delocal-counter-max)))
        ;;Worst case we want dont want our delocal counter to exceed the max incase we get adjusted
        (setf delocal-counter (min (- delocal-counter 1) delocal-counter-max)))

      (cl-mpm:iterate-over-mps
       mps
       (lambda (mp)
         (when (typep mp 'cl-mpm/particle:particle-damage)
           (damage-model-calculate-y mp dt))))

      (if non-local-damage
          (progn
            (when (sim-enable-length-localisation sim)
              (update-localisation-lengths sim))
            (if (cl-mpm/damage::sim-enable-ekl sim)
                (delocalise-damage-ekl sim)
                (delocalise-damage sim)))
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
 (notinline diff-damaged)
 (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle cl-mpm/particle:particle) double-float) diff-damaged))
(defun diff-damaged (mesh mp-a mp-b)
  (flet ((node-dam (node)
           (if (cl-mpm/mesh::node-active node)
               (min 1d0 (max 0d0 (the double-float (cl-mpm/mesh::node-damage node))))
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
                               (let ((damage 0d0))
                                 (declare (double-float damage))
                                 (cl-mpm::iterate-over-neighbours-point-linear
                                  mesh
                                  step-point
                                  (lambda (m node weight grads)
                                    (declare (double-float damage weight))
                                    ;; (pprint (cl-mpm/mesh::node-damage node))
                                    (incf damage
                                          (* (node-dam node)
                                             weight))))

                                 ;; (print damage)
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

                                        ; (print damage)
                      (incf final-distance (* step-size
                                              (/ 1d0
                                                 (max epsilon
                                                      (the double-float (sqrt (max 0d0 (min 1d0 (- 1d0 damage))))))))))
                    )
                  )
                )
              )
            )
          (max 0d0 (* final-distance final-distance)))
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
                        (p cl-mpm/particle:mp-position))
           mp-other
         (declare (double-float d m length))
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
                   (* (the double-float (cl-mpm/particle::mp-damage-y-local mp-other))
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
                                   (* (the double-float (cl-mpm/particle::mp-damage-y-local mp-other))
                                      weight m))))))
             (reflect-axis 0 *enable-reflect-x*)
             (reflect-axis 1 *enable-reflect-y*)
             (reflect-axis 2 *enable-reflect-z*)
             )
           ))
       ))
    (when (> mass-total 0d0)
      (setf damage-inc (/ damage-inc mass-total)))
    damage-inc))

(declaim (notinline length-localisation))
(defun length-localisation (local-length local-length-damaged damage)
  ;; (+ (* local-length (- 1d0 damage)) (* local-length-damaged damage))
  (* local-length (max (sqrt (- 1d0 damage)) 1d-10))
  ;; local-length
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

(defgeneric g2p-damage (sim)
  (:documentation
   "Map damage to grid"))

(defmethod g2p-damage ((sim cl-mpm/damage::mpm-sim-damage))
  (cl-mpm:iterate-over-nodes
   (cl-mpm:sim-mesh sim)
   (lambda (n)
     (when (cl-mpm/mesh::node-active n)
       (setf (cl-mpm/mesh::node-damage n) 0d0
             (cl-mpm/mesh::node-volume n) 0d0
             ))))
  (cl-mpm::iterate-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (when (typep mp 'cl-mpm/particle::particle-damage)
       (cl-mpm:iterate-over-neighbours
        (cl-mpm:sim-mesh sim)
        mp
        (lambda (mesh mp node svp grads fsvp fgrads)
          (sb-thread:with-mutex ((cl-mpm/mesh:node-lock node))
            (incf (cl-mpm/mesh::node-volume node)
                  (* svp
                     (cl-mpm/particle::mp-volume mp)))
            (incf (cl-mpm/mesh::node-damage node)
                  (* svp
                     (cl-mpm/particle::mp-volume mp)
                     (cl-mpm/particle::mp-damage mp)))))))))
  (cl-mpm:iterate-over-nodes
   (cl-mpm:sim-mesh sim)
   (lambda (n)
     (if (> (cl-mpm/mesh::node-volume n) 0d0);(cl-mpm/mesh::node-active n)
         (setf (cl-mpm/mesh::node-damage n)
               (min 1d0
                    (max 0d0 (/ (cl-mpm/mesh::node-damage n)
                                (cl-mpm/mesh::node-volume n)))))
       (setf (cl-mpm/mesh::node-damage n) 0d0)
       ;; (setf (cl-mpm/mesh::node-damage n) (/ (cl-mpm/mesh::node-damage n) (cl-mpm/mesh::node-svp-sum n)))
       ))))

(defgeneric update-localisation-lengths (sim))
(defmethod update-localisation-lengths ((sim cl-mpm/damage::mpm-sim-damage))
  (with-accessors ((mps cl-mpm::sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (g2p-damage sim)
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
         (with-accessors ((damage-ybar cl-mpm/particle::mp-damage-ybar)
                          (damage cl-mpm/particle::mp-damage)
                          (local-length-t cl-mpm/particle::mp-local-length))
             mp
           (setf damage-ybar (calculate-delocalised-damage mesh mp local-length-t)))))))
  (values))

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

(defun localise-damage (mesh mps dt)
  "Apply local damage model"
  (cl-mpm:iterate-over-mps
   mps
    (lambda (mp)
      (when (typep mp 'cl-mpm/particle:particle-damage)
        (with-accessors ((damage-ybar cl-mpm/particle::mp-damage-ybar)
                         (damage-y cl-mpm/particle::mp-damage-y-local))
            mp
          (setf damage-ybar damage-y)))))
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
            (let ((cauchy-undamaged (cl-mpm/fastmaths:fast-scale-voigt stress (/ 1d0 (cl-mpm/fastmaths:det-3x3 def)))))
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





(defgeneric compute-damage (mp))
(defmethod compute-damage ((mp cl-mpm/particle::particle))
  )
(defmethod compute-damage ((mp cl-mpm/particle::particle-damage))
  )

(defgeneric update-damage (mp dt))
(defmethod update-damage ((mp cl-mpm/particle::particle) dt))

(defmethod update-damage ((mp cl-mpm/particle::particle-damage) dt)
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
         (when t;(< (the double-float d) 1d0)
           (let (
                 ;;Nodally averaged local funcj
                 ;; (weight (weight-func-mps mesh mp mp-other (* 0.5d0 (+ length ll))))
                 ;; (weight (weight-func-mps mesh mp mp-other (sqrt (* length ll))))
                 ;;
                 (weight
                   ;; (weight-func-mps mesh mp mp-other (sqrt (* length ll)))
                   ;; (weight-func-mps mesh mp mp-other length)
                   ;; (weight-func-mps mesh mp mp-other (* 0.5d0 (+ length ll)))
                   ;; (weight-func-mps mesh mp mp-other ll)
                   ;; (weight-func-mps mesh mp mp-other ll)
                   (weight-func-mps-damaged mesh mp mp-other
                                            (cl-mpm/particle::mp-local-length mp)
                                            )
                   ))
             (declare (double-float weight m d mass-total damage-inc))
             (incf mass-total (* weight m))
             (incf damage-inc
                   (* (the double-float (cl-mpm/particle::mp-damage-y-local mp-other))
                      weight m)))
           ;; (macrolet ((reflect-axis (axis enable)
           ;;              (declare (fixnum axis))
           ;;              `(when (and ,enable
           ;;                          (< (magicl:tref (cl-mpm/particle::mp-position mp) ,axis 0) (* 0.5d0 length)))
           ;;                 (let ((weight (weight-func-pos
           ;;                                mesh
           ;;                                (cl-mpm/particle::mp-position mp)
           ;;                                (cl-mpm/fastmaths::fast-.* (cl-mpm/particle:mp-position mp-other)
           ;;                                                           (cl-mpm/utils::vector-from-list
           ;;                                                            (list ,(if (= axis 0) -1d0 0d0)
           ;;                                                                  ,(if (= axis 1) -1d0 0d0)
           ;;                                                                  ,(if (= axis 2) -1d0 0d0)
           ;;                                                                  )))
           ;;                                (sqrt (* length ll))
           ;;                                )))
           ;;                   (declare (double-float weight m d mass-total damage-inc))
           ;;                   (incf mass-total (* weight m))
           ;;                   (incf damage-inc
           ;;                         (* (the double-float (cl-mpm/particle::mp-damage-y-local mp-other))
           ;;                            weight m))))))
           ;;   (reflect-axis 0 *enable-reflect-x*)
           ;;   (reflect-axis 1 *enable-reflect-y*)
           ;;   (reflect-axis 2 *enable-reflect-z*)
           ;;   )
           ))
       ))
    (when (> mass-total 0d0)
      (setf damage-inc (/ damage-inc mass-total)))
    damage-inc))

(defun calculate-debug-ekl (mesh mp length)
  )

(defun get-nonlocal-interactions (sim mp)
  (g2p-damage sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps)) sim
    (let ((data-positions (list))
          (data-weights (list))
          (damage-inc 0d0)
          (mass-total 0d0)
          (length (* 1d0 (cl-mpm/particle::mp-local-length mp))))
      (declare (double-float damage-inc mass-total))
      (iterate-over-neighour-mps
       mesh mp length
       (lambda (mesh mp mp-other dist)
         (with-accessors ((d cl-mpm/particle::mp-damage)
                          (m cl-mpm/particle:mp-volume)
                          (ll cl-mpm/particle::mp-local-length)
                          (p cl-mpm/particle:mp-position))
             mp-other
           (when t;(< (the double-float d) 1d0)
             (push (cl-mpm/fastmaths::fast-.-
                    (cl-mpm/particle::mp-position mp-other)
                    (cl-mpm/particle::mp-position mp)
                    ) data-positions)
             (let ((weight
                     ;; (weight-func-mps mesh mp mp-other length)
                     (weight-func-mps-damaged mesh mp mp-other

                                              length
                                              )
                     ))
               (declare (double-float weight m d mass-total damage-inc))
               (incf mass-total (* weight m))
               (incf damage-inc
                     (* (the double-float (cl-mpm/particle::mp-damage-y-local mp-other))
                        weight m))
               (push (* weight m) data-weights)
               )
             ))
         ))
      (when (> mass-total 0d0)
        (setf damage-inc (/ damage-inc mass-total))
        (setf data-weights (mapcar (lambda (w) (/ w mass-total)) data-weights)))
      (values data-positions data-weights))))
