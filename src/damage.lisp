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
      (* (expt (max 0d0 (/ (- stress init-stress) init-stress)) 0.5d0) rate)
      0d0))

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
    (with-accessors ((stress cl-mpm/particle::mp-stress)
                     ;(stress cl-mpm/particle::mp-undamaged-stress)
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
                ;; (setf damage-increment (* s_1 (- 1d0 damage)))
                ;; (setf damage-increment s_1)
                (if (< damage 1d0)
                    (setf damage-increment (/ s_1 (expt (- 1d0 damage) 0.5d0)))
                  (setf damage-increment s_1)
                  )
          ;;       (let* ((omega (matrix-to-voigt (magicl:inv (magicl:.- (magicl:from-diag '(1d0 1d0)) (voigt-to-matrix damage-tensor)))))
          ;;              (identity (magicl:from-list '(1d0 1d0 0d0) '(3 1) :type 'double-float))
          ;;              (M (magicl:.+ (magicl:@ omega identity)
          ;;                            (magicl:transpose (magicl:@ identity omega))))
          ;;              (su (magicl:@ M stress))
          ;;              )
          ;; (multiple-value-bind (l v) (magicl:eig
          ;;                             su)
          ;;   (let* ()
          ;;     (loop for i from 0 to 1
          ;;           do (let* ((sii (nth i l)))
          ;;                (declare (double-float sii))
          ;;                (when (< sii 0d0)
          ;;                  (setf (nth i l)
          ;;                        0d0))
          ;;                ))
          ;;     (setf ybar-tensor (matrix-to-voight (magicl:@ v
          ;;                                       (magicl:from-diag l :type 'double-float)
          ;;                                       (magicl:transpose v)))))))
                ))
            (when (>= damage 1d0)
              (setf damage-increment 0d0))

            ;; (when (= damage 1d0)
            ;;   (setf damage-increment 0d0))
            ;;(setf ybar damage-increment)

            ;; (setf damage-increment (* dt (damage-rate-profile damage-increment damage damage-rate init-stress)))
            (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
            ;(setf local-length-t (length-localisation local-length local-length-damaged damage))
            ;; (setf local-length-t local-length)
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
          (setf damage-inc (* dt (damage-rate-profile damage-inc damage damage-rate init-stress)))
          (when (>= damage 1d0)
            (setf damage-inc 0d0)
            (setf ybar 0d0))
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
(defparameter *delocal-counter-max* 4)
(defun calculate-damage (mesh mps dt len non-local-damage)
  (when non-local-damage
    (when (<= *delocal-counter* 0)
      (create-delocalisation-list mesh mps len)
      (setf *delocal-counter* *delocal-counter-max*))
    (decf *delocal-counter*))
  (lparallel:pdotimes (i (length mps))
                      (when (typep (aref mps i) 'cl-mpm/particle:particle-damage)
                        (find-nodal-local-length mesh (aref mps i))
                        (calculate-damage-increment (aref mps i) dt)))
  (if non-local-damage
    (delocalise-damage mesh mps dt len)
    (localise-damage mesh mps dt))
  (lparallel:pdotimes (i (length mps))
                      (when (typep (aref mps i) 'cl-mpm/particle:particle-damage)
                        ;; (find-nodal-local-length mesh (aref mps i))
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

(defun local-list-add-particle (mesh mp)
  (let ((node-id (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp))))
    (when (cl-mpm/mesh:in-bounds mesh node-id)
      (let ((node (cl-mpm/mesh:get-node mesh node-id)))
        (sb-thread:with-mutex ((cl-mpm/mesh:node-lock node))
          (setf (cl-mpm/particle::mp-damage-position mp) (magicl:copy-tensor (cl-mpm/particle:mp-position mp)))
          (vector-push-extend mp (cl-mpm/mesh::node-local-list node)))))))

(defun local-list-remove-particle (mesh mp)
  (let ((node-id (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp))))
    (when (cl-mpm/mesh:in-bounds mesh node-id)
      (let ((node (cl-mpm/mesh:get-node mesh node-id)))
        (with-accessors ((ll cl-mpm/mesh::node-local-list)
                         (lock cl-mpm/mesh:node-lock)
                         ) node
          (sb-thread:with-mutex (lock)
            (setf ll (delete mp ll))))))))

(defun update-delocalisation-list (mesh mps length)
  (with-accessors ((nodes cl-mpm/mesh:mesh-nodes)
                   (h cl-mpm/mesh:mesh-resolution))
        mesh
      (lparallel:pdotimes (i (array-total-size nodes))
        (let ((node (row-major-aref nodes i)))
          (setf (fill-pointer (cl-mpm/mesh::node-local-list node)) 0)))
    (lparallel:pdotimes (i (length mps))
      (let ((mp (aref mps i)))
        (when (typep mp 'cl-mpm/particle:particle-damage)
          (let* ((delta (simd-accumulate (cl-mpm/particle:mp-position mp)
                                         (cl-mpm/particle::mp-damage-position mp)
                                         )))
            (when (> delta (/ h 4d0))
              (let ((node-id (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp))))
                (when (cl-mpm/mesh:in-bounds mesh node-id)
                  (let ((node (cl-mpm/mesh:get-node mesh node-id)))
                    (sb-thread:with-mutex ((cl-mpm/mesh:node-lock node))
                      (vector-push-extend mp (cl-mpm/mesh::node-local-list node)))))))))))
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

(defun diff-damaged (mesh mp-a mp-b)
  (flet ((node-dam (node)
           (if (cl-mpm/mesh::node-active node)
               (/ (cl-mpm/mesh::node-damage node)
                  (cl-mpm/mesh::node-svp-sum node)
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
        (magicl:.- pos-b pos-a diff)
        (let* ((length (sqrt (cl-mpm/fastmath::dot diff diff))))

          ;; Sample damage at midpoint and integrated with constant damage assumption
          ;; Can utilise linear shape function damage assumption but produced a nasty looking intergral
          (declare (double-float length final-distance))
          (when (> length 0d0)
            (let* (
                   (step-norm (magicl:scale diff (/ 1d0 length)))
                   ;;Start at point a and step through to b
                   (step-point (magicl:scale pos-a 1d0))
                   ;;Resolution of our midpoint integration
                   (step-size (/ h 1d0))
                   )
              (multiple-value-bind (steps remainder) (floor length step-size)
                (when (> steps 0)
                  (let* ((dhstep (magicl:scale step-norm (* 0.5d0 step-size)))
                         )
                    (loop for i from 0 below steps
                          do
                             (progn
                               (magicl:.+ step-point dhstep step-point)
                               (let ((damage 0d0))
                                 (cl-mpm::iterate-over-neighbours-point-linear-simd
                                  mesh
                                  step-point
                                  (lambda (m node weight grads)
                                    (incf damage
                                          (* (node-dam node) weight))))
                                 ;; (print damage)
                                 (incf final-distance (* step-size
                                                         (/ 1d0 (max epsilon (sqrt (- 1d0 damage)))))))
                               (magicl:.+ step-point dhstep step-point)
                               )
                          )))
                (let* ((step-size remainder)
                       (dhstep (magicl:scale step-norm (* 0.5d0 step-size)))
                       )
                  (progn
                    (magicl:.+ step-point dhstep step-point)
                    (let ((damage 0d0))
                      (cl-mpm::iterate-over-neighbours-point-linear-simd
                       mesh
                       step-point
                       (lambda (m node weight grads)
                         (incf damage
                               (* (node-dam node) weight))))
                                        ; (print damage)
                      (incf final-distance (* step-size (/ 1d0
                                                           (max epsilon
                                                                (the double-float (sqrt (- 1d0 damage))))))))
                    )
                  )
                )))
                                        ;(setf final-distance length)
          (* final-distance final-distance))
        ))))



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
 (ftype (function (cl-mpm/mesh::mesh
                   cl-mpm/particle:particle
                   cl-mpm/particle:particle
                   double-float
                   ) double-float)
        weight-func-mps
        ))
(defun weight-func-mps (mesh mp-a mp-b length)
  (weight-func (diff-squared mp-a mp-b) length))
(defun weight-func-mps-damaged (mesh mp-a mp-b length)
  (weight-func (diff-damaged mesh mp-a mp-b) length))


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
      (loop for dx fixnum from (- node-reach) to node-reach
            do (loop for dy fixnum from (- node-reach) to node-reach
                     do
                        (let ((idx (mapcar #'+ node-id (list dx dy))))
                          (declare (dynamic-extent idx))
                          (when (cl-mpm/mesh:in-bounds mesh idx)
                            (let ((node (cl-mpm/mesh:get-node mesh idx)))
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
                                                 (weight (weight-func-mps mesh mp mp-other (* 0.5d0 (+ length ll))))
                                                 ;;
                                                 ;; (weight (weight-func-mps-damaged mesh mp mp-other (* 0.5d0 (+ length ll))))
                                                 )
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
               )
                sim
    (declare (type double-float mass-filter))
                (progn
                    (cl-mpm::reset-grid mesh)
                    (cl-mpm::p2g mesh mps)
                    (when (> mass-filter 0d0)
                      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
                      ;; (cl-mpm::filter-grid-volume mesh 1d-8)
                    (cl-mpm::update-node-kinematics mesh dt )
                    (cl-mpm::apply-bcs mesh bcs dt)
                    (cl-mpm::update-stress mesh mps dt)
                    (when enable-damage
                     (cl-mpm/damage::calculate-damage mesh
                                                      mps
                                                      dt
                                                      50d0
                                                      nonlocal-damage
                                                      ))
                    ;Map forces onto nodes
                    (cl-mpm::p2g-force mesh mps)
                    ;(cl-mpm::apply-bcs mesh bcs-force dt)
                    (loop for bcs-f in bcs-force-list
                          do
                             (cl-mpm::apply-bcs mesh bcs-f dt))
                    (cl-mpm::update-node-forces mesh (cl-mpm::sim-damping-factor sim) dt (cl-mpm::sim-mass-scale sim))
                    ;Reapply velocity BCs
                    (cl-mpm::apply-bcs mesh bcs dt)
                    ;Also updates mps inline
                    (cl-mpm::g2p mesh mps dt)

                    (when remove-damage
                      (cl-mpm::remove-material-damaged sim))
                    (when split
                      (cl-mpm::split-mps sim))
                    (cl-mpm::check-mps mps)
                    )))


(defstruct damage-model
  name)

(defstruct (damage-model-creep (:include damage-model))
           (initiation-stress
            0d0
            :type DOUBLE-FLOAT)
           (damage-rate
            0d0
            :type DOUBLE-FLOAT))

(defstruct (damage-model-rankine (:include damage-model))
  (initiation-stress
   0d0
   :type DOUBLE-FLOAT))

(defgeneric damage-model-calculate-y (mp damage-model dt))

(defmethod damage-model-calculate-y (mp (dm damage-model) dt)
  0d0)
(defmethod damage-model-calculate-y (mp (dm damage-model-creep) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-stress)
                     (damage cl-mpm/particle:mp-damage)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (pressure cl-mpm/particle::mp-pressure)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     (damage-tensor cl-mpm/particle::mp-damage-tensor)
                     (ybar-tensor cl-mpm/particle::mp-damage-ybar-tensor)
                     (local-length cl-mpm/particle::mp-local-length)
                     (local-length-damaged cl-mpm/particle::mp-local-length-damaged)
                     (local-length-t cl-mpm/particle::mp-true-local-length)
                     ) mp
      (declare (double-float pressure damage))
        (progn
          (progn
            (let* ((j 1d0)
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
                   (vm (- s_1 s_2))
                   ;(s_1 vm)
                   )
              (when (> s_1 0d0)
                (if (< damage 1d0)
                    (setf damage-increment (/ s_1 (expt (- 1d0 damage) 0.5d0)))
                  (setf damage-increment s_1))))
            (when (>= damage 1d0)
              (setf damage-increment 0d0))
            ;; (setf local-length-t local-length)
            ;; (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
            damage-increment
            ))))
  )

(defmethod damage-model-calculate-y (mp (dm damage-model-rankine) dt) 
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-stress)
                     ;(stress cl-mpm/particle::mp-undamaged-stress)
                     (strain cl-mpm/particle::mp-strain)
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
          (multiple-value-bind (l v) (magicl:eig  (voight-to-matrix strain))
            (let* ((l (sort l #'>))
                   (s_1 (nth 0 l))
                   (s_2 (nth 1 l))
                   (pressure-effective (* pressure 1d0))
                   (s_1 (- s_1 pressure-effective))
                   (s_2 (- s_2 pressure-effective))
                   (s_1 (max 0d0 s_1))
                   (s_2 (max 0d0 s_2))
                   (vm (- s_1 s_2))
                   )
              (when (> s_1 0d0)
                (if (< damage 1d0)
                    (setf damage-increment (/ s_1 (expt (- 1d0 damage) 0.5d0)))
                  (setf damage-increment s_1))))
            (when (>= damage 1d0)
              (setf damage-increment 0d0))
            (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
            (setf local-length-t local-length)
            ))))
  )
(defgeneric damage-model-update-damage (mp damage-model dt))

(defmethod damage-model-update-damage (mp damage-model dt))

(defmethod damage-model-update-damage (mp (dm damage-model-creep) dt)
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
          (setf damage-inc (* dt (damage-rate-profile damage-inc damage damage-rate init-stress)))
          (when (>= damage 1d0)
            (setf damage-inc 0d0)
            (setf ybar 0d0))
          ;;Transform to log damage
          (incf damage damage-inc)
          ;;Transform to linear damage
          (setf damage (max 0d0 (min 1d0 damage)))
          (when (> damage critical-damage)
            (setf damage 1d0)
            (setf damage-inc 0d0)))
  (values)))
