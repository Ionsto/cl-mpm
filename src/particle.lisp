(defpackage :cl-mpm/particle
  (:use :cl
        :cl-mpm/utils)
  (:export
    #:make-particle
    #:make-particle-elastic
    #:make-particle-elastic-damage
    #:mp-mass
    #:mp-nd
    #:mp-volume
    #:mp-position
    #:mp-velocity
    #:mp-stress
    #:mp-strain
    #:mp-strain-rate
    #:mp-vorticity
    #:mp-gravity
    #:mp-body-force
    #:mp-damage
    #:mp-temperature
    #:mp-heat-capacity
    #:mp-critical-stress
    #:mp-deformation-gradient
    #:constitutive-model
    #:particle
    #:particle-damage
    #:post-stress-step
    )
  )
(in-package :cl-mpm/particle)
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(defstruct node-cache
  node
  weight
  grads
  weight-fbar
  grads-fbar)

(defclass particle ()
  ((mass
     :accessor mp-mass
     :type double-float
     :initarg :mass
     :initform 1d0)
   (nD
     :accessor mp-nd
     :type integer
     :initarg :nD)
   (index
     :accessor mp-index
     :type integer
     :initform -1
     :initarg :index)
   (mpi-index
    :accessor mp-mpi-index
    :type integer
    :initform -1
    :initarg :mpi-index)
   (volume
     :accessor mp-volume
     :type double-float
     :initarg :volume
     :initform 1d0)
   (volume-0
    :accessor mp-volume-0
    :type double-float
    :initarg :volume
    :initform 1d0)
   (size-0
    :accessor mp-domain-size-0
    :type magicl:matrix/double-float
    :initarg :size-0
    :initform (cl-mpm/utils:vector-zeros))
   (size
     :accessor mp-domain-size
     :type magicl:matrix/double-float
     :initarg :size
     :initform (cl-mpm/utils:vector-zeros))
   (position
     :accessor mp-position
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :initarg :position)
   (corners
    :accessor mp-corners
    :type (array MAGICL:MATRIX/DOUBLE-FLOAT 4))
   (velocity
     :accessor mp-velocity
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :initarg :velocity
     :initform (cl-mpm/utils:vector-zeros))
   (acceleration
    :accessor mp-acceleration
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils:vector-zeros))
   (stress
     :accessor mp-stress
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :initarg :stress
     :initform (cl-mpm/utils::voigt-zeros))
   (stress-kirchoff
    :accessor mp-stress-kirchoff
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initarg :stress
    :initform (cl-mpm/utils:voigt-zeros))
   (int-force
    :accessor mp-int-force
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils:vector-zeros))
   (strain
     :accessor mp-strain
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :initarg :strain
     :initform (cl-mpm/utils:voigt-zeros))
   (stretch-tensor
    :accessor mp-stretch-tensor
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils::matrix-zeros))
   (stretch-tensor-fbar
    :accessor mp-stretch-tensor-fbar
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils::matrix-zeros))
   (strain-rate-tensor
    :accessor mp-strain-rate-tensor
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils::matrix-zeros))
   (strain-rate
     :accessor mp-strain-rate
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :accessor mp-strain-rate
     :initarg :strain-rate
     :initform (cl-mpm/utils:voigt-zeros))
   (velocity-rate
    :accessor mp-velocity-rate
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :accessor mp-velocity-rate
    :initform (cl-mpm/utils::voigt-zeros))
   (eng-strain-rate
    :accessor mp-eng-strain-rate
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils:voigt-zeros))
   (vorticity
    :accessor mp-vorticity
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils::voigt-zeros))
   (deformation-gradient
     :accessor mp-deformation-gradient
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :initarg :deformation-gradient
     :initform (magicl:eye 3))
   (gravity
     :type double-float
     :accessor mp-gravity
     :initform 0d0;-9.8d0
     :initarg :gravity
     )
   (pressure
    :type double-float
    :accessor mp-pressure
    :initform 0d0
    )
   (pressure-datum
    :type double-float
    :accessor mp-pressure-datum
    :initform 0d0)
   (pressure-head
    :type double-float
    :accessor mp-pressure-head
    :initform 0d0)

   (boundary
    :type double-float
    :accessor mp-boundary
    :initform 0d0
    )
   (gravity-axis
    :type magicl:matrix/double-float
    :accessor mp-gravity-axis
    :initform (vector-from-list '(0d0 1d0 0d0))
    :initarg :gravity-axis)
   (body-force
     :accessor mp-body-force
     :type MAGICL:MATRIX/DOUBLE-FLOAT
     :initarg :body-force
     :initform (cl-mpm/utils::vector-zeros))
   (displacement
    :accessor mp-displacement
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils::vector-zeros))
   (cached-nodes
    :accessor mp-cached-nodes
    :initarg :nc
    :initform (make-array 8 :fill-pointer 0 :element-type 'node-cache))
   (p-modulus
    :accessor mp-p-modulus
    :initform 1d0
    )
   (damage
    :accessor mp-damage
    :type DOUBLE-FLOAT
    :initarg :damage
    :initform 0d0)
   (fixed-velocity
    :accessor mp-fixed-velocity
    :type list
    :initarg :fixed-velocity
    :initform nil)
   (j-fbar
    :accessor mp-j-fbar
    :initform 1d0)
   (split-depth
    :accessor mp-split-depth
    :type integer
    :initarg :split-depth
    :initform 0)
   (single-particle
    :accessor mp-single-particle
    :type boolean
    :initform nil)
   (fbar-j
    :accessor mp-debug-j
    :type double-float
    :initform 0d0)
   (fbar-j-gather
    :accessor mp-debug-j-gather
    :type double-float
    :initform 0d0)
   )
  (:documentation "A single material point"))

;; (defun mp-mass (mp)
;;   (sb-mop:standard-instance-access mp 0))
;; (defun mp-volume (mp)
;;   (sb-mop:standard-instance-access mp 4))
;; (defun mp-volume-0 (mp)
;;   (sb-mop:standard-instance-access mp 5))
;; (defun mp-deformation-gradient (mp)
;;   (sb-mop:standard-instance-access mp 22))

(defclass particle-elastic (particle)
  ((E
     :accessor mp-E
     :initarg :E
     )
   (nu
     :accessor mp-nu
     :initarg :nu)
   (elastic-matrix
    :accessor mp-elastic-matrix
    :type magicl:matrix/double-float)
   (2d-approximation
    :accessor mp-elastic-approximation
    :initarg :elastic-approxmation
    :initform :plane-strain
    :type t
    )
   )
  (:documentation "A linear-elastic material point"))

(defclass particle-vm (particle-elastic)
  ((rho
    :accessor mp-rho
    :initarg :rho
    )
   (strain-plastic-vm
    :accessor mp-strain-plastic-vm
    :type DOUBLE-FLOAT
    :initform 0d0)
   (strain-plastic
    :accessor mp-strain-plastic
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initarg :strain-plastic
    :initform (cl-mpm/utils:voigt-zeros))
   (yield-func
    :accessor mp-yield-func
    :type double-float
    :initform 0d0)
   )
  (:documentation "A vm perfectly plastic material point"))

(defclass particle-mc (particle-elastic)
  ((phi
    :accessor mp-phi
    :initarg :phi
    )
   (psi
    :accessor mp-psi
    :initarg :psi
    )
   (c
    :accessor mp-c
    :initarg :c
    )
   (strain-plastic-vm
    :accessor mp-strain-plastic-vm
    :type DOUBLE-FLOAT
    :initform 0d0)
   (strain-plastic
    :accessor mp-strain-plastic
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initarg :strain-plastic
    :initform (cl-mpm/utils:voigt-zeros))
   (yield-func
    :accessor mp-yield-func
    :type double-float
    :initform 0d0)
   )
  (:documentation "A vm perfectly plastic material point"))

(defun update-elastic-matrix (particle)
  (with-accessors ((de mp-elastic-matrix)
                   (E  mp-E)
                   (nu mp-nu)
                   (p mp-p-modulus)
                   )
      particle
    (setf p (/ E (* (+ 1d0 nu) (- 1d0 nu))))
    (setf de (cl-mpm/constitutive::linear-elastic-matrix E nu))
    ;; (if (eq (mp-elastic-approximation particle) :plane-strain)
    ;;     (setf de (cl-mpm/constitutive::linear-elastic-matrix E nu))
    ;;     (setf de (cl-mpm/constitutive::linear-elastic-matrix-ps E nu))
    ;;     )
    ))

(defmethod (setf mp-E) :after (value (p particle-elastic))
  (update-elastic-matrix p))
(defmethod (setf mp-nu) :after (value (p particle-elastic))
  (update-elastic-matrix p))
(defmethod initialize-instance :after ((p particle-elastic) &key)
  (update-elastic-matrix p))
(defmethod (setf mp-elastic-approximation) :after (value (p particle-elastic))
  (update-elastic-matrix p))

(defclass particle-fluid (particle)
  ((rest-density
    :accessor mp-rest-density
    :initarg :rest-density
    )
   (stiffness
    :accessor mp-stiffness
    :initarg :stiffness)
   (adiabatic-index
    :accessor mp-adiabatic-index
    :initarg :adiabatic-index)
   (viscosity
    :accessor mp-viscosity
    :initarg :viscosity)
   )
  (:documentation "A fluid material point"))

(defclass particle-viscoelastic (particle-elastic)
  (
   (viscosity
    :accessor mp-viscosity
    :initarg :viscosity)
   )
  (:documentation "A visco-elastic material point"))

(defclass particle-glen (particle-elastic)
  ((visc-factor
    :accessor mp-visc-factor
    :initarg :visc-factor)
   (visc-power
    :accessor mp-visc-power
    :initarg :visc-power)
   (true-visc
    :accessor mp-true-visc
    :initform 0d0)
   )
  (:documentation "A glen flow law material point"))

(defclass particle-viscoplastic (particle-elastic)
  ((visc-factor
    :accessor mp-visc-factor
    :initarg :visc-factor)
   (visc-power
    :accessor mp-visc-power
    :initarg :visc-power)
   (true-visc
    :accessor mp-true-visc
    :initform 0d0)
   (time-averaged-visc
    :accessor mp-time-averaged-visc
    :initform 0d0)
   )
  (:documentation "A visco-plastic material point"))

(defclass particle-m-ice (particle-viscoplastic)
  ((plastic-stress
    :accessor mp-plastic-stress
    :initarg :plastic-stress)
   (visc-plastic
    :accessor mp-visc-plastic
    :initform 0d0)
   (visc-glen
    :accessor mp-visc-glen
    :initform 0d0)
   )
  (:documentation "A visco-plastic material point"))

(defclass particle-damage (particle)
  (
   (log-damage
    :accessor mp-log-damage
    :type DOUBLE-FLOAT
    :initform 0d0)
   ;; (damage
   ;;  :accessor mp-damage
   ;;  :type DOUBLE-FLOAT
   ;;  :initarg :damage
   ;;  :initform 0d0)
   (local-damage
    :accessor mp-local-damage
    :type DOUBLE-FLOAT
    :initform 0d0)
   (local-damage-increment
    :accessor mp-local-damage-increment
    :type DOUBLE-FLOAT
    :initarg :damage-y
    :initform 0d0)
   (damage-increment
    :accessor mp-damage-increment
    :type DOUBLE-FLOAT
    :initform 0d0)
   (undamaged-stress
    :accessor mp-undamaged-stress
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils:voigt-zeros))
   (initiation-stress
    :accessor mp-initiation-stress
    :type DOUBLE-FLOAT
    :initarg :initiation-stress
    :initform 0d0)
   (critical-stress
    :accessor mp-critical-stress
    :type DOUBLE-FLOAT
    :initarg :critical-stress
    :initform 0d0)
   (damage-rate
    :accessor mp-damage-rate
    :type DOUBLE-FLOAT
    :initarg :damage-rate
    :initform 0d0)
   (critical-damage
    :accessor mp-critical-damage
    :type DOUBLE-FLOAT
    :initarg :critical-damage
    :initform 1d0)
   (damage-ybar
    :accessor mp-damage-ybar
    :type DOUBLE-FLOAT
    :initform 0d0
    :initarg :damage-ybar
    )
   (damage-y-local
    :accessor mp-damage-y-local
    :type DOUBLE-FLOAT
    :initform 0d0
    :initarg :damage-y
    )
   (time-averaged-ybar
    :accessor mp-time-averaged-ybar
    :initform 1d0)
   (time-averaged-damage-inc
    :accessor mp-time-averaged-damage-inc
    :initform 1d0)
   (time-averaged-counter
    :accessor mp-time-averaged-counter
    :initform 0d0)
   (local-length
    :accessor mp-local-length
    :type DOUBLE-FLOAT
    :initarg :local-length
    :initform 1d0)
   (local-length-damaged
    :accessor mp-local-length-damaged
    :type DOUBLE-FLOAT
    :initarg :local-length-damaged
    :initform 1d0)
   (true-local-length
    :accessor mp-true-local-length
    :type DOUBLE-FLOAT
    :initarg :local-length-t
    :initform 1d0)
   (damage-position
    :accessor mp-damage-position
    ;:type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform nil
    ;:initform (cl-mpm/utils::vector-zeros)
    )
   (damage-tensor
    :accessor mp-damage-tensor
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils::matrix-zeros)
    )
   (damage-ybar-tensor
    :accessor mp-damage-ybar-tensor
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils::matrix-zeros))
   (damage-model
    :accessor mp-damage-model
    :initarg :damage-model
    )
   (damage-domain-update-rate
    :accessor mp-damage-domain-update-rate
    :initarg :damage-domain-rate
    :initform 0d0))
  (:documentation "A material point with a damage tensor"))
(defclass particle-thermal (particle)
  (
   (temperature
    :accessor mp-temperature
    :type DOUBLE-FLOAT
    :initarg :temperature
    :initform 0d0)
   (heat-capacity
    :accessor mp-heat-capacity
    :type DOUBLE-FLOAT
    :initarg :heat-capacity
    :initform 0d0)
   (thermal-conductivity
    :accessor mp-thermal-conductivity
    :type DOUBLE-FLOAT
    :initarg :thermal-conductivity
    :initform 1d0)
   )
  (:documentation "A material point with a thermal properties"))
(defclass particle-fracture (particle-damage)
  (
   (strain-energy-density
    :accessor mp-strain-energy-density
    :type DOUBLE-FLOAT
    :initform 0d0)
   (strain-energy-density-local
    :accessor mp-strain-energy-density-local
    :type DOUBLE-FLOAT
    :initform 0d0)
   (fracture-toughness
    :accessor mp-fracture-toughness
    :type DOUBLE-FLOAT
    :initform 0d0
    :initarg :fracture-toughness)
   )
  (:documentation "A material point with fracture mechanics"))

(defclass particle-elastic-damage (particle-elastic particle-damage)
  ()
  (:documentation "A mp with damage influanced elastic model"))

(defclass particle-creep-damage (particle-elastic-damage)
  ()
  (:documentation "A mp with isotropic damage - reduced for use in creep test (higher performance)"))



(defclass particle-elastic-fracture (particle-elastic particle-fracture)
  ()
  (:documentation "A mp with fracture mechanics"))

(defclass particle-viscoelastic-damage (particle-viscoelastic particle-damage)
  ()
  (:documentation "A visco-elastic material point with damage"))

(defclass particle-viscoplastic-damage (particle-viscoplastic particle-damage)
  ()
  (:documentation "A mp with damage mechanics"))

(defclass particle-glen-damage (particle-glen particle-damage)
  ()
  (:documentation "A weakly compressible glen flow mp with damage mechanics"))

(defclass particle-thermoelastic-damage (particle-elastic particle-damage particle-thermal)
  ()
  (:documentation "A mp with elastic mechanics with variable thermal fields"))
(defclass particle-thermoviscoplastic-damage (particle-viscoplastic particle-damage particle-thermal)
  ()
  (:documentation "A mp with viscoplastic mechanics with variable thermal fields"))

(defclass particle-thermofluid-damage (particle-fluid particle-damage particle-thermal)
  ()
  (:documentation "A mp with viscoplastic mechanics with variable thermal fields"))
(defclass particle-viscoelastic-fracture (particle-viscoelastic particle-fracture)
  ()
  (:documentation "A viscoelastic mp with fracture mechanics"))

(defun make-particle (nD &optional (constructor 'particle) &rest args &key  (position nil) (volume 1) (mass 1) &allow-other-keys)
  (progn
    (if (eq position nil)
        (setf position (magicl:zeros (list nD 1)))
        (setf position (magicl:from-list (mapcar (lambda (x) (coerce x 'double-float)) position) (list nD 1))))
    (let ((stress-size 3))
      (let ((mp (apply #'make-instance constructor
                      :nD nD
                      :volume (coerce volume 'double-float)
                      :mass (coerce mass 'double-float)
                      :position position
                      args)))
        (progn
          (setf (mp-domain-size-0 mp) (magicl:scale (mp-domain-size-0 mp) 1d0))
          mp)))))
;; (defun make-particle (nD &rest args &key (constructor 'particle) (pos nil) (volume 1) (mass 1))
;;   (progn
;;     (if (eq pos nil)
;;         (setf pos (magicl:zeros (list nD 1)))
;;         (setf pos (magicl:from-list (mapcar (lambda (x) (coerce x 'double-float)) pos) (list nD 1))))
;;     (let ((stress-size 3))
;;       (make-instance constructor
;;                      :nD nD
;;                      :volume (coerce volume 'double-float)
;;                      :mass (coerce mass 'double-float)
;;                      :position pos))))
(defun make-particle-elastic (nD E nu &key (pos nil) (volume 1) (mass 1))
    (let ((p (make-particle nD 'particle-elastic :position pos :volume volume :mass mass) ))
        (progn
            (setf (mp-E p) E)
            (setf (mp-nu p) nu)
        p)))

(defun make-particle-elastic-damage (nD E nu &key (pos nil) (volume 1) (mass 1))
    (let ((p (make-particle nD 'particle-elastic-damage :position pos :volume volume :mass mass
                            ) ))
        (progn
            (setf (mp-E p) E)
            (setf (mp-nu p) nu)
        p)))

(defgeneric constitutive-model (mp elastic-trial-strain dt)
    (:documentation "Compute new stress state given elastic strain")
    (:method (mp strain dt)
        (magicl:scale strain 0)))

(defmethod constitutive-model :before ((mp particle-damage) strain dt)
  "For rate-based equations we want them to see the undamaged stress"
  (with-slots ((stress stress)
               (stress-u undamaged-stress))
      mp
    ;; (setf stress (magicl:scale stress-u 1d0))
    ))


(defmethod constitutive-model ((mp particle-elastic) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((de elastic-matrix))
      mp
    (cl-mpm/constitutive::linear-elastic-mat strain de)))

(defmethod constitutive-model ((mp particle-elastic-damage) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (stress stress)
               (stress-undamaged undamaged-stress)
               (strain-rate strain-rate)
               (D stretch-tensor)
               (velocity-rate velocity-rate)
               (vorticity vorticity)
               (def deformation-gradient)
               (damage damage)
               (pressure pressure)
               ;; (datum pressure-datum)
               ;; (rho pressure-head)
               (pos position)
               (calc-pressure pressure-func)
               )
      mp
    (declare (double-float pressure damage)
             (function calc-pressure))
    ;; Non-objective stress intergration

    ;; (setf stress-undamaged (cl-mpm/constitutive::linear-elastic-mat strain de))
    ;; (setf pressure (funcall calc-pressure pos))
    (magicl.simd::.+-simd
     stress-undamaged
     (objectify-stress-logspin
      (cl-mpm/constitutive::linear-elastic-mat strain-rate de)
      stress-undamaged
      def
      vorticity
      D)
     stress-undamaged)

    (setf stress (magicl:scale stress-undamaged 1d0))

    (when (> damage 0.0d0)
      (let* ((j 1d0)
             (p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
             (s (cl-mpm/constitutive::deviatoric-voigt stress)))
        (setf stress (magicl:.+ (cl-mpm/constitutive::voight-eye p)
                                (magicl:scale! s (max 1d-3 (expt (- 1d0 damage) 2d0)))
                                ))
        (multiple-value-bind (l v) (cl-mpm/utils::eig
                                    (magicl:scale! (voight-to-matrix stress) (/ 1d0 j)))
          (let* (;(tp (funcall calc-pressure (magicl:tref pos 1 0) datum rho))
                 (tp (funcall calc-pressure pos))
                 (driving-pressure (* tp 1d0 (expt (min 1.00d0 damage) 1)))
                 (degredation (expt (- 1d0 damage) 2d0)))
            (setf pressure tp)
            (loop for i from 0 to 2
                  do
                     (let* ((sii (nth i l))

                              (esii (- sii driving-pressure)))
                         (when (> esii 0d0)
                           ;;Tensile damage -> unbounded
                           (setf (nth i l) (* esii (max 1d-6 degredation)))
                           (setf (nth i l) (+ (nth i l) driving-pressure))
                           )
                         (when (< esii 1d0)
                           ;;Bounded compressive damage
                           (setf (nth i l) (* esii (max 1d0 degredation)))
                           (setf (nth i l) (+ (nth i l) driving-pressure))
                           )
                         ;; (setf (nth i l) (* sii (max 0d0 (- 1d0 damage))))
                         )
                  )
              (setf stress (magicl:scale! (matrix-to-voight (magicl:@ v
                                                                      (magicl:from-diag l :type 'double-float)
                                                                      (magicl:transpose v))) j))
              ))
        ))
    stress
    ))

(defmethod constitutive-model ((mp particle-creep-damage) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (stress stress)
               (stress-undamaged undamaged-stress)
               (strain-rate strain-rate)
               (def deformation-gradient)
               (damage damage)
               )
      mp
    (declare (double-float damage))
    (setf stress-undamaged (cl-mpm/constitutive::linear-elastic-mat strain de))
    (setf stress (magicl:scale stress-undamaged (- 1d0 damage)))
    stress
    ))

(defclass particle-elastic-inc (particle-elastic)
  ()
  (:documentation "A linear-elastic material point"))
(defmethod constitutive-model ((mp particle-elastic-inc) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (stress stress)
               (strain-rate strain-rate)
               (velocity-rate velocity-rate)
               (vorticity vorticity)
               (def deformation-gradient)
               )
      mp
    ;; Non-objective stress intergration
    ;; (magicl.simd::.+-simd
    ;;  stress
    ;;  (cl-mpm/constitutive::linear-elastic-mat strain-rate de))
    (magicl:.+ stress (magicl:@ de strain-rate))
    ))
(defclass particle-elastic-jaumann (particle-elastic)
  ()
  (:documentation "A linear-elastic material point"))
(defmethod constitutive-model ((mp particle-elastic-jaumann) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (stress stress)
               (strain-rate strain-rate)
               (velocity-rate velocity-rate)
               (vorticity vorticity)
               (def deformation-gradient)
               )
      mp
    ;;Jaumann rate equation
    (magicl.simd::.+-simd
     stress
     (objectify-stress-jaumann
      (cl-mpm/constitutive::linear-elastic-mat strain-rate de)
      stress
      vorticity))
    ))
(defclass particle-elastic-truesdale (particle-elastic)
  ()
  (:documentation "A linear-elastic material point"))
(defmethod constitutive-model ((mp particle-elastic-truesdale) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (stress stress)
               (strain-rate strain-rate)
               (velocity-rate velocity-rate)
               (vorticity vorticity)
               (def deformation-gradient)
               )
      mp
    ;; Truesdale rate
    (magicl.simd::.+-simd
     stress
     (objectify-stress-kirchoff-truesdale
      (cl-mpm/constitutive::linear-elastic-mat strain-rate de)
      stress
      strain-rate))
    ))
(defclass particle-elastic-logspin (particle-elastic)
  ()
  (:documentation "A linear-elastic material point"))
(defmethod constitutive-model ((mp particle-elastic-logspin) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (stress stress-kirchoff)
               (strain-rate strain-rate)
               (D stretch-tensor)
               (vorticity vorticity)
               (def deformation-gradient))
      mp
    (let (
          ;; (strain-rate (magicl:scale strain-rate (/ 1d0 dt)))
          ;; (vorticity (magicl:scale vorticity (/ 1d0 dt)))
          )
      (magicl.simd::.+-simd
       stress
       (objectify-stress-logspin
        (cl-mpm/constitutive::linear-elastic-mat strain-rate de)
        stress
        def
        vorticity
        ;; strain-rate
        D
        )))
    ))

(defmethod constitutive-model ((mp particle-fluid) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((viscosity viscosity)
               (mass mass)
               (volume volume)
               (rest-density rest-density)
               (stiffness stiffness)
               (adiabatic-index adiabatic-index)
               )
      mp
    (let* ((density (/ mass volume))
           (pressure (* stiffness (expt (- (/ density rest-density) 1) adiabatic-index))))
      (cl-mpm/constitutive:newtonian-fluid strain pressure viscosity))))

(defmethod constitutive-model ((mp particle-viscoelastic) strain dt)
  "Function for modelling stress intergrated viscoelastic maxwell material"
  (with-slots ((E E)
               (nu nu)
               (viscosity viscosity)
               (strain-rate strain-rate) ;Note strain rate is actually strain increment through dt
               (vorticity vorticity)
               (stress stress))
      mp
    (cl-mpm/constitutive:maxwell strain-rate stress E nu viscosity dt vorticity)))

(defmethod constitutive-model ((mp particle-glen) strain dt)
  "Function for modeling stress intergrated viscoplastic norton-hoff material"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (visc-factor visc-factor)
               (visc-power visc-power)
               (strain-rate strain-rate) ;Note strain rate is actually strain increment through dt
               (velocity-rate velocity-rate) ;Note strain rate is actually strain increment through dt
               (stress stress)
               (strain strain)
               (stretch stretch-tensor)
               (eng-strain-rate eng-strain-rate)
               )
      mp
    (declare (double-float E visc-factor visc-power))
    (let* (;(viscosity (cl-mpm/constitutive::glen-viscosity-strain velocity-rate visc-factor visc-power))
           (viscosity (cl-mpm/constitutive::glen-viscosity-strain eng-strain-rate visc-factor visc-power)))
      (cl-mpm/constitutive::elasto-glen strain-rate stress E nu de viscosity dt strain)
      )))
(defmethod constitutive-model ((mp particle-viscoplastic) strain dt)
  "Function for modeling stress intergrated viscoplastic norton-hoff material"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (visc-factor visc-factor)
               (visc-power visc-power)
               (strain-rate strain-rate) ;Note strain rate is actually strain increment through dt
               (velocity-rate velocity-rate) ;Note strain rate is actually strain increment through dt
               (strain-plastic strain-plastic)
               (def deformation-gradient)
               (vorticity vorticity)
               (D stretch-tensor)
               (stress stress)
               (temp true-visc)
               (eng-strain-rate eng-strain-rate)
               )
      mp
    (declare (double-float E visc-factor visc-power))
    (let* ((eng-strain-rate (magicl.simd::.*-simd (magicl:map (lambda (x) (* x (exp x))) strain) velocity-rate
                                                  (cl-mpm/utils:voigt-from-list '(1d0 1d0 1d0 0.5d0 0.5d0 0.5d0))))
          (viscosity (cl-mpm/constitutive::glen-viscosity-strain eng-strain-rate visc-factor visc-power))
          ;(viscosity (cl-mpm/constitutive::glen-viscosity-stress stress visc-factor visc-power))
          )
      ;; stress
      (setf temp viscosity)
      (magicl.simd::.+-simd
       stress
       (objectify-stress-logspin
        (if (> viscosity 0d0)
            (cl-mpm/constitutive::maxwell strain-rate stress E nu de viscosity dt)
            (cl-mpm/constitutive::linear-elastic-mat strain-rate de))
        stress
        def
        vorticity
        D
        ))
      )
    ))

(defmethod constitutive-model ((mp particle-m-ice) strain dt)
  "Function for modeling stress intergrated viscoplastic norton-hoff material"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (visc-factor visc-factor)
               (visc-power visc-power)
               (strain-rate strain-rate) ;Note strain rate is actually strain increment through dt
               (velocity-rate velocity-rate) ;Note strain rate is actually strain increment through dt
               (strain-plastic strain-plastic)
               (def deformation-gradient)
               (vorticity vorticity)
               (D stretch-tensor)
               (stress stress)
               (temp true-visc)
               (plastic-stress plastic-stress)
               (eng-strain-rate eng-strain-rate)
               (visc-plastic visc-plastic)
               (visc-glen visc-glen)
               (time-averaged-visc time-averaged-visc)
               (p p-modulus)
               )
      mp
    (declare (double-float E visc-factor visc-power))
    (flet ((inv (x) (/ 1d0 (the double-float x))))
      (let* (;(eng-strain-rate (magicl.simd::.*-simd (magicl:map (lambda (x) (* x (exp x))) strain) velocity-rate
                                        ;                            (cl-mpm/utils:stress-from-list '(1d0 1d0 0.5d0))))
             (estrain (cl-mpm/constitutive::effective-strain-rate eng-strain-rate))
             (viscosity-diff (/ 6.5d8 (* 2d0 (expt 10d-3 2))))
             (viscosity-glen (cl-mpm/constitutive::glen-viscosity-strain eng-strain-rate visc-factor visc-power))
             (viscosity-plastic (/ plastic-stress (+ 1d-40 (* 2d0 estrain))))
             (viscosity (+ 1d-40 (inv (+ ;(inv viscosity-diff)
                                         (inv viscosity-glen)
                                         (inv viscosity-plastic)))))
             ;; (viscosity viscosity-glen)
             )
        ;; stress
        ;; (print viscosity-glen)
        ;; ;; (print viscosity-plastic)
        ;; (setf p
        ;;       (* (/ E (* (+ 1 nu) (- 1 nu))) (/ 1d12 viscosity)))
        (setf temp viscosity)
        (setf visc-plastic viscosity-plastic
              visc-glen viscosity-glen)
        (incf time-averaged-visc viscosity)
        ;; (setf stress-u
        ;;       (cl-mpm/constitutive:maxwell-exp-v strain-rate stress E nu de visc-u dt))
        ;; (setf stress
        ;;       (cl-mpm/constitutive:maxwell-exp-v strain-rate stress E nu de viscosity dt))

        (magicl.simd::.+-simd
         stress
         (objectify-stress-logspin
          (if (> viscosity 0d0)
              (cl-mpm/constitutive::maxwell-exp-inc strain-rate stress E nu de viscosity dt)
              (cl-mpm/constitutive::linear-elastic-mat strain-rate de))
          stress
          def
          vorticity
          D
          ))
        ))
    ))

(defun inplace-stress-update (mp strain dt stress)
  (declare (cl-mpm/particle::particle-viscoplastic mp)
           (magicl:matrix/double-float strain stress)
           (double-float dt))
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (visc-factor visc-factor)
               (visc-power visc-power)
               (strain-rate strain-rate) ;Note strain rate is actually strain increment through dt
               (velocity-rate velocity-rate) ;Note strain rate is actually strain increment through dt
               (strain-plastic strain-plastic)
               (def deformation-gradient)
               (vorticity vorticity)
               (D stretch-tensor)
               )
      mp
    (declare (double-float E visc-factor visc-power))
    (let ((viscosity (cl-mpm/constitutive::glen-viscosity-strain velocity-rate visc-factor visc-power)))
      (magicl.simd::.+-simd
       stress
       (objectify-stress-logspin
        ;; (if (> viscosity 0d0)
            (cl-mpm/constitutive::maxwell strain-rate stress E nu de viscosity dt)
            ;; (cl-mpm/constitutive::linear-elastic-mat strain-rate de))
        stress
        def
        vorticity
        D
        )
       stress))))

(defmethod constitutive-model ((mp particle-viscoplastic-damage) strain dt)
  "Function for modeling stress intergrated viscoplastic norton-hoff material"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (visc-factor visc-factor)
               (visc-power visc-power)
               (strain-rate strain-rate) ;Note strain rate is actually strain increment through dt
               (velocity-rate velocity-rate) ;Note strain rate is actually strain increment through dt
               (strain-plastic strain-plastic)
               (damage damage)
               (critical-damage critical-damage)
               (def deformation-gradient)
               (D stretch-tensor)
               (vorticity vorticity)
               (stress stress)
               (stress-u undamaged-stress)
               (pressure pressure)
               (pos position)
               (eng-strain-rate eng-strain-rate)
               (time-averaged-visc time-averaged-visc)
               (p p-modulus)

               (calc-pressure pressure-func)

               ;; (stress undamaged-stress)
               ;; (stress-damaged stress)
               )
      mp
    (declare (double-float E visc-factor visc-power))
    (let* (;; (eng-strain-rate (magicl.simd::.*-simd (magicl:map (lambda (x) (* x (exp x))) strain) velocity-rate
           ;;                                        (cl-mpm/utils:voigt-from-list '(1d0 1d0 1d0 0.5d0 0.5d0 0.5d0))))
           (viscosity (cl-mpm/constitutive::glen-viscosity-strain eng-strain-rate visc-factor visc-power))
           ;; (viscosity 1d-20)
           (visc-u viscosity)
           (viscosity (* viscosity (max 1d-10 (expt (- 1d0 damage) 2))))
           )
      ;; (setf true-visc viscosity)
      ;; (setf p
      ;;       (/ (/ E (* (+ 1 nu) (- 1 nu)))
      ;;          (expt (max 1d-2 (- 1d0 damage)) 2)
      ;;          ;(expt (/ 1d13 viscosity) 1)
      ;;          ))
      ;; (setf stress (magicl:scale stress-undamaged (- 1d0 damage)))

      (incf time-averaged-visc viscosity)
      ;; (setf stress-u
      ;;       (magicl.simd::.+-simd
      ;;        stress-u
      ;;        (objectify-stress-logspin
      ;;         (if (> viscosity 0d0)
      ;;             ;; (cl-mpm/constitutive::maxwell-damage strain-rate stress E nu de viscosity dt damage)
      ;;             ;; (cl-mpm/constitutive:maxwell-exp-v strain-rate stress E nu de visc-u dt))
      ;;             (cl-mpm/constitutive::maxwell strain-rate stress E nu de viscosity dt)
      ;;             (cl-mpm/constitutive::linear-elastic-mat strain-rate de))
      ;;         stress-u
      ;;         def
      ;;         vorticity
      ;;         D)))

      ;; (magicl:.+
      ;;  stress-u
      ;;  (cl-mpm/constitutive:maxwell strain-rate stress E nu de viscosity dt)
      ;;  stress-u)
      (setf stress-u

            (cl-mpm/constitutive:maxwell-exp-v
             strain-rate
             stress-u
             E nu de
             viscosity dt))
      ;; (when (> (abs(magicl:tref stress-u 3 0)) 1d-10)
      ;;   (break)
      ;;   )
      ;; (setf (magicl:tref stress-u 3 0) 0d0
      ;;       (magicl:tref stress-u 4 0) 0d0
      ;;       )


      ;; (setf stress-u (cl-mpm/constitutive:maxwell strain-rate stress E nu de viscosity dt))
      ;; (setf stress
      ;;       (cl-mpm/constitutive:maxwell-exp-v strain-rate stress E nu de viscosity dt))
      ;; (setf stress-u (cl-mpm/constitutive::linear-elastic-mat strain de))
      (setf stress (magicl:scale stress-u 1d0))
      (when (> damage 0.0d0)
        (let ((j (magicl:det def)))
          (multiple-value-bind (l v) (cl-mpm/utils::eig
                                      (magicl:scale! (voight-to-matrix stress) (/ 1d0 j)))

            (let* ((tp (funcall calc-pressure pos))
                  (driving-pressure (* tp (expt (min 0.90d0 damage) (magicl:det def))))
                  ;; (driving-pressure 0d0) ;;(* pressure (min 1.00d0 damage)))
                  (degredation (expt (- 1d0 damage) 2d0)))

              (loop for i from 0 to 2
                    do (let* ((sii (nth i l))
                              (esii (- sii driving-pressure)))
                         (when (> esii 0d0)
                           ;;Tensile damage -> unbounded
                           (setf (nth i l) (* esii (max 1d-8 degredation)))
                           (setf (nth i l) (+ (nth i l) driving-pressure))
                           )
                         (when (< esii 0d0)
                           ;;Bounded compressive damage
                           (setf (nth i l) (* esii (max 1d-1 degredation)))
                           (setf (nth i l) (+ (nth i l) driving-pressure))
                           )
                         ;; (setf (nth i l) (* sii (max 0d0 (- 1d0 damage))))
                         ))
              (setf stress (magicl:scale! (matrix-to-voight (magicl:@ v
                                                                      (magicl:from-diag l :type 'double-float)
                                                                      (magicl:transpose v))) j))
              ))))
      ;; (when nil;(> damage 0.0d0)
      ;;   (let ((j 1d0))
      ;;     (multiple-value-bind (l v) (cl-mpm/utils::eig
      ;;                                 (magicl:scale! (voight-to-matrix stress) (/ 1d0 j)))
      ;;       (let ((driving-pressure (* pressure (min 1.00d0 damage)))
      ;;             (degredation (expt (- 1d0 damage) 2d0)))
      ;;         (loop for i from 0 to 2
      ;;               do (let* ((sii (nth i l))
      ;;                         (esii (- sii driving-pressure)))
      ;;                    (when (> esii 0d0)
      ;;                      ;;Tensile damage -> unbounded
      ;;                      (setf (nth i l) (* esii (max 0d-8 degredation)))
      ;;                      (setf (nth i l) (+ (nth i l) driving-pressure))
      ;;                      )
      ;;                    (when (< esii 0d0)
      ;;                      ;;Bounded compressive damage
      ;;                      (setf (nth i l) (* esii (max 1d-2 degredation)))
      ;;                      (setf (nth i l) (+ (nth i l) driving-pressure))
      ;;                      )
      ;;                    ;; (setf (nth i l) (* sii (max 0d0 (- 1d0 damage))))
      ;;                    ))
      ;;         (setf stress (magicl:scale! (matrix-to-voight (magicl:@ v
      ;;                                                                 (magicl:from-diag l :type 'double-float)
      ;;                                                                 (magicl:transpose v))) j))
      ;;         ))))

      stress
      ;; (cl-mpm/constitutive::elasto-glen strain-rate stress E nu de viscosity dt)
      )))

(defmethod constitutive-model ((mp particle-glen-damage) strain dt)
  "Function for modeling stress intergrated viscoplastic norton-hoff material"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (visc-factor visc-factor)
               (visc-power visc-power)
               (strain-rate strain-rate) ;Note strain rate is actually strain increment through dt
               (velocity-rate velocity-rate) ;Note strain rate is actually strain increment through dt
               (stress stress)
               (strain strain)
               (damage damage)
               )
      mp
    (declare (double-float E visc-factor visc-power))
    (let* ((viscosity (cl-mpm/constitutive::glen-viscosity-strain velocity-rate visc-factor visc-power))
          ;(viscosity (* viscosity (- 1 (* damage (- 1 0.1)))))
           (viscosity (* viscosity (max 1d-3 (expt (- 1d0 damage) visc-power))))
          )
      ;; (cl-mpm/constitutive::elasto-glen-damage strain-rate stress E nu de viscosity dt strain damage)
      (cl-mpm/constitutive::elasto-glen strain-rate stress E nu de viscosity dt strain)
      )
    ))

(defmethod constitutive-model ((mp particle-viscoelastic-damage) strain dt)
  "Function for modelling stress intergrated viscoelastic maxwell material"
  (with-slots ((E E)
               (nu nu)
               (viscosity viscosity)
               (strain-rate strain-rate) ;Note strain rate is actually strain increment through dt
               (vorticity vorticity)
               (def deformation-gradient)
               (stress undamaged-stress)
               ;; (stress stress)
               )
      mp
    ;; (magicl.simd::.+-simd stress (cl-mpm/constitutive:linear-elastic strain-rate E nu) (objectify-stress-jaumann stress vorticity))
   ;(cl-mpm/constitutive:maxwell strain-rate (magicl:scale stress (magicl:det def)) E nu viscosity dt vorticity)
    (cl-mpm/constitutive:maxwell-exp strain-rate stress E nu viscosity dt vorticity)
    ))

(defmethod constitutive-model ((mp particle-thermoviscoplastic-damage) strain dt)
  "Function for modelling stress intergrated viscoplastic norton-hoff material"
  (with-slots ((E E)
               (nu nu)
               (visc-factor visc-factor)
               (visc-power visc-power)
               (strain-rate strain-rate) ;Note strain rate is actually strain increment through dt
               (strain-plastic strain-plastic)
               (deformation-gradient deformation-gradient)
               (vorticity vorticity)
               (temperature temperature)
               (stress-u undamaged-stress))
      mp
    (let ((viscosity (cl-mpm/constitutive::glen-viscosity stress-u visc-factor visc-power)))
            (cl-mpm/constitutive::maxwell strain-rate stress-u E nu viscosity dt vorticity))
    ))

(defmethod constitutive-model ((mp particle-thermofluid-damage) strain dt)
  "Function for modelling stress intergrated viscoplastic norton-hoff material"
  (with-slots ((viscosity viscosity)
               (mass mass)
               (volume volume)
               (rest-density rest-density)
               (stiffness stiffness)
               (temp temperature)
               (adiabatic-index adiabatic-index)
               )
      mp
    (let* ((density (/ mass volume))
           (effective-rest-density rest-density)
           (pressure (* stiffness
                        (expt
                         (- (/ density effective-rest-density)
                            1) adiabatic-index))))
      (cl-mpm/constitutive:newtonian-fluid strain pressure viscosity)))
  )

(defgeneric post-stress-step (mesh mp dt)
  (:documentation "This step gets called after full stress state resolved and allows for other processing"))
(defmethod post-stress-step (mesh mp dt)
  ())

(declaim (inline assemble-vorticity-matrix)
         (ftype (function (magicl:matrix/double-float)
                          magicl:matrix/double-float
                          ) assemble-vorticity-matrix))
(defun assemble-vorticity-matrix (vorticity)
  (let ((dx (magicl:tref vorticity 0 0))
        (dy (magicl:tref vorticity 1 0))
        (dxdy (magicl:tref vorticity 2 0))
        )
    (declare (double-float dx dy dxdy))
    (cl-mpm/utils::matrix-from-list (list dx dxdy (- dxdy) dy))))

(defun objectify-stress (mp)
  (cl-mpm/particle:mp-stress mp))

(defun objectify-stress-jaumann (stress-inc stress vorticity)
  (magicl:.-
   stress-inc
   (matrix-to-voight
    (magicl::.- (magicl:@ (voight-to-matrix stress) (assemble-vorticity-matrix vorticity))
                (magicl:@ (assemble-vorticity-matrix vorticity) (voight-to-matrix stress)))))
             )

(defun objectify-stress-kirchoff-truesdale (stress-inc stress velocity-rate)
  (magicl:.-
   stress-inc
   (matrix-to-voight
    (magicl::.- (magicl:@ (voight-to-matrix stress) (voight-to-matrix velocity-rate))
                (magicl:@ (voight-to-matrix velocity-rate) (voight-to-matrix stress))
                ))))

(declaim (inline objectify-stress-logspin)
         (ftype (function (magicl:matrix/double-float
                           magicl:matrix/double-float
                           magicl:matrix/double-float
                           magicl:matrix/double-float
                           magicl:matrix/double-float
                           )
                          magicl:matrix/double-float
                          )
                objectify-stress-logspin))
(defun objectify-stress-logspin (stress-inc stress def vorticity D)
    (let ((b (magicl:@ def (magicl:transpose def)))
          ;; (omega (assemble-vorticity-matrix vorticity))
          (omega (magicl:scale! (magicl:.- D (magicl:transpose D)) 0.5d0))
          (D (magicl:scale! (magicl.simd::.+-simd D (magicl:transpose D)) 0.5d0))
          ;; (D (cl-mpm/utils::voigt-to-matrix D))
          )
        (multiple-value-bind (l v) (cl-mpm/utils::eig b)
          (loop for i from 0 to 2
                do (loop for j from 0 to 2
                         do
                            ;;For all the pairs of eigenvalues
                            (when (not (= i j))
                              (let ((l_i (nth i l))
                                    (l_j (nth j l))
                                    (v_i (magicl:column v i))
                                    (v_j (magicl:column v j)))
                                (declare (double-float l_i l_j)
                                         (magicl:matrix/double-float v_i v_j))
                                ;;When the eigenvalues are distinct
                                (when (and
                                       ;;When they are nonzero
                                       (> (abs (- l_i l_j)) 1d-6)
                                       )
                                  ;; When we have pairs of unique nonzero eigenvalues
                                  (setf omega
                                        (magicl.simd::.+-simd omega
                                                   (magicl:scale!
                                                    (magicl:@
                                                     (magicl:@
                                                      v_i
                                                      (magicl:transpose v_i))
                                                     D
                                                     (magicl:@
                                                      v_j
                                                      (magicl:transpose v_j)))
                                                    (+
                                                     (/ (+ 1d0 (/ l_i l_j)) (- 1d0 (/ l_i l_j)))
                                                     (/ 2d0 (the double-float (log (/ l_i l_j)))))
                                                    )))
                                  ))))))
      (magicl:.-
       stress-inc
       (cl-mpm/utils::matrix-to-voight
        (magicl::.- (magicl:@ (cl-mpm/utils::voight-to-matrix stress) omega)
                    (magicl:@ omega (cl-mpm/utils::voight-to-matrix stress)))))
      ))


(defmethod constitutive-model ((mp particle-vm) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-accessors ((de mp-elastic-matrix)
                   (stress mp-stress)
                   (rho mp-rho)
                   (plastic-strain mp-strain-plastic)
                   (ps-vm mp-strain-plastic-vm)
                   (strain mp-strain)
                   (yield-func mp-yield-func)
                   )
      mp
    ;;Train elastic strain - plus trail kirchoff stress
    (setf stress
          (cl-mpm/constitutive::linear-elastic-mat strain de))
    (multiple-value-bind (sig eps-e f) (cl-mpm/constitutive::vm-plastic stress de strain rho)
      (setf stress
            sig
            plastic-strain (magicl:.- strain eps-e)
            strain eps-e
            yield-func f
            ))
    (incf ps-vm
          (multiple-value-bind (l v)
                     (cl-mpm/utils:eig (cl-mpm/utils:voigt-to-matrix (cl-mpm/particle::mp-strain-plastic mp)))
                   (destructuring-bind (s1 s2 s3) l
                     (sqrt
                      (/ (+ (expt (- s1 s2) 2d0)
                            (expt (- s2 s3) 2d0)
                            (expt (- s3 s1) 2d0)
                            ) 2d0)))))
    stress
    ))

(defmethod constitutive-model ((mp particle-mc) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-accessors ((de mp-elastic-matrix)
                   (stress mp-stress)
                   (E mp-E)
                   (nu mp-nu)
                   (phi mp-phi)
                   (psi mp-psi)
                   (c mp-c)
                   (plastic-strain mp-strain-plastic)
                   (ps-vm mp-strain-plastic-vm)
                   (strain mp-strain)
                   (yield-func mp-yield-func)
                   )
      mp
    ;;Train elastic strain - plus trail kirchoff stress
    (setf stress
          (cl-mpm/constitutive::linear-elastic-mat strain de))
    (multiple-value-bind (sig eps-e f) (cl-mpm/constitutive::mc-plastic stress de strain E nu phi psi c)
      (setf stress
            sig
            plastic-strain (magicl:.- strain eps-e)
            strain eps-e
            yield-func f
            ))
    (incf ps-vm
          (multiple-value-bind (l v)
                     (cl-mpm/utils:eig (cl-mpm/utils:voigt-to-matrix (cl-mpm/particle::mp-strain-plastic mp)))
                   (destructuring-bind (s1 s2 s3) l
                     (sqrt
                      (/ (+ (expt (- s1 s2) 2d0)
                            (expt (- s2 s3) 2d0)
                            (expt (- s3 s1) 2d0)
                            ) 2d0)))))
    stress
    ))


(defclass particle-chalk (particle-elastic-damage)
  (
   (coheasion
    :accessor mp-coheasion
    :initarg :coheasion
    :initform 0d0)
   (friction-angle
    :accessor mp-friction-angle
    :initarg :friction-angle
    :initform 30d0)
   (max-strain
    :initform 0d0
    )
   (ductility
    :accessor mp-ductility
    :initarg :ductility
    :initform 1d0)

   (shear-residual-ratio
    :accessor mp-shear-residual-ratio
    :initarg :g-res-ratio
    :initform 1d-3
    )
   (k-tensile-residual-ratio
    :accessor mp-k-tensile-residual-ratio
    :initarg :kt-res-ratio
    :initform 1d-6
    )
   (k-compressive-residual-ratio
    :accessor mp-k-compressive-residual-ratio
    :initarg :kc-res-ratio
    :initform 1d-2
    )
   )
  (:documentation "A chalk damage model"))

(defclass particle-chalk-creep (particle-chalk)
  ()
  (:documentation "A chalk damage model"))

(defclass particle-chalk-brittle (particle-chalk particle-mc)
  (
   (enable-plasticity
    :accessor mp-enable-plasticity
    :initarg :enable-plasticity
    :initform t
    )
   (fracture-energy
    :accessor mp-gf
    :initarg :fracture-energy
    :initform 1d0)
   (history-stress
    :accessor mp-history-stress
    :initform 0d0)
   (delay-time
    :accessor mp-delay-time
    :initform 1d0
    :initarg :delay-time
    )
   (ft
    :accessor mp-ft
    :initform 200d3
    :initarg :ft)
   (fc
    :accessor mp-fc
    :initform 200d3
    :initarg :fc)
   )
  (:documentation "A chalk damage model"))

(defclass particle-chalk-delayed (particle-chalk-brittle)
  (
   (delay-time
    :accessor mp-delay-time
    :initform 1d0
    :initarg :delay-time
    )
   )
  (:documentation "A chalk damage model"))
(defclass particle-chalk-anisotropic (particle-chalk-delayed)
  ()
  (:documentation "A chalk damage model"))


(defmethod constitutive-model ((mp particle-chalk) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-accessors ((de mp-elastic-matrix)
                   (stress mp-stress)
                   (stress-u mp-undamaged-stress)
                   (strain mp-strain)
                   (damage mp-damage)
                   (pressure mp-pressure)
                   (def mp-deformation-gradient)
                   (pos mp-position)
                   (calc-pressure mp-pressure-func)
                   (rho mp-rho)
                   (ps-vm mp-strain-plastic-vm)
                   (plastic-strain mp-strain-plastic)
                   (yield-func mp-yield-func)
                   (enable-plasticity mp-enable-plasticity)
                   )
      mp
    (declare (function calc-pressure))
    ;;Train elastic strain - plus trail kirchoff stress
    (setf stress-u
          (cl-mpm/constitutive::linear-elastic-mat strain de))
    ;; (call-next-method)

    (if enable-plasticity
        (let* ((rho_0 rho)
               (rho_1 rho_0)
               (rho-d (+ rho_1 (* (- rho_0 rho_1) (- 1d0 damage))))
               ;; (rho-d rho_0)
               )
          (multiple-value-bind (sig eps-e f) (cl-mpm/constitutive::vm-plastic stress-u de strain rho)
            (setf stress
                  sig
                  plastic-strain (magicl:.- strain eps-e)
                  yield-func f
                  )
            (setf strain eps-e)
            )
          (incf ps-vm
                (multiple-value-bind (l v)
                    (cl-mpm/utils:eig (cl-mpm/utils:voigt-to-matrix (cl-mpm/particle::mp-strain-plastic mp)))
                  (destructuring-bind (s1 s2 s3) l
                    (sqrt
                     (/ (+ (expt (- s1 s2) 2d0)
                           (expt (- s2 s3) 2d0)
                           (expt (- s3 s1) 2d0)
                           ) 2d0))))))
        (setf stress (magicl:scale stress-u 1d0))
        )
    (when (> damage 0.0d0)
      (let* ((j (magicl:det def))
             (degredation (expt (- 1d0 damage) 2d0))
             )
        ;; (setf stress (magicl:.+ (cl-mpm/constitutive::voight-eye p)
        ;;                         (magicl:scale! s (max 1d-9 degredation))
        ;;                         ))
        ;; (setf stress (magicl:scale stress (max 1d-9 degredation)))
        ;; (multiple-value-bind (l v) (cl-mpm/utils::eig
        ;;                             (magicl:scale! (voight-to-matrix stress) (/ 1d0 j)))
        ;;   (let* ((tp 0d0)
        ;;          ;(tp (funcall calc-pressure pos))
        ;;          (driving-pressure (* tp 1d0 (expt (min 1.00d0 damage) 1)))
        ;;          )
        ;;     (setf pressure tp)
        ;;     (loop for i from 0 to 2
        ;;           do
        ;;              (let* ((sii (nth i l))
        ;;                       (esii (- sii driving-pressure)))
        ;;                  (when (> esii 0d0)
        ;;                    ;;tensile damage -> unbounded
        ;;                    (setf (nth i l) (* sii (max 1d-8 degredation)))
        ;;                    ;; (setf (nth i l) (+ (nth i l) driving-pressure))
        ;;                    ;; (setf (nth i l) (* sii degredation))
        ;;                    )
        ;;                ;; (setf (nth i l) (* sii (max 1d-5 degredation)))
        ;;                  ;; (when (< esii 1d0)
        ;;                  ;;   ;;bounded compressive damage
        ;;                  ;;   (setf (nth i l) (* (nth i l) (max 1d-1 degredation)))
        ;;                  ;;   ;; (setf (nth i l) (+ (nth i l) driving-pressure))
        ;;                  ;;   )
        ;;                ;; (setf (nth i l) 0)
        ;;                  ;; (setf (nth i l) (* sii (max 0d0 (- 1d0 damage))))
        ;;                  )
        ;;           )
        ;;       (setf stress (magicl:scale! (matrix-to-voight (magicl:@ v
        ;;                                                               (magicl:from-diag l :type 'double-float)
        ;;                                                               (magicl:transpose v))) j))
        ;;       ))
        (let ((p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
              (s (cl-mpm/constitutive::deviatoric-voigt stress)))
          (setf stress (magicl:.+ (cl-mpm/constitutive::voight-eye p)
                                  (magicl:scale! s (max 0d-3 degredation))
                                  )))
        ))
    stress
    ))

(defmethod constitutive-model ((mp particle-chalk-delayed) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-accessors ((de mp-elastic-matrix)
                   (stress mp-stress)
                   (stress-u mp-undamaged-stress)
                   (strain mp-strain)
                   (strain-rate mp-strain-rate)
                   (damage mp-damage)
                   (pressure mp-pressure)
                   (def mp-deformation-gradient)
                   (pos mp-position)
                   (calc-pressure mp-pressure-func)
                   (coheasion mp-c)
                   (ps-vm mp-strain-plastic-vm)
                   (plastic-strain mp-strain-plastic)
                   (yield-func mp-yield-func)
                   (enable-plasticity mp-enable-plasticity)
                   (D mp-stretch-tensor)
                   (vorticity mp-vorticity)
                   (E mp-E)
                   (nu mp-nu)
                   (phi mp-phi)
                   (psi mp-psi)
                   (kc-r mp-k-compressive-residual-ratio)
                   (kt-r mp-k-tensile-residual-ratio)
                   (g-r mp-shear-residual-ratio)
                   )
      mp
    (declare (function calc-pressure))
    ;;Train elastic strain - plus trail kirchoff stress
    (setf stress-u
          (cl-mpm/constitutive::linear-elastic-mat strain de))

    (if enable-plasticity
        (progn
          (multiple-value-bind (sig eps-e f)
              (cl-mpm/constitutive::mc-plastic stress-u
                                               de
                                               strain
                                               E
                                               nu
                                               phi
                                               psi
                                               coheasion)
            (setf stress
                  sig
                  plastic-strain (magicl:.- strain eps-e)
                  yield-func f)
            (setf strain eps-e))
          (incf ps-vm
                (multiple-value-bind (l v)
                    (cl-mpm/utils:eig (cl-mpm/utils:voigt-to-matrix (cl-mpm/particle::mp-strain-plastic mp)))
                  (destructuring-bind (s1 s2 s3) l
                    (sqrt
                     (/ (+ (expt (- s1 s2) 2d0)
                           (expt (- s2 s3) 2d0)
                           (expt (- s3 s1) 2d0)
                           ) 2d0)))))))
        (setf stress (magicl:scale stress-u 1d0))
    (when (> damage 0.0d0)
      (let* ((j (magicl:det def))
             (exponential 1)
             (degredation (expt (- 1d0 damage) 1d0))
             )
        (declare (double-float damage))
        (let ((p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
              (s (cl-mpm/constitutive::deviatoric-voigt stress)))
          (setf p
                (if (> p 0d0)
                    (* (- 1d0 (* (- 1d0 kt-r) damage)) p)
                    (* (- 1d0 (* (- 1d0 kc-r) damage)) p)))
          (setf stress (magicl:.+ (cl-mpm/constitutive::voight-eye p)
                                  (magicl:scale! s (- 1d0 (* (- 1d0 g-r) damage))))))
      ))
    stress
    ))

(defun matrix-sqrt (mat)
  (multiple-value-bind (l v) (cl-mpm/utils::eig mat)
    (magicl:@ v
              (cl-mpm/utils::matrix-from-list
               (list
                (the double-float (sqrt (the double-float (nth 0 l)))) 0d0 0d0
                0d0 (the double-float (sqrt (the double-float (nth 1 l)))) 0d0
                0d0 0d0 (the double-float (sqrt (the double-float (nth 2 l))))
                ))
              (magicl:transpose v)
              )))

(defmethod constitutive-model ((mp particle-chalk-anisotropic) strain dt)
  (declare (optimize (speed 0) (debug 3)))
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-accessors ((de mp-elastic-matrix)
                   (stress mp-stress)
                   (stress-u mp-undamaged-stress)
                   (strain mp-strain)
                   (damage mp-damage)
                   (pressure mp-pressure)
                   (def mp-deformation-gradient)
                   (pos mp-position)
                   (calc-pressure mp-pressure-func)
                   (rho mp-rho)
                   (ps-vm mp-strain-plastic-vm)
                   (plastic-strain mp-strain-plastic)
                   (yield-func mp-yield-func)
                   (enable-plasticity mp-enable-plasticity)
                   (damage-tensor mp-damage-tensor)
                   (E mp-E)
                   (nu mp-nu)
                   (coheasion mp-c)
                   (phi mp-phi)
                   (psi mp-psi)
                   )
      mp
    (declare (function calc-pressure))
    ;;Train elastic strain - plus trail kirchoff stress
    (setf stress-u
          (cl-mpm/constitutive::linear-elastic-mat strain de))
    ;; (call-next-method)

    (if enable-plasticity
        (progn
          (multiple-value-bind (sig eps-e f)
              ;; (cl-mpm/constitutive::vm-plastic stress-u de strain coheasion)
            (cl-mpm/constitutive::mc-plastic stress-u de strain
                                             E
                                             nu
                                             phi
                                             psi
                                             coheasion)
            (setf stress
                  sig
                  plastic-strain (magicl:.- strain eps-e)
                  yield-func f
                  )
            (setf strain eps-e)
            )
          (incf ps-vm
                (multiple-value-bind (l v)
                    (cl-mpm/utils:eig (cl-mpm/utils:voigt-to-matrix (cl-mpm/particle::mp-strain-plastic mp)))
                  (destructuring-bind (s1 s2 s3) l
                    (sqrt
                     (/ (+ (expt (- s1 s2) 2d0)
                           (expt (- s2 s3) 2d0)
                           (expt (- s3 s1) 2d0)
                           ) 2d0))))))
        (setf stress (magicl:scale stress-u 1d0))
        )
    ;; (setf stress (magicl:scale stress-u 1d0))
    (when (> damage 0.0d0)
      (let* ((j (magicl:det def))
             (degredation (expt (- 1d0 damage) 1d0))
             (min-strength 1d-6)
             )

        (let* ((min-damage (- 1d0 0d-6))
               (d1 (* min-damage (magicl:tref damage-tensor 0 0)))
               (d2 (* min-damage (magicl:tref damage-tensor 1 1)))
               (d3 (* min-damage (magicl:tref damage-tensor 2 2)))
               (d11 (- 1d0 d1))
               (d22 (- 1d0 d2))
               (d33 (- 1d0 d3))
               (d12 (sqrt (* d11 d22)))
               (d13 (sqrt (* d11 d33)))
               (d23 (sqrt (* d22 d33)))
               (g12 d12)
               (g13 d13)
               (g23 d23)
               (damage-effect (cl-mpm/utils:tensor-voigt-4th-from-list
                               (list
                                d11 d12 d13 0d0 0d0 0d0
                                d12 d22 d23 0d0 0d0 0d0
                                d13 d23 d33 0d0 0d0 0d0
                                0d0 0d0 0d0 g23 0d0 0d0
                                0d0 0d0 0d0 0d0 g13 0d0
                                0d0 0d0 0d0 0d0 0d0 g12
                                )))
               (new-de (magicl:.* damage-effect de))
               )
          ;; (pprint de)
          ;; (pprint new-de)
          ;; (pprint stress)
          (setf stress (magicl:@ new-de strain))
          ;; (pprint stress)
          ;; (break)
          )

        ;; (let* ((b (magicl:.-
        ;;            (magicl:from-diag (list 1d0 1d0 1d0) '(3 3))
        ;;           damage-tensor))
        ;;        (bsqrt (matrix-sqrt b)))
        ;;   )
        ;; (multiple-value-bind (l v)
        ;;     (cl-mpm/utils:eig damage-tensor)
        ;;   (let* ((d1 (- 1d0 (nth 0 l)))
        ;;          (d2 (- 1d0 (nth 1 l)))
        ;;          (d3 (- 1d0 (nth 2 l)))
        ;;          (m (magicl:from-diag
        ;;              (list d1
        ;;                    d2
        ;;                    d3
        ;;                    (sqrt (* d2 d3))
        ;;                    (sqrt (* d3 d1))
        ;;                    (sqrt (* d1 d2)))))
        ;;          (Q
        ;;            (magicl:transpose!
        ;;             (magicl:block-matrix (list
        ;;                                   (magicl:.* (magicl:column v 0)
        ;;                                              (magicl:column v 0))

        ;;                                   (magicl:.* (magicl:column v 1)
        ;;                                              (magicl:column v 1))
        ;;                                   (magicl:.* (magicl:column v 2)
        ;;                                              (magicl:column v 2))

        ;;                                   (magicl:scale!
        ;;                                    (magicl:.* (magicl:column v 0)
        ;;                                               (magicl:column v 1)) 2d0)
        ;;                                   (magicl:scale!
        ;;                                    (magicl:.* (magicl:column v 1)
        ;;                                               (magicl:column v 2)) 2d0)
        ;;                                   (magicl:scale!
        ;;                                    (magicl:.* (magicl:column v 2)
        ;;                                               (magicl:column v 0)) 2d0)

        ;;                                   (magicl:.* (magicl:column v 0)
        ;;                                              (cl-mpm/constitutive::rotate-vector (magicl:column v 0)))
        ;;                                   (magicl:.* (magicl:column v 1)
        ;;                                              (cl-mpm/constitutive::rotate-vector (magicl:column v 1)))
        ;;                                   (magicl:.* (magicl:column v 2)
        ;;                                              (cl-mpm/constitutive::rotate-vector (magicl:column v 2)))

        ;;                                   (magicl:scale! (magicl:.+
        ;;                                                   (magicl:.* (magicl:column v 0)
        ;;                                                              (cl-mpm/constitutive::rotate-vector (magicl:column v 1)))
        ;;                                                   (magicl:.* (magicl:column v 1)
        ;;                                                              (cl-mpm/constitutive::rotate-vector (magicl:column v 0)))) 1d0)

        ;;                                   (magicl:scale! (magicl:.+
        ;;                                                   (magicl:.* (magicl:column v 1)
        ;;                                                              (cl-mpm/constitutive::rotate-vector (magicl:column v 2)))
        ;;                                                   (magicl:.* (magicl:column v 2)
        ;;                                                              (cl-mpm/constitutive::rotate-vector (magicl:column v 1)))) 1d0)
        ;;                                   (magicl:scale! (magicl:.+
        ;;                                                   (magicl:.* (magicl:column v 2)
        ;;                                                              (cl-mpm/constitutive::rotate-vector (magicl:column v 0)))
        ;;                                                   (magicl:.* (magicl:column v 0)
        ;;                                                              (cl-mpm/constitutive::rotate-vector (magicl:column v 2)))) 1d0)
        ;;                                   ) '(2 6))))
        ;;          )
        ;;     (setf stress (magicl:@
        ;;                   (magicl:transpose Q)
        ;;                   m
        ;;                   Q
        ;;                   stress
        ;;                   )
        ;;           )

        ;;     ))
        
        ;; (magicl:scale! stress degredation)
        ;; (multiple-value-bind (l v) (cl-mpm/utils::eig
        ;;                             (magicl:scale! (voight-to-matrix stress) (/ 1d0 j)))
        ;;   (let* ()
        ;;     (loop for i from 0 to 2
        ;;           do
        ;;              (let* ((sii (nth i l))
        ;;                       (esii sii)
        ;;                     (vii (magicl::column v i))
        ;;                     )
        ;;                  (when t;(> esii 0d0)
        ;;                    ;;tensile damage -> unbounded
        ;;                    (let* (
        ;;                           (vsi (magicl:@ vii (magicl:transpose vii)))
        ;;                           (dii (magicl::trace (magicl:@ damage-tensor vsi)))
        ;;                           (degredation (expt (- 1d0 dii) 2d0)))
        ;;                      (setf (nth i l) (* sii (max 1d-12) degredation)))
        ;;                    ;; (setf (nth i l) (* sii (max 1d-8 (- 1d0 damage))))
        ;;                    )
        ;;                  )
        ;;           )
        ;;       (setf stress (magicl:scale! (matrix-to-voight (magicl:@ v
        ;;                                                               (magicl:from-diag l :type 'double-float)
        ;;                                                               (magicl:transpose v))) j))
        ;;       ))
        ))
    stress
    ))

(defclass particle-concrete (particle-elastic-damage)
  ((fracture-energy
    :accessor mp-gf
    :initarg :fracture-energy
    :initform 1d0)
   (ductility
    :accessor mp-ductility
    :initarg :ductility
    :initform 0d0)
   (history-stress
    :accessor mp-history-stress
    :initform 0d0)
   (dissipated-energy
    :accessor mp-dissipated-energy
    :initform 0d0)
   (dissipated-energy-inc
    :accessor mp-dissipated-energy-inc
    :initform 0d0)
   )
  (:documentation "A concrete damage model"))

;; (defclass particle-elastic-damage-ideal (particle-concrete)
;;   )

(defclass particle-limestone (particle-concrete)
  (
   (compression-ratio
    :accessor mp-compression-ratio
    :initarg :compression-ratio
    :initform 1d0)
   (interal-length
    :accessor mp-internal-length
    :type DOUBLE-FLOAT
    :initarg :internal-length
    :initform 1d0)
   )
  (:documentation "A concrete damage model"))

(defclass particle-limestone-delayed (particle-limestone)
  ((delay-time
    :accessor mp-delay-time
    :initform 1d0
    :initarg :delay-time
    ))
  (:documentation "A time dependant limestone elastic damage model"))

(defmethod constitutive-model ((mp particle-limestone) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (stress stress)
               (stress-undamaged undamaged-stress)
               (strain-rate strain-rate)
               (D stretch-tensor)
               (velocity-rate velocity-rate)
               (vorticity vorticity)
               (def deformation-gradient)
               (damage damage)
               (pressure pressure)
               ;; (datum pressure-datum)
               ;; (rho pressure-head)
               (pos position)
               (calc-pressure pressure-func)
               )
      mp
    (declare (double-float pressure damage)
             (function calc-pressure))
    ;; Non-objective stress intergration
    (setf stress-undamaged (cl-mpm/constitutive::linear-elastic-mat strain de))

    (setf stress (magicl:scale stress-undamaged 1d0))
    (when (> damage 0.0d0)
      (let ((degredation (max 1d-9 (expt (- 1d0 damage) 1d0))))
        (magicl:scale! stress degredation)))
    stress
    ))

(defmethod constitutive-model ((mp particle-concrete) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((E E)
               (nu nu)
               (de elastic-matrix)
               (stress stress)
               (stress-undamaged undamaged-stress)
               (strain-rate strain-rate)
               (D stretch-tensor)
               (velocity-rate velocity-rate)
               (vorticity vorticity)
               (def deformation-gradient)
               (damage damage)
               (pressure pressure)
               ;; (datum pressure-datum)
               ;; (rho pressure-head)
               (pos position)
               (calc-pressure pressure-func)
               )
      mp
    (declare (double-float pressure damage)
             (function calc-pressure))
    ;; Non-objective stress intergration
    (setf stress-undamaged (cl-mpm/constitutive::linear-elastic-mat strain de))

    (setf stress (magicl:scale stress-undamaged 1d0))
    (when (> damage 0.0d0)
      (let ((degredation (expt (- 1d0 damage) 1d0)))
        (magicl:scale! stress (max 0d-9 degredation)))
      ;; (let* ((j 1d0)
      ;;        (p (/ (cl-mpm/constitutive::voight-trace stress) 3d0))
      ;;        (s (cl-mpm/constitutive::deviatoric-voigt stress)))
      ;;   ;; (setf stress (magicl:.+ (cl-mpm/constitutive::voight-eye p)
      ;;   ;;                         (magicl:scale! s (max 1d-3 (expt (- 1d0 damage) 2d0)))
      ;;   ;;                         ))
      ;;   (multiple-value-bind (l v) (cl-mpm/utils::eig
      ;;                               (magicl:scale! (voight-to-matrix stress) (/ 1d0 j)))
      ;;     (let* (;(tp (funcall calc-pressure (magicl:tref pos 1 0) datum rho))
      ;;            (tp (funcall calc-pressure pos))
      ;;            (driving-pressure (* tp 1d0 (expt (min 1.00d0 damage) 1)))
      ;;            (degredation (expt (- 1d0 damage) 2d0)))
      ;;       ;; (setf stress (magicl:scale stress-undamaged (max 1d-3 degredation)))
      ;;       (setf pressure tp)
      ;;       (loop for i from 0 to 2
      ;;             do
      ;;                (let* ((sii (nth i l))
      ;;                         (esii (- sii driving-pressure)))
      ;;                    (when (> esii 0d0)
      ;;                      ;;Tensile damage -> unbounded
      ;;                      (setf (nth i l) (* esii (max 0d-6 degredation)))
      ;;                      (setf (nth i l) (+ (nth i l) driving-pressure))
      ;;                      )
      ;;                    ;; (when (< esii 1d0)
      ;;                    ;;   ;;Bounded compressive damage
      ;;                    ;;   (setf (nth i l) (* esii (max 1d0 degredation)))
      ;;                    ;;   (setf (nth i l) (+ (nth i l) driving-pressure))
      ;;                    ;;   )
      ;;                    ;; (setf (nth i l) (* sii (max 0d0 (- 1d0 damage))))
      ;;                    )
      ;;             )
      ;;         (setf stress (magicl:scale! (matrix-to-voight (magicl:@ v
      ;;                                                                 (magicl:from-diag l :type 'double-float)
      ;;                                                                 (magicl:transpose v))) j))
      ;;         ))
      ;;   )
      )
    stress
    ))
