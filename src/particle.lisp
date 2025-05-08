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
    #:mp-acceleration
    #:mp-stress
    #:mp-strain
    #:mp-strain-n
    #:mp-strain-rate
    #:mp-vorticity
    #:mp-gravity
    #:mp-body-force
    #:mp-damage
    #:mp-temperature
    #:mp-heat-capacity
    #:mp-critical-stress
    #:mp-deformation-gradient
    #:mp-penalty-contact
    #:mp-penalty-frictional-force
    #:constitutive-model
    #:particle
    #:particle-damage
    #:post-stress-step
    ))

(in-package :cl-mpm/particle)

(declaim (optimize (debug 0) (safety 0) (speed 3)))

(declaim (inline make-node-cache))
(defstruct (node-cache
            (:constructor make-node-cache
                (node
                 weight
                 grad-x
                 grad-y
                 grad-z
                 weight-fbar
                 grad-fbar-x
                 grad-fbar-y
                 grad-fbar-z)))
  (node nil :type (or null cl-mpm/mesh::node))
  (weight 0d0 :type double-float)
  (grad-x 0d0 :type double-float)
  (grad-y 0d0 :type double-float)
  (grad-z 0d0 :type double-float)
  (weight-fbar 0d0 :type double-float)
  (grad-fbar-x 0d0 :type double-float)
  (grad-fbar-y 0d0 :type double-float)
  (grad-fbar-z 0d0 :type double-float)
  )

;(defun make-node-cache (&key node weight grads weight-fbar grads-fbar))

;; (defun push-mp-cache-node (mp &key node weight grads weight-fbar grads-fbar)
;;   (declare (particle mp))
;;   (with-slots ((nc cached-nodes))
;;       mp
;;     (when (< (fill-pointer nc) (array-total-size nc))
;;       ()
;;       ))

;;   )

(defclass particle ()
  ((position
    :accessor mp-position
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils:vector-zeros)
    :initarg :position)
   (size
    :accessor mp-domain-size
    :type magicl:matrix/double-float
    :initarg :size
    :initform (cl-mpm/utils:vector-zeros))
   (true-domain
    :accessor mp-true-domain
    :type magicl:matrix/double-float
    :initarg :true-domain
    :initform (cl-mpm/utils::matrix-eye 1d0))
   (mass
     :accessor mp-mass
     :type double-float
     :initarg :mass
     :initform 1d0)
   (velocity
    :accessor mp-velocity
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initarg :velocity
    :initform (cl-mpm/utils:vector-zeros))
   (volume
    :accessor mp-volume
    :type double-float
    :initarg :volume
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
   (unique-index
    :accessor mp-unique-index
    :type integer
    :initform -1
    :initarg :unique-index)
   (mpi-index
    :accessor mp-mpi-index
    :type integer
    :initform -1
    :initarg :mpi-index)
   (volume-0
    :accessor mp-volume-0
    :type double-float
    ;;By specifying both :volume and :volume-0 as initargs, it will capture volume by default
    ;;but it also allows for direct specification
    :initarg :volume
    :initarg :volume-0
    :initform 1d0)
   (size-0
    :accessor mp-domain-size-0
    :type magicl:matrix/double-float
    :initarg :size-0
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
   (strain-n
    :accessor mp-strain-n
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils:voigt-zeros))
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
     :initform (cl-mpm/utils:matrix-eye 1d0))
   (deformation-gradient-0
       :accessor mp-deformation-gradient-0
       :type MAGICL:MATRIX/DOUBLE-FLOAT
       :initarg :deformation-gradient
       :initform (cl-mpm/utils:matrix-eye 1d0))
   (deformation-gradient-increment
       :accessor mp-deformation-gradient-increment
       :type MAGICL:MATRIX/DOUBLE-FLOAT
       :initarg :deformation-gradient-increment
       :initform (cl-mpm/utils:matrix-eye 1d0))
   (gravity
     :type double-float
     :accessor mp-gravity
     :initform 0d0;-9.8d0
     :initarg :gravity
     )
   (viscous-damping
    :type double-float
    :accessor mp-viscous-damping
    :initform 0d0
    :initarg :damping
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
   (displacement-increment
    :accessor mp-displacement-increment
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils::vector-zeros))
   (position-trial
    :accessor mp-position-trial
    ;; :initarg :position
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils:vector-zeros))
   (displacement
    :accessor mp-displacement
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils::vector-zeros))
   (cached-nodes
    :accessor mp-cached-nodes
    :initarg :nc
    :initform (make-array 8 :fill-pointer 0 :element-type 'node-cache :initial-element (make-node-cache nil 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0)))
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
   (penalty-frictional-force
    :accessor mp-penalty-frictional-force
    :type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform (cl-mpm/utils::vector-zeros))
   (penalty-normal-force
    :accessor mp-penalty-normal-force
    :type double-float
    :initform 0d0)
   (penalty-friction-stick
    :accessor mp-penalty-friction-stick
    :type boolean
    :initform nil)
   ;;This reports whether we had contact over the course of the last timestep
   (penalty-contact-step
    :accessor mp-penalty-contact-step
    :type boolean
    :initform nil)
   ;;This is for keeping track of whether we will be having contact coming up
   (penalty-contact
    :accessor mp-penalty-contact
    :type boolean
    :initform nil)
   (penalty-energy
    :accessor mp-penalty-energy
    :type double-float
    :initform 0d0
    )
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
    :type double-float
    :accessor mp-E
    :initarg :E)
   (nu
    :type double-float
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

(defun update-elastic-matrix (particle)
  (with-accessors ((de mp-elastic-matrix)
                   (E  mp-E)
                   (nu mp-nu)
                   (p mp-p-modulus))
      particle
    (setf p (/ E (* (+ 1d0 nu) (- 1d0 nu))))
    (setf de (cl-mpm/constitutive::linear-elastic-matrix E nu))))

(defmethod (setf mp-E) :after (value (p particle-elastic))
  (update-elastic-matrix p))
(defmethod (setf mp-nu) :after (value (p particle-elastic))
  (update-elastic-matrix p))

(defmethod (setf mp-domain-size) :after (value (p particle))
  ;; (with-accessors ((domain-true mp-true-domain))
  ;;     p
  ;;   (break)
  ;;   (setf
  ;;    (magicl:tref domain-true 0 0) (varef value 0)
  ;;    (magicl:tref domain-true 1 1) (varef value 1)
  ;;    (magicl:tref domain-true 2 2) (varef value 2)))
  )

(defmethod initialize-instance :after ((p particle) &key)
  (with-accessors ((domain mp-domain-size)
                   (domain-true mp-true-domain)
                   (position mp-position)
                   (position-trial mp-position-trial)
                   )
      p
    (cl-mpm/utils:matrix-copy-into position position-trial)
    (setf
     (magicl:tref domain-true 0 0) (varef domain 0)
     (magicl:tref domain-true 1 1) (varef domain 1)
     (magicl:tref domain-true 2 2) (varef domain 2)))
  )

(defmethod reinitialize-instance :after ((p particle) &key)
  ;; (with-accessors ((domain mp-domain-size)
  ;;                  (domain-true mp-true-domain))
  ;;     p
  ;;   (setf
  ;;    (magicl:tref domain-true 0 0) (varef domain 0)
  ;;    (magicl:tref domain-true 1 1) (varef domain 1)
  ;;    (magicl:tref domain-true 2 2) (varef domain 2)))
  )

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
    :initarg :viscosity))
  (:documentation "A fluid material point"))

(defclass particle-viscoelastic (particle-elastic)
  (
   (viscosity
    :accessor mp-viscosity
    :initarg :viscosity)
   )
  (:documentation "A visco-elastic material point"))

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



;; (defclass particle-thermoelastic-damage (particle-elastic particle-damage particle-thermal)
;;   ()
;;   (:documentation "A mp with elastic mechanics with variable thermal fields"))
;; (defclass particle-thermoviscoplastic-damage (particle-viscoplastic particle-damage particle-thermal)
;;   ()
;;   (:documentation "A mp with viscoplastic mechanics with variable thermal fields"))

;; (defclass particle-thermofluid-damage (particle-fluid particle-damage particle-thermal)
;;   ()
;;   (:documentation "A mp with viscoplastic mechanics with variable thermal fields"))
;; (defclass particle-viscoelastic-fracture (particle-viscoelastic particle-fracture)
;;   ()
;;   (:documentation "A viscoelastic mp with fracture mechanics"))


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

;; (defun make-particle-elastic-damage (nD E nu &key (pos nil) (volume 1) (mass 1))
;;     (let ((p (make-particle nD 'particle-elastic-damage :position pos :volume volume :mass mass
;;                             ) ))
;;         (progn
;;             (setf (mp-E p) E)
;;             (setf (mp-nu p) nu)
;;         p)))

(defgeneric constitutive-model (mp elastic-trial-strain dt)
    (:documentation "Compute new stress state given elastic strain")
    (:method (mp strain dt)
        (magicl:scale strain 0)))

(defmethod constitutive-model ((mp particle-elastic) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((de elastic-matrix))
      mp
    (cl-mpm/constitutive::linear-elastic-mat strain de)))


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
    ;; (cl-mpm/fastmaths::fast-.+
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
    (cl-mpm/fastmaths::fast-.+
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
    (cl-mpm/fastmaths::fast-.+
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
      (cl-mpm/fastmaths::fast-.+
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
    (let* ((eng-strain-rate (cl-mpm/fastmaths::fast-.* (magicl:map (lambda (x) (* x (exp x))) strain) velocity-rate
                                                  (cl-mpm/utils:voigt-from-list '(1d0 1d0 1d0 0.5d0 0.5d0 0.5d0))))
          (viscosity (cl-mpm/constitutive::glen-viscosity-strain eng-strain-rate visc-factor visc-power))
          ;(viscosity (cl-mpm/constitutive::glen-viscosity-stress stress visc-factor visc-power))
          )
      ;; stress
      (setf temp viscosity)
      (cl-mpm/fastmaths::fast-.+
       stress
       (objectify-stress-logspin
        (if (> viscosity 0d0)
            (cl-mpm/constitutive::maxwell strain-rate stress E nu de viscosity dt)
            (cl-mpm/constitutive::linear-elastic-mat strain-rate de))
        stress
        def
        vorticity
        D
        )))))

(defgeneric post-stress-step (mesh mp dt)
  (:documentation "This step gets called after full stress state resolved and allows for other processing"))

(defmethod post-stress-step (mesh mp dt)
  ())

;; (declaim (inline assemble-vorticity-matrix)
;;          (ftype (function (magicl:matrix/double-float)
;;                           magicl:matrix/double-float
;;                           ) assemble-vorticity-matrix))
(defun assemble-vorticity-matrix (vorticity)
  (let ((dx (magicl:tref vorticity 0 0))
        (dy (magicl:tref vorticity 1 0))
        (dxdy (magicl:tref vorticity 2 0))
        )
    (declare (double-float dx dy dxdy))
    (cl-mpm/utils::matrix-from-list (list dx dxdy 0d0
                                          (- dxdy) dy 0d0
                                          0d0 0d0 0d0))))

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

;; (declaim (inline objectify-stress-logspin)
;;          (ftype (function (magicl:matrix/double-float
;;                            magicl:matrix/double-float
;;                            magicl:matrix/double-float
;;                            magicl:matrix/double-float
;;                            magicl:matrix/double-float
;;                            )
;;                           magicl:matrix/double-float
;;                           )
;;                 objectify-stress-logspin))
(defun objectify-stress-logspin (stress-inc stress def vorticity D)
    (let ((b (magicl:@ def (magicl:transpose def)))
          ;; (omega (assemble-vorticity-matrix vorticity))
          (omega (magicl:scale! (magicl:.- D (magicl:transpose D)) 0.5d0))
          (D (magicl:scale! (cl-mpm/fastmaths::fast-.+ D (magicl:transpose D)) 0.5d0))
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
                                ;; (declare (double-float l_i l_j)
                                ;;          (magicl:matrix/double-float v_i v_j))
                                ;;When the eigenvalues are distinct
                                (when (and
                                       ;;When they are nonzero
                                       (> (abs (- l_i l_j)) 1d-6)
                                       )
                                  ;; When we have pairs of unique nonzero eigenvalues
                                  (setf omega
                                        (cl-mpm/fastmaths::fast-.+ omega
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


;; (defclass particle-elastic-damage-ideal (particle-concrete)
;;   )




(defclass particle-damage-fundemental ()
  ())
(defclass particle-damage (particle particle-damage-fundemental)
  ((local-damage
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
   (damage-ybar-prev
    :accessor mp-damage-ybar-prev
    :type DOUBLE-FLOAT
    :initform 0d0
    :initarg :damage-ybar-prev)
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
    :initarg :local-length
    :initform 1d0)
   (damage-position
    :accessor mp-damage-position
    ;:type MAGICL:MATRIX/DOUBLE-FLOAT
    :initform nil
    ;:initform (cl-mpm/utils::vector-zeros)
    )
   (damage-model
    :accessor mp-damage-model
    :initarg :damage-model)
   (enable-damage
    :accessor mp-enable-damage
    :initarg :enable-damage
    :initform t)
   (damage-domain-update-rate
    :accessor mp-damage-domain-update-rate
    :initarg :damage-domain-rate
    :initform 0d0))
  (:documentation "A material point with a damage tensor"))



(defgeneric new-loadstep-mp (mp)
  (:documentation "Finalise a dynamic relaxation loadstep, setting converged values to previous values"))

(defmethod new-loadstep-mp ((mp cl-mpm::particle))
  (with-accessors ((strain cl-mpm/particle:mp-strain)
                   (strain-n cl-mpm/particle:mp-strain-n)
                   (disp cl-mpm/particle::mp-displacement-increment)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (def-0 cl-mpm/particle::mp-deformation-gradient-0)
                   (df-inc    cl-mpm/particle::mp-deformation-gradient-increment)
                   )
      mp
    (cl-mpm/utils:matrix-copy-into def def-0)
    (cl-mpm/utils:matrix-copy-into (cl-mpm/utils:matrix-eye 1d0) df-inc)
    (cl-mpm/utils:voigt-copy-into strain strain-n)
    ;; (cl-mpm/fastmaths:fast-zero disp)
    ))

(defmethod new-loadstep-mp ((mp cl-mpm::particle-damage))
  (with-accessors ((strain cl-mpm/particle:mp-strain)
                   (strain-n cl-mpm/particle:mp-strain-n)
                   (disp cl-mpm/particle::mp-displacement-increment)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (def-0 cl-mpm/particle::mp-deformation-gradient-0)
                   (df-inc    cl-mpm/particle::mp-deformation-gradient-increment)
                   (ybar    cl-mpm/particle::mp-damage-ybar)
                   (ybar-prev    cl-mpm/particle::mp-damage-ybar-prev)
                   )
      mp
    (cl-mpm/utils:matrix-copy-into def def-0)
    (cl-mpm/utils:matrix-copy-into (cl-mpm/utils:matrix-eye 1d0) df-inc)
    (cl-mpm/utils:voigt-copy-into strain strain-n)
    ;; (cl-mpm/fastmaths:fast-zero disp)
    (setf ybar-prev ybar)
    (call-next-method)
    ))



