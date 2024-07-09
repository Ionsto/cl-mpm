(defpackage :cl-mpm/models/limestone
  (:use
   :cl
   :cl-mpm/utils
   :cl-mpm/particle)
  (:export))
(in-package :cl-mpm/models/limestone)

(defclass particle-limestone (particle-concrete particle-mc)
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
   (enable-plasticity
    :accessor mp-enable-plasticity
    :initarg :enable-plasticity
    :initform nil))
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
    (with-accessors
          ((phi mp-phi)
           (psi mp-psi)
           (E mp-E)
           (nu mp-nu)
           (coheasion mp-c)
           (plastic-strain mp-strain-plastic)
           (yield-func mp-yield-func)
           (enable-plasticity mp-enable-plasticity)
           (ps-vm mp-strain-plastic-vm)
           )
        mp
      (declare (double-float pressure damage)
               (function calc-pressure))
      ;; Non-objective stress intergration
      (cl-mpm/constitutive::linear-elastic-mat strain de stress-undamaged)
      ;; (setf stress-undamaged (cl-mpm/constitutive::linear-elastic-mat strain de))
      (if enable-plasticity
          (progn
            (multiple-value-bind (sig eps-e f)
                (cl-mpm/constitutive::mc-plastic stress-undamaged
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
                             ) 2d0))))))
          (setf stress (cl-mpm/utils:voigt-copy stress-undamaged)))
      (when (> damage 0.0d0)
        (cl-mpm/fastmath::fast-scale! stress (- 1d0 (* (- 1d0 1d-9) damage)))))
    stress))
