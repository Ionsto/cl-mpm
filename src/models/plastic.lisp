;; (defpackage :cl-mpm/models/plastic
;;   (:use :cl
;;    :cl-mpm/utils)
;;   (:export
;;    )
;;   )
;; (in-package :cl-mpm/models/plastic)
(in-package :cl-mpm/particle)

(defclass particle-plastic (particle)
  ((enable-plasticity
    :accessor mp-enable-plasticity
    :initarg :enable-plasticity
    :initform t
    )
   (strain-plastic-vm
    :accessor mp-strain-plastic-vm
    :type DOUBLE-FLOAT
    :initform 0d0)
   (strain-plastic-vm-inc
    :accessor mp-strain-plastic-vm-inc
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
   ))

(defclass particle-vm (particle-elastic particle-plastic)
  ((rho
    :accessor mp-rho
    :initarg :rho
    )
   (yield-func
    :accessor mp-yield-func
    :type double-float
    :initform 0d0))
  (:documentation "A vm perfectly plastic material point"))

(defclass particle-mc (particle-elastic particle-plastic)
  ((phi
    :accessor mp-phi
    :initarg :phi
    :initform 0d0
    )
   (psi
    :accessor mp-psi
    :initarg :psi
    :initform 0d0
    )
   (c
    :accessor mp-c
    :initarg :c
    :initform 0d0
    )
   (phi-0
    :accessor mp-phi-0
    :initarg :phi
    :initform 0d0
    )
   (psi-0
    :accessor mp-psi-0
    :initarg :psi
    :initform 0d0
    )
   (c-0
    :accessor mp-c-0
    :initarg :c
    :initform 0d0
    )
   (phi-r
    :accessor mp-phi-r
    :initarg :phi-r
    :initform 0d0
    )
   (psi-r
    :accessor mp-psi-r
    :initarg :psi-r
    :initform 0d0
    )
   (c-r
    :accessor mp-c-r
    :initarg :c-r
    :initform 0d0
    )
   (softening
    :accessor mp-softening
    :initarg :softening
    :initform 0d0))

  (:documentation "A mohr-coloumb perfectly plastic material point"))

(defclass particle-dp (particle-mc)
  ())

(defmethod constitutive-model ((mp particle-vm) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-accessors ((de mp-elastic-matrix)
                   (stress mp-stress)
                   (rho mp-rho)
                   (plastic-strain mp-strain-plastic)
                   (ps-vm mp-strain-plastic-vm)
                   (strain mp-strain)
                   (yield-func mp-yield-func)
                   (enable-plasticity mp-enable-plasticity))
      mp
    ;;Train elastic strain - plus trail kirchoff stress
    (setf stress
          (cl-mpm/constitutive::linear-elastic-mat strain de))
    (when enable-plasticity
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
                       ) 2d0))))))
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
                   (soft mp-softening)
                   (enabled mp-enable-plasticity)
                   )
      mp
    (declare (double-float soft ps-vm E nu phi psi c))
    ;;Train elastic strain - plus trail kirchoff stress
    (setf stress (cl-mpm/constitutive::linear-elastic-mat strain de stress))
    (when enabled
      (multiple-value-bind (sig eps-e f) (cl-mpm/constitutive::mc-plastic stress de strain E nu phi psi c)
        (setf stress
              sig
              plastic-strain (cl-mpm/fastmath:fast-.+ plastic-strain (magicl:.- strain eps-e) plastic-strain)
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
                       ) 2d0))))))
    (when (> soft 0d0)
      (with-accessors ((c-0 mp-c-0)
                       (phi-0 mp-phi-0)
                       (psi-0 mp-psi-0)
                       (c-r mp-c-r)
                       (phi-r mp-phi-r)
                       (psi-r mp-psi-r))
          mp
        (declare (double-float c-0 c-r phi-0 phi-r psi-0 psi-r))
          (setf
           c (+ c-r (* (- c-0 c-r) (exp (- (* soft ps-vm)))))
           phi (+ phi-r (* (- phi-0 phi-r) (exp (- (* soft ps-vm)))))))
      )
    stress
    ))

(defmethod constitutive-model ((mp particle-dp) strain dt)
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
                   (soft mp-softening)
                   (enabled mp-enable-plasticity)
                   )
      mp
    (declare (double-float soft ps-vm E nu phi psi c))
    ;;Train elastic strain - plus trail kirchoff stress
    (cl-mpm/constitutive::linear-elastic-mat strain de stress)
    (when enabled
      (multiple-value-bind (sig eps-e f) (cl-mpm/constitutive::plastic-dp strain de E nu phi psi c)
        (setf stress
              sig
              plastic-strain (cl-mpm/fastmath:fast-.+ plastic-strain (magicl:.- strain eps-e) plastic-strain)
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
                       ) 2d0))))))
    (when (> soft 0d0)
      (with-accessors ((c-0 mp-c-0)
                       (phi-0 mp-phi-0)
                       (psi-0 mp-psi-0)
                       (c-r mp-c-r)
                       (phi-r mp-phi-r)
                       (psi-r mp-psi-r))
          mp
        (declare (double-float c-0 c-r phi-0 phi-r psi-0 psi-r))
          (setf
           c (+ c-r (* (- c-0 c-r) (exp (- (* soft ps-vm)))))
           phi (+ phi-r (* (- phi-0 phi-r) (exp (- (* soft ps-vm)))))))
      )
    stress
    ))
