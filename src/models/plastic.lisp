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
   (strain-plastic-vm-1
    :accessor mp-strain-plastic-vm-1
    :type DOUBLE-FLOAT
    :initform 0d0)
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
   (softening
    :accessor mp-softening
    :initarg :softening
    :initform 0d0)
   (plastic-iterations
    :accessor mp-plastic-iterations
    :initform 0)))

(defclass particle-vm (particle-elastic particle-plastic)
  ((rho
    :accessor mp-rho
    :initarg :rho
    )
   (rho-0
    :accessor mp-rho-0
    :initarg :rho
    )
   (rho-r
    :accessor mp-rho-r
    :initarg :rho-r
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
   )

  (:documentation "A mohr-coloumb perfectly plastic material point"))

(defclass particle-dp (particle-mc)
  ())

(defmethod constitutive-model ((mp particle-vm) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-accessors ((de mp-elastic-matrix)
                   (stress mp-stress)
                   (rho mp-rho)
                   (soft mp-softening)
                   (plastic-strain mp-strain-plastic)
                   (ps-vm mp-strain-plastic-vm)
                   (ps-vm-1 mp-strain-plastic-vm-1)
                   (ps-vm-inc mp-strain-plastic-vm-inc)
                   (strain mp-strain)
                   (strain-n mp-strain-n)
                   (yield-func mp-yield-func)
                   (enable-plasticity mp-enable-plasticity))
      mp
    ;;Train elastic strain - plus trail kirchoff stress
    (cl-mpm/constitutive::linear-elastic-mat strain de stress)
    (when enable-plasticity
      (multiple-value-bind (sig eps-e f) (cl-mpm/constitutive::vm-plastic stress de strain rho)
        (setf stress
              sig
              plastic-strain (cl-mpm/fastmaths:fast-.- strain-n eps-e)
              strain eps-e
              yield-func f
              ))
      (setf ps-vm-inc (sqrt (cl-mpm/constitutive::voigt-j2 (cl-mpm/utils:deviatoric-voigt (cl-mpm/particle::mp-strain-plastic mp))))
            ;; (multiple-value-bind (l v)
            ;;     (cl-mpm/utils:eig (cl-mpm/utils:voigt-to-matrix ))
            ;;   (destructuring-bind (s1 s2 s3) l
            ;;     (sqrt
            ;;      (/ (+ (expt (- s1 s2) 2d0)
            ;;            (expt (- s2 s3) 2d0)
            ;;            (expt (- s3 s1) 2d0)
            ;;            ) 2d0))))
            )
      (setf ps-vm (+ ps-vm-1 ps-vm-inc))
      )
    (when (> soft 0d0)
      (with-accessors ((rho-r mp-rho-r)
                       (rho-0 mp-rho-0))
          mp
        (declare (double-float rho-0 rho-r))
        (setf
         rho (+ rho-r (* (- rho-0 rho-r) (exp (- (* soft ps-vm))))))))
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
                   (ps-vm-1 mp-strain-plastic-vm-1)
                   (ps-vm-inc mp-strain-plastic-vm-inc)
                   (strain mp-strain)
                   (yield-func mp-yield-func)
                   (soft mp-softening)
                   (enabled mp-enable-plasticity)
                   )
      mp
    (declare (double-float soft ps-vm ps-vm-1 ps-vm-inc E nu phi psi c))
    ;;Train elastic strain - plus trail kirchoff stress
    (setf stress (cl-mpm/constitutive::linear-elastic-mat strain de stress))
    (when enabled
      (let ((f-r t))
        (loop for i from 0 to 0;to 50
              while f-r
              do
                 (progn
                   (multiple-value-bind (sig eps-e f inc)
                       (cl-mpm/ext::constitutive-mohr-coulomb stress
                                                              de
                                                              strain
                                                              E
                                                              nu
                                                              phi
                                                              psi
                                                              c)
                     (declare (double-float f inc))
                     (setf f-r (> f 1d-5))
                     (setf
                      stress sig
                      strain eps-e
                      yield-func f)
                     (setf ps-vm-inc inc)
                     (setf ps-vm (+ ps-vm-1 inc)))
                   (setf (mp-plastic-iterations mp) i)
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
                        phi (atan (+ (tan phi-r) (* (- (tan phi-0) (tan phi-r)) (exp (- (* soft ps-vm))))))))
                     )))))
    stress))

(defmethod cl-mpm/particle::reset-loadstep-mp ((mp particle-plastic))
  (with-accessors ((ps-vm-1 mp-strain-plastic-vm-1)
                   (ps-vm mp-strain-plastic-vm)
                   )
      mp
    (declare (double-float ps-vm ps-vm-1))
    (setf ps-vm ps-vm-1))
  (call-next-method))
(defmethod cl-mpm/particle::new-loadstep-mp ((mp particle-plastic))
  (with-accessors ((ps-vm-1 mp-strain-plastic-vm-1)
                   (ps-vm-inc mp-strain-plastic-vm-inc)
                   (ps-vm mp-strain-plastic-vm)
                   )
      mp
    (declare (double-float ps-vm ps-vm-1 ps-vm-inc))
    (setf ps-vm-1 ps-vm))
  (call-next-method))


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
      (multiple-value-bind (sig eps-e f yield)
          (cl-mpm/constitutive::plastic-dp stress de strain E nu phi psi c)
          ;; (cl-mpm/ext::constitutive-drucker-prager strain de E nu phi psi c)
        (if yield
            (progn
              (setf stress sig
                    plastic-strain (cl-mpm/fastmaths:fast-.- strain eps-e)
                    yield-func f
                    )
              (setf strain eps-e))
          (progn
            (cl-mpm/fastmaths:fast-zero plastic-strain)
            (setf yield-func 0d0))))
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
           phi (atan (+ (tan phi-r) (* (- (tan phi-0) (tan phi-r)) (exp (- (* soft ps-vm)))))))))
    stress))

;; (defmethod new-loadstep-mp ((mp particle-plastic))
;;   (with-accessors ((strain cl-mpm/particle:mp-strain)
;;                    (strain-n cl-mpm/particle:mp-strain-n)
;;                    (disp cl-mpm/particle::mp-displacement)
;;                    (def    cl-mpm/particle:mp-deformation-gradient)
;;                    (def-0 cl-mpm/particle::mp-deformation-gradient-0)
;;                    (df-inc    cl-mpm/particle::mp-deformation-gradient-increment)
;;                    )
;;       mp
;;     (cl-mpm/utils:matrix-copy-into def def-0)
;;     (cl-mpm/utils:matrix-copy-into (cl-mpm/utils:matrix-eye 1d0) df-inc)
;;     (cl-mpm/utils:voigt-copy-into strain strain-n)
;;     (cl-mpm/fastmaths:fast-zero disp)
;;     (break)
;;     ))
