(in-package :cl-mpm/particle)
(declaim #.cl-mpm/settings:*optimise-setting*)
(defclass particle-ssa (particle)
  ((height
    :initarg :height
    :initform 1d0
    :accessor mp-height)))

(defclass particle-ssa-elastic (particle-ssa particle-elastic)
  ())


(defclass particle-ssa-glen (particle-ssa-elastic)
  ((visc-factor
    :initform 111d6
    :accessor mp-visc-factor
    :initarg :visc-factor)
   (visc-power
    :initform 3d0
    :accessor mp-visc-power
    :initarg :visc-power)
   (true-visc
    :accessor mp-true-visc
    :initform 0d0)))

(defmethod cl-mpm::update-stress-mp (mesh (mp particle-ssa) dt fbar)
  (cl-mpm::update-stress-linear mesh mp dt fbar))

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-ssa) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  (cl-mpm::update-domain-corner mesh mp dt)
  ;; (cl-mpm::update-domain-stretch mesh mp dt)
  ;; (cl-mpm::scale-domain-size mesh mp)
  )



(defun glen-visco (stress visc-factor visc-power)
  (let* ((dev (cl-mpm/utils:deviatoric-voigt stress))
         (effective-stress (sqrt (* 1/2 (cl-mpm/fastmaths::voigt-j2 dev))))
         (visc-factor (expt visc-factor (- visc-power)))
         )
    (/ 1d0
       (+ 1d-20
          (* visc-factor 2 (expt effective-stress (- visc-power 1)))))))

(defun ssa-maxwell-glen (strain-inc stress E nu visc dt)
  (let* ((p (+ (varef stress 0)
               (varef stress 1)))
         (K (/ E (* 3 (- 1 (* 2 nu)))))
         (G (/ E (* 2 (+ 1 nu))))
         (relaxation-const (/ (* dt G) visc))
         (dp (* K (* 0.5d0
                     (varef strain-inc 0)
                     (varef strain-inc 1)))))
    (cl-mpm/fastmaths:fast-.+
     stress
     (cl-mpm/utils:voigt-from-list (list dp dp dp 0d0 0d0 0d0))
     stress)
    (cl-mpm/fastmaths:fast-.-
     stress
     stress
     )
    ))

(defmethod cl-mpm::constitutive-model ((mp cl-mpm/particle::particle-ssa-glen) strain dt)
  "Function for modeling stress intergrated viscoplastic norton-hoff material"
  (with-accessors ((E mp-E)
                   (nu mp-nu)
                   ;; (de elastic-matrix)
                   (visc-factor mp-visc-factor)
                   (visc-power mp-visc-power)
                   ;; (visc-power visc-power)
                   ;; (strain-rate strain-rate) ;Note strain rate is actually strain increment through dt
                   ;; (velocity-rate velocity-rate) ;Note strain rate is actually strain increment through dt
                   (stress mp-stress)
                   (def mp-deformation-gradient)
                   ;; (strain strain)
                   ;; (stretch stretch-tensor)
                   (de mp-elastic-matrix)
                   (strain-inc mp-strain-rate)
                   )
      mp
    (let ((visc
            (glen-visco (cl-mpm/fastmaths:fast-scale-voigt
                         stress
                         (/ 1d0 (magicl:det def)))
                        visc-factor
                        visc-power)))
      (ssa-maxwell-glen strain-inc stress E nu visc dt))))
