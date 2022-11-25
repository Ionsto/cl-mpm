(defpackage :cl-mpm/constitutive
  (:use :cl
        :cl-mpm/utils)
  (:export
    #:linear-elastic
    #:maxwell
    #:maxwell-exp
    )
  )
(in-package :cl-mpm/constitutive)
(declaim (optimize (debug 3) (safety 3) (speed 2)))

(defun linear-elastic-matrix (E nu)
  "Create an isotropic linear elastic matrix"
  (magicl:scale  
    (magicl:from-list (list 
                        (- 1d0 nu) nu 0d0 
                        nu (- 1d0 nu) 0d0 
                        0d0 0d0 (- 1d0 (* 2 nu)))
                      '(3 3) :type 'double-float)
    (/ E (* (+ 1 nu) (- 1 nu nu)))))

(defun linear-elastic (strain E nu)
  "Isotropic linear-elastic constitutive model"
   (magicl:@ (linear-elastic-matrix E nu) strain))

(defun eos-cole (density rest-density stiffness power)
  "Cole equation of state for pressure"
  (* stiffness (- (expt (/ density rest-density) power) 1)))

(defun eos-ideal-gas (density rest-density rest-pressure adiabatic-index)
  "Ideal gas law equation of state for pressure"
  (* rest-pressure (expt (/ density rest-density) adiabatic-index)))


(defun newtonian-fluid (strain pressure viscosity)
  "A newtonian fluid model"
  (magicl:.+ (magicl:from-list (list pressure pressure 0) '(3 1))
             (magicl:scale strain viscosity)))

(defun maxwell-linear (strain-increment stress elasticity viscosity dt)
  (magicl:.+ stress (magicl:@ (linear-elastic-matrix elasticity 0d0) strain-increment)))

(defun maxwell (strain-increment stress elasticity viscosity dt vorticity)
  "A stress increment form of a viscoelastic maxwell material"
  (let* ((order 3)
         (strain-matrix (voight-to-matrix strain-increment))
         (pressure (/ (magicl:trace strain-matrix) 3d0))
         (pressure-matrix (magicl:eye 2 :value pressure))
         (viscosity-matrix (magicl:eye 3 :value (/ elasticity viscosity)))
         (dev-stress (magicl:.- strain-matrix pressure-matrix))
         )
    (magicl:.+ stress
               ;;This is the jaumann stress-rate relationship fo 
               (magicl:.- (magicl:@ (linear-elastic-matrix elasticity 0.33d0) strain-increment)
                                        ;(magicl:scale stress (/ (* dt elasticity) viscosity))
                                ;; (magicl:scale (magicl:@ viscosity-matrix (matrix-to-voight dev-stress)) dt)
                                ;; (magicl:scale stress (/ (* dt elasticity) viscosity))
                                (magicl:zeros '(3 1))
                                )
               ;;To translate our jaumann stress increment to cauchy stress increment we need this corrector
               ;; Note this is also a dt adjusted vorticity factor
               (matrix-to-voight
                (magicl::.- (magicl:@ (voight-to-matrix stress) (voight-to-matrix vorticity))
                            (magicl:@ (voight-to-matrix vorticity) (voight-to-matrix stress))
                            ))))
    )
(defun maxwell-exp (strain-increment stress elasticity viscosity dt)
  "A stress increment form of a viscoelastic maxwell material"
  (let* ((order 3)
         ;; (strain-matrix (voight-to-matrix strain-increment))
         ;; (pressure (/ (magicl:trace strain-matrix) 3d0))
         ;; (pressure-matrix (magicl:eye 2 :value pressure))
         ;; (viscosity-matrix (magicl:eye 3 :value (/ elasticity viscosity)))
         ;; (dev-stress (magicl:.- strain-matrix pressure-matrix))
         (rho (/ viscosity elasticity))
         (exp-rho (exp (- (/ dt rho))))
         (lam (* (- 1 exp-rho) (/ rho dt)))
         )
    (magicl:.+ (magicl:scale stress exp-rho)
              (magicl:scale (magicl:@ (linear-elastic-matrix elasticity 0.33d0) strain-increment) lam)
              )
    ;; (magicl:.- (magicl:@ (linear-elastic-matrix elasticity 0.33d0) strain-increment)
    ;;                                     ;(magicl:scale stress (/ (* dt elasticity) viscosity))
    ;;            ;; (magicl:scale (magicl:@ viscosity-matrix (matrix-to-voight dev-stress)) dt)
    ;;            ;; (magicl:scale stress (/ (* dt elasticity) viscosity))
    ;;            (magicl:zeros '(3 1))
    ;;            ))
  ))
