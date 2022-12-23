(defpackage :cl-mpm/constitutive
  (:use :cl
        :cl-mpm/utils)
  (:export
    #:linear-elastic
    #:newtonian-fluid
    #:maxwell
    #:maxwell-exp
    #:norton-hoff
    )
  )
(in-package :cl-mpm/constitutive)
(declaim (optimize (debug 0) (safety 0) (speed 3)))

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
  (let* ((strain-matrix (voight-to-matrix strain))
         (dev-strain (matrix-to-voight (magicl:.- strain-matrix (magicl:eye 2 :value (/ (magicl:trace strain-matrix) 2))))))
    (magicl:.+ (magicl:from-list (list (- pressure) (- pressure) 0) '(3 1))
               (magicl:scale dev-strain viscosity))))

(defun maxwell-linear (strain-increment stress elasticity viscosity dt)
  (magicl:.+ stress (magicl:@ (linear-elastic-matrix elasticity 0d0) strain-increment)))

(defun assemble-vorticity-matrix (vorticity)
  (let ((dx (magicl:tref vorticity 0 0))
        (dy (magicl:tref vorticity 1 0))
        (dxdy (magicl:tref vorticity 2 0))
        )
    (magicl:from-list (list dx dxdy (- dxdy) dy) '(2 2))))

(defun maxwell (strain-increment stress elasticity viscosity dt vorticity)
  "A stress increment form of a viscoelastic maxwell material"
  (let* ((order 3)
         (strain-matrix (voight-to-matrix strain-increment))
         (pressure (/ (magicl:trace strain-matrix) 2d0))
         (pressure-matrix (magicl:eye 2 :value pressure))
         (viscosity-matrix (magicl:eye 3 :value (/ elasticity viscosity)))
         (dev-stress (magicl:.- strain-matrix pressure-matrix))
         (relaxation-const (/ (* dt elasticity) viscosity))
         )
    (magicl:.+ stress
               (magicl:.-
                ;;I think this is correct but not sure
                (magicl:@ (linear-elastic-matrix elasticity 0.30d0) strain-increment)
                ;; (magicl:scale (magicl:@ viscosity-matrix (matrix-to-voight dev-stress)) dt)
                ;; (magicl:scale (magicl:from-list (list pressure pressure 0) '(3 1)) elasticity)
                (magicl:from-list (list 0d0 0d0 (* (magicl:tref stress 2 0) relaxation-const)) '(3 1))
                ;; (magicl:scale stress (/ (* dt elasticity) viscosity))
                                        ;(magicl:scale stress (/ (* dt elasticity) viscosity))
                          ;; (magicl:scale stress (/ (* dt elasticity) viscosity))
                          ;; (magicl:zeros '(3 1))
                          ;;To translate our jaumann stress increment to cauchy stress increment we need this corrector
                          ;; Note this is also a dt adjusted vorticity factor
                          (matrix-to-voight
                           (magicl::.- (magicl:@ (voight-to-matrix stress) (assemble-vorticity-matrix vorticity))
                                       (magicl:@ (assemble-vorticity-matrix vorticity) (voight-to-matrix stress))
                                       )))))
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
               (magicl:scale (magicl:@ (linear-elastic-matrix elasticity 0.33d0) strain-increment) lam))
    ;; (magicl:.- (magicl:@ (linear-elastic-matrix elasticity 0.33d0) strain-increment)
    ;;                                     ;(magicl:scale stress (/ (* dt elasticity) viscosity))
    ;;            ;; (magicl:scale (magicl:@ viscosity-matrix (matrix-to-voight dev-stress)) dt)
    ;;            ;; (magicl:scale stress (/ (* dt elasticity) viscosity))
    ;;            (magicl:zeros '(3 1))
    ;;            ))
  ))
(defun norton-hoff (strain-increment stress youngs-modulus poisson-ratio visc-factor visc-power dt vorticity)
  "A stress of a viscoplastic norton-off material"
  (let* ((order 3)
         ;(strain-matrix (voight-to-matrix strain-increment))
         (pressure (/ (magicl:trace (voight-to-matrix stress)) 2d0))
         (pressure-matrix (magicl:eye 2 :value pressure))
         (dev-stress (matrix-to-voight (magicl:.- (voight-to-matrix stress) pressure-matrix)))
         (glenn-strain-rate (magicl:scale dev-stress (* dt
                                                        visc-factor
                                                        (expt (magicl::sum (magicl:.* dev-stress dev-stress
                                                                                        (magicl:from-list
                                                                                         '(0.5d0 0.5d0 1d0) '(3 1))))
                                                              (* 0.5 (- visc-power 1))))))
         )
    (magicl:.+ stress
               (magicl:.-
                ;;I think this is correct but not sure
                (magicl:@ (linear-elastic-matrix youngs-modulus poisson-ratio)
                          (magicl:.- strain-increment glenn-strain-rate))
                (magicl:scale (matrix-to-voight
                               (magicl::.- (magicl:@ (voight-to-matrix stress) (assemble-vorticity-matrix vorticity))
                                           (magicl:@ (assemble-vorticity-matrix vorticity) (voight-to-matrix stress)))
                               ) 0)
                ))
    ))
(defun glen-flow (strain-increment stress bulk-modulus visc-factor visc-power dt vorticity)
  "A stress of a viscoplastic glen flow law material"
  (let* ((order 3)
         ;(strain-matrix (voight-to-matrix strain-increment))
         (pressure (/ (magicl:trace (voight-to-matrix stress)) 2d0))
         (pressure-increment (* bulk-modulus (magicl:trace (voight-to-matrix strain-increment))))
         (pressure-matrix (magicl:eye 2 :value (+ pressure pressure-increment)))
         (dev-stress (glen-stress strain-increment visc-factor visc-power dt)))
    (magicl:.+ (matrix-to-voight pressure-matrix) dev-stress)))
(defun norton-hoff-plastic-strain (stress visc-factor visc-power dt)
  (let* ((order 3)
         ;(strain-matrix (voight-to-matrix strain-increment))
         (pressure (/ (magicl:trace (voight-to-matrix stress)) 2d0))
         (pressure-matrix (magicl:eye 2 :value pressure))
         (dev-stress (matrix-to-voight (magicl:.- (voight-to-matrix stress) pressure-matrix)))
         (glenn-strain-rate (magicl:scale dev-stress (* dt
                                                        visc-factor
                                                        (expt (magicl::sum (magicl:.* dev-stress dev-stress
                                                                                        (magicl:from-list
                                                                                         '(0.5d0 0.5d0 1d0) '(3 1))))
                                                              (- visc-power 1)))))
         )
    glenn-strain-rate
    ))

(defun glen-stress (strain visc-factor visc-power dt)
  (let* ((strain-trace (/ (magicl:trace (voight-to-matrix strain)) 3d0))
         (dev-strain (matrix-to-voight (magicl:.- (voight-to-matrix strain) (magicl:eye 2 :value strain-trace))))
         (second-invar (magicl:from-list '(0.5d0 0.5d0 1d0) '(3 1) :type 'double-float))
         (effective-strain (magicl::sum (magicl:.* dev-strain dev-strain second-invar)))
         )
    (if (> effective-strain 0d0)
        (magicl:scale dev-strain (* visc-factor (expt effective-strain
                                                 (* 0.5 (- (/ 1 visc-power)  1d0)))))
        (magicl:scale dev-strain 0d0))))
