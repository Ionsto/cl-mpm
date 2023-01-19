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
  (magicl:scale!
    (magicl:from-list (list
                        (- 1d0 nu) nu 0d0
                        nu (- 1d0 nu) 0d0
                        0d0 0d0 (- 1d0 (* 2 nu)))
                      '(3 3) :type 'double-float)
    (/ E (* (+ 1 nu) (- 1 nu nu)))))

(defun linear-elastic-mat (strain elastic-matrix)
  "Isotropic linear-elastic constitutive model"
  (magicl:@ elastic-matrix strain))
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
    (magicl:.- (magicl:from-list (list (- pressure) (- pressure) 0d0) '(3 1))
               (magicl:scale dev-strain viscosity))))

(defun maxwell-linear (strain-increment stress elasticity viscosity dt)
  (magicl:.+ stress (magicl:@ (linear-elastic-matrix elasticity 0d0) strain-increment)))


(defun maxwell (strain-increment stress elasticity poisson-ratio de viscosity dt)
  "A stress increment form of a viscoelastic maxwell material"
  (let* ((order 2)
         (stress-matrix (voight-to-matrix stress))
         (pressure (/ (magicl:trace stress-matrix) 2d0))
         (pressure-matrix (magicl:eye 2 :value pressure))
         (dev-stress (magicl:.- stress-matrix pressure-matrix))
         (relaxation-const (/ (* dt elasticity) (* 2d0 (- 1d0 poisson-ratio) viscosity)))
         )
    (declare (type double-float relaxation-const))
    (magicl:.+ stress
               ;; (matrix-to-voight (magicl:eye 2 :value (* (/ elasticity (* 2 (+ 1 poisson-ratio) (- 1 poisson-ratio poisson-ratio))) (magicl:trace (voight-to-matrix strain-increment))) :type 'double-float))
               (magicl:.-
                (magicl:@ de strain-increment)
                (magicl:scale! (matrix-to-voight dev-stress) relaxation-const))
               )))
(declaim (inline voight-eye)
         (ftype (function (double-float) magicl:matrix/double-float) voight-eye))
(defun voight-eye (val)
  ;(let ((arr (make-array 3 :element-type 'double-float :initial-contents (list val val 0d0))))
  ;  (magicl:from-array arr '(3 1) :layout :column-major))
  (magicl:from-array (make-array 3 :element-type 'double-float :initial-contents (list val val 0d0))
                     '(3 1) :layout :column-major)
  )
(declaim (inline voight-trace)
         (ftype (function (magicl:matrix/double-float) double-float) voight-trace))
(defun voight-trace (m)
  (+ (the double-float (magicl:tref m 0 0)) (the double-float (magicl:tref m 1 0))))
(defun maxwell-exp (strain-increment stress elasticity nu viscosity dt)
  "A stress increment form of a viscoelastic maxwell material"
  (let* (
         (stress-matrix (voight-to-matrix stress))
         (pressure (/ (magicl:trace stress-matrix) 2d0))
         (pressure-matrix (magicl:eye 2 :value pressure))

         (dev-stress (magicl:.- stress-matrix pressure-matrix))
         (rho (/ (* 2 (- 1 nu) viscosity) elasticity))
         (exp-rho (exp (- (/ dt rho))))
         (lam (* (- 1 exp-rho) (/ rho dt)))
         (stress-inc (voight-to-matrix (magicl:@ (linear-elastic-matrix elasticity nu) strain-increment)))
         (stress-inc-pressure (magicl:eye 2 :value (/ (magicl:trace stress-inc) 2d0)))
         (stress-inc-dev (magicl:.- stress-inc stress-inc-pressure)
         ))
    (magicl:.+
     (matrix-to-voight (magicl:.+ pressure-matrix stress-inc-pressure))
               (magicl:.+ (magicl:scale! (matrix-to-voight dev-stress) exp-rho)
                          (magicl:scale! (matrix-to-voight stress-inc-dev) lam)))))

(defun maxwell-exp-v (strain-increment stress elasticity nu de viscosity dt &key (result-stress))
  "A stress increment form of a viscoelastic maxwell material"
  (let* (
         (pressure (/ (voight-trace stress) 2d0))
         (pressure-matrix (voight-eye pressure))
         (dev-stress (magicl:.- stress pressure-matrix))
         (rho (/ (* 2d0 (- 1d0 nu) viscosity) elasticity))
         (exp-rho (exp (- (/ dt rho))))
         (lam (* (- 1 exp-rho) (/ rho dt)))
         (stress-inc (magicl:@ (linear-elastic-matrix elasticity nu) strain-increment))
         (stress-inc-pressure (voight-eye (/ (voight-trace stress-inc) 2d0)))
         (stress-inc-dev (magicl:.- stress-inc stress-inc-pressure))
         )
    (declare (double-float rho exp-rho lam viscosity elasticity dt nu))
    ;(declare (dynamic-extent pressure-matrix dev-stress stress-inc stress-inc-pressure stress-inc-dev))
     (magicl:.+
      (magicl:.+ pressure-matrix stress-inc-pressure result-stress)
      (magicl:.+ (magicl:scale! dev-stress exp-rho)
                 (magicl:scale! stress-inc-dev lam) result-stress))
  ))

(defun maxwell-exp-v-simd (strain-increment stress elasticity nu de viscosity dt &key (result-stress))
  "A stress increment form of a viscoelastic maxwell material"
  (let* (
         (pressure (/ (voight-trace stress) 2d0))
         (pressure-matrix (voight-eye pressure))
         (dev-stress (magicl:.- stress pressure-matrix))
         (rho (/ (* 2d0 (- 1d0 nu) viscosity) elasticity))
         (exp-rho (exp (- (/ dt rho))))
         (lam (* (- 1 exp-rho) (/ rho dt)))
         (stress-inc (magicl:@ de strain-increment))
         (stress-inc-pressure (voight-eye (/ (voight-trace stress-inc) 2d0)))
         (stress-inc-dev (magicl:.- stress-inc stress-inc-pressure))
         (result-stress (unless result-stress
                          (magicl:from-array (make-array 3 :element-type 'double-float :initial-element 0d0) '(3 1) :layout :column-major))))
    (declare (double-float rho exp-rho lam viscosity elasticity dt nu))
    (declare (dynamic-extent pressure-matrix dev-stress stress-inc stress-inc-pressure stress-inc-dev))
    (cl-mpm/fastmath:fast-add result-stress pressure-matrix)
    (cl-mpm/fastmath:fast-add result-stress stress-inc-pressure)
    (cl-mpm/fastmath:fast-add result-stress stress-inc-pressure)
    (cl-mpm/fastmath:fast-add result-stress stress-inc-pressure)
    (cl-mpm/fastmath:fast-fmacc result-stress dev-stress exp-rho)
    (cl-mpm/fastmath:fast-fmacc result-stress stress-inc-dev lam)
    result-stress
    ;; (magicl:.+
    ;;  (magicl:.+ pressure-matrix stress-inc-pressure result-stress)
    ;;  (magicl:.+ (magicl:scale! dev-stress exp-rho)
    ;;             (magicl:scale! stress-inc-dev lam) result-stress))
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
                                                              (* 0.5 (- visc-power 1)))))))
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

(defun glen-flow (strain-increment stress bulk-modulus visc-factor visc-power dt vorticity)
  "A stress of a viscoplastic glen flow law material"
  (let* ((order 3)
         ;(strain-matrix (voight-to-matrix strain-increment))
         (pressure (/ (magicl:trace (voight-to-matrix stress)) 2d0))
         (pressure-increment (* bulk-modulus (magicl:trace (voight-to-matrix strain-increment))))
         (pressure-matrix (magicl:eye 2 :value (+ pressure pressure-increment)))
         (dev-stress (glen-stress strain-increment visc-factor visc-power dt)))
    (magicl:.+ (matrix-to-voight pressure-matrix) dev-stress)))

(defun glen-stress (strain visc-factor visc-power dt)
  (let* ((strain-trace (/ (magicl:trace (voight-to-matrix strain)) 2d0))
         (dev-strain (matrix-to-voight (magicl:.- (voight-to-matrix strain) (magicl:eye 2 :value strain-trace))))
         (second-invar (magicl:from-list '(0.5d0 0.5d0 1d0) '(3 1) :type 'double-float))
         (effective-strain (magicl::sum (magicl:.* dev-strain dev-strain second-invar)))
         )
    (if (> effective-strain 0d0)
        (magicl:scale dev-strain (* visc-factor (expt effective-strain
                                                 (* 0.5d0 (- (/ 1d0 visc-power) 1d0)))))
        (magicl:scale dev-strain 0d0))))

(declaim (ftype (function (magicl:matrix/double-float double-float double-float) double-float) glen-viscosity))
(defun glen-viscosity (stress visc-factor visc-power)
  "Get the viscosity for a given stress state"
  (let* ((stress-matrix (voight-to-matrix stress))
         (stress-trace (/ (the double-float (magicl:trace stress-matrix)) 2d0))
         (dev-stress (matrix-to-voight (magicl:.- stress-matrix (magicl:eye 2 :value (the double-float stress-trace)))))
         (second-invar (magicl:from-list '(0.5d0 0.5d0 1d0) '(3 1) :type 'double-float))
         (effective-stress (+ 1d-14 (magicl::sum (magicl:.* dev-stress dev-stress second-invar))))
         )
    (declare (type double-float effective-stress))
    (if (> effective-stress 0d0)
        (/ 1d0 (* 2d0 visc-factor (expt effective-stress (* 0.5d0 (- visc-power 1)))))
        0d0)))

(declaim (ftype (function (magicl:matrix/double-float double-float double-float) double-float) glen-viscosity-strain))
(defun glen-viscosity-strain (strain visc-factor visc-power)
  "Get the viscosity for a given strain state"
  (let* ((strain-matrix (voight-to-matrix strain))
         (strain-trace (/ (the double-float (magicl:trace strain-matrix)) 2d0))
         (dev-strain (matrix-to-voight (magicl:.- strain-matrix (magicl:eye 2 :value (the double-float strain-trace)))))
         (second-invar (magicl:from-list '(0.5d0 0.5d0 1d0) '(3 1) :type 'double-float))
         (effective-strain (+ 1d-14 (magicl::sum (magicl:.* strain strain second-invar))))
         )
    (declare (type double-float effective-strain))
    (if (> effective-strain 0d0)
        (* 0.5d0 visc-factor (expt effective-strain (* 0.5d0 (- (/ 1d0 visc-power) 1d0))))
        0d0)))
