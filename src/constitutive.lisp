(defpackage :cl-mpm/constitutive
  (:use :cl
        :cl-mpm/utils)
  (:export
    #:linear-elastic
    #:maxwell
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

(defun maxwell (strain-increment elasticity viscosity)
  "A stress increment form of a viscoelastic maxwell material"
  (let* ((order 2)
         (strain-matrix (voight-to-matrix strain-increment))
         (pressure (/ (magicl:trace strain-matrix) 3d0))
         (pressure-matrix (magicl:eye order :value pressure))
         (dev-stress (magicl:.- strain-matrix pressure-matrix)))
    (matrix-to-voight (magicl:.-
                      (magicl:scale pressure-matrix elasticity)
                      (magicl:scale dev-stress (/ 1d0 viscosity))
                      ))))
