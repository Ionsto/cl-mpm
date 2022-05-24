(ql:quickload "magicl")
(defpackage :cl-mpm
  (:use :cl))
(in-package :cl-mpm)

(defun linear-elastic (E nu)
  "Create an isotropic linear elastic matrix"
  (magicl:scale  
    (magicl:from-list (list 
                        (- 1d0 nu) nu 0d0 
                        nu (- nu) 0d0 
                        0d0 0d0 (- 1d0 (* 2 nu)))
                      '(3 3) :type 'double-float)
    (/ E (* (+ 1 nu) (- 1 nu nu)))))

(defun linear-constitutive-model (mp strain)
  (with-slots ((E youngs-modulus)
               (nu poissons-ratio))
               mp
      (magicl:@ (linear-elastic E nu) strain)))
