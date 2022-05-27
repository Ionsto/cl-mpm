(defpackage :cl-mpm
  (:use :cl))
(in-package :cl-mpm)

(defun linear-elastic (E nu)
  "Create an isotropic linear elastic matrix"
  (magicl:scale  
    (magicl:from-list (list 
                        (- 1d0 nu) nu 0d0 
                        nu (- 1d0 nu) 0d0 
                        0d0 0d0 (- 1d0 (* 2 nu)))
                      '(3 3) :type 'double-float)
    (/ E (* (+ 1 nu) (- 1 nu nu)))))

(defmethod constitutive-model ((mp particle-elastic) strain)
  (with-slots ((E E)
               (nu nu))
               mp
               (magicl:@ (linear-elastic E nu) strain)))

