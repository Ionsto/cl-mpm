(defpackage :cl-mpm/constitutive
  (:use :cl
        :cl-mpm/utils)
  (:export
    #:linear-elastic
    #:newtonian-fluid
    #:maxwell
    #:maxwell-exp
    #:maxwell-exp-v
    #:maxwell-exp-v-simd
    #:norton-hoff
    )
  )
(in-package :cl-mpm/constitutive)
;; (declaim (optimize (debug 3) (safety 3) (speed 0)))
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(defun linear-elastic-matrix (E nu)
  "Create an isotropic linear elastic matrix"
  (let ((inv-nu (- 1d0 nu)))
    (magicl:scale!
     (magicl:from-list (list
                        (- 1d0 nu) nu nu 0d0 0d0 0d0
                        nu (- 1d0 nu) nu 0d0 0d0 0d0
                        nu nu (- 1d0 nu) 0d0 0d0 0d0
                        0d0 0d0 0d0 (* 0.5d0 (- 1d0 (* 2d0 nu))) 0d0 0d0
                        0d0 0d0 0d0 0d0 (* 0.5d0 (- 1d0 (* 2d0 nu))) 0d0
                        0d0 0d0 0d0 0d0 0d0 (* 0.5d0 (- 1d0 (* 2d0 nu))))
                       '(6 6) :type 'double-float)
     (/ E (* (+ 1d0 nu) (- 1d0 (* 2d0 nu)))))))

(defun linear-elastic-matrix-ps (E nu)
  "Create an isotropic linear elastic matrix"
  (let ((inv-nu (/ 1d0 (- 1d0 (expt nu 2)))))
    (magicl:scale!
     (magicl:from-list (list
                        1d0 nu 0d0 0d0 0d0 0d0
                        nu 1d0 0d0 0d0 0d0 0d0
                        0d0 0d0 0d0 0d0 0d0 0d0
                        0d0 0d0 0d0 0d0 0d0 0d0
                        0d0 0d0 0d0 0d0 0d0 0d0
                        0d0 0d0 0d0 0d0 0d0 (* 0.5d0 (- 1d0 nu)))
                       '(6 6) :type 'double-float)
     (/ E (- 1d0 (* nu nu))))))

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
         (dev-strain (matrix-to-voight (magicl:.- strain-matrix (magicl:eye 2 :value (/ (magicl:trace strain-matrix) 3))))))
    (magicl:.- (magicl:from-list (list (- pressure) (- pressure) 0d0) '(3 1))
               (magicl:scale dev-strain viscosity))))

(defun maxwell-linear (strain-increment stress elasticity viscosity dt)
  (magicl.simd::.+-simd
   stress
   (magicl:@
    (linear-elastic-matrix elasticity 0d0) strain-increment)))


(defun maxwell (strain-increment stress elasticity poisson-ratio de viscosity dt)
  "A stress increment form of a viscoelastic maxwell material"
  (let* ((relaxation-const (/ (* dt elasticity) (* 2d0 (- 1d0 poisson-ratio) viscosity))))
    (declare (type double-float relaxation-const))
    (magicl:.-
     (magicl:@ de strain-increment)
     (magicl:scale! (deviatoric-voigt stress) relaxation-const))
    ))

(defun maxwell-exp-inc (strain-increment stress elasticity poisson-ratio de viscosity dt)
  "A stress increment form of a viscoelastic maxwell material"
  (magicl:.-
   stress
   (maxwell-exp-v strain-increment stress elasticity poisson-ratio de viscosity dt)))
(defun maxwell-damage (strain-increment stress elasticity poisson-ratio de viscosity dt damage strain)
  "A stress increment form of a viscoelastic maxwell material"
  (let* ((relaxation-const (/ (* dt elasticity) (* 2d0 (- 1d0 poisson-ratio) viscosity)))
         (elastic-increment (magicl:@ de strain-increment))
         )
    (declare (type double-float relaxation-const))
    (let* ((stress-matrix (voight-to-matrix elastic-increment))
           (p (/ (magicl:trace stress-matrix) 3d0))
           (pressure-matrix (magicl:eye 2 :value p))
           (dev-stress (magicl:.- stress-matrix pressure-matrix)))
      (setf elastic-increment (matrix-to-voight
                    (magicl.simd::.+-simd pressure-matrix
                               (magicl:scale! dev-stress (max 0d0 (- 1d0 damage)))))))
    ;; (magicl:scale! elastic-increment (max 1d-3 (- 1d0 damage)))
    (magicl:.-
     elastic-increment
     ;; (magicl:scale! (magicl:@ de strain-increment) (max 1d-3 (expt (- 1d0 damage) 1d0)))
     ;; (magicl:@ (apply-damage-constitutive de strain damage) strain-increment)
     (magicl:scale! (deviatoric-voigt stress) relaxation-const))
    ))
;; (defun maxwell-damage (strain-increment stress elasticity poisson-ratio de viscosity dt damage)
;;   "A stress increment form of a viscoelastic maxwell material"
;;   (let* ((order 2)
;;          (stress-matrix (voight-to-matrix stress))
;;          (pressure (/ (magicl:trace stress-matrix) 3d0))
;;          (pressure-matrix (magicl:eye 2 :value pressure))
;;          (dev-stress (magicl:.- stress-matrix pressure-matrix))
;;          (relaxation-const (/ (* dt elasticity)
;;                               (* 2d0 (- 1d0 poisson-ratio) viscosity)))
;;          )
;;     (declare (type double-float relaxation-const))
;;     (magicl:.-
;;      (magicl:scale! (magicl:@ de strain-increment) (max 1d-3 (- 1d0 damage)))
;;                                         ;(magicl:@ (linear-elastic-matrix (* (- 1d0 damage) elasticity) poisson-ratio) strain-increment)
;;      (magicl:scale! (matrix-to-voight dev-stress) relaxation-const))
;;     ;; (matrix-to-voight (tensile-project (voight-to-matrix (magicl:@ de strain-increment)) stress damage))
;;     ))

(defun tensile-projection-eig (stress)
  (flet ((H (x) (if (> x 0d0) 1d0 0d0)))
    (multiple-value-bind (l v) (cl-mpm/utils::eig (voight-to-matrix stress))
      (loop for i from 0 to 1
            do (let* ((sii (nth i l)))
                 (setf (nth i l) (* sii (H sii)))))
      (matrix-to-voight (magicl:@ v
                                  (magicl:from-diag l :type 'double-float)
                                  (magicl:transpose v))))))

(defun tensile-projection-Q-mandel (stress)
  (flet ((H (x) (if (> x 0d0) 1d0 0d0)))
    (let ((Q (magicl:zeros '(3 3) :type 'double-float)))
      (multiple-value-bind (l v) (cl-mpm/utils::eig (voight-to-matrix stress))
        (loop for i from 0 to 1
              do (let* ((sii (nth i l))
                        (vii (magicl::column v i))

                        (vsi (magicl:@ vii (magicl:transpose vii)))
                        (comp-prod
                          (magicl:scale (magicl:@
                                         (matrix-to-mandel vsi)
                                         (magicl:transpose (matrix-to-mandel vsi)))
                                        (H sii)))

                        (v1 (magicl:tref vii 0 0))
                        (v2 (magicl:tref vii 1 0))

                        (v1111 (* v1 v1 v1 v1))
                        (v2222 (* v2 v2 v2 v2))
                        (v1122 (* v1 v1 v2 v2))
                        (v1112 (* v1 v1 v1 v2))
                        (v1222 (* v1 v2 v2 v2))
                        (v1212 (* v1 v2 v1 v2))

                        (comp
                          (magicl:scale!
                           (magicl:from-list (list (* 1 v1111) (* 1 v1122) (* (sqrt 2) v1112)
                                                   (* 1 v1122) (* 1 v2222) (* (sqrt 2) v1222)
                                                   (* (sqrt 2) v1112) (* (sqrt 2) v1222) (* 2 v1212)) '(3 3))
                           (H sii))
                          )
                        )
                   (magicl.simd::.+-simd Q comp Q)
                   )
              )
      Q))))

(defun tensile-projection-Q-cw-mandel (strain)
  (flet ((H (x) (if (> x 0d0) 1d0 0d0)))
    (let ((Q (magicl:zeros '(3 3) :type 'double-float)))
      (multiple-value-bind (l v) (cl-mpm/utils::eig (voight-to-matrix (magicl.simd::.*-simd strain
                                                                          (magicl:from-list '(1d0 1d0 0.5d0) '(3 1))
                                                                          )))
        (loop for i from 0 to 1
              do (let* ((si (nth i l))
                        (vi (magicl::column v i))
                        (vii (magicl:@ vi (magicl:transpose vi)))
                        (comp-prod
                          (magicl:scale (magicl:@
                                         (matrix-to-mandel vii)
                                         (magicl:transpose (matrix-to-mandel vii)))
                                        (H si))))
                   (magicl.simd::.+-simd Q comp-prod Q)
                   )
              )
        (let* ((si (nth 0 l))
               (sj (nth 1 l))
               (vi (magicl::column v 0))
               (vj (magicl::column v 1))
               (vij (magicl:@ vi (magicl:transpose vj)))
               (vji (magicl:@ vj (magicl:transpose vi)))
               (pij (magicl:scale!
                     (magicl.simd::.+-simd vij vji)
                     0.5d0))
               (comp-prod
                 (magicl:scale! (magicl:@
                                (matrix-to-mandel pij)
                                (magicl:transpose (matrix-to-mandel pij)))
                               (+ (H si) (H sj)))))
          (magicl.simd::.+-simd Q comp-prod Q))
          )
      Q)))

(defun tensile-project-q (stress)
  (mandel-to-voigt
   (magicl:@
    (tensile-projection-q-mandel stress)
    (voigt-to-mandel stress))))

(defun tensile-projection-A (strain damage)
  "Generate a mandel form tensile projection matrix A* from stress"
  (let ((Q (tensile-projection-q-cw-mandel strain))
        (I (magicl:from-diag '(1d0 1d0 1d0))))
    (magicl.simd::.+-simd I (magicl:scale Q (- (sqrt (- 1d0 damage)) 1d0)))))

;; (defun test-tensile ()
;;   (loop for stress in (list
;;                        (magicl:from-list '(1 0 0) '(3 1) :type 'double-float)
;;                        (magicl:from-list '(0 1 0) '(3 1) :type 'double-float)
;;                        (magicl:from-list '(1 1 0) '(3 1) :type 'double-float)
;;                        (magicl:from-list '(-1 0 0) '(3 1) :type 'double-float)
;;                        (magicl:from-list '(0 -1 0) '(3 1) :type 'double-float)
;;                        (magicl:from-list '(-1 -1 0) '(3 1) :type 'double-float)
;;                        (magicl:from-list '(1 -1 0) '(3 1) :type 'double-float)
;;                        (magicl:from-list '(-1 1 0) '(3 1) :type 'double-float)
;;                        (magicl:from-list '(0 0 1) '(3 1) :type 'double-float)
;;                        (magicl:from-list '(1 0 1) '(3 1) :type 'double-float)
;;                        (magicl:from-list '(0 1 1) '(3 1) :type 'double-float)
;;                        (magicl:from-list '(-1 0 1) '(3 1) :type 'double-float)
;;                        (magicl:from-list '(0 -1 1) '(3 1) :type 'double-float)
;;                         )
;;         do (let* ((stress-eig (tensile-projection-eig stress 1d0))
;;                   (stress-project (tensile-project-q stress))
;;                  )
;;              (format t "Stress~%")
;;              (loop for i from 0 to 2
;;                    do (format t "~A ~A ~A ~%"
;;                               (magicl:tref stress i 0)
;;                               (magicl:tref stress-eig i 0)
;;                               (magicl:tref stress-project i 0)
;;                               )))
;;     ))

(defun tensile-project (stiffness stress damage)
  (if (> damage 0.0d0)
    (let ((damaged-stiffness (magicl:zeros '(2 2) :type 'double-float)))
      (multiple-value-bind (l v) (cl-mpm/utils::eig (voight-to-matrix stress))
        (loop for i from 0 to 1
              do (let* ((sii (nth i l))
                        (vii (magicl::column v i))
                        (scale 1d0)
                        (A (magicl:@ vii (magicl:transpose vii)))
                        )
                   (when (> sii 0d0)
                     (setf scale (sqrt (- 1d0 damage)))
                     (magicl:scale! A scale))
                   (magicl.simd::.+-simd damaged-stiffness
                              (magicl:@ A stiffness A)
                              damaged-stiffness)
                     )))
      damaged-stiffness)
    stiffness))

(let ((mandel-constant
        (magicl:from-list (list 1d0 1d0 (sqrt 2)
                                1d0 1d0 (sqrt 2)
                                (sqrt 2) (sqrt 2) 2d0) '(3 3))))
  (defun voigt-4th-to-mandel (tensor)
    (magicl.simd::.*-simd tensor mandel-constant))
  (defun mandel-to-voigt-4th (tensor)
    (magicl:./ tensor mandel-constant))
  )
(defun apply-damage-constitutive (de strain damage)
  (let ((A (tensile-projection-a strain damage)))
    (mandel-to-voigt-4th (magicl:@ A (voigt-4th-to-mandel de) A))))


(defun elasto-glen (strain-increment stress E nu de viscosity dt strain)
  "A absolute stress form of a nonlinear glen flow with elastic volume"
  (let* ((order 2)
         (strain-dev (deviatoric strain-increment))
         (bulk-modulus (/ E (* (+ 1 nu) (- 1 nu nu))))
         ;; (pressure-inc (* bulk-modulus (voight-trace strain-increment)))
         ;; (pressure (+ pressure-inc (* 0.5 (voight-trace stress))))
         (pressure (* bulk-modulus 0.5 (voight-trace strain)))
         )
    (declare (type double-float pressure))
    (magicl.simd::.+-simd
     (magicl:from-list (list pressure pressure 0d0) '(3 1))
     (magicl:scale! strain-dev (/ (* 2d0 viscosity) dt))
    )))

(defun elasto-glen-damage (strain-increment stress E nu de viscosity dt strain damage)
  "A absolute stress form of a nonlinear glen flow with elastic volume"
  (let* ((order 2)
         (strain-dev (deviatoric strain-increment))
         (bulk-modulus (/ E (* (+ 1 nu) (- 1 nu nu))))
         (pressure (* bulk-modulus 0.5 (voight-trace strain)))
         )
    (declare (type double-float pressure))
    (when (< pressure 0)
      (setf pressure (* pressure (- 1d0 damage))))
    (magicl.simd::.+-simd
     (magicl:from-list (list pressure pressure 0d0) '(3 1))
     (magicl:scale! strain-dev (/ (* 2d0 viscosity) dt))
     )))

(declaim (inline voight-eye)
         (ftype (function (double-float) magicl:matrix/double-float) voight-eye))
(defun voight-eye (val)
  ;(let ((arr (make-array 3 :element-type 'double-float :initial-contents (list val val 0d0))))
  ;  (magicl:from-array arr '(3 1) :layout :column-major))
  ;; (magicl:from-array (make-array 3 :element-type 'double-float :initial-contents (list val val 0d0))
  ;;                    '(3 1) :layout :column-major)
  (cl-mpm/utils:voigt-from-list (list val val val 0d0 0d0 0d0)))
(declaim (inline voight-trace)
         (ftype (function (magicl:matrix/double-float) double-float) voight-trace))
(defun voight-trace (m)
  "Take the trace of a voigt notation stress vector"
  (let ((s (magicl::matrix/double-float-storage m)))
    (+ (the double-float (aref s 0))
       (the double-float (aref s 1))
       (the double-float (aref s 2)))))

(declaim (ftype (function (magicl:matrix/double-float
                           magicl:matrix/double-float
                           double-float
                           double-float
                           magicl:matrix/double-float
                           double-float
                           double-float
                           )) maxwell-exp))
(defun maxwell-exp (strain-increment stress elasticity nu de viscosity dt)
  ;; (declare (optimize (safety 3) (speed 0)))
  "A stress increment form of a viscoelastic maxwell material"
  (let* (
         (stress-matrix (voight-to-matrix stress))
         (pressure (/ (magicl:trace stress-matrix) 3d0))
         (pressure-matrix (magicl:eye 3 :value pressure))
         (dev-stress (magicl:.- stress-matrix pressure-matrix))
         (rho (/ (* 2 (- 1 nu) viscosity) elasticity))
         (exp-rho (exp (- (/ dt rho))))
         (lam (* (- 1 exp-rho) (/ rho dt)))
         ;; (exp-rho 1d0)
         ;; (lam 1d0)
         (stress-inc (voight-to-matrix (magicl:@ de strain-increment)))
         (stress-inc-pressure (magicl:eye 3 :value (/ (magicl:trace stress-inc) 3d0)))
         (stress-inc-dev (magicl:.- stress-inc stress-inc-pressure)
         ))
    (magicl.simd::.+-simd
     (matrix-to-voight (magicl.simd::.+-simd pressure-matrix stress-inc-pressure))
               (magicl.simd::.+-simd
                (magicl:scale! (matrix-to-voight dev-stress) exp-rho)
                (magicl:scale! (matrix-to-voight stress-inc-dev) lam)))))

(defun maxwell-exp-v (strain-increment stress elasticity nu de viscosity dt &key (result-stress))
  "A stress increment form of a viscoelastic maxwell material"
  (let* (
         (pressure (/ (voight-trace stress) 3d0))
         (pressure-matrix (voight-eye pressure))
         (dev-stress (magicl:.- stress pressure-matrix))
         (rho (/ (* 2d0 (- 1d0 nu) viscosity) elasticity))
         (exp-rho (exp (- (/ dt rho))))
         (lam (* (- 1 exp-rho) (/ rho dt)))
         (stress-inc (magicl:@ (linear-elastic-matrix elasticity nu) strain-increment))
         (stress-inc-pressure (voight-eye (/ (voight-trace stress-inc) 3d0)))
         (stress-inc-dev (magicl:.- stress-inc stress-inc-pressure))
         )
    (declare (double-float rho exp-rho lam viscosity elasticity dt nu))
    ;(declare (dynamic-extent pressure-matrix dev-stress stress-inc stress-inc-pressure stress-inc-dev))
     (magicl.simd::.+-simd
      (magicl.simd::.+-simd pressure-matrix stress-inc-pressure)
      (magicl.simd::.+-simd (magicl:scale! dev-stress exp-rho)
                 (magicl:scale! stress-inc-dev lam)))))

(declaim (ftype (function (magicl:matrix/double-float
                           magicl:matrix/double-float
                           double-float
                           double-float
                           magicl:matrix/double-float
                           double-float
                           double-float
                           &optional
                           (magicl:matrix/double-float nil)
                           )) maxwell-exp-v-simd))

(defun maxwell-exp-v-simd (strain-increment stress elasticity nu de viscosity dt &optional (result-stress nil))
  "A stress increment form of a viscoelastic maxwell material"
  (let* (
         (pressure (/ (voight-trace stress) 3d0))
         (pressure-matrix (voight-eye pressure))
         (dev-stress (magicl:.- stress pressure-matrix))
         (rho (/ (* 2d0 (- 1d0 nu) viscosity) elasticity))
         (exp-rho (exp (- (/ dt rho))))
         (lam (* (- 1 exp-rho) (/ rho dt)))
         (stress-inc (magicl:@ de strain-increment))
         (stress-inc-pressure (voight-eye (/ (voight-trace stress-inc) 3d0)))
         (stress-inc-dev (magicl:.- stress-inc stress-inc-pressure))
         (result-stress (if result-stress
                            result-stress
                            (magicl:zeros '(3 1)))))
    (declare (double-float rho exp-rho lam viscosity elasticity dt nu))
    ;; (declare (dynamic-extent pressure-matrix dev-stress stress-inc stress-inc-pressure stress-inc-dev))
    (cl-mpm/fastmath:fast-add result-stress pressure-matrix)
    (cl-mpm/fastmath:fast-add result-stress stress-inc-pressure)
    (cl-mpm/fastmath:fast-add result-stress stress-inc-pressure)
    (cl-mpm/fastmath:fast-add result-stress stress-inc-pressure)
    (cl-mpm/fastmath:fast-fmacc result-stress dev-stress exp-rho)
    (cl-mpm/fastmath:fast-fmacc result-stress stress-inc-dev lam)
    result-stress
    ;; (magicl.simd::.+-simd
    ;;  (magicl.simd::.+-simd pressure-matrix stress-inc-pressure result-stress)
    ;;  (magicl.simd::.+-simd (magicl:scale! dev-stress exp-rho)
    ;;             (magicl:scale! stress-inc-dev lam) result-stress))
  ))

(defun norton-hoff (strain-increment stress youngs-modulus poisson-ratio visc-factor visc-power dt vorticity)
  "A stress of a viscoplastic norton-off material"
  (let* ((order 3)
         ;(strain-matrix (voight-to-matrix strain-increment))
         (pressure (/ (magicl:trace (voight-to-matrix stress)) 3d0))
         (pressure-matrix (magicl:eye 2 :value pressure))
         (dev-stress (matrix-to-voight (magicl:.- (voight-to-matrix stress) pressure-matrix)))
         (glenn-strain-rate (magicl:scale dev-stress (* dt
                                                        visc-factor
                                                        (expt (magicl::sum (magicl.simd::.*-simd dev-stress dev-stress
                                                                                        (magicl:from-list
                                                                                         '(0.5d0 0.5d0 1d0) '(3 1))))
                                                              (* 0.5 (- visc-power 1)))))))
    (magicl.simd::.+-simd stress
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
         (pressure (/ (magicl:trace (voight-to-matrix stress)) 3d0))
         (pressure-matrix (magicl:eye 2 :value pressure))
         (dev-stress (matrix-to-voight (magicl:.- (voight-to-matrix stress) pressure-matrix)))
         (glenn-strain-rate (magicl:scale dev-stress (* dt
                                                        visc-factor
                                                        (expt (magicl::sum (magicl.simd::.*-simd dev-stress dev-stress
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
         (pressure (/ (magicl:trace (voight-to-matrix stress)) 3d0))
         (pressure-increment (* bulk-modulus (magicl:trace (voight-to-matrix strain-increment))))
         (pressure-matrix (magicl:eye 2 :value (+ pressure pressure-increment)))
         (dev-stress (glen-stress strain-increment visc-factor visc-power dt)))
    (magicl.simd::.+-simd (matrix-to-voight pressure-matrix) dev-stress)))

(defun glen-stress (strain visc-factor visc-power dt)
  (let* ((strain-trace (/ (magicl:trace (voight-to-matrix strain)) 3d0))
         (dev-strain (matrix-to-voight (magicl:.- (voight-to-matrix strain) (magicl:eye 2 :value strain-trace))))
         (second-invar (magicl:from-list '(0.5d0 0.5d0 1d0) '(3 1) :type 'double-float))
         (effective-strain (magicl::sum (magicl.simd::.*-simd dev-strain dev-strain second-invar)))
         )
    (if (> effective-strain 0d0)
        (magicl:scale dev-strain (* visc-factor (expt effective-strain
                                                 (* 0.5d0 (- (/ 1d0 visc-power) 1d0)))))
        (magicl:scale dev-strain 0d0))))

(defun deviatoric (voigt)
  (let* ((mat (voight-to-matrix voigt))
         (trace (/ (magicl:trace mat) 3d0)))
    (matrix-to-voight (magicl:.- mat (magicl:eye 3 :value (the double-float trace))))))

(defun deviatoric-mat (mat)
  (let* ((trace (/ (magicl:trace mat) 3d0)))
    (matrix-to-voight (magicl:.- mat (magicl:eye 3 :value (the double-float trace))))))


(declaim (ftype (function (magicl:matrix/double-float double-float double-float) double-float) glen-viscosity))
(defun glen-viscosity (stress visc-factor visc-power)
  "Get the viscosity for a given stress state"
  (let* ((dev-stress (deviatoric stress))
         (second-invar (magicl:from-list '(0.5d0 0.5d0 1d0) '(3 1) :type 'double-float))
         (effective-stress (+ 1d-30 (magicl::sum (magicl.simd::.*-simd dev-stress dev-stress second-invar))))
         )
    (declare (type double-float effective-stress))
    (if (> effective-stress 0d0)
        (/ 1d0 (* 2d0 visc-factor (expt effective-stress (* 0.5d0 (- visc-power 1)))))
        0d0)))
(defun effective-strain-rate (strain)
  (sqrt (* 0.5d0 (cl-mpm/fastmath::voigt-tensor-reduce-simd (deviatoric strain)))))

(declaim (ftype (function (magicl:matrix/double-float double-float double-float) double-float) glen-viscosity-strain))
(defun glen-viscosity-strain (strain visc-factor visc-power)
  "Get the viscosity for a given strain state"
  (let* ((effective-strain (+ 1d-20 (effective-strain-rate strain))))
    (declare (type double-float effective-strain))
    (* visc-factor (the double-float (expt effective-strain (* 1d0 (- (/ 1d0 visc-power) 1d0)))))))

(defun glen-viscosity-stress (stress visc-factor visc-power)
  (declare (double-float visc-factor visc-power))
  "Get the viscosity for a given strain state"
  (let* ((A (expt visc-factor (- visc-power)))
         (effective-strain (sqrt (* 0.5d0 (expt (cl-mpm/utils::trace-voigt (deviatoric stress)) 2d0)))))
    (declare (type double-float effective-strain A))
    (if (> effective-strain 0d0)
        (/ 1d0 (* 2d0 A (the double-float (expt effective-strain (- visc-power 1d0)))))
        0d0)
    ))


(declaim (ftype (function (magicl:matrix/double-float double-float double-float) (values double-float)) glen-viscosity-strain))
(defun glen-viscosity-stretch (stretch visc-factor visc-power)
  "Get the viscosity for a given strain state"
  (let* (;(effective-strain (+ 1d-15 (cl-mpm/fastmath::voigt-tensor-reduce-simd (deviatoric strain))))
         (dev-stretch (deviatoric-voigt stretch))
         (effective-strain (+ 1d-15 (* 0.5d0 (magicl::sum (magicl.simd::.*-simd dev-stretch dev-stretch)))))
         )
    (declare (type double-float effective-strain))
    (if (> effective-strain 0d0)
        (values (* 0.5d0 visc-factor (expt effective-strain (* 0.5d0 (- (/ 1d0 visc-power) 1d0)))))
        (values 0d0))))

(defun voigt-j2 (s)
  "Calculate j2 invarient from deviatoric stress"
  (let ((storage (magicl::matrix/double-float-storage s)))
    (/ (+ (the double-float (cl-mpm/fastmath:dot s s))
          (the double-float (expt (aref storage 3) 2))
          (the double-float (expt (aref storage 4) 2))
          (the double-float (expt (aref storage 5) 2))
          ) 2d0)))

(defun vm-yield-func (j2 rho)
  (declare (double-float j2 rho))
  (- (/ (the double-float (sqrt (* 2d0 j2))) rho) 1d0))

(defun vm-derivatives (sig rho)
  (declare (double-float rho))
  (let* ((s (deviatoric-voigt sig))
         (j2 (the double-float (voigt-j2 s)))
         (dj2 (magicl:.* s (cl-mpm/utils::voigt-from-list '(1d0 1d0 1d0 2d0 2d0 2d0))))
         (ddj2 (magicl:block-matrix (list (magicl:.- (magicl:eye 3) (magicl:const (/ 1d0 3d0) '(3 3)))
                                          (magicl:zeros '(3 3))
                                          (magicl:zeros '(3 3))
                                          (magicl:eye 3 :value 2d0)) '(2 2)))
         (df (magicl:scale dj2 (the double-float (/ 1d0 (the double-float (* rho (the double-float (sqrt (the double-float (* 2d0 j2))))))))))
         (ddf (magicl:scale!
               (magicl:.-
                (magicl:scale! ddj2 (/ 1d0 (the double-float (sqrt (* 2d0 j2)))))
                (magicl:scale! (magicl:@ dj2 (magicl:transpose dj2)) (/ 1d0 (the double-float (expt (* 2d0 j2) 3/2))))
                )
               (/ 1d0 rho))))
    (values df ddf)))

(defun b-norm (b)
  (sqrt (loop for i from 0 below 6
             sum (expt (magicl:tref b i 0) 2))))
(defun vm-plastic (stress de trial-elastic-strain rho)
  (let* ((tol 1d-9)
        (max-iter 5)
        ;; (sig stress)
        (sig (cl-mpm/utils::voigt-copy stress))
        (s (deviatoric-voigt stress))
        (j2 (voigt-j2 s))
        (f (vm-yield-func j2 rho))
        )
    (declare (dynamic-extent s ))
    (if (> f tol)
      (let ((eps-e (cl-mpm/utils::voigt-copy trial-elastic-strain))
            (df (cl-mpm/utils:voigt-zeros))
            (ddf (cl-mpm/utils:matrix-zeros)))
        (declare (dynamic-extent df ddf))
        (multiple-value-bind (ndf nddf) (vm-derivatives sig rho)
          (setf df ndf
                ddf nddf))
        ;;Plastic deformation is occuring
        (let* (
               (b (magicl:from-list (list 0d0 0d0 0d0 0d0 0d0 0d0 f) '(7 1) :type 'double-float))
               (dgam 0d0))
          (loop for iters from 0 to max-iter
                when (or (> (b-norm b) tol)
                         (> (abs (magicl:tref b 6 0)) tol))
                  do
                     (progn
                       ;; (format t "it:~D f:~F~%" iters (magicl:tref b 6 0))
                       ;;Calculate backstress
                       (let* ((A (magicl:block-matrix
                                  (list
                                   (magicl:.+ (magicl:eye 6)
                                              (magicl:scale! (magicl:@ ddf de) dgam))
                                   df
                                   (magicl:@ (magicl:transpose df) de)
                                   (magicl:zeros '(1 1))
                                   )
                                  '(2 2)
                                  ))
                              (dx (magicl:scale! (magicl:@ (magicl:inv A) b) -1d0)))
                         ;;Add just the 6 components of stress
                         (loop for i from 0 below 6
                               do (incf (magicl:tref eps-e i 0) (magicl:tref dx i 0)))
                         (incf dgam (magicl:tref dx 6 0))
                         (setf sig (magicl:@ de eps-e))
                         (setf s (deviatoric-voigt sig))
                         (setf j2 (voigt-j2 s))
                         (setf f (vm-yield-func j2 rho))
                         (multiple-value-bind (ndf nddf) (vm-derivatives sig rho)
                           (setf df ndf
                                 ddf nddf))
                         (let ((b-eps (magicl:.+ eps-e
                                                 (magicl:scale trial-elastic-strain -1d0)
                                                 (magicl:scale! df dgam)
                                                 )
                                      )
                               (b-f (vm-yield-func j2 rho)))
                           (loop for i from 0 below 6
                                 do (setf (magicl:tref b i 0) (magicl:tref b-eps i 0)))
                           (setf (magicl:tref b 6 0) b-f)))
                       ))
          ;; (when (or (> (b-norm b) tol)
          ;;           (> (abs (magicl:tref b 6 0)) tol))
          ;;   (format t "Bad VM solve~%"))
          (values sig eps-e f)
          ))
      (values sig trial-elastic-strain f))))



(defun mc-yield-func (stress phi c)
  (let* ((k (/ (+ 1d0 (sin phi)) (- 1d0 (sin phi))))
         (sigc (* 2d0 c (sqrt k)))
         (f (- (* k (magicl:tref stress 0 0)) (magicl:tref stress 2 0) sigc)))
    f))

(defun rotate-vector (vec)
  (vector-from-list (list (magicl:tref vec 1 0)
                          (magicl:tref vec 2 0)
                          (magicl:tref vec 0 0))))
(defun swizzle-voigt->coombs (vec)
  (voigt-from-list (list
                    (magicl:tref vec 0 0)
                    (magicl:tref vec 1 0)
                    (magicl:tref vec 2 0)
                    (magicl:tref vec 4 0)
                    (magicl:tref vec 5 0)
                    (magicl:tref vec 3 0)
                    )))
(defun swizzle-coombs->voigt (vec)
  (voigt-from-list (list
                    (magicl:tref vec 0 0)
                    (magicl:tref vec 1 0)
                    (magicl:tref vec 2 0)
                    (magicl:tref vec 4 0)
                    (magicl:tref vec 5 0)
                    (magicl:tref vec 3 0)
                    )))

(declaim (notinline  mc-plastic))
(defun mc-plastic (stress de trial-elastic-strain E nu phi psi c)
  (declare (optimize (speed 0) (safety 3) (debug 3)))
  (let* ((tol 1d-9)
         (sig (cl-mpm/utils::voigt-copy (swizzle-voigt->coombs stress)))
         (eps-e (cl-mpm/utils:vector-zeros))
        )
    (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix (swizzle-voigt->coombs trial-elastic-strain)))
      (let* ((l-sort (sort (mapcar #'cons l (list 0 1 2)) #'> :key #'car))
             (l (mapcar #'car l-sort))
             (v (magicl:block-matrix (list
                                      (magicl:column v (cdr (nth 0 l-sort)))
                                      (magicl:column v (cdr (nth 1 l-sort)))
                                      (magicl:column v (cdr (nth 2 l-sort)))) '(1 3))))
        (let* ((De3
                 (magicl:scale!
                  (magicl:from-list (list
                                     (- 1d0 nu) nu nu
                                     nu (- 1d0 nu) nu
                                     nu nu (- 1d0 nu))
                                    '(3 3) :type 'double-float)
                  (/ E (* (+ 1d0 nu) (- 1d0 (* 2d0 nu))))))
               (epsTr (cl-mpm/utils:vector-from-list l))
               (Ce (magicl:inv De3))
               (sig (magicl:@ De3 epsTr))
               (eps-e (cl-mpm/utils::vector-copy epsTr))

               (k (/ (+ 1 (sin phi)) (- 1d0 (sin phi))))
               (sigc (* 2d0 c (sqrt k)))
               (f (mc-yield-func sig phi c)))

          (if (> f tol)
              (let* ((m (/ (+ 1 (sin psi)) (- 1d0 (sin psi))))
                     (siga (magicl:scale! (cl-mpm/utils:vector-from-list (list 1d0 1d0 1d0))
                                          (/ sigc (- k 1d0))
                                          ))
                     (r1 (vector-from-list (list 1d0 1d0 k)))
                     (r2 (vector-from-list (list 1d0 k k)))
                     (rg1 (vector-from-list (list 1d0 1d0 m)))
                     (rg2 (vector-from-list (list 1d0 m m)))
                     (df (vector-from-list (list k 0d0 -1d0)))
                     (dg (vector-from-list (list m 0d0 -1d0)))
                     (rp (magicl:scale! (magicl:@ De3 dg)
                                        (/ 1d0 (magicl:tref (magicl:@ (magicl:transpose dg) De3 df) 0 0))))
                     (t1 (/ (magicl:tref (magicl:@ (magicl:transpose rg1) Ce (magicl:.- sig siga)) 0 0)
                            (magicl:tref (magicl:@ (magicl:transpose rg1) Ce r1) 0 0)))
                     (t2 (/ (magicl:tref (magicl:@ (magicl:transpose rg2) Ce (magicl:.- sig siga)) 0 0)
                            (magicl:tref (magicl:@ (magicl:transpose rg2) Ce r2) 0 0)))
                     (f12 (magicl:tref (magicl:@ (magicl:transpose! (cl-mpm/fastmath::cross-product rp r1))
                                                 (magicl:.- sig siga)) 0 0))
                     (f13 (magicl:tref (magicl:@ (magicl:transpose! (cl-mpm/fastmath::cross-product rp r2))
                                                 (magicl:.- sig siga)) 0 0))
                     (path :no-return)
                     (Q
                       (magicl:transpose!
                        (magicl:block-matrix (list
                                              (magicl:.* (magicl:column v 0)
                                                         (magicl:column v 0))

                                              (magicl:.* (magicl:column v 1)
                                                         (magicl:column v 1))
                                              (magicl:.* (magicl:column v 2)
                                                         (magicl:column v 2))

                                              (magicl:scale!
                                               (magicl:.* (magicl:column v 0)
                                                          (magicl:column v 1)) 2d0)
                                              (magicl:scale!
                                               (magicl:.* (magicl:column v 1)
                                                          (magicl:column v 2)) 2d0)
                                              (magicl:scale!
                                               (magicl:.* (magicl:column v 2)
                                                          (magicl:column v 0)) 2d0)

                                              (magicl:.* (magicl:column v 0)
                                                         (rotate-vector (magicl:column v 0)))
                                              (magicl:.* (magicl:column v 1)
                                                         (rotate-vector (magicl:column v 1)))
                                              (magicl:.* (magicl:column v 2)
                                                         (rotate-vector (magicl:column v 2)))

                                              (magicl:scale! (magicl:.+
                                                             (magicl:.* (magicl:column v 0)
                                                                        (rotate-vector (magicl:column v 1)))
                                                             (magicl:.* (magicl:column v 1)
                                                                        (rotate-vector (magicl:column v 0)))) 1d0)

                                              (magicl:scale! (magicl:.+
                                                              (magicl:.* (magicl:column v 1)
                                                                         (rotate-vector (magicl:column v 2)))
                                                              (magicl:.* (magicl:column v 2)
                                                                         (rotate-vector (magicl:column v 1)))) 1d0)
                                              (magicl:scale! (magicl:.+
                                                              (magicl:.* (magicl:column v 2)
                                                                         (rotate-vector (magicl:column v 0)))
                                                              (magicl:.* (magicl:column v 0)
                                                                         (rotate-vector (magicl:column v 2)))) 1d0)
                                              ) '(2 6)))))
                (cond
                  ((and
                    (> t1 tol)
                    (> t2 tol)
                    )
                   ;;Apex stress return
                   (setf sig siga)
                   (setf path :apex)
                   )
                  ((and
                    (< f12 tol)
                    (< f13 tol)
                    )
                   (setf sig (magicl:.+ siga (magicl:scale! r1 t1)))
                   (setf path :line-1)
                   ;;line 1
                   )
                  ((and
                    (> f12 tol)
                    (> f13 tol)
                    )
                   (setf sig (magicl:.+ siga (magicl:scale! r2 t2)))
                   ;;line 2
                   (setf path :line-2)
                   )
                  (t
                   ;(break)
                   (setf sig (magicl:.- sig (magicl:scale! rp f)))
                   (setf path :main)
                   ;main
                   )
                  )
                (setf eps-e (magicl:@ Ce sig))

                (setf f (mc-yield-func sig phi c))
                ;; (when (> f (* 2 tol))
                ;;   (error "Mohr-coloumb return misscalculated on path: ~A with an error of f: ~F" path f))

                (let ((pad-eps (magicl:block-matrix (list eps-e
                                                          (cl-mpm/utils:vector-zeros))
                                                    '(2 1))))
                  ;; (setf eps-e (magicl:@ (magicl:inv Q) pad-eps))
                  (setf eps-e (magicl:@ (magicl:linear-solve Q pad-eps)))
                  )

                (setf sig (magicl:@ (magicl:transpose! Q)
                                    (cl-mpm/utils:voigt-from-list
                                     (list
                                      (magicl:tref sig 0 0)
                                      (magicl:tref sig 1 0)
                                      (magicl:tref sig 2 0)
                                      0d0 0d0 0d0))))
                (values
                 sig
                        (swizzle-coombs->voigt eps-e) f)
                )
              (values stress
                      trial-elastic-strain f)))))))

