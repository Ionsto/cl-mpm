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
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(defun linear-elastic-matrix (E nu)
  "Create an isotropic linear elastic matrix"
  (magicl:scale!
    (magicl:from-list (list
                        (- 1d0 nu) nu 0d0
                        nu (- 1d0 nu) 0d0
                        0d0 0d0 (* 0.5 (- 1d0 (* 2 nu))))
                      '(3 3) :type 'double-float)
    (/ E (* (+ 1 nu) (- 1 (* 2 nu))))))

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
    (magicl:.-
     (magicl:@ de strain-increment)
     (magicl:scale! (matrix-to-voight dev-stress) relaxation-const))
    ))

(defun tensile-projection-eig (stress)
  (flet ((H (x) (if (> x 0d0) 1d0 0d0)))
    (multiple-value-bind (l v) (magicl:eig (voight-to-matrix stress))
      (loop for i from 0 to 1
            do (let* ((sii (nth i l)))
                 (setf (nth i l) (* sii (H sii)))))
      (matrix-to-voight (magicl:@ v
                                  (magicl:from-diag l :type 'double-float)
                                  (magicl:transpose v))))))

(defun tensile-projection-Q-mandel (stress)
  (flet ((H (x) (if (> x 0d0) 1d0 0d0)))
    (let ((Q (magicl:zeros '(3 3) :type 'double-float)))
      (multiple-value-bind (l v) (magicl:eig (voight-to-matrix stress))
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
                   ;; (format t "l: ~A ~%" sii)
                   ;; (format t "v: ~A ~%" vii)
                   ;; (format t "~A ~%" comp-prod)
                   ;; (format t "~A ~%" comp)
                   (magicl:.+ Q comp Q)
                   )
              )
      Q))))

(defun tensile-project-q (stress)
  (matrix-to-voight (mandel-to-matrix
                     (magicl:@
                      (tensile-projection-q-mandel stress)
                      (matrix-to-mandel (voight-to-matrix stress))))))

(defun tensile-projection-A (stress damage)
  "Generate a mandel form tensile projection matrix A* from stress"
  (let ((Q (tensile-projection-q stress))
        (I (magicl:from-diag '(1d0 1d0 1d0))))
    (magicl:.- I (magicl:scale Q (- 1 (sqrt (- 1 damage)))))))

(defun test-tensile ()
  (loop for stress in (list
                       (magicl:from-list '(1 0 0) '(3 1) :type 'double-float)
                       (magicl:from-list '(0 1 0) '(3 1) :type 'double-float)
                       (magicl:from-list '(1 1 0) '(3 1) :type 'double-float)
                       (magicl:from-list '(-1 0 0) '(3 1) :type 'double-float)
                       (magicl:from-list '(0 -1 0) '(3 1) :type 'double-float)
                       (magicl:from-list '(-1 -1 0) '(3 1) :type 'double-float)
                       (magicl:from-list '(1 -1 0) '(3 1) :type 'double-float)
                       (magicl:from-list '(-1 1 0) '(3 1) :type 'double-float)
                       (magicl:from-list '(0 0 1) '(3 1) :type 'double-float)
                       (magicl:from-list '(1 0 1) '(3 1) :type 'double-float)
                       (magicl:from-list '(0 1 1) '(3 1) :type 'double-float)
                       (magicl:from-list '(-1 0 1) '(3 1) :type 'double-float)
                       (magicl:from-list '(0 -1 1) '(3 1) :type 'double-float)
                        )
        do (let* ((stress-eig (tensile-projection-eig stress 1d0))
                  (stress-project (tensile-project-q stress))
                 )
             (format t "Stress~%")
             (loop for i from 0 to 2
                   do (format t "~A ~A ~A ~%"
                              (magicl:tref stress i 0)
                              (magicl:tref stress-eig i 0)
                              (magicl:tref stress-project i 0)
                              )))
    ))

(defun tensile-project (stiffness stress damage)
  (if (> damage 0.0d0)
    (let ((damaged-stiffness (magicl:zeros '(2 2) :type 'double-float)))
      (multiple-value-bind (l v) (magicl:eig (voight-to-matrix stress))
        (loop for i from 0 to 1
              do (let* ((sii (nth i l))
                        (vii (magicl::column v i))
                        (scale 1d0)
                        (A (magicl:@ vii (magicl:transpose vii)))
                        )
                   (when (> sii 0d0)
                     (setf scale (sqrt (- 1d0 damage)))
                     (magicl:scale! A scale))
                   (magicl:.+ damaged-stiffness
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
    (magicl:.* tensor mandel-constant)))

(defun maxwell-damage (strain-increment stress elasticity poisson-ratio de viscosity dt damage)
  "A stress increment form of a viscoelastic maxwell material"
  (let* ((order 2)
         (stress-matrix (voight-to-matrix stress))
         (pressure (/ (magicl:trace stress-matrix) 2d0))
         (pressure-matrix (magicl:eye 2 :value pressure))
         (dev-stress (magicl:.- stress-matrix pressure-matrix))
         (relaxation-const (/ (* dt elasticity)
                              (* 2d0 (- 1d0 poisson-ratio) viscosity)))
         )
    (declare (type double-float relaxation-const))
    (magicl:.-
     (magicl:@ de strain-increment)
     (magicl:scale! (matrix-to-voight dev-stress) relaxation-const))
    ;; (matrix-to-voight (tensile-project (voight-to-matrix (magicl:@ de strain-increment)) stress damage))
    ))

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
    (magicl:.+
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
    (magicl:.+
     (magicl:from-list (list pressure pressure 0d0) '(3 1))
     (magicl:scale! strain-dev (/ (* 2d0 viscosity) dt))
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
(declaim (ftype (function (magicl:matrix/double-float
                           magicl:matrix/double-float
                           double-float
                           double-float
                           magicl:matrix/double-float
                           double-float
                           double-float
                           )) maxwell-exp))
(defun maxwell-exp (strain-increment stress elasticity nu de viscosity dt)
  "A stress increment form of a viscoelastic maxwell material"
  (let* (
         (stress-matrix (voight-to-matrix stress))
         (pressure (/ (magicl:trace stress-matrix) 2d0))
         (pressure-matrix (magicl:eye 2 :value pressure))

         (dev-stress (magicl:.- stress-matrix pressure-matrix))
         (rho (/ (* 2 (- 1 nu) viscosity) elasticity))
         (exp-rho (exp (- (/ dt rho))))
         (lam (* (- 1 exp-rho) (/ rho dt)))
         (stress-inc (voight-to-matrix (magicl:@ de strain-increment)))
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
      (magicl:.+ pressure-matrix stress-inc-pressure)
      (magicl:.+ (magicl:scale! dev-stress exp-rho)
                 (magicl:scale! stress-inc-dev lam)))))

(declaim (ftype (function (magicl:matrix/double-float
                           magicl:matrix/double-float
                           double-float
                           double-float
                           magicl:matrix/double-float
                           double-float
                           double-float
                           &optional
                           result-stress
                           )) maxwell-exp-v-simd))
(defun maxwell-exp-v-simd (strain-increment stress elasticity nu de viscosity dt &optional (result-stress nil))
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

(defun deviatoric (voigt)
  (let* ((mat (voight-to-matrix voigt))
         (trace (/ (magicl:trace mat) 2)))
    (matrix-to-voight (magicl:.- mat (magicl:eye 2 :value (the double-float trace))))))
(defun deviatoric-mat (mat)
  (let* ((trace (/ (magicl:trace mat) 2)))
    (matrix-to-voight (magicl:.- mat (magicl:eye 2 :value (the double-float trace))))))


(declaim (ftype (function (magicl:matrix/double-float double-float double-float) double-float) glen-viscosity))
(defun glen-viscosity (stress visc-factor visc-power)
  "Get the viscosity for a given stress state"
  (let* ((dev-stress (deviatoric stress))
         (second-invar (magicl:from-list '(0.5d0 0.5d0 1d0) '(3 1) :type 'double-float))
         (effective-stress (+ 1d-30 (magicl::sum (magicl:.* dev-stress dev-stress second-invar))))
         )
    (declare (type double-float effective-stress))
    (if (> effective-stress 0d0)
        (/ 1d0 (* 2d0 visc-factor (expt effective-stress (* 0.5d0 (- visc-power 1)))))
        0d0)))

(declaim (ftype (function (magicl:matrix/double-float double-float double-float) (values double-float)) glen-viscosity-strain))
(defun glen-viscosity-strain (strain visc-factor visc-power)
  "Get the viscosity for a given strain state"
  (let* (;(effective-strain (+ 1d-15 (cl-mpm/fastmath::voigt-tensor-reduce-simd (deviatoric strain))))
         (effective-strain (+ 1d-15 (cl-mpm/fastmath::voigt-tensor-reduce-simd (deviatoric strain))))
         )
    (declare (type double-float effective-strain))
    (if (> effective-strain 0d0)
        (values (* 0.5d0 visc-factor (expt effective-strain (* 0.5d0 (- (/ 1d0 visc-power) 1d0)))))
        (values 0d0))))


(declaim (ftype (function (magicl:matrix/double-float double-float double-float) (values double-float)) glen-viscosity-strain))
(defun glen-viscosity-stretch (stretch visc-factor visc-power)
  "Get the viscosity for a given strain state"
  (let* (;(effective-strain (+ 1d-15 (cl-mpm/fastmath::voigt-tensor-reduce-simd (deviatoric strain))))
         (dev-stretch (deviatoric-mat stretch))
         (effective-strain (+ 1d-15 (* 0.5d0 (magicl::sum (magicl:.* dev-stretch dev-stretch)))))
         )
    (declare (type double-float effective-strain))
    (if (> effective-strain 0d0)
        (values (* 0.5d0 visc-factor (expt effective-strain (* 0.5d0 (- (/ 1d0 visc-power) 1d0)))))
        (values 0d0))))
