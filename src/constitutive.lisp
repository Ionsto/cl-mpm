(defpackage :cl-mpm/constitutive
  (:use :cl
        :cl-mpm/utils)
  (:export
   #:linear-elastic
   #:linear-elastic-mat
   #:linear-elastic-matrix
   #:newtonian-fluid
   #:maxwell
   #:maxwell-exp
   #:maxwell-exp-v
   #:maxwell-exp-v-simd
   #:norton-hoff
   ))
(in-package :cl-mpm/constitutive)
;; (declaim (optimize (debug 3) (safety 3) (speed 0)))
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(defun linear-elastic-matrix (E nu)
  "Create an isotropic linear elastic matrix"
  (declare (double-float E nu))
    (cl-mpm/fastmaths:fast-scale!
     (magicl:from-list (list
                        (- 1d0 nu) nu nu 0d0 0d0 0d0
                        nu (- 1d0 nu) nu 0d0 0d0 0d0
                        nu nu (- 1d0 nu) 0d0 0d0 0d0
                        0d0 0d0 0d0 (* 0.5d0 (- 1d0 (* 2d0 nu))) 0d0 0d0
                        0d0 0d0 0d0 0d0 (* 0.5d0 (- 1d0 (* 2d0 nu))) 0d0
                        0d0 0d0 0d0 0d0 0d0 (* 0.5d0 (- 1d0 (* 2d0 nu))))
                       '(6 6) :type 'double-float)
     (/ E (* (+ 1d0 nu) (- 1d0 (* 2d0 nu))))))

(defun linear-elastic-principal-matrix (E nu)
  "Create an isotropic linear elastic matrix"
  (declare (double-float E nu))
  (cl-mpm/fastmaths:fast-scale!
   (matrix-from-list (list
                      (- 1d0 nu) nu nu
                      nu (- 1d0 nu) nu
                      nu nu (- 1d0 nu)))
   (/ E (* (+ 1d0 nu) (- 1d0 (* 2d0 nu))))))

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

(defun linear-elastic-mat (strain elastic-matrix &optional result)
  "Isotropic linear-elastic constitutive model"
  ;(cl-mpm/utils:@-tensor-voigt elastic-matrix strain)
  ;; (let ((result (if result
  ;;                   result
  ;;                   (cl-mpm/utils:voigt-zeros))))
  ;   (magicl:mult elastic-matrix strain :target result)
  ;;   result)
  (let ((result (if result
                    result
                    (cl-mpm/utils:voigt-zeros))))
    (cl-mpm/fastmaths::fast-@-tensor-voigt elastic-matrix strain result)
    result))

(declaim (ftype (function (magicl:matrix/double-float double-float double-float)
                          magicl:matrix/double-float
                          ) linear-elastic))
(defun linear-elastic (strain E nu)
  "Isotropic linear-elastic constitutive model"
   (cl-mpm/fastmaths::fast-@-tensor-voigt (linear-elastic-matrix E nu) strain))

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
  (cl-mpm/fastmaths::fast-.+
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
                    (cl-mpm/fastmaths::fast-.+ pressure-matrix
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
                   (cl-mpm/fastmaths::fast-.+ Q comp Q)
                   )
              )
      Q))))

(defun tensile-projection-Q-cw-mandel (strain)
  (flet ((H (x) (if (> x 0d0) 1d0 0d0)))
    (let ((Q (magicl:zeros '(3 3) :type 'double-float)))
      (multiple-value-bind (l v) (cl-mpm/utils::eig (voight-to-matrix (cl-mpm/fastmaths::fast-.* strain
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
                   (cl-mpm/fastmaths::fast-.+ Q comp-prod Q)
                   )
              )
        (let* ((si (nth 0 l))
               (sj (nth 1 l))
               (vi (magicl::column v 0))
               (vj (magicl::column v 1))
               (vij (magicl:@ vi (magicl:transpose vj)))
               (vji (magicl:@ vj (magicl:transpose vi)))
               (pij (magicl:scale!
                     (cl-mpm/fastmaths::fast-.+ vij vji)
                     0.5d0))
               (comp-prod
                 (magicl:scale! (magicl:@
                                (matrix-to-mandel pij)
                                (magicl:transpose (matrix-to-mandel pij)))
                               (+ (H si) (H sj)))))
          (cl-mpm/fastmaths::fast-.+ Q comp-prod Q))
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
    (cl-mpm/fastmaths::fast-.+ I (magicl:scale Q (- (sqrt (- 1d0 damage)) 1d0)))))

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
                   (cl-mpm/fastmaths::fast-.+ damaged-stiffness
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
    (cl-mpm/fastmaths::fast-.* tensor mandel-constant))
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
    (cl-mpm/fastmaths::fast-.+
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
    (cl-mpm/fastmaths::fast-.+
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
  (let* ((stress-matrix (voight-to-matrix stress))
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
    (cl-mpm/fastmaths::fast-.+
     (matrix-to-voight (cl-mpm/fastmaths::fast-.+ pressure-matrix stress-inc-pressure))
               (cl-mpm/fastmaths::fast-.+
                (magicl:scale! (matrix-to-voight dev-stress) exp-rho)
                (magicl:scale! (matrix-to-voight stress-inc-dev) lam)))))


(defun maxwell-exp-v (strain-increment stress elasticity nu de viscosity dt &key (result-stress))
  "A stress increment form of a viscoelastic maxwell material"
  (let* ((pressure (/ (voight-trace stress) 3d0))
         (pressure-matrix (voight-eye pressure))
         (dev-stress (magicl:.- stress pressure-matrix))
         (rho (/ (* 2d0 (- 1d0 nu) viscosity) elasticity))
         (exp-rho (exp (- (/ dt rho))))
         (lam (* (- 1 exp-rho) (/ rho dt)))
         (stress-inc (magicl:@ (linear-elastic-matrix elasticity nu) strain-increment))
         (stress-inc-pressure (voight-eye (/ (voight-trace stress-inc) 3d0)))
         (stress-inc-dev (magicl:.- stress-inc stress-inc-pressure)))
    (declare (double-float rho exp-rho lam viscosity elasticity dt nu))
    ;(declare (dynamic-extent pressure-matrix dev-stress stress-inc stress-inc-pressure stress-inc-dev))
    (cl-mpm/fastmaths::fast-.+
     (cl-mpm/fastmaths::fast-.+ pressure-matrix stress-inc-pressure)
     (cl-mpm/fastmaths::fast-.+ (cl-mpm/fastmaths:fast-scale! dev-stress exp-rho)
                               (cl-mpm/fastmaths:fast-scale! stress-inc-dev lam)))))

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
    (cl-mpm/fastmaths:fast-add result-stress pressure-matrix)
    (cl-mpm/fastmaths:fast-add result-stress stress-inc-pressure)
    (cl-mpm/fastmaths:fast-add result-stress stress-inc-pressure)
    (cl-mpm/fastmaths:fast-add result-stress stress-inc-pressure)
    (cl-mpm/fastmaths:fast-fmacc result-stress dev-stress exp-rho)
    (cl-mpm/fastmaths:fast-fmacc result-stress stress-inc-dev lam)
    result-stress
    ;; (cl-mpm/fastmaths::fast-.+
    ;;  (cl-mpm/fastmaths::fast-.+ pressure-matrix stress-inc-pressure result-stress)
    ;;  (cl-mpm/fastmaths::fast-.+ (magicl:scale! dev-stress exp-rho)
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
                                                        (expt (magicl::sum (cl-mpm/fastmaths::fast-.* dev-stress dev-stress
                                                                                        (magicl:from-list
                                                                                         '(0.5d0 0.5d0 1d0) '(3 1))))
                                                              (* 0.5 (- visc-power 1)))))))
    (cl-mpm/fastmaths::fast-.+ stress
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
                                                        (expt (magicl::sum (cl-mpm/fastmaths::fast-.* dev-stress dev-stress
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
    (cl-mpm/fastmaths::fast-.+ (matrix-to-voight pressure-matrix) dev-stress)))

(defun glen-stress (strain visc-factor visc-power dt)
  (let* ((strain-trace (/ (magicl:trace (voight-to-matrix strain)) 3d0))
         (dev-strain (matrix-to-voight (magicl:.- (voight-to-matrix strain) (magicl:eye 2 :value strain-trace))))
         (second-invar (magicl:from-list '(0.5d0 0.5d0 1d0) '(3 1) :type 'double-float))
         (effective-strain (magicl::sum (cl-mpm/fastmaths::fast-.* dev-strain dev-strain second-invar)))
         )
    (if (> effective-strain 0d0)
        (magicl:scale dev-strain (* visc-factor (expt effective-strain
                                                 (* 0.5d0 (- (/ 1d0 visc-power) 1d0)))))
        (magicl:scale dev-strain 0d0))))

(defun deviatoric (voigt)
  (deviatoric-voigt voigt)
  ;; (let* ((mat (voight-to-matrix voigt))
  ;;        (trace (/ (magicl:trace mat) 3d0)))
  ;;   (matrix-to-voight (magicl:.- mat (magicl:eye 3 :value (the double-float trace)))))
  )

(defun deviatoric-mat (mat)
  (let* ((trace (/ (magicl:trace mat) 3d0)))
    (matrix-to-voight (magicl:.- mat (magicl:eye 3 :value (the double-float trace))))))


(declaim (ftype (function (magicl:matrix/double-float double-float double-float) double-float) glen-viscosity))
(defun glen-viscosity (stress visc-factor visc-power)
  "Get the viscosity for a given stress state"
  (let* ((s-dev (deviatoric-voigt stress))
         (second-invar (cl-mpm/utils:voigt-from-list (list 1d0 1d0 1d0 2d0 2d0 2d0)))
         (visc-factor (expt visc-factor (- visc-power)))
         (effective-stress (+ 1d-20 (* 0.5d0
                                       (cl-mpm/fastmaths::fast-sum
                                        (cl-mpm/fastmaths:fast-.*
                                         (cl-mpm/fastmaths:fast-.* s-dev s-dev)
                                         second-invar))))))
    (declare (type double-float effective-stress))
    (if (> effective-stress 0d0)
        (/ 1d0 (* 2d0 visc-factor (expt effective-stress (* 0.5d0 (- visc-power 1)))))
        0d0)))

(defun effective-strain-rate (strain)
  (sqrt
   (* 0.5d0
      (cl-mpm/fastmaths::voigt-tensor-reduce-simd
       (deviatoric strain)))))

(declaim (ftype (function (magicl:matrix/double-float double-float double-float) double-float) glen-viscosity-strain))
(defun glen-viscosity-strain (strain visc-factor visc-power)
  "Get the viscosity for a given strain state"
  (let* ((effective-strain (+ 1d-20 (effective-strain-rate strain))))
    (declare (type double-float effective-strain))
    (* visc-factor (the double-float (expt effective-strain (* 1d0 (- (/ 1d0 visc-power) 1d0)))))))

(defun glen-viscosity-stress (stress visc-factor visc-power)
  (declare (double-float visc-factor visc-power))
  "Get the viscosity for a given strain state"
  (let* ((A (expt visc-factor (/ -1d0 visc-power)))
         (effective-strain (+ 1d-50 (sqrt (* 0.5d0 (expt (cl-mpm/utils::trace-voigt (deviatoric-voigt stress)) 2d0))))))
    (declare (type double-float effective-strain A))
    (if (> effective-strain 0d0)
        (/ 1d0 (* 2d0 A (the double-float (expt effective-strain (- visc-power 1d0)))))
        0d0)))


(declaim (ftype (function (magicl:matrix/double-float double-float double-float) (values double-float)) glen-viscosity-strain))
(defun glen-viscosity-stretch (stretch visc-factor visc-power)
  "Get the viscosity for a given strain state"
  (let* (;(effective-strain (+ 1d-15 (cl-mpm/fastmaths::voigt-tensor-reduce-simd (deviatoric strain))))
         (dev-stretch (deviatoric-voigt stretch))
         (effective-strain (+ 1d-15 (* 0.5d0 (magicl::sum (cl-mpm/fastmaths::fast-.* dev-stretch dev-stretch))))))
    (declare (type double-float effective-strain))
    (if (> effective-strain 0d0)
        (values (* 0.5d0 visc-factor (expt effective-strain (* 0.5d0 (- (/ 1d0 visc-power) 1d0)))))
        (values 0d0))))

(defun voigt-j2 (s)
  "Calculate j2 invarient from deviatoric stress"
  (let ((storage (magicl::matrix/double-float-storage s)))
    (/ (+ (the double-float (cl-mpm/fastmaths:dot s s))
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
         (dj2 (cl-mpm/fastmaths::fast-.* s (cl-mpm/utils::voigt-from-list '(1d0 1d0 1d0 2d0 2d0 2d0))))
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
         (inc 0d0)
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
                                   (cl-mpm/fastmaths::fast-.+ (magicl:eye 6)
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
                         (let ((b-eps (cl-mpm/fastmaths::fast-.+ eps-e
                                                 (magicl:scale trial-elastic-strain -1d0)
                                                 (magicl:scale! df dgam)
                                                 )
                                      )
                               (b-f (vm-yield-func j2 rho)))
                           (loop for i from 0 below 6
                                 do (setf (magicl:tref b i 0) (magicl:tref b-eps i 0)))
                           (setf (magicl:tref b 6 0) b-f)))
                       ))
          (setf inc (the double-float
                           (cl-mpm/fastmaths::voigt-von-mises
                            (cl-mpm/fastmaths::fast-.--voigt
                             eps-e
                             trial-elastic-strain))))
          (values sig eps-e f inc)
          ))
      (values sig trial-elastic-strain f inc))))



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
                    (magicl:tref vec 5 0)
                    (magicl:tref vec 3 0)
                    (magicl:tref vec 4 0)
                    )))
(defun swizzle-coombs->voigt (vec &optional res)
  (let ((res (if res
                 res
                 (cl-mpm/utils:voigt-zeros))))
    (setf
     (varef res 0) (varef vec 0)
     (varef res 1) (varef vec 1)
     (varef res 2) (varef vec 2)
     (varef res 3) (varef vec 4)
     (varef res 4) (varef vec 5)
     (varef res 5) (varef vec 3))
    res))

(defun dp-yield-mc-circumscribe (stress phi c)
  (declare (double-float c phi))
  (let* ((factor (sqrt (+ 3d0 (expt (sin phi) 2))))
         (a (/ (* (sqrt 3) c (cos phi)) factor))
         (b (/ (* (sqrt 3) (sin phi)) factor))
         (p (/ (cl-mpm/utils::trace-voigt stress) 3d0))
         (j2 (cl-mpm/constitutive::voigt-j2
              (cl-mpm/utils::deviatoric-voigt stress))))
    (declare (double-float a b j2 p))
    (+
     (the double-float (sqrt j2))
     (- a)
     (* b p))))

(defun dp-yield-mc-dp2 (stress phi c)
  (declare (double-float c phi))
  (let* ((a (* c (cos phi)))
         (b (sin phi))
         (p (cl-mpm/utils::trace-voigt stress))
         (j2 (cl-mpm/constitutive::voigt-j2
              (cl-mpm/utils::deviatoric-voigt stress))))
    (declare (double-float a b j2 p))
    (+
     (the double-float (sqrt j2))
     (- a)
     (* b p))))

(defun fast-mc (stress phi c)
  (multiple-value-bind (s1 s2 s3) (cl-mpm/utils::principal-stresses-3d stress)
      (mc-yield-func (cl-mpm/utils:vector-from-list (list s1 s2 s3)) phi c)))


;; (defun plastic-dp (trial-elastic-strain de E nu phi psi c)
;;   (cl-mpm/ext::constitutive-drucker-prager trial-elastic-strain de E nu phi psi c))

(declaim (notinline  mc-plastic))
(defun mc-plastic (stress de trial-elastic-strain E nu phi psi c)
  (declare (optimize (speed 3) (safety 0) (debug 0)))
  (declare (double-float E nu phi psi c)
           (magicl:matrix/double-float stress de trial-elastic-strain))
  (let* ((tol 1d-9)
         (initial-f 0d0))
    (let (;(f-dp (dp-yield-mc-circumscribe stress phi c))
          (f-dp (fast-mc stress phi c))
          )
      ;;Early check for if we should yield - DP eval is much faster?
      (if (> f-dp tol)
          (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix trial-elastic-strain))
            (let* ((l-sort (sort (mapcar #'cons l (list 0 1 2)) #'> :key #'car))
                   (l (mapcar #'car l-sort))
                   (v (magicl:block-matrix (list
                                            (magicl:column v (cdr (nth 0 l-sort)))
                                            (magicl:column v (cdr (nth 1 l-sort)))
                                            (magicl:column v (cdr (nth 2 l-sort)))) '(1 3))))
              (let* ((De3
                       (cl-mpm/fastmaths::fast-scale!
                        (cl-mpm/utils:matrix-from-list (list
                                                        (- 1d0 nu) nu nu
                                                        nu (- 1d0 nu) nu
                                                        nu nu (- 1d0 nu)))
                        (/ E (* (+ 1d0 nu) (- 1d0 (* 2d0 nu))))))
                     (epsTr (cl-mpm/utils:vector-from-list l))
                     (sig (cl-mpm/utils:@-mat-vec De3 epsTr))
                     (f (mc-yield-func sig phi c))
                     (initial-f f))
                ;; (when (and
                ;;        (> f tol)
                ;;        (not (> f-dp tol)))
                ;;   (error "DP misses mc case ~E ~E" f f-dp))
                (if (> f tol)
                    (let* (
                           (Ce (magicl:inv De3))
                           (eps-e (cl-mpm/utils::vector-copy epsTr))
                           ;; (eps-e epsTr)
                           (k (/ (+ 1 (sin phi)) (- 1d0 (sin phi))))
                           (sigc (* 2d0 c (sqrt k)))
                           (m (/ (+ 1 (sin psi)) (- 1d0 (sin psi))))
                           (siga (magicl:scale!
                                  (cl-mpm/utils:vector-from-list (list 1d0 1d0 1d0))
                                  (/ sigc (- k 1d0))
                                  ))
                           (r1 (vector-from-list (list 1d0 1d0 k)))
                           (r2 (vector-from-list (list 1d0 k k)))
                           (rg1 (vector-from-list (list 1d0 1d0 m)))
                           (rg2 (vector-from-list (list 1d0 m m)))
                           (df (vector-from-list (list k 0d0 -1d0)))
                           (dg (vector-from-list (list m 0d0 -1d0)))
                           (rp (magicl:scale! (cl-mpm/utils:@-mat-vec De3 dg)
                                              (/ 1d0 (the double-float
                                                          (magicl:tref (magicl:@ (magicl:transpose dg) De3 df) 0 0)
                                                          ))))
                           (t1 (/ (magicl:tref (magicl:@ (magicl:transpose rg1) Ce (magicl:.- sig siga)) 0 0)
                                  (magicl:tref (magicl:@ (magicl:transpose rg1) Ce r1) 0 0)))
                           (t2 (/ (magicl:tref (magicl:@
                                                (magicl:transpose rg2)
                                                Ce (magicl:.- sig siga)) 0 0)
                                  (magicl:tref (magicl:@ (magicl:transpose rg2) Ce r2) 0 0)))
                           (f12 (magicl:tref (magicl:@ (magicl:transpose! (cl-mpm/fastmaths::cross-product rp r1))
                                                       (magicl:.- sig siga)) 0 0))
                           (f13 (magicl:tref (magicl:@ (magicl:transpose! (cl-mpm/fastmaths::cross-product rp r2))
                                                       (magicl:.- sig siga)) 0 0))
                           (path :no-return)
                           (Q
                             (magicl:transpose!
                              (magicl:block-matrix
                               (list
                                (cl-mpm/fastmaths::fast-.* (magicl:column v 0)
                                                           (magicl:column v 0))

                                (cl-mpm/fastmaths::fast-.* (magicl:column v 1)
                                                           (magicl:column v 1))
                                (cl-mpm/fastmaths::fast-.* (magicl:column v 2)
                                                           (magicl:column v 2))

                                (magicl:scale!
                                 (cl-mpm/fastmaths::fast-.* (magicl:column v 0)
                                                            (magicl:column v 1)) 2d0)
                                (magicl:scale!
                                 (cl-mpm/fastmaths::fast-.* (magicl:column v 1)
                                                            (magicl:column v 2)) 2d0)
                                (magicl:scale!
                                 (cl-mpm/fastmaths::fast-.* (magicl:column v 2)
                                                            (magicl:column v 0)) 2d0)

                                (cl-mpm/fastmaths::fast-.* (magicl:column v 0)
                                                           (rotate-vector (magicl:column v 0)))
                                (cl-mpm/fastmaths::fast-.* (magicl:column v 1)
                                                           (rotate-vector (magicl:column v 1)))
                                (cl-mpm/fastmaths::fast-.* (magicl:column v 2)
                                                           (rotate-vector (magicl:column v 2)))

                                (magicl:scale! (cl-mpm/fastmaths::fast-.+
                                                (cl-mpm/fastmaths::fast-.*
                                                 (magicl:column v 0)
                                                 (rotate-vector (magicl:column v 1)))
                                                (cl-mpm/fastmaths::fast-.*
                                                 (magicl:column v 1)
                                                 (rotate-vector (magicl:column v 0)))) 1d0)

                                (magicl:scale! (cl-mpm/fastmaths::fast-.+
                                                (cl-mpm/fastmaths::fast-.*
                                                 (magicl:column v 1)
                                                 (rotate-vector (magicl:column v 2)))
                                                (cl-mpm/fastmaths::fast-.*
                                                 (magicl:column v 2)
                                                 (rotate-vector (magicl:column v 1)))) 1d0)
                                (magicl:scale! (cl-mpm/fastmaths::fast-.+
                                                (cl-mpm/fastmaths::fast-.*
                                                 (magicl:column v 2)
                                                 (rotate-vector (magicl:column v 0)))
                                                (cl-mpm/fastmaths::fast-.*
                                                 (magicl:column v 0)
                                                 (rotate-vector (magicl:column v 2)))) 1d0)
                                ) '(2 6))))
                           (psinc 0d0)
                           )
                      (declare (double-float t1 t2 f12 f13))
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
                         (setf sig (cl-mpm/fastmaths::fast-.+ siga (magicl:scale! r1 t1)))
                         (setf path :line-1)
                         ;;line 1
                         )
                        ((and
                          (> f12 tol)
                          (> f13 tol)
                          )
                         (setf sig (cl-mpm/fastmaths::fast-.+ siga (magicl:scale! r2 t2)))
                         ;;line 2
                         (setf path :line-2)
                         )
                        (t
                         (setf sig (magicl:.- sig (magicl:scale! rp f)))
                         (setf path :main)
                                        ;main
                         )
                        )
                      (setf eps-e (cl-mpm/fastmaths::fast-@-matrix-vector Ce sig eps-e))

                      (setf f (mc-yield-func sig phi c))
                      (when (> f (* 10000d0 tol))
                        (error "Mohr-coloumb return misscalculated on path: ~A with an error of f: ~F" path f))
                      (setf psinc (the double-float
                                       (cl-mpm/fastmaths::vector-principal-von-mises
                                        (cl-mpm/fastmaths::fast-.--vector
                                         eps-e
                                         epstr))))
                      ;; (break)
                      (let ((pad-eps
                              (cl-mpm/utils:voigt-from-list
                               (list
                                (varef eps-e 0)
                                (varef eps-e 1)
                                (varef eps-e 2)
                                0d0
                                0d0
                                0d0))))
                        (setf eps-e (magicl:linear-solve Q pad-eps)))
                      (setf sig (cl-mpm/fastmaths::fast-@-tensor-voigt de eps-e))
                      (values
                       (swizzle-coombs->voigt sig stress)
                       (swizzle-coombs->voigt eps-e trial-elastic-strain)
                       ;; sig
                       ;; eps-e
                       initial-f
                       ;; 0d0
                       psinc
                       ))
                    ;;No MC yield - just return
                    (values stress
                            trial-elastic-strain
                            initial-f
                            0d0
                            )))))
          ;;No DP yield - just return
          (values stress
                  trial-elastic-strain
                  f-dp
                  0d0)))))
(defun mc-plastic-terzaghi (stress de trial-elastic-strain E nu phi psi c pore-pressure)
  (declare (optimize (speed 3) (safety 0) (debug 0)))
  (declare (double-float E nu phi psi c)
           (magicl:matrix/double-float stress de trial-elastic-strain))
  (let* ((tol 1d-9)
         (initial-f 0d0))
    (let ((f-dp (dp-yield-mc-circumscribe
                 (cl-mpm/fastmaths:fast-.+ stress (cl-mpm/utils::voigt-eye pore-pressure))
                 phi c)))
      ;;Early check for if we should yield - DP eval is much faster?
      (if (> f-dp tol)
          (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix trial-elastic-strain))
            (let* ((l-sort (sort (mapcar #'cons l (list 0 1 2)) #'> :key #'car))
                   (l (mapcar #'car l-sort))
                   (v (magicl:block-matrix (list
                                            (magicl:column v (cdr (nth 0 l-sort)))
                                            (magicl:column v (cdr (nth 1 l-sort)))
                                            (magicl:column v (cdr (nth 2 l-sort)))) '(1 3))))
              (let* ((De3
                       (cl-mpm/fastmaths::fast-scale!
                        (cl-mpm/utils:matrix-from-list (list
                                                        (- 1d0 nu) nu nu
                                                        nu (- 1d0 nu) nu
                                                        nu nu (- 1d0 nu)))
                        (/ E (* (+ 1d0 nu) (- 1d0 (* 2d0 nu))))))
                     (epsTr (cl-mpm/utils:vector-from-list l))
                     (sig (cl-mpm/utils:@-mat-vec De3 epsTr))
                     (sig (cl-mpm/fastmaths:fast-.+ sig (cl-mpm/utils::voigt-eye pore-pressure) sig))
                     (f (mc-yield-func sig phi c))
                     (initial-f f))
                ;; (format t "f: ~E - Fdp: ~E~%" f f-dp)
                ;; (when (and
                ;;        (> f tol)
                ;;        (not (> f-dp tol)))
                ;;   (error "DP misses mc case ~E ~E" f f-dp))
                (if (> f tol)
                    (let* (
                           (Ce (magicl:inv De3))
                           (eps-e (cl-mpm/utils::vector-copy epsTr))
                           ;; (eps-e epsTr)
                           (k (/ (+ 1 (sin phi)) (- 1d0 (sin phi))))
                           (sigc (* 2d0 c (sqrt k)))
                           (m (/ (+ 1 (sin psi)) (- 1d0 (sin psi))))
                           (siga (magicl:scale! (cl-mpm/utils:vector-from-list (list 1d0 1d0 1d0))
                                                (/ sigc (- k 1d0))
                                                ))
                           (r1 (vector-from-list (list 1d0 1d0 k)))
                           (r2 (vector-from-list (list 1d0 k k)))
                           (rg1 (vector-from-list (list 1d0 1d0 m)))
                           (rg2 (vector-from-list (list 1d0 m m)))
                           (df (vector-from-list (list k 0d0 -1d0)))
                           (dg (vector-from-list (list m 0d0 -1d0)))
                           (rp (magicl:scale! (cl-mpm/utils:@-mat-vec De3 dg)
                                              (/ 1d0 (the double-float
                                                          (magicl:tref (magicl:@ (magicl:transpose dg) De3 df) 0 0)
                                                          ))))
                           (t1 (/ (magicl:tref (magicl:@ (magicl:transpose rg1) Ce (magicl:.- sig siga)) 0 0)
                                  (magicl:tref (magicl:@ (magicl:transpose rg1) Ce r1) 0 0)))
                           (t2 (/ (magicl:tref (magicl:@ (magicl:transpose rg2) Ce (magicl:.- sig siga)) 0 0)
                                  (magicl:tref (magicl:@ (magicl:transpose rg2) Ce r2) 0 0)))
                           (f12 (magicl:tref (magicl:@ (magicl:transpose! (cl-mpm/fastmaths::cross-product rp r1))
                                                       (magicl:.- sig siga)) 0 0))
                           (f13 (magicl:tref (magicl:@ (magicl:transpose! (cl-mpm/fastmaths::cross-product rp r2))
                                                       (magicl:.- sig siga)) 0 0))
                           (path :no-return)
                           (Q
                             (magicl:transpose!
                              (magicl:block-matrix (list
                                                    (cl-mpm/fastmaths::fast-.* (magicl:column v 0)
                                                                              (magicl:column v 0))

                                                    (cl-mpm/fastmaths::fast-.* (magicl:column v 1)
                                                                              (magicl:column v 1))
                                                    (cl-mpm/fastmaths::fast-.* (magicl:column v 2)
                                                                              (magicl:column v 2))

                                                    (magicl:scale!
                                                     (cl-mpm/fastmaths::fast-.* (magicl:column v 0)
                                                                               (magicl:column v 1)) 2d0)
                                                    (magicl:scale!
                                                     (cl-mpm/fastmaths::fast-.* (magicl:column v 1)
                                                                               (magicl:column v 2)) 2d0)
                                                    (magicl:scale!
                                                     (cl-mpm/fastmaths::fast-.* (magicl:column v 2)
                                                                               (magicl:column v 0)) 2d0)

                                                    (cl-mpm/fastmaths::fast-.* (magicl:column v 0)
                                                                              (rotate-vector (magicl:column v 0)))
                                                    (cl-mpm/fastmaths::fast-.* (magicl:column v 1)
                                                                              (rotate-vector (magicl:column v 1)))
                                                    (cl-mpm/fastmaths::fast-.* (magicl:column v 2)
                                                                              (rotate-vector (magicl:column v 2)))

                                                    (magicl:scale! (cl-mpm/fastmaths::fast-.+
                                                                    (cl-mpm/fastmaths::fast-.* (magicl:column v 0)
                                                                                              (rotate-vector (magicl:column v 1)))
                                                                    (cl-mpm/fastmaths::fast-.* (magicl:column v 1)
                                                                                              (rotate-vector (magicl:column v 0)))) 1d0)

                                                    (magicl:scale! (cl-mpm/fastmaths::fast-.+
                                                                    (cl-mpm/fastmaths::fast-.* (magicl:column v 1)
                                                                                              (rotate-vector (magicl:column v 2)))
                                                                    (cl-mpm/fastmaths::fast-.* (magicl:column v 2)
                                                                                              (rotate-vector (magicl:column v 1)))) 1d0)
                                                    (magicl:scale! (cl-mpm/fastmaths::fast-.+
                                                                    (cl-mpm/fastmaths::fast-.* (magicl:column v 2)
                                                                                              (rotate-vector (magicl:column v 0)))
                                                                    (cl-mpm/fastmaths::fast-.* (magicl:column v 0)
                                                                                              (rotate-vector (magicl:column v 2)))) 1d0)
                                                    ) '(2 6)))))
                      (declare (double-float t1 t2 f12 f13))
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
                         (setf sig (cl-mpm/fastmaths::fast-.+ siga (magicl:scale! r1 t1)))
                         (setf path :line-1)
                         ;;line 1
                         )
                        ((and
                          (> f12 tol)
                          (> f13 tol)
                          )
                         (setf sig (cl-mpm/fastmaths::fast-.+ siga (magicl:scale! r2 t2)))
                         ;;line 2
                         (setf path :line-2)
                         )
                        (t
                         (setf sig (magicl:.- sig (magicl:scale! rp f)))
                         (setf path :main)
                                        ;main
                         )
                        )


                      (setf f (mc-yield-func sig phi c))
                      (when (> f (* 10000d0 tol))
                        (error "Mohr-coloumb return misscalculated on path: ~A with an error of f: ~F" path f))

                      ;; (cl-mpm/fastmaths:fast-.- sig (cl-mpm/utils::voigt-eye pore-pressure) sig)
                      (setf eps-e (magicl:@ Ce sig))

                      ;; (break)
                      (let ((pad-eps (magicl:block-matrix (list eps-e
                                                                (cl-mpm/utils:vector-zeros))
                                                          '(2 1))))
                        ;; (setf eps-e (magicl:@ (magicl:inv Q) pad-eps))
                        (setf eps-e (magicl:linear-solve Q pad-eps))
                        )

                      (setf sig (magicl:@ (magicl:transpose! Q)
                                          (cl-mpm/utils:voigt-from-list
                                           (list
                                            (magicl:tref sig 0 0)
                                            (magicl:tref sig 1 0)
                                            (magicl:tref sig 2 0)
                                            0d0 0d0 0d0))))

                      (values
                       ;; sig
                       (swizzle-coombs->voigt sig)
                       (swizzle-coombs->voigt eps-e) initial-f)
                      )
                    ;;No MC yield - just return
                    (values stress trial-elastic-strain initial-f)))))
          ;;No DP yield - just return
          (values stress
                  trial-elastic-strain f-dp)))))



(defun Q-matrix (v)
  (magicl:transpose!
   (magicl:block-matrix
    (list
     (cl-mpm/fastmaths::fast-.* (magicl:column v 0) (magicl:column v 0))

     (cl-mpm/fastmaths::fast-.* (magicl:column v 1) (magicl:column v 1))
     (cl-mpm/fastmaths::fast-.* (magicl:column v 2) (magicl:column v 2))

     (magicl:scale!
      (cl-mpm/fastmaths::fast-.* (magicl:column v 0) (magicl:column v 1)) 2d0)
     (magicl:scale!
      (cl-mpm/fastmaths::fast-.* (magicl:column v 1) (magicl:column v 2)) 2d0)
     (magicl:scale!
      (cl-mpm/fastmaths::fast-.* (magicl:column v 2) (magicl:column v 0)) 2d0)

     (cl-mpm/fastmaths::fast-.* (magicl:column v 0) (rotate-vector (magicl:column v 0)))
     (cl-mpm/fastmaths::fast-.* (magicl:column v 1) (rotate-vector (magicl:column v 1)))
     (cl-mpm/fastmaths::fast-.* (magicl:column v 2) (rotate-vector (magicl:column v 2)))

     (magicl:scale! (cl-mpm/fastmaths::fast-.+
                     (cl-mpm/fastmaths::fast-.* (magicl:column v 0)
                                               (rotate-vector (magicl:column v 1)))
                     (cl-mpm/fastmaths::fast-.* (magicl:column v 1)
                                               (rotate-vector (magicl:column v 0)))) 1d0)

     (magicl:scale! (cl-mpm/fastmaths::fast-.+
                     (cl-mpm/fastmaths::fast-.* (magicl:column v 1)
                                               (rotate-vector (magicl:column v 2)))
                     (cl-mpm/fastmaths::fast-.* (magicl:column v 2)
                                               (rotate-vector (magicl:column v 1)))) 1d0)
     (magicl:scale! (cl-mpm/fastmaths::fast-.+
                     (cl-mpm/fastmaths::fast-.* (magicl:column v 2)
                                               (rotate-vector (magicl:column v 0)))
                     (cl-mpm/fastmaths::fast-.* (magicl:column v 0)
                                               (rotate-vector (magicl:column v 2)))) 1d0)
     ) '(2 6))))

(declaim (notinline plastic-dp))
(defun plastic-dp (stress de trial-elastic-strain E nu phi psi c)
  ;; (declare (optimize (speed 3) (safety 0) (debug 0)))
  (declare (optimize (speed 0) (safety 3) (debug 3)))
  (declare (double-float E nu phi psi c)
           (magicl:matrix/double-float stress de trial-elastic-strain))
  (let* ((tol 1d-9)
         (initial-f 0d0))

    (let ((f-dp (dp-yield-mc-circumscribe stress phi c)))
      ;;Early check for if we should yield - DP eval is much faster?
      (if t;(> f-dp tol)
          (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix trial-elastic-strain))
            (let* ((l-sort (sort (mapcar #'cons l (list 0 1 2)) #'> :key #'car))
                   (l (mapcar #'car l-sort))
                   (v (magicl:block-matrix (list
                                            (magicl:column v (cdr (nth 0 l-sort)))
                                            (magicl:column v (cdr (nth 1 l-sort)))
                                            (magicl:column v (cdr (nth 2 l-sort)))) '(1 3))))
              (let* ((alfa (- (tan phi)))
                     (bta (- (tan psi)))
                     (xsic (* (sqrt 3) (/ 1d0 (tan phi)) c))
                     (De3
                       (cl-mpm/fastmaths::fast-scale!
                        (cl-mpm/utils:matrix-from-list (list
                                                        (- 1d0 nu) nu nu
                                                        nu (- 1d0 nu) nu
                                                        nu nu (- 1d0 nu)))
                        (/ E (* (+ 1d0 nu) (- 1d0 (* 2d0 nu))))))
                     ;; (Ce (magicl:inv De3))
                     (Ce (cl-mpm/fastmaths:fast-scale!
                          (cl-mpm/utils:matrix-from-list
                                                       (list
                                                        1d0 (- nu) (- nu)
                                                        (- nu) 1d0 (- nu)
                                                        (- nu) (- nu) 1d0))
                          (/ 1d0 E)))
                     (epsTr (cl-mpm/utils:vector-from-list l))
                     (sig (magicl:.-
                           (cl-mpm/utils:@-mat-vec De3 epsTr)
                           (/ xsic (sqrt 3))))
                     (epsE (magicl:@ Ce sig))
                     (xi (/ (cl-mpm/fastmaths:fast-sum sig) (sqrt 3)))
                     (s (magicl:.-
                         sig
                         (/ xi (sqrt 3))))
                     (rho (sqrt (cl-mpm/fastmaths:dot s s)))
                     (f (- rho (* alfa xi)))
                     (initial-f f))
                (if (> f tol)
                    (let* ((epsTr (cl-mpm/utils:vector-copy epsE))
                           (fap (if (not (= bta 0d0))
                                     (+
                                      (* rho (sqrt (+ 1d0 nu)))
                                      (/
                                       (* xi (sqrt (- 1d0 (* 2d0 nu))))
                                       (* bta
                                          (/
                                           (sqrt (+ 1d0 nu))
                                           (sqrt (- 1d0 (* 2d0 nu)))))))
                                    (* xi sb-ext:double-float-negative-infinity)))
                           (path :no-path)
                           (Q (Q-matrix v)))
                      (cond
                        ((< fap tol)
                         (setf path :apex-return)
                         (setf sig (cl-mpm/fastmaths:fast-scale!
                                    (cl-mpm/utils:vector-from-list (list 1d0 1d0 1d0))
                                    (/ xsic (sqrt 3))))
                         (setf epsE (magicl:@ Ce sig)))
                        (t
                         (setf path :surface-return)
                         (let ((tolf 1d-6)
                               (b (magicl:zeros '(4 1)))
                               (dgam 0d0)
                               (df (magicl:.-
                                    (cl-mpm/fastmaths:fast-scale-vector s (/ 1d0 rho))
                                    (/ alfa (sqrt 3))))
                               (dg (magicl:.-
                                    (cl-mpm/fastmaths:fast-scale-vector s (/ 1d0 rho))
                                    (/ bta (sqrt 3))))
                               (dgg
                                 (magicl:.-
                                  (cl-mpm/fastmaths:fast-scale!
                                   (magicl:.-
                                    (magicl:eye 3 :value 3)
                                    (magicl:ones (list 3 3))
                                    )
                                   (/ 1d0 (* 3d0 rho)))
                                  (cl-mpm/fastmaths:fast-scale!
                                   (magicl:@ s (magicl:transpose s))
                                   (/ 1d0 (expt rho 3))))))
                           (setf (varef b 3) f)
                           (loop for i from 0 to 4
                                 while (or (> (sqrt
                                                (+
                                                 (expt (varef b 0) 2)
                                                 (expt (varef b 1) 2)
                                                 (expt (varef b 2) 2)))
                                               tol)
                                            (> (abs (varef b 3)) tolf)
                                            )
                                 do
                                    (let* ((A (magicl:block-matrix (list
                                                                    (magicl:.+
                                                                     (magicl:eye 3)
                                                                     (cl-mpm/fastmaths:fast-scale! (magicl:@ dgg De3) dgam)
                                                                     )
                                                                    dg
                                                                    (magicl:@ (magicl:transpose df) De3)
                                                                    (magicl:zeros '(1 1)))
                                                                   (list 2 2)
                                                ))
                                           (dx (cl-mpm/fastmaths:fast-scale! (magicl:linear-solve A b) -1d0)))
                                      ;; (pprint A)
                                      ;; (pprint dx)
                                      (loop for i from 0 to 2 do 
                                        (incf (varef epsE i) (varef dx i)))
                                      (incf dgam (varef dx 3))
                                      (setf sig (magicl:@ De3 epsE))
                                      (setf xi (/ (cl-mpm/fastmaths:fast-sum sig) (sqrt 3)))
                                      (setf s (magicl:.-
                                               sig
                                               (/ xi (sqrt 3))))
                                      (setf rho (sqrt (cl-mpm/fastmaths:dot s s)))
                                      (setf df (magicl:.-
                                                (cl-mpm/fastmaths:fast-scale-vector s (/ 1d0 rho))
                                                (/ alfa (sqrt 3))))
                                      (setf dg (magicl:.-
                                           (cl-mpm/fastmaths:fast-scale-vector s (/ 1d0 rho))
                                           (/ bta (sqrt 3))))
                                      (setf dgg
                                       (magicl:.-
                                        (cl-mpm/fastmaths:fast-scale!
                                         (magicl:.-
                                          (magicl:eye 3 :value 3)
                                          (magicl:ones (list 3 3))
                                          )
                                         (/ 1d0 (* 3d0 rho)))
                                        (cl-mpm/fastmaths:fast-scale!
                                         (magicl:@ s (magicl:transpose s))
                                         (/ 1d0 (expt rho 3)))))

                                      (loop for i from 0 to 2 do
                                        (setf (varef b i) (+ (- (varef epsE i) (varef epstr i))
                                                             (* dgam (varef dg i)))))
                                      ;; (pprint (varef b 3))
                                      (setf (varef b 3) (- rho (* alfa xi))))))
                         (loop for i from 0 to 2 do 
                           (incf (varef sig i) (/ xsic (sqrt 3))))
                         ))

                      (setf epsE (magicl:@ Ce sig))
                      (let ((pad-eps (magicl:block-matrix (list epsE
                                                                (cl-mpm/utils:vector-zeros))
                                                          '(2 1))))
                        (setf epsE (magicl:linear-solve Q pad-eps))
                        )

                      (setf sig (magicl:@ (magicl:transpose! Q)
                                          (cl-mpm/utils:voigt-from-list
                                           (list
                                            (magicl:tref sig 0 0)
                                            (magicl:tref sig 1 0)
                                            (magicl:tref sig 2 0)
                                            0d0 0d0 0d0))))
                      (values
                       ;; sig
                       (swizzle-coombs->voigt sig)
                       (swizzle-coombs->voigt epsE)
                       initial-f
                       t))
                    ;;No MC yield - just return
                    (values stress
                            trial-elastic-strain
                            initial-f
                            0d0
                            nil)))))
          ;;No DP yield - just return
          (values stress
                  trial-elastic-strain
                  f-dp
                  0d0
                  nil)))))


(defun maxwell-forwards-fs ()

  )

