(defpackage :cl-mpm/examples/dr/collapse-mass-matrix
  (:use :cl
   :cl-mpm/example
   :cl-mpm/utils)
  )
(in-package :cl-mpm/examples/dr/collapse-mass-matrix)

;; (setf cl-mpm/settings:*optimise-setting* cl-mpm/settings::*optimise-speed*)
;; (sb-ext:restrict-compiler-policy 'speed  3 3)
;; (sb-ext:restrict-compiler-policy 'debug  0 0)
;; (sb-ext:restrict-compiler-policy 'safety 0 0)
(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)
(setf cl-mpm/settings:*optimise-setting* cl-mpm/settings::*optimise-debug*)


;;Can be :K-0 :K-UPDATED :P-ELASTIC P-PLASTIC-SCALAR P-PLASTIC-TANGENT

(defparameter *mass-matrix-types* (list :K-UPDATED :K-0 :P-ELASTIC :P-PLASTIC-SCALAR :P-PLASTIC-TANGENT))
(defvar *mass-matrix* :K-UPDATED)

(defclass mpm-sim-test (cl-mpm/dynamic-relaxation::mpm-sim-dr-ul)
  ((mass-matrix-init
    :initform nil
    :accessor mpm-sim-mass-matrix-init)))

(defclass particle-vm-adjusted-pwaves (cl-mpm/particle::particle-vm)
  ())

(defmethod cl-mpm/particle::constitutive-model ((mp particle-vm-adjusted-pwaves) strain-in dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-accessors ((de cl-mpm/particle::mp-elastic-matrix)
                   (dep cl-mpm/particle::mp-tangent-stiffness)
                   (stress cl-mpm/particle::mp-stress-kirchoff)
                   (rho cl-mpm/particle::mp-rho)
                   (plastic-strain cl-mpm/particle::mp-strain-plastic)
                   (ps-vm  cl-mpm/particle::mp-strain-plastic-vm)
                   (ps-vm-1  cl-mpm/particle::mp-strain-plastic-vm-1)
                   (ps-vm-inc  cl-mpm/particle::mp-strain-plastic-vm-inc)
                   (strain  cl-mpm/particle::mp-strain)
                   (yield-func cl-mpm/particle::mp-yield-func)
                   (p-mod cl-mpm/particle::mp-p-modulus-0)
                   (E cl-mpm/particle::mp-e)
                   (nu cl-mpm/particle::mp-nu)
                   (enable-plasticity cl-mpm/particle::mp-enable-plasticity))
      mp
    ;;Train elastic strain - plus trail kirchoff stress
    (cl-mpm/constitutive::linear-elastic-mat strain de stress)
    (if (equal dep de)
      (setf dep (cl-mpm/utils::deep-copy de))
      (cl-mpm/utils::copy-into de dep))
    (setf p-mod (cl-mpm/particle::compute-p-modulus mp))
    (when enable-plasticity
      (multiple-value-bind (sig eps-e f inc pmod) (cl-mpm/ext::constitutive-vm-tangent
                                                   stress
                                                   strain
                                                   de
                                                   e nu
                                                   rho
                                                   dep)
        (setf stress sig
              yield-func f)
        (setf ps-vm-inc inc)
        (setf strain eps-e)
        (setf ps-vm (+ ps-vm-1 ps-vm-inc))
        (case *mass-matrix*
          (:P-ELASTIC
           ;;Done
           (setf p-mod (cl-mpm/particle::compute-p-modulus mp)))
          (:P-PLASTIC-SCALAR
           (when (> inc 0d0)
             (setf p-mod (cl-mpm/utils::calculate-bulk-modulus e nu))))
          (:P-PLASTIC-TANGENT
           (let ((p-mod-elastic (* 1d-9 (cl-mpm/particle::compute-p-modulus mp))))
             (setf p-mod p-mod-elastic)
             (loop for d from 0 below 2
                   do (let ((probe-vec (cl-mpm/utils::voigt-zeros)))
                        (setf (varef probe-vec d) 1d0)
                        (setf
                         p-mod
                         (max
                          p-mod
                          (cl-mpm/fastmaths:dot
                           probe-vec
                           (cl-mpm/fastmaths::fast-@-tensor-voigt dep probe-vec)))))))))
        )
      )
    stress))

(defmethod cl-mpm/dynamic-relaxation::pre-step :before ((sim mpm-sim-test))
  (setf (mpm-sim-mass-matrix-init sim) nil))


(defmethod cl-mpm/dynamic-relaxation::update-node-fictious-mass ((sim mpm-sim-test))
  (case *mass-matrix*
    (:K-0
     ;;Do nothing on inner solve
     (unless (mpm-sim-mass-matrix-init sim)
       (setf (mpm-sim-mass-matrix-init sim) t)
       (cl-mpm/dynamic-relaxation::implicit-assemble-stiffness sim)
       (cl-mpm/aggregate::update-mass-matrix sim)))
    (:K-UPDATED
     ;;Do nothing
     (cl-mpm/dynamic-relaxation::implicit-assemble-stiffness sim)
     (cl-mpm/aggregate::update-mass-matrix sim)
     )
    ;; (:P-PLASTIC-SCALAR
    ;;  ;;Do nothing
    ;;  (cl-mpm/dynamic-relaxation::map-stiffness-quasi-static sim)
    ;;  (cl-mpm/aggregate::update-mass-matrix sim)
    ;;  )
    ;; (:P-PLASTIC-TANGENT
    ;;  ;;Do nothing
    ;;  ;; (pprint "elasto plastic")
    ;;  (cl-mpm/dynamic-relaxation::map-stiffness-quasi-static sim)
    ;;  (cl-mpm/aggregate::update-mass-matrix sim))
    (t
     ;;Do nothing
     (cl-mpm/dynamic-relaxation::map-stiffness-quasi-static sim)
     (cl-mpm/aggregate::update-mass-matrix sim)
     ))
  )

(defmethod cl-mpm::update-particle (mesh (mp particle-vm-adjusted-pwaves) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  (cl-mpm::update-domain-stretch mesh mp dt))

(defun setup (&key (refine 1) (mps 2))
  (defparameter *sim* nil)
  (let ((mps-per-dim mps)
        (size '(16 16))
        (block-size '(8 8))
        )
    ;; (setf *sim* (setup-test-column  sim-type refine mps-per-dim multigrid-refine))
    (let* ((E 1d6)
           (density 1d3)
           (sim (cl-mpm/setup::make-simple-sim
                 (/ 1d0 refine)
                 (mapcar (lambda (x) (* x refine)) size)
                 :sim-type 'mpm-sim-test
                 ;; :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static
                 :args-list
                 (list
                  :enable-fbar nil
                  :enable-aggregate t
                  :ghost-factor nil
                  :mass-update-count 1
                  :damping-update-count 1
                  :max-split-depth 6
                  :mp-removal-size nil
                  ;; :split-factor nil
                  :enable-split nil
                  :gravity -10d0)))
           (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim))))
      (declare (double-float h density))
      (setf *sim* sim)
      (progn
        (cl-mpm::add-mps
         sim
         (cl-mpm/setup::make-block-mps
          (list 0d0 0d0 0d0)
          block-size
          (mapcar (lambda (e) (* (/ e h) mps)) block-size)
          density
          'particle-vm-adjusted-pwaves
          :E E
          :nu 0.3d0
          :rho 20d3))
        (cl-mpm::domain-sort-mps sim)
        (defparameter *density* density)
        (cl-mpm/setup::set-mass-filter sim density :proportion 1d-9))))
  (cl-mpm/setup::setup-bcs
   *sim*
   :left (list 0 nil nil)
   :bottom (list nil 0 nil)
   :top (list nil nil nil))
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defun run (&key (output-dir "./output/")
              (damping-scale (sqrt 2d0))
              (dt-scale 1d0)
              )
  (let* ((lstps 10))
    (assert (member *mass-matrix* *mass-matrix-types*))
    (cl-mpm/dynamic-relaxation::run-load-control
     *sim*
     :output-dir output-dir
     :plotter (lambda (sim) (plot-domain))
     :load-steps lstps
     :damping damping-scale
     :substeps 20
     :conv-steps 1000
     :criteria 1d-9
     :save-vtk-dr nil
     :save-vtk-loadstep nil
     :dt-scale dt-scale)))

(defun test ()
  (setup :refine 1 :mps 6)
  ;; (setf *mass-matrix* :P-PLASTIC-SCALAR)
  (run)
  )

(defun full-test ()
  (cl-mpm/utils:set-workers 8)
  (dolist (mm (list
               ;; :P-ELASTIC
               :P-PLASTIC-SCALAR
               :K-0
               :K-UPDATED
               ;; :P-PLASTIC-TANGENT
               ))
    (format t "Testing ~A~%" *mass-matrix*)
    (setup :refine 1d0 :mps 3)
    (setf *mass-matrix* mm)
    ;; (ignore-errors)
    (run :output-dir (format nil "./output-~a/" mm))))

(defun damping-test ()
  (cl-mpm/utils:set-workers 8)
  (dolist (mm (list
               ;; :VISC
               ;; :VISC-ADJUST-0.8
               :VISC-UNDER-0.8
               ;; :VISC-OVER
               ;; :KINETIC
               ))
    (setf *mass-matrix* :P-PLASTIC-TANGENT)
    ;; (format t "Testing ~A~%" *mass-matrix*)
    (setup :refine 1d0 :mps 6)
    (let ((damping-scale 0d0)
          (dt-scale 1d0)
          )
      (case mm
        (:VISC
         (setf damping-scale (sqrt 2d0)))
        ;; (:VISC-ADJUST-0.8
        ;;  (setf damping-scale (* (sqrt 2) 0.75d0)))
        (:VISC-UNDER-0.8
         (setf damping-scale (* (sqrt 2d0) 0.8d0)))
        ;; (:VISC-OVER
        ;;  (setf damping-scale (sqrt 4d0)))
        (:KINETIC
         (setf dt-scale 0.5d0)
         (setf (cl-mpm/dynamic-relaxation::sim-kinetic-damping *sim*) t)))
      ;; (ignore-errors)
      (run :output-dir (format nil "./output-~a/" mm)
           :damping-scale (* damping-scale)
           :dt-scale dt-scale
           ))))

(defmacro time-form (it form)
  `(progn
     (declaim (optimize speed))
     (let* ((iterations ,it)
            (start (get-internal-real-time))
            (gc-start sb-ext:*gc-real-time*)
            )
       (time
        (dotimes (i ,it)
          ,form))
       (let* ((end (get-internal-real-time))
              (gc-end sb-ext:*gc-real-time*)
              (units internal-time-units-per-second)
              (dt (/ (- end start) (* iterations units))))
         (format t "Total time: ~f ~%" (/ (- end start) units))
         (format t "Time per iteration: ~f~%" (/ (- end start) (* iterations units)))
         (format t "Total gc time: ~f ~%" (/ (- gc-end gc-start) units))
         (format t "Throughput: ~f~%" (/ 1 dt))
         (format t "Time per MP: ~E~%" (/ dt (length (cl-mpm:sim-mps *sim*))))
         (format t "MP Throughput: ~E~%" (/ (length (cl-mpm:sim-mps *sim*)) dt))
         dt))))
