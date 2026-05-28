(defpackage :cl-mpm/examples/dr/collapse-mass-matrix
  (:use :cl
   :cl-mpm/example
   :cl-mpm/utils)
  )
(in-package :cl-mpm/examples/dr/collapse-mass-matrix)


;;Can be :K-0 :K-UPDATED :P-ELASTIC P-ELASTO-PLASTIC

(defparameter *mass-matrix* :K-0)

(defclass mpm-sim-test (cl-mpm/dynamic-relaxation::mpm-sim-quasi-static)
  ((mass-matrix-init
    :initform nil
    :accessor mpm-sim-mass-matrix-init)))

(defmethod cl-mpm/dynamic-relaxation::pre-step :before ((sim mpm-sim-test))
  (setf (mpm-sim-mass-matrix-init sim) nil))

(defun k-0-mass (sim))

(defmethod cl-mpm/dynamic-relaxation::update-node-fictious-mass ((sim mpm-sim-test))
  (case *mass-matrix*
    (:K-0
     ;;Do nothing on inner solve
     (unless (mpm-sim-mass-matrix-init sim)
       (cl-mpm/dynamic-relaxation::implicit-assemble-stiffness sim)
       (cl-mpm/aggregate::update-mass-matrix sim)))
    (:K-UPDATED
     ;;Do nothing
     (cl-mpm/dynamic-relaxation::implicit-assemble-stiffness sim)
     (cl-mpm/aggregate::update-mass-matrix sim)
     )
    (:P-ELASTIC
     ;;Do nothing
     (cl-mpm/dynamic-relaxation::map-stiffness-quasi-static sim)
     (cl-mpm/aggregate::update-mass-matrix sim)
     )
    (:P-ELASTO-PLASTIC
     ;;Do nothing
     (cl-mpm/dynamic-relaxation::map-stiffness-quasi-static sim)
     (cl-mpm/aggregate::update-mass-matrix sim))))

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-elastic) dt)
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
                 :args-list
                 (list
                  :enable-fbar t
                  :enable-aggregate t
                  :ghost-factor nil
                  :mass-update-count 1
                  :damping-update-count 1
                  :max-split-depth 6
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
          'cl-mpm/particle::particle-vm
          :E E
          :nu 0.3d0
          :rho 20d3
          ))
        (cl-mpm::domain-sort-mps sim)
        (defparameter *density* density)
        (cl-mpm/setup::set-mass-filter sim density :proportion 1d-15))))
  (cl-mpm/setup::setup-bcs
   *sim*
   :left (list 0 nil nil)
   :bottom (list nil 0 nil)
   :top (list nil nil nil)
   )
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defun run (&key (output-dir "./output/") )
  (let* ((lstps 10))
    (cl-mpm/dynamic-relaxation::run-load-control
     *sim*
     :output-dir output-dir 
     :plotter (lambda (sim) (plot-domain))
     :load-steps lstps
     :damping (sqrt 2d0)
     :substeps 50
     :conv-steps 1000
     :criteria 1d-9
     :save-vtk-dr nil
     :save-vtk-loadstep t
     :dt-scale 0.9d0)))

(defun test ()
  (cl-mpm/utils:set-workers 8)
  (setup :mps 3)
  (run)
  )
