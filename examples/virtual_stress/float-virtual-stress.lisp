(defpackage :cl-mpm/examples/float-virtual-stress
  (:use :cl
   :cl-mpm/example
   :cl-mpm/utils)
  )
(in-package :cl-mpm/examples/float-virtual-stress)
(declaim (notinline plot-domain))
(defun plot-domain ()
  (when *sim*
    ;; (cl-mpm::g2p (cl-mpm:sim-mesh *sim*)
    ;;              (cl-mpm:sim-mps *sim*)
    ;;              (cl-mpm:sim-dt *sim*)
    ;;              0d0
    ;;              :QUASI-STATIC)
    (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
           (h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
           (ms-x (first ms))
           (ms-y (second ms)))
      (vgplot:format-plot t "set object 1 rect from 0,0 to ~f,~f fc rgb 'blue' fs transparent solid 0.5 noborder behind" ms-x *datum*)
      (cl-mpm/plotter:simple-plot
       *sim*
       :plot :deformed
       :trial t
       :colour-func (lambda (mp) (cl-mpm/particle::mp-index mp))))))

(defun setup (&key (refine 1) (mps 3))
  (let* ((L 10d0)
         (d 1d0)
         (domain-width 30d0)
         (h (/ 1d0 refine))
         (offset h)
         (height 10d0)
         (domain-height 20d0)
         (density 1000d0)
         (E 1d6)
         (domain-size (list domain-width domain-height))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list height height)))

    (setf *sim* (cl-mpm/setup::make-simple-sim
                 h
                 element-count
                 :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static
                 :args-list (list :enable-aggregate t
                                  :enable-fbar nil
                                  :enable-split t
                                  :vel-algo :QUASI-STATIC
                                  )))
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0 offset)
      block-size
      (mapcar (lambda (e) (* (/ e h) mps)) block-size)
      density
      'cl-mpm/particle::particle-elastic
      :E E
      :nu 0.2d0
      ))
    ;; (setf (cl-mpm::sim-ghost-factor *sim*) (* E 1d-4))
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)
    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 nil nil))

    (let ((datum (+ offset height))
          (water-density (* density 1.5d0))
          (pressure-condition t)
          )
      (defparameter *datum* datum)
      (defparameter *water-bc*
        (if pressure-condition
            (cl-mpm/buoyancy::make-bc-buoyancy-clip
             *sim*
             datum
             water-density
             (lambda (pos datum) t)
             :visc-damping 1d0)
            (cl-mpm/buoyancy::make-bc-buoyancy-body
             *sim*
             datum
             water-density
             (lambda (pos) t)))))

    (cl-mpm:add-bcs-force-list
     *sim*
     *water-bc*)

    (defparameter *current-inc* 0d0)
    (setf (cl-mpm:sim-dt *sim*)
          (* 0.5d0 (cl-mpm/setup:estimate-elastic-dt *sim*)))
    (setf *run-sim* t))
  (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
  )

(defun run ()
  (let* ((lstps 20)
         (total-disp -5d0)
         (delta (/ total-disp lstps)))
    (cl-mpm/dynamic-relaxation::run-load-control
     *sim*
     :output-dir (format nil "./output/")
     :plotter (lambda (sim) (plot-domain))
     :load-steps lstps
     :kinetic-damping nil
     :damping 1d0
     :substeps 50
     :criteria 1d-3
     :save-vtk-dr t
     :save-vtk-loadstep t
     :dt-scale 0.9d0
     )))

(defun test ()
  (setup :mps 3 :refine 2)
  (run))
