(defpackage :cl-mpm/examples/virtual-stress/float
  (:use :cl
   :cl-mpm/example
   :cl-mpm/utils)
  )
(in-package :cl-mpm/examples/virtual-stress/float)

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-elastic) dt fbar)
  (cl-mpm::update-stress-kirchoff-dynamic-relaxation mesh mp dt fbar)
  )
(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-elastic) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  (cl-mpm::update-domain-polar-2d mesh mp dt)
  ;; (cl-mpm::update-domain-deformation mesh mp dt)
  )

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
       :colour-func (lambda (mp) (cl-mpm/utils:varef (cl-mpm/particle::mp-stress mp) 1))))))

(defun setup (&key (refine 1) (mps 3))
  (let* ((L 1d0)
         (d 1d0)
         (domain-width 5d0)
         (h (/ 0.25d0 refine))
         (offset l)
         (height 1d0)
         (domain-height (* 3 L))
         (density 1000d0)
         (E 1d7)
         (domain-size (list domain-width domain-height))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list height height)))

    (setf *sim* (cl-mpm/setup::make-simple-sim
                 h
                 element-count
                 ;:sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static
                 :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
                 :args-list (list :enable-aggregate t
                                  :enable-fbar nil
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
      :nu 0.2d0))

    (setf (cl-mpm::sim-gravity *sim*) -1d1)
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)
    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 nil nil))

    (let ((datum (+ offset (* 1 height)))
          (water-density (* density 2d0))
          (pressure-condition t))
      (defparameter *datum* datum)
      (defparameter *water-bc*
        (if pressure-condition
            (cl-mpm/buoyancy::make-bc-buoyancy-clip
             *sim*
             datum
             water-density
             (lambda (pos datum) t)
             :visc-damping 0d0)
            (cl-mpm/buoyancy::make-bc-buoyancy-body
             *sim*
             datum
             water-density
             (lambda (pos) t)))))
    (cl-mpm:add-bcs-force-list
     *sim*
     *water-bc*)
    (defparameter *current-inc* 0d0)
    ;; (setf
    ;;  (cl-mpm:sim-dt *sim*)
    ;;  (* 0.5d0 (cl-mpm/setup:estimate-elastic-dt *sim*)))
    (setf *run-sim* t))
  (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*))))

(defun get-height ()
  )

(defun run ()
  (let* ((lstps 1))
    (cl-mpm/dynamic-relaxation::run-load-control
     *sim*
     :output-dir (format nil "./output/")
     :plotter (lambda (sim) (plot-domain))
     :load-steps lstps
     :damping (sqrt 2d0)
     :substeps 100
     :criteria 1d-9
     :save-vtk-dr t
     :save-vtk-loadstep t
     :dt-scale 1d0
     )))

(defun run-explicit ()
  (change-class *sim* 'cl-mpm/aggregate::mpm-sim-agg-usf)
  (setf (cl-mpm::sim-velocity-algorithm *sim*) :TBLEND)
  (cl-mpm/dynamic-relaxation::run-time
   *sim*
   :output-dir (format nil "./output/")
   :plotter (lambda (sim) (plot-domain))
   :damping 1d-2
   :dt-scale 0.9d0
   :dt 0.1d0
   :total-time 100d0
   :initial-quasi-static nil
   ))

(defun test ()
  (cl-mpm/utils::set-workers 8)
  (setup :mps 3 :refine 2)
  (run-explicit)
  )
