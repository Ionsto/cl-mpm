(defpackage :cl-mpm/examples/sliding
  (:use :cl
   :cl-mpm/example
        :cl-mpm/utils))
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(in-package :cl-mpm/examples/sliding)
(declaim (optimize (debug 3) (safety 3) (speed 0)))


(defun setup (&key (refine 1) (mps 3)
                (angle 0d0)
                (mu 0.25d0)
                )
  (let* ((d 1d0)
         (dsize 50d0)
         (density 1d3)
         (mesh-resolution (/ 0.25d0 refine))
         (offset (* 2d0 mesh-resolution))
         (domain-size (list dsize (+ offset (* 2 d))))
         (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
         (block-size (list d d))
         (E 1d6))
    (setf *sim* (cl-mpm/setup::make-simple-sim
                 mesh-resolution
                 element-count
                 :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
                 :args-list (list :enable-aggregate t
                                  :enable-fbar nil)))
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list offset offset)
      block-size
      (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
      density
      'cl-mpm/particle::particle-elastic
      :E E
      :nu 0.2d0
      :index 0
      ;; :gravity-axis
      ;; 0d0
      ;; (cl-mpm/utils:vector-from-list
                    ;;  (list
                    ;;   (sin (* angle (/ pi 180d0)))
                    ;;   (cos (* angle (/ pi 180d0)))
                    ;;   0d0))
      ))

    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)

    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 0 nil))

    (let ((eps-scale 1d1))
      (defparameter *loader*
        (cl-mpm/penalty::make-bc-penalty-distance-point
         *sim*
         (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0))
         (cl-mpm/utils:vector-from-list (list offset (+ offset (* d 0.5d0)) 0d0))
         (* d 0.5d0)
         (* E 1d2)
         0d0
         0d0
         ))
      (defparameter *follower*
        (cl-mpm/penalty::make-bc-penalty-distance-point
         *sim*
         (cl-mpm/utils:vector-from-list (list -1d0 0d0 0d0))
         (cl-mpm/utils:vector-from-list (list (+ d offset) (+ offset (* d 0.5d0)) 0d0))
         (* d 0.5d0)
         (* E 1d2)
         0d0
         0d0
         ))
      (defparameter *floor*
        (cl-mpm/penalty::make-bc-penalty-distance-point 
         *sim*
         (cl-mpm/utils:vector-from-list (list 0d0 1d0 0d0))
         (cl-mpm/utils:vector-from-list (list (* dsize 0.5 ) offset 0d0))
         (* dsize 0.5d0)
         (* E 1d2)
         mu
         0d0
         )))
    (cl-mpm::add-bcs-force-list
     *sim*
     *floor*)
    (cl-mpm::add-bcs-force-list
     *sim*
     *loader*)
    (cl-mpm::add-bcs-force-list
     *sim*
     *follower*)
    (setf *run-sim* t))
  (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))

  )

(defun plot-vel ()
  (let ((times (reverse *data-time*))
        (vel (reverse *data-vel*))
        )
  (when (> (length *data-time*) 1)
    (let ((dt (- (second times) (first times))))
      (let ((acc (mapcar (lambda (x) (/ x dt)) (mapcar #'- (rest vel) (butlast vel)))))
        (vgplot:plot
         ;; times
         ;; vel
         (rest times)
         acc
         "Velocity"))))))

(defun run-time ()
  ;; (change-class *sim* 'cl-mpm/aggregate::mpm-sim-agg-usf)
  (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic)
  (setf (cl-mpm::sim-velocity-algorithm *sim*) :FLIP)
  (defparameter *data-time* (list 0d0))
  (defparameter *data-vel* (list 0d0))
  (let ((grav-angle (cl-mpm/particle::mp-gravity-axis (aref (cl-mpm::sim-mps *sim*) 0))))
    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (setf (cl-mpm/particle::mp-gravity-axis mp) (cl-mpm/utils:vector-from-list (list 0d0 1d0 0d0)))))
    (cl-mpm/dynamic-relaxation::run-time
     *sim*
     :output-dir (format nil "./output/")
     :plotter (lambda (sim) (plot-vel))
     :damping 0d-4
     :dt 0.01d0
     :total-time 10d0
     ;; :dt-scale 1d0
     :dt-scale 1d0
     :post-iter-step
     (lambda (sim)
       (push (cl-mpm::sim-time *sim*) *data-time*)
       (push (cl-mpm::reduce-over-mps
              (cl-mpm:sim-mps *sim*)
              (lambda (mp)
                (* 
                 (cl-mpm/particle::mp-volume mp)
                 (cl-mpm/utils::varef (cl-mpm/particle::mp-velocity mp) 0)))
              #'+)
             *data-vel*))
     :post-conv-step
     (lambda (sim)
       (cl-mpm::iterate-over-mps
        (cl-mpm:sim-mps *sim*)
        (lambda (mp)
          (setf (cl-mpm/particle::mp-gravity-axis mp) grav-angle)))))))

(defun run ()
  (let ((total-disp 1d-1)
        (current-disp 0d0))
    (defparameter *data-disp* (list 0d0))
    (defparameter *data-load* (list 0d0))
    (vgplot:close-all-plots)
    (cl-mpm/dynamic-relaxation::run-load-control
     *sim*
     :output-dir (format nil "./output/")
     :plotter (lambda (sim)
                (vgplot:plot
                 *data-disp*
                 *data-load*
                 "Load")
                ;; (plot-domain)
                )
     :load-steps 10
     :kinetic-damping nil
     :damping 1d0
     :substeps 50
     :criteria 1d-3
     :save-vtk-dr t
     :save-vtk-loadstep t
     :loading-function (lambda (percent)
                         (setf current-disp (* total-disp percent))
                         (cl-mpm/penalty::bc-set-displacement *loader*
                                                              (cl-mpm/utils:vector-from-list
                                                               (list
                                                                current-disp
                                                                0d0
                                                                0d0)))
                         (cl-mpm/penalty::bc-set-displacement *follower*
                                                              (cl-mpm/utils:vector-from-list
                                                               (list
                                                                current-disp
                                                                0d0
                                                                0d0))))
     :post-conv-step (lambda (sim)
                       (push current-disp *data-disp*)
                       (push (cl-mpm/penalty::resolve-load *loader*) *data-load*))
     :dt-scale 1d0))
  )

(defun test ()
  (let* ((angle -30d0)
         (mu 0.5d0))
    (setup :refine 2
           :mps 3
           :mu mu)
    (setf (cl-mpm::sim-gravity *sim*) -10d0)
    (format t "Maximum force ~E~%" (* mu 1d3))
    ;; (angle-calc angle mu)
    (run)
    ))

(defun angle-calc (angle mu)
  (let* ((psi (* (- 90d0 (abs angle)) (/ pi 180)))
         (gravity 10d0))
    (format t "~%Estimated acceleration ~F~%" (* gravity (max 0d0 (- (cos psi) (* mu (sin psi)))) ))))
