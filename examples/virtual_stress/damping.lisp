
(defpackage :cl-mpm/examples/virtual-stress/damping
  (:use :cl))
(in-package :cl-mpm/examples/virtual-stress/damping)
;; (sb-ext:restrict-compiler-policy 'speed  3 3)
;; (sb-ext:restrict-compiler-policy 'debug  0 0)
;; (sb-ext:restrict-compiler-policy 'safety 0 0)
(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)

(declaim (optimize (debug 3) (safety 3) (speed 2)))



;; (defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-elastic) dt fbar)
;;   (cl-mpm::update-stress-linear mesh mp dt fbar))

(declaim (notinline plot))
(defun plot (sim &optional (plot :deformed))
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:axis (list 0 ms-x
                       0 ms-y))
    (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
  (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim))))
    (vgplot:format-plot t "set ytics ~f" h)
    (vgplot:format-plot t "set xtics ~f" h))
  (cl-mpm/plotter::simple-plot
   sim
   :plot :deformed
   ;:colour-func (lambda (mp) (compute-radial-stress mp))
   :colour-func (lambda (mp) (abs (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :xy)))
   ;; :contact-bcs *penalty-bc*
   ))

(declaim (notinline setup-test-column))
(defun setup (&key (refine 1d0) (mps 2d0))
  (let* ((L 2d0)
         (d 2d0)
         (dsize 10d0)
         (E 1d6)
         (density 900d0)
         (mesh-resolution (/ 1d0 refine))
         (domain-size (list dsize dsize))
         (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
         (depth 1d0)
         (offset-y 5d0)
         (block-size (list L d)))
    (defparameter *sim* (cl-mpm/setup::make-simple-sim
                 mesh-resolution
                 element-count
                 ;; :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-multi grid
                 ;; :sim-type 'cl-mpm/aggregate::mpm-sim-agg-usf
                 :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic
                 :args-list (list :enable-aggregate t
                                  ;; :ghost-factor (* E 1d-4)
                                  :vel-algo :FLIP)))
    (setf mesh-resolution (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0 offset-y)
      block-size
      (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
      density
      'cl-mpm/particle::particle-elastic
      :E E
      :nu 0.2d0
      :index 0))
    (setf (cl-mpm:sim-gravity *sim*) -10d0)
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)
    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 nil nil))
    (setf (cl-mpm:sim-dt *sim*)
          (* 0.5d0 (cl-mpm/setup:estimate-elastic-dt *sim*)))
    (setf *run-sim* t)
    (let ((datum (+ offset-y (* d 0.8d0)))
          (water-density 1000d0))
      (defparameter *water*
        (cl-mpm/buoyancy::make-bc-buoyancy-clip
         *sim*
         datum
         water-density
         (lambda (pos datum) t)
         :visc-damping 1d0)
        )
      (cl-mpm:add-bcs-force-list
       *sim*
       *water*)
      )
    )
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  (defparameter *sim-step* 0)
  )



(defparameter *run-sim* nil)

(defun stop ()
  (setf (cl-mpm::sim-run-sim *sim*) nil)
  (setf *run-sim* nil))

(defun plot-disp ()
  (vgplot:plot
   *data-time*
   *data-disp*)
  )

(defun get-disp ()
  (let ((mps (cl-mpm:sim-mps *sim*)))
    (/
     (cl-mpm::reduce-over-mps
      mps
      (lambda (mp)
        (cl-mpm/utils:get-vector (cl-mpm/particle::mp-velocity mp) :y))
      #'+)
     (length mps))))


(defun output-disp-header (output-dir)
  (with-open-file (stream (merge-pathnames "disp.csv" output-dir) :direction :output :if-exists :supersede)
    (format stream "time,disp~%")))

(defun output-disp-data (output-dir)
  (with-open-file (stream (merge-pathnames "disp.csv" output-dir) :direction :output :if-exists :append)
    (format stream "~f,~f~%"
            (cl-mpm::sim-time *sim*)
            (get-disp))))

(defun run (&key (output-dir "./output/"))
  (defparameter *data-time* '())
  (defparameter *data-disp* '())
  (vgplot:close-all-plots)
  (cl-mpm/dynamic-relaxation::run-time
   *sim*
   :output-dir output-dir
   :total-time 5d0
   :dt 0.1d0
   :dt-scale 5d0;0.25d0
   :initial-quasi-static nil
   :damping 1d-6
   :post-conv-step (lambda (sim)
                     (output-disp-header output-dir))
   :post-iter-step (lambda (sim)
                     (output-disp-data output-dir)
                     (push (cl-mpm::sim-time sim) *data-time*)
                     (push
                      (get-disp)
                       *data-disp*))
   :plotter (lambda (sim) (plot-disp))))

(defun test ()
  (let ((r 2)
        (mps 4))
    (setup :refine r :mps mps)
    (run :output-dir (format nil "./output-~D-~D/" r mps))))

(defun test-refine ()
  (dolist (r (list 2 3 4))
    (let ((mps 3))
      (setup :refine r :mps mps)
      (run :output-dir (format nil "./output-~D-~D/" r mps)))))
