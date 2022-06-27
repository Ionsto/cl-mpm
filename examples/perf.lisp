(declaim (optimize (debug 0) (safety 0) (speed 3)))
(ql:quickload "cl-mpm")
(ql:quickload "cl-mpm/setup")
(ql:quickload "vgplot")
(ql:quickload "swank.live")
(ql:quickload "cl-mpm/output")
(ql:quickload "magicl")

(defmacro time-form-noprint (form it)
  `(progn
    ;(declaim (optimize speed))
    (let* ((iterations ,it)
           (start (get-internal-real-time)))
      (dotimes (i iterations) ,form)
      (let* ((end (get-internal-real-time))
             (units internal-time-units-per-second)
             )
        (/ (- end start) (* iterations units))))))

(defun make-test-column (size)
  (let* ((sim (cl-mpm/setup::make-column 1 size)))
    (let* ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
           (e-scale 1)
           (h-x h)
           (h-y (/ h e-scale))
           (elements (* e-scale (- size 1)))
           )
      (setf (cl-mpm:sim-mps sim)
            (cl-mpm/setup::make-column-mps-elastic
              elements
              (list h-x h-y)
              1e4 0d0))) 
    (setf (cl-mpm:sim-damping-factor sim) 0)
    (setf (cl-mpm:sim-mass-filter sim) 0)
    sim
    ))


(defun test-throughput ()
  (progn
    (let* ((sizes (loop for i from 1 to 12 collect (expt 2 i)))
           (times (loop for size in sizes collect 
                        (let ((sim (make-test-column size)))
                          (time-form-noprint (cl-mpm::update-sim sim) 10))))
           (throughput (loop for dt in times collect (/ 1 dt)))
           )
      (values sizes times throughput)
      )))

(defun test-plot-throughput ()
  (multiple-value-bind (sizes times throughput) (test-throughput)
      (vgplot:close-all-plots)
      (vgplot:figure)
      (vgplot:xlabel "Size")
      (vgplot:ylabel "Throughput")
      (vgplot:plot sizes throughput)
      (vgplot:figure)
      (vgplot:loglog sizes times)
      (vgplot:xlabel "Size")
      (vgplot:ylabel "Time per step")))

(progn
  (defparameter *perf-sizes* '())
  (defparameter *perf-times* '())
  (loop for cores in '(1 4 8)
        do (progn
             (setf lparallel:*kernel* (lparallel:make-kernel cores :name "custom-kernel"))
             (multiple-value-bind (sizes times throughput) (test-throughput)
               (push sizes *perf-sizes*)
               (push times *perf-times*))))

  (vgplot:close-all-plots)
  (vgplot:figure)
  (vgplot:xlabel "Size")
  (vgplot:ylabel "Throughput")
  (vgplot:loglog    (nth 0 *perf-sizes*) (nth 0 *perf-times*) "0"
                    (nth 1 *perf-sizes*) (nth 1 *perf-times*) "4"
                    (nth 2 *perf-sizes*) (nth 2 *perf-times*) "8"
                 )
  (vgplot:legend)
  )
;(setf lparallel:*kernel* (lparallel:make-kernel 1 :name "custom-kernel"))
;(test-throughput)
;(vgplot:figure)
