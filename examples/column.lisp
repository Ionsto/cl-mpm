(declaim (optimize (debug 0) (safety 0) (speed 3)))
(ql:quickload "cl-mpm")
(ql:quickload "cl-mpm/setup")
(ql:quickload "vgplot")
(ql:quickload "swank.live")
(ql:quickload "cl-mpm/output")
(ql:quickload "magicl")


(defun plot (sim)
  (multiple-value-bind (x y)
    (loop for mp in (cl-mpm:sim-mps sim)
          collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
          collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
          finally (return (values x y)))
    (vgplot:plot x y ";;with points pt 7"))
  (vgplot:replot))

(setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))

(defparameter *run-sim* nil)
(defparameter *run-sim* t)

(defun setup-test-column (size)
  (let* ((sim (cl-mpm/setup::make-column 1 size)) 
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (e-scale 1)
         (h-x h)
         (h-y (/ h e-scale))
         (elements (* e-scale (- size 1))))
    (progn
      (setf (cl-mpm:sim-mps sim) 
            (cl-mpm/setup::make-column-mps-elastic
              elements
              (list h-x h-y)
              1e4 0d0))
      (setf (cl-mpm:sim-damping-factor sim) 0)
      (setf (cl-mpm:sim-mass-filter sim) 0)
      sim)))
(progn
  (defparameter *sim* (setup-test-column 4))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  )

(progn 
    (vgplot:close-all-plots)
    (vgplot:figure)
    (vgplot:axis (list 0 (nth 0 (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*))) 
                       0 (nth 1 (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))))
    (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
      (vgplot:format-plot t "set ytics ~f" h)
      (vgplot:format-plot t "set xtics ~f" h))
    (time (loop for steps from 0 to 10
                while *run-sim*
                do
                (progn
                  (format t "Step ~d ~%" steps)
                  (dotimes (i 100)
                    (cl-mpm::update-sim *sim*)
                    (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))
                    (let ((h (/ (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)) 2)))
                      (setf *velocity* (cons (magicl:tref (cl-mpm/output::sample-point-velocity *sim* (list h (* h 2))) 1 0) *velocity*)))
                    (setf *time*     (cons *t* *time*))
                    )
                  (plot *sim*)
                  (swank.live:update-swank)
                  (sleep .10)
                  )))
    (vgplot:figure)
    (vgplot:title "Velocity over time")
    (vgplot:plot *time* *velocity*)
    )
