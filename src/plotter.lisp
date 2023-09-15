(defpackage :cl-mpm/plotter
  (:use
   :cl
   :cl-mpm/particle
   :cl-mpm/utils
   :cl-mpm/mesh)
  (:import-from
   :magicl tref .+ .-
   )
  (:export
   simple-plot
   ))
(declaim (optimize (debug 0) (safety 0) (speed 3)))
                                        ;    #:make-shape-function
(in-package :cl-mpm/plotter)

(defun simple-plot (sim &key (plot :point) (colour-func (lambda (mp) 0d0)))
  (declare (function colour-func))
  "A simple GIMP plot that display only the position and size of the MPs in a sim"
  (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  (multiple-value-bind (x y lx ly c)
    (loop for mp across (cl-mpm:sim-mps sim)
          collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
          collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
          collect (magicl:tref (cl-mpm/particle::mp-domain-size mp) 0 0) into lx
          collect (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0) into ly
          collect (funcall colour-func mp) into c
          finally (return (values x y lx ly c)))
    (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 1d-30 (apply #'max c)))
    (cond
      ((eq plot :point)
       (vgplot:plot x y c ";;with points pt 7 lc palette"))
      ((eq plot :deformed)
       (vgplot:plot x y lx ly c ";;with ellipses lc palette"))))
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms)))
    (vgplot:axis (list 0 ms-x
                       0 ms-y))
    (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
    (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim))))
      (vgplot:format-plot t "set ytics ~f" h)
      (vgplot:format-plot t "set xtics ~f" h))
  (vgplot:replot))
