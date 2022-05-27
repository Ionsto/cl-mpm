(ql:quickload "cl-mpm")
(ql:quickload "vgplot")
(ql:quickload "swank.live")
(ql:quickload "magicl")


(defparameter *nD* 2)

(defun dostep (sim n)
  (progn
    (dotimes (i n)
      (cl-mpm::update-sim sim)
      )
      (when  (cl-mpm:sim-mps sim)
      (let* ( 
             (pos (loop for mp in (cl-mpm:sim-mps sim) collect (slot-value mp 'cl-mpm::position)))
             (sig (loop for mp in (cl-mpm:sim-mps sim) collect (slot-value mp 'cl-mpm::stress)))
             (x (loop for p in pos collect (magicl:tref p 0 0)))
             (sigy (loop for s in sig collect (magicl:tref s 1 0)))
             (y (loop for p in pos collect (magicl:tref p 1 0))))
        (vgplot:plot x y ";;with points pt 7")
        ))
      (vgplot:replot)
      (swank.live:update-swank)
      (sleep .05)
      ))
(progn 
  (defparameter *mp-spacing* 1d0)
  (defparameter *sim* (cl-mpm:make-mpm-sim '(1 15) 1 1e-3 
                                           (cl-mpm::make-shape-function-linear *nD*)))
  (setf (cl-mpm:sim-mps *sim*) 
        (list (cl-mpm::make-particle-elastic *nD* 1e2 0 :pos '(0.5 0.5))))
  (setf (cl-mpm:sim-mps *sim*) 
        (loop for i from 0 to 5 collect 
              (cl-mpm::make-particle-elastic *nD*
                                             1e2
                                             0
                                             :pos (list 0.5 (+ 0.5d0 (* *mp-spacing* i)))
                                             :volume *mp-spacing*)))
  (setf (cl-mpm:sim-bcs *sim*) (cl-mpm/bc:make-outside-bc (cl-mpm:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))))
(vgplot:figure)
(vgplot:axis (list 0 (nth 0 (cl-mpm:mesh-mesh-size (cl-mpm:sim-mesh *sim*))) 
                   0 (nth 1 (cl-mpm:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))))
(inspect (cl-mpm:sim-mps *sim*))
(dostep *sim* 1)
(loop for mp in (cl-mpm:sim-mps *sim*) do
      (setf (cl-mpm::mp-E mp) 1e1))
(dotimes (k 10)
  (format t "Frame ~d ~%" k)
  (dostep *sim* 50))
