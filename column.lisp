(ql:quickload "cl-mpm")
(ql:quickload "cl-mpm/setup")
(ql:quickload "vgplot")
(ql:quickload "swank.live")

(setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))
(defparameter *sim* (cl-mpm/setup::make-column 16))

(defun plot (sim)
      (let* ( 
             (pos (loop for mp in (cl-mpm:sim-mps sim) collect (slot-value mp 'cl-mpm::position)))
             (sig (loop for mp in (cl-mpm:sim-mps sim) collect (slot-value mp 'cl-mpm::stress)))
             (x (loop for p in pos collect (magicl:tref p 0 0)))
             (sigy (loop for s in sig collect (magicl:tref s 1 0)))
             (y (loop for p in pos collect (magicl:tref p 1 0))))
        (vgplot:plot x y ";;with points pt 7")
        ) 
      (vgplot:replot)
  )

(defparameter *run-sim* nil)
(loop for steps from 0 to 100
      while *run-sim*
      do
      (progn
        (format t "Step ~d ~%" steps)
        (dotimes (i 100)
          (cl-mpm::update-sim *sim*))
        (plot *sim*)
        (swank.live:update-swank)
        (sleep .05)))
