(ql:quickload :cl-mpm)
(ql:quickload :cl-mpm/setup)
;;Create a simulation of 10x10 nodes 1x1 meter wide
(defparameter *sim*
  (cl-mpm/setup::make-block
   0.1d0 ;; Resolution h=h_x=h_y
   (list 10 10) ;; Elements 
   ))
;;Create a block of elastic material points 0.5x0.5 meters across with 10x10 material points total
(setf (cl-mpm:sim-mps *sim*)
      (cl-mpm/setup::make-block-mps
       (list 0d0 0d0);; Offset
       (list 0.5d0 0.5d0) ;; Size
       (list 20 20);; Mp count
       1000d0 ;Density 1kg/m^3
       ;; 'cl-mpm/particle::particle-elastic
       'cl-mpm/particle::particle-vm
       :E 1d5 ;; Young's modulus 1GPa
       :nu 0.35d0 ;; Poission's ratio
       :rho 1d4
       :gravity -9.8d0 ;;Gravity acting in the y direction
       ))

;;Apply a damping n/m^2
(setf (cl-mpm:sim-damping-factor *sim*) 10d0)
;;Estimate a dt value from the p wave modulus
(setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup:estimate-elastic-dt *sim*))
(format t "~&Estimated dt ~f~%" (cl-mpm:sim-dt *sim*))

(ql:quickload :cl-mpm/plotter)
(cl-mpm/plotter:simple-plot *sim* :plot :deformed)

;;Setup how many threads
(setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))

(time
 (dotimes (i 500)
   (format t "Steps ~D~%" i)
   (cl-mpm:update-sim *sim*)
   (format t "Velocity norm : ~F~%"
           (/
            (loop for mp across (cl-mpm:sim-mps *sim*)
                  sum (cl-mpm/fastmath::mag (cl-mpm/particle:mp-velocity mp)))
            (length (cl-mpm:sim-mps *sim*))))

   ))

;;Replot the xx stress
(cl-mpm/plotter:simple-plot
 *sim*
 :plot :deformed
 :colour-func (lambda (mp) (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0))
 )

(vgplot:print-plot
 (merge-pathnames (format nil "tutorial.png"))
 :terminal "png size 1920,1080")

