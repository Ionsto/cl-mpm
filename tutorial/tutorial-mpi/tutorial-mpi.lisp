(in-package :cl-mpm-worker)
;;Setup out local threads
(setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))

;;Work out which node this instance is
(defparameter *rank* (cl-mpi:mpi-comm-rank))

;;Check out output folder exists - but only do this on one node
(when (= *rank* 0)
  (uiop:ensure-all-directories-exist (list (uiop:merge-pathnames* "./output/"))))

;;Create a simulation of 10x10 nodes 1x1 meter wide
(defparameter *sim*
  (cl-mpm/setup::make-block
   0.1d0 ;; Resolution h=h_x=h_y
   (list 10 10) ;; Elements 
   :sim-type 'cl-mpm/mpi::mpm-sim-mpi-nodes
   ))

;;Create a full block of elastic material points 0.5x0.5 meters across with 10x10 material points total
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
       :rho 1.5d3
       :gravity -9.8d0 ;;Gravity acting in the y direction
       ))

;;Apply a damping n/m^2
(setf (cl-mpm:sim-damping-factor *sim*) 10d0)
;;Estimate a dt value from the p wave modulus
(setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup:estimate-elastic-dt *sim*))

;;Only output this information once
(when (= *rank* 0)
  (format t "~&Estimated dt ~f~%" (cl-mpm:sim-dt *sim*)))

;;Save the mesh out once
(when (= *rank* 0)
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk") *sim*))


;;The amount of nodes we've run with mpirun
(let ((node-count (cl-mpi:mpi-comm-size)))
  ;;Set the amount of nodes in the x-y-z axis
  ;;This setup will make each node have a slice along the x axis
  ;;n0-n1-n2-n3
  (setf (cl-mpm/mpi::mpm-sim-mpi-domain-count *sim*) (list node-count 1 1))
  ;;If the domain count is (2 2 1) the problem will be decomposed like:
  ;; n2-n3
  ;; n0-n1
  )

;;Output the amount of mps that the total problem has
(when (= *rank* 0)
  (format t "Sim MPs: ~a~%" (length (cl-mpm:sim-mps *sim*)))
  (format t "Decompose~%"))

;;Split up our domains
(cl-mpm/mpi:domain-decompose *sim*)

;;Check how many MPs per node we have after decomposition
(format t "Rank ~D - Sim MPs: ~a~%" *rank* (length (cl-mpm:sim-mps *sim*)))

;;Run through 100 steps - outputting the mps at each step
(time
 (dotimes (i 100)
   (format t "Steps ~D~%" i)
   ;;Normal update sim step
   (cl-mpm:update-sim *sim*)
   ;;Aim to estimate a velocity norm of sorts however we need to reduce this over all nodes
   (let* ((nodal-norm (loop for mp across (cl-mpm:sim-mps *sim*)
                            sum (cl-mpm/fastmath::mag (cl-mpm/particle:mp-velocity mp))))
          (overall-average (cl-mpm/mpi:mpi-average nodal-norm (length (cl-mpm:sim-mps *sim*)))))
     ;;Output our reduced norm
     (when (= *rank* 0)
       (format t "Velocity norm : ~E~%" overall-average)))

   ;;Output vtks for every node at every timestep
   (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~2,'0d_~5,'0d.vtk" *rank* i)) *sim*)
   ))

;;Wait for everybody to finish
(cl-mpi:mpi-waitall)
;;Tear down threads
(lparallel:end-kernel)
;;Return back to the mpi-worker image where it will close everything for us
