# CL-MPM
 MPM/GIMP/B-Spline MPM implementation in common lisp, only SBCL conforming.
 Uses an explicit dynamic formulation with mass scaling, and adaptive timestepping.
 Used for simulating viscoelastic glacier dynamics.
 
 
# Features
 - GIMP shape functions
 - Extensible design using CLOS
 - Non-local continuum damage model
 - Nonconforming boundary conditions (virtual stress method)
 - Releasing arbrbitary penalty boundary conditions
 - Particle addition and deletion
 - Eigen errosion post-process step
 - MPI capabilities
 - Multithreaded with lparallel
 
Output
 - VTK mesh & mp output
 - Live GNUplot plotter
 - CSV format to conform with cb-geo

# Installation
cl-mpm only runs with SBCL due to its use of SIMD intrinsics, install sbcl from https://www.sbcl.org/
In the console we can start a lisp REPL with:
```sbcl --dynamic-space-size 4000```

Most of the dependancies are on quicklisp but some need to be installed manually:
After installing quicklisp (https://www.quicklisp.org/beta/#installation) you can install cl-mpm by cloneing the following repositories to ~/quicklisp/local-projects/
 - https://github.com/Ionsto/cl-mpm.git
 - https://github.com/Ionsto/magicl.git
 - https://github.com/Ionsto/vgplot.git

To include the accelerated c++ kirchoff strain update; run build.sh in the libs/ folder
This will require gcc & an install of eigen3.

Then to get started you can quickload cl-mpm and other packages
```
(ql:quickload :cl-mpm)
```

# Tutorial - quick start
Lets run a quasi-static elastic 2D problem where a square of material loaded under gravity squishes a little.
We start by loading the cl-mpm and setup package.
This can either be done in an interactive REPL, or
```
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
       'cl-mpm/particle::particle-elastic
       :E 1d5 ;; Young's modulus 1GPa
       :nu 0.35d0 ;; Poission's ratio
       :gravity -9.8d0 ;;Gravity acting in the y direction
       ))

;;Apply a damping n/m^2
(setf (cl-mpm:sim-damping-factor *sim*) 10d0)
;;Estimate a dt value from the p wave modulus
(setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup:estimate-elastic-dt *sim*))
;;Output this to the terminal
(format t "Estimated dt ~f~%" (cl-mpm:sim-dt *sim*))

```
You can then visulise the domain with gnuplot (if installed)
```
(ql:quickload :cl-mpm/plotter)
(cl-mpm/plotter:simple-plot *sim*)
```
We need to set a kernel count for the parallel implementation (I have a 4 core machine).
We can then run through 500 steps to let the system settle
```
;;Set the kernel thread count
(setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))
;;Do 500 updates
(dotimes (i 500)
  (format t "Steps ~D~%" i)
  (cl-mpm:update-sim *sim*))
```
We can then check the displacement of the system, using a lambda function to extract the x displacement from the (3x1) sized vector
```
(cl-mpm/plotter:simple-plot
 *sim*
 :plot :deformed
 :colour-func (lambda (mp) (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0)))
```
We get a displacement variation as expected:
![tutorial](https://github.com/Ionsto/cl-mpm/assets/117826225/95163c5a-35a8-4312-a71c-c9a4941ee388)

# MPI
To use mpi you need to compile an image worker - check the ham-mpi repo 
