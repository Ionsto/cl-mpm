# CL-MPM
 Fast and extensible material point method (MPM) code implementation in common lisp, only SBCL conforming.
 Uses an explicit dynamic formulation with mass scaling, and adaptive timestepping.
 Uses threading and MPI for large scale simulations.
 Nonlocal intergral continuum damage working within the MPI framework for effective damage simulations.
 An extensible constitutive model system for easy development.
 
 Used for simulating viscoelastic glacier dynamics.

Tested on Ubuntu & WSL.
Theoretically should work on windows however MAGICL library is required which is broken as of right now.
 
 
# Features
 - GIMP shape functions
 - Extensible design using CLOS
 - Non-local continuum damage model for quasi-static and dynamic damage
 - Nonconforming boundary conditions (virtual stress method)
 - Releasing arbrbitary penalty boundary conditions
 - Particle addition and deletion
 - MPI capabilities
 - Multithreaded with lparallel
 - Common plastic models: Von-Mises, Mohr-Coloumb
 
Output
 - VTK mesh & mp output
 - Live GNUplot plotter
 - CSV format to conform with cb-geo

# Installation
cl-mpm only runs with SBCL due to its use of SIMD intrinsics, install sbcl from https://www.sbcl.org/

```sudo apt install sbcl```

In the console we can start a lisp REPL with:
```sbcl --dynamic-space-size 4000```

Most of the dependancies are on quicklisp but some need to be installed manually:
After installing quicklisp (https://www.quicklisp.org/beta/#installation) you can install cl-mpm by cloning the following repositories to ~/quicklisp/local-projects/
 - https://github.com/Ionsto/cl-mpm.git
 - https://github.com/Ionsto/magicl.git
 - https://github.com/Ionsto/vgplot.git

To include the accelerated c++ kirchoff strain update; run build.sh in the libs/ folder
This will require gcc & an install of eigen3.

```sudo apt install build-essential libeigen3-dev```

Then to get started you can quickload cl-mpm and other packages
```
(ql:quickload :cl-mpm)
```

For live plotting gnuplot is used:
```sudo apt install gnuplot```

# Tutorial - quick start
Lets run a quasi-static elastic 2D problem where a square of material loaded under gravity squishes a little.
We start by loading the cl-mpm and setup package.
This can either be done in an interactive REPL, or putting this in a file and loading it with sbcl --load.
```lisp
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
```lisp
(ql:quickload :cl-mpm/plotter)
(cl-mpm/plotter:simple-plot *sim*)
```
We need to set a kernel count for the parallel implementation (I have a 4 core machine).
We can then run through 500 steps to let the system settle
```lisp
;;Set the kernel thread count
(setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))
;;Do 500 updates
(dotimes (i 500)
  (format t "Steps ~D~%" i)
  (cl-mpm:update-sim *sim*))
```
We can then check the displacement of the system, using a lambda function to extract the x displacement from the (3x1) sized vector
```lisp
(cl-mpm/plotter:simple-plot
 *sim*
 :plot :deformed
 :colour-func (lambda (mp) (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0)))
```
We get a displacement variation as expected:
![tutorial](https://github.com/Ionsto/cl-mpm/assets/117826225/95163c5a-35a8-4312-a71c-c9a4941ee388)

To change consitutitve model to a von-mises perfectly plastic model, we just change our setup code from 
```lisp
...
'cl-mpm/particle::particle-elastic
:E 1d5 ;; Young's modulus 1GPa
:nu 0.35d0 ;; Poission's ratio
...
```
To a von-mises particle with an extra rho field.
```lisp
...
'cl-mpm/particle::particle-vm
:E 1d5 ;; Young's modulus 1GPa
:nu 0.35d0 ;; Poission's ratio
:rho 1.5d3 ;; Yield stress
...
```
This results in a much greater displacement under plastic deformation
![tutorial](https://github.com/Ionsto/cl-mpm/assets/117826225/e3b55ef9-f5aa-419f-8abf-84be209ab2c0)

# MPI
To use mpi you need to compile an image worker with a guide here: https://github.com/samjvsutcliffe/cl-mpm-worker
You can then use an MPI-ready ```mpm-sim``` class like ```mpm-sim-mpi-nodes```.
[MPI tutorial file](./tutorial/tutorial-mpi/tutorial-mpi.lisp)

When we use MPI domain decomposition, we end up with multiple part-files that may be visualised in paraview.

The material points get distributed across the mpi domains spacially - here with the mpi index of the domains being shown:

![image](https://github.com/Ionsto/cl-mpm/assets/117826225/35ef92d6-4eae-4e29-9ce8-22105ad86136)

In a larger example many domains are used, however the occupancy of each sub domain is not well distributed.

![image](https://github.com/Ionsto/cl-mpm/assets/117826225/a9488de4-d5bc-411e-90aa-9dda1524d773)

