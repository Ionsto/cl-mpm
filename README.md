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
