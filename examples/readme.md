# Example files for cl-mpm
A base package is stored in core.lisp, defining a 2d plotting function ```(plot)``` which plots the simulation stored in the global variable ```*sim*```.
In general, each example will have 2 main functions, ```(setup)``` and ```(run)```.
Setup takes two keyword arguments ```:refine N``` which refines the problem by a factor of N, and ```:mps N``` which sets the MPs per cell to NxN.
Run may have additional keyword arguments.

## AMPLE examples
The three quasi-static examples from the AMPLE paper [AMPLE: A Material Point Learning Environment][https://www.sciencedirect.com/science/article/pii/S0965997819304648] are implemented in the examples/ample folder.
This includes:
- 1D self weight compression
- quasi-static elastic beam bending
- 2d von-Mises elasto-plastic collapse

## Chalk
The problem setup for the chalk cliff domain from [Modelling cliff collapse and run-out with the material point method][https://www.sciencedirect.com/science/article/pii/S0266352X25005579] is in examples/chalk/joss.

## Penalty contact
Test cases for the penalty contact formulation are in the examples/penalty folder, including a frictional contact test and the difficult billet upsetting problem.
