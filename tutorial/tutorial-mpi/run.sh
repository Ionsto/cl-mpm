sbcl --dynamic-space-size 4000 --load "build_step.lisp" --quit
mpirun -N 2 ./mpi-worker --dynamic-space-size 8000
