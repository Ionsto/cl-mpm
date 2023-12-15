(restrict-compiler-policy 'speed  0 0)
(restrict-compiler-policy 'debug  3 3)
(restrict-compiler-policy 'safety 3 3)
(ql:quickload :cl-mpm-worker)
(in-package :cl-mpm-worker)

(ql:quickload :cl-mpm)
(ql:quickload :cl-mpm/setup)
(ql:quickload :cl-mpm/particle)
(ql:quickload :cl-mpm/mpi)

(asdf:compile-system :cl-mpm/mpi :force t)

(defun primary-main ()
  ;;By default the compiled lisp interperater will load a file called "test_file.lisp"
  (load "tutorial-mpi.lisp"))
(sb-ext:save-lisp-and-die
 "mpi-worker"
 :executable t
 :toplevel #'main
 :save-runtime-options t)
(uiop:quit)
