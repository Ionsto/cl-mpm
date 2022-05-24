(ql:quickload "cl-mpm")
(cl-mpm)
(cl-mpm:make-shape-function  2)
(cl-mpm:make-mpm-sim 2 '(1 1) 1 (cl-mpm:sha) )
