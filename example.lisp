(ql:quickload "cl-mpm")


(defparameter *nD* 2)
(defparameter *mesh* (cl-mpm::make-mesh *nD* '(1 1) 1 (cl-mpm::make-shape-function-linear *nD*)))
(defparameter *mps* (list (cl-mpm::make-particle *nD* :pos '(0.5 0.5))))
(describe (get-node *mesh* '(0)))
(describe (first *mps*))
(defparameter *sim* (cl-mpm:make-mpm-sim *nD* '(1 1) 1 (cl-mpm::make-shape-function-linear *nD*)))
(cl-mpm::update-sim *sim*)
