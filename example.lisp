(ql:quickload "cl-mpm")


(defparameter *nD* 2)
(defparameter *mesh* (cl-mpm::make-mesh *nD* '(1) 1 (cl-mpm::make-shape-function-linear *nD*)))
(defparameter *mps* (list (cl-mpm::make-particle *nD* :pos '(0.5))))
(cl-mpm::p2g *mesh* *mps*)
(describe *mesh*)
(describe (slot-value *mesh* 'shape-func))
(with-slots ((shape-func shape-func))
  *mesh*
  (with-slots ((svp svp)
               (dsvp dsvp)
               )
    shape-func

    (funcall svp 1)
    (funcall dsvp 1)
    ))
(describe (get-node *mesh* '(0)))
(describe (first *mps*))
(defparameter *sim* (cl-mpm:make-mpm-sim 2 '(1 1) 1 (cl-mpm::make-shape-function-linear 2)))
(cl-mpm::update-sim *sim*)
