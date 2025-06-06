(defpackage :cl-mpm/examples/aggregate-condition
  (:use :cl
   :cl-mpm/example))
(in-package :cl-mpm/examples/aggregate-condition)
(defun setup (&key (refine 1) (mps 2))
  (let* ((density 1d0)
         (h (/ 1d0 refine))
         (offset (* h 1.5))
         (element-count (list 16 1))
         (block-size (list (* 4 h) h)))
    (defparameter *sim* (cl-mpm/setup::make-simple-sim h
                                               element-count
                                               :sim-type
                                               'cl-mpm/aggregate::mpm-sim-agg-usf))
    (let* ()
      (format t "Mesh size ~F~%" h)
      (cl-mpm:add-mps
       *sim*
       (cl-mpm/setup:make-block-mps
        (list offset 0 0)
        block-size
        (mapcar (lambda (e) (* (/ e h) mps)) block-size)
        density
        'cl-mpm/particle::particle-elastic
        :E 1d9
        :nu 0.325d0
        )))
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 0d-15)
    (cl-mpm/setup::setup-bcs
     *sim*
     :top '(nil nil nil)
     :bottom '(nil 0 nil))
    *sim*))

(defun condition-number (ma)
  (multiple-value-bind (u s v) (magicl:svd ma)
    (let ((sv (magicl::diag s)))
      (/ (reduce #'max sv)
         (reduce #'min sv)))))

(defun test-condition ()
  (setup)
  (let* ((h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
         (resolution 10)
         (x (coerce (loop for x from 0d0 below (* 1 h) by (/ h 100) collect x) 'vector))
         (condition-ma (make-array (length x)))
         (condition-mii (make-array (length x))))
    (loop for p across x
          for i from 0
          while (cl-mpm::sim-run-sim *sim*)
          do
             (setup)
             (cl-mpm::iterate-over-mps
              (cl-mpm:sim-mps *sim*)
              (lambda (mp)
                (incf (cl-mpm/utils::varef (cl-mpm/particle::mp-position mp) 0) p)))
             (cl-mpm:update-sim *sim*)
             (plot-domain)
             ;; (plot-condition)
             (multiple-value-bind (ma mii) (cl-mpm/aggregate::assemble-global-full-mass *sim*)
               (setf (aref condition-ma i) (condition-number ma)
                     (aref condition-mii i) (condition-number mii))))
    (defparameter *cond-x* x)
    (defparameter *cond-ma* condition-ma)
    (defparameter *cond-mii* condition-mii)

    (vgplot:figure)
    (plot-condition)
    )
  )
(defun plot-condition ()
    (vgplot:semilogy *cond-x* *cond-ma* "m^a"
     *cond-x* *cond-mii* "m lumped"))(list h)
