(defpackage :cl-mpm/examples/aggregate-condition
  (:use :cl
   :cl-mpm/example))
(in-package :cl-mpm/examples/aggregate-condition)

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-elastic) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  ;; (cl-mpm::update-domain-polar-2d mesh mp dt)
  ;; (cl-mpm::scale-domain-size mesh mp)
  )

(defun setup (&key (refine 1) (mps 2))
  (let* ((density 1d0)
         (h (/ 1d0 refine))
         (offset (* h 2.5))
         (element-count (list 16 1))
         ;; (block-size (list (* 4 h) h))
         (block-size (list (* 4 h) (* 1 h)))
         )
    (defparameter *sim* (cl-mpm/setup::make-simple-sim h
                                               element-count
                                               :sim-type
                                               'cl-mpm/aggregate::mpm-sim-agg-usf))
    (let* ()
      (format t "Mesh size ~F~%" h)
      (cl-mpm:add-mps
       *sim*
       (cl-mpm/setup:make-block-mps
        (list offset 0d0 0)
        block-size
        (mapcar (lambda (e) (* (/ e h) mps)) block-size)
        density
        'cl-mpm/particle::particle-elastic
        :E 1d9
        :nu 0.325d0
        )))
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 0d-15)
    (setf (cl-mpm:sim-mass-filter *sim*) 0d0)
    (cl-mpm/setup::setup-bcs
     *sim*
     :top '(nil nil nil)
     :bottom '(nil 0 nil))
    *sim*))

(defun condition-number (ma)
  (multiple-value-bind (u s v) (magicl:svd ma)
    (let ((sv (magicl::diag s)))
      (let ((smin (reduce #'min sv))
            (smax (reduce #'max sv))
            )
        (if (> smin 0d0)
            (/ smax smin)
            sb-ext::most-positive-double-float)))))

(defun test-condition (&key (resolution 10)
                         (save-vtk nil)
                         )
  (ensure-directories-exist (merge-pathnames "./output/"))
  (setup)
  (let* ((h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
         (dx (/ h resolution))
         (x (coerce (loop for x from 0d0 below (* 1 h) by dx collect x) 'vector))
         (condition-ma (make-array (length x)))
         (condition-mii (make-array (length x))))
    (defparameter *cond-x* x)
    (defparameter *cond-ma* condition-ma)
    (defparameter *cond-mii* condition-mii)
    (setup)
    (loop for p across x
          for i from 0
          while (cl-mpm::sim-run-sim *sim*)
          do
             (cl-mpm::iterate-over-mps
              (cl-mpm:sim-mps *sim*)
              (lambda (mp)
                (incf (cl-mpm/utils::varef (cl-mpm/particle::mp-position mp) 0) dx)))
             (cl-mpm:update-sim *sim*)
             ;; (plot-domain)
             (plot-condition)
             (let* ((mii (cl-mpm/aggregate::assemble-global-mass-matrix *sim*))
                    (E (cl-mpm/aggregate::assemble-e *sim*))
                    (ma (magicl:@  (magicl:transpose E) mii E)))
               (setf (aref condition-ma i) (condition-number ma)
                     (aref condition-mii i) (condition-number mii)))
             (when save-vtk
               (save-vtks "./output/" i))
             ;; (multiple-value-bind (ma mii) (cl-mpm/aggregate::assemble-global-full-mass *sim*)
             ;;   (setf (aref condition-ma i) (condition-number ma)
             ;;         (aref condition-mii i) (condition-number mii)))
          )
    (vgplot:figure)
    (plot-condition)
    ;; (save-data)
    ))

(defun test-cfl (&key (resolution 10)
                         (save-vtk nil)
                         )
  (ensure-directories-exist (merge-pathnames "./output/"))
  (setup)
  (let* ((h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
         (dx (/ h resolution))
         (x (coerce (loop for x from 0d0 below (* 1 h) by dx collect x) 'vector))
         (condition-ma (make-array (length x)))
         (condition-mii (make-array (length x))))
    (defparameter *cond-x* x)
    (defparameter *cond-ma* condition-ma)
    (defparameter *cond-mii* condition-mii)
    (setup)
    (loop for p across x
          for i from 0
          while (cl-mpm::sim-run-sim *sim*)
          do
             (cl-mpm::iterate-over-mps
              (cl-mpm:sim-mps *sim*)
              (lambda (mp)
                (incf (cl-mpm/utils::varef (cl-mpm/particle::mp-position mp) 0) dx)))
             (cl-mpm:update-sim *sim*)
             ;(plot-condition)
             (plot-domain)
             (setf (aref condition-ma i) (cl-mpm::calculate-min-dt *sim*))

             (cl-mpm::iterate-over-nodes
              (cl-mpm:sim-mesh *sim*)
              (lambda (n)
                (setf (cl-mpm/mesh::node-agg n) nil)))

             (setf (aref condition-mii i) (cl-mpm::calculate-min-dt *sim*))
             ;; (let* ((mii (cl-mpm/aggregate::assemble-global-mass *sim*))
             ;;        (E (cl-mpm/aggregate::assemble-global-e *sim*))
             ;;        (ma (magicl:@  (magicl:transpose E) mii E)))
             ;;   (setf (aref condition-ma i) (condition-number ma)
             ;;         (aref condition-mii i) (condition-number mii)))
             (when save-vtk
               (save-vtks "./output/" i)))
    (vgplot:figure)
    (plot-condition)
    ))




(defun vec-to-single-float (vec)
  (map '(vector single-float) (lambda (x) (coerce x 'single-float)) vec))

(defun save-data ()
  (let* ((output-file (merge-pathnames "./analysis_scripts/aggregate/1d_rigid/data.csv"))
         (df (lisp-stat:make-df '(:x :ma :m-lumped)
                                (list
                                 (vec-to-single-float *cond-x*)
                                 (vec-to-single-float *cond-ma*)
                                 (vec-to-single-float *cond-mii*)))))
    (lisp-stat:write-csv df output-file :add-first-row t))
  )
(defun save-vtks (output-dir step)
  (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" step)) *sim*)
  (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) *sim*)
  (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d.vtk" step)) *sim*))

(defun plot-condition ()
  (vgplot:semilogy *cond-x* *cond-ma* "m^a"
                   *cond-x* *cond-mii* "m lumped"))
