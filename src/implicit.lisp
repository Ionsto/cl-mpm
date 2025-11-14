(defpackage :cl-mpm/implicit
  (:use :cl
   :cl-mpm
   :cl-mpm/particle
   :cl-mpm/mesh
   :cl-mpm/utils
   :cl-mpm/fastmaths
   )
  (:export
   #:mpm-sim-agg-usf)
  )
(in-package :cl-mpm/implicit)

(defmacro project-global-vec (sim vector accessor)
  `(progn
     (let* ((active-nodes (cl-mpm/aggregate::sim-agg-nodes-fd ,sim))
            (proj-val ,vector)
            (nd 3))
       (cl-mpm::iterate-over-nodes-array
        active-nodes
        (lambda (n)
          (loop for d from 0 below nd
                do (setf (cl-mpm/utils:varef (funcall ,accessor n) d)
                         (cl-mpm/utils:varef proj-val (+ (* nd (cl-mpm/mesh::node-agg-fd n)) d))))))
       )))

(defun assemble-global-vec (sim accessor)
  (declare (function accessor))
  (let* ((active-nodes (cl-mpm/aggregate::sim-agg-nodes-fd sim))
         (nd 3)
         (ndof (* nd (length active-nodes)))
         (v (cl-mpm/utils::arb-matrix ndof 1)))
    (cl-mpm::iterate-over-nodes-array
     active-nodes
     (lambda (n)
       (loop for d from 0 below nd
             do (setf (cl-mpm/utils:varef v (+ (* nd (cl-mpm/mesh::node-agg-fd n)) d))
                      (cl-mpm/utils:varef (funcall accessor n) d)))))
    (values v)))

(defun assemble-global-bcs (sim)
  (let* ((active-nodes (cl-mpm/aggregate::sim-agg-nodes-fd sim))
         (nd 3)
         (ndof (* nd (length active-nodes)))
         (v (cl-mpm/utils::arb-matrix ndof 1)))
    (cl-mpm::iterate-over-nodes-array
     active-nodes
     (lambda (n)
       (loop for d from 0 below nd
             do
                (let ((index (+ d (* nd (cl-mpm/mesh::node-agg-fd n)))))
                  (let ((bcs (cl-mpm/mesh::node-bcs n)))
                    (if bcs
                        (setf (cl-mpm/utils:varef v index) (cl-mpm/utils:varef bcs d))
                        (setf (cl-mpm/utils:varef v index) 1d0))))
                )))
    (values v)))

(defun setup-implicit (sim)
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (mass-filter cl-mpm::mass-filter)
               (initial-setup initial-setup))
      sim
    (cl-mpm::reset-grid mesh)
    (cl-mpm::reset-node-displacement sim)
    (cl-mpm::p2g mesh mps)
    (when (> mass-filter 0d0)
      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
    (cl-mpm::filter-cells sim)
    (setf
     (cl-mpm/aggregate::sim-agg-nodes-fd sim)
     (cl-mpm/aggregate::filter-nodes sim #'cl-mpm/mesh::node-active))
    (let ((fd 0))
      (loop for n across (cl-mpm/aggregate::sim-agg-nodes-fd sim)
            do (progn
                 (setf (cl-mpm/mesh::node-agg-fd n) fd)
                 (incf fd))))
    )
  )
(defun reduce-with-bcs (v bcs)
  (let ((bc-map (make-array (magicl:nrows bcs) :fill-pointer 0)))
    ;; (cl-mpm/fastmaths:fast-zero target-vi)
    (loop for i from 0
          for bc across (magicl::storage bcs)
          do (when (> bc 0d0)
               (vector-push-extend i bc-map)))
    (let ((reduced-size (length bc-map)))
      (if (> reduced-size 0)
        (let* ((v-r (cl-mpm/utils::arb-matrix reduced-size 1)))
          (lparallel:pdotimes (i reduced-size)
            (setf (cl-mpm/utils:varef v-r i) (cl-mpm/utils:varef v (aref bc-map i))))
          ;; (let ((vs (magicl::linear-solve A-r v-r)))
          ;;   (lparallel:pdotimes (i reduced-size)
          ;;     (setf (varef target-vi (aref bc-map i)) (varef vs i))))
          v-r
          )
        (error "no degrees of freedom")))))

(defun expand-with-bcs (v bcs)
  (let* ((expanded-size (magicl:nrows bcs))
         (target (cl-mpm/utils::arb-matrix expanded-size 1))
         (bc-map (make-array (magicl:nrows bcs) :fill-pointer 0)))
    (loop for i from 0
          for bc across (magicl::storage bcs)
          do (when (> bc 0d0)
               (vector-push-extend i bc-map)))
    (lparallel:pdotimes (i (magicl:nrows v))
      (setf (cl-mpm/utils:varef target (aref bc-map i)) (cl-mpm/utils:varef v i)))
    target))

(defun linear-solve-with-bcs (ma v bcs &optional (target-vi nil))
  (let ((target-vi (if target-vi
                       target-vi
                       (cl-mpm/utils::arb-matrix (magicl:nrows v) 1)
                       ))
        (bc-map (make-array (magicl:nrows bcs) :fill-pointer 0)))
    (cl-mpm/fastmaths:fast-zero target-vi)
    (loop for i from 0
          for bc across (magicl::storage bcs)
          do (when (> bc 0d0)
               (vector-push-extend i bc-map)))
    (let ((reduced-size (length bc-map)))
      (when (> reduced-size 0)
        (let* ((A-r (cl-mpm/utils::arb-matrix reduced-size reduced-size))
               (v-r (cl-mpm/utils::arb-matrix reduced-size 1)))
          (lparallel:pdotimes (i reduced-size)
            (dotimes (j reduced-size)
              (setf (mtref a-r i j) (mtref ma (aref bc-map i) (aref bc-map j))))
            (setf (varef v-r i) (varef v (aref bc-map i))))
          (let ((vs (magicl::linear-solve A-r v-r)))
            (lparallel:pdotimes (i reduced-size)
              (setf (varef target-vi (aref bc-map i)) (varef vs i))))
          )))
    target-vi
    ))

(defun test-imp ()
  (setup :refine 0.5)
  ;; (let ((d  (assemble-global-vec *sim* #'cl-mpm/mesh::node-displacment)))
  ;;   (setf (cl-mpm/utils:varef d 0) 1d0)
  ;;   (project-global-vec *sim* d cl-mpm/mesh::node-displacment)
  ;;   )
  (let ((lstps 50))
    (loop for l from 1 to lstps
          do
             (let ((sim *sim*)
                   (iter 0))

               (setf (cl-mpm::sim-ghost-factor *sim*) (* 1d6 1d-3))
               (setup-implicit *sim*)
               (setf (cl-mpm/aggregate::sim-enable-aggregate sim) nil)
               (with-slots ((mesh cl-mpm::mesh)
                            (mps cl-mpm::mps)
                            (bcs cl-mpm::bcs)
                            (dt-loadstep cl-mpm/dynamic-relaxation::dt-loadstep)
                            (dt cl-mpm::dt)
                            (ghost-factor cl-mpm::ghost-factor)
                            (fbar cl-mpm::enable-fbar)
                            (initial-setup initial-setup))
                   *sim*
                 (cl-mpm::p2g-force-fs *sim*)
                 (cl-mpm::apply-bcs mesh bcs dt)
                 (let* ((bcs-vec (assemble-global-bcs sim))
                        (f-ext-full (assemble-global-vec sim #'cl-mpm/mesh::node-external-force))
                        (f-ext
                          (reduce-with-bcs
                           f-ext-full
                           bcs-vec)))
                   (cl-mpm/fastmaths:fast-scale! f-ext (* -1d0 (/ (float l) lstps)))
                   (cl-mpm/linear-solver::solve-richardson
                    (lambda (d)
                      (project-global-vec sim
                                          (expand-with-bcs
                                           d
                                           bcs-vec)
                                          #'cl-mpm/mesh::node-displacment)
                      (cl-mpm::apply-bcs mesh bcs dt)
                      (cl-mpm::reset-nodes-force sim)
                      (cl-mpm::update-stress mesh mps dt-loadstep fbar)
                      (cl-mpm::p2g-force-fs sim)
                      (when ghost-factor
                        (cl-mpm/ghost::apply-ghost sim ghost-factor)
                        (cl-mpm::apply-bcs mesh bcs dt))
                      (cl-mpm::apply-bcs mesh bcs dt)
                      (incf iter)
                      (cl-mpm:iterate-over-nodes
                       mesh
                       (lambda (n)
                         (when (cl-mpm/mesh::node-active n)
                           (cl-mpm/fastmaths:fast-.+ (cl-mpm/mesh::node-internal-force n)
                                                     (cl-mpm/mesh::node-ghost-force n)
                                                     (cl-mpm/mesh::node-force n)))))
                      (reduce-with-bcs
                       (assemble-global-vec sim #'cl-mpm/mesh::node-force)
                       bcs-vec))
                    f-ext
                    :tol 1d-3
                    )))
               (cl-mpm::finalise-loadstep sim)
               (cl-mpm/output:save-vtk (merge-pathnames "./output/" (format nil "sim_~5,'0d.vtk" l)) sim)
               (cl-mpm/output:save-vtk-nodes (merge-pathnames "./output/" (format nil "sim_n_~5,'0d.vtk" l)) sim)
               (plot *sim*)
               (pprint iter)))))

(defun test-nr-imp ()
  (setup :refine 0.5)
  (setup-implicit *sim*)
  ;; (let ((d  (assemble-global-vec *sim* #'cl-mpm/mesh::node-displacment)))
  ;;   (setf (cl-mpm/utils:varef d 0) 1d0)
  ;;   (project-global-vec *sim* d cl-mpm/mesh::node-displacment)
  ;;   )
  (let ((sim *sim*)
        (iter 0))
    (setf (cl-mpm/aggregate::sim-enable-aggregate sim) nil)

    (let* ((bcs-vec (assemble-global-bcs sim))
           (d-0
             (reduce-with-bcs
              (assemble-global-vec sim #'cl-mpm/mesh::node-displacment)
              bcs-vec))
           (nr-crit 1d-9)
           (nr-error nr-crit))
      ;;NR iter
      (loop for niter from 0 to 10
            while (>= nr-error nr-crit)
            do
               (progn
                 (with-slots ((mesh cl-mpm::mesh)
                              (mps cl-mpm::mps)
                              (bcs cl-mpm::bcs)
                              (dt-loadstep cl-mpm/dynamic-relaxation::dt-loadstep)
                              (dt cl-mpm::dt)
                              (fbar cl-mpm::enable-fbar)
                              (initial-setup initial-setup))
                     *sim*
                   (cl-mpm::p2g-force-fs *sim*)
                   (cl-mpm::apply-bcs mesh bcs dt)
                   (let* ((f-ext
                            (reduce-with-bcs
                             (assemble-global-vec sim #'cl-mpm/mesh::node-external-force)
                             bcs-vec)))
                     (cl-mpm::reset-node-displacement sim)
                     (cl-mpm/fastmaths:fast-scale! f-ext 1d-3)
                     (cl-mpm/linear-solver::solve-conjugant-gradients
                      (lambda (d)
                        (project-global-vec sim
                                            (expand-with-bcs
                                             d
                                             bcs-vec)
                                            #'cl-mpm/mesh::node-displacment)
                        ;; (pprint d)
                        ;; (cl-mpm::reset-nodes-force sim)
                        (cl-mpm::apply-bcs mesh bcs dt)
                        ;; (cl-mpm::update-nodes sim)
                        ;; (cl-mpm::update-cells sim)
                        (cl-mpm::reset-nodes-force sim)
                        (cl-mpm::update-stress mesh mps dt-loadstep fbar)
                        (cl-mpm::p2g-force-fs sim)
                        ;; (cl-mpm::update-node-forces sim)
                        (cl-mpm::apply-bcs mesh bcs dt)
                        ;; (plot sim)
                        ;; (cl-mpm/output:save-vtk (merge-pathnames "./output/" (format nil "sim_~5,'0d.vtk" iter)) sim)
                        ;; (cl-mpm/output:save-vtk-nodes (merge-pathnames "./output/" (format nil "sim_n_~5,'0d.vtk" iter)) sim)
                        (incf iter)
                        (reduce-with-bcs
                         (assemble-global-vec sim #'cl-mpm/mesh::node-internal-force)
                         bcs-vec))
                      f-ext)))
                 (let ((r
                         (reduce-with-bcs
                          (cl-mpm/fastmaths::fast-.+
                           (assemble-global-vec sim #'cl-mpm/mesh::node-external-force)
                           (assemble-global-vec sim #'cl-mpm/mesh::node-internal-force))
                          bcs-vec)))
                   (setf nr-crit (cl-mpm/fastmaths::mag-squared r)))
                 (format t "NR iter ~D - ~E ~%" niter nr-crit)
                 )))))
