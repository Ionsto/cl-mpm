(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)

(defun combi-stats (sim)
  (destructuring-bind (mass
                       energy
                       oobf-num
                       oobf-denom
                       power)
      (cl-mpm::reduce-over-nodes
       (cl-mpm:sim-mesh sim)
       (lambda (node)
         (if (and (cl-mpm/mesh:node-active node)
                  ;; (or
                  ;;  (not (cl-mpm/mesh::node-agg node))
                  ;;  ;; (cl-mpm/mesh::node-interior node)
                  ;;  )
                  )
             (with-accessors ((active cl-mpm/mesh::node-active)
                              (f-ext cl-mpm/mesh::node-external-force)
                              ;; (f-int cl-mpm/mesh::node-residual)
                              (res cl-mpm/mesh::node-residual)
                              (node-oobf cl-mpm/mesh::node-oobf)
                              (mass cl-mpm/mesh::node-mass)
                              (volume cl-mpm/mesh::node-volume)
                              (volume-t cl-mpm/mesh::node-volume-true)
                              (vel cl-mpm/mesh::node-velocity)
                              (disp cl-mpm/mesh::node-displacment)
                              )
                 node
               (declare (double-float mass))
               (let (;(mass 1d0)
                     ;; (scale-factor (expt mass 1))
                     (scale-factor 1d0)
                     ;; (scale-factor (/ volume volume-t))
                     )
                 (list
                  scale-factor
                  (* scale-factor (* 0.5d0 mass (cl-mpm/fastmaths::mag-squared vel)))
                  ;(* scale-factor (cl-mpm/fastmaths::mag (cl-mpm/fastmaths::fast-.+-vector f-ext f-int)))
                  (* scale-factor (cl-mpm/fastmaths::mag-squared res))
                  (* scale-factor (cl-mpm/fastmaths::mag-squared f-ext))
                  (* scale-factor
                     (cl-mpm/fastmaths:dot
                      disp f-ext)))))
             (list 0d0 0d0 0d0 0d0 0d0)))
       (lambda (a b) (mapcar (lambda (x y) (declare (double-float x y)) (+ x y)) a b)))
    (declare (double-float mass energy oobf-num oobf-denom power))
    (let ((oobf 0d0))
      ;; (format t "~E ~E~%" oobf-num oobf-denom)
      (if (> oobf-denom 0d0)
          (setf oobf (sqrt (/ oobf-num oobf-denom)))
          (setf oobf (if (> oobf-num 0d0) sb-ext:double-float-positive-infinity 0d0)))
      (when (> oobf-denom 0d0)
        (cl-mpm::iterate-over-nodes
         (cl-mpm:sim-mesh sim)
         (lambda (node)
           (with-accessors ((active cl-mpm/mesh::node-active)
                            (agg cl-mpm/mesh::node-agg)
                            ;; (f-ext cl-mpm/mesh::node-external-force)
                            (res cl-mpm/mesh::node-residual)
                            (n-mass cl-mpm/mesh::node-mass)
                            (node-oobf cl-mpm/mesh::node-oobf))
               node
             (if (and active
                        (or
                         (not agg) (cl-mpm/mesh::node-interior node)))
               (when t
                 (setf node-oobf
                       (if (> oobf 0d0)
                           (/ (* n-mass (cl-mpm/fastmaths::mag res))
                               oobf-denom)
                           0d0
                           )))
               (setf node-oobf 0d0))))))
      (if (> mass 0d0)
          ;(values (/ energy mass) oobf (/ power mass))
          (values energy oobf power)
          (values 0d0 0d0 0d0)))))

(defun combi-stats-aggregated (sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (sim-agg cl-mpm/aggregate::sim-enable-aggregate))
      sim
    (destructuring-bind (energy
                         oobf-num
                         oobf-denom
                         power)
        (cl-mpm::reduce-over-nodes
         (cl-mpm:sim-mesh sim)
         (lambda (node)
           (if (and (cl-mpm/mesh:node-active node)
                    (not (cl-mpm/mesh::node-agg node)))
               (with-accessors ((active cl-mpm/mesh::node-active)
                                (f-ext cl-mpm/mesh::node-external-force)
                                (res cl-mpm/mesh::node-residual)
                                (f-int cl-mpm/mesh::node-internal-force)
                                (node-oobf cl-mpm/mesh::node-oobf)
                                (mass cl-mpm/mesh::node-mass)
                                (volume cl-mpm/mesh::node-volume)
                                (volume-t cl-mpm/mesh::node-volume-true)
                                (vel cl-mpm/mesh::node-velocity)
                                (disp cl-mpm/mesh::node-displacment)
                                )
                   node
                 (declare (double-float mass))
                 (let ()
                   (list
                    (* 0.5d0 mass (cl-mpm/fastmaths::mag-squared vel))
                    (cl-mpm/fastmaths::mag-squared res)
                    (cl-mpm/fastmaths::mag-squared f-ext)
                    (cl-mpm/fastmaths:dot disp f-ext))))
               (list 0d0 0d0 0d0 0d0)))
         (lambda (a b) (mapcar (lambda (x y) (declare (double-float x y)) (+ x y)) a b)))
      (declare (double-float energy oobf-num oobf-denom power))
      (when sim-agg
        (loop for d from 0 below (cl-mpm/mesh::mesh-nd mesh)
              do (let* ((res (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-residual d))
                        (f-ext (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-external-force d))
                        (E (cl-mpm/aggregate::sim-global-e sim))
                        (ma (cl-mpm/aggregate::sim-global-ma sim))
                        (mv (cl-mpm/aggregate::assemble-global-scalar sim #'cl-mpm/mesh::node-mass))
                        (vi (magicl:@
                             (magicl:transpose E)
                             (cl-mpm/fastmaths:fast-.*
                              mv
                              (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-velocity d))))
                        (disp (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-displacment d)))

                   ;; (setf vi (cl-mpm/aggregate::apply-internal-bcs sim vi d))
                   (setf
                    vi
                    (cl-mpm/aggregate::linear-solve-with-bcs ma vi (cl-mpm/aggregate::assemble-internal-bcs sim d)))
                   ;; (setf res (cl-mpm/aggregate::apply-internal-bcs sim res d))
                   ;; (setf f-ext (cl-mpm/aggregate::apply-internal-bcs sim f-ext d))
                   ;; (setf disp (cl-mpm/aggregate::apply-internal-bcs sim disp d))
                   ;; (setf vi (magicl:@ (magicl:transpose E) vi))

                   (incf oobf-num (cl-mpm/fastmaths::mag-squared
                                   (cl-mpm/aggregate::apply-internal-bcs
                                    sim
                                    (magicl:@ (magicl:transpose E) res)
                                    d
                                    )))
                   (incf oobf-denom (cl-mpm/fastmaths::mag-squared
                                     (cl-mpm/aggregate::apply-internal-bcs
                                      sim
                                      (magicl:@ (magicl:transpose E) f-ext)
                                      d)))
                   (incf power (cl-mpm/fastmaths:dot
                                disp f-ext))
                   (incf energy (* 0.5d0 (cl-mpm/utils::mtref (magicl:@ (magicl:transpose vi) ma vi) 0 0)))
                   )))
      (let ((oobf 0d0)
            (oobf-num (sqrt oobf-num))
            (oobf-denom (sqrt oobf-denom)))
        (if (> oobf-denom 0d0)
            (setf oobf (/ oobf-num oobf-denom))
            (setf oobf (if (> oobf-num 0d0) sb-ext:double-float-positive-infinity 0d0)))
        (when (> oobf-denom 0d0)
          (cl-mpm::iterate-over-nodes
           (cl-mpm:sim-mesh sim)
           (lambda (node)
             (with-accessors ((active cl-mpm/mesh::node-active)
                              (agg cl-mpm/mesh::node-agg)
                              (f-ext cl-mpm/mesh::node-external-force)
                              ;; (f-int cl-mpm/mesh::node-internal-force)
                              (f-int cl-mpm/mesh::node-residual)
                              (n-mass cl-mpm/mesh::node-mass)
                              (node-oobf cl-mpm/mesh::node-oobf))
                 node
               (if (and active
                          (not agg))
                 (when t
                   (setf node-oobf
                         (if (> oobf 0d0)
                             (/ (* n-mass (cl-mpm/fastmaths::mag (cl-mpm/fastmaths::fast-.+-vector f-ext f-int)))
                                oobf-denom)
                             0d0
                             )))
                 (setf node-oobf 0d0)))))

          (when sim-agg
            ;; (print "Hello")
            (loop for d from 0 below (cl-mpm/mesh::mesh-nd mesh)
                  do (cl-mpm/aggregate::zero-global-scalar sim cl-mpm/mesh::node-oobf))
            (loop for d from 0 below (cl-mpm/mesh::mesh-nd mesh)
                  do (let* ((E (cl-mpm/aggregate::sim-global-e sim))
                            ;; (f-int (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-residual d))
                            (f-int (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-internal-force d))
                            (f-ext (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-external-force d)))

                       (let* ((fagg (cl-mpm/aggregate::apply-internal-bcs
                                     sim
                                     (magicl:@ (magicl:transpose E)
                                               (cl-mpm/fastmaths::fast-.+
                                                f-int
                                                f-ext))
                                     d))
                              (fnorm (cl-mpm/fastmaths:fast-.* fagg fagg)))
                         ;; (pprint fnorm)
                         (cl-mpm/aggregate::increment-global-scalar sim (magicl:@ E fnorm) cl-mpm/mesh::node-oobf)))))
          (cl-mpm::iterate-over-nodes
           (cl-mpm:sim-mesh sim)
           (lambda (node)
             (with-accessors ((active cl-mpm/mesh::node-active)
                              (agg cl-mpm/mesh::node-agg)
                              (node-oobf cl-mpm/mesh::node-oobf))
                 node
               (when (and active
                        agg)
                   ;; (setf node-oobf (/ (sqrt (max node-oobf 0d0)) oobf-denom))
                   ))))
          )
        (values energy oobf power)
        ))))

(defgeneric estimate-energy-norm (sim))
(defmethod estimate-energy-norm ((sim cl-mpm::mpm-sim))
  (let ((energy 0d0)
        (mass 0d0)
        (lock (sb-thread:make-mutex))
        )
    (cl-mpm:iterate-over-nodes
     (cl-mpm:sim-mesh sim)
     (lambda (n)
       (when (cl-mpm/mesh:node-active n)
         (sb-thread:with-mutex (lock)
           (incf mass (cl-mpm/mesh:node-mass n))
           (incf energy
                 (*
                  ;; (/ (cl-mpm/mesh::node-volume n) (cl-mpm/mesh::node-volume-true n))
                  (cl-mpm/mesh:node-mass n)
                  (cl-mpm/mesh::node-mass n)
                  (cl-mpm/fastmaths::mag-squared (cl-mpm/mesh::node-velocity n))
                  ))))))
    ;; (setf
    ;;  energy
    ;;  (cl-mpm::reduce-over-nodes
    ;;   (cl-mpm:sim-mesh sim)
    ;;   (lambda (n)
    ;;     (if (cl-mpm/mesh:node-active n)
    ;;         (*
    ;;          (cl-mpm/mesh::node-mass n)
    ;;          (cl-mpm/fastmaths::mag-squared (cl-mpm/mesh::node-velocity n)))
    ;;         0d0
    ;;         ))
    ;;   #'+))
    ;; (setf
    ;;  mass
    ;;  (cl-mpm::reduce-over-nodes
    ;;   (cl-mpm:sim-mesh sim)
    ;;   (lambda (n)
    ;;     (if (cl-mpm/mesh:node-active n)
    ;;         (cl-mpm/mesh::node-mass n)
    ;;         0d0))
    ;;   #'+))
    (if (= mass 0d0)
        0d0
        (/ energy mass))))
(defmethod estimate-energy-norm ((sim cl-mpm/mpi::mpm-sim-mpi))
  (cl-mpm/mpi::mpi-sum
   (let ((energy 0d0)
         (mass 0d0)
         (lock (sb-thread:make-mutex)))
     (cl-mpm:iterate-over-nodes
      (cl-mpm:sim-mesh sim)
      (lambda (n)
        (when (cl-mpm/mesh:node-active n)
          (when (cl-mpm/mpi::in-computational-domain sim (cl-mpm/mesh::node-position n))
            (sb-thread:with-mutex (lock)
              (incf mass (cl-mpm/mesh:node-mass n))
              (incf energy
                    (*
                     ;; (/ (cl-mpm/mesh::node-volume n) (cl-mpm/mesh::node-volume-true n))
                     (cl-mpm/mesh:node-mass n)
                     (cl-mpm/mesh::node-mass n)
                     (cl-mpm/fastmaths::mag-squared (cl-mpm/mesh::node-velocity n))
                     )))))))
     (if (= mass 0d0)
         0d0
         (/ energy mass)))))

(defgeneric estimate-power-norm (sim))
(defmethod estimate-power-norm ((sim cl-mpm::mpm-sim))
  (let ((dt (cl-mpm:sim-dt sim)))
    (let ((energy 0d0)
          (mass 0d0)
          (lock (sb-thread:make-mutex)))
      (cl-mpm:iterate-over-nodes
       (cl-mpm:sim-mesh sim)
       (lambda (n)
         (when (cl-mpm/mesh:node-active n)
           (sb-thread:with-mutex (lock)
             (incf mass (cl-mpm/mesh:node-mass n))
             (incf energy
                   (*
                    dt
                    ;; (cl-mpm/mesh::node-volume n)
                    ;; (/ (cl-mpm/mesh::node-volume n) (cl-mpm/mesh::node-volume-true n))
                    (cl-mpm/mesh:node-mass n)
                    (cl-mpm/fastmaths:dot
                     ;(cl-mpm/mesh::node-displacment n)
                     (cl-mpm/mesh::node-velocity n)
                     (cl-mpm/mesh::node-external-force n))
                    ))))))
      (if (= mass 0d0)
          0d0
          ;; energy
          (/ energy mass)
          ))))

(defmethod estimate-power-norm ((sim cl-mpm/mpi::mpm-sim-mpi))
  (let ((dt (cl-mpm:sim-dt sim)))
    (cl-mpm/mpi::mpi-sum
     (let ((energy 0d0)
           (mass 0d0)
           (lock (sb-thread:make-mutex)))
       (cl-mpm:iterate-over-nodes
        (cl-mpm:sim-mesh sim)
        (lambda (n)
          (when (cl-mpm/mesh:node-active n)
            (when (cl-mpm/mpi::in-computational-domain sim (cl-mpm/mesh::node-position n))
              (sb-thread:with-mutex (lock)
                (incf mass (cl-mpm/mesh:node-mass n))
                (incf energy
                      (*
                       dt
                       ;; (cl-mpm/mesh::node-volume n)
                       (cl-mpm/mesh:node-mass n)
                       ;; (/ (cl-mpm/mesh::node-volume n) (cl-mpm/mesh::node-volume-true n))
                       (cl-mpm/fastmaths:dot
                        (cl-mpm/mesh::node-velocity n)
                        (cl-mpm/mesh::node-external-force n))
                       )))))))
       (if (= mass 0d0)
           0d0
           (/ energy mass))
       ))))

(defun estimate-oobf-debug (sim)
  (let ((oobf 0d0)
        (nmax 0d0)
        (dmax 0d0)
        (oobf-norm 0d0)
        (lock (sb-thread:make-mutex))
        )
    (cl-mpm::iterate-over-nodes
     (cl-mpm:sim-mesh sim)
     (lambda (node)
       (with-accessors ((active cl-mpm/mesh::node-active)
                        (agg cl-mpm/mesh::node-agg)
                        (f-ext cl-mpm/mesh::node-external-force)
                        (f-int cl-mpm/mesh::node-internal-force)
                        (f-damp cl-mpm/mesh::node-damping-force)
                        (f-ghost cl-mpm/mesh::node-ghost-force)
                        (node-oobf cl-mpm/mesh::node-oobf)
                        )
           node
         (when (and active (not agg))
           (when t;(> (cl-mpm/fastmaths::mag-squared f-ext) 0d0)
             (sb-thread:with-mutex (lock)
               (let ((inc (*
                           (expt (cl-mpm/mesh:node-mass node) 1)
                           ;; (expt (/ (cl-mpm/mesh::node-volume node) (cl-mpm/mesh::node-volume-true node)) 2)
                           (cl-mpm/fastmaths::mag
                            (reduce #'cl-mpm/fastmaths::fast-.+-vector
                                    (list
                                     ;; residual
                                     f-ext
                                     f-int
                                     ;; ;; f-damp
                                     ;; f-ghost
                                     )
                                    )))))
                 (incf node-oobf inc)
                 (setf nmax (+
                             nmax
                             inc)
                       dmax (+
                             dmax
                             (*
                              (expt (cl-mpm/mesh:node-mass node) 1)
                              ;; (expt (/ (cl-mpm/mesh::node-volume node) (cl-mpm/mesh::node-volume-true node)) 2)
                              (cl-mpm/fastmaths::mag
                               f-ext)))))))))))
    (cl-mpm::iterate-over-nodes
     (cl-mpm:sim-mesh sim)
     (lambda (node)
       (with-accessors ((active cl-mpm/mesh::node-active)
                        (agg cl-mpm/mesh::node-agg)
                        (node-oobf cl-mpm/mesh::node-oobf))
           node
         (when (and active (not agg))
           (when t;(> (cl-mpm/fastmaths::mag-squared f-ext) 0d0)
             (sb-thread:with-mutex (lock)
               (setf node-oobf
                     (if (> dmax 0d0)
                         (/ node-oobf dmax)
                         0d0
                         ))
               ))))))
    (if (> dmax 0d0)
      (setf oobf (/ nmax dmax))
      ;;Very odd case, we have external force but no internal forces
      (setf oobf (if (> nmax 0d0) sb-ext:double-float-positive-infinity 0d0)))
    oobf))
(defun estimate-oobf-prod (sim)
(let ((oobf 0d0)
        (nmax 0d0)
        (dmax 0d0)
        (oobf-norm 0d0)
        (lock (sb-thread:make-mutex))
        )
    (cl-mpm::iterate-over-nodes
     (cl-mpm:sim-mesh sim)
     (lambda (node)
       (with-accessors ((active cl-mpm/mesh::node-active)
                        (f-ext cl-mpm/mesh::node-external-force)
                        (f-int cl-mpm/mesh::node-internal-force)
                        (f-damp cl-mpm/mesh::node-damping-force)
                        )
           node
         (when active
           (when t;(> (cl-mpm/fastmaths::mag-squared f-ext) 0d0)
             (sb-thread:with-mutex (lock)
               (setf nmax (+
                           nmax
                           (*
                            (expt (cl-mpm/mesh:node-mass node) 2)
                            ;; (expt (/ (cl-mpm/mesh::node-volume node) (cl-mpm/mesh::node-volume-true node)) 2)
                            (cl-mpm/fastmaths::mag-squared
                             (reduce #'cl-mpm/fastmaths::fast-.+-vector
                                     (list
                                      f-ext
                                      f-int
                                      ;; f-damp
                                      )
                                     ))))
                     dmax (+
                           dmax
                           (*
                            (expt (cl-mpm/mesh:node-mass node) 2)
                            ;; (expt (/ (cl-mpm/mesh::node-volume node) (cl-mpm/mesh::node-volume-true node)) 2)
                            (cl-mpm/fastmaths::mag-squared
                             f-ext))))))))))
    (if (> dmax 0d0)
      (setf oobf (sqrt (/ nmax dmax)))
      ;;Very odd case, we have external force but no internal forces
      (setf oobf (if (> nmax 0d0) sb-ext:double-float-positive-infinity 0d0))
      )
    oobf)
  )


(defun calculate-ke (sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (dt cl-mpm:sim-dt))
      sim
    (let ((ke 0d0))
      (setf
       ke
       (cl-mpm::reduce-over-nodes
        mesh
        (lambda (node)
          (if (and (cl-mpm/mesh:node-active node)
                   (not (cl-mpm/mesh::node-agg node)))
              (cl-mpm/fastmaths::mag-squared
               (cl-mpm/mesh:node-velocity node))
              0d0))
        #'+))

      (when (cl-mpm/aggregate::sim-enable-aggregate sim)
        (loop for d from 0 below (cl-mpm/mesh::mesh-nd mesh)
              do
                 (let* ((ma (cl-mpm/aggregate::sim-global-ma sim))
                        (vi (cl-mpm/aggregate::assemble-internal-vec sim #'cl-mpm/mesh::node-velocity d)))
                   (incf ke (* 0.5d0 (cl-mpm/utils::mtref (magicl:@ (magicl:transpose vi) ma vi) 0 0))))))
      ke)))

(defgeneric estimate-oobf (sim))
(defmethod estimate-oobf (sim)
  "With better reporting"
  )
(defmethod estimate-oobf (sim)
  (estimate-oobf-debug sim))


(defgeneric estimate-static-oobf (sim))

(defmethod estimate-static-oobf ((sim cl-mpm::mpm-sim))
  (estimate-oobf sim))


(defmethod cl-mpm::update-dynamic-stats ((sim cl-mpm:mpm-sim))
  (with-accessors ((stats-energy cl-mpm::sim-stats-energy)
                   (stats-oobf cl-mpm::sim-stats-oobf)
                   (stats-power cl-mpm::sim-stats-power)
                   (stats-work cl-mpm::sim-stats-work))
      sim
    (multiple-value-bind (e o p) (combi-stats sim)
      (setf stats-energy e
            stats-oobf o
            stats-power p)
      (incf stats-work p))))


(defmethod cl-mpm::update-dynamic-stats ((sim cl-mpm/aggregate::mpm-sim-aggregated))
  (with-accessors ((stats-energy cl-mpm::sim-stats-energy)
                   (stats-oobf cl-mpm::sim-stats-oobf)
                   (stats-power cl-mpm::sim-stats-power)
                   (stats-work cl-mpm::sim-stats-work))
      sim
    (if (cl-mpm/aggregate::sim-enable-aggregate sim)
      (multiple-value-bind (e o p) (cl-mpm/dynamic-relaxation::combi-stats-aggregated sim)
        (setf stats-energy e
              stats-oobf o
              stats-power p)
        (incf stats-work p))
      (multiple-value-bind (e o p) (cl-mpm/dynamic-relaxation::combi-stats sim)
        (setf stats-energy e
              stats-oobf o
              stats-power p)
        (incf stats-work p)))
    ;; (setf stats-oobf (res-norm-aggregated sim))
    ))

(defun df-to-jacobian (sim df)
  (let* ((h (the double-float (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim))))
         (ddf (cl-mpm/fastmaths::fast-scale! (cl-mpm/fastmaths:fast-.- df (cl-mpm/utils:matrix-eye 1d0)) (/ 1d0 h))))
    (cl-mpm/fastmaths:fast-.+ (cl-mpm/utils:matrix-eye 1d0) ddf)))


(defun compute-max-deformation-2d (sim)
  (cl-mpm::reduce-over-cells
   (cl-mpm:sim-mesh sim)
   (lambda (c)
     (if (and (cl-mpm/mesh::cell-active c)
              (not (cl-mpm/mesh::cell-partial c)))
         (multiple-value-bind (l v) (cl-mpm::eig (cl-mpm/utils::slice-matrix-2d (df-to-jacobian sim (cl-mpm/mesh::cell-deformation-gradient c))))
           (abs (/ (the double-float (reduce #'max l))
                   (the double-float (reduce #'min l)))))
         0d0))
   #'max)
  )
(defun compute-max-deformation-3d (sim)
  (cl-mpm::reduce-over-cells
   (cl-mpm:sim-mesh sim)
   (lambda (c)
     (if (and (cl-mpm/mesh::cell-active c)
              (not (cl-mpm/mesh::cell-partial c)))
         (multiple-value-bind (l v) (cl-mpm::eig (df-to-jacobian sim (cl-mpm/mesh::cell-deformation-gradient c)))
           (/ (the double-float (reduce #'max l))
              (the double-float (reduce #'min l))))
         0d0))
   #'max)
  )

(defmethod compute-max-deformation (sim)
  (/ (if (= (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)) 2)
         (compute-max-deformation-2d sim)
         (compute-max-deformation-3d sim))
     (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim))
     ))

(defmethod print-max-deformation (sim)
  (let ((cmax 0d0)
        (fmax nil))
    (cl-mpm::iterate-over-cells-serial
     (cl-mpm:sim-mesh sim)
     (lambda (c)
       (if (and (cl-mpm/mesh::cell-active c)
                (not (cl-mpm/mesh::cell-partial c)))
           (let ((ct (multiple-value-bind (l v) (cl-mpm::eig (df-to-jacobian sim (cl-mpm/mesh::cell-deformation-gradient c)))
                      (/ (the double-float (reduce #'max l))
                         (the double-float (reduce #'min l))))))
             (when (> ct cmax)
               (setf cmax ct
                     fmax (df-to-jacobian sim (cl-mpm/mesh::cell-deformation-gradient c))))))))
    (format t "cmax ~E ~% f: ~A ~%" cmax fmax)
    ))

(defun criteria-deformation-gradient (sim &key (criteria 10d0))
  (%criteria-deformation-gradient sim criteria))
(defgeneric %criteria-deformation-gradient (sim criteria))
(defmethod %criteria-deformation-gradient (sim criteria)
  (>
   (compute-max-deformation sim)
   criteria))

(defun res-norm-aggregated (sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (sim-agg cl-mpm/aggregate::sim-enable-aggregate))
      sim
    (let ((res-norm 0d0))
      (incf res-norm
            (cl-mpm::reduce-over-nodes
             mesh
             (lambda (node)
               (with-accessors ((active cl-mpm/mesh::node-active)
                                (agg cl-mpm/mesh::node-agg)
                                (internal cl-mpm/mesh::node-interior)
                                (res cl-mpm/mesh::node-residual)
                                (mass cl-mpm/mesh::node-mass)
                                (res-prev cl-mpm/mesh::node-residual-prev)
                                (vel cl-mpm/mesh::node-velocity)
                                )
                   node
                 (if (and
                      (or
                       internal
                       (not agg))
                      active)
                     ;; (loop for r across (cl-mpm/utils:fast-storage res)
                     ;;       for v across (cl-mpm/utils:fast-storage vel)
                     ;;       sum (* (abs r) (abs v)))
                     ;; (abs (cl-mpm/fastmaths:dot res vel))
                     ;; (/
                     ;;  (cl-mpm/fastmaths:mag (cl-mpm/fastmaths:fast-.- res res-prev))
                     ;;  (+ 1d-10 (cl-mpm/fastmaths:mag res)))
                     (* ;mass
                        (cl-mpm/fastmaths:mag vel))
                     0d0)))
             #'+))
      ;; (when sim-agg
      ;;   (loop for d from 0 below (cl-mpm/mesh::mesh-nd mesh)
      ;;         do (let* ((res (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-residual d))
      ;;                   (vel (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-velocity d))
      ;;                   )
      ;;              (setf res (cl-mpm/aggregate::apply-internal-bcs sim res d))
      ;;              (setf vel (cl-mpm/aggregate::apply-internal-bcs sim res d))
      ;;              (incf res-norm (abs (cl-mpm/fastmaths:mag res)))
      ;;              ;; (incf res-norm (abs (cl-mpm/fastmaths:dot res vel)))
      ;;              )))
    res-norm)))

