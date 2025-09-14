(defpackage :cl-mpm/dynamic-relaxation
  (:use :cl)
  (:export
   #:estimate-energy-norm
   #:estimate-power-norm
   #:estimate-oobf
   #:converge-quasi-static))
(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)

(defclass mpm-sim-dr (cl-mpm::mpm-sim)
  ((dt-loadstep
    :initform 1d0
    :initarg :dt-loadstep
    :accessor sim-dt-loadstep)
   (damping-scale
    :initform 1d0
    :initarg :damping-scale
    :accessor sim-damping-scale))
  (:default-initargs
   :vel-algo :QUASI-STATIC)
  (:documentation "DR psudo-linear step with update stress last update"))

(defclass mpm-sim-dr-ul (mpm-sim-dr cl-mpm/aggregate::mpm-sim-aggregated)
  ((initial-setup
    :initform nil
    :accessor sim-initial-setup))
  (:documentation "DR solved quasi-static implicit mpm"))

(defclass mpm-sim-implict-dynamic (mpm-sim-dr-ul)
  ()
  (:documentation "DR implicit-dynamic mpm"))

(defclass mpm-sim-dr-usf (mpm-sim-dr)
  ()
  (:documentation "DR psudo-linear step with update stress last update"))

(defclass mpm-sim-dr-damage-usf (mpm-sim-dr cl-mpm/damage:mpm-sim-damage)
  ()
  (:documentation "DR psudo-linear step with update stress last update"))

(defclass mpm-sim-dr-damage-ul (mpm-sim-dr-ul cl-mpm/damage:mpm-sim-damage)
  ()
  (:documentation "DR psudo-linear step with update stress last update"))


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
                  (or
                   (not (cl-mpm/mesh::node-agg node))
                   (cl-mpm/mesh::node-interior node)))
             (with-accessors ((active cl-mpm/mesh::node-active)
                              (f-ext cl-mpm/mesh::node-external-force)
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
               (let (;(mass 1d0)
                     ;; (scale-factor (expt mass 1))
                     (scale-factor 1d0)
                     ;; (scale-factor (/ volume volume-t))
                     )
                 (list
                  scale-factor
                  (* scale-factor (* 0.5d0 mass (cl-mpm/fastmaths::mag-squared vel)))
                  (* scale-factor (cl-mpm/fastmaths::mag (cl-mpm/fastmaths::fast-.+-vector f-ext f-int)))
                  (* scale-factor (cl-mpm/fastmaths::mag f-ext))
                  (* scale-factor
                     (cl-mpm/fastmaths:dot
                      disp f-ext)))))
             (list 0d0 0d0 0d0 0d0 0d0)))
       (lambda (a b) (mapcar (lambda (x y) (declare (double-float x y)) (+ x y)) a b)))
    (declare (double-float mass energy oobf-num oobf-denom power))
    ;; (break)
    (let ((oobf 0d0))
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
                            (f-int cl-mpm/mesh::node-internal-force)
                            (n-mass cl-mpm/mesh::node-mass)
                            (node-oobf cl-mpm/mesh::node-oobf))
               node
             (when (and active
                        (or
                         (not agg) (cl-mpm/mesh::node-interior node)))
               (when t
                 (setf node-oobf
                       (if (> oobf 0d0)
                            (/ (* n-mass (cl-mpm/fastmaths::mag (cl-mpm/fastmaths::fast-.+-vector f-ext f-int)))
                               oobf-denom)
                           0d0
                           ))))))))
      (values (/ energy mass) oobf (/ power mass)))))

(defun combi-stats-aggregated (sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (sim-agg cl-mpm/aggregate::sim-enable-aggregate))
      sim
    (destructuring-bind (mass
                         energy
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
                 (let ((scale-factor 1d0))
                   (list
                    scale-factor
                    (* scale-factor (* 0.5d0 mass (cl-mpm/fastmaths::mag-squared vel)))
                    (* scale-factor (cl-mpm/fastmaths::mag-squared (cl-mpm/fastmaths::fast-.+-vector f-ext f-int)))
                    (* scale-factor (cl-mpm/fastmaths::mag-squared f-ext))
                    (* scale-factor
                       (cl-mpm/fastmaths:dot
                        disp f-ext)))))
               (list 0d0 0d0 0d0 0d0 0d0)))
         (lambda (a b) (mapcar (lambda (x y) (declare (double-float x y)) (+ x y)) a b)))
      (declare (double-float mass energy oobf-num oobf-denom power))
      (setf mass 1d0)
      (when sim-agg
        (loop for d from 0 below (cl-mpm/mesh::mesh-nd mesh)
              do (let* ((f-int (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-internal-force d))
                        (f-ext (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-external-force d))
                        (E (cl-mpm/aggregate::sim-global-e sim))
                        (ma (cl-mpm/aggregate::sim-global-ma sim))
                        (vi (magicl:@ (magicl:transpose E) (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-velocity d)))
                        (disp (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-displacment d)))
                   (incf oobf-num (cl-mpm/fastmaths::mag-squared
                                   (cl-mpm/aggregate::apply-internal-bcs
                                    sim
                                    (magicl:@ (magicl:transpose E)
                                              (cl-mpm/fastmaths::fast-.+
                                               f-int
                                               f-ext))
                                    d
                                    )))
                   (incf oobf-denom (cl-mpm/fastmaths::mag-squared
                                     (cl-mpm/aggregate::apply-internal-bcs
                                      sim
                                      (magicl:@ (magicl:transpose E) f-ext)
                                      d)))
                   (incf power (cl-mpm/fastmaths:dot
                                disp f-ext
                                ;; (magicl:@ (magicl:transpose E) disp)
                                ;; (magicl:@ (magicl:transpose E) f-ext)
                                ))
                   (incf energy (* 0.5d0 (cl-mpm/utils::mtref (magicl:@ (magicl:transpose vi) ma vi) 0 0)))
                   )))
      (let ((oobf 0d0)
            (oobf-num (sqrt oobf-num))
            (oobf-denom (sqrt oobf-denom))
            )
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
                              (f-int cl-mpm/mesh::node-internal-force)
                              (n-mass cl-mpm/mesh::node-mass)
                              (node-oobf cl-mpm/mesh::node-oobf))
                 node
               (when (and active
                          (not agg))
                 (when t
                   (setf node-oobf
                         (if (> oobf 0d0)
                             (/ (* n-mass (cl-mpm/fastmaths::mag (cl-mpm/fastmaths::fast-.+-vector f-ext f-int)))
                                oobf-denom)
                             0d0
                             ))))))))
        ;; (values 0d0 oobf 1d0)
        (values energy oobf power)
        ))))

(define-condition non-convergence-error (error)
  ((text :initarg :text :reader text)
   (ke-norm :initarg :ke-norm :reader ke-norm)
   (oobf-norm :initarg :oobf-norm :reader oobf-norm)))


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

(defgeneric estimate-oobf (sim))
(defmethod estimate-oobf (sim)
  "With better reporting"
  )
(defmethod estimate-oobf (sim)
  (estimate-oobf-debug sim))


(declaim (notinline plot-time-disp))
;; (defun plot-time-disp (full-time full-load full-energy)
;;   (vgplot:xlabel "Displacement (mm)")
;;   (vgplot:ylabel "Load (N)")
;;   (vgplot:plot
;;    (mapcar (lambda (x) (* x -1d0)) full-time) (mapcar (lambda (x) (* x 0.1)) full-load) "mps"
;;    (mapcar (lambda (x) (* x -1d0)) full-time) (mapcar (lambda (x) (* x 0.1)) (mapcar (lambda (x) (* x 1d4)) full-energy)) "energy"
;;    )
;;   )
;; (plot-time-disp *full-step* *full-load* *full-energy*)

;; (defparameter *full-load* nil)
;; (defparameter *full-step* nil)
;; (defparameter *full-energy* nil)
(defparameter *run-convergance* t)
;; (defparameter *run-convergance* nil)
(declaim (notinline converge-quasi-static))
(defun converge-quasi-static (sim &key
                                    (energy-crit 1d-8)
                                    (oobf-crit 1d-8)
                                    (live-plot nil)
                                    (dt-scale 0.5d0)
                                    (substeps 50)
                                    (conv-steps 50)
                                    (post-iter-step nil)
                                    (convergance-criteria nil)
                                    (pic-update t)
                                    (kinetic-damping t)
                                    (damping-factor 1d-1)
                                    )
  "Converge a simulation to a quasi-static solution via dynamic relaxation
   This is controlled by a kinetic energy norm"
  (let ((current-vel (cl-mpm::sim-velocity-algorithm sim)))
    (when pic-update
      (setf (cl-mpm::sim-velocity-algorithm sim) :BLEND))
    (when (typep sim 'mpm-sim-dr)
      (setf (cl-mpm/dynamic-relaxation::sim-damping-scale sim)
            (if damping-factor
                damping-factor
                0d0)))
    (restart-case
        (%converge-quasi-static sim energy-crit oobf-crit live-plot dt-scale substeps conv-steps post-iter-step convergance-criteria kinetic-damping damping-factor)
      (continue ())
      (retry-convergence ()
        (%converge-quasi-static sim energy-crit oobf-crit live-plot dt-scale substeps conv-steps post-iter-step convergance-criteria kinetic-damping damping-factor)
        ))
    (setf (cl-mpm::sim-velocity-algorithm sim) current-vel)))

(defgeneric %converge-quasi-static (sim
                                    energy-crit
                                    oobf-crit
                                    live-plot
                                    dt-scale
                                    substeps
                                    conv-steps
                                    post-iter-step
                                    convergance-criteria
                                    kinetic-damping
                                    damping-factor
                                    )
 )
(defparameter *work* 0d0)
(defmethod %converge-quasi-static (sim
                                   energy-crit
                                   oobf-crit
                                   live-plot
                                   dt-scale
                                   substeps
                                   conv-steps
                                   post-iter-step
                                   convergance-criteria
                                   kinetic-damping
                                   damping-factor
                                   )
  (declare (double-float dt-scale))
  (setf *run-convergance* t)
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (let* ((fnorm 0d0)
           (energy-total 0d0)
           (oobf 0d0)
           (target-time 1d-4)
           (estimated-t 1d-5)
           (total-step 0)
           (load 0d0)
           (converged nil))

      (setf (cl-mpm:sim-dt sim) (cl-mpm/setup::estimate-elastic-dt sim :dt-scale dt-scale))
      (if damping-factor
          (setf (cl-mpm:sim-damping-factor sim)
              (* damping-factor (cl-mpm/setup::estimate-critical-damping sim)))
          (setf (cl-mpm:sim-damping-factor sim) 0d0))

      (format t "Substeps ~D~%" substeps)
      (let ((full-load (list))
            (full-step (list))
            (full-enegy (list))
            (energy-list (list))
            (power-last 0d0)
            (power-current 0d0)
            (energy-first 0d0)
            (energy-last 0d0)
            ;; (dt-0 (cl-mpm/setup::estimate-elastic-dt sim :dt-scale dt-scale))
            )
        (setf *work* 0d0)
        (setf (cl-mpm::sim-stats-work sim) 0d0)
        (loop for i from 0 to conv-steps
              while (and *run-convergance*
                         (cl-mpm::sim-run-sim sim)
                         (not converged))
              do
                 (progn
                   (setf fnorm 0d0
                         load 0d0)
                   (time
                    (dotimes (j substeps)
                      (setf cl-mpm/penalty::*debug-force* 0d0)
                      (cl-mpm:update-sim sim)
                      ;; (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm::calculate-min-dt sim)))
                      ;; (format t "dt adjustment ~E~%" (/ (cl-mpm:sim-dt sim) dt-0))
                      (when damping-factor
                        (setf (cl-mpm:sim-damping-factor sim) (* damping-factor (dr-estimate-damping sim))))
                      (let ((power (cl-mpm::sim-stats-power sim))
                            (energy (cl-mpm::sim-stats-energy sim)))
                        (incf *work* power)
                        (when kinetic-damping
                          (if (and
                               ;; (< (* power power-last) 0d0)
                               (> energy-last energy-first)
                               (> energy-last energy))
                              (progn
                                ;; (format t "Peak found resetting KE - ~E ~E ~E~%" energy-first energy-last energy)
                                (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
                                (cl-mpm:iterate-over-mps
                                 mps
                                 (lambda (mp)
                                   (cl-mpm/fastmaths:fast-zero (cl-mpm/particle:mp-velocity mp))))
                                (setf power-last 0d0
                                      energy-first 0d0
                                      energy-last 0d0
                                      energy 0d0))
                              (progn
                                (setf energy-first energy-last)
                                (setf power-last power
                                      energy-last energy)))))))
                   (setf load cl-mpm/penalty::*debug-force*)
                   (setf energy-total (cl-mpm::sim-stats-energy sim))
                   (if (= *work* 0d0)
                       (setf fnorm 0d0)
                       (setf fnorm (abs (/ energy-total *work*))))
                   (setf oobf (cl-mpm::sim-stats-oobf sim))
                   (format t "Estimated dt ~E~%" (cl-mpm:sim-dt sim))
                   (format t "Conv step ~D - KE norm: ~E - Work: ~E - OOBF: ~E - Load: ~E~%" i fnorm *work* oobf
                           load)
                   (when (if convergance-criteria
                             (funcall convergance-criteria sim fnorm oobf)
                             (and
                              (< fnorm energy-crit)
                              (< oobf oobf-crit)))
                     (format t "Took ~D steps to converge~%" i)
                     (setf converged t))
                   (when post-iter-step
                     (funcall post-iter-step i fnorm oobf))
                   (swank.live:update-swank))))
      (when (not converged)
        (error (make-instance 'non-convergence-error
                              :text "System failed to converge"
                              :ke-norm fnorm
                              :oobf-norm 0d0)))
      (values load fnorm oobf))))
(defmethod %converge-quasi-static ((sim cl-mpm/mpi:mpm-sim-mpi)
                                   energy-crit
                                   oobf-crit
                                   live-plot
                                   dt-scale
                                   substeps
                                   conv-steps
                                   post-iter-step
                                   convergance-criteria
                                   kinetic-damping
                                    damping-factor
                                   )
  (setf *run-convergance* t)
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (let* ((fnorm 0d0)
           (oobf 0d0)
           (rank (cl-mpi:mpi-comm-rank))
           (load 0d0)
           (target-time 1d-4)
           (converged nil))
      (setf (cl-mpm:sim-dt sim)
            (cl-mpm/setup::estimate-elastic-dt sim :dt-scale dt-scale))
      (when (= rank 0)
        (format t "Substeps ~D~%" substeps))
      (setf *work* 0d0)
      (let ((full-load (list))
            (full-step (list))
            (energy-list (list))
            (energy-total 0d0)
            (power-last 0d0))
        (loop for i from 0 to conv-steps
              while (and *run-convergance*
                         (not converged))
              do
                 (progn
                   (dotimes (j substeps)
                     (setf cl-mpm/penalty::*debug-force* 0d0)
                     (cl-mpm:update-sim sim)
                     (let ((power (cl-mpm::sim-stats-power sim)))
                       (incf *work* power)
                       (when kinetic-damping
                         (if (< (* power-last power) 0d0)
                             (progn
                               (when (= rank 0)
                                 (format t "Peak found resetting KE~%"))
                               (cl-mpm:iterate-over-mps
                                mps
                                (lambda (mp)
                                  (cl-mpm/fastmaths:fast-zero (cl-mpm/particle:mp-velocity mp))))
                               (setf power-last 0d0))
                             (setf power-last power)))))

                   (setf load (cl-mpm/mpi::mpi-sum cl-mpm/penalty::*debug-force*))
                   (setf energy-total (cl-mpm::sim-stats-energy sim))
                   (if (= *work* 0d0)
                       (setf fnorm 0d0)
                       (setf fnorm (abs (/ energy-total *work*))))

                   (setf oobf (cl-mpm::sim-stats-oobf sim))
                   (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm::calculate-min-dt sim)))

                   (when (= 0 rank)
                     (format t "Estimated dt ~E~%" (cl-mpm:sim-dt sim)))


                   (when (= 0 rank)
                     (format t "Conv step ~D - KE norm: ~E - Work: ~E - OOBF: ~E - Load: ~E~%" i fnorm *work* oobf load))

                   (push energy-total energy-list)
                   (when (and (< fnorm energy-crit)
                              (< oobf oobf-crit)
                              (if convergance-criteria (funcall convergance-criteria sim) t))
                     (when (= 0 rank)
                       (format t "Took ~D steps to converge~%" i))
                     (setf converged t))
                   (when post-iter-step
                     (funcall post-iter-step i fnorm oobf))
                   (swank.live:update-swank))))
      (when (not converged)
        (error (make-instance 'non-convergence-error
                              :text "System failed to converge"
                              :ke-norm fnorm
                              :oobf-norm 0d0))
        (when (= 0 rank)
          (format t "System didn't converge~%"))
        )
      (values load fnorm oobf))))

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
    (multiple-value-bind (e o p) (cl-mpm/dynamic-relaxation::combi-stats-aggregated sim)
      (setf stats-energy e
            stats-oobf o
            stats-power p)
      (incf stats-work p))))


;; (defun dr-estimate-damping (sim)
;;   (with-accessors ((mesh cl-mpm:sim-mesh)
;;                    (dt cl-mpm:sim-dt))
;;       sim
;;     (let ((num 0d0)
;;           (denom 0d0))
;;       (setf num
;;             (cl-mpm::reduce-over-nodes
;;              mesh
;;              (lambda (node)
;;                (if (and (cl-mpm/mesh:node-active node)
;;                         (not (cl-mpm/mesh::node-agg node)))
;;                    (/ 
;;                     (cl-mpm/fastmaths:dot
;;                      (cl-mpm/mesh:node-velocity node)
;;                      ;; (cl-mpm/mesh::node-internal-force node)
;;                      (cl-mpm/fastmaths:fast-.-
;;                       (cl-mpm/mesh::node-residual node)
;;                       (cl-mpm/mesh::node-residual-prev node)
;;                       )
;;                      )
;;                     dt
;;                     )
;;                    0d0))
;;              #'+))
;;       (setf denom
;;             (* ;dt
;;                (cl-mpm::reduce-over-nodes
;;                 mesh
;;                 (lambda (node)
;;                   (if (and (cl-mpm/mesh:node-active node)
;;                            (not (cl-mpm/mesh::node-agg node)))
;;                       (* (cl-mpm/mesh:node-mass node)
;;                          (cl-mpm/fastmaths::dot
;;                           (cl-mpm/mesh::node-displacment node)
;;                           (cl-mpm/mesh::node-displacment node)))
;;                       0d0))
;;                 #'+)))
;;       (if (= denom 0d0)
;;           (cl-mpm/setup::estimate-critical-damping sim)
;;           (* 2d0 (sqrt (/ (max 0d0 num) denom)))))))

(defun dr-estimate-damping (sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (dt cl-mpm:sim-dt))
      sim
    (let ((num 0d0)
          (denom 0d0))
      (setf
       num
       (cl-mpm::reduce-over-nodes
        mesh
        (lambda (node)
          (if (and (cl-mpm/mesh:node-active node)
                   (not (cl-mpm/mesh::node-agg node)))
              (cl-mpm/fastmaths:dot
               (cl-mpm/mesh:node-velocity node)
               (cl-mpm/fastmaths:fast-.-
                (cl-mpm/mesh::node-residual-prev node)
                (cl-mpm/mesh::node-residual node)))
              0d0))
        #'+))
      (setf
       denom
       (* dt
          (cl-mpm::reduce-over-nodes
           mesh
           (lambda (node)
             (if (and (cl-mpm/mesh:node-active node)
                      (not (cl-mpm/mesh::node-agg node)))
                 (* (cl-mpm/mesh:node-mass node)
                    (cl-mpm/fastmaths::dot
                     (cl-mpm/mesh::node-velocity node)
                     (cl-mpm/mesh::node-velocity node)))
                 0d0))
           #'+)))

      (when (cl-mpm/aggregate::sim-enable-aggregate sim)
        (loop for d from 0 below (cl-mpm/mesh::mesh-nd mesh)
              do
                 (let* ((res (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-residual d))
                        (res-prev (cl-mpm/aggregate::assemble-global-vec sim #'cl-mpm/mesh::node-residual-prev d))
                        (E (cl-mpm/aggregate::sim-global-e sim))
                        (ma (cl-mpm/aggregate::sim-global-ma sim))
                        (vi (cl-mpm/aggregate::assemble-internal-vec sim #'cl-mpm/mesh::node-velocity d)))
                   (incf num (cl-mpm/fastmaths:dot
                              vi
                              (magicl:@ (magicl:transpose E)
                                        (cl-mpm/fastmaths::fast-.-
                                         res-prev
                                         res))))
                   (incf denom (* 0.5d0 (cl-mpm/utils::mtref (magicl:@ (magicl:transpose vi) ma vi) 0 0))))))
      (if (> num 0d0)
          (if (= denom 0d0)
              ;; (cl-mpm/setup::estimate-critical-damping sim)
              0d0
              (* 2 (sqrt (/ num denom))))
          0d0))))


(defun reset-mp-velocity (sim)
  (cl-mpm:iterate-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (cl-mpm/fastmaths::fast-zero (cl-mpm/particle:mp-velocity mp)))))


(defmethod cl-mpm::finalise-loadstep ((sim mpm-sim-dr))
  ;;DR algorithm requires that finalisation is called once
  (incf (cl-mpm::sim-time sim) (sim-dt-loadstep sim))
  (cl-mpm::new-loadstep sim))
