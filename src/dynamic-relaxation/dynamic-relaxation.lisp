(defpackage :cl-mpm/dynamic-relaxation
  (:use :cl)
  (:export
   #:estimate-energy-norm
   #:estimate-power-norm
   #:estimate-oobf
   #:converge-quasi-static))
(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)

(define-condition non-convergence-error (error)
  ((text :initarg :text :reader text)
   (ke-norm :initarg :ke-norm :reader ke-norm)
   (oobf-norm :initarg :oobf-norm :reader oobf-norm)))

(defgeneric estimate-energy-norm (sim))
(defmethod estimate-energy-norm ((sim cl-mpm::mpm-sim))
  ;; (loop for mp across (cl-mpm:sim-mps sim)
  ;;           sum (* (cl-mpm/particle:mp-mass mp)
  ;;                  (cl-mpm/fastmaths::mag (cl-mpm/particle:mp-velocity mp))))
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
                  ;; (/ (cl-mpm/mesh::node-volume n) (cl-mpm/mesh::node-volume-true n))
                  (cl-mpm/mesh:node-mass n)
                  (cl-mpm/mesh::node-mass n)
                  (cl-mpm/fastmaths::mag-squared (cl-mpm/mesh::node-velocity n))
                  ))))))
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
                        (f-ext cl-mpm/mesh::node-external-force)
                        (f-int cl-mpm/mesh::node-internal-force)
                        (f-damp cl-mpm/mesh::node-damping-force)
                        (f-ghost cl-mpm/mesh::node-ghost-force)
                        (node-oobf cl-mpm/mesh::node-oobf)
                        )
           node
         (when active
           (when t;(> (cl-mpm/fastmaths::mag-squared f-ext) 0d0)
             (sb-thread:with-mutex (lock)
               (let ((inc (*
                           (expt (cl-mpm/mesh:node-mass node) 2)
                           ;; (expt (/ (cl-mpm/mesh::node-volume node) (cl-mpm/mesh::node-volume-true node)) 2)
                           (cl-mpm/fastmaths::mag-squared
                            (reduce #'cl-mpm/fastmaths::fast-.+-vector
                                    (list
                                     f-ext
                                     f-int
                                     ;; f-damp
                                     f-ghost
                                     )
                                    )))))
                 (incf node-oobf inc)
                 (setf nmax (+
                             nmax
                             inc)
                       dmax (+
                             dmax
                             (*
                              (expt (cl-mpm/mesh:node-mass node) 2)
                              ;; (expt (/ (cl-mpm/mesh::node-volume node) (cl-mpm/mesh::node-volume-true node)) 2)
                              (cl-mpm/fastmaths::mag-squared
                               f-ext)))))))))))
    (cl-mpm::iterate-over-nodes
     (cl-mpm:sim-mesh sim)
     (lambda (node)
       (with-accessors ((active cl-mpm/mesh::node-active)
                        (f-ext cl-mpm/mesh::node-external-force)
                        (f-int cl-mpm/mesh::node-internal-force)
                        (f-damp cl-mpm/mesh::node-damping-force)
                        (node-oobf cl-mpm/mesh::node-oobf)
                        )
           node
         (when active
           (when t;(> (cl-mpm/fastmaths::mag-squared f-ext) 0d0)
             (sb-thread:with-mutex (lock)
               (setf node-oobf
                     (if (> dmax 0d0)
                         (sqrt (/ node-oobf dmax))
                         ;;Very odd case, we have external force but no internal forces
                         (if (> nmax 0d0) sb-ext:double-float-positive-infinity 0d0)))
               ))))))
    (if (> dmax 0d0)
      (setf oobf (sqrt (/ nmax dmax)))
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
      (format t "Substeps ~D~%" substeps)
      (let ((full-load (list))
            (full-step (list))
            (full-energy (list))
            (energy-list (list))
            (power-last 0d0)
            (power-current 0d0)
            (energy-first 0d0)
            (energy-last 0d0)
            )
        (setf *work* 0d0)
        (loop for i from 0 to conv-steps
              while (and *run-convergance*
                         (not converged))
              do
                 (progn
                   (setf fnorm 0d0
                         load 0d0)
                   (dotimes (j substeps)
                     (setf cl-mpm/penalty::*debug-force* 0d0)
                     (cl-mpm:update-sim sim)
                     ;; (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm::calculate-min-dt sim)))
                     (setf (cl-mpm:sim-damping-factor sim) (* damping-factor (dr-estimate-damping sim)))
                     (let ((power (cl-mpm::sim-stats-power sim))
                           (energy (cl-mpm::sim-stats-energy sim)))
                       (incf *work* power)
                       (when kinetic-damping
                         (if (and
                              (> energy-last energy-first)
                              (> energy-last energy))
                                        ;(< (* power-last power) 0d0)
                             (progn
                               (format t "Peak found resetting KE~%")
                               ;; (format t "~E ~E ~E ~%" energy-first energy-last energy)
                               (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
                               (cl-mpm:iterate-over-mps
                                mps
                                (lambda (mp)
                                  (cl-mpm/fastmaths:fast-zero (cl-mpm/particle:mp-velocity mp))
                                  (cl-mpm/fastmaths:fast-zero (cl-mpm/particle::mp-acceleration mp))))
                               (setf power-last 0d0
                                     energy-first 0d0
                                     energy-last 0d0))
                             (progn
                               (setf energy-first energy-last)
                               (setf power-last power
                                     energy-last energy))))
                       ))
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
                              :oobf-norm 0d0))
        ;; (error "System didn't converge")
        ;; (format t "System didn't converge~%")
        )
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
           ;; (estimated-t 0.5d0)
           (target-time 1d-4)
           ;; (work 0d0)
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
                                  (cl-mpm/fastmaths:fast-zero (cl-mpm/particle:mp-velocity mp))
                                  (cl-mpm/fastmaths:fast-zero (cl-mpm/particle::mp-acceleration mp))))
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
                   (stats-power cl-mpm::sim-stats-power))
      sim
    (setf stats-energy (cl-mpm/dynamic-relaxation:estimate-energy-norm sim)
          stats-oobf (cl-mpm/dynamic-relaxation:estimate-oobf sim)
          stats-power (cl-mpm/dynamic-relaxation:estimate-power-norm sim))))

(defclass mpm-sim-dr (cl-mpm::mpm-sim)
  ((dt-loadstep
    :initform 1d0
    :initarg :dt-loadstep
    :accessor sim-dt-loadstep))
  (:documentation "DR psudo-linear step with update stress last update"))

(defclass mpm-sim-dr-ul (mpm-sim-dr)
  ((initial-setup
    :initform nil
    :accessor sim-initial-setup))
  (:documentation "DR psudo-linear using one displacement step"))

(defclass mpm-sim-dr-ul-usl (mpm-sim-dr)
  ()
  (:documentation "DR psudo-linear using one displacement step"))

(defclass mpm-sim-dr-usl (mpm-sim-dr)
  ()
  (:documentation "DR psudo-linear step with update stress last update"))

(defclass mpm-sim-dr-usf (mpm-sim-dr)
  ()
  (:documentation "DR psudo-linear step with update stress last update"))

(defmethod cl-mpm::finalise-loadstep ((sim mpm-sim-dr))
  ;;DR algorithm requires that finalisation is called once
  (incf (cl-mpm::sim-time sim) (sim-dt-loadstep sim))
  (cl-mpm::new-loadstep sim))

(defmethod cl-mpm::update-sim ((sim mpm-sim-dr-usl))
  "Update stress last algorithm"
  (declare (cl-mpm::mpm-sim sim))
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (dt cl-mpm::dt)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (nonlocal-damage cl-mpm::nonlocal-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (vel-algo cl-mpm::velocity-algorithm))
                sim
    (declare (type double-float mass-filter))
                (progn
                    (cl-mpm::reset-grid mesh)
                    (cl-mpm::p2g mesh mps)
                    (when (> mass-filter 0d0)
                      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
                    (cl-mpm::update-node-kinematics sim)
                    (cl-mpm::p2g-force mesh mps)
                    (loop for bcs-f in bcs-force-list
                          do (cl-mpm::apply-bcs mesh bcs-f dt))
                    ;; (apply-bcs mesh bcs-force dt)
                    (cl-mpm::update-node-forces sim)
                    (cl-mpm::apply-bcs mesh bcs dt)
                    (cl-mpm::g2p mesh mps dt vel-algo)

                    ;; (cl-mpm::update-nodes sim)
                    ;;2nd round of momentum mapping
                    (cl-mpm::reset-grid-velocity mesh)
                    (cl-mpm::p2g mesh mps)
                    (when (> mass-filter 0d0)
                      (cl-mpm::filter-grid-velocity mesh (cl-mpm::sim-mass-filter sim)))
                    (cl-mpm::update-node-kinematics sim)
                    (cl-mpm::update-nodes sim)
                    (cl-mpm::update-cells sim)
                    (cl-mpm::apply-bcs mesh bcs dt)
                    ;; (cl-mpm::update-nodes sim)
                    ;;Update stress last
                    (cl-mpm::update-stress mesh mps dt fbar)
                    (cl-mpm::update-dynamic-stats sim)
                    ;; (update-particles sim)

                    (when remove-damage
                      (cl-mpm::remove-material-damaged sim))
                    (when split
                      (cl-mpm::split-mps sim))
                    (cl-mpm::check-mps sim)
                    )))


(defun dr-estimate-damping (sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (dt cl-mpm:sim-dt))
      sim
    (let ((num 0d0)
          (denom 0d0))
      (setf num
            (cl-mpm::reduce-over-nodes
             mesh
             (lambda (node)
               (if (cl-mpm/mesh:node-active node)
                   (cl-mpm/fastmaths:dot
                    (cl-mpm/mesh:node-velocity node)
                    (cl-mpm/fastmaths:fast-.-
                     (cl-mpm/mesh::node-residual-prev node)
                     (cl-mpm/mesh::node-residual node)))
                   0d0))
             #'+))
      (setf denom
            (* dt
               (cl-mpm::reduce-over-nodes
                mesh
                (lambda (node)
                  (if (cl-mpm/mesh:node-active node)
                      (* (cl-mpm/mesh:node-mass node)
                         (cl-mpm/fastmaths:mag
                          (cl-mpm/mesh:node-velocity node)))
                      0d0))
                #'+)))
      (if (= denom 0d0)
          (cl-mpm:sim-damping-factor sim)
          (* 2d0 (sqrt (abs (/ num denom))))))))

(defmethod cl-mpm::update-sim ((sim mpm-sim-dr-usf))
  "Update stress last algorithm"
  (declare (cl-mpm::mpm-sim sim))
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (dt cl-mpm::dt)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (nonlocal-damage cl-mpm::nonlocal-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (vel-algo cl-mpm::velocity-algorithm))
                sim
    (declare (type double-float mass-filter))
                (progn
                    (cl-mpm::reset-grid mesh)
                    (cl-mpm::p2g mesh mps)
                    (when (> mass-filter 0d0)
                      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
                    (cl-mpm::update-node-kinematics sim)
                    (cl-mpm::apply-bcs mesh bcs dt)

                    (cl-mpm::update-nodes sim)

                    (cl-mpm::update-stress mesh mps dt fbar)
                    (cl-mpm::p2g-force mesh mps)
                    (cl-mpm::apply-bcs mesh bcs-force dt)
                    (loop for bcs-f in bcs-force-list
                          do (cl-mpm::apply-bcs mesh bcs-f dt))
                    ;;Update our nodes after force mapping
                    (cl-mpm::update-node-forces sim)
                    (cl-mpm::apply-bcs mesh bcs dt)
                    ;; (cl-mpm::update-nodes sim)
                    (cl-mpm::update-dynamic-stats sim)
                    (cl-mpm::g2p mesh mps dt vel-algo)
                    ;; (when remove-damage
                    ;;   (cl-mpm::remove-material-damaged sim))
                    ;; (when split
                    ;;   (cl-mpm::split-mps sim))
                    ;; (cl-mpm::check-mps sim)
                    )))
(defmethod cl-mpm::finalise-loadstep ((sim mpm-sim-dr-ul))
  ;;DR algorithm requires that finalisation is called once
  (setf (sim-initial-setup sim) nil)
  ;; (incf (cl-mpm::sim-time sim) (sim-dt-loadstep sim))
  ;; (cl-mpm::new-loadstep sim)
  (call-next-method)
  )
(defmethod cl-mpm::update-sim ((sim mpm-sim-dr-ul))
  "Update stress last algorithm"
  (declare (cl-mpm::mpm-sim sim))
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (dt cl-mpm::dt)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (nonlocal-damage cl-mpm::nonlocal-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (ghost-factor cl-mpm::ghost-factor)
               (initial-setup initial-setup)
               (vel-algo cl-mpm::velocity-algorithm))
                sim
    (declare (type double-float mass-filter))
    (progn
      (unless initial-setup
        (cl-mpm::reset-grid mesh)
        (cl-mpm::p2g mesh mps)
        (when (> mass-filter 0d0)
          (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
        (cl-mpm::update-node-kinematics sim)
        (cl-mpm::apply-bcs mesh bcs dt)
        (setf initial-setup t))
      ;; (cl-mpm::reset-grid mesh)
      (cl-mpm::apply-bcs mesh bcs dt)
      (cl-mpm::update-nodes sim)
      (cl-mpm::update-cells sim)

      ;; (cl-mpm::zero-grid-velocity mesh)
      (cl-mpm::iterate-over-nodes
       mesh
       (lambda (n)
         (when (cl-mpm/mesh::node-active n)
           (cl-mpm/mesh::reset-node-force n))))

      (when ghost-factor
        (cl-mpm/ghost::apply-ghost sim ghost-factor))
      ;; (update-node-fictious-mass sim)
      (cl-mpm::update-stress mesh mps dt fbar)
      (cl-mpm::p2g-force mesh mps)
      (cl-mpm::apply-bcs mesh bcs-force dt)
      (loop for bcs-f in bcs-force-list
            do (cl-mpm::apply-bcs mesh bcs-f dt))
      ;;Update our nodes after force mapping
      (cl-mpm::update-node-forces sim)
      (cl-mpm::apply-bcs mesh bcs dt)
      (cl-mpm::update-dynamic-stats sim)
      (cl-mpm::g2p mesh mps dt vel-algo)
      )))
(defmethod cl-mpm::update-sim ((sim mpm-sim-dr-ul-usl))
  "Update stress last algorithm"
  (declare (cl-mpm::mpm-sim sim))
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (dt cl-mpm::dt)
               (dt-loadstep dt-loadstep)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (nonlocal-damage cl-mpm::nonlocal-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (ghost-factor cl-mpm::ghost-factor)
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (vel-algo cl-mpm::velocity-algorithm))
                sim
    (declare (type double-float mass-filter))
                (progn
                    (cl-mpm::reset-grid mesh)
                    ;; Map momentum to grid
                    (cl-mpm::p2g mesh mps)
                    ;;Reset nodes below our mass-filter
                    (when (> mass-filter 0d0)
                      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
                    ;;Turn momentum into velocity
                    (cl-mpm::update-node-kinematics sim)
                    (cl-mpm::p2g-force mesh mps)
                    (loop for bcs-f in bcs-force-list
                          do (cl-mpm::apply-bcs mesh bcs-f dt))
                    (when ghost-factor
                      (cl-mpm/ghost::apply-ghost sim ghost-factor))
                    ;;Update our nodes after force mapping
                    (cl-mpm::update-node-forces sim)
                    (cl-mpm::apply-bcs mesh bcs dt)
                    (cl-mpm::g2p mesh mps dt vel-algo)
                    (cl-mpm::update-dynamic-stats sim)

                    ;;2nd round of momentum mapping

                    (cl-mpm::reset-grid-velocity mesh)
                    (cl-mpm::p2g mesh mps)
                    (when (> mass-filter 0d0)
                      (cl-mpm::filter-grid-velocity mesh (cl-mpm::sim-mass-filter sim)))
                    (cl-mpm::update-node-kinematics sim)
                    (cl-mpm::update-nodes sim)
                    (cl-mpm::apply-bcs mesh bcs dt)
                    (cl-mpm::update-cells sim)
                    ;;Update stress last
                    (cl-mpm::update-stress mesh mps dt-loadstep fbar)

                    )))

(defclass mpm-sim-dr-damage-usf (mpm-sim-dr cl-mpm/damage:mpm-sim-damage)
  ()
  (:documentation "DR psudo-linear step with update stress last update"))

(defclass mpm-sim-dr-damage-usl (mpm-sim-dr cl-mpm/damage::mpm-sim-usl-damage)
  ()
  (:documentation "DR psudo-linear step with update stress last update"))


(defun update-node-fictious-mass (sim)
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (dt cl-mpm::sim-dt))
      sim
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (setf (cl-mpm/mesh:node-mass node) 1d0))))))

(defclass mpm-sim-dr-damage-ul (mpm-sim-dr-damage-usf)
  ()
  (:documentation "DR psudo-linear step with update stress last update"))

(defclass mpm-sim-dr-damage-ul-usl (mpm-sim-dr-damage-ul)
  ()
  (:documentation "DR psudo-linear step with update stress last update"))

(defmethod cl-mpm::update-sim ((sim mpm-sim-dr-damage-UL))
  "Update stress last algorithm"
  (declare (cl-mpm::mpm-sim sim))
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (dt cl-mpm::dt)
               (dt-loadstep dt-loadstep)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (nonlocal-damage cl-mpm::nonlocal-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (ghost-factor cl-mpm::ghost-factor)
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (vel-algo cl-mpm::velocity-algorithm))
      sim
    (declare (type double-float mass-filter))
    (progn
      (cl-mpm::reset-grid mesh)
      (cl-mpm::p2g mesh mps)
      (when (> mass-filter 0d0)
        (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
      (cl-mpm::update-node-kinematics sim)
      (cl-mpm::apply-bcs mesh bcs dt)
      (cl-mpm::update-nodes sim)
      (cl-mpm::update-cells sim)
      (cl-mpm::update-stress mesh mps dt-loadstep fbar)
      (cl-mpm/damage::calculate-damage sim dt-loadstep)
      (cl-mpm::p2g-force mesh mps)
      (cl-mpm::apply-bcs mesh bcs-force dt)
      (loop for bcs-f in bcs-force-list
            do (cl-mpm::apply-bcs mesh bcs-f dt))
      (when ghost-factor
        (cl-mpm/ghost::apply-ghost sim ghost-factor))
      ;; ;;Update our nodes after force mapping
      (cl-mpm::update-node-forces sim)
      (cl-mpm::apply-bcs mesh bcs dt)

      (cl-mpm::update-dynamic-stats sim)
      (cl-mpm::g2p mesh mps dt vel-algo)
      )))

(defmethod cl-mpm::update-sim ((sim mpm-sim-dr-damage-ul-usl))
  "Update stress last algorithm"
  (declare (cl-mpm::mpm-sim sim))
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (dt cl-mpm::dt)
               (dt-loadstep dt-loadstep)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (nonlocal-damage cl-mpm::nonlocal-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (ghost-factor cl-mpm::ghost-factor)
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (vel-algo cl-mpm::velocity-algorithm))
      sim
    (declare (type double-float mass-filter))
    (progn
      (cl-mpm::reset-grid mesh)
      ;; Map momentum to grid
      (cl-mpm::p2g mesh mps)
      ;;Reset nodes below our mass-filter
      (when (> mass-filter 0d0)
        (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
      ;;Turn momentum into velocity
      (cl-mpm::update-node-kinematics sim)
      (cl-mpm::p2g-force mesh mps)
      (loop for bcs-f in bcs-force-list
            do (cl-mpm::apply-bcs mesh bcs-f dt))
      ;; (when ghost-factor
      ;;   (cl-mpm/ghost::apply-ghost sim ghost-factor))
      ;;Update our nodes after force mapping
      (cl-mpm::update-node-forces sim)
      (cl-mpm::apply-bcs mesh bcs dt)
      (cl-mpm::g2p mesh mps dt vel-algo)

      ;;2nd round of momentum mapping

      (cl-mpm::reset-grid-velocity mesh)
      (cl-mpm::p2g mesh mps)
      (when (> mass-filter 0d0)
        (cl-mpm::filter-grid-velocity mesh (cl-mpm::sim-mass-filter sim)))
      (cl-mpm::update-node-kinematics sim)
      (cl-mpm::update-nodes sim)
      (cl-mpm::apply-bcs mesh bcs dt)
      (cl-mpm::update-cells sim)
      (cl-mpm::update-dynamic-stats sim)
      ;;Update stress last
      (cl-mpm::update-stress mesh mps dt-loadstep fbar)
      (cl-mpm/damage::calculate-damage sim dt-loadstep)
      )))


(defmethod cl-mpm::finalise-loadstep ((sim mpm-sim-dr-damage-UL))
  (incf (cl-mpm::sim-time sim) (sim-dt-loadstep sim))
  (cl-mpm::new-loadstep sim))

(defmethod cl-mpm::update-sim ((sim mpm-sim-dr-damage-usf))
  "Update stress last algorithm"
  (declare (cl-mpm::mpm-sim sim))
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (dt cl-mpm::dt)
               (dt-loadstep dt-loadstep)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (nonlocal-damage cl-mpm::nonlocal-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (vel-algo cl-mpm::velocity-algorithm))
      sim
    (declare (type double-float mass-filter))
    (progn
      (cl-mpm::reset-grid mesh)
      (when (> (length mps) 0)
        (cl-mpm::p2g mesh mps)
        (when (> mass-filter 0d0)
          (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
        (cl-mpm::update-node-kinematics sim)
        ;; (update-node-fictious-mass sim)
        (cl-mpm::apply-bcs mesh bcs dt)
        (cl-mpm::update-nodes sim)
        (cl-mpm::update-stress mesh mps dt-loadstep fbar)
        (cl-mpm/damage::calculate-damage sim dt-loadstep)

        (cl-mpm::p2g-force mesh mps)
        ;; (cl-mpm::apply-bcs mesh bcs-force dt)
        (loop for bcs-f in bcs-force-list
              do (cl-mpm::apply-bcs mesh bcs-f dt))
        ;;Update our nodes after force mapping
        (cl-mpm::update-node-forces sim)
        (cl-mpm::apply-bcs mesh bcs dt)
        (cl-mpm::update-dynamic-stats sim)
        (cl-mpm::g2p mesh mps dt vel-algo)
        ;;
        (cl-mpm::update-particles sim)
        (cl-mpm::reset-node-displacement sim)

        (when remove-damage
          (cl-mpm::remove-material-damaged sim))
        (when split
          (cl-mpm::split-mps sim))
        (cl-mpm::check-mps sim))
      )))

;; (defmethod cl-mpm::finalise-loadstep ((sim mpm-sim-dr-damage-usf))
;;   (incf (cl-mpm::sim-time sim) (sim-dt-loadstep sim))
;;   (cl-mpm:iterate-over-mps
;;    (cl-mpm:sim-mps sim)
;;    (lambda (mp)
;;      (cl-mpm/particle::new-loadstep-mp mp))))


(defun reset-mp-velocity (sim)
  (cl-mpm:iterate-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (cl-mpm/fastmaths::fast-zero (cl-mpm/particle:mp-velocity mp)))))
