(defpackage :cl-mpm/dynamic-relaxation
  (:use :cl)
  (:export
   #:estimate-energy-norm
   #:estimate-power-norm
   #:estimate-oobf
   #:converge-quasi-static))
(in-package :cl-mpm/dynamic-relaxation)

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

(defgeneric estimate-oobf (sim))
(defmethod estimate-oobf (sim)
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
                        (f-int cl-mpm/mesh::node-internal-force))
           node
         (when active
           (when t;(> (cl-mpm/fastmaths::mag-squared f-ext) 0d0)
             (sb-thread:with-mutex (lock)
               ;(setf oobf-norm
               ;      (+
               ;       oobf-norm
               ;       (*
               ;        (cl-mpm/mesh:node-mass node)
               ;        ;; (/ (cl-mpm/mesh::node-volume node) (cl-mpm/mesh::node-volume-true node))
               ;        (/
               ;         (cl-mpm/fastmaths::mag-squared
               ;          (magicl:.+ f-ext f-int))
               ;         (cl-mpm/fastmaths::mag-squared f-ext)))))

               (setf nmax (+
                           nmax
                           (*
                            (expt (cl-mpm/mesh:node-mass node) 2)
                            ;; (expt (/ (cl-mpm/mesh::node-volume node) (cl-mpm/mesh::node-volume-true node)) 2)
                            (cl-mpm/fastmaths::mag-squared
                             (cl-mpm/fastmaths::fast-.+-vector f-ext f-int))))
                     dmax (+
                           dmax
                           (*
                            (expt (cl-mpm/mesh:node-mass node) 2)
                            ;; (expt (/ (cl-mpm/mesh::node-volume node) (cl-mpm/mesh::node-volume-true node)) 2)
                            (cl-mpm/fastmaths::mag-squared
                             f-ext))))))))))
    (when (> dmax 0d0)
      (setf oobf (sqrt (/ nmax dmax)))
      ;; (setf oobf (/ nmax dmax))
      )
    ;; (setf oobf (/ oobf-norm (lparallel:pmap-reduce #'cl-mpm/particle:mp-mass #'+ (cl-mpm:sim-mps sim))))
    ;; (setf oobf oobf-norm)
    oobf))
(defmethod estimate-oobf ((sim cl-mpm/mpi::mpm-sim-mpi))
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
                        (f-int cl-mpm/mesh::node-internal-force))
           node
         (when active
           (when (cl-mpm/mpi::in-computational-domain sim (cl-mpm/mesh::node-position node))
             (when t;(> (cl-mpm/fastmaths::mag-squared f-ext) 0d0)
               (sb-thread:with-mutex (lock)
                 (setf nmax (+
                             nmax
                             (*
                              (expt (cl-mpm/mesh:node-mass node) 2)
                              ;; (/ (cl-mpm/mesh::node-volume node) (cl-mpm/mesh::node-volume-true node))
                              (cl-mpm/fastmaths::mag-squared
                               (cl-mpm/fastmaths::fast-.+-vector f-ext f-int))))
                       dmax (+
                             dmax
                             (*
                              ;; (/ (cl-mpm/mesh::node-volume node) (cl-mpm/mesh::node-volume-true node))
                              (expt (cl-mpm/mesh:node-mass node) 2)
                              (cl-mpm/fastmaths::mag-squared f-ext)))))))))))

    (setf nmax (cl-mpm/mpi::mpi-sum nmax))
    (setf dmax (cl-mpm/mpi::mpi-sum dmax))
    (if (> dmax 0d0)
        (setf oobf (sqrt (/ nmax dmax)))
        ;; (setf oobf (/ nmax dmax))
      ;;Odd case where we have no forces?
      (setf oobf sb-ext:double-float-negative-infinity))

    ;; (let ((mass-total (cl-mpm/mpi::mpi-sum
    ;;                    (lparallel:pmap-reduce #'cl-mpm/particle:mp-mass #'+ (cl-mpm:sim-mps sim))))
    ;;       (oobf-norm (cl-mpm/mpi::mpi-sum oobf-norm)))
    ;;   (setf oobf (/ oobf-norm mass-total)))

    ;;Sanity check the floating point errors
    (setf oobf (cl-mpm/mpi::mpi-max oobf))
    ;; (setf oobf (cl-mpm/mpi::mpi-max oobf))
    oobf))

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
                                    )
  "Converge a simulation to a quasi-static solution via dynamic relaxation
   This is controlled by a kinetic energy norm"
  (let ((current-vel (cl-mpm::sim-velocity-algorithm sim)))
    (when pic-update
      (setf (cl-mpm::sim-velocity-algorithm sim) :BLEND))
    (restart-case (%converge-quasi-static sim energy-crit oobf-crit live-plot dt-scale substeps conv-steps post-iter-step convergance-criteria kinetic-damping)
      (continue ())
      (retry-convergence ()
        (%converge-quasi-static sim energy-crit oobf-crit live-plot dt-scale substeps conv-steps post-iter-step convergance-criteria kinetic-damping)
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
                     (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm::calculate-min-dt sim)))
                     (let ((power (estimate-power-norm sim)))
                       (incf *work* power)
                       (when kinetic-damping
                         (if (< (* power-last power) 0d0)
                             (progn
                               (format t "Peak found resetting KE~%")
                               (cl-mpm:iterate-over-mps
                                mps
                                (lambda (mp)
                                  (cl-mpm/fastmaths:fast-zero (cl-mpm/particle:mp-velocity mp))
                                  (cl-mpm/fastmaths:fast-zero (cl-mpm/particle::mp-acceleration mp))))
                               ;; (setf *work* 0d0)
                               (setf power-last 0d0))
                             (setf power-last power)))
                       ))
                   (setf load cl-mpm/penalty::*debug-force*)
                   (setf energy-total (estimate-energy-norm sim))
                   (if (= *work* 0d0)
                       (setf fnorm 0d0)
                       (setf fnorm (abs (/ energy-total *work*))))
                   (setf oobf (estimate-oobf sim))
                   (format t "Estimated dt ~E~%" (cl-mpm:sim-dt sim))
                   (format t "Conv step ~D - KE norm: ~E - Work: ~E - OOBF: ~E - Load: ~E~%" i fnorm *work* oobf
                           load)
                   ;; (when kinetic-damping
                   ;;   (push energy-total energy-list)
                   ;;   (when (> (length energy-list) 2)
                   ;;     (when (and
                   ;;            (< (nth 0 energy-list) (nth 1 energy-list))
                   ;;            (> (nth 1 energy-list) (nth 2 energy-list))
                   ;;                      ;(> (nth 2 energy-list) (nth 3 energy-list))
                   ;;            )
                   ;;       (format t "Peak found resetting KE~%")
                   ;;       (cl-mpm:iterate-over-mps
                   ;;        mps
                   ;;        (lambda (mp)
                   ;;          (cl-mpm/fastmaths:fast-zero (cl-mpm/particle:mp-velocity mp))
                   ;;          (cl-mpm/fastmaths:fast-zero (cl-mpm/particle::mp-acceleration mp))))
                   ;;       )))
                   (when (and
                          (if convergance-criteria (funcall convergance-criteria sim) t)
                          (< fnorm energy-crit)
                          (< oobf oobf-crit))
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
                                   kinetic-damping)
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
                     (let ((power (estimate-power-norm sim)))
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
                   (setf energy-total (estimate-energy-norm sim))
                   (if (= *work* 0d0)
                       (setf fnorm 0d0)
                       (setf fnorm (abs (/ energy-total *work*))))

                   (setf oobf (estimate-oobf sim))
                   (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm::calculate-min-dt sim)))

                   (when (= 0 rank)
                     (format t "Estimated dt ~E~%" (cl-mpm:sim-dt sim)))


                   (when (= 0 rank)
                     (format t "Conv step ~D - KE norm: ~E - Work: ~E - OOBF: ~E - Load: ~E~%" i fnorm *work* oobf load))

                   (push energy-total energy-list)
                   ;; (when (> (length energy-list) 2)
                   ;;   (when (and
                   ;;          (< (nth 0 energy-list) (nth 1 energy-list))
                   ;;          (> (nth 1 energy-list) (nth 2 energy-list))
                   ;;                      ;(> (nth 2 energy-list) (nth 3 energy-list))
                   ;;          )
                   ;;     (when (= rank 0)
                   ;;       (format t "Peak found resetting KE~%"))
                   ;;     (cl-mpm:iterate-over-mps
                   ;;      mps
                   ;;      (lambda (mp)
                   ;;        (cl-mpm/fastmaths:fast-zero (cl-mpm/particle:mp-velocity mp))
                   ;;        (cl-mpm/fastmaths:fast-zero (cl-mpm/particle::mp-acceleration mp))))))

                   (when (and (< fnorm energy-crit)
                              (< oobf oobf-crit)
                              (if convergance-criteria (funcall convergance-criteria sim) t)
                              )
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
