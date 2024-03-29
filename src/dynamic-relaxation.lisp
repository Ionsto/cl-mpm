(defpackage :cl-mpm/dynamic-relaxation
  (:use :cl))
(in-package :cl-mpm/dynamic-relaxation)
(defgeneric estimate-energy-norm (sim))
(defmethod estimate-energy-norm ((sim cl-mpm::mpm-sim))
  ;; (loop for mp across (cl-mpm:sim-mps sim)
  ;;           sum (* (cl-mpm/particle:mp-mass mp)
  ;;                  (cl-mpm/fastmath::mag (cl-mpm/particle:mp-velocity mp))))
  (let ((energy 0d0))
    (cl-mpm:iterate-over-nodes-serial
     (cl-mpm:sim-mesh sim)
     (lambda (n)
       (when (cl-mpm/mesh:node-active n)
         (incf energy
               (*
                (cl-mpm/mesh::node-mass n)
                (cl-mpm/fastmath::mag (cl-mpm/mesh::node-velocity n))
                )))))
    energy))
(defmethod estimate-energy-norm ((sim cl-mpm/mpi::mpm-sim-mpi))
  (cl-mpm/mpi::mpi-sum
   (let ((energy 0d0))
     (cl-mpm:iterate-over-nodes-serial
      (cl-mpm:sim-mesh sim)
      (lambda (n)
        (when (cl-mpm/mesh:node-active n)
          (when (cl-mpm/mpi::in-computational-domain sim (cl-mpm/mesh::node-position n))
            (incf energy
                  (*
                   (cl-mpm/mesh::node-mass n)
                   (cl-mpm/fastmath::mag (cl-mpm/mesh::node-velocity n))
                   ))))))
     energy)))
(defmethod estimate-oobf (sim))

(declaim (notinline plot-time-disp))
(defun plot-time-disp (full-time full-load full-energy)
  (vgplot:xlabel "Displacement (mm)")
  (vgplot:ylabel "Load (N)")
  (vgplot:plot
   (mapcar (lambda (x) (* x -1d0)) full-time) (mapcar (lambda (x) (* x 0.1)) full-load) "mps"
   (mapcar (lambda (x) (* x -1d0)) full-time) (mapcar (lambda (x) (* x 0.1)) (mapcar (lambda (x) (* x 1d4)) full-energy)) "energy"
   )
  )
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
                                     )
  "Converge a simulation to a quasi-static solution via dynamic relaxation
   This is controlled by a kinetic energy norm"
  (%converge-quasi-static sim energy-crit oobf-crit live-plot dt-scale substeps conv-steps post-iter-step))

(defgeneric %converge-quasi-static (sim
                                   energy-crit
                                   oobf-crit
                                   live-plot
                                   dt-scale
                                   substeps
                                   conv-steps
                                   post-iter-step)
 )
(defmethod %converge-quasi-static (sim
                                    energy-crit
                                    oobf-crit
                                    live-plot
                                    dt-scale
                                    substeps
                                    conv-steps
                                    post-iter-step)
  (setf *run-convergance* t)
  (let* ((fnorm 0d0)
         (oobf 0d0)
         ;; (estimated-t 0.5d0)
         (target-time 1d-4)
         ;; (substeps 40)
         ;; (substeps (floor estimated-t (cl-mpm:sim-dt sim)))
         (estimated-t 1d-5)
         ;; (substeps (floor estimated-t (cl-mpm:sim-dt sim)))
         (total-step 0)
        (converged nil))
    (format t "Substeps ~D~%" substeps)
    ;; (format t "dt ~D~%" dt)
    ;; (setf *full-load* (list)
    ;;       *full-step* (list)
    ;;       *full-energy* (list)
    ;;       )
    (let ((full-load (list))
          (full-step (list))
          (full-energy (list)))
      (loop for i from 0 to conv-steps
            while (and *run-convergance*
                   (not converged))
            do
               (progn
                 (dotimes (j substeps)
                   ;; (push
                   ;;  ;; (get-reaction-force *fixed-nodes*)
                   ;;  cl-mpm/penalty::*debug-force*
                   ;;  *full-load*)
                   ;; (push
                   ;;  (estimate-energy-norm sim)
                   ;;  *full-energy*)
                   ;; (push
                   ;;  total-step
                   ;;  *full-step*)
                   ;; (incf total-step)

                   ;; (push
                   ;;  (get-reaction-force *fixed-nodes*)
                   ;;  ;; cl-mpm/penalty::*debug-force*
                   ;;  *data-full-reaction*)
                   ;; (push
                   ;;  (estimate-energy-crit *sim*)
                   ;;  *data-full-energy*)
                   ;; (push
                   ;;  *t*
                   ;;  *data-full-time*)
                   ;; (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))
                   (setf cl-mpm/penalty::*debug-force* 0d0)
                   (cl-mpm:update-sim sim)
                   )
                 ;; (when t;live-plot
                 ;;   (plot-time-disp full-step full-load full-energy))
                 ;; (plot-time-disp)
                 (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time sim target-time :dt-scale dt-scale)
                   (format t "CFL dt estimate: ~f~%" dt-e)
                   (format t "CFL step count estimate: ~D~%" substeps-e)
                   ;; (setf substeps substeps-e)
                   )
                 (setf fnorm (estimate-energy-norm sim))

                 (setf oobf 0d0)
                 (let ((nmax 0d0)
                       (dmax 0d0)
                       (imax 0)
                       (iter 0))
                   ;; (cl-mpm::iterate-over-nodes-serial
                   ;;  (cl-mpm:sim-mesh sim)
                   ;;  (lambda (node)
                   ;;    (with-accessors ((active cl-mpm/mesh::node-active)
                   ;;                     (f-ext cl-mpm/mesh::node-external-force)
                   ;;                     (f-int cl-mpm/mesh::node-internal-force))
                   ;;        node
                   ;;      (when active
                   ;;        (setf imax iter)
                   ;;        (setf nmax (+ nmax
                   ;;                      (cl-mpm/fastmath::mag-squared
                   ;;                       (magicl:.- f-ext f-int)))
                   ;;              dmax (+ dmax (cl-mpm/fastmath::mag-squared f-ext))))
                   ;;      )
                   ;;    (incf iter)
                   ;;    ))
                   (when (> dmax 0d0)
                     (setf oobf (/ nmax dmax)))
                   (format t "Conv step ~D - KE norm: ~E - OOBF: ~E - Load: ~E~%" i fnorm oobf
                           cl-mpm/penalty::*debug-force*)
                   (when (and (< fnorm energy-crit)
                              ;(< oobf oobf-crit)
                              )
                     (format t "Took ~D steps to converge~%" i)
                     (setf converged t))
                   (when post-iter-step
                     (funcall post-iter-step)))
                 (swank.live:update-swank))))
    (when (not converged)
      (error "System didn't converge")
      ;; (format t "System didn't converge~%")
      )
    ))
(defmethod %converge-quasi-static ((sim cl-mpm/mpi:mpm-sim-mpi)
                                    energy-crit
                                    oobf-crit
                                    live-plot
                                    dt-scale
                                    substeps
                                    conv-steps
                                    post-iter-step)
  (setf *run-convergance* t)
  (let* ((fnorm 0d0)
         (oobf 0d0)
         (rank (cl-mpi:mpi-comm-rank))
         ;; (estimated-t 0.5d0)
         (target-time 1d-4)
         ;; (substeps 40)
         ;; (substeps (floor estimated-t (cl-mpm:sim-dt sim)))
         (estimated-t 1d-5)
         ;; (substeps (floor estimated-t (cl-mpm:sim-dt sim)))
         (total-step 0)
        (converged nil))
    (format t "Substeps ~D~%" substeps)
    ;; (format t "dt ~D~%" dt)
    ;; (setf *full-load* (list)
    ;;       *full-step* (list)
    ;;       *full-energy* (list)
    ;;       )
    (let ((full-load (list))
          (full-step (list))
          (full-energy (list)))
      (loop for i from 0 to conv-steps
            while (and *run-convergance*
                   (not converged))
            do
               (progn
                 (dotimes (j substeps)
                   (setf cl-mpm/penalty::*debug-force* 0d0)
                   (cl-mpm:update-sim sim)
                   )
                 (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time sim target-time :dt-scale dt-scale)
                   (when (= 0 rank)
                     (format t "CFL dt estimate: ~f~%" dt-e)
                     (format t "CFL step count estimate: ~D~%" substeps-e))
                   )
                 (setf fnorm (estimate-energy-norm sim))

                 (setf oobf 0d0)
                 (let ((nmax 0d0)
                       (dmax 0d0)
                       (imax 0)
                       (iter 0)
                       (force (cl-mpm/mpi:mpi-sum cl-mpm/penalty::*debug-force*))
                       )
                   (when (> dmax 0d0)
                     (setf oobf (/ nmax dmax)))

                   (when (= 0 rank)
                     (format t "Conv step ~D - KE norm: ~E - OOBF: ~E - Load: ~E~%" i fnorm oobf
                             force))
                   (when (and (< fnorm energy-crit)
                              ;(< oobf oobf-crit)
                              )
                     (when (= 0 rank)
                       (format t "Took ~D steps to converge~%" i))
                     (setf converged t))
                   (when post-iter-step
                     (funcall post-iter-step)))
                 (swank.live:update-swank))))
    (when (not converged)
      ;; (error "System didn't converge")
      (when (= 0 rank)
        (format t "System didn't converge~%"))
      )
    ))

;; (defun converge-quasi-static (sim &key
;;                                     (energy-crit 1d-8)
;;                                     (oobf-crit 1d-8)
;;                                     (live-plot nil)
;;                                     (dt-scale 0.5d0)
;;                                     (substeps 50)
;;                                     (conv-steps 50)
;;                                     (post-iter-step nil)
;;                                     )
;;   (dotimes (i substeps)
;;     (cl-mpm:update-sim sim))
;;   (setf cl-mpm/penalty::*debug-force* 0d0)
;;   (cl-mpm:update-sim sim)
;;   )

;; (defun quasi-test ()
;;   ;; (setup)
;;   ;; (setf *run-sim* t)
;;   ;; (incf *target-displacement* (* -0.002d-3 4))
;;   (setf cl-mpm/penalty::*debug-force* 0d0)
;;   (setf cl-mpm/penalty::*debug-force-count* 0d0)
;;   (setf cl-mpm/damage::*enable-reflect-x* t)
;;   (defparameter *data-full-time* '(0d0))
;;   (defparameter *data-full-load* '(0d0))
;;   (defparameter *data-full-reaction* '(0d0))
;;   (time
;;    (converge-quasi-static *sim*))
;;   )

