(defpackage :cl-mpm/examples/collapse
  (:use :cl))
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
;; (setf *block-compile-default* t)
;(sb-int:set-floating-point-modes :traps '(:overflow :invalid :inexact :divide-by-zero :underflow))
;; (sb-int:set-floating-point-modes :traps '(:overflow :divide-by-zero :underflow))

(in-package :cl-mpm/examples/collapse)
(declaim (optimize (debug 3) (safety 3) (speed 0)))

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle) dt fbar)
  (cl-mpm::update-stress-kirchoff mesh mp dt fbar)
  ;; (cl-mpm::update-domain-polar-2d mesh mp dt)
  ;; (cl-mpm::scale-domain-size mesh mp)
  )

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-finite-viscoelastic) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  (cl-mpm::update-domain-polar-2d mesh mp dt)
  (cl-mpm::scale-domain-size mesh mp)
  )
(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-vm) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  (cl-mpm::update-domain-polar-2d mesh mp dt)
  ;; (cl-mpm::update-domain-midpoint mesh mp dt)
  ;(cl-mpm::update-domain-max-corner-2d mesh mp dt)
  ;; (cl-mpm::scale-domain-size mesh mp)
  )


(defmethod cl-mpm::update-dynamic-stats ((sim cl-mpm::mpm-sim-usf))
  (with-accessors ((stats-energy cl-mpm::sim-stats-energy)
                   (stats-oobf cl-mpm::sim-stats-oobf)
                   (stats-power cl-mpm::sim-stats-power))
      sim
    (setf stats-energy (cl-mpm/dynamic-relaxation:estimate-energy-norm sim)
          stats-oobf (cl-mpm/dynamic-relaxation:estimate-oobf sim)
          stats-power (cl-mpm/dynamic-relaxation:estimate-power-norm sim))))

(defun plot-load-disp ()
  (vgplot:plot *data-steps* *data-energy*))
(defun plot (sim)
  ;; (plot-load-disp)
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :deformed
   ;; :colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :xy))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-damage-ybar mp))
   :colour-func (lambda (mp) (cl-mpm/particle::mp-strain-plastic-vm mp))
   )
  )
(defmethod cl-mpm/output::save-vtk (filename (sim cl-mpm/damage::mpm-sim-damage))
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)) sim
    (with-open-file (fs filename :direction :output :if-exists :supersede)
      (format fs "# vtk DataFile Version 2.0~%")
      (format fs "Lisp generated vtk file, WMC~%")
      (format fs "ASCII~%")
      (format fs "DATASET UNSTRUCTURED_GRID~%")
      (format fs "POINTS ~d double~%" (length mps))
      (loop for mp across mps
            do (format fs "~E ~E ~E ~%"
                       (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 1 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 2 0) 'single-float)))
      (format fs "~%")

      ;; (with-parameter-list fs mps
      ;;   ("mass" (lambda (mp) (cl-mpm/particle:mp-mass mp))))

      (let ((id 1)
            (nd (cl-mpm/mesh:mesh-nd mesh)))
        (declare (special id))
        (format fs "POINT_DATA ~d~%" (length mps))

        (cl-mpm/output::save-parameter "mass" (cl-mpm/particle:mp-mass mp))
        (cl-mpm/output::save-parameter "density" (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)))
        (cl-mpm/output::save-parameter "index" (cl-mpm/particle::mp-index mp))
        (cl-mpm/output::save-parameter "mpi-domain" (cl-mpm/particle::mp-mpi-index mp))
        (cl-mpm/output::save-parameter "split-depth" (cl-mpm/particle::mp-split-depth mp))
        (cl-mpm/output::save-parameter "vel_x" (magicl:tref (cl-mpm/particle:mp-velocity mp) 0 0))
        (cl-mpm/output::save-parameter "vel_y" (magicl:tref (cl-mpm/particle:mp-velocity mp) 1 0))
        (cl-mpm/output::save-parameter "vel_z" (magicl:tref (cl-mpm/particle:mp-velocity mp) 2 0))

        ;; (save-parameter "acc_x" (magicl:tref (cl-mpm/particle::mp-acceleration mp) 0 0))
        ;; (save-parameter "acc_y" (magicl:tref (cl-mpm/particle::mp-acceleration mp) 1 0))

        (cl-mpm/output::save-parameter "disp_x" (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0))
        (cl-mpm/output::save-parameter "disp_y" (magicl:tref (cl-mpm/particle::mp-displacement mp) 1 0))
        (cl-mpm/output::save-parameter "disp_z" (magicl:tref (cl-mpm/particle::mp-displacement mp) 2 0))

        (cl-mpm/output::save-parameter "damage"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage)
                                           (cl-mpm/particle:mp-damage mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-k"
                                       (if (slot-exists-p mp 'cl-mpm/particle::history-stress)
                                           (cl-mpm/particle::mp-history-stress mp)
                                           0d0))
        (cl-mpm/output::save-parameter "local-length"
                                       (if (slot-exists-p mp 'cl-mpm/particle::true-local-length)
                                           (cl-mpm/particle::mp-true-local-length mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-ybar"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-ybar)
                                           (cl-mpm/particle::mp-damage-ybar mp)
                                           0d0))
        (cl-mpm/output::save-parameter "sig_xx" (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0))
        (cl-mpm/output::save-parameter "sig_yy" (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0))
        (cl-mpm/output::save-parameter "sig_xy" (magicl:tref (cl-mpm/particle:mp-stress mp) 5 0))
        (when (= nd 3)
          (cl-mpm/output::save-parameter "sig_zz" (magicl:tref (cl-mpm/particle:mp-stress mp) 2 0))
          (cl-mpm/output::save-parameter "sig_yz" (magicl:tref (cl-mpm/particle:mp-stress mp) 3 0))
          (cl-mpm/output::save-parameter "sig_zx" (magicl:tref (cl-mpm/particle:mp-stress mp) 4 0)))

        (cl-mpm/output::save-parameter "eps_xx" (magicl:tref (cl-mpm/particle:mp-strain mp) 0 0))
        (cl-mpm/output::save-parameter "eps_yy" (magicl:tref (cl-mpm/particle:mp-strain mp) 1 0))
        (cl-mpm/output::save-parameter "eps_xy" (magicl:tref (cl-mpm/particle:mp-strain mp) 5 0))

        (when (= nd 3)
          (cl-mpm/output::save-parameter "eps_zz" (magicl:tref (cl-mpm/particle:mp-strain mp) 2 0))
          (cl-mpm/output::save-parameter "eps_yz" (magicl:tref (cl-mpm/particle:mp-strain mp) 3 0))
          (cl-mpm/output::save-parameter "eps_zx" (magicl:tref (cl-mpm/particle:mp-strain mp) 4 0)))

        (cl-mpm/output::save-parameter "size_x" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 0 0))
        (cl-mpm/output::save-parameter "size_y" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0))
        ))))
(defparameter *eta* 1d0)
;; (defmethod cl-mpm/particle::constitutive-model ((mp cl-mpm/particle::particle-elastic-damage-delayed) strain dt)
;;   "Strain intergrated elsewhere, just using elastic tensor"
;;   (with-accessors ((E                cl-mpm/particle::mp-E)
;;                    (nu               cl-mpm/particle::mp-nu)
;;                    (de               cl-mpm/particle::mp-elastic-matrix)
;;                    (stress           cl-mpm/particle::mp-stress)
;;                    (stress-undamaged cl-mpm/particle::mp-undamaged-stress)
;;                    (L                cl-mpm/particle::mp-stretch-tensor)
;;                    (damage           cl-mpm/particle::mp-damage)
;;                )
;;       mp
;;     (declare (double-float damage))
;;     (setf stress-undamaged (cl-mpm/constitutive::linear-elastic-mat strain de stress-undamaged))
;;     (cl-mpm/utils::voigt-copy-into stress-undamaged stress)
;;     (when (> damage 0d0)
;;       (cl-mpm/fastmaths::fast-scale! stress (- 1d0 (* (- 1d0 1d-9) damage))))
;;     (let ((D (cl-mpm/utils:matrix-to-voigt
;;               (cl-mpm/fastmaths:fast-scale!
;;                (cl-mpm/fastmaths:fast-.+ L (magicl:transpose L))
;;                0.5d0))))
;;     (cl-mpm/fastmaths:fast-.+ stress
;;                               (cl-mpm/fastmaths:fast-scale!
;;                                D
;;                                (/ *eta* dt))
;;                               stress)
;;     )
;;     stress))

(defun stop ()
  (setf *run-sim* nil))

(defun setup-test-column (size block-size sim-type &optional (e-scale 1) (mp-scale 1))
  (let* ((sim (cl-mpm/setup::make-simple-sim;-sd
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               ;; :sim-type
               ;; 'cl-mpm::mpm-sim-sd
               :sim-type sim-type
               ;; 'cl-mpm/damage::mpm-sim-damage
               ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 1d3)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
         )
    (declare (double-float h density))
    (progn
      (let ()
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-mps-from-list
               (cl-mpm/setup::make-block-mps-list
                (mapcar (lambda (x) 0) size)
                block-size
                (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                density
                ;'cl-mpm/particle::particle-finite-viscoelastic-ice
                ;; 'cl-mpm/particle::particle-finite-viscoelastic
                ;; 'cl-mpm/particle::particle-elastic-damage-delayed
                'cl-mpm/particle::particle-elastic
                ;; 'cl-mpm/particle::particle-vm
<<<<<<< HEAD
                :E 1d6
=======
                :E 0.5d6
>>>>>>> 0076062c928785d0442066c9942d57fb0ae47e24
                :nu 0.24d0
                ;:viscosity 1.11d6
                ;; :viscosity 1d08
                ;; :visc-power 3d0
                ;; :rho 30d3
<<<<<<< HEAD
                ;; :rho-r 1d3
                ;; :softening 1d0
=======
>>>>>>> 0076062c928785d0442066c9942d57fb0ae47e24
                ;; :initiation-stress 1d4
                ;; :delay-time 1d0
                ;; :delay-exponent 1d0
                ;; :local-length h
                ;; :ductility 10d0
                :gravity -10.0d0
                :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
                ))))
      ;; (format t "Charictoristic time ~E~%" (/ ))
      (setf (cl-mpm:sim-allow-mp-split sim) t)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-enable-fbar sim) t)
      ;; (setf (cl-mpm::sim-mass-filter sim) 0d0)
      ;; (cl-mpm/setup::set-mass-filter sim density :proportion 1d-4)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      ;; (setf (cl-mpm/damage::sim-enable-length-localisation sim) t)
      (let ((ms 1d0))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        (setf (cl-mpm:sim-damping-factor sim)
              (* 0d-1
                 (cl-mpm/setup:estimate-critical-damping sim))))

      (setf *eta* (* 0.5d0
                     (cl-mpm/setup::estimate-stiffness-critical-damping sim 1d6 density)))

      (let ((dt-scale 1d0))
        (setf
         (cl-mpm:sim-dt sim)
         (* dt-scale h
            (sqrt (cl-mpm::sim-mass-scale sim))
            (sqrt (/ density (cl-mpm/particle::mp-p-modulus (aref (cl-mpm:sim-mps sim) 0)))))))

      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))
      ;; (setf (cl-mpm:sim-bcs sim)
      ;;       (cl-mpm/bc::make-outside-bc-var
      ;;        (cl-mpm:sim-mesh sim)
      ;;        (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
      ;;        (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
      ;;        (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil 0 nil)))
      ;;        (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil 0 nil)))
      ;;        ;; (lambda (i) nil)
      ;;        ;; (lambda (i) nil)
      ;;        (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
      ;;        (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
      ;;       ))
      (cl-mpm/setup::setup-bcs
       sim
       :bottom '(0 0 nil))
      ;; (setf
      ;;  (cl-mpm:sim-bcs sim)
      ;;  (cl-mpm/bc::make-outside-bc-varfix
      ;;   (cl-mpm:sim-mesh sim)
      ;;   '(0 nil 0)
      ;;   '(0 nil 0)
      ;;   '(nil 0 nil)
      ;;   '(nil 0 nil)
      ;;   '(nil nil 0)
      ;;   '(nil nil 0)))
      sim)))


(defparameter *sim* nil)
(defparameter *run-sim* t)
(defparameter *t* 0)
(defparameter *sim-step* 0)
(defparameter *refine* (/ 1d0 1d0))
(let ((refine (uiop:getenv "REFINE")))
  (when refine
    (setf *refine* (parse-integer (uiop:getenv "REFINE")))
    ))

(defun setup (&key (refine 1) (mps 2)
              (sim-type 'cl-mpm:mpm-sim-usf)
                )
  (let ((mps-per-dim mps))
    ;(defparameter *sim* (setup-test-column '(16 16 8) '(8 8 8) *refine* mps-per-dim))
    (defparameter *sim* (setup-test-column '(32 16) '(8 8) sim-type refine mps-per-dim))
    )
  ;; (cl-mpm/setup::initialise-stress-self-weight
  ;;  *sim*
  ;;  8d0)
  ;; (defparameter *sim* (setup-test-column '(1 1 1) '(1 1 1) 1 1))
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defun run (&key (output-dir "./output/")
              (ms 1d0))
  (uiop:ensure-all-directories-exist (list output-dir))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (cl-mpm/output:save-vtk-mesh (merge-pathnames output-dir "mesh.vtk")
                          *sim*)

  (defparameter *data-steps* (list))
  (defparameter *data-oobf* (list))
  (defparameter *data-energy* (list))
  (setf (cl-mpm:sim-damping-factor *sim*)
        (* 1d-2
           (cl-mpm/setup:estimate-critical-damping *sim*)))
  (cl-mpm:iterate-over-mps
   (cl-mpm:sim-mps *sim*)
   (lambda (mp)
     (setf (cl-mpm/particle::mp-enable-plasticity mp) nil)))
  (cl-mpm/dynamic-relaxation:converge-quasi-static
   *sim*
   :dt-scale 0.5d0
   :oobf-crit 1d-2
   :energy-crit 1d-2)
  (cl-mpm:iterate-over-mps
   (cl-mpm:sim-mps *sim*)
   (lambda (mp)
     (setf (cl-mpm/particle::mp-enable-plasticity mp) t)))
  (let* ((target-time 0.1d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 0.5d0)
         (dt-min (cl-mpm:sim-dt *sim*))
         )
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 1d-3
             (cl-mpm/setup:estimate-critical-damping *sim*)))
    (setf (cl-mpm:sim-mass-scale *sim*) ms)
    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup:estimate-elastic-dt *sim*))
    ;; (setf (cl-mpm:sim-damping-factor *sim*) 0d0)
    (cl-mpm::update-sim *sim*)
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                    (format t "CFL dt estimate: ~f~%" dt-e)
                    (format t "CFL step count estimate: ~D~%" substeps-e)
                    (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
    (time (loop for steps from 0 to 100
                while *run-sim*
                do
                   (let ((energy 0d0)
                         (work 0d0)
                         (oobf 0d0))
                     ;; (when t;(> steps 20)
                     ;;   (setf (cl-mpm::sim-damping-factor *sim*) 1d0)
                     ;;   ;; (setf (cl-mpm::sim-enable-damage *sim*) t)
                     ;;   )
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     ;; (cl-mpm/output::save-vtk-mesh-nodes (merge-pathnames output-dir (format nil "sim_nodes_p_~5,'0d.vtk" *sim-step*)) (cl-mpm::sim-mesh-p *sim*))
                     (setf dt-min (cl-mpm::calculate-min-dt *sim*))
                     (time
                      (dotimes (i substeps)
                        (cl-mpm::update-sim *sim*)
                        ;; (setf dt-min (min (cl-mpm::calculate-min-dt *sim*) dt-min))
                        ;; (setf (cl-mpm:sim-dt *sim*) (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
                        (incf oobf (cl-mpm::sim-stats-oobf *sim*))
                        (incf energy (cl-mpm::sim-stats-energy *sim*))
                        (incf work (cl-mpm::sim-stats-power *sim*))
                        (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))

                     (setf
                      energy (/ energy substeps)
                      oobf (/ oobf substeps))

                     (if (= work 0d0)
                         (setf energy 0d0)
                         (setf energy (abs (/ energy work))))

                     (format t "Min dt est ~f~%" dt-min)
                     (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       (format t "CFL dt estimate: ~f~%" dt-e)
                       (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf substeps substeps-e))
                     ;; (setf (cl-mpm:sim-dt *sim*) (* dt-min dt-scale))
                     ;; (setf substeps (floor target-time (cl-mpm:sim-dt *sim*)))
                     (format t "Substeps ~D~%" substeps)
                     ;; (let ((energy 0d0))
                     ;;   (setf energy (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
                     ;;   (push energy *data-energy*)
                     ;;   (push steps *data-steps*)
                     ;;   (format t "Energy ~F~%" energy))
                     (format t "~E ~E~%" energy oobf)

                     (incf *sim-step*)
                     (plot *sim*)
                     (vgplot:title (format nil "Time:~F - KE ~E - OOBF ~E - Work ~E"  (cl-mpm::sim-time *sim*) energy oobf work))
                     (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                        :terminal "png size 1920,1080"
                                        )

                     (swank.live:update-swank)
                     ))))
  (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" *sim-step*)) *sim*))


;; (setf lparallel:*kernel* (lparallel:make-kernel 8 :name "custom-kernel"))
;; (push (lambda ()
;;         (format t "Closing kernel~%")
;;         (lparallel:end-kernel))
;;       sb-ext:*exit-hooks*)
;; (setup)
;; (run)

;; (time
;;  (dotimes (i 1)
;;    (cl-mpm::update-stress (cl-mpm:sim-mesh *sim*)
;;                           (cl-mpm:sim-mps *sim*)
;;                           (cl-mpm:sim-dt *sim*))))

(defun save-data (name)
  (with-open-file (stream (merge-pathnames name) :direction :output :if-exists :supersede)
    (format stream "time,oobf,energy,damage~%")
    (loop for time in (reverse *data-t*)
          for oobf in (reverse *data-oobf*)
          for energy in (reverse *data-energy*)
          for damage in (reverse *data-damage*)
          do (format stream "~f,~f,~f,~f~%" time oobf energy damage))))


(defmacro time-form (it form)
  `(progn
     (declaim (optimize speed))
     (let* ((iterations ,it)
            (start (get-internal-real-time)))
       (time
        (dotimes (i ,it)
              ,form))
       (let* ((end (get-internal-real-time))
              (units internal-time-units-per-second)
              (dt (/ (- end start) (* iterations units)))
              )
         (format t "Total time: ~f ~%" (/ (- end start) units)) (format t "Time per iteration: ~f~%" (/ (- end start) (* iterations units)))
         (format t "Throughput: ~f~%" (/ 1 dt))
         (format t "Time per MP: ~E~%" (/ dt (length (cl-mpm:sim-mps *sim*))))
         dt))))
(defun test-mult ()
  (time
   (cl-mpm:iterate-over-mps
    (cl-mpm:sim-mps *sim*)
    (lambda (mp)
      (with-accessors ((strain cl-mpm/particle:mp-strain)
                       (de cl-mpm/particle::mp-elastic-matrix)
                       (stress cl-mpm/particle::mp-stress-kirchoff))
          mp
        (cl-mpm/constitutive::linear-elastic-mat strain de stress))))))
(defun profile ()
  (setup :refine 16)
  ;; (sb-profile:unprofile)
  ;; (sb-profile:profile "CL-MPM")
  ;; ;; (sb-profile:profile "CL-MPM/PARTICLE")
  ;; ;; (sb-profile:profile "CL-MPM/MESH")
  ;; ;; (sb-profile:profile "CL-MPM/SHAPE-FUNCTION")
  ;; (sb-profile:reset)
  (time-form 100
             (progn
               (format t "~D~%" i)
                  (cl-mpm::update-sim *sim*)))
  ;; (time-form 1000000
  ;;            (progn
  ;;              ;; (cl-mpm:iterate-over-neighbours (cl-mpm:sim-mesh *sim*) (aref (cl-mpm:sim-mps *sim*) 0) (lambda (mesh mp &rest args) (setf (fill-pointer (cl-mpm/particle::mp-cached-nodes mp)) 0)))
  ;;              (cl-mpm::iterate-over-neighbours-shape-gimp-simd (cl-mpm:sim-mesh *sim*) (aref (cl-mpm:sim-mps *sim*) 0) (lambda (&rest args)))
  ;;              ))
  ;; (time
  ;;  (dotimes (i 100)
  ;;    (format t "~D~%" i)
  ;;    (cl-mpm::update-sim *sim*)))
  (format t "MPS ~D~%" (length (cl-mpm:sim-mps *sim*)))
  ;; (sb-profile:report)
  )

(defun profile ()
  (setup :refine 16)
  (time-form 100
             (progn
               (format t "~D~%" i)
               (cl-mpm::update-sim *sim*)))
  ;; (time-form
  ;;  100
  ;;  (progn
  ;;    (cl-mpm::check-mps *sim*)))
  ;; (time
  ;;  (dotimes (i 100)
  ;;    (format t "~D~%" i)
  ;;    (cl-mpm::update-sim *sim*)))
  (format t "MPS ~D~%" (length (cl-mpm:sim-mps *sim*)))
  ;; (sb-profile:report)
  )

;; (lparallel:end-kernel)
;; (sb-ext::exit)
;; (uiop:quit)


;; (declaim (optimize (debug 0) (safety 0) (speed 3)))
;; (defun test ()
;;     (let ((iters 10000000))
;;       (let ((a (cl-mpm/utils:vector-zeros)))
;;         (time
;;          (lparallel:pdotimes (i iters)
;;            (magicl:.+ a (cl-mpm/utils:vector-zeros) a))))
;;       (let ((a (make-array 2 :element-type 'double-float)))
;;         (time
;;          (lparallel:pdotimes (i iters)
;;            (let ((b (make-array 2 :element-type 'double-float)))
;;              (loop for i fixnum from 0 to 1
;;                    do (incf (aref a i) (aref b i))))
;;            )))))


;; (cl-mpm::iterate-over-nodes-serial
;;  (cl-mpm:sim-mesh *sim*)
;;  (lambda (node)
;;    (loop for v across (magicl::storage (cl-mpm/mesh::node-velocity node))
;;          do
;;             (when (sb-ext::float-nan-p v)
;;               (format t "Found nan vel~%")
;;               (pprint node)
;;               ))))

(defun find-nans ()
  (cl-mpm::iterate-over-nodes
   (cl-mpm:sim-mesh *sim*)
   (lambda (mp)
     (loop for v across (magicl::storage (cl-mpm/mesh::node-velocity mp))
           do
              ;; (when (> v 0d0)
              ;;   (format t "~E~%" v))
              (when (or (sb-ext::float-nan-p v)
                        ;; (> v 1d10)
                        ;; (< v -1d10)
                        )
                (format t "Found nan vel~%")
                (pprint mp)
                )
           ))))



;;Old
;; Evaluation took:
;; 7.969 seconds of real time
;; 71.424494 seconds of total run time (53.687638 user, 17.736856 system)
;; [ Real times consist of 2.410 seconds GC time, and 5.559 seconds non-GC time. ]
;; [ Run times consist of 2.467 seconds GC time, and 68.958 seconds non-GC time. ]
;; 896.27% CPU
;; 33,467,579,605 processor cycles
;; 22,937,955,888 bytes consed

;; Total time: 7.969999 
;; Time per iteration: 0.07969999
;; Throughput: 12.547053
;; Time per MP: 4.8645017e-6
;; MPS 16384

;;New


(defun stop ()
  (setf *run-sim* nil)
  (setf cl-mpm/dynamic-relaxation::*run-convergance* nil))
(defun test-conv ()
  (defparameter *data-steps* (list))
  (defparameter *data-oobf* (list))
  (defparameter *data-energy* (list))
  (defparameter *data-name* (list))

  (vgplot:figure)
  (let ((path (merge-pathnames "./analysis_scripts/vel_algo/data/")))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* path)) do (uiop:delete-file-if-exists f))
    (dolist (algo (list :FLIP))
      (dolist (dt-scale (list 0.5d0))
        (dolist (refine (list 1 2 3))
          (dolist (fbar (list nil))
            (dolist (sim-type (list
                               'cl-mpm::mpm-sim-usf
                               'cl-mpm::mpm-sim-usl
                               ))
              (let* ((steps (list))
                     (energy (list))
                     (oobf (list))
                     ;; (dt-scale 0.5d0)
                     (substeps (round 50 dt-scale))
                     )
                (setup :refine refine :mps 3
                       :sim-type sim-type
                       )
                ;; (change-class *sim* sim-type)
                (setf (cl-mpm:sim-damping-factor *sim*)
                      (* 1d-1 (cl-mpm/setup::estimate-critical-damping *sim*)))
                (setf
                 (cl-mpm::sim-velocity-algorithm *sim*)
                 algo
                 (cl-mpm::sim-enable-fbar *sim*) fbar)
                (let ((name (format nil "~A-~A-~E-~A-~A"
                                    sim-type
                                    (cl-mpm::sim-velocity-algorithm *sim*)
                                    dt-scale
                                    refine
                                    fbar
                                    )))
                  (vgplot:title name)
                  (handler-case
                      (cl-mpm/dynamic-relaxation:converge-quasi-static
                       *sim*
                       :dt-scale dt-scale
                       :energy-crit 1d-9
                       :oobf-crit 1d-9
                       :substeps substeps
                       :conv-steps 100
                       :dt-scale dt-scale
                       :pic-update nil
                       :post-iter-step
                       (lambda (i e o)
                         ;; (plot *sim*)
                         (push i steps)
                         (push e energy)
                         (push o oobf)
                         (apply #'vgplot:semilogy
                                (reduce #'append
                                        (append
                                         (mapcar #'list
                                                 (append *data-steps* (list steps))
                                                 (append *data-oobf*  (list oobf))
                                                 (append *data-name*  (list name))
                                                 )))))
                       )
                    (cl-mpm/dynamic-relaxation::non-convergence-error (c)
                      (format t "Sim failed to converge~%"))
                    )
                  (with-open-file (stream (merge-pathnames path (format nil "~A.csv" name)) :direction :output :if-exists :supersede)
                    (format stream "step,energy,oobf~%")
                    (loop for s in steps
                          for e in energy
                          for o in oobf
                          do (format stream "~D,~F,~F~%" s e o)))
                  (cl-mpm/output:save-vtk (merge-pathnames path (format nil "sim_~A.vtk" name)) *sim*)
                  (cl-mpm/output:save-vtk-nodes (merge-pathnames path (format nil "sim_nodes_~A.vtk" name)) *sim*)

                  (push steps *data-steps*)
                  (push energy *data-energy*)
                  (push oobf *data-oobf*)
                  (push name *data-name*)
                  )))))))))

(defun plot-all ()
  (apply #'vgplot:semilogy
         (reduce #'append
                 (mapcar #'list
                         *data-steps*
                         *data-oobf*
                         *data-name*
                         ))))

(defun plot-some ()
  (let* ((indexes (list 0 1 2
                        9 10 11
                        ))
         (steps (loop for i in indexes collect (nth i *data-steps*)))
         (energy (loop for i in indexes collect (nth i *data-energy*)))
         (oobf (loop for i in indexes collect (nth i *data-oobf*)))
         (name (loop for i in indexes collect (nth i *data-name*)))
        )
    (apply #'vgplot:semilogy
           (reduce #'append
                   (mapcar #'list
                           steps
                           oobf
                           name
                           )))))

;; (let ((mesh (cl-mpm:sim-mesh *sim*))
;;       (mp0 (aref (cl-mpm:sim-mps *sim*) 0)))
;;   (pprint (cl-mpm/particle:mp-position mp0))
;;   (cl-mpm::iterate-over-neighbours
;;    mesh mp0
;;    (lambda (mesh mp node svp grads fsvp fgrads)
;;      (with-accessors ((node-active cl-mpm/mesh:node-active)
;;                       (node-j-inc cl-mpm/mesh::node-jacobian-inc)
;;                       (node-volume cl-mpm/mesh::node-volume))
;;          node
;;        (pprint svp)))))

(defun test-ms ()
  (setf *run-sim* t)
  (loop for  ms in (list 1d0 1d1 1d2 1d3)
        while *run-sim*
        do
           (progn
             (setup :refine 2)
             (run :output-dir (format nil "./output-~F/" ms) :ms ms))))

(defun test-usl-usf ()
  (setf *run-sim* t)
  (loop for sim-type in (list
                         'cl-mpm::mpm-sim-usf
                         'cl-mpm::mpm-sim-usl)
        while *run-sim*
        do
           (progn
             (setup :refine 8 :mps 2 :sim-type sim-type)
             (setf (cl-mpm:sim-enable-fbar *sim*) t)
             (run :output-dir (format nil "./output-~A/" sim-type)))))


(defun p2g (mesh x y)
  )
(defclass particle-component ()
  (
   (pos-x
    :accessor particle-component-pos-x
    ))
  )

(defun test ()
  (setup)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)
                   )
      *sim*
      (let* ((mps-count (length mps))
             (pos-x (make-array   mps-count :initial-element 0d0 :element-type 'double-float))
             (pos-y (make-array   mps-count :initial-element 0d0 :element-type 'double-float))
             (size-x (make-array  mps-count :initial-element 0d0 :element-type 'double-float))
             (size-y (make-array  mps-count :initial-element 0d0 :element-type 'double-float))
             ;; (index-x (make-array mps-count :initial-element 0 :element-type 'fixnum))
             ;; (index-y (make-array mps-count :initial-element 0 :element-type 'fixnum))
             (h (cl-mpm/mesh:mesh-resolution mesh))
             )
        (loop for x across pos-x
              for y across pos-y
              do
                 (let ((index-x (floor x h))
                       (index-y (floor y h)))

                   )
                 ;; (setf index-x (floor x h)
                 ;;       index-y (floor y h))
              ))))

(defun test-max-conv ()
  (defparameter *data-steps* (list))
  (defparameter *data-oobf* (list))
  (defparameter *data-energy* (list))
  (defparameter *data-name* (list))

  (vgplot:figure)
  (let ((path (merge-pathnames "./analysis_scripts/vel_algo/data/")))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* path)) do (uiop:delete-file-if-exists f))
    (dolist (algo (list :FLIP))
      (dolist (dt-scale (list 1d0 0.5d0 0.1d0))
        (dolist (refine (list 1 2 3))
          (dolist (fbar (list nil))
            (dolist (sim-type (list
                               'cl-mpm::mpm-sim-usf
                               'cl-mpm::mpm-sim-usl
                               ))
              (let* ((steps (list))
                     (energy (list))
                     (oobf (list))
                     ;; (dt-scale 0.5d0)
                     (substeps (round (* refine 50) dt-scale))
                     )
                (setup :refine refine :mps 3
                       :sim-type sim-type
                       )
                ;; (change-class *sim* sim-type)
                (setf (cl-mpm:sim-damping-factor *sim*)
                      (* 1d-1 (cl-mpm/setup::estimate-critical-damping *sim*)))
                (setf
                 (cl-mpm::sim-velocity-algorithm *sim*)
                 algo
                 (cl-mpm::sim-enable-fbar *sim*) fbar)
                (let ((name (format nil "~A-~A-~E-~A-~A"
                                    sim-type
                                    (cl-mpm::sim-velocity-algorithm *sim*)
                                    dt-scale
                                    refine
                                    fbar
                                    )))
                  (vgplot:title name)
                  (let ((work 0d0)
                        (e 0d0)
                        (o 0d0)
                        (sim *sim*)
                        (max-steps 200))
                      (loop for i from 0 to max-steps
                            while *run-sim*
                            do
                               (progn
                                 (dotimes (j substeps)
                                   (cl-mpm:update-sim *sim*))
                                 (setf (cl-mpm:sim-dt sim) (* dt-scale (cl-mpm::calculate-min-dt sim)))
                                 (setf e (cl-mpm::sim-stats-energy sim))
                                 (if (= work 0d0)
                                     (setf e 0d0)
                                     (setf e (abs (/ e work))))
                                 (setf o (cl-mpm::sim-stats-oobf sim))
                                 (push i steps)
                                 (push e energy)
                                 (push o oobf)
                                 (apply #'vgplot:semilogy
                                        (reduce #'append
                                                (append
                                                 (mapcar #'list
                                                         (append *data-steps* (list steps))
                                                         (append *data-oobf*  (list oobf))
                                                         (append *data-name*  (list name))
                                                         )))))
                               (swank.live:update-swank)
                       ))
                  (with-open-file (stream (merge-pathnames path (format nil "~A.csv" name)) :direction :output :if-exists :supersede)
                    (format stream "step,energy,oobf~%")
                    (loop for s in steps
                          for e in energy
                          for o in oobf
                          do (format stream "~D,~F,~F~%" s e o)))
                  (cl-mpm/output:save-vtk (merge-pathnames path (format nil "sim_~A.vtk" name)) *sim*)
                  (cl-mpm/output:save-vtk-nodes (merge-pathnames path (format nil "sim_nodes_~A.vtk" name)) *sim*)

                  (push steps *data-steps*)
                  (push energy *data-energy*)
                  (push oobf *data-oobf*)
                  (push name *data-name*)
                  )))))))))
