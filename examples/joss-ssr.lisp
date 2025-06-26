(defpackage :cl-mpm/examples/joss
  (:use :cl))

(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)

(setf *block-compile-default* t)

(in-package :cl-mpm/examples/joss)
;(declaim (optimize (debug 3) (safety 3) (speed 0)))
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(defun full-recompile ()
  ;; (asdf:compile-system :cl-mpm/utils :force t)
  ;; (asdf:compile-system :cl-mpm/fastmaths :force t)
  ;; (asdf:compile-system :cl-mpm/forces :force t)
  ;; (asdf:compile-system :cl-mpm/constitutive :force t)
  ;; (asdf:compile-system :cl-mpm/mesh :force t)
  ;; (asdf:compile-system :cl-mpm/particle :force t)
  ;; (asdf:compile-system :cl-mpm/bc :force t)
  ;; (asdf:compile-system :cl-mpm :force t)
  (asdf:compile-system :cl-mpm/all :force t)
  (ql:quickload :cl-mpm/examples/joss)
  )

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-elastic) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  ;; (cl-mpm::update-domain-det mesh mp dt)
  ;; (cl-mpm::update-domain-corner mesh mp dt)
  ;; (cl-mpm::co-domain-corner-2d mesh mp dt)
  (cl-mpm::update-domain-polar-2d mesh mp dt)
  (cl-mpm::scale-domain-size mesh mp)
  )

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-elastic) dt fbar)
  (cl-mpm::update-stress-kirchoff mesh mp dt fbar)
  )



(declaim (notinline plot))
(defun plot (sim)
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :deformed
   ;; :colour-func #'cl-mpm/particle::mp-damage
   :colour-func #'cl-mpm/particle::mp-strain-plastic-vm
   )
  )


(declaim (notinline setup-test-column))
(defun setup-test-column (size block-size offset &optional (e-scale 1) (mp-scale 1))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               :sim-type
               'cl-mpm::mpm-sim-usf
               ;; 'cl-mpm/damage::mpm-sim-usl-damage
               ;; 'cl-mpm/damage::mpm-sim-damage
               ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 1.7d3)
         )
    (declare (double-float h density))
    (progn
      (let* (
             ;; (init-stress (* 26d3 0.63d0))
             (E 1d9)
             (angle 50d0)
             (init-c 26d3)
             )
        (cl-mpm:add-mps
         sim
         (cl-mpm/setup::make-mps-from-list
          (cl-mpm/setup::make-block-mps-list
           offset
           block-size
           (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
           density
           'cl-mpm/particle::particle-mc
           :E 1d9
           :nu 0.24d0
           :psi (* 0d0 (/ pi 180))
           :phi (* angle (/ pi 180))
           :c init-c
           :psi-r (* 0d0 (/ pi 180))
           :phi-r (* 30d0 (/ pi 180))
           :c-r 0d0
           :softening 0d0

           :gravity (* -9.8d0 1d0)
           :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
           ))))

      (setf (cl-mpm:sim-allow-mp-split sim) t)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-velocity-algorithm sim) :BLEND)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-enable-fbar sim) t)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (cl-mpm/setup::set-mass-filter sim density :proportion 1d-3)

      (let ((ms 1d0))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        (setf (cl-mpm:sim-damping-factor sim)
              (* 0.1d0
                 (sqrt ms)
                 (cl-mpm/setup::estimate-critical-damping sim)))
        )
      (format t "Set damping: ~F~%" (cl-mpm:sim-damping-factor sim))

      (setf
       (cl-mpm:sim-dt sim)
       (cl-mpm/setup:estimate-elastic-dt sim :dt-scale 0.1d0))

      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))
      (setf
       (cl-mpm:sim-bcs sim)
       (cl-mpm/bc::make-outside-bc-varfix
        (cl-mpm:sim-mesh sim)
        '(0 nil 0)
        '(0 nil 0)
        '(0 0 0)
        '(0 0 0)
        '(nil nil 0)
        '(nil nil 0)))


      sim)))


(defun estimate-gf (eta ft lc &optional (E 1d9))
  (let* (
         (gf (* eta (/ (expt ft 2) (* 2 E))))
         )
    (* gf lc)))




(defparameter *sim* nil)
(defparameter *run-sim* t)
(defparameter *t* 0)
(defparameter *sim-step* 0)
(defparameter *refine* (/ 1d0 2d0))

(defun plot-collapse-curves ()
  ;; (vgplot:figure)
  (when (> (length *data-t*) 0)
    (vgplot:plot
     *data-t* (mapcar (lambda (x) (/ x (reduce #'max *data-e*))) *data-e*) "KE"
     *data-t* (mapcar (lambda (x) (/ x (reduce #'max *data-oobf*))) *data-oobf*) "oobf"
     *data-t* (mapcar (lambda (x) (/ x (reduce #'max *data-d*))) *data-d*) "Damage"
     ))
  )

(defun estimate-total-energy (sim)
  (loop for mp across (cl-mpm:sim-mps sim)
        sum (* (cl-mpm/particle::mp-mass mp)
                      (cl-mpm/fastmaths::mag-squared (cl-mpm/particle::mp-velocity mp)))))

(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk") *sim*)
  ;; (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
  (let* ((target-time 1d0)
         (target-time-original target-time)
         (mass-scale (cl-mpm::sim-mass-scale *sim*))
         (accelerate-target-time 1d1)
         (accelerate-mass-scale 1d4)
         (collapse-target-time 0.1d0)
         (collapse-mass-scale 1d0)
         (plasticity-enabled (cl-mpm/particle::mp-enable-plasticity (aref (cl-mpm:sim-mps *sim*) 0)))
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 0.8d0)
         (settle-steps 0)
         (damp-steps 0)
         (sim-state :settle)
         (dt-0 0d0)
         (last-e 0d0)
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
         (enable-damage t)
         (criteria-energy 1d-1)
         (criteria-oobf 1d-1)
         (damping-0
           (* 1d-4
              (cl-mpm/setup::estimate-critical-damping *sim*)))
         (damage-0
           (lparallel:pmap-reduce (lambda (mp)
                                    (*
                                     (cl-mpm/particle::mp-damage mp)))
                                  #'+ (cl-mpm:sim-mps *sim*)
                                  :initial-value 0d0))
         )

    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (cl-mpm/output::save-simulation-parameters
     #p"output/settings.json"
     *sim*
     (list :dt target-time))
    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp) (setf (cl-mpm/particle::mp-enable-plasticity mp) nil)))

    (with-open-file (stream (merge-pathnames "output/timesteps.csv") :direction :output :if-exists :supersede)
      (format stream "steps,time,damage,plastic,energy,oobf,step-type,mass~%"))

    ;; (cl-mpm::update-sim *sim*)
    ;; (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
    ;;   (format t "CFL dt estimate: ~f~%" dt-e)
    ;;   (format t "CFL step count estimate: ~D~%" substeps-e)
    ;;   (setf substeps substeps-e))

    (setf dt-0 (/ (cl-mpm:sim-dt *sim*) (sqrt (cl-mpm::sim-mass-scale *sim*))))
    (defparameter *data-damage* 0d0)
    (defparameter *data-plastic* 0d0)
    (defparameter *data-energy* 0d0)
    (defparameter *inst-data-oobf* 0d0)
    (defparameter *data-mass* (lparallel:pmap-reduce
                               #'cl-mpm/particle::mp-mass
                               #'+
                               (cl-mpm:sim-mps *sim*)
                               :initial-value 0d0))
    (setf *data-d* (list))


    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (let ((ms 1d0))
      (setf (cl-mpm:sim-mass-scale *sim*) ms)
      (setf (cl-mpm:sim-damping-factor *sim*)
            (* 0.1d0
               (cl-mpm/setup::estimate-critical-damping *sim*))))

    (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_conv_~5,'0d.vtk" 0)) *sim*)

    (defparameter *collapse* nil)
    (cl-mpm/dynamic-relaxation:converge-quasi-static
     *sim*
     :dt-scale dt-scale
     :energy-crit 1d-1
     :oobf-crit 1d-1
     :substeps 50
     :conv-steps 1000
     :dt-scale dt-scale
     :post-iter-step
     (lambda (i e o)
       ;; (cl-mpm/damage::calculate-damage *sim*)
       (plot *sim*)
       (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_conv_~5,'0d.vtk" (1+ i))) *sim*)
       ;; (let ((dy (lparallel:pmap-reduce (lambda (mp)
       ;;                                    (cl-mpm/particle::mp-damage-ybar mp))
       ;;                                  #'max (cl-mpm:sim-mps *sim*))))
       ;;   (format t "Max damage-ybar ~E~%" dy)
       ;;   )
       ;; (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" i))
       ;;                    :terminal "png size 1920,1080")
       ))
    (cl-mpm:iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (cl-mpm/fastmaths::fast-zero (cl-mpm/particle:mp-velocity mp))
       (cl-mpm/fastmaths::fast-zero (cl-mpm/particle::mp-acceleration mp))))

    (when t
      (setf sim-state :settle)
      (setf (cl-mpm:sim-damping-factor *sim*)
            damping-0
            ))
    (when t
      (setf (cl-mpm::sim-enable-damage *sim*) enable-damage)
      (cl-mpm::iterate-over-mps
       (cl-mpm:sim-mps *sim*)
       (lambda (mp) (setf (cl-mpm/particle::mp-enable-plasticity mp) plasticity-enabled))))

    (setf (cl-mpm::sim-mass-scale *sim*) accelerate-mass-scale)
    (setf target-time accelerate-target-time)
    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (setf substeps (floor target-time (cl-mpm:sim-dt *sim*)))
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 0d-3
             (cl-mpm/setup::estimate-critical-damping *sim*)))

    (format t "Substeps ~D~%" substeps)
    (loop for vals in *to-damage-mps*
          do
          (destructuring-bind (mp k) vals
            (setf (cl-mpm/particle::mp-history-stress mp) k)
            (cl-mpm/damage::update-damage mp 1d-3)))

    (let ((work 0d0)
          )
      (time (loop for steps from 0 to 500
                  while *run-sim*
                  do
                     (progn
                       (format t "Step ~d ~%" steps)
                       (format t "MPs ~d ~%" (length (cl-mpm:sim-mps *sim*)))
                       (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                       (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                       (with-open-file (stream (merge-pathnames "output/timesteps.csv") :direction :output :if-exists :append)
                         (format stream "~D,~f,~f,~f,~f,~f,~A,~f~%"
                                 steps
                                 *t*
                                 *data-damage*
                                 *data-plastic*
                                 *data-energy*
                                 *inst-data-oobf*
                                 sim-state
                                 *data-mass*))
                       (let ((energy-estimate 0d0)
                             (oobf 0d0)
                             (total-energy 0d0)
                             (total-strain-energy 0d0)
                             (total-gpe 0d0)
                             (damage-inc 0d0)
                             (plastic-inc 0d0)
                             ;; (work 0d0)
                             
                             )
                         (time
                          (let ((current-damage (cl-mpm::sim-enable-damage *sim*)))
                            (loop for i from 1 to substeps
                                  while *run-sim*
                                  do (progn
                                       (cl-mpm::update-sim *sim*)
                                       (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))
                                       (incf work (cl-mpm/dynamic-relaxation::estimate-power-norm *sim*))
                                       (incf oobf (cl-mpm/dynamic-relaxation::estimate-oobf *sim*))
                                       (incf energy-estimate (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
                                       ))))


                         ;; (setf work (cl-mpm/dynamic-relaxation::estimate-power-norm *sim*))
                         ;; (setf energy-estimate (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
                         ;; (setf oobf (cl-mpm/dynamic-relaxation::estimate-oobf *sim*))
                         (setf
                          energy-estimate (/ energy-estimate substeps)
                          oobf (/ oobf substeps))
                         (if (= work 0d0)
                             (setf energy-estimate 0d0)
                             (setf energy-estimate (abs (/ energy-estimate work))))

                         ;; (setf
                         ;;  energy-estimate (/ energy-estimate substeps)
                         ;;  oobf (/ oobf substeps))
                         ;; (setf energy-estimate
                         ;;       (/
                         ;;        energy-estimate
                         ;;        (lparallel:pmap-reduce
                         ;;         #'cl-mpm/particle::mp-mass
                         ;;         #'+
                         ;;         (cl-mpm:sim-mps *sim*)
                         ;;         :initial-value 0d0)))


                         ;; (setf energy-estimate (/
                         ;;                        (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*)
                         ;;                        (lparallel:pmap-reduce
                         ;;                         #'cl-mpm/particle::mp-mass
                         ;;                         #'+
                         ;;                         (cl-mpm:sim-mps *sim*)
                         ;;                         :initial-value 0d0)))

                         ;; (setf total-energy (abs (/ (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*) work)))
                         (setf *data-energy* energy-estimate)
                         (setf *inst-data-oobf* oobf)
                         (let ((damage-est
                                 (-
                                  (lparallel:pmap-reduce (lambda (mp)
                                                           (*
                                                            1d0
                                                            (cl-mpm/particle::mp-damage mp)))
                                                         #'+ (cl-mpm:sim-mps *sim*)
                                                         :initial-value 0d0)
                                  damage-0)))
                           (setf *data-damage* damage-est))
                         (let ((damage-est
                                 (lparallel:pmap-reduce (lambda (mp)
                                                          (*
                                                           1d0
                                                           (cl-mpm/particle::mp-strain-plastic-vm mp)))
                                                        #'+ (cl-mpm:sim-mps *sim*)
                                                        :initial-value 0d0)
                                 ))
                           (setf *data-plastic* damage-est))
                         ;; (setf
                         ;;  energy-estimate
                         ;;  (/
                         ;;   (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*)
                         ;;   (lparallel:pmap-reduce
                         ;;    #'cl-mpm/particle::mp-mass
                         ;;    #'+
                         ;;    (cl-mpm:sim-mps *sim*)
                         ;;    :initial-value 0d0)))
                         ;; (setf work (/ work (* target-time (cl-mpm::sim-mass-scale *sim*))))
                         (setf *data-mass* (lparallel:pmap-reduce
                                            #'cl-mpm/particle::mp-mass
                                            #'+
                                            (cl-mpm:sim-mps *sim*)
                                            :initial-value 0d0))

                         (push *t* *data-t*)
                         (push total-energy *data-e*)
                         (push oobf *data-oobf*)

                         ;; (push *data-damage* *data-d*)

                         ;; (push total-strain-energy *data-se*)
                         ;; (push total-gpe *data-gpe*)
                         ;; (push damage-inc *data-di*)
                         ;; (push plastic-inc *data-pi*)
                         ;; (format t "Energy estimate: ~E~%" energy-estimate)
                         (format t "Energy: ~E~%" energy-estimate)
                         (format t "OOBF : ~E~%" oobf)
                         (format t "Work : ~E~%" work)

                         (setf last-e total-energy)
                         (defparameter *oobf* oobf)
                         (defparameter *energy* energy-estimate)
                         (defparameter *total-energy* total-energy)
                         (defparameter *total-strain-energy* total-strain-energy)
                         (defparameter *damage-inc* damage-inc)
                         (defparameter *plastic-inc* plastic-inc)

                         (when (= steps damp-steps)
                           (setf sim-state :settle)
                           (setf (cl-mpm:sim-damping-factor *sim*) damping-0))
                         (when (= steps settle-steps)
                           (setf (cl-mpm::sim-enable-damage *sim*) enable-damage)
                           ;; (setf (cl-mpm/buoyancy::bc-enable *bc-erode*) t)
                           (cl-mpm::iterate-over-mps
                            (cl-mpm:sim-mps *sim*)
                            (lambda (mp) (setf (cl-mpm/particle::mp-enable-plasticity mp) plasticity-enabled))))
                         (when (>= steps settle-steps)
                           (if (or
                                ;; t
                                (> energy-estimate criteria-energy)
                                (> oobf criteria-oobf)
                                ;; t
                                ;; nil
                                ;; (> work 1d6)
                                )
                               (when (not (eq sim-state :collapse))
                                 (setf sim-state :collapse)
                                 (format t "Changed to collapse~%")
                                 ;; (setf work 0d0)
                                 )
                               (progn
                                 (when (not (eq sim-state :accelerate))
                                   (format t "Changed to accelerate~%")
                                   (setf work 0d0)
                                   (setf sim-state :accelerate)
                                   ;; (cl-mpm::remove-mps-func
                                   ;;  *sim*
                                   ;;  (lambda (p)
                                   ;;    (and (> (cl-mpm::mp-damage p) 0.99d0)
                                   ;;         (= (cl-mpm/particle::mp-index p) 0))
                                   ;;    ))
                                   (cl-mpm:iterate-over-mps
                                    (cl-mpm:sim-mps *sim*)
                                    (lambda (mp)
                                      ;; (cl-mpm/fastmaths::fast-zero (cl-mpm/particle:mp-acceleration mp))
                                      (cl-mpm/fastmaths::fast-zero (cl-mpm/particle:mp-velocity mp))
                                      ))
                                   )))
                           (case sim-state
                             (:accelerate
                              (format t "Accelerate timestep~%")
                              (setf
                               target-time accelerate-target-time
                               (cl-mpm::sim-mass-scale *sim*) accelerate-mass-scale))
                             (:collapse
                              (format t "Collapse timestep~%")
                              (setf
                               work 0d0
                               target-time collapse-target-time
                               (cl-mpm::sim-mass-scale *sim*) collapse-mass-scale))))

                         ;; (when (> steps 20)
                         ;;   (setf target-time 1d2)
                         ;;   ;; (setf *collapse* t)
                         ;;   )
                         (format t "Sim state - ~A~%" sim-state)


                         (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                           (format t "CFL dt estimate: ~f~%" dt-e)
                           (format t "CFL step count estimate: ~D~%" substeps-e)
                           (setf substeps substeps-e))

                         ;; (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
                         ;; (setf substeps (floor target-time (cl-mpm:sim-dt *sim*)))

                         ;; (setf (cl-mpm:sim-dt *sim*) (* min-dt dt-scale))
                         ;; (setf substeps (floor target-time (cl-mpm:sim-dt *sim*)))
                         ;; (format t "CFL dt estimate: ~f~%" dt)
                         ;; (format t "CFL step count estimate: ~D~%" substeps)
                         ;; (setf (cl-mpm:sim-damping-factor *sim*)
                         ;;       (* (cl-mpm:sim-damping-factor *sim*) (expt 1d-3 1/40)))

                         ;; (incf *sim-step*)
                         (plot *sim*)
                         (vgplot:title (format nil "Time:~F - KE ~E - OOBF ~E - Work ~E - ~A"  time energy oobf work sim-state)))
                       (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                          :terminal "png size 1920,1080"
                                          )

                       (swank.live:update-swank)
                       )))))
  (values))


;; (setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))
;; (push (lambda ()
;;         (format t "Closing kernel~%")
;;         (lparallel:end-kernel))
;;       sb-ext:*exit-hooks*)
;; (setup)
;; (run)

(defun test-update-stress ()
  (setup-3d)
  (let* ((mesh (cl-mpm:sim-mesh *sim*))
         (mps (cl-mpm:sim-mps *sim*))
         (mp (aref mps 0))
         (dt (cl-mpm:sim-dt *sim*))
         (iters 10))
    (cl-mpm::update-sim *sim*)
    ;; (cl-mpm::update-stress-mp mesh mp dt nil)
    ;; (time-form
    ;;  iters
    ;;  (progn
    ;;    (setf
    ;;     (cl-mpm/particle:mp-damage mp)
    ;;     0.5d0
    ;;     (cl-mpm/particle::mp-strain mp)
    ;;     (cl-mpm/utils:voigt-from-list (list 0d0 0d0 0d0 0d0 0d0 0d0)))
    ;;    (cl-mpm/particle:constitutive-model mp (cl-mpm/particle::mp-strain mp) dt)
    ;;    ))
    ;; (time-form
    ;;  iters
    ;;  (progn
    ;;    (setf
    ;;     (cl-mpm/particle:mp-damage mp)
    ;;     0.5d0
    ;;     (cl-mpm/particle:mp-strain mp) (cl-mpm/utils:voigt-from-list (list 10d0 0d0 0d0 0d0 0d0 0d0)))
    ;;    (cl-mpm/particle:constitutive-model mp (cl-mpm/particle::mp-strain mp) dt)
    ;;   ))
    ;; (time-form
    ;;  iters
    ;;  (progn
    ;;    (setf (cl-mpm/particle:mp-strain mp) (cl-mpm/utils:voigt-from-list (list 0d0 0d0 0d0 0d0 0d0 0d0)))
    ;;    (cl-mpm/particle:constitutive-model mp (cl-mpm/particle::mp-strain mp) dt)))

    (format t "Elastic:~%")
    (time-form
     iters
     (progn
       (loop for mp across mps
             do (setf (cl-mpm/particle:mp-strain mp) (cl-mpm/utils:voigt-from-list (list 0d0 0d0 0d0 0d0 0d0 0d0))))
       (cl-mpm::update-stress mesh mps dt nil)
      ))

    (format t "Plastic:~%")
    (time-form
     iters
     (progn
       (loop for mp across mps
             do (setf (cl-mpm/particle:mp-strain mp) (cl-mpm/utils:voigt-from-list (list 1d0 0d0 0d0 0d0 0d0 0d0))))
       (cl-mpm::update-stress mesh mps dt nil)
       ))

    (format t "Disable plastic:~%")
    (loop for mp across mps
          do (setf (cl-mpm/particle::mp-enable-plasticity mp) nil))
    (time-form
     iters
     (progn
       (loop for mp across mps
             do (setf (cl-mpm/particle:mp-strain mp) (cl-mpm/utils:voigt-from-list (list 0d0 0d0 0d0 0d0 0d0 0d0))))
       (cl-mpm::update-stress mesh mps dt nil)
       ))
    ;; (time
    ;;  (dotimes (i 10000)
    ;;    ;; (cl-mpm::calculate-strain-rate mesh mp dt)
    ;;    ;; (cl-mpm/particle:constitutive-model mp (cl-mpm::mp-strain mp) dt)
    ;;    ;; (cl-mpm::update-strain-kirchoff mesh mp dt nil)
    ;;    ;; (cl-mpm::update-domain-midpoint mesh mp dt)
    ;;    ;; (cl-mpm::update-stress-mp mesh mp dt nil)
    ;;    ;; (cl-mpm::post-stress-step mesh mp dt)
    ;;    ))
    ))

(defun test ()
  (setup-3d)
  ;; (sb-profile:profile "CL-MPM")
  ;; (sb-profile:profile "CL-MPM/PARTICLE")
  ;; (sb-profile:profile "CL-MPM/CONSTITUTIVE")
  ;; ;; (sb-profile:profile "CL-MPM/MESH")
  ;; ;; (sb-profile:profile "CL-MPM/SHAPE-FUNCTION")
  ;; (sb-profile:reset)
  (setf (cl-mpm::sim-enable-damage *sim*) nil)
  (time
   (dotimes (i 10)
         (cl-mpm::update-sim *sim*)))
  ;; (sb-profile:report)
  )
;; (defun test-undercut ()
;;   (cl-mpm/output:save-vtk-mesh (merge-pathnames "output_chalk/mesh.vtk")
;;                                *sim*)
;;   (loop for c in (list 0d0 10d0 20d0 30d0 40d0 50d0)
;;         while *run-sim*
;;         do
;;            (progn
;;              (setup :undercut (- c))
;;              (run)
;;              (cl-mpm/output:save-vtk (merge-pathnames (format nil "output_chalk/chalk_~5,'0d.vtk" c)) *sim*)
;;              )))

(defun profile ()
  (sb-profile:unprofile)
  (sb-profile:reset)
  ;; (sb-profile:profile "CL-MPM")
  ;; (sb-profile:profile "CL-MPM/PARTICLE")
  ;; (sb-profile:profile "CL-MPM/MESH")
  (loop repeat 100
        do (progn
             (cl-mpm::update-sim *sim*)
             ))
  (sb-profile:report))

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


(defmacro time-form (it form)
  `(progn
     (declaim (optimize speed))
     (let* ((iterations ,it)
            (start (get-internal-real-time)))
       (dotimes (i ,it)
         ,form)
       (let* ((end (get-internal-real-time))
              (units internal-time-units-per-second)
              (dt (/ (- end start) (* iterations units)))
              )
         (format t "Total time: ~f ~%" (/ (- end start) units)) (format t "Time per iteration: ~f~%" (/ (- end start) (* iterations units)))
         (format t "Throughput: ~f~%" (/ 1 dt))
         dt))))

;; (defun test ()
;;   (let ((data-cores '())
;;         (data-dt '()))
;;     (loop for k from 1 to 8
;;           do
;;              (progn
;;                (setf lparallel:*kernel* (lparallel:make-kernel k :name "custom-kernel"))
;;                (format t "Kernel: ~D~%" k)
;;                (let* ((iterations 100000)
;;                       (start (get-internal-real-time)))
;;                  (let ((a (cl-mpm/utils:vector-zeros)))
;;                    (time
;;                     (lparallel:pdotimes (i iterations)
;;                       (cl-mpm::update-strain-kirchoff (cl-mpm:sim-mesh *sim*) (aref (cl-mpm:sim-mps *sim*) 0) 0d0 nil)
;;                       )))
;;                  (let* ((end (get-internal-real-time))
;;                         (units internal-time-units-per-second)
;;                         (dt (/ (- end start) (* iterations units)))
;;                         )
;;                    (format t "Total time: ~f ~%" (/ (- end start) units)) (format t "Time per iteration: ~f~%" (/ (- end start) (* iterations units)))
;;                    (format t "Throughput: ~f~%" (/ 1 dt))
;;                    (push (/ 1d0 dt) data-dt)
;;                    (push k data-cores)
;;                    ))

;;                ;; (let ((iters 100000))
;;                ;;   (let ((a (cl-mpm/utils:vector-zeros)))
;;                ;;     (time
;;                ;;      (lparallel:pdotimes (i iters)
;;                ;;        (magicl:@ (magicl:eye 3) (cl-mpm/utils:vector-zeros))))))
;;                (lparallel:end-kernel)
;;                ))
;;     (vgplot:close-all-plots)
;;     (vgplot:plot data-cores data-dt))
;;   )


(defun test-mc (exx eyy ezz eyz ezx exy)
  (let* ((eps (cl-mpm/constitutive::swizzle-coombs->voigt
               (cl-mpm/utils:voigt-from-list (list exx eyy ezz eyz ezx exy))))
         (E 1d0)

         (nu 0.0440d0)
         (angle 1.22d0)
         (c 0.1101d0)
         (de (cl-mpm/constitutive::linear-elastic-matrix E nu))
         (sig (magicl:@ de eps)))
    (multiple-value-bind (sig eps f) (cl-mpm/constitutive::mc-plastic sig de eps E nu angle angle c)
      (pprint (cl-mpm/constitutive::swizzle-voigt->coombs sig))
      (pprint (cl-mpm/constitutive::swizzle-voigt->coombs eps))
      ;; (pprint (cl-mpm/utils:vector-from-list (multiple-value-list (cl-mpm/damage::principal-stresses-3d sig))))
      ;; (format t "Recalculated f ~E~%" (cl-mpm/constitutive::mc-yield-func (cl-mpm/utils:vector-from-list (multiple-value-list (cl-mpm/damage::principal-stresses-3d sig))) angle c))
      )
    ))

(defun test-mc-random ()
  (with-open-file (stream (uiop:merge-pathnames* "random.csv")
                          :direction :output
                          :if-exists :supersede)
    (dotimes (i 1000)
      (let* ((eps (cl-mpm/utils:voigt-from-list (loop for i from 0 to 5
                                                      collect (- (random 2d0) 1d0))))
             (E 1d0)
             (nu (random 0.45d0))
             (angle (* (random 80d0) (/ pi 180d0)))
             (c (random 1.0d0))
             (de (cl-mpm/constitutive::linear-elastic-matrix E nu))
             (sig-i (magicl:@ de eps)))
        (format t "Iter ~D~%" i)
        (pprint eps)
        (multiple-value-bind (sig eps-e f) (cl-mpm/constitutive::mc-plastic sig-i de eps E nu angle angle c)
          ;; (let ((actual-f (cl-mpm/constitutive::mc-yield-func (cl-mpm/utils:vector-from-list (multiple-value-list (cl-mpm/damage::principal-stresses-3d sig))) angle c)))
          ;;   (when (> actual-f 1d-6)
          ;;     (format t "MC f: ~E - actual f: ~E ~%" f actual-f)
          ;;     ;; (pprint )
          ;;     ;; (break)
          ;;     (cl-mpm/constitutive::mc-plastic sig de eps E nu angle angle c)
          ;;     ))
          (cl-mpm/constitutive::swizzle-voigt->coombs sig)
          (loop for comp across (magicl::storage (cl-mpm/constitutive::swizzle-voigt->coombs eps))
                do (format stream "~F," comp))
          (loop for comp across (magicl::storage (cl-mpm/constitutive::swizzle-voigt->coombs eps-e))
                do (format stream "~F," comp))
          (format stream "~F,~F,~F,~F" E nu angle c)
          (format stream "~%")
          )
        ;; (pprint (cl-mpm/constitutive::mc-plastic sig de eps E nu angle angle c))
        ))))




(defun plot-stress-damage ()
  (vgplot:close-all-plots)
  (let* ((mp (aref (cl-mpm:sim-mps *sim*))))
    (with-accessors ((damage cl-mpm/particle:mp-damage)
                     (E cl-mpm/particle::mp-e)
                     (Gf cl-mpm/particle::mp-Gf)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (length cl-mpm/particle::mp-local-length)
                     ) mp
      (let* ((stress (loop for x from 1d5 to 1d6 by 1d4
                           collect x))
             (damage (mapcar (lambda (stress)
                               ;(cl-mpm/damage::brittle-concrete-linear-d stress E Gf length init-stress)
                               (cl-mpm/damage::test-damage stress 2d5 1d5)
                               ;; (cl-mpm/damage::damage-response-exponential stress E Gf length init-stress)
                               )
                             stress)))
        (vgplot:plot stress damage)))))



(defun setup (&key (notch-length 1d0) (refine 1.0d0) (mps 2)
                (mp-refine 0)
                )
  (let* ((mesh-size (/ 1d0 refine))
         (mps-per-cell mps)
         (shelf-height 15)
         ;(soil-boundary (floor (* 15 1)))
         (soil-boundary 2)
         (shelf-aspect 1)
         (runout-aspect 1)
         (shelf-length (* shelf-height shelf-aspect))
         (domain-length (+ shelf-length (* runout-aspect shelf-height)))
         (shelf-height-true shelf-height)
         (shelf-height (+ shelf-height soil-boundary))
         (depth 400)
         (offset (list (* 0 mesh-size)
                       (* 0 mesh-size)))
         )
    (format t "Mesh size ~E~%" mesh-size)
    (defparameter *sim*
      (setup-test-column (list domain-length
                               (+ shelf-height (* 5 mesh-size))
                               ;; depth
                               )
                         (list domain-length shelf-height
                               ;; depth
                               )
                         offset
                         (/ 1d0 mesh-size) mps-per-cell))
    (cl-mpm:iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (when (<= (cl-mpm/utils:varef (cl-mpm/particle:mp-position mp) 1) soil-boundary)
         (setf (cl-mpm/particle::mp-index mp) 1))))


    (cl-mpm/setup::initialise-stress-self-weight
     *sim*
     shelf-height
     )
    (when nil
      (let ((transition-point (- shelf-length (* shelf-height-true 1d0))))
        (cl-mpm/setup::initialise-stress-self-weight
         *sim*
         soil-boundary
         :scaler
         (lambda (pos)
           (let ((start transition-point)
                 (end shelf-length))
             (max 0d0
                  (min 1d0
                       (/ (- (cl-mpm/utils:varef pos 0) start)
                          (- end start)))))))
        (cl-mpm/setup::initialise-stress-self-weight
         *sim*
         (* 0.8d0 shelf-height)
         :scaler
         (lambda (pos)
           (let ((end shelf-length)
                 (start transition-point))
             (- 1d0
                (max 0d0
                     (min 1d0
                          (/ (- (cl-mpm/utils:varef pos 0) start)
                             (- end start))))))))))

    ;;Refine around tip
    (dotimes (i mp-refine)
      (dolist (dir (list :x
                         :y
                         ))
        (cl-mpm::split-mps-criteria
         *sim*
         (lambda (mp h)
           (when
               (and
                (> (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                   (- shelf-length
                      15d0
                      ))
                (> (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                   (- shelf-length 8d0))
                (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                   (+ shelf-length 2d0))
                (> (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                   (- soil-boundary 0d0)
                   ))
             dir)))))
    (loop for mp across (cl-mpm::sim-mps *sim*)
          do (setf (cl-mpm/particle::mp-split-depth mp) 0))

    (let* ((sloped-height (- (- shelf-height soil-boundary) 6.8d0))
           (measured-angle 78d0)
           (undercut-angle              ;(- 82.5d0 90d0)
             (- measured-angle 90d0))
           (normal (magicl:from-list (list
                                      (cos (- (* pi (/ undercut-angle 180d0))))
                                      (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1)))
           (sloped-inflex-point
             (magicl:from-list (list (- shelf-length (* (tan (- (* pi (/ undercut-angle 180d0)))) sloped-height))
                                     soil-boundary)
                               '(2 1) :type 'double-float))
           (edge-refine 1))
      (flet ((cutout (p)
               (if (and
                    (> (magicl:tref p 1 0) soil-boundary))
                   (if (< (magicl:tref p 1 0) (+ soil-boundary sloped-height))
                       (cl-mpm/setup::plane-point-sdf
                        (magicl:from-list (list (magicl:tref p 0 0)
                                                (magicl:tref p 1 0)) '(2 1))
                        normal
                        (magicl:from-list (list shelf-length soil-boundary)
                                          '(2 1) :type 'double-float))

                       (cl-mpm/setup::plane-point-sdf
                        (magicl:from-list (list (magicl:tref p 0 0)
                                                (magicl:tref p 1 0)) '(2 1))
                        (magicl:from-list (list 1d0 0d0) '(2 1)  :type 'double-float)
                        sloped-inflex-point)
                       )
                   1d0)))
        (dotimes (i edge-refine)
          (dolist (dir (list :x :y))
            (cl-mpm::split-mps-criteria
             *sim*
             (lambda (mp h)
               (when (and
                      (or
                       (<= (cutout
                            (cl-mpm/fastmaths:fast-.+
                             (cl-mpm/particle:mp-position mp)
                             (cl-mpm/fastmaths:fast-.*
                              (cl-mpm/particle::mp-domain-size mp)
                              (cl-mpm/utils:vector-from-list (list 0.5d0 0.5d0 0d0))))
                            ) (* mesh-size 0d0))
                       (<= (cutout
                            (cl-mpm/fastmaths:fast-.+
                             (cl-mpm/particle:mp-position mp)
                             (cl-mpm/fastmaths:fast-.*
                              (cl-mpm/particle::mp-domain-size mp)
                              (cl-mpm/utils:vector-from-list (list 0.5d0 -0.5d0 0d0))))
                            ) (* mesh-size 0d0))
                       )
                          (> (cutout (cl-mpm/particle:mp-position mp)) (- (* mesh-size 1d0))))
                 dir)))))
        (cl-mpm/setup::remove-sdf *sim*
                                  #'cutout
                                  )
        )
      (let* ((notched-depth notch-length)
           ;; (undercut-angle 45d0)
           (undercut-angle 45d0)
           (normal (magicl:from-list (list
                                      (cos (- (* pi (/ undercut-angle 180d0))))
                                      (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1)))
           (sloped-inflex-point
             (magicl:from-list (list
                                (- shelf-length notched-depth)
                                soil-boundary)
                               '(2 1) :type 'double-float)))

      (flet ((cutout (p)
               (if (and
                    (> (magicl:tref p 1 0) soil-boundary))
                   (cl-mpm/setup::plane-point-sdf
                    (magicl:from-list (list (magicl:tref p 0 0)
                                            (magicl:tref p 1 0)) '(2 1))
                    normal
                    sloped-inflex-point)
                   1d0)))
        (when (> notched-depth 0d0)
          (dotimes (i edge-refine)
            (dolist (dir (list :x :y))
              (cl-mpm::split-mps-criteria
               *sim*
               (lambda (mp h)
                 (when (and
                        (or
                         (<= (cutout
                              (cl-mpm/fastmaths:fast-.+
                               (cl-mpm/particle:mp-position mp)
                               (cl-mpm/fastmaths:fast-.*
                                (cl-mpm/particle::mp-domain-size mp)
                                (cl-mpm/utils:vector-from-list (list 0.5d0 0.5d0 0d0))))
                              ) (* mesh-size 0d0))
                         (<= (cutout
                              (cl-mpm/fastmaths:fast-.+
                               (cl-mpm/particle:mp-position mp)
                               (cl-mpm/fastmaths:fast-.*
                                (cl-mpm/particle::mp-domain-size mp)
                                (cl-mpm/utils:vector-from-list (list 0.5d0 -0.5d0 0d0))))
                              ) (* mesh-size 0d0))
                         )
                        (> (cutout (cl-mpm/particle:mp-position mp)) (- (* mesh-size 1d0))))
                   dir)))))
          (cl-mpm/setup:remove-sdf *sim*
                                   #'cutout
                                   ))))

      (when nil
        (let ((cut-height (* 0.5d0 shelf-height-true))
              (cut-back-distance 0.15d0)
              (width                    ;(* 2d0 mesh-size)
                ;; 0.5d0
                (* 1d0 (cl-mpm/particle::mp-local-length (aref (cl-mpm:sim-mps *sim*) 0)))
                ))
          (cl-mpm/setup::apply-sdf *sim* (lambda (p) (cl-mpm/setup::line-sdf
                                                      (cl-mpm/utils:vector-from-list (list (magicl:tref p 0 0)
                                                                              (magicl:tref p 1 0)
                                                                              0d0))
                                                      (list (- (magicl:tref sloped-inflex-point 0 0)
                                                               (* cut-back-distance shelf-height-true))
                                                            (float shelf-height 0d0)
                                                            0d0)
                                                      (list (- (magicl:tref sloped-inflex-point 0 0)
                                                               (* cut-back-distance shelf-height-true))
                                                            (float (- shelf-height cut-height) 0d0)
                                                            0d0)
                                                      width))
                                   (lambda (mp v)
                                     (let ((d ;(* 0.99d0 (exp (- (expt (/ (+ width v) width) 1))))
                                             (* 1d0 (cl-mpm/damage::weight-func (expt (+ v width) 2) width)))
                                           )
                                       (setf (cl-mpm/particle:mp-damage mp)
                                             d)
                                       (let ((k (cl-mpm/damage::find-k-damage-mp mp d)))
                                         (setf (cl-mpm/particle::mp-history-stress mp)
                                               k)))
                                     (cl-mpm/damage::update-damage mp 1d-3)
                                     ))
          ))
      (defparameter *to-damage-mps* (list))
      (when nil
        (let ((cut-height (* 0.5d0 shelf-height-true))
              (cut-back-distance 0.5d0)
              (width                    ;(* 2d0 mesh-size)
                ;; 0.5d0
                (* 1d0 (cl-mpm/particle::mp-local-length (aref (cl-mpm:sim-mps *sim*) 0)))
                ))
          (cl-mpm/setup::apply-sdf *sim* (lambda (p) (cl-mpm/setup::line-sdf
                                                      (cl-mpm/utils:vector-from-list (list (magicl:tref p 0 0)
                                                                              (magicl:tref p 1 0)
                                                                              0d0))
                                                      (list (- (magicl:tref sloped-inflex-point 0 0)
                                                               (* cut-back-distance shelf-height-true))
                                                            (float shelf-height 0d0)
                                                            0d0)
                                                      (list (magicl:tref sloped-inflex-point 0 0)
                                                            (float soil-boundary 0d0)
                                                            0d0)
                                                      width))
                                   (lambda (mp v)
                                     (let ((d (* (- 1d0 1d-1)
                                                 (exp
                                                  (- (expt (/ (+ width v) width) 2))))
                                             ;; (* (- 1d0 1d-5) (cl-mpm/damage::weight-func (expt (+ v width) 2) width))
                                             )
                                           )
                                       (when (< v 0d0)
                                         (setf d (- 1d0 1d-9))
                                         )
                                       ;; (setf (cl-mpm/particle:mp-damage mp)
                                       ;;       d)
                                       (let ((k (cl-mpm/damage::find-k-damage-mp mp d)))
                                         (push (list mp k) *to-damage-mps*)
                                         ;; (setf (cl-mpm/particle::mp-history-stress mp)
                                         ;;       k)
                                         ))
                                     ;; (cl-mpm/damage::update-damage mp 1d-3)
                                     ))
          ))
      )

    ;; (let ((notch-height 0.0d0))
    ;;   (cl-mpm/setup:remove-sdf
    ;;    *sim*
    ;;    (lambda (p)
    ;;      (if t
    ;;          (funcall
    ;;           (cl-mpm/setup::rectangle-sdf
    ;;            (list shelf-length (+ notch-height soil-boundary))
    ;;            (list (* 2d0 notch-height) notch-height)
    ;;            ) p)
    ;;          1d0))))

    (setf cl-mpm::*max-split-depth* 4))

    (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
    (defparameter *run-sim* t)
    (defparameter *t* 0)
    (defparameter *oobf* 0)
    (defparameter *energy* 0)
    (defparameter *sim-step* 0)

  (defparameter *data-t* (list))
  (defparameter *data-e* (list))
  (defparameter *data-se* (list))
  (defparameter *data-gpe* (list))
  (defparameter *data-di* (list))
  (defparameter *data-pi* (list))
  (defparameter *data-oobf* (list))

  (defparameter *data-d* (list))
  )
(defparameter *data-t* (list))
(defparameter *data-e* (list))
(defparameter *data-se* (list))
(defparameter *data-gpe* (list))
(defparameter *data-di* (list))
(defparameter *data-pi* (list))
(defparameter *data-oobf* (list))

(defun stop ()
  (setf *run-sim* nil)
  (setf cl-mpm/dynamic-relaxation::*run-convergance* nil))


(defun test-ssr ()
  (let ((phi 50d0)
        (phi-r 30d0)
        (c (* 2d0 26d3))
        (c-r 0d0))
    (setup :refine 2 :mps 3)
    (cl-mpm/dynamic-relaxation::run-load-control
     *sim*
     :output-dir (format nil "./output/")
     :plotter #'plot
     :load-steps 100
     :damping 0.5d0
     :substeps 50
     :criteria 1d-2
     :adaptive-damping nil
     :kinetic-damping nil
     :save-vtk-dr t
     :save-vtk-loadstep t
     :dt-scale 0.8d0
     :loading-function
     (lambda (percent)
       (cl-mpm:iterate-over-mps
        (cl-mpm:sim-mps *sim*)
        (lambda (mp)
          (setf (cl-mpm/particle::mp-c mp)
                (+ c-r (* (- 1d0 percent) (- c c-r))))))))
    ))
