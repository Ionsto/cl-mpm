(defpackage :cl-mpm/examples/chalk
  (:use :cl))
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
;; (setf *block-compile-default* t)
(in-package :cl-mpm/examples/chalk)
(declaim (optimize (debug 3) (safety 3) (speed 0)))

;; (defmethod print-object ((object magicl:matrix) stream)
;;   (pprint object stream))

(declaim (notinline plot))
(defun plot (sim)
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :deformed
   ;; :colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :xx))
   :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-damage-ybar mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-yield-func mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-strain-plastic-vm mp))
   ;; :colour-func #'cl-mpm/particle::mp-strain-plastic-vm
   ;; :colour-func #'cl-mpm/particle::mp-boundary
   ))

(defclass bc-weathering (cl-mpm/buoyancy::bc-scalar)
  ((weathering-rate
    :initform 1d0
    )))
(defun make-bc-weathering (sim datum)
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (with-accessors ((h cl-mpm/mesh:mesh-resolution))
        mesh
      (make-instance 'bc-weathering
                     :index '(0 0 0)
                     :sim sim))))
(defmethod cl-mpm/bc::apply-bc ((bc bc-weathering) node mesh dt)
  (call-next-method)
  (with-accessors ((sim cl-mpm/buoyancy::bc-buoyancy-sim))
      bc
    (when (cl-mpm::sim-enable-damage sim)
      (loop for mp across (cl-mpm:sim-mps sim)
            do
               (let ((weathering 0d0))
                 (cl-mpm:iterate-over-neighbours
                  mesh
                  mp
                  (lambda (mesh mp node svp grads fsvp fgrads)
                    (incf weathering (* svp (cl-mpm/mesh::node-boundary-scalar node)))))
                 (setf (cl-mpm/particle::mp-boundary mp)
                       weathering)
                 (incf
                  (cl-mpm/particle::mp-history-stress mp)
                  (abs (*
                        10d0
                        weathering dt))))
            ))))

(declaim (notinline setup-test-column))
(defun setup-test-column (size block-size offset &optional (e-scale 1) (mp-scale 1))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               :sim-type
               ;; 'cl-mpm::mpm-sim-usf
               'cl-mpm/damage::mpm-sim-damage
               ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 1.7d3)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
         )
    (declare (double-float h density))
    (progn
      (let ()
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-mps-from-list
               (cl-mpm/setup::make-block-mps-list
                ;; (mapcar (lambda (x) 0) size)
                offset
                ;; '(0 0)
                block-size
                (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                density
                ;; ;; ;; 'cl-mpm/particle::particle-chalk-brittle
                'cl-mpm/particle::particle-chalk-delayed
                ;; 'cl-mpm/particle::particle-chalk-anisotropic
                ;; ;; 'cl-mpm/particle::particle-limestone
                :E 1d9
                :nu 0.35d0
                ;; :rho 1d6
                :enable-plasticity t

                :ft 200d3
                :fc 500d3
                :friction-angle 40d0

                :fracture-energy 3000d0
                :initiation-stress 200d3
                :delay-time 1d1
                :ductility 10d0
                ;; :compression-ratio 8d0

                :critical-damage 1.0d0;0.999d0
                :damage-domain-rate 0.5d0;This slider changes how GIMP update turns to uGIMP under damage
                :local-length (* 3.0d0 (sqrt 7))
                ;; :local-length-damaged (* 0.1d0 (sqrt 7))
                ;; :local-length-damaged 1d0
                :local-length-damaged 10d-10

                ;; 'cl-mpm/particle::particle-mc
                ;; :E 1d9
                ;; :nu 0.35d0
                :psi (* 00d0 (/ pi 180))
                :phi (* 40d0 (/ pi 180))
                :c 1000d3
                ;; :c 1d6

                ;; 'cl-mpm/particle::particle-vm
                ;; :E 1d9
                ;; :nu 0.35d0
                ;; :rho 1d6

                :gravity -9.8d0
                :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
                ))))
      ;; (let* ((mp-0 (aref (cl-mpm:sim-mps *sim*) 0))
      ;;        (fc (cl-mpm/particle::mp-fc mp-0))
      ;;        )
      ;;   (format t "Chalk damage angle: ~F~%"
      ;;           (atan (* 3 (/ (- fc ft) (+ fc ft))))))
      ;; (cl-mpm/examples/tpb::calculate-ductility-param 1d9 200d0 1d0 200d3)
      (setf (cl-mpm:sim-allow-mp-split sim) t)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-enable-fbar sim) t)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (setf (cl-mpm::sim-mass-filter sim) 1d-15)
      (let ((ms 1d4))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        (setf (cl-mpm:sim-damping-factor sim) (* 1d-1 ms))
        ;; (setf (cl-mpm:sim-damping-factor sim) 10.0d0)
        ;; (setf (cl-mpm:sim-damping-factor sim) (* 1d-1 ms))
        )

      ;; (dotimes (i 2)
      ;;   (dolist (dir (list :x :y))
      ;;     (cl-mpm::split-mps-criteria
      ;;      sim
      ;;      (lambda (mp h)
      ;;        (when
      ;;            (and
      ;;             (> (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
      ;;                80)
      ;;             (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
      ;;                200
      ;;                )
      ;;             (> (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
      ;;                50
      ;;                )
      ;;             )
      ;;          dir
      ;;          )))))

      (let ((dt-scale 1d0))
        (setf
         (cl-mpm:sim-dt sim)
         (* dt-scale h
            (sqrt (cl-mpm::sim-mass-scale sim))
            (sqrt (/ density (cl-mpm/particle::mp-p-modulus (aref (cl-mpm:sim-mps sim) 0)))))))

      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))
      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-var
             (cl-mpm:sim-mesh sim)
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 0 nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 0 nil)))
             ;; (lambda (i) nil)
             ;; (lambda (i) nil)
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
            ))

      (defparameter *floor-bc*
        (cl-mpm/penalty::make-bc-penalty-point-normal
         sim
         (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
         (cl-mpm/utils:vector-from-list (list 00d0 (+ h-y) 0d0))
         (* density 1d5)
         0.9d0
         ;; 1d1
         ))

      ;; (setf (cl-mpm::sim-bcs-force-list sim)
      ;;       (list
      ;;        (cl-mpm/bc:make-bcs-from-list
      ;;         (list
      ;;          (make-bc-weathering sim 0d0)
      ;;          ;; *floor-bc*
      ;;          ))))

      ;; (let* ((terminus-size (second block-size))
      ;;        (ocean-y (* terminus-size 1.0d0)))
      ;;   (setf (cl-mpm::sim-bcs-force-list sim)
      ;;         (list
      ;;          (cl-mpm/bc:make-bcs-from-list
      ;;           (list
      ;;            (cl-mpm/buoyancy::make-bc-buoyancy-clip
      ;;             sim
      ;;             ocean-y
      ;;             1000d0
      ;;             (lambda (pos datum)
      ;;               t)
      ;;             ))))))

      sim)))

(defmethod cl-mpm::post-stress-step (mesh (mp cl-mpm/particle::particle-mc) dt)
  (with-accessors ((ps cl-mpm/particle::mp-strain-plastic-vm)
                   (c cl-mpm/particle::mp-c)
                   (phi cl-mpm/particle::mp-phi)
                   )
      mp
    ;; (break)
    (let ((rho_0 1000d3)
          (rho_1 0d0)
          (phi_0 (* 30d0 (/ pi 180)))
          (phi_1 (* 20d0 (/ pi 180)))
          (soft 10d0))
      ;; (setf c (+ rho_1
      ;;            (* (- rho_0 rho_1) (exp (- (* soft ps))))))
      ;; (setf
      ;;  phi
      ;;  (+ phi_1
      ;;     (* (- phi_0 phi_1) (exp (- (* soft ps))))))
      )
    )
  )

(defmethod cl-mpm::post-stress-step (mesh (mp cl-mpm/particle::particle-vm) dt)
  (with-accessors ((ps cl-mpm/particle::mp-strain-plastic-vm)
                   (c cl-mpm/particle::mp-rho))
      mp
    ;; (let ((rho_0 1d5)
    ;;       (rho_1 0.0100d6)
    ;;       (soft 5d1))
    ;;   (setf c (max rho_1
    ;;                (* rho_0 (exp (- (* soft ps)))))))
    )
  )

(defparameter *sim* nil)
(defparameter *run-sim* t)
(defparameter *t* 0)
(defparameter *sim-step* 0)
(defparameter *refine* (/ 1d0 2d0))
(let ((refine (uiop:getenv "REFINE")))
  (when refine
    (setf *refine* (parse-integer (uiop:getenv "REFINE")))
    ))

(defun setup (&key (undercut 0d0))
  (let* ((mesh-size 5)
         (mps-per-cell 2)
         (shelf-height 100)
         (soil-boundary 20)
         (shelf-aspect 2.0)
         (runout-aspect 2.00)
         (shelf-length (* shelf-height shelf-aspect))
         (domain-length (+ shelf-length (* runout-aspect shelf-height)))
         (shelf-height (+ shelf-height soil-boundary))
         (depth 400)
         (offset (list 0 (* 0 mesh-size)
                       ;; 0
                       ))
         )
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

    ;;Refine around tip
    (dotimes (i 0)
      (dolist (dir (list :x
                         :y
                         ))
        (cl-mpm::split-mps-criteria
         *sim*
         (lambda (mp h)
           (when
               (and
                (> (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                   (- shelf-length (* 1d0 shelf-height)))
                (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                   shelf-length)
                (> (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                   soil-boundary
                   )
                )
             dir
             )))))
    (loop for mp across (cl-mpm::sim-mps *sim*)
          do (setf (cl-mpm/particle::mp-split-depth mp) 0))

    (let* ((undercut-angle undercut)
           (normal (magicl:from-list (list
                                      (cos (- (* pi (/ undercut-angle 180d0))))
                                      (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1))))
      (cl-mpm/setup:remove-sdf *sim*
                               (lambda (p)
                                 (if (> (magicl:tref p 1 0) soil-boundary)
                                     (cl-mpm/setup::plane-point-sdf
                                      (magicl:from-list (list (magicl:tref p 0 0)
                                                              (magicl:tref p 1 0)) '(2 1))
                                      normal
                                      (magicl:from-list (list shelf-length soil-boundary)
                                                        '(2 1) :type 'double-float))
                                     1d0)
                                 ))
      )

   
     (let ((notch-height 0d0))
       (cl-mpm/setup:remove-sdf *sim*
                                (lambda (p)
                                  (if t ;(< (abs (- (magicl:tref p 2 0) (* 0.5d0 depth))) 20d0)
                                      (funcall
                                       (cl-mpm/setup::rectangle-sdf
                                        (list shelf-length (+ notch-height soil-boundary))
                                        (list (* 2d0 notch-height) notch-height)
                                        ) p)
                                      1d0)
                                  )))

     (setf cl-mpm::*max-split-depth* 6)

     ;; (let ((ratio 1.0d0))
     ;;   (cl-mpm/setup::damage-sdf *sim* (lambda (p) (cl-mpm/setup::line-sdf
     ;;                                                (magicl:from-list (list (magicl:tref p 0 0)
     ;;                                                                        (magicl:tref p 1 0)) '(2 1))
     ;;                                                (list (- shelf-length (* shelf-height ratio)) shelf-height)
     ;;                                                (list shelf-length soil-boundary)
     ;;                                                10d0
     ;;                                                )) 1.0d0))
     )
    ;; (let ((upper-random-bound 0.5d0))
    ;;   (loop for mp across (cl-mpm:sim-mps *sim*)
    ;;         do (setf (cl-mpm/particle::mp-damage mp)
    ;;                  (reduce #'*
    ;;                          (loop for i from 0 to 2
    ;;                                collect (random upper-random-bound))))))
    (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
    (defparameter *run-sim* t)
    (defparameter *t* 0)
    (defparameter *oobf* 0)
    (defparameter *energy* 0)
    (defparameter *sim-step* 0))

(defgeneric estimate-energy-crit (sim))

(defmethod estimate-energy-crit ((sim cl-mpm::mpm-sim))
  ;; (/)
  (loop for mp across (cl-mpm:sim-mps sim)
            maximizing (* (cl-mpm/particle::mp-mass mp)
                   ;; (cl-mpm/particle::mp-damage mp)
                   (cl-mpm/fastmaths::mag-squared (cl-mpm/particle::mp-velocity mp))))
  ;; (cl-mpm::sim-mass-scale sim)
  ;; (* (loop for mp across (cl-mpm:sim-mps *sim*)
  ;;          summing (cl-mpm/particle::mp-mass mp))
  ;;    ;; (cl-mpm::sim-mass-scale sim)
  ;;    )
  
   )
(defmethod estimate-energy-crit ((sim cl-mpm/mpi::mpm-sim-mpi))
  (cl-mpm/mpi::mpi-sum
   (loop for mp across (cl-mpm:sim-mps sim)
         maximizing (* (cl-mpm/particle::mp-mass mp)
                (cl-mpm/fastmaths::mag-squared (cl-mpm/particle::mp-velocity mp)))))
  ;; (/
  ;;  (cl-mpm/mpi::mpi-sum
  ;;   (loop for mp across (cl-mpm:sim-mps *sim*)
  ;;                            summing (* (cl-mpm/particle::mp-mass mp)
  ;;                                       (cl-mpm/fastmaths::mag-squared (cl-mpm/particle::mp-velocity mp)))))
  ;;  (* (cl-mpm/mpi:mpi-sum
  ;;      (loop for mp across (cl-mpm:sim-mps *sim*)
  ;;            summing (cl-mpm/particle::mp-mass mp)))
  ;;     (cl-mpm::sim-mass-scale sim)))
  )

(defun estimate-oobf (sim)
  (let ((oobf 0d0)
        (nmax 0d0)
        (dmax 0d0)
        (iter 0))
    (cl-mpm::iterate-over-nodes-serial
     (cl-mpm:sim-mesh sim)
     (lambda (node)
       (with-accessors ((active cl-mpm/mesh::node-active)
                        (f-ext cl-mpm/mesh::node-external-force)
                        (f-int cl-mpm/mesh::node-internal-force))
           node
         (when active
           ;; (setf imax iter)
           (setf nmax (+ nmax
                         (cl-mpm/fastmaths::mag-squared
                          (magicl:.- f-ext f-int)))
                 dmax (+ dmax (cl-mpm/fastmaths::mag-squared f-ext))))
         )
       (incf iter)
       ))
    (when (> dmax 0d0)
      ;; (pprint (row-major-aref (cl-mpm/mesh::mesh-nodes (cl-mpm:sim-mesh *sim*)) imax))
      ;; (break)
      (setf oobf (/ nmax dmax)))
    oobf
    ;; (when (and ;(< fnorm 1d-6) (< oobf 5d0)
    ;;        (< (abs (- penalty reaction)) 100d0)
    ;;        )
    ;;   (setf converged t))
    ))

(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)

  (let* ((target-time 1d1)
         (target-time-original target-time)
         (mass-scale (cl-mpm::sim-mass-scale *sim*))
         (collapse-target-time 1d0)
         (collapse-mass-scale 1d0)
         (plasticity-enabled (cl-mpm/particle::mp-enable-plasticity (aref (cl-mpm:sim-mps *sim*) 0)))
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 1.0d0)
         (settle-steps 10)
         (damp-steps 5)
         )
    (cl-mpm/output::save-simulation-parameters #p"output/settings.json"
                                               *sim*
                               

                (list :dt target-time))
    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp) (setf (cl-mpm/particle::mp-enable-plasticity mp) nil)))

    (cl-mpm::update-sim *sim*)
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                    (format t "CFL dt estimate: ~f~%" dt-e)
                    (format t "CFL step count estimate: ~D~%" substeps-e)
                    (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
    (time (loop for steps from 0 to 1000
                while *run-sim*
                do
                   (progn
                     ;; (when (= steps settle-steps)
                     ;;   (setf (cl-mpm::sim-enable-damage *sim*) t)
                     ;;   (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                     ;;     (format t "CFL dt estimate: ~f~%" dt-e)
                     ;;     (format t "CFL step count estimate: ~D~%" substeps-e)
                     ;;     (setf substeps substeps-e))
                     ;;   )
                     (format t "Step ~d ~%" steps)
                     (format t "MPs ~d ~%" (length (cl-mpm:sim-mps *sim*)))
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     ;; (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     (let ((energy-estimate 0d0)
                           (oobf 0d0)
                           )
                       (time
                        (let ((current-damage (cl-mpm::sim-enable-damage *sim*)))
                          (dotimes (i substeps)
                            (cl-mpm::update-sim *sim*)
                            (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)));)
                            (incf energy-estimate (estimate-energy-crit *sim*))
                            (incf oobf (estimate-oobf *sim*))
                            )

                          ))
                       (setf energy-estimate (/ energy-estimate substeps)
                             oobf (/ oobf substeps)
                             )
                       (format t "Energy estimate: ~E~%" energy-estimate)
                       (defparameter *oobf* oobf)
                       (defparameter *energy* energy-estimate)
                       (when (>= steps damp-steps)
                         (let ((ms (cl-mpm::sim-mass-scale *sim*)))
                           (setf (cl-mpm:sim-damping-factor *sim*)
                                 0d-2
                                 ;; 0d0
                                 ;(* 1d-1 ms)
                                 )))
                       (when (>= steps settle-steps)
                           (setf (cl-mpm::sim-enable-damage *sim*) t)
                           (cl-mpm::iterate-over-mps
                            (cl-mpm:sim-mps *sim*)
                            (lambda (mp) (setf (cl-mpm/particle::mp-enable-plasticity mp) plasticity-enabled)))

                         (if (and
                              ;; t
                              (> energy-estimate 1d0)
                              ;; nil
                              )
                             (progn
                               (format t "Collapse timestep~%")
                               (setf
                                target-time collapse-target-time
                                (cl-mpm::sim-mass-scale *sim*) collapse-mass-scale
                                ))
                             (progn
                               (format t "Accelerated timestep~%")
                               ;; (setf
                               ;;  target-time target-time-original
                               ;;  ;; (cl-mpm::sim-mass-scale *sim*) 1d3
                               ;;  (cl-mpm::sim-mass-scale *sim*) mass-scale

                               ;;  )
                               ))))

                     (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       (format t "CFL dt estimate: ~f~%" dt-e)
                       (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf substeps substeps-e))
                     ;; (setf (cl-mpm:sim-damping-factor *sim*)
                     ;;       (* (cl-mpm:sim-damping-factor *sim*) (expt 1d-3 1/40)))

                     (incf *sim-step*)
                     (plot *sim*)
                     (vgplot:title (format nil "Time:~F - OOBF ~E - Energy ~E"  *t* *oobf* *energy*))
                     (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                        :terminal "png size 1920,1080"
                                        )

                     (swank.live:update-swank)
                     ))))
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*))


;; (setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))
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

(defun test ()
  (setup)
  (sb-profile:profile "CL-MPM")
  (sb-profile:profile "CL-MPM/PARTICLE")
  (sb-profile:profile "CL-MPM/MESH")
  (sb-profile:profile "CL-MPM/SHAPE-FUNCTION")
  (sb-profile:reset)
  (time
   (dotimes (i 100)
         (cl-mpm::update-sim *sim*)))
  (sb-profile:report)
  )
(defun test-undercut ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output_chalk/mesh.vtk")
                               *sim*)
  (loop for c in (list 0d0 10d0 20d0 30d0 40d0 50d0)
        while *run-sim*
        do
           (progn
             (setup :undercut (- c))
             (run)
             (cl-mpm/output:save-vtk (merge-pathnames (format nil "output_chalk/chalk_~5,'0d.vtk" c)) *sim*)
             )))

(defun profile ()
  (sb-profile:unprofile)
  (sb-profile:reset)
  (sb-profile:profile "CL-MPM")
  (sb-profile:profile "CL-MPM/PARTICLE")
  (sb-profile:profile "CL-MPM/MESH")
  (loop repeat 100
        do (progn
             (cl-mpm::update-sim *sim*)
             ;; (cl-mpm/damage::calculate-damage (cl-mpm:sim-mesh *sim*)
             ;;                                  (cl-mpm:sim-mps *sim*)
             ;;                                  (cl-mpm:sim-dt *sim*)
             ;;                                  25d0)
             ;; (cl-mpm/eigenerosion:update-fracture *sim*)
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


(defun plot-interaction ()
  (vgplot:close-all-plots)
  (let* ((length 1d0)
         (x (loop for x from -5d0 to 5d0 by 0.1d0 collect x)))
    (vgplot:plot x (mapcar (lambda (x) (cl-mpm/damage::weight-func (expt x 2) length)) x))
    ))




(defun merge-mps (sim length)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (let (main-merge-lock (sb-thread:make-mutex))
      (cl-mpm:iterate-over-mps
       mps
       (lambda (mp)
         (let ((merge-list (list)))
           (cl-mpm/damage::iterate-over-neighour-mps
            mesh
            mp
            length
            (lambda (mesh mp mp-other)
              (push mp-other merge-list)
              ))))))))



(defun setup-joss ()
  (let* ((mesh-size 1.0)
         (mps-per-cell 4)
         (shelf-height 15.4)
         (soil-boundary 0)
         (shelf-aspect 0.5)
         (runout-aspect 1.00)
         (shelf-length (* shelf-height shelf-aspect))
         (domain-length (+ shelf-length (* runout-aspect shelf-height)))
         (shelf-height (+ shelf-height soil-boundary))
         (depth 400)
         (offset (list 0 (* 0 mesh-size)
                       ;; 0
                       ))
         )
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

    ;;Refine around tip
    (dotimes (i 0)
      (dolist (dir (list :x
                         :y
                         ))
        (cl-mpm::split-mps-criteria
         *sim*
         (lambda (mp h)
           (when
               (and
                (> (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                   (- shelf-length (* 1d0 shelf-height)))
                (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                   shelf-length)
                (> (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                   soil-boundary
                   )
                )
             dir
             )))))
    (loop for mp across (cl-mpm::sim-mps *sim*)
          do (setf (cl-mpm/particle::mp-split-depth mp) 0))

    (let* ((sloped-height (- shelf-height 7.8d0))
           (undercut-angle (- 82.5d0 90d0))
           (normal (magicl:from-list (list
                                      (cos (- (* pi (/ undercut-angle 180d0))))
                                      (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1)))
           (sloped-inflex-point
             (magicl:from-list (list (- shelf-length (* (tan (- (* pi (/ undercut-angle 180d0)))) sloped-height))
                                     soil-boundary)
                               '(2 1) :type 'double-float)

             )
           )
      (cl-mpm/setup:remove-sdf *sim*
                               (lambda (p)
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
                                     1d0)
                                 ))
      )

   
     (let ((notch-height 0d0))
       (cl-mpm/setup:remove-sdf *sim*
                                (lambda (p)
                                  (if t ;(< (abs (- (magicl:tref p 2 0) (* 0.5d0 depth))) 20d0)
                                      (funcall
                                       (cl-mpm/setup::rectangle-sdf
                                        (list shelf-length (+ notch-height soil-boundary))
                                        (list (* 2d0 notch-height) notch-height)
                                        ) p)
                                      1d0)
                                  )))

     (setf cl-mpm::*max-split-depth* 6)

     ;; (let ((ratio 1.0d0))
     ;;   (cl-mpm/setup::damage-sdf *sim* (lambda (p) (cl-mpm/setup::line-sdf
     ;;                                                (magicl:from-list (list (magicl:tref p 0 0)
     ;;                                                                        (magicl:tref p 1 0)) '(2 1))
     ;;                                                (list (- shelf-length (* shelf-height ratio)) shelf-height)
     ;;                                                (list shelf-length soil-boundary)
     ;;                                                10d0
     ;;                                                )) 1.0d0))
     )
    ;; (let ((upper-random-bound 0.5d0))
    ;;   (loop for mp across (cl-mpm:sim-mps *sim*)
    ;;         do (setf (cl-mpm/particle::mp-damage mp)
    ;;                  (reduce #'*
    ;;                          (loop for i from 0 to 2
    ;;                                collect (random upper-random-bound))))))
    (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
    (defparameter *run-sim* t)
    (defparameter *t* 0)
    (defparameter *oobf* 0)
    (defparameter *energy* 0)
    (defparameter *sim-step* 0))

(defun plot-mohr-c (stress)
  (let ((angle 30)
        (coheasion 1))
    (vgplot:figure)
    (vgplot:plot '(0) '(0))
    (vgplot:format-plot t "replot (~f*x + ~f)~%" (sin (* angle (/ pi 180))) coheasion)
    (vgplot:replot))
  )
