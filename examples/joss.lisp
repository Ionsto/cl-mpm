(defpackage :cl-mpm/examples/joss
  (:use :cl))
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)
(in-package :cl-mpm/examples/joss)
(declaim (optimize (debug 3) (safety 3) (speed 0)))

(defun cl-mpm/damage::length-localisation (local-length local-length-damaged damage)
  ;; (+ (* local-length (- 1d0 damage)) (* local-length-damaged damage))
  ;; (* local-length (max (sqrt (- 1d0 damage)) 1d-10))
  local-length
  )
;; (defmethod print-object ((object magicl:matrix) stream)
;;   (pprint object stream))

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-chalk-delayed) dt fbar)
  ;; (incf (cl-mpm/particle::mp-c mp) 1d-3)
  (cl-mpm::update-stress-kirchoff-damaged mesh mp dt fbar)
  ;; (cl-mpm::update-stress-kirchoff mesh mp dt fbar)
  )
(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-chalk-brittle) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
                     (strain cl-mpm/particle::mp-strain)
                     (damage cl-mpm/particle:mp-damage)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (pressure cl-mpm/particle::mp-pressure)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     (angle cl-mpm/particle::mp-friction-angle)
                     (c cl-mpm/particle::mp-coheasion)
                     (nu cl-mpm/particle::mp-nu)
                     (ft cl-mpm/particle::mp-ft)
                     (fc cl-mpm/particle::mp-fc)
                     (E cl-mpm/particle::mp-e)
                     (de cl-mpm/particle::mp-elastic-matrix)
                     (kc-r cl-mpm/particle::mp-k-compressive-residual-ratio)
                     (kt-r cl-mpm/particle::mp-k-tensile-residual-ratio)
                     (g-r cl-mpm/particle::mp-shear-residual-ratio)
                     ) mp
      (declare (double-float pressure damage))
      (progn
        (when (< damage 1d0)
          (setf damage-increment (cl-mpm/damage::tensile-energy-norm strain E de))
          ;; (multiple-value-bind (s_1 s_2 s_3) (cl-mpm/damage::principal-stresses-3d (magicl:scale stress (/ 1d0 (magicl:det def))))
          ;;   (setf damage-increment (max damage-increment (max s_1 0d0)))
          ;;   )
          ;(setf damage-increment (criterion-mc strain (* angle (/ pi 180d0)) E nu))
          ;; (setf damage-increment (criterion-dp strain (* angle (/ pi 180d0)) E nu))
          ;; (setf damage-increment (cl-mpm/damage::smooth-rankine-criterion (magicl:scale stress (/ 1d0 (magicl:det def)))))
          ;; (setf damage-increment
          ;;       (max 0d0
          ;;            (cl-mpm/damage::drucker-prager-criterion
          ;;             (magicl:scale stress (/ 1d0 (magicl:det def))) (* angle (/ pi 180d0)))))
          )
        (when (>= damage 1d0)
          (setf damage-increment 0d0))
        ;;Delocalisation switch
        (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
        (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
        ))))

;; (defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-chalk-brittle) dt)
;;   (let ((damage-increment 0d0))
;;     (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
;;                      (strain cl-mpm/particle::mp-strain)
;;                      (damage cl-mpm/particle:mp-damage)
;;                      (init-stress cl-mpm/particle::mp-initiation-stress)
;;                      (critical-damage cl-mpm/particle::mp-critical-damage)
;;                      (damage-rate cl-mpm/particle::mp-damage-rate)
;;                      (pressure cl-mpm/particle::mp-pressure)
;;                      (ybar cl-mpm/particle::mp-damage-ybar)
;;                      (def cl-mpm/particle::mp-deformation-gradient)
;;                      (angle cl-mpm/particle::mp-friction-angle)
;;                      (c cl-mpm/particle::mp-coheasion)
;;                      (nu cl-mpm/particle::mp-nu)
;;                      (ft cl-mpm/particle::mp-ft)
;;                      (fc cl-mpm/particle::mp-fc)
;;                      (E cl-mpm/particle::mp-e)
;;                      (de cl-mpm/particle::mp-elastic-matrix)
;;                      ) mp
;;       (declare (double-float pressure damage))
;;         (progn
;;           (when (< damage 1d0)
;;             (let* ((strain+
;;                      (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix strain))
;;                        (loop for i from 0 to 2
;;                              do
;;                                 (setf (nth i l) (max (nth i l) 0d0)))
;;                        (cl-mpm/utils:matrix-to-voigt (magicl:@ v
;;                                                                (magicl:from-diag l :type 'double-float)
;;                                                                (magicl:transpose v)))))
;;                    (strain- (magicl:.- strain strain+))
;;                    (k (/ fc ft))
;;                    (e+ (sqrt (max 0d0 (* E (cl-mpm/fastmath::dot strain+ (magicl:@ de strain+))))))
;;                    (e- (sqrt (max 0d0 (* E (cl-mpm/fastmath::dot strain- (magicl:@ de strain-)))))))

;;               (setf damage-increment
;;                     (if (= *sim-value* 0)
;;                         e+
;;                         (/
;;                           (+ (* k e+) e-)
;;                           (+ k 1d0)
;;                           ))
;;                     ))

;;             )
;;           (when (>= damage 1d0)
;;             (setf damage-increment 0d0))
;;           ;;Delocalisation switch
;;           (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
;;           (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)))))


(declaim (notinline plot))
(defun plot (sim)
  ;; (plot-collapse-curves)
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :deformed
   ;; :colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :xx))
   :colour-func #'cl-mpm/particle::mp-damage
   ;; :colour-func #'cl-mpm/particle::mp-damage-ybar
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-yield-func mp))
   ;; :colour-func (lambda (mp) (magicl:tref (cl-mpm/particle::mp-stress mp) 0 0))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-strain-plastic-vm mp))
   ;; :colour-func #'cl-mpm/particle::mp-strain-plastic-vm
   ;; :colour-func #'cl-mpm/particle::mp-boundary
   ;; :colour-func (lambda (mp) (*
   ;;                            ;; (- 1d0 (cl-mpm/particle::mp-damage mp))
   ;;                            (cl-mpm/particle::mp-boundary mp)
   ;;                            ))
   )
  )

(defclass bc-weathering (cl-mpm/buoyancy::bc-scalar)
  ((weathering-rate
    :initform 1d0
    )))
(declaim (notinline make-bc-weathering))
(defun make-bc-weathering (sim datum)
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (with-accessors ((h cl-mpm/mesh:mesh-resolution))
        mesh
      (make-instance 'bc-weathering
                     :index '(0 0 0)
                     :sim sim
                     :scalar-func (lambda (pos)
                                    1d0
                                    ;; (+ 1.0d0
                                    ;;    (* 1d0 (exp (- (expt (- 3.5d0 (magicl:tref pos 1 0)) 2))))
                                    ;;    )
                                    )
                     ))))
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
                    (when (cl-mpm/mesh:node-active node)
                      (incf weathering (* svp (cl-mpm/mesh::node-boundary-scalar node))))))
                 (setf (cl-mpm/particle::mp-boundary mp)
                       weathering)
                 ;; (setf weathering (cl-mpm/particle::mp-boundary mp))
                 (incf
                  (cl-mpm/particle::mp-history-stress mp)
                  (abs (*
                        10d3
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
         )
    (declare (double-float h density))
    (progn
      (let ()
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-mps-from-list
               (cl-mpm/setup::make-block-mps-list
                offset
                block-size
                (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                density
                'cl-mpm/particle::particle-chalk-delayed
                :E 1d9
                :nu 0.24d0
                ;; :nu 0.15d0

                :enable-plasticity t

                :ft 1d0
                :fc 10d0

                :friction-angle 60.0d0
                :kt-res-ratio 1d-10
                :kc-res-ratio 1d-2
                :g-res-ratio 6.5d-3

                ;; :kt-res-ratio 1d-10
                ;; :kc-res-ratio 1d-2
                ;; :g-res-ratio  1d-3

                ;; :g-res-ratio 1.0d-3
                ;; :g-res-ratio 1.0d-4
                :fracture-energy 3000d0
                :initiation-stress 30d3;18d3
                :delay-time 10d0
                :delay-exponent 3d0
                ;:ductility 6.7d0
                :ductility 10d0
                ;; :ductility 6.7d0

                :critical-damage 1d0;(- 1.0d0 1d-3)
                :damage-domain-rate 0.9d0;This slider changes how GIMP update turns to uGIMP under damage

                :local-length 1.0d0;(* 0.20d0 (sqrt 7))
                :local-length-damaged 10d-10

                :psi (* 00d0 (/ pi 180))
                :phi (* 42d0 (/ pi 180))
                ;; :phi 0d0;(* 30d0 (/ pi 180))
                ;; :phi (* 70d0 (/ pi 180))
                :c 131d3
                ;; :c 1000d3
                ;; :c 500d3

                :gravity -9.8d0
                :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
                ))))
      (let* ((mp-0 (aref (cl-mpm:sim-mps sim) 0))
             (fc (cl-mpm/particle::mp-fc mp-0))
             (ft (cl-mpm/particle::mp-ft mp-0))
             (angle-d (* (/ 180 pi) (atan (* 3 (/ (- fc ft) (+ fc ft))))))
             (rc (cl-mpm/particle::mp-k-compressive-residual-ratio mp-0))
             (rs (cl-mpm/particle::mp-shear-residual-ratio mp-0))
             (angle-plastic (cl-mpm/particle::mp-phi mp-0))
             (angle-plastic-damaged (atan (* (/ rs rc) (tan angle-plastic))))
             )
        (format t "Estimated Gf ~F~%" (estimate-gf (cl-mpm/particle::mp-ductility mp-0)
                                                   (cl-mpm/particle::mp-initiation-stress mp-0)
                                                   (cl-mpm/particle::mp-local-length mp-0)))
        (format t "Chalk damage growth angle: ~F~%"
                angle-d)
        (format t "Chalk plastic virgin angle: ~F~%"
                (* (/ 180 pi) angle-plastic))
        (format t "Chalk plastic residual angle: ~F~%"
                (* (/ 180 pi) angle-plastic-damaged)))

      ;; (cl-mpm/examples/tpb::calculate-ductility-param 1d9 47.5d0 0.05d0 20d3)

      (setf (cl-mpm:sim-allow-mp-split sim) t)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-enable-fbar sim) t)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (setf (cl-mpm::sim-mass-filter sim) 0d0)
      (let ((ms 1d6))
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
      ;;          (make-bc-weathering sim 10d0)
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


(defun estimate-gf (eta ft lc &optional (E 1d9))
  (let* (
         (gf (* eta (/ (expt ft 2) (* 2 E))))
         )
    (* gf lc)))


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

(defun plot-collapse-curves ()
  (vgplot:figure)
  (vgplot:plot
   ;; *data-t* *data-e*   "Energy"
   *data-t* (mapcar (lambda (x) (/ x 1d0 ;(reduce #'max *data-e*)
                                   )) *data-e*  ) "Kinetic Energy"
   *data-t* (mapcar (lambda (x) (/ x 1d0
                                   ;(reduce #'max *data-se*)
                                   )) *data-se*  ) "Strain Energy"
   ;; *data-t* (mapcar (lambda (x) (/ x (reduce #'max *data-di*))) *data-di* ) "Damage inc"
   ;; *data-t* (mapcar (lambda (x) (/ x (reduce #'max *data-pi*))) *data-pi* ) "plastic inc"
   )
  )

(let ((refine (uiop:getenv "REFINE")))
  (when refine
    (setf *refine* (parse-integer (uiop:getenv "REFINE")))
    ))


(defun estimate-total-energy (sim)
  (loop for mp across (cl-mpm:sim-mps sim)
        sum (* (cl-mpm/particle::mp-mass mp)
                      (cl-mpm/fastmath::mag-squared (cl-mpm/particle::mp-velocity mp)))))

(defgeneric estimate-energy-crit (sim))

(defmethod estimate-energy-crit ((sim cl-mpm::mpm-sim))
  ;; (loop for mp across (cl-mpm:sim-mps sim)
  ;;           sum (*
  ;;                0.5d0
  ;;                (cl-mpm/particle::mp-mass mp)
  ;;                  ;; (cl-mpm/particle::mp-damage mp)
  ;;                  (cl-mpm/fastmath::mag-squared (cl-mpm/particle::mp-velocity mp))))
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
    energy)
  )
(defmethod estimate-energy-crit ((sim cl-mpm/mpi::mpm-sim-mpi))
  ;; (cl-mpm/mpi::mpi-sum
  ;;  (loop for mp across (cl-mpm:sim-mps sim)
  ;;        sum (* (cl-mpm/particle::mp-mass mp)
  ;;               (cl-mpm/fastmath::mag-squared (cl-mpm/particle::mp-velocity mp)))))
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
                         (cl-mpm/fastmath::mag-squared
                          (magicl:.- f-ext f-int)))
                 dmax (+ dmax (cl-mpm/fastmath::mag-squared f-ext))))
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

  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
  (let* ((target-time 1d1)
         (target-time-original target-time)
         (mass-scale (cl-mpm::sim-mass-scale *sim*))
         (collapse-target-time 0.1d0)
         (collapse-mass-scale 1d0)
         (plasticity-enabled (cl-mpm/particle::mp-enable-plasticity (aref (cl-mpm:sim-mps *sim*) 0)))
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 1.0d0)
         (settle-steps 6)
         (damp-steps 5)
         (dt-0 0d0)
         (last-e 0d0)
         )
    (cl-mpm/output::save-simulation-parameters
     #p"output/settings.json"
     *sim*
     (list :dt target-time))
    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp) (setf (cl-mpm/particle::mp-enable-plasticity mp) nil)))
    (with-open-file (stream (merge-pathnames "output/timesteps.csv") :direction :output :if-exists :supersede)
      (format stream "steps,time,energy~%"))

    (cl-mpm::update-sim *sim*)
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
      (format t "CFL dt estimate: ~f~%" dt-e)
      (format t "CFL step count estimate: ~D~%" substeps-e)
      (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
    (setf dt-0 (/ (cl-mpm:sim-dt *sim*) (sqrt (cl-mpm::sim-mass-scale *sim*))))
    (time (loop for steps from 0 to 500
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
                     (with-open-file (stream (merge-pathnames "output/timesteps.csv") :direction :output :if-exists :append)
                       (format stream "~D,~f~%"
                               steps
                               *t*
                               last-e
                               ))
                     ;; (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     (let ((energy-estimate 0d0)
                           (oobf 0d0)
                           (total-energy 0d0)
                           (total-strain-energy 0d0)
                           (total-gpe 0d0)
                           (damage-inc 0d0)
                           (plastic-inc 0d0)
                           )
                       (time
                        (let ((current-damage (cl-mpm::sim-enable-damage *sim*)))
                          (loop for i from 1 to substeps
                                while *run-sim*
                                do (progn
                                     (cl-mpm::update-sim *sim*)
                                     (setf *t* (+ *t* (cl-mpm::sim-dt *sim*))) ;)
                                     ;; (incf energy-estimate (estimate-energy-crit *sim*))
                                     ;(incf total-energy (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
                                     (incf total-energy (estimate-energy-crit *sim*))
                                     (incf total-strain-energy
                                           (lparallel:pmap-reduce (lambda (mp)
                                                                    (* 0.5d0
                                                                       (cl-mpm/particle::mp-volume mp)
                                                                       (cl-mpm/fastmath:dot
                                                                        (cl-mpm/particle::mp-stress mp)
                                                                        (cl-mpm/particle::mp-strain mp)
                                                                        )))
                                                                  #'+ (cl-mpm:sim-mps *sim*))
                                           
                                           )
                                     (incf total-gpe
                                           (lparallel:pmap-reduce (lambda (mp)
                                                                    (* 
                                                                     (cl-mpm/particle:mp-mass mp)
                                                                     (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                                                                     (cl-mpm/particle::mp-gravity mp)))
                                                                  #'+ (cl-mpm:sim-mps *sim*))
                                           )
                                     ;; (incf damage-inc
                                     ;;       (loop for mp across (cl-mpm:sim-mps *sim*)
                                     ;;             sum
                                     ;;             (cl-mpm/particle::mp-damage-increment mp)))
                                     ;; (incf plastic-inc
                                     ;;       (loop for mp across (cl-mpm:sim-mps *sim*)
                                     ;;             sum
                                     ;;             (cl-mpm/particle::mp-strain-plastic-vm-inc mp)))
                                     (incf oobf (estimate-oobf *sim*))))))
                       (setf energy-estimate (/ energy-estimate substeps)
                             oobf (/ oobf substeps)
                             total-energy (/ total-energy substeps)
                             total-strain-energy (/ total-strain-energy substeps)
                             total-gpe (/ total-gpe substeps)
                             damage-inc (/ damage-inc substeps)
                             plastic-inc (/ plastic-inc substeps)
                             )
                       (push *t* *data-t*)
                       (push total-energy *data-e*)
                       (push total-strain-energy *data-se*)
                       (push total-gpe *data-gpe*)
                       (push damage-inc *data-di*)
                       (push plastic-inc *data-pi*)
                       ;; (format t "Energy estimate: ~E~%" energy-estimate)
                       (format t "Total energy: ~E~%" total-energy)
                       (setf last-e total-energy)
                       (defparameter *oobf* oobf)
                       (defparameter *energy* energy-estimate)
                       (defparameter *total-energy* total-energy)
                       (defparameter *total-strain-energy* total-strain-energy)
                       (defparameter *damage-inc* damage-inc)
                       (defparameter *plastic-inc* plastic-inc)
                       (when (>= steps settle-steps)
                         (setf (cl-mpm::sim-enable-damage *sim*) t)
                         (cl-mpm::iterate-over-mps
                          (cl-mpm:sim-mps *sim*)
                          (lambda (mp) (setf (cl-mpm/particle::mp-enable-plasticity mp) plasticity-enabled)))

                         (if (or
                              ;; t
                              ;; (> energy-estimate 1d-3)
                              ;; (> total-energy 1d2)
                              ;; (> total-energy 1d2)
                              ;; (> plastic-inc  1d-3)
                              (> total-energy 1d2)
                              ;; (< 0.5d0
                              ;;    (loop for mp across (cl-mpm:sim-mps *sim*)
                              ;;          maximizing (cl-mpm/particle::mp-damage mp)))
                              ;; nil
                              ;; t
                              )
                             (progn
                               (format t "Collapse timestep~%")
                               (setf
                                target-time collapse-target-time
                                (cl-mpm::sim-mass-scale *sim*) collapse-mass-scale
                                ;; (cl-mpm::sim-mass-scale *sim*) 1d0
                                ;; (cl-mpm:sim-damping-factor *sim*) 10d0
                                ))
                             (progn
                               (format t "Accelerated timestep~%")
                               (setf
                                ;; dt-scale 1.0d0
                                target-time 1d0
                                (cl-mpm::sim-mass-scale *sim*) 1d4
                                )
                               ;; (setf
                               ;;  target-time target-time-original
                               ;;  ;; (cl-mpm::sim-mass-scale *sim*) 1d3
                               ;;  (cl-mpm::sim-mass-scale *sim*) mass-scale
                               ;;  )
                               )))

                       (when (>= steps damp-steps)
                         (let ((ms (cl-mpm::sim-mass-scale *sim*)))
                           (setf (cl-mpm:sim-damping-factor *sim*)
                                 ;; (* 0d-3 ms)
                                 0.0d0
                                 )))
                       )
                     


                     (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       (format t "CFL dt estimate: ~f~%" dt-e)
                       (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf substeps substeps-e))
                     (setf (cl-mpm:sim-dt *sim*) (* dt-0 (sqrt (cl-mpm::sim-mass-scale *sim*))))
                     (setf substeps (floor target-time (cl-mpm:sim-dt *sim*)))
                     (format t "CFL dt estimate: ~f~%" dt)
                     (format t "CFL step count estimate: ~D~%" substeps)
                     ;; (setf (cl-mpm:sim-damping-factor *sim*)
                     ;;       (* (cl-mpm:sim-damping-factor *sim*) (expt 1d-3 1/40)))

                     (incf *sim-step*)
                     (plot *sim*)
                     (vgplot:title (format nil "Time:~F - KE ~E - SE ~E - dD ~E - dP ~E"  *t* *total-energy* *total-strain-energy* *damage-inc* *plastic-inc*))
                     (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                        :terminal "png size 1920,1080"
                                        )

                     (swank.live:update-swank)
                     ))))
  (values))


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

(let* ((s (cl-mpm/utils::stress-from-list (list 500d3 0d0 0d0 0d0 0d0 0d0)))
       (j2 (cl-mpm/constitutive::voigt-j2
            (cl-mpm/utils::deviatoric-voigt s)))
       (p (cl-mpm/constitutive::voight-trace s))
       (ft 300d3)
       (fc 500d3)
       (angle (atan (* 3 (/ (- fc ft) (+ fc ft)))))
       (k (* (/ 2d0 (sqrt 3)) (/ (* ft fc) (+ ft fc))))
       (s_1 (-
             (* (/ 3d0 (+ 3 (tan angle)))
                (+ (sqrt (* 3 j2)) (* 1/3 (tan angle) p)))
             k
             )))
  (format t "~A~%" s_1))

;; (defun test (s1 s2 s3)
;;   (let ((damage-inc-mat (cl-mpm/utils:matrix-zeros))
;;         (damage-tensor (cl-mpm/utils:voight-to-matrix (cl-mpm/utils:voigt-from-list (list 0.1d0 0.5d0 0d0 0d0 0d0 0d0))))
;;         (ybar-tensor (cl-mpm/utils:matrix-zeros))
;;         (E 1d9)
;;         (length 10d0)
;;         (Gf 1000d0)
;;         (init-stress 100d3)
;;         (cauchy-undamaged (cl-mpm/utils:voight-to-matrix (cl-mpm/utils:voigt-from-list (list s1 s2 s3 0d0 0d0 0d0)))))
;;     (multiple-value-bind (ls v) (cl-mpm/utils:eig cauchy-undamaged)
;;       (loop for i from 0 to 2
;;             do
;;                (let* ((sii (nth i ls))
;;                       (vii (magicl::column v i))
;;                       (vsi (magicl:@ vii (magicl:transpose vii)))
;;                       (dii (magicl::trace (magicl:@ damage-tensor vsi)))
;;                       ;; (dii 0d0)
;;                       (new-damage (cl-mpm/damage::damage-response-exponential sii E Gf length init-stress))
;;                       (damage-increment (- (max dii new-damage) dii))
;;                       )
;;                  ;; (when (> damage-increment 0d0)
;;                  ;;   (break))
;;                  (format t "Sii ~F : new-damage ~F ~%" sii new-damage)
;;                  (format t "Damage  dim ~D tensor ~A~%" i (magicl:@ damage-tensor vsi))
;;                  (format t "Damage  dim ~D vsi ~A~%" i vsi)
;;                  ;; (magicl:.+ ybar-tensor
;;                  ;;            (magicl:scale! vsi sii)
;;                  ;;            ybar-tensor)
;;                  (magicl:.+ damage-inc-mat
;;                             (magicl:scale! vsi damage-increment)
;;                             ;; (magicl:scale! vsi (* (/ dt tau) (- 1d0 (exp (- (* 1d0 (abs (- new-damage dii))))))))
;;                             damage-inc-mat))))
;;     (format t "Damage inc ~A~%" damage-inc-mat)
;;     (magicl:.+ damage-tensor
;;                damage-inc-mat
;;                damage-tensor)

;;     (format t "Damage tensor ~A~%" damage-tensor)
;;     )
;;   )



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




(defun setup ()
  (let* ((mesh-size 0.5d0)
         (mps-per-cell 2)
         (shelf-height 15.5)
         (soil-boundary 2)
         (shelf-aspect 1.0)
         (runout-aspect 2.0)
         (shelf-length (* shelf-height shelf-aspect))
         (domain-length (+ shelf-length (* runout-aspect shelf-height)))
         (shelf-height-true shelf-height)
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

    (when nil
      (let ((density 1.7d3)
            (g 9.8d0))
        (loop for mp across (cl-mpm:sim-mps *sim*)
              do
                 (let* ((depth (- (magicl:tref (cl-mpm/particle:mp-position mp) 1 0) shelf-height ))
                        (p (* density depth g 1/3))
                        (E (cl-mpm/particle::mp-e mp))
                        (nu (cl-mpm/particle::mp-e mp))
                        (K (/ E (* 2 (+ 1 nu) (- 1 nu nu))))
                        (p-e (/ p K))
                        (stress (cl-mpm/utils:voigt-from-list (list p p p 0d0 0d0 0d0)))
                        )
                   (setf
                    (cl-mpm/particle::mp-strain mp)
                    ;(cl-mpm/utils:voigt-from-list (list p-e p-e 0d0 0d0 0d0 0d0))
                    (magicl:linear-solve
                     (cl-mpm/particle::mp-elastic-matrix mp)
                     stress)
                    (cl-mpm/particle::mp-undamaged-stress mp)
                    stress
                    (cl-mpm/particle::mp-stress mp)
                    stress
                    )))
              ));)
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
                   (- shelf-length
                      15d0
                   ;   (* 0.30d0 shelf-height)
                      ))
                ;; (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                ;;    shelf-length)
                (< (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                   (+ soil-boundary 15d0))
                (> (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                   soil-boundary
                   )
                ;; (< (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                ;;    (- shelf-height 13d0)
                ;;    )
                )
             dir
             )))))
    (loop for mp across (cl-mpm::sim-mps *sim*)
          do (setf (cl-mpm/particle::mp-split-depth mp) 0))

    (let* ((sloped-height (- (- shelf-height soil-boundary) 7d0))
           (measured-angle 78d0)
           (undercut-angle ;(- 82.5d0 90d0)
             (- measured-angle 90d0)
                           )
           ;; (undercut-angle 0d0)
           (normal (magicl:from-list (list
                                      (cos (- (* pi (/ undercut-angle 180d0))))
                                      (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1)))
           (sloped-inflex-point
             (magicl:from-list (list (- shelf-length (* (tan (- (* pi (/ undercut-angle 180d0)))) sloped-height))
                                     soil-boundary)
                               '(2 1) :type 'double-float)))
      (cl-mpm/setup::remove-sdf *sim*
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
      (when nil
        (let ((cut-height (* 0.5d0 shelf-height-true))
              (cut-back-distance 0.15d0)
              (width ;(* 2d0 mesh-size)
                ;; 0.5d0
                (* 1.5d0 (cl-mpm/particle::mp-local-length (aref (cl-mpm:sim-mps *sim*) 0)))
                ))
          (cl-mpm/setup::apply-sdf *sim* (lambda (p) (cl-mpm/setup::line-sdf
                                                      (magicl:from-list (list (magicl:tref p 0 0)
                                                                              (magicl:tref p 1 0)) '(2 1))
                                                      (list (- (magicl:tref sloped-inflex-point 0 0)
                                                               (* cut-back-distance shelf-height-true))
                                                            shelf-height)
                                                      (list (- (magicl:tref sloped-inflex-point 0 0)
                                                               (* cut-back-distance shelf-height-true))
                                                            (- shelf-height cut-height))
                                        ;10d0
                                                      width
                                                      ))
                                   (lambda (mp v)
                                     (setf (cl-mpm/particle:mp-damage mp)
                                           (exp (- (expt (/ (+ width v) width) 4)))
                                           )
                                     ))
          )))

    (let* ((notched-depth 1d0)
           (undercut-angle 30d0)
           (normal (magicl:from-list (list
                                      (cos (- (* pi (/ undercut-angle 180d0))))
                                      (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1)))
           (sloped-inflex-point
             (magicl:from-list (list
                                (- shelf-length notched-depth)
                                soil-boundary)
                               '(2 1) :type 'double-float)

             )
           )

      (cl-mpm/setup:remove-sdf *sim*
                               (lambda (p)
                                 (if (and
                                      (> (magicl:tref p 1 0) soil-boundary))
                                     (cl-mpm/setup::plane-point-sdf
                                      (magicl:from-list (list (magicl:tref p 0 0)
                                                              (magicl:tref p 1 0)) '(2 1))
                                      normal
                                      sloped-inflex-point)
                                     1d0)
                                 ))
      )

   
     (let ((notch-height 0.0d0))
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
    ;;                  (max (cl-mpm/particle::mp-damage mp) (random 1.0d0)))))
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
    (defparameter *sim-step* 0)

  (defparameter *data-t* (list))
  (defparameter *data-e* (list))
  (defparameter *data-se* (list))
  (defparameter *data-gpe* (list))
  (defparameter *data-di* (list))
  (defparameter *data-pi* (list))
  )
(defparameter *data-t* (list))
(defparameter *data-e* (list))
(defparameter *data-se* (list))
(defparameter *data-gpe* (list))
(defparameter *data-di* (list))
(defparameter *data-pi* (list))


(defun setup-3d ()
  (let* ((mesh-size 0.25)
         (mps-per-cell 2)
         (shelf-height 15.0)
         (soil-boundary 2)
         (shelf-aspect 0.25)
         (runout-aspect 0.5)
         (shelf-length (* shelf-height shelf-aspect))
         (domain-length (+ shelf-length (* runout-aspect shelf-height)))
         (shelf-height-true shelf-height)
         (shelf-height (+ shelf-height soil-boundary))
         (depth (* 30d0))
         (offset (list 0 (* 0 mesh-size)
                       0
                       ))
         )
    (defparameter *sim*
      (setup-test-column (list domain-length
                               (+ shelf-height (* 5 mesh-size))
                               depth
                               )
                         (list domain-length shelf-height
                               depth
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
                   (- shelf-length
                      (* 0.2d0 shelf-height)))

                ;; (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                ;;    shelf-length)

                (< (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                   (+ soil-boundary 2d0))

                (> (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                   soil-boundary)

                ;; (< (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                ;;    (- shelf-height 13d0)
                ;;    )
                )
             dir
             )))))
    (loop for mp across (cl-mpm::sim-mps *sim*)
          do (setf (cl-mpm/particle::mp-split-depth mp) 0))

    (let* ((sloped-height (- (- shelf-height soil-boundary) 7d0))
           (measured-angle 60d0)
           (undercut-angle ;(- 82.5d0 90d0)
             (- measured-angle 90d0)
                           )
           ;; (undercut-angle 0d0)
           (normal (magicl:from-list (list
                                      (cos (- (* pi (/ undercut-angle 180d0))))
                                      (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1)))
           (sloped-inflex-point
             (magicl:from-list (list (- shelf-length (* (tan (- (* pi (/ undercut-angle 180d0)))) sloped-height))
                                     soil-boundary)
                               '(2 1) :type 'double-float)))
      (cl-mpm/setup::remove-sdf *sim*
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
      (let ((cut-height (* 0.0d0 shelf-height-true)))
        (cl-mpm/setup::damage-sdf
         *sim*
         (lambda (p)
           (if t ;(< (abs (- (magicl:tref p 2 0) (* 0.5d0 depth))) 20d0)
               (funcall
                (cl-mpm/setup::rectangle-sdf
                 (list (- (magicl:tref sloped-inflex-point 0 0)
                          (* 0.15d0 shelf-height-true))
                       shelf-height)
                 (list (* 1.0d0 mesh-size) cut-height)
                 ) p)
               0.99d0)
           ))))
    (let* ((notched-depth 10.0d0)
           (undercut-angle 30d0)
           (normal (magicl:from-list (list
                                      (cos (- (* pi (/ undercut-angle 180d0))))
                                      (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1)))
           (sloped-inflex-point
             (magicl:from-list (list
                                (- shelf-length notched-depth)
                                soil-boundary)
                               '(2 1) :type 'double-float)

             )
           )
      (cl-mpm/setup:remove-sdf *sim*
                               (lambda (p)
                                 (let* ((d (* 0.3d0 (abs (- 15d0 (magicl:tref p 2 0)))))
                                       (sloped-inflex-point
                                         (magicl:from-list (list
                                                            (- shelf-length
                                                               (* notched-depth
                                                                  (- 
                                                                     (/ 1d0 (+ 1d0 d))
                                                                     0d0
                                                                     )))
                                                            soil-boundary)
                                                           '(2 1) :type 'double-float)))
                                   (if (and
                                        (> (magicl:tref p 1 0) soil-boundary))
                                       (cl-mpm/setup::plane-point-sdf
                                        (magicl:from-list (list (magicl:tref p 0 0)
                                                                (magicl:tref p 1 0)) '(2 1))
                                        normal
                                        sloped-inflex-point)
                                       1d0))
                                 ))
      )

     (let ((notch-height 0.0d0))
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

(defun setup-under ()
  (let* ((mesh-size 0.5)
         (mps-per-cell 2)
         (shelf-height 15.0)
         (soil-boundary 2)
         (shelf-aspect 1.0)
         (runout-aspect 1.0)
         (shelf-length (* shelf-height shelf-aspect))
         (domain-length (+ shelf-length (* runout-aspect shelf-height)))
         (shelf-height-true shelf-height)
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
                   (- shelf-length
                      ;; 3d0
                      (* 0.5d0 shelf-height)
                      ))
                (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                   shelf-length)
                (> (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                   soil-boundary
                   )
                ;; (< (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                ;;    (- shelf-height 13d0)
                ;;    )
                )
             dir
             )))))
    (loop for mp across (cl-mpm::sim-mps *sim*)
          do (setf (cl-mpm/particle::mp-split-depth mp) 0))

    (let* ((sloped-height (- (- shelf-height soil-boundary) 7d0))
           (measured-angle 90d0)
           (undercut-angle ;(- 82.5d0 90d0)
             (- measured-angle 90d0)
                           )
           ;; (undercut-angle 0d0)
           (normal (magicl:from-list (list
                                      (cos (- (* pi (/ undercut-angle 180d0))))
                                      (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1)))
           (sloped-inflex-point
             (magicl:from-list (list (- shelf-length (* (tan (- (* pi (/ undercut-angle 180d0)))) sloped-height))
                                     soil-boundary)
                               '(2 1) :type 'double-float)))
      (cl-mpm/setup::remove-sdf *sim*
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
      (let ((cut-height (* 0.0d0 shelf-height-true)))
        (cl-mpm/setup::damage-sdf
         *sim*
         (lambda (p)
           (if t ;(< (abs (- (magicl:tref p 2 0) (* 0.5d0 depth))) 20d0)
               (funcall
                (cl-mpm/setup::rectangle-sdf
                 (list (- (magicl:tref sloped-inflex-point 0 0)
                          (* 0.15d0 shelf-height-true))
                       shelf-height)
                 (list (* 1.0d0 mesh-size) cut-height)
                 ) p)
               0.99d0)
           ))))
    ;; (let ((sand-layer 3d0))
    ;;   (cl-mpm/setup::apply-sdf
    ;;    *sim*
    ;;    (lambda (p)
    ;;      (if (and
    ;;           (> (magicl:tref p 1 0) (- soil-boundary sand-layer))
    ;;           (> (magicl:tref p 0 0) shelf-length))
    ;;          0d0
    ;;          1d0))
    ;;    (lambda (mp)
    ;;      (setf (cl-mpm/particle::mp-damage mp) 0.9d0)))
    ;;   )
    (let* ((notched-depth 0.5d0)
           (undercut-angle 40d0)
           (normal (magicl:from-list (list
                                      (cos (- (* pi (/ undercut-angle 180d0))))
                                      (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1)))
           (sloped-inflex-point
             (magicl:from-list (list
                                (- shelf-length notched-depth)
                                soil-boundary)
                               '(2 1) :type 'double-float)

             )
           )
      (cl-mpm/setup:remove-sdf *sim*
                               (lambda (p)
                                 (if (and
                                      (> (magicl:tref p 1 0) soil-boundary))
                                     (cl-mpm/setup::plane-point-sdf
                                      (magicl:from-list (list (magicl:tref p 0 0)
                                                              (magicl:tref p 1 0)) '(2 1))
                                      normal
                                      sloped-inflex-point)
                                     1d0)
                                 ))
      )

     (let ((notch-height 0.0d0))
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

(defun test-energy (strain k)
  (let* ((E 1d0)
         (de (cl-mpm/constitutive::linear-elastic-matrix E 0.1d0)))
    (let* ((strain+
             (multiple-value-bind (l v) (cl-mpm/utils::eig
                                         (cl-mpm/utils:voight-to-matrix strain))
               (loop for i from 0 to 2
                     do
                        (setf (nth i l) (max (nth i l) 0d0)))
               (cl-mpm/utils:matrix-to-voight (magicl:@ v
                                                        (magicl:from-diag l :type 'double-float)
                                                        (magicl:transpose v)))))
           (strain- (magicl:.- strain strain+))
           (e+ (sqrt (max 0d0 (* E (cl-mpm/fastmath::dot strain+ (magicl:@ de strain+))))))
           (e- (sqrt (max 0d0 (* E (cl-mpm/fastmath::dot strain- (magicl:@ de strain-))))))
           )
      (/
       (+ (* k e+) e-)
       (+ k 1d0))
      )))



(defun estimate-ductility-jirsek2004 (GF R ft E &optional (k 1d0))
 (let* ((e0 (/ ft E))
        (ef (+ (/ GF (* k R E e0)) (/ e0 2)))
        )
   (format t "Ductility aprox: ~F ~%" (- (* 2 (/ ef e0)) 1d0))))
