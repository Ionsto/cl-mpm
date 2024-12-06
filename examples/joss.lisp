(defpackage :cl-mpm/examples/joss
  (:use :cl))
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  0 0)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)

(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)

(setf *block-compile-default* t)

(in-package :cl-mpm/examples/joss)
;(declaim (optimize (debug 3) (safety 3) (speed 0)))
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(defun full-recompile ()
  (asdf:compile-system :cl-mpm/utils :force t)
  (asdf:compile-system :cl-mpm/fastmaths :force t)
  (asdf:compile-system :cl-mpm/forces :force t)
  (asdf:compile-system :cl-mpm/constitutive :force t)
  (asdf:compile-system :cl-mpm/mesh :force t)
  (asdf:compile-system :cl-mpm/particle :force t)
  (asdf:compile-system :cl-mpm/bc :force t)
  (asdf:compile-system :cl-mpm :force t)
  (ql:quickload :cl-mpm/examples/joss)
  )

(defun cl-mpm/damage::length-localisation (local-length local-length-damaged damage)
  ;; (+ (* local-length (- 1d0 damage)) (* local-length-damaged damage))
  (* local-length (max (sqrt (- 1d0 damage)) 1d-10))
  ;; local-length
  )

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-chalk-delayed) dt fbar)
  (cl-mpm::update-stress-kirchoff-damaged mesh mp dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-ugimp mesh mp dt fbar)
  ;; (cl-mpm::update-stress-kirchoff mesh mp dt fbar)
  ;; (cl-mpm::update-stress-linear mesh mp dt fbar)
  )

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-chalk-delayed-grassl) dt fbar)
  (cl-mpm::update-stress-kirchoff-damaged mesh mp dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-ugimp mesh mp dt fbar)
  ;; (cl-mpm::update-stress-kirchoff mesh mp dt fbar)
  ;; (cl-mpm::update-stress-linear mesh mp dt fbar)
  )
(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-mc) dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-damaged mesh mp dt fbar)
  (cl-mpm::update-stress-kirchoff-ugimp mesh mp dt fbar)
  ;; (cl-mpm::update-stress-kirchoff mesh mp dt fbar)
  ;; (cl-mpm::update-stress-linear mesh mp dt fbar)
  )

(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-chalk-delayed) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
                     (strain cl-mpm/particle::mp-strain)
                     (plastic-strain cl-mpm/particle::mp-strain-plastic)
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
        ;; (setf damage-increment (cl-mpm/damage::tensile-energy-norm
        ;;                         strain
        ;;                         ;; (cl-mpm/fastmaths:fast-.+
        ;;                         ;;  (magicl:scale plastic-strain (- 1d0 damage)))
        ;;                         E de))
        ;; (setf damage-increment (cl-mpm/damage::tensile-energy-norm strain E de))
        (setf damage-increment
              (max 0d0
                   (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile
                    stress
                    (* angle (/ pi 180d0)))))
        ;; (setf damage-increment
        ;;       (max 0d0
        ;;            (cl-mpm/damage::criterion-y
        ;;             stress
        ;;             strain
        ;;             damage
        ;;             E
        ;;             kc-r
        ;;             kt-r
        ;;             g-r)))

        ;; (setf damage-increment
        ;;       (max 0d0
        ;;            (cl-mpm/damage::criterion-dp-coheasion
        ;;             ;; cl-mpm/damage::criterion-dp-tensile
        ;;             ;; stress
        ;;             (cl-mpm/fastmaths:fast-scale-voigt stress (/ 1d0 (magicl:det def)))
        ;;             (* angle (/ pi 180d0)))))
        ;; (setf damage-increment
        ;;       (sqrt (* 3/2
        ;;                (cl-mpm/constitutive::voigt-j2
        ;;                 (cl-mpm/utils::deviatoric-voigt (magicl:scale stress (/ 1d0 (magicl:det def))))))))
        ;; (setf damage-increment (* E (cl-mpm/damage::modified-vm-criterion strain nu (/ fc ft))))

        ;; (multiple-value-bind (s_1 s_2 s_3) (cl-mpm/damage::principal-stresses-3d (magicl:scale stress (/ 1d0 (magicl:det def))))
        ;;   (setf damage-increment (max damage-increment
        ;;                               (max s_1
        ;;                                    s_2
        ;;                                    s_3
        ;;                                    0d0
        ;;                                    )))
        ;;   )
                                        ;(setf damage-increment (criterion-mc strain (* angle (/ pi 180d0)) E nu))
        ;; (setf damage-increment (criterion-dp strain (* angle (/ pi 180d0)) E nu))
        ;; (setf damage-increment (cl-mpm/damage::smooth-rankine-criterion (magicl:scale stress (/ 1d0 (magicl:det def)))))
        ;; (setf damage-increment
        ;;       (max 0d0
        ;;            (cl-mpm/damage::drucker-prager-criterion
        ;;             (magicl:scale stress (/ 1d0 (magicl:det def))) (* angle (/ pi 180d0)))))
        ;; (when (>= damage 1d0)
        ;;   (setf damage-increment 0d0))
        ;;Delocalisation switch
        (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
        (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
        ))))


(declaim (notinline plot))
(defun plot (sim)
  ;; (plot-collapse-curves)
  ;; (vgplot:plot *time-data* *damage-data*)
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :deformed
   ;; :colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :xy))
   :colour-func #'cl-mpm/particle::mp-damage
   ;; :colour-func #'cl-mpm/particle::mp-damage-shear
   ;; :colour-func #'cl-mpm/particle::mp-damage-ybar
   ;; :colour-func #'cl-mpm/particle::mp-strain-plastic-vm
   ;; :colour-func (lambda (mp)
   ;;                (let ((drive 
   ;;                        (*
   ;;                         (cl-mpm/particle::mp-damage-ybar mp)
   ;;                           ;; (- 1d0 (cl-mpm/particle::mp-damage mp))
   ;;                           )))
   ;;                  (if (> drive (cl-mpm/particle::mp-initiation-stress mp))
   ;;                      drive
   ;;                      0d0
   ;;                      )
   ;;                  ))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-yield-func mp))
   ;; :colour-func (lambda (mp) (magicl:tref (cl-mpm/particle::mp-stress mp) 0 0))
   ;; :colour-func (lambda (mp) (magicl:tref (cl-mpm/particle::mp-displacement mp) 1 0))
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
                                    1d0)))))

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
               ;; 'cl-mpm/damage::mpm-sim-usl-damage
               'cl-mpm/damage::mpm-sim-damage
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
             (angle 50d0)
             (init-c 26d3)
             (init-stress (cl-mpm/damage::mohr-coloumb-coheasion-to-tensile init-c (* angle (/ pi 180))))
             (downscale (/ 1d0 1d0))
             ;(gf (/ (expt (/ init-stress 6.88d0) 2) 1d9))
             (gf 5d0)
             ;; (gf 10d0)
             (length-scale (* h 3d0))
             ;; (length-scale 1.0d0)
             ;; (length-scale (/ (* 1d9 gf) (expt init-stress 2)))
             (ductility (estimate-ductility-jirsek2004 gf length-scale init-stress 1d9))
             ;; (ductility 10d0)
             )
        (format t "Estimated ductility ~E~%" ductility)
        (format t "Estimated lc ~E~%" length-scale)
        (format t "Estimated init stress ~E~%" init-stress)
        ;; (when (< ductility 1d0)
        ;;   (error "Ductility too low ~A" ductility))
        (cl-mpm:add-mps
         sim
         (cl-mpm/setup::make-mps-from-list
          (cl-mpm/setup::make-block-mps-list
           offset
           block-size
           (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
           density

           ;; 'cl-mpm/particle::particle-mc
           ;; ;; 'cl-mpm/particle::particle-dp
           ;; :E 1d9
           ;; :nu 0.24d0
           ;; :psi 0d0
           ;; ;; :psi (* 42d0 (/ pi 180))
           ;; :phi (* 42d0 (/ pi 180))
           ;; :c (* 131d3 1d0)
           ;; :phi-r (* 30d0 (/ pi 180))
           ;; :c-r 0d0
           ;; :softening 10d0
           'cl-mpm/particle::particle-chalk-delayed
           :E 1d9
           :nu 0.24d0

           :enable-plasticity t

           :ft 1d0
           :fc 10d0

           :friction-angle angle

           :kt-res-ratio 1d0
           :kc-res-ratio 0d0
           :g-res-ratio 0.5d0
           :peerlings-damage t

           :fracture-energy 3000d0
           :initiation-stress init-stress;18d3
           :delay-time 1d0
           :delay-exponent 3d0

           ;; :ductility 5d0
           :ductility ductility

           :critical-damage 1d0;(- 1.0d0 1d-3)
           :damage-domain-rate 0.9d0;This slider changes how GIMP update turns to uGIMP under damage

           :local-length length-scale
           :local-length-damaged 10d-10

           ;; :psi 0d0
           ;; :phi (* 42d0 (/ pi 180))
           ;; :c 131d3
           :psi (* 0d0 (/ pi 180))
           :phi (* angle (/ pi 180))
           :c (* init-c 1d2)

           :gravity -9.8d0
           :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
           ))))
      ;; (let ((mp-0 (aref (cl-mpm:sim-mps sim) 0)))
      ;;   (when (typep mp-0 'cl-mpm/particle::particle-chalk-brittle)
      ;;     (let* (
      ;;            (fc (cl-mpm/particle::mp-fc mp-0))
      ;;            (ft (cl-mpm/particle::mp-ft mp-0))
      ;;            (angle-d (* (/ 180 pi) (atan (* 3 (/ (- fc ft) (+ fc ft))))))
      ;;            (rc (cl-mpm/particle::mp-k-compressive-residual-ratio mp-0))
      ;;            (rs (cl-mpm/particle::mp-shear-residual-ratio mp-0))
      ;;            (angle-plastic (cl-mpm/particle::mp-phi mp-0))
      ;;            (angle-plastic-damaged (atan (* (/ rs rc) (tan angle-plastic))))
      ;;            )
      ;;       (format t "Estimated Gf ~F~%" (estimate-gf (cl-mpm/particle::mp-ductility mp-0)
      ;;                                                  (cl-mpm/particle::mp-initiation-stress mp-0)
      ;;                                                  (cl-mpm/particle::mp-local-length mp-0)))
      ;;       (format t "Chalk damage growth angle: ~F~%"
      ;;               angle-d)
      ;;       (format t "Chalk plastic virgin angle: ~F~%"
      ;;               (* (/ 180 pi) angle-plastic))
      ;;       (format t "Chalk plastic residual angle: ~F~%"
      ;;               (* (/ 180 pi) angle-plastic-damaged)))))

      ;; (cl-mpm/examples/tpb::calculate-ductility-param 1d9 47.5d0 0.05d0 20d3)

      (setf (cl-mpm:sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-enable-fbar sim) t)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (setf (cl-mpm::sim-mass-filter sim) 1d-10)
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
      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-varfix
             (cl-mpm:sim-mesh sim)
             '(0 nil nil)
             '(0 nil nil)
             '(0 0 nil)
             '(0 0 nil)
             '(nil nil 0)
             '(nil nil 0)))

      (defparameter *floor-bc*
        (cl-mpm/penalty::make-bc-penalty-point-normal
         sim
         (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
         (cl-mpm/utils:vector-from-list (list 0d0 (+ h-y) 0d0))
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
    ;; (let ((phi_0 (* 30d0;42d0
    ;;                 (/ pi 180)))
    ;;       (phi_1 (* 30d0 (/ pi 180)))
    ;;       (c_0
    ;;         0d0
    ;;         ;(* 131d3 0.25d0)
    ;;         )
    ;;       (soft 10000d0))
    ;;   ;; (setf c (+ rho_1
    ;;   ;;            (* (- rho_0 rho_1) (exp (- (* soft ps))))))
    ;;   (setf
    ;;    c (* c_0 (exp (- (* soft ps))))
    ;;    phi (+ phi_1 (* (- phi_0 phi_1) (exp (- (* soft ps))))))
    ;;   )
    ))

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
  ;; (vgplot:figure)
  (when (> (length *data-t*) 0)
    (vgplot:plot
     *data-t* (mapcar (lambda (x) (/ x (reduce #'max *data-e*))) *data-e*) "KE"
     *data-t* (mapcar (lambda (x) (/ x (reduce #'max *data-oobf*))) *data-oobf*) "oobf"
     *data-t* (mapcar (lambda (x) (/ x (reduce #'max *data-d*))) *data-d*) "Damage"
     ))
  )

;; (let ((refine (uiop:getenv "REFINE")))
;;   (when refine
;;     (setf *refine* (parse-integer (uiop:getenv "REFINE")))
;;     ))


(defun estimate-total-energy (sim)
  (loop for mp across (cl-mpm:sim-mps sim)
        sum (* (cl-mpm/particle::mp-mass mp)
                      (cl-mpm/fastmaths::mag-squared (cl-mpm/particle::mp-velocity mp)))))

(defgeneric estimate-energy-crit (sim))

(defmethod estimate-energy-crit ((sim cl-mpm::mpm-sim))
  ;; (loop for mp across (cl-mpm:sim-mps sim)
  ;;           sum (*
  ;;                0.5d0
  ;;                (cl-mpm/particle::mp-mass mp)
  ;;                  ;; (cl-mpm/particle::mp-damage mp)
  ;;                  (cl-mpm/fastmaths::mag-squared (cl-mpm/particle::mp-velocity mp))))
  (let ((energy 0d0)
        (mut (sb-thread:make-mutex)))
    (cl-mpm:iterate-over-nodes
     (cl-mpm:sim-mesh sim)
     (lambda (n)
       (when (cl-mpm/mesh:node-active n)
         (let ((e-inc (*
                       (cl-mpm/mesh::node-mass n)
                       (cl-mpm/fastmaths::mag (cl-mpm/mesh::node-velocity n))
                       )))
           (sb-thread:with-mutex (mut)
             (incf energy e-inc))))))
    energy)
  )
(defmethod estimate-energy-crit ((sim cl-mpm/mpi::mpm-sim-mpi))
  ;; (cl-mpm/mpi::mpi-sum
  ;;  (loop for mp across (cl-mpm:sim-mps sim)
  ;;        sum (* (cl-mpm/particle::mp-mass mp)
  ;;               (cl-mpm/fastmaths::mag-squared (cl-mpm/particle::mp-velocity mp)))))
  (cl-mpm/mpi::mpi-sum
   (let ((energy 0d0)
         (mut (sb-thread:make-mutex)))
     (cl-mpm:iterate-over-nodes
      (cl-mpm:sim-mesh sim)
      (lambda (n)
        (when (cl-mpm/mesh:node-active n)
          (when (cl-mpm/mpi::in-computational-domain sim (cl-mpm/mesh::node-position n))
            (let ((e-inc (*
                          (cl-mpm/mesh::node-mass n)
                          (cl-mpm/fastmaths::mag (cl-mpm/mesh::node-velocity n))
                          )))
              (sb-thread:with-mutex (mut)
                (incf energy e-inc)))))))
     energy)))


(defun estimate-energy-gradient (data h)
  (*
   (cond
       ((= (length data) 1)
        0d0)
       ((= (length data) 2)
        (/ (- (first data) (second data)) h))
       (t
        (/ (- (+ (first data) (third data)) (* 2 (second data))) (* h h))
        ))
   (cl-mpm:sim-mass-scale *sim*)
   ))

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
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk") *sim*)
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
  (let* ((target-time 1d0)
         (target-time-original target-time)
         (mass-scale (cl-mpm::sim-mass-scale *sim*))
         (accelerate-target-time 1d0)
         (accelerate-mass-scale 1d4)
         (collapse-target-time 0.1d0)
         (collapse-mass-scale 1d0)
         (plasticity-enabled (cl-mpm/particle::mp-enable-plasticity (aref (cl-mpm:sim-mps *sim*) 0)))
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 0.5d0)
         (settle-steps 0)
         (damp-steps 0)
         (sim-state :settle)
         (dt-0 0d0)
         (last-e 0d0)
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
         (damping-0
           (* 0d-2
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
      (format stream "steps,time,damage,plastic,energy,oobf,step-type~%"))

    (cl-mpm::update-sim *sim*)
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
      (format t "CFL dt estimate: ~f~%" dt-e)
      (format t "CFL step count estimate: ~D~%" substeps-e)
      (setf substeps substeps-e))

    (format t "Substeps ~D~%" substeps)
    (setf dt-0 (/ (cl-mpm:sim-dt *sim*) (sqrt (cl-mpm::sim-mass-scale *sim*))))
    (defparameter *data-damage* 0d0)
    (defparameter *data-plastic* 0d0)
    (defparameter *data-energy* 0d0)
    (defparameter *inst-data-oobf* 0d0)
    (setf *data-d* (list))


    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (setf (cl-mpm:sim-damping-factor *sim*) 
          (* 0.1d0 (cl-mpm/setup::estimate-critical-damping *sim*)))

    (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_conv_~5,'0d.vtk" 0)) *sim*)
    (cl-mpm/dynamic-relaxation:converge-quasi-static
     *sim*
     :dt-scale dt-scale
     :energy-crit 1d-2
     :oobf-crit 1d-2
     :substeps 10
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
       (cl-mpm/fastmaths::fast-zero (cl-mpm/particle:mp-velocity mp))))

    (setf (cl-mpm::sim-mass-scale *sim*) accelerate-mass-scale)
    (setf target-time accelerate-target-time)
    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (setf substeps (floor target-time (cl-mpm:sim-dt *sim*)))
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 1d-4
             (cl-mpm/setup::estimate-critical-damping *sim*)))

    (loop for vals in *to-damage-mps*
          do
          (destructuring-bind (mp k) vals
            (setf (cl-mpm/particle::mp-history-stress mp) k)
            (cl-mpm/damage::update-damage mp 1d-3)))

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
                       (format stream "~D,~f,~f,~f,~f,~f,~A~%"
                               steps
                               *t*
                               *data-damage*
                               *data-plastic*
                               *data-energy*
                               *inst-data-oobf*
                               sim-state))
                     (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     (let ((energy-estimate 0d0)
                           (oobf 0d0)
                           (total-energy 0d0)
                           (total-strain-energy 0d0)
                           (total-gpe 0d0)
                           (damage-inc 0d0)
                           (plastic-inc 0d0)
                           (work 0d0)
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


                       ;; (setf energy-estimate (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
                       (setf oobf (cl-mpm/dynamic-relaxation::estimate-oobf *sim*))
                       (setf
                        energy-estimate (/ energy-estimate substeps)
                        oobf (/ oobf substeps))
                       (when (> work 0d0)
                         (setf energy-estimate (abs (/ energy-estimate work))))
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
                         (setf *data-plastic damage-est))
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
                       (format t "Energy gradient: ~E~%" (estimate-energy-gradient *data-e* target-time))

                       (setf last-e total-energy)
                       (defparameter *oobf* oobf)
                       (defparameter *energy* energy-estimate)
                       (defparameter *total-energy* total-energy)
                       (defparameter *total-strain-energy* total-strain-energy)
                       (defparameter *damage-inc* damage-inc)
                       (defparameter *plastic-inc* plastic-inc)

                       (when (= steps damp-steps)
                          (setf sim-state :settle)
                          (setf (cl-mpm:sim-damping-factor *sim*)
                                0d0
                                ;; damping-0
                                ))
                       (when (= steps settle-steps)
                         (setf (cl-mpm::sim-enable-damage *sim*) t)
                         (cl-mpm::iterate-over-mps
                          (cl-mpm:sim-mps *sim*)
                          (lambda (mp) (setf (cl-mpm/particle::mp-enable-plasticity mp) plasticity-enabled))))
                        (when (>= steps settle-steps)
                          (if (or
                               ;; t
                               (> energy-estimate 1d-2)
                               (> oobf 1d-2)
                               ;; t
                               ;; nil
                               ;; (> work 1d6)
                               )
                              (when (not (eq sim-state :collapse))
                                (setf sim-state :collapse)
                                (format t "Changed to collapse~%"))
                              (progn
                                (when (not (eq sim-state :accelerate))
                                  (format t "Changed to accelerate~%")
                                  (setf sim-state :accelerate)
                                  (cl-mpm:iterate-over-mps
                                   (cl-mpm:sim-mps *sim*)
                                   (lambda (mp)
                                     (cl-mpm/fastmaths::fast-zero (cl-mpm/particle:mp-velocity mp)))))))
                          (case sim-state
                            (:accelerate
                             (format t "Accelerate timestep~%")
                             (setf
                              target-time accelerate-target-time
                              (cl-mpm::sim-mass-scale *sim*) accelerate-mass-scale))
                            (:collapse
                             (format t "Collapse timestep~%")
                             (setf
                              target-time collapse-target-time
                              (cl-mpm::sim-mass-scale *sim*) collapse-mass-scale))))
                       (format t "Sim state - ~A~%" sim-state)


                       (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                         (format t "CFL dt estimate: ~f~%" dt-e)
                         (format t "CFL step count estimate: ~D~%" substeps-e)
                         (setf substeps substeps-e))

                       (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
                       (setf substeps (floor target-time (cl-mpm:sim-dt *sim*)))
                       (format t "CFL dt estimate: ~f~%" dt)
                       (format t "CFL step count estimate: ~D~%" substeps)
                       ;; (setf (cl-mpm:sim-damping-factor *sim*)
                       ;;       (* (cl-mpm:sim-damping-factor *sim*) (expt 1d-3 1/40)))

                       (incf *sim-step*)
                       (plot *sim*)
                       (vgplot:title (format nil "Time:~F - KE ~E - OOBF ~E - Work ~E - ~A"  *t* energy-estimate *oobf* work sim-state)))
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


(defun est-angle ()
  (let* ((rc 0d0)
         (rs 0.5d0)
         (ratio (/ (- 1d0 rs) (- 1d0 rc))
                )
         (angle-plastic (* 50d0 (/ pi 180)))
         (angle-plastic-damaged (atan (* ratio (tan angle-plastic))))
         )
     (format t "Chalk plastic virgin angle: ~F~%"
             (* (/ 180 pi) angle-plastic))
     (format t "Chalk plastic residual angle: ~F~%"
             (* (/ 180 pi) angle-plastic-damaged))))

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
              (push mp-other merge-list)))))))))




(defun setup (&key (notch-length 1d0) (refine 1.0d0) (mps 2))
  (let* ((mesh-size (/ 1d0 refine))
         (mps-per-cell mps)
         (shelf-height 15.0)
         ;(soil-boundary (floor (* 15 1)))
         (soil-boundary 5)
         (shelf-aspect 1)
         (runout-aspect 1)
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

    ;; (cl-mpm/setup::initialise-stress-self-weight
    ;;  *sim*
    ;;  ;; shelf-height
    ;;  soil-boundary
    ;;  ;; :clipping-func
    ;;  ;; (lambda (pos)
    ;;  ;;   (< (cl-mpm/utils:varef pos 0) shelf-length))
    ;;  )
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
                      ))
                (> (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                   (- shelf-length 6d0))
                (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                   (+ shelf-length 2d0))
                (> (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                   (- soil-boundary 2d0)
                   ))
             dir)))))
    (loop for mp across (cl-mpm::sim-mps *sim*)
          do (setf (cl-mpm/particle::mp-split-depth mp) 0))

    (let* ((sloped-height (- (- shelf-height soil-boundary) 7d0))
           (measured-angle 78d0)
           (undercut-angle              ;(- 82.5d0 90d0)
             (- measured-angle 90d0))
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
                (* 2d0 (cl-mpm/particle::mp-local-length (aref (cl-mpm:sim-mps *sim*) 0)))
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
                                     (let ((d ;(* 0.99d0 (exp (- (expt (/ (+ width v) width) 1))))
                                             (* (- 1d0 1d-5) (cl-mpm/damage::weight-func (expt (+ v width) 2) width)))
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
      (cl-mpm/setup:remove-sdf
       *sim*
       (lambda (p)
         (if t
             (funcall
              (cl-mpm/setup::rectangle-sdf
               (list shelf-length (+ notch-height soil-boundary))
               (list (* 2d0 notch-height) notch-height)
               ) p)
             1d0))))

    (setf cl-mpm::*max-split-depth* 6))

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


(defun setup-3d ()
  (let* ((mesh-size 1.0)
         (mps-per-cell 3)
         (shelf-height 15.0)
         (soil-boundary 2)
         (shelf-aspect 1.0)
         (runout-aspect 1.0)
         (shelf-length (* shelf-height shelf-aspect))
         (domain-length (+ shelf-length (* runout-aspect shelf-height)))
         (shelf-height-true shelf-height)
         (shelf-height (+ shelf-height soil-boundary))
         (depth-aspect 1.00d0)
         (depth 5d0
           ;; (* depth-aspect shelf-height)
           )
         (offset (list 0 (* 0 mesh-size)
                       0
                       ))
         )
    (setf cl-mpm/damage::*enable-reflect-z* t)
    (defparameter *sim*
      (setup-test-column (list domain-length
                               (+ shelf-height (* 5 mesh-size))
                               depth
                               ;; mesh-size
                               )
                         (list domain-length shelf-height
                               depth
                               ;; mesh-size
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
           )))
      (let* ((notched-depth 2.0d0)
             (undercut-angle 20d0)
             (notched-width (* notched-depth))
             (d-mid 0d0
                                        ;(/ depth 2)
                    )
             (normal (magicl:from-list (list
                                        (cos (- (* pi (/ undercut-angle 180d0))))
                                        (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1)))

             )
        (cl-mpm/setup:remove-sdf *sim*
                                 (lambda (p)
                                   (if (and
                                        (> (magicl:tref p 1 0) soil-boundary)
                                        (> (magicl:tref p 2 0) (- d-mid notched-width))
                                        (< (magicl:tref p 2 0) (+ d-mid notched-width))
                                        )
                                       (let
                                           (
                                            (sloped-inflex-point
                                              (magicl:from-list (list
                                                                 (- shelf-length (* notched-depth
                                                                                    (- 1d0 (/ (abs (- (magicl:tref p 2 0) d-mid)) notched-width))
                                                                                    ))
                                                                 soil-boundary)
                                                                '(2 1) :type 'double-float))
                                            )
                                         (cl-mpm/setup::plane-point-sdf
                                          (magicl:from-list (list (magicl:tref p 0 0)
                                                                  (magicl:tref p 1 0)) '(2 1))
                                          normal
                                          sloped-inflex-point)
                                         )
                                       1d0)
                                   ))
  )

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

(defun setup-under (&key (mesh-size 0.5d0) (mps-per-cell 2) (height 15d0))
  (let* (
         ;; (mesh-size 0.5)
         ;; (mps-per-cell 2)
         (shelf-height ;15.0
           height
                       )
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
           (measured-angle 78d0)
           (undercut-angle ;(- 82.5d0 90d0)
             0d0
             ;; (- measured-angle 90d0)
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
    (let* ((notched-depth 0.0d0)
       (undercut-angle 45d0)
       (notched-width (* 0.5 notched-depth))
       (d-mid 0d0
              ;(/ depth 2)
              )
       (normal (magicl:from-list (list
                                  (cos (- (* pi (/ undercut-angle 180d0))))
                                  (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1)))

       )
      (cl-mpm/setup:remove-sdf *sim*
                               (lambda (p)
                                 (if (and
                                      (> (magicl:tref p 1 0) soil-boundary)
                                      (> (magicl:tref p 2 0) (- d-mid notched-width))
                                      (< (magicl:tref p 2 0) (+ d-mid notched-width))
                                      )
                                     (let
                                         (
                                          (sloped-inflex-point
                                            (magicl:from-list (list
                                                               (- shelf-length (* notched-depth
                                                                                  (- 1d0 (/ (abs (- (magicl:tref p 2 0) d-mid)) notched-width))
                                                                                  ))
                                                               soil-boundary)
                                                              '(2 1) :type 'double-float))
                                          )
                                       (cl-mpm/setup::plane-point-sdf
                                        (magicl:from-list (list (magicl:tref p 0 0)
                                                                (magicl:tref p 1 0)) '(2 1))
                                        normal
                                        sloped-inflex-point)
                                       )
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
           (e+ (sqrt (max 0d0 (* E (cl-mpm/fastmaths::dot strain+ (magicl:@ de strain+))))))
           (e- (sqrt (max 0d0 (* E (cl-mpm/fastmaths::dot strain- (magicl:@ de strain-))))))
           )
      (/
       (+ (* k e+) e-)
       (+ k 1d0))
      )))



(defun estimate-ductility-jirsek2004 (GF R ft E &optional (k 1d0))
 (let* ((e0 (/ ft E))
        (ef (+ (/ GF (* k R E e0)) (/ e0 2)))
        )
   (- (* 2d0 (/ ef e0)) 1d0)))



;; (let* ((mpi-count 16)
;;        (x (loop for x from 0 upto 1 by (/ 1d0 mpi-count) collect x))
;;        (z (loop for i from 0 upto 1 by (/ 1d0 mpi-count) collect 1d0))
;;        (density-ratio 0.5)
;;        (inflection (/ 0.4 density-ratio))
;;       (di (* density-ratio inflection)))
;;   (let ((y (mapcar
;;             (lambda (p)
;;                                         (expt p 1.6)
;;               ;; (if (> p inflection)
;;               ;;     (+
;;               ;;      (* (/ (- 1d0 di)
;;               ;;            (- 1d0 inflection)) p)
;;               ;;      (-
;;               ;;       di
;;               ;;       (*
;;               ;;        (/ (- 1d0 di)
;;               ;;           (- 1d0 inflection))
;;               ;;        inflection)
;;               ;;       ))
;;               ;;     (* density-ratio p)
;;               ;;     )
;;               ) x))
;;         (inf (mapcar
;;             (lambda (p)
;;               ;; (expt p 1.6)
;;               (if (> p inflection)
;;                   (+
;;                    (* (/ (- 1d0 di)
;;                          (- 1d0 inflection)) p)
;;                    (-
;;                     di
;;                     (*
;;                      (/ (- 1d0 di)
;;                         (- 1d0 inflection))
;;                      inflection)
;;                     ))
;;                   (* density-ratio p)
;;                   )
;;               ) x))
;;         )
;;     ;; (vgplot:figure)
;;     (vgplot:plot y (mapcar (lambda (i) (* 0.25d0 i)) z) ";;with points"
;;                  inf (mapcar (lambda (i) (* 0.5d0 i)) z) ";;with points"
;;                  x (mapcar (lambda (i) (* 0.75d0 i)) z) ";;with points"
;;                  )
;;     ;; (vgplot:axis '())
;;     (vgplot:axis (list 0d0 1d0 0d0 1d0))
;;     ))



;; (defun test ()
;;   (let ((iters 100000000)
;;         (eps (magicl:zeros '(3 1)))
;;         (res (magicl:zeros '(3 1))))
;;     (declare (magicl:matrix/double-float eps res))
;;     ;; (magicl:.+ eps eps res)
;;     ;; (cl-mpm/fastmaths::fast-add res eps)
;;     ;; (magicl.simd::simd-add (magicl::matrix/double-float-storage eps) 
;;     ;;           (magicl::matrix/double-float-storage eps)
;;     ;;           (magicl::matrix/double-float-storage res))
;;     ;; (magicl.simd::.+-unsafe eps eps res)

;;     (format t "Test magicl~%")
;;     (time (lparallel:pdotimes (i iters)
;;             (magicl:.+ eps eps res)))
;;     (format t "Test magicl.simd~%")
;;     (time (lparallel:pdotimes (i iters)
;;             (magicl.simd::.+-simd eps eps res)))
;;     (format t "Test magicl.simd~%")
;;     (time (lparallel:pdotimes (i iters)
;;             (cl-mpm/fastmaths:fast-add res eps)))))
;; (setup)
;; (let* ((eps (cl-mpm/constitutive::swizzle-coombs->voigt
;;               (cl-mpm/utils:voigt-from-list (list 10d0 0d0 5d0 1d0 0d0 0d0))))
;;         (E 1d0)
;;         (nu 0.0440d0)
;;         (angle 0.10d0)
;;         (c 0.1000d0)
;;         (de (cl-mpm/constitutive::linear-elastic-matrix E nu))
;;        (iters 100000))

;;   (time (lparallel:pdotimes (i iters)
;;           (declare (ignore i))
;;           (let ((sig (cl-mpm/constitutive::linear-elastic-mat eps de)))
;;             (cl-mpm/constitutive::mc-plastic sig de eps
;;                                              E
;;                                              nu
;;                                              angle
;;                                              0d0
;;                                              c)
;;             )
;;           ))
;;   )
 


(defun test-add () 
  (let ((a (cl-mpm/utils::stretch-dsvp-3d-zeros))
        (b (cl-mpm/utils::stretch-dsvp-3d-zeros))
        (res (cl-mpm/utils::stretch-dsvp-3d-zeros))
        (iter 100000000))
    ;; (time
    ;;  (dotimes (i iter)
    ;;    (magicl:.+ a b res)))
    (time
     (dotimes (i iter)
       (cl-mpm/fastmaths::fast-.+ a b res)))
    (time
     (dotimes (i iter)
       (cl-mpm/fastmaths::fast-.+-stretch a b res)))
    ;; (time
    ;;  (lparallel:pdotimes (i iter)
    ;;    (magicl:scale! a 0d0)))
    ;; (time
    ;;  (lparallel:pdotimes (i iter)
    ;;    (cl-mpm/fastmaths::fast-scale! a 0d0)))
    ;; (time
    ;;  (lparallel:pdotimes (i iter)
    ;;    (cl-mpm/fastmaths::fast-zero a)))
    ))


;; (defun test-strain-rate ()
;;   (setup)
;;   (let* ((mesh (cl-mpm:sim-mesh *sim*))
;;          (mps (cl-mpm:sim-mps *sim*))
;;          (mp (aref mps 0))
;;          (dt (cl-mpm:sim-dt *sim*))
;;         )
;;     (cl-mpm:iterate-over-nodes
;;      mesh
;;      (lambda (n) (setf (cl-mpm/mesh:node-active n) t)))
;;     (time
;;      (dotimes (i 100000)
;;        (cl-mpm::calculate-strain-rate mesh mp dt))))
;;   )
;; (defun test-stretch ()
;;   (let ((res (cl-mpm/utils::stretch-dsvp-3d-zeros)))
;;     (time
;;      (dotimes (i 1000000)
;;        (cl-mpm/shape-function::assemble-dstretch-3d-prealloc (list 1d0 2d0 3d0) res)))
;;     (time
;;      (dotimes (i 1000000)
;;        (cl-mpm/shape-function::assemble-dstretch-3d-prealloc-fast (list 1d0 2d0 3d0) res)))
;;     )
;;   )


(defun shit-mult (ma mb mc)
  (declare (magicl:matrix/double-float ma mb mc)
           (optimize (speed 3) (safety 0) (debug 0)))
  (let ((a (magicl::matrix/double-float-storage ma))
        (b (magicl::matrix/double-float-storage mb))
        (c (magicl::matrix/double-float-storage mc))
        )
    (declare ((simple-array double-float (6)) a)
             ((simple-array double-float (3)) b c))
    (flet ((tref (m x y)
             (aref m (+ x (* 3 y)))))
      (loop for i fixnum from 0 below 3
            do
               (setf
                (aref c i)
                (+
                 (* (aref b 0) (tref a i 0))
                 (* (aref b 1) (tref a i 1))
                 (* (aref b 2) (tref a i 2))
                 )
                )
            )))
  mc)

;; (let ((a (cl-mpm/utils::matrix-from-list (list
;;                                           1d0 5d0 8d0
;;                                           5d0 2d0 3d0
;;                                           9d0 1.2d0 3d0)))
;;       (b (cl-mpm/utils:vector-from-list (list 2d0 3d0 4d0)))
;;       (res (cl-mpm/utils:vector-zeros))
;;       (iters 10000000)
;;       )
;;   (magicl:mult a b :target res)
;;   (pprint res)
;;   (setf res (cl-mpm/utils:vector-zeros))
;;   (cl-mpm/fastmaths::@-m-v a b res)
;;   (pprint res)
;;   ;; (time
;;   ;;  (dotimes (i iters)
;;   ;;    (magicl:mult a b :target res)))
;;   ;; (time
;;   ;;  (dotimes (i iters)
;;   ;;    (shit-mult a b res)))
;;   )

;; (let* ((grads (list 1d0 2d0 -1d0))
;;        (stretch (cl-mpm/shape-function::assemble-dstretch-3d grads))
;;       (vel (cl-mpm/utils::vector-from-list (list 0.1d0 1d0 2d0)))
;;       (res (cl-mpm/utils::stretch-dsvp-voigt-zeros))
;;       (iters 10000000)
;;       )
;;   ;; (magicl:mult stretch vel :target res)
;;   (time
;;    (dotimes (i iters)
;;      (let ((res (cl-mpm/utils::stretch-dsvp-voigt-zeros))
;;            (st (cl-mpm/utils::matrix-zeros))
;;            (temp (cl-mpm/utils::matrix-zeros))
;;            )
;;        (cl-mpm/fastmaths::@-stretch-vec stretch vel res)
;;        (cl-mpm/utils::voight-to-stretch-prealloc res temp)
;;        (cl-mpm/fastmaths::fast-.+-matrix st temp st)
;;        )))
;;   (time
;;    (dotimes (i iters)
;;     (let ((st (cl-mpm/utils::matrix-zeros)))
;;       (cl-mpm/shape-function::@-combi-assemble-dstretch-3d grads vel st))))
;;   ;; (setf res (cl-mpm/utils::stretch-dsvp-voigt-zeros))
;;   ;; (cl-mpm/fastmaths::@-stretch-vec stretch vel res)
;;   ;; (pprint res)
;;   ;; (setf res (cl-mpm/utils::stretch-dsvp-voigt-zeros))
;;   ;; (cl-mpm/fastmaths::@-stretch-vec-simd stretch vel res)
;;   ;; (pprint res)
;;   ;; (time
;;   ;;  (dotimes (i iters)
;;   ;;    (magicl:mult stretch vel :target res)))
;;   ;; (time
;;   ;;  (dotimes (i iters)
;;   ;;    (cl-mpm/fastmaths::@-stretch-vec stretch vel res)))
;;   ;; (time
;;   ;;  (dotimes (i iters)
;;   ;;    (cl-mpm/fastmaths::@-stretch-vec-simd stretch vel res)))
;;   )


(defun setup-dsvp ()
  (let ((dsvp (cl-mpm/shape-function::assemble-dsvp-3d (list 1d0 2d0 1d0)))
        (stress (cl-mpm/utils:voigt-from-list (list 1d0 2d0 1d0 2d0 5d0 8d0)))
        (res (cl-mpm/utils:vector-zeros))
        (grads (list 1d0 2d0 1d0))
        (iters 1000000))
    ;; (pprint (magicl:@ (magicl:transpose dsvp) stress))
    ;; (pprint (progn (cl-mpm/forces::mult-force dsvp stress 1d0 res) res))
    ;; (setf res (cl-mpm/utils:vector-zeros))
    ;; (pprint (progn (cl-mpm/forces::@-dsvp-vec dsvp stress 1d0 res) res))
    ;; (setf res (cl-mpm/utils:vector-zeros))
    ;; (pprint (progn (cl-mpm/forces::@-dsvp-vec-simd dsvp stress 1d0 res) res))
    ;; (time
    ;;  (dotimes (i iters)
    ;;    (magicl:@ (magicl:transpose dsvp) stress)))
    ;; (time
    ;;  (dotimes (i iters)
    ;;    (cl-mpm/forces::mult-force dsvp stress -1d0 res)))
    (time
     (dotimes (i iters)
       (cl-mpm/shape-function::assemble-dsvp-3d-prealloc grads dsvp)
       (cl-mpm/forces::@-dsvp-vec-simd dsvp stress -1d0 res)))
    (time
     (dotimes (i iters)
       (cl-mpm/shape-function::@-assemble-dsvp-3d-prealloc grads stress res)))
    )
  )


;; (let* ((ddF (cl-mpm/utils::voight-to-stretch (magicl:@ (cl-mpm/shape-function::assemble-dstretch-3d (list 1d0 2d0 0d0)) (cl-mpm/utils:vector-from-list (list 5d0 7d0 1d0)))))
;;        (df (cl-mpm/utils::matrix-from-list '(1d0 0d0 0d0
;;                                              0d0 1d0 0d0
;;                                              0d0 0d0 1d0)))
;;        (strain (cl-mpm/utils:voigt-zeros))
;;       )

;;   (cl-mpm/fastmaths::fast-.+-matrix df ddF df)
;;   (cl-mpm/ext:kirchoff-update strain df)
;;   (pprint (cl-mpm/constitutive::swizzle-voigt->coombs strain))
;;   )


(defun test-simd-3d ()
  (let ((mp (aref (cl-mpm:sim-mps *sim*) 0))
        (mesh (cl-mpm:sim-mesh *sim*))
        (iter 100000))
    (time
     (dotimes (i iter)
       (cl-mpm::iterate-over-neighbours-shape-gimp-3d
        mesh
        mp
        (lambda (mesh mp node weight grads &rest args)
          (* weight (nth 0 grads))
          ;(pprint weight)
          ;(pprint grads)
          ))))
    (time
     (dotimes (i iter)
       (cl-mpm::iterate-over-neighbours-shape-gimp-simd-3d
        mesh
        mp
        (lambda (mesh mp node weight grads &rest args)
          (* weight (nth 0 grads))
          ;(pprint weight)
          ;(pprint grads)
          )
        )))
    (cl-mpm::iterate-over-neighbours-shape-gimp-3d
     mesh
     mp
     (lambda (mesh mp node weight grads &rest args)
       (pprint weight)
       (pprint grads)
       ))
    (cl-mpm::iterate-over-neighbours-shape-gimp-simd-3d
     mesh
     mp
     (lambda (mesh mp node weight grads &rest args)
       (pprint weight)
       (pprint grads)
       )
     )
    )
  )


;; (let ((iter 10000000))
;;    (let ((def (cl-mpm/utils:matrix-zeros))
;;          (df (cl-mpm/utils:matrix-zeros)))
;;      (time
;;       (lparallel:pdotimes (i iter)
;;         (setf def (magicl:@ df def))
;;         )))
;;   (let ((def (cl-mpm/utils:matrix-zeros))
;;         (df (cl-mpm/utils:matrix-zeros)))
;;     (time
;;      (lparallel:pdotimes (i iter)
;;        (let ((temp (cl-mpm/utils:matrix-zeros)))
;;          (magicl:mult df def :target temp)
;;          (cl-mpm/utils:matrix-copy-into temp def))
;;        )))
;;   )

;; (let* ((length-scale 0.25d0)
;;        (init-stress 1.5d6)
;;        (gf (/ (expt (/ init-stress 6.88d0) 2) 1d9))
;;        (ductility (estimate-ductility-jirsek2004 gf length-scale init-stress 1d9)))
;;   (format t "Init-stress ~E~%" init-stress)
;;   (format t "Gf ~E~%" gf)
;;   (format t "Ductility ~E~%" ductility))

(defun test-mc-dp ()
  (dotimes (i 100)
    (let* ((p 1d0)
           (eps
             (cl-mpm/utils:voigt-from-list (loop for i from 0 to 5
                                                    collect (- (random 2d0) 1d0)))
             ;; (cl-mpm/utils:voigt-from-list (list 10d0 0d0 0d0 0d0 0d0 0d0))
             )
           (E 1d0)
           (nu (random 0.45d0))
           (angle (* (random 80d0) (/ pi 180d0)))
           (c (random 1.0d0))
           (de (cl-mpm/constitutive::linear-elastic-matrix E nu))
           (sig-i (magicl:@ de eps)))
      (format t "Iter ~D~%" i)
      ;; (pprint eps)
      (multiple-value-bind (sig eps-e f) (cl-mpm/constitutive::mc-plastic sig-i de eps E nu angle angle c)
        (cl-mpm/constitutive::swizzle-voigt->coombs sig))
      )))


;; (let* (
;;        (init-stress 60d3)
;;        (gf (* 3d0))
;;        (length-scale 0.25d0)
;;        (ductility (estimate-ductility-jirsek2004 gf length-scale init-stress 1d9)))
;;   (format t "~%Estimated ductility ~E~%" ductility)
;;   (format t "~%Estimated Gf ~F~%" gf)
;;   (format t "~%Estimated lc ~E~%" length-scale)
;;   (format t "~%Estimated init stress ~F kPa~%" (* 1d-3 init-stress)))


;; (let* ((E 1d9)
;;        (nu 0.24d0)
;;        (K (/ E (* 3 (- 1d0 (* 2 nu)))))
;;        (G (/ E (* 2 (+ 1d0 nu))))
;;        )
;;   (format t "~%E: ~E nu: ~E~%" E nu)
;;   (format t "~%K: ~E G: ~E~%" K G)
;;   (let* ((kr 1d-2)
;;          (gr 1d-3)
;;          (K (* K kr))
;;          (G (* G gr))
;;          (Ed (/ (* 9 K G)
;;                 (+ (* 3 K) G)))
;;          (nud (/ (- (* 3 K) (* 2 G))
;;                  (+ (* 6 K) (* 2 G)))))
;;     (format t "~%Kd: ~E Gd: ~E~%" K G)
;;     (format t "~%Ed: ~E nud: ~E~%" Ed nud)))


;;1-0.5-6.47
;;0.5-0.5-7.41
;;0.25-0.5-7.41

;; (defmethod cl-mpm::update-node-forces ((sim cl-mpm::mpm-sim))
;;   (with-accessors ((damping cl-mpm::sim-damping-factor)
;;                    (mass-scale cl-mpm::sim-mass-scale)
;;                    (mesh cl-mpm::sim-mesh)
;;                    (dt cl-mpm::sim-dt))
;;       sim
;;     (cl-mpm::iterate-over-nodes
;;      mesh
;;      (lambda (node)
;;        (cl-mpm::calculate-forces-cundall node damping dt mass-scale)))))
;; (setf *run-sim* nil)

;;h=1
;;0 - 1.83379480274271d+4
;;1 - 3.694570161832935d+4
;;h=0.5
;;0 - 1.83379480274271d+4
;;1 - 3.


(defun make-plot-data ()
  )

(defun estimate-max-stress-good (&key (refine 1)
                                   (mps 2)
                                   (notch-length 0d0))
  (vgplot:close-all-plots)
  (let ((dt-scale 0.50d0))
    (setup :notch-length notch-length
           :refine refine
           :mps mps)
    (cl-mpm/output::save-simulation-parameters
     #p"output/settings.json"
     *sim*
     (list :dt 1d0))
    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (setf (cl-mpm/particle::mp-initiation-stress mp) 1d10))

    (setf (cl-mpm:sim-mass-scale *sim*) 1d0)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 0.1d0 (cl-mpm/setup:estimate-critical-damping *sim*)))
    (setf (cl-mpm:sim-dt *sim*) (* dt-scale (cl-mpm/setup:estimate-elastic-dt *sim*)))

    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (setf
              (cl-mpm/particle::mp-initiation-stress mp) 1d10
              (cl-mpm/particle::mp-enable-plasticity mp) nil))

    (defparameter *damage-data* (list))
    (defparameter *energy-data* (list))
    (defparameter *time-data* (list))

    (cl-mpm/output::save-simulation-parameters
     #p"output/settings.json"
     *sim*
     (list :dt 1d0))
    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (setf
              (cl-mpm/particle::mp-initiation-stress mp) 1d10
              (cl-mpm/particle::mp-enable-plasticity mp) nil
              ))
    (setf (cl-mpm::sim-enable-damage *sim*) nil)
    (let ((step 0))
      (cl-mpm/dynamic-relaxation:converge-quasi-static
       *sim*
       :dt-scale dt-scale
       :energy-crit 1d-2
       :oobf-crit 1d-2
       :substeps 50
       :conv-steps 200
       :dt-scale dt-scale
       :post-iter-step
       (lambda (i e o)
         (cl-mpm/damage::calculate-damage *sim*)
         (setf step i)
         (plot *sim*)
         (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" i)) *sim*)
         ;; (let ((dy (lparallel:pmap-reduce (lambda (mp)
         ;;                                  (cl-mpm/particle::mp-damage-ybar mp))
         ;;                                #'max (cl-mpm:sim-mps *sim*))))
         ;;              (format t "Max damage-ybar ~E~%" dy)
         ;;              (push dy *damage-data*)
         ;;              (push i *time-data*))
         ;; (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" i))
         ;;                    :terminal "png size 1920,1080")
         ))
      (cl-mpm/damage::calculate-damage *sim*)
      (plot *sim*)
      (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" (1+ step))) *sim*))
    (let ((dy (lparallel:pmap-reduce (lambda (mp)
                                       (cl-mpm/particle::mp-damage-ybar mp))
                                     #'max (cl-mpm:sim-mps *sim*))))
      (format t "Converged max damage-ybar ~E~%" dy)
      ;; (push dy *damage-data*)
      ;; (push i *time-data*)
      dy)))

(ql:quickload :str)
(defun sweep-stress ()
  (setf *run-sim* t)
  (uiop:ensure-all-directories-exist (list "./max-stress/"))
  ;; (lisp-stat:defdf *sweep-df* (lisp-stat:make-df '(:refine :mps :length :max-stress)
  ;;                                             '(#() #() #() #())))
  ;; (lisp-stat:matrix-df (lisp-stat:keys *sweep-df*)
  ;;                      (lisp-stat:stack-rows
  ;;                       (lisp-stat:rows *sweep-df*))
  ;;                      #(0 0 0 0))
  (defparameter *sweep-data* (list))
  (let (;(mps 2)
        (name "LS1")
        (notch-length 0d0))
    (dolist (mps '(2 3 4))
      (loop for refine in (list 1 2 3 4)
            while *run-sim*
            do
               (let ((ybar
                       (estimate-max-stress-good :refine refine
                                                 :mps mps
                                                 :notch-length notch-length
                                                 )))
                 (format t "Converged REFINE: ~F - MPS: ~D - YBAR: ~E~%" refine mps ybar)
                 (str:to-file
                  (merge-pathnames (format nil "./max-stress/data_~A_~D_~D_~F.json" name refine mps notch-length))
                  (jonathan:to-json
                   (list :refine refine :mps mps :notch-length notch-length :max-stress ybar :name name
                         :length-scale (cl-mpm/particle::mp-local-length (cl-mpm::get-mp *sim* 0)))))
                 )
               (cl-mpm/output:save-vtk-mesh (merge-pathnames (format nil "./max-stress/mesh_~D.vtk" refine)) *sim*)
               (cl-mpm/output:save-vtk (merge-pathnames (format nil "./max-stress/final_~A_~D_~D_~F.vtk" name refine mps notch-length)) *sim*)))))
;; (defmethod cl-mpm::update-node-forces ((sim cl-mpm::mpm-sim))

;;   (with-accessors ((damping cl-mpm::sim-damping-factor)
;;                    (mass-scale cl-mpm::sim-mass-scale)
;;                    (mesh cl-mpm::sim-mesh)
;;                    (dt cl-mpm::sim-dt))
;;       sim
;;     (cl-mpm::iterate-over-nodes
;;      mesh
;;      (lambda (node)
;;        (cl-mpm::calculate-forces-cundall node damping dt mass-scale)))))



(defun stop ()
  (setf *run-sim* nil)
  (setf cl-mpm/dynamic-relaxation::*run-convergance* nil))


(defun plot-mpi-distr ()
  (let* ((density-ratio 0.50d0)
         (inflection (/ 0.25d0 density-ratio))
         (di (* density-ratio inflection))
         ;; (density-ratio (/ di inflection))
         )
    (let* ((x (loop for i from 0 upto 1d0 by 0.01d0 collect i))
           (y (mapcar
               (lambda (p)
                 (float
                  (if (> p inflection)
                      (+
                       (* (/ (- 1d0 di)
                             (- 1d0 inflection)) p)
                       (-
                        di
                        (*
                         (/ (- 1d0 di)
                            (- 1d0 inflection))
                         inflection)
                        ))
                      (* density-ratio p)
                      )
                  0d0
                  )) x)))
      ;; (vgplot:figure)
      (vgplot:plot x y ";;x")
      )
    ))
