(defpackage :cl-mpm/examples/joss
  (:use :cl))
(defclass cl-mpm/particle::particle-chalk-erodable (cl-mpm/particle::particle-chalk-delayed
                                                     cl-mpm/particle::particle-erosion)
  ())

(defmethod cl-mpm/erosion::mp-erosion-enhancment ((mp cl-mpm/particle::particle-chalk-erodable))
  (+ 1d0 (* 10 (cl-mpm/particle::mp-damage mp)))
  ;(setf weathering (* weathering (+ 1d0 (* 10 (cl-mpm/particle::mp-strain-plastic-vm mp)))))
  )
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
(defclass bc-water-damage (cl-mpm/buoyancy::bc-scalar)
  ((damage-rate
    :initform 1d0
    :initarg :damage-rate
    :accessor bc-water-damage-damage-rate)))

;; (defclass bc-erode (cl-mpm/buoyancy::bc-scalar)
;;   ((damage-rate
;;     :initform 1d0
;;     :initarg :damage-rate
;;     :accessor bc-water-damage-damage-rate)))

;; (defun make-bc-erode (sim datum rate lower-datum upper-datum)
;;   (with-accessors ((mesh cl-mpm:sim-mesh))
;;       sim
;;     (with-accessors ((h cl-mpm/mesh:mesh-resolution))
;;         mesh
;;       (make-instance 'bc-erode
;;                      :index nil
;;                      :damage-rate rate
;;                      :damage-volume nil
;;                      :scalar-func (lambda (pos) (expt (- 1d0 (/ (abs (- (cl-mpm/utils:varef pos 1) datum))
;;                                                                 (- upper-datum datum)
;;                                                                 ))
;;                                                       12))
;;                      :clip-func (lambda (pos) (and (>= (cl-mpm/utils:varef pos 1) lower-datum)
;;                                                    (<= (cl-mpm/utils:varef pos 1) upper-datum)
;;                                                    ))
;;                      :sim sim))))
;; (defmethod cl-mpm/bc::apply-bc ((bc bc-erode) node mesh dt)
;;   (call-next-method)
;;   ;; (break)
;;   (with-accessors ((sim cl-mpm/buoyancy::bc-buoyancy-sim))
;;       bc
;;     (when (cl-mpm::sim-enable-damage sim)
;;       (loop for mp across (cl-mpm:sim-mps sim)
;;             do
;;                (let ((weathering 0d0))
;;                  (cl-mpm:iterate-over-neighbours
;;                   mesh
;;                   mp
;;                   (lambda (mesh mp node svp grads fsvp fgrads)
;;                     (incf weathering (* svp (cl-mpm/mesh::node-boundary-scalar node)))))
;;                  (setf weathering (* weathering (+ 1d0 (* 8 (cl-mpm/particle:mp-damage mp)))))
;;                  (setf (cl-mpm/particle::mp-boundary mp) weathering)
;;                  (setf weathering (min weathering 0d0))
;;                  (let ((density (/ (cl-mpm/particle::mp-mass mp) (cl-mpm/particle::mp-volume mp))))
;;                    (setf
;;                     (cl-mpm/particle::mp-volume mp)
;;                     (max
;;                      0d0
;;                      (-
;;                       (cl-mpm/particle::mp-volume mp)
;;                       (abs (*
;;                             (bc-water-damage-damage-rate bc)
;;                             weathering dt)))))
;;                    (setf (cl-mpm/particle::mp-mass mp) (* density (cl-mpm/particle::mp-volume mp))))))
;;       (cl-mpm::remove-mps-func sim (lambda (mp) (= 0d0 (cl-mpm/particle::mp-mass mp)))))))

(defun make-bc-water-damage (sim datum rate)
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (with-accessors ((h cl-mpm/mesh:mesh-resolution))
        mesh
      (make-instance 'bc-water-damage
                     :index nil
                     :damage-rate rate
                     :damage-volume t
                     :sim sim))))
(defmethod cl-mpm/bc::apply-bc ((bc bc-water-damage) node mesh dt)
  (call-next-method)
  ;; (break)
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
                  (cl-mpm/particle::mp-damage mp)
                  (abs (*
                        (bc-water-damage-damage-rate bc)
                        weathering dt)))
                 (setf
                  (cl-mpm/particle::mp-damage mp)
                  (min
                   (cl-mpm/particle::mp-damage mp)
                   1d0)))))))

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



(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-chalk-delayed) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  ;; (cl-mpm::update-domain-det mesh mp dt)
  ;; (cl-mpm::update-domain-stretch mesh mp dt)
  ;; (cl-mpm::update-domain-corner mesh mp dt)
  ;; (cl-mpm::co-domain-corner-2d mesh mp dt)
  (cl-mpm::update-domain-polar-2d mesh mp dt)
  ;; (cl-mpm::scale-domain-size mesh mp)
  )

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-chalk-delayed) dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-damaged mesh mp dt fbar)
  ;; (cl-mpm::update-stress-kirchoff mesh mp dt fbar)
  (cl-mpm::update-stress-kirchoff-dynamic-relaxation mesh mp dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-mapped-jacobian mesh mp dt fbar)
  ;; (cl-mpm::update-domain-det mesh mp)
  ;; (cl-mpm::update-stress-kirchoff-noscale mesh mp dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-det mesh mp dt fbar)
  ;; (cl-mpm::update-stress-kirchoff mesh mp dt fbar)
  ;; (cl-mpm::update-domain-corner mesh mp)
  ;; (cl-mpm::co-domain-corner-2d mesh mp dt)
  ;; (cl-mpm::update-domain-corner-2d mesh mp dt)
  ;; (cl-mpm::update-domain-polar-2d mesh mp dt)
  ;; (cl-mpm::scale-domain-size mesh mp)
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
  ;; (cl-mpm::update-stress-kirchoff-ugimp mesh mp dt fbar)
  (cl-mpm::update-stress-kirchoff mesh mp dt fbar)
  ;; (cl-mpm::update-stress-linear mesh mp dt fbar)
  )
(in-package :cl-mpm/particle)
;; (defmethod constitutive-model ((mp cl-mpm/particle::particle-chalk-delayed) strain dt)
;;   "Strain intergrated elsewhere, just using elastic tensor"
;;   (with-accessors ((de mp-elastic-matrix)
;;                    (stress mp-stress)
;;                    (E mp-E)
;;                    (nu mp-nu)
;;                    (phi mp-phi)
;;                    (psi mp-psi)
;;                    (c mp-c)
;;                    (plastic-strain mp-strain-plastic)
;;                    (ps-vm mp-strain-plastic-vm)
;;                    (strain mp-strain)
;;                    (yield-func mp-yield-func)
;;                    (soft mp-softening)
;;                    (enabled mp-enable-plasticity)
;;                    )
;;       mp
;;     (declare (double-float soft ps-vm E nu phi psi c))
;;     ;;Train elastic strain - plus trail kirchoff stress
;;     (setf stress (cl-mpm/constitutive::linear-elastic-mat strain de stress))
;;     (when enabled
;;       (let ((f-r t))
;;         (loop for i from 0 to 50
;;               while f-r
;;               do
;;                  (progn
;;                    (multiple-value-bind (sig eps-e f inc)
;;                        (cl-mpm/ext::constitutive-mohr-coulomb stress
;;                                                               de
;;                                                               strain
;;                                                               E
;;                                                               nu
;;                                                               phi
;;                                                               psi
;;                                                               c)
;;                      (setf f-r (> f 1d-5))
;;                      (setf
;;                       stress sig
;;                       strain eps-e
;;                       yield-func f)
;;                      (incf ps-vm inc))
;;                    (setf (mp-plastic-iterations mp) i)
;;                    (when (> soft 0d0)
;;                      (with-accessors ((c-0 mp-c-0)
;;                                       (phi-0 mp-phi-0)
;;                                       (psi-0 mp-psi-0)
;;                                       (c-r mp-c-r)
;;                                       (phi-r mp-phi-r)
;;                                       (psi-r mp-psi-r))
;;                          mp
;;                        (declare (double-float c-0 c-r phi-0 phi-r psi-0 psi-r))
;;                        (setf
;;                         c (+ c-r (* (- c-0 c-r) (exp (- (* soft ps-vm)))))
;;                         phi (atan (+ (tan phi-r) (* (- (tan phi-0) (tan phi-r)) (exp (- (* soft ps-vm))))))))
;;                      )))))
;;     stress))

(in-package :cl-mpm/examples/joss)

;; (let ((phi-r (* 15d0 (/ pi 180)))
;;       (phi-0 (* 50d0 (/ pi 180)))
;;       (c-0 26d3)
;;       (c-r 0d0)
;;       (soft (* 1d-2 1000d0))
;;       )
;;   (pprint (* (atan (+ (tan phi-r) (* (- (tan phi-0) (tan phi-r)) (exp (- soft ))))) (/ 180 pi)))
;;   (pprint (+ c-r (* (- c-0 c-r) (exp (- soft)))))
;;   )

(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-chalk-delayed) dt)
  (let ((damage-increment 0d0))
    (with-accessors (;(stress cl-mpm/particle::mp-undamaged-stress)
                     (trial-strain cl-mpm/particle::mp-trial-strain)
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
                     )
        mp
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
                    (cl-mpm/constitutive:linear-elastic-mat strain de)
                    ;; stress
                    (* angle (/ pi 180d0)))))
        ;;Delocalisation switch
        ;; (setf damage-increment
        ;;       (max 0d0
        ;;            (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile
        ;;             stress
        ;;             (* angle (/ pi 180d0)))))
        (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
        (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
        ))))


(declaim (notinline plot))
(defun plot (sim)
  ;; (plot-collapse-curves)
  ;; (vgplot:plot *time-data* *damage-data*)
  (cl-mpm::g2p (cl-mpm:sim-mesh *sim*)
               (cl-mpm:sim-mps *sim*)
               (cl-mpm:sim-dt *sim*)
               0d0
               :QUASI-STATIC)
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :deformed
   ;; :colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :xy))
   :trial t
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
                     :damage-volume t
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
  (let* ((sim (cl-mpm/setup:make-simple-sim
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               :sim-type
               ;; 'cl-mpm::mpm-sim-usf
               ;; 'cl-mpm/damage::mpm-sim-usl-damage
               ;; 'cl-mpm/damage::mpm-sim-damage
               'cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul
               :args-list (list
                           ;; :split-factor 0.5d0
                           :enable-fbar t
                           :enable-aggregate t
                           :vel-algo :QUASI-STATIC
                           :max-split-depth 4
                           :enable-length-localisation nil
                           :enable-split t))

              )
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 1.7d3)
         )
    (declare (double-float h density))
    (progn
      (let* ((E 1d9)
             (angle 50d0)
             (init-c 26d3)
             (init-stress (cl-mpm/damage::mohr-coloumb-coheasion-to-tensile init-c (* angle (/ pi 180))))
             (gf (* 4.8d0 1d0))
             (length-scale (* h 1d0))
             (ductility (estimate-ductility-jirsek2004 gf length-scale init-stress E))
             (oversize (cl-mpm/damage::compute-oversize-factor 0.999d0 ductility)))
        (format t "Estimated oversize ~F~%" oversize)
        (format t "Estimated D-r ~F~%" (cl-mpm/damage::damage-response-exponential (* oversize init-stress)
                                                                                   E
                                                                                   init-stress
                                                                                   ductility
                                                                                   ))
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

           ;'cl-mpm/particle::particle-chalk-delayed
           'cl-mpm/particle::particle-chalk-erodable
           :E E
           :nu 0.24d0
           :enable-damage t
           :enable-plasticity t
           :friction-angle angle
           :kt-res-ratio 1d0
           :kc-res-ratio 0d0
           ;; :g-res-ratio 0.51d0
           :g-res-ratio 0.90d0
           :peerlings-damage t
           :fracture-energy 3000d0
           :initiation-stress init-stress;18d3
           :delay-time 1d2
           :delay-exponent 2d0
           :ductility ductility
           :critical-damage 1d0;(- 1.0d0 1d-3)
           ;; :damage-domain-rate 0.9d0;This slider changes how GIMP update turns to uGIMP under damage

           :local-length length-scale
           :local-length-damaged 10d-10
           ;; :psi 0d0
           ;; :phi (* 42d0 (/ pi 180))
           ;; :c 131d3
           :psi (* 0d0 (/ pi 180))
           :phi (* angle (/ pi 180))
           :c (* init-c oversize)
           :psi-r (* 0d0 (/ pi 180))
           :phi-r (* 30d0 (/ pi 180))
           :c-r 0d0
           :softening 0d0

           :index 0

           :gravity -9.8d0
           :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
           ))))

      (setf (cl-mpm:sim-allow-mp-split sim) t)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      ;; (setf (cl-mpm::sim-velocity-algorithm sim) :PIC)
      ;; (setf (cl-mpm::sim-velocity-algorithm sim) :FLIP)
                                        ;(setf (cl-mpm::sim-velocity-algorithm sim) :BLEND-2ND-ORDER)
      (setf (cl-mpm::sim-velocity-algorithm sim) :BLEND)
      ;; (setf (cl-mpm::sim-velocity-algorithm sim) :FLIP)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-enable-fbar sim) t)
      (setf (cl-mpm/damage::sim-enable-length-localisation sim) t)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      ;; (let ((mass-filter (* density (expt h 2) 1d-2)))
      ;;   (format t "Mass filter: ~F~%" mass-filter)
      ;;   (setf (cl-mpm::sim-mass-filter sim) mass-filter))
      (cl-mpm/setup::set-mass-filter sim density :proportion 0d-3)
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
         (accelerate-target-time 1d0)
         (accelerate-mass-scale 1d4)
         (collapse-target-time 0.1d0)
         (collapse-mass-scale 1d0)
         (damage-enabled nil)
         (plasticity-enabled nil)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 0.5d0)
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
           (* 1d-1
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


    (let ((ms 1d0))
      (setf (cl-mpm:sim-mass-scale *sim*) ms)
      (setf (cl-mpm:sim-damping-factor *sim*)
            (* 0.1d0
               (cl-mpm/setup::estimate-critical-damping *sim*))))

    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_conv_~5,'0d.vtk" 0)) *sim*)

    (defparameter *collapse* nil)
    (cl-mpm/dynamic-relaxation:converge-quasi-static
     *sim*
     :dt-scale dt-scale
     :energy-crit 1d-4
     :oobf-crit 1d-4
     :kinetic-damping t
     :damping-factor 1d-4
     :substeps 50
     :conv-steps 1000
     :post-iter-step
     (lambda (i e o)
       ;; (cl-mpm/damage::calculate-damage *sim*)
       (plot *sim*)
       (vgplot:title (format nil "step:~D - KE ~E - OOBF ~E" i e o))
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
          (* 1d-4
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
                                       ;; (incf work (cl-mpm/dynamic-relaxation::estimate-power-norm *sim*))
                                       ;; (incf oobf (cl-mpm/dynamic-relaxation::estimate-oobf *sim*))
                                       ;; (incf energy-estimate (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
                                       ;; (incf oobf (cl-mpm::sim-stats-oobf *sim*))
                                       ;; (incf energy-estimate (cl-mpm::sim-stats-energy *sim*))
                                       ;; (incf work (cl-mpm::sim-stats-power *sim*))
                                       ))))


                         (setf oobf (cl-mpm::sim-stats-oobf *sim*))
                         (setf energy-estimate (cl-mpm::sim-stats-energy *sim*))
                         (setf work (cl-mpm::sim-stats-power *sim*))
                         ;; (setf work (cl-mpm/dynamic-relaxation::estimate-power-norm *sim*))
                         ;; (setf energy-estimate (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
                         ;; (setf oobf (cl-mpm/dynamic-relaxation::estimate-oobf *sim*))
                         ;; (setf
                         ;;  energy-estimate (/ energy-estimate substeps)
                         ;;  oobf (/ oobf substeps))
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

                         (incf *sim-step*)
                         (plot *sim*)
                         (vgplot:title (format nil "Time:~F - KE ~E - OOBF ~E - Work ~E - ~A"  (cl-mpm::sim-time *sim*) energy-estimate oobf work sim-state)))
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
  (let* ((rc 0.9d0)
         ;(rs 0.359d0)
         (rs 0.950d0)
         (ratio (/ (- 1d0 rs) (- 1d0 rc)))
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

    (let ((datum (+ soil-boundary 1d0) )
          (lower-datum soil-boundary)
          (upper-datum (+ 16d0 soil-boundary))
          )
      (cl-mpm:iterate-over-mps
       (cl-mpm:sim-mps *sim*)
       (lambda (mp)
         (when (<= (cl-mpm/utils:varef (cl-mpm/particle:mp-position mp) 1) soil-boundary)
           (setf (cl-mpm/particle::mp-erosion-modulus mp) 1000d0))))
      (defparameter *bc-erode*
        (cl-mpm/erosion::make-bc-erode
         *sim*
         :rate 10d0
         :scalar-func (lambda (pos) (expt (- 1d0 (/ (abs (- (cl-mpm/utils:varef pos 1) datum))
                                                    (- upper-datum datum)
                                                    ))
                                          8))
         :clip-func (lambda (pos) (and (>= (cl-mpm/utils:varef pos 1) lower-datum)
                                       (<= (cl-mpm/utils:varef pos 1) upper-datum)
                                       )))))
    (cl-mpm:add-bcs-force-list
     *sim*
     *bc-erode*
     )
    (setf (cl-mpm/buoyancy::bc-enable *bc-erode*) nil)
    ;; (cl-mpm:add-bcs-force-list
    ;;  *sim*
    ;;  (make-bc-water-damage *sim*
    ;;                        (+ soil-boundary 1d0)
    ;;                        (/ 1d0 100d0)))

    ;; (cl-mpm/setup::initialise-stress-self-weight
    ;;  *sim*
    ;;  shelf-height
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
      ;;Pre-cut failure plane
      (when nil
        (let ((cut-height (* 0.5d0 shelf-height-true))
              (cut-back-distance 0.5d0)
              (start-point (+ 0.5d0 (magicl:tref sloped-inflex-point 0 0)))
              (width                    ;(* 2d0 mesh-size)
                ;; 0.5d0
                (* 0.5d0 (cl-mpm/particle::mp-local-length (aref (cl-mpm:sim-mps *sim*) 0)))
                ))
          (cl-mpm/setup::apply-sdf *sim* (lambda (p) (cl-mpm/setup::line-sdf
                                                      (cl-mpm/utils:vector-from-list (list (magicl:tref p 0 0)
                                                                              (magicl:tref p 1 0)
                                                                              0d0))
                                                      (list (- start-point
                                                               (* cut-back-distance shelf-height-true))
                                                            (float shelf-height 0d0)
                                                            0d0)
                                                      (list start-point
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
                                       (setf (cl-mpm/particle:mp-damage mp)
                                             d)
                                       (let ((k (cl-mpm/damage::find-k-damage-mp mp d)))
                                         ;; (push (list mp k) *to-damage-mps*)
                                         (setf (cl-mpm/particle::mp-history-stress mp)
                                               k)
                                         ))
                                     (cl-mpm/damage::update-damage mp 1d-3)
                                     )))))

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


(defun setup-3d ()
  (let* ((mesh-size 0.5)
         (mps-per-cell 2)
         (shelf-height 15.0)
         (soil-boundary 2)
         (shelf-aspect 1)
         (runout-aspect 0.5)
         (shelf-length (* shelf-height shelf-aspect))
         (domain-length (+ shelf-length (* runout-aspect shelf-height)))
         (shelf-height-true shelf-height)
         (shelf-height (+ shelf-height soil-boundary))
         (depth-aspect 2.00d0)
         (depth 15d0)
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
                                  )
                                :refine 1)
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
                                   )
                                 :refine 1
                                 )
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

;; (defun setup-under (&key (mesh-size 0.5d0) (mps-per-cell 2) (height 15d0))
;;   (let* (
;;          ;; (mesh-size 0.5)
;;          ;; (mps-per-cell 2)
;;          (shelf-height ;15.0
;;            height
;;                        )
;;          (soil-boundary 2)
;;          (shelf-aspect 1.0)
;;          (runout-aspect 2.0)
;;          (shelf-length (* shelf-height shelf-aspect))
;;          (domain-length (+ shelf-length (* runout-aspect shelf-height)))
;;          (shelf-height-true shelf-height)
;;          (shelf-height (+ shelf-height soil-boundary))
;;          (depth 400)
;;          (offset (list 0 (* 0 mesh-size)
;;                        ;; 0
;;                        ))
;;          )
;;     (defparameter *sim*
;;       (setup-test-column (list domain-length
;;                                (+ shelf-height (* 5 mesh-size))
;;                                ;; depth
;;                                )
;;                          (list domain-length shelf-height
;;                                ;; depth
;;                                )
;;                          offset
;;                          (/ 1d0 mesh-size) mps-per-cell))

;;     ;;Refine around tip
;;     (dotimes (i 0)
;;       (dolist (dir (list :x
;;                          :y
;;                          ))
;;         (cl-mpm::split-mps-criteria
;;          *sim*
;;          (lambda (mp h)
;;            (when
;;                (and
;;                 (> (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
;;                    (- shelf-length
;;                       ;; 3d0
;;                       (* 0.5d0 shelf-height)
;;                       ))
;;                 (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
;;                    shelf-length)
;;                 (> (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
;;                    soil-boundary
;;                    )
;;                 ;; (< (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
;;                 ;;    (- shelf-height 13d0)
;;                 ;;    )
;;                 )
;;              dir
;;              )))))
;;     (loop for mp across (cl-mpm::sim-mps *sim*)
;;           do (setf (cl-mpm/particle::mp-split-depth mp) 0))

;;     (let* ((sloped-height (- (- shelf-height soil-boundary) 7d0))
;;            (measured-angle 78d0)
;;            (undercut-angle ;(- 82.5d0 90d0)
;;              0d0
;;              ;; (- measured-angle 90d0)
;;                            )
;;            ;; (undercut-angle 0d0)
;;            (normal (magicl:from-list (list
;;                                       (cos (- (* pi (/ undercut-angle 180d0))))
;;                                       (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1)))
;;            (sloped-inflex-point
;;              (magicl:from-list (list (- shelf-length (* (tan (- (* pi (/ undercut-angle 180d0)))) sloped-height))
;;                                      soil-boundary)
;;                                '(2 1) :type 'double-float)))
;;       (cl-mpm/setup::remove-sdf *sim*
;;                                (lambda (p)
;;                                  (if (and
;;                                       (> (magicl:tref p 1 0) soil-boundary))
;;                                      (if (< (magicl:tref p 1 0) (+ soil-boundary sloped-height))
;;                                          (cl-mpm/setup::plane-point-sdf
;;                                           (magicl:from-list (list (magicl:tref p 0 0)
;;                                                                   (magicl:tref p 1 0)) '(2 1))
;;                                           normal
;;                                           (magicl:from-list (list shelf-length soil-boundary)
;;                                                             '(2 1) :type 'double-float))

;;                                          (cl-mpm/setup::plane-point-sdf
;;                                           (magicl:from-list (list (magicl:tref p 0 0)
;;                                                                   (magicl:tref p 1 0)) '(2 1))
;;                                           (magicl:from-list (list 1d0 0d0) '(2 1)  :type 'double-float)
;;                                           sloped-inflex-point)
;;                                          )
;;                                      1d0)
;;                                  ))
;;       (let ((cut-height (* 0.0d0 shelf-height-true)))
;;         (cl-mpm/setup::damage-sdf
;;          *sim*
;;          (lambda (p)
;;            (if t ;(< (abs (- (magicl:tref p 2 0) (* 0.5d0 depth))) 20d0)
;;                (funcall
;;                 (cl-mpm/setup::rectangle-sdf
;;                  (list (- (magicl:tref sloped-inflex-point 0 0)
;;                           (* 0.15d0 shelf-height-true))
;;                        shelf-height)
;;                  (list (* 1.0d0 mesh-size) cut-height)
;;                  ) p)
;;                0.99d0)
;;            ))))
;;     ;; (let ((sand-layer 3d0))
;;     ;;   (cl-mpm/setup::apply-sdf
;;     ;;    *sim*
;;     ;;    (lambda (p)
;;     ;;      (if (and
;;     ;;           (> (magicl:tref p 1 0) (- soil-boundary sand-layer))
;;     ;;           (> (magicl:tref p 0 0) shelf-length))
;;     ;;          0d0
;;     ;;          1d0))
;;     ;;    (lambda (mp)
;;     ;;      (setf (cl-mpm/particle::mp-damage mp) 0.9d0)))
;;     ;;   )
;;     (let* ((notched-depth 0.0d0)
;;        (undercut-angle 45d0)
;;        (notched-width (* 0.5 notched-depth))
;;        (d-mid 0d0
;;               ;(/ depth 2)
;;               )
;;        (normal (magicl:from-list (list
;;                                   (cos (- (* pi (/ undercut-angle 180d0))))
;;                                   (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1)))

;;        )
;;       (cl-mpm/setup:remove-sdf *sim*
;;                                (lambda (p)
;;                                  (if (and
;;                                       (> (magicl:tref p 1 0) soil-boundary)
;;                                       (> (magicl:tref p 2 0) (- d-mid notched-width))
;;                                       (< (magicl:tref p 2 0) (+ d-mid notched-width))
;;                                       )
;;                                      (let
;;                                          (
;;                                           (sloped-inflex-point
;;                                             (magicl:from-list (list
;;                                                                (- shelf-length (* notched-depth
;;                                                                                   (- 1d0 (/ (abs (- (magicl:tref p 2 0) d-mid)) notched-width))
;;                                                                                   ))
;;                                                                soil-boundary)
;;                                                               '(2 1) :type 'double-float))
;;                                           )
;;                                        (cl-mpm/setup::plane-point-sdf
;;                                         (magicl:from-list (list (magicl:tref p 0 0)
;;                                                                 (magicl:tref p 1 0)) '(2 1))
;;                                         normal
;;                                         sloped-inflex-point)
;;                                        )
;;                                      1d0)
;;                                  ))
;;   )

;;      (let ((notch-height 0.0d0))
;;        (cl-mpm/setup:remove-sdf *sim*
;;                                 (lambda (p)
;;                                   (if t ;(< (abs (- (magicl:tref p 2 0) (* 0.5d0 depth))) 20d0)
;;                                       (funcall
;;                                        (cl-mpm/setup::rectangle-sdf
;;                                         (list shelf-length (+ notch-height soil-boundary))
;;                                         (list (* 2d0 notch-height) notch-height)
;;                                         ) p)
;;                                       1d0)
;;                                   )))
;;      (setf cl-mpm::*max-split-depth* 6)

;;      ;; (let ((ratio 1.0d0))
;;      ;;   (cl-mpm/setup::damage-sdf *sim* (lambda (p) (cl-mpm/setup::line-sdf
;;      ;;                                                (magicl:from-list (list (magicl:tref p 0 0)
;;      ;;                                                                        (magicl:tref p 1 0)) '(2 1))
;;      ;;                                                (list (- shelf-length (* shelf-height ratio)) shelf-height)
;;      ;;                                                (list shelf-length soil-boundary)
;;      ;;                                                10d0
;;      ;;                                                )) 1.0d0))
;;      )
;;     ;; (let ((upper-random-bound 0.5d0))
;;     ;;   (loop for mp across (cl-mpm:sim-mps *sim*)
;;     ;;         do (setf (cl-mpm/particle::mp-damage mp)
;;     ;;                  (reduce #'*
;;     ;;                          (loop for i from 0 to 2
;;     ;;                                collect (random upper-random-bound))))))
;;     (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
;;     (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
;;     (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
;;     (defparameter *run-sim* t)
;;     (defparameter *t* 0)
;;     (defparameter *oobf* 0)
;;     (defparameter *energy* 0)
;;     (defparameter *sim-step* 0))

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

(defun test-mc ()
  (let* ((E 1d0)
         (nu 0.2d0)
         (phi 0.1d0)
         (psi 0d0)
         (c 1000d0)
         (de (cl-mpm/constitutive:linear-elastic-matrix E nu))
         (eps (cl-mpm/utils:voigt-from-list (list -10d0
                                                  -10d0
                                                  3d0
                                                  10d0
                                                  0d0
                                                  100d0)))
         (stress (magicl:@ de eps))
         (iters 1000)
         (sub-iters 1000)
         )
    (multiple-value-bind (sig eps f) (cl-mpm/constitutive::mc-plastic stress de eps E nu phi psi c)
      (print f)
      (pprint eps))
    (multiple-value-bind (sig eps f) (cl-mpm/ext::constitutive-mohr-coulomb stress de eps E nu phi psi c)
      (print f)
      (pprint eps))
    ;; (time
    ;;  (lparallel:pdotimes (i iters)
    ;;    (dotimes (j sub-iters)
    ;;      (cl-mpm/constitutive::mc-plastic stress de eps E nu phi psi c))))
    (time
     (lparallel:pdotimes (i iters)
       (dotimes (j sub-iters)
         (cl-mpm/ext::constitutive-mohr-coulomb stress de eps E nu phi psi c))))
    ))

;; (dotimes (i 1000)
;;   (let* ((E (random 100d0))
;;          (nu (random 0.4999d0))
;;          (phi  (random 60d0))
;;          (phi-rad (tan (* (/ pi 180) phi)))
;;          (psi 0d0)
;;          (c (random 200d0))
;;          (de (cl-mpm/constitutive:linear-elastic-matrix E nu))
;;          (scale (random 100d0))
;;          (eps (cl-mpm/utils:voigt-from-list (list (- (random scale) (* 2 scale))
;;                                                   (- (random scale) (* 2 scale))
;;                                                   (- (random scale) (* 2 scale))
;;                                                   (- (random scale) (* 2 scale))
;;                                                   (- (random scale) (* 2 scale))
;;                                                   (- (random scale) (* 2 scale))
;;                                                   )))
;;          (stress (magicl:@ de eps))
;;          (iters 100000))
;;     (let ((fdp
;;             (cl-mpm/constitutive::dp-yield-mc-circumscribe stress phi c)
;;             ;; (cl-mpm/constitutive::dp-yield-mc-dp2 stress phi-rad c)
;;             )
;;           (fmc (cl-mpm/constitutive::fast-mc stress phi-rad c)))
;;       ;; (format t "DP : ~E~%" fdp)
;;       ;; (format t "MC : ~E~%" fmc)
;;       (when (and (> fdp 0d0)
;;                  (< fmc 0d0))
;;         (format t "~%")
;;         (format t "DP : ~E~%" fdp)
;;         (format t "MC : ~E~%" fmc)
;;         ;; (format t "E: ~A~%" E)
;;         ;; (format t "nu: ~A~%" nu)
;;         (format t "phi: ~A~%" phi)
;;         ;; (format t "c: ~A~%" c)
;;         ;; (format t "eps: ~A~%" eps)
;;         (pprint "Odd failure")))))


;; (let* ((E 2d0)
;;        (nu 0.49d0)
;;        (De3
;;          (cl-mpm/fastmaths::fast-scale!
;;           (cl-mpm/utils:matrix-from-list (list
;;                                           (- 1d0 nu) nu nu
;;                                           nu (- 1d0 nu) nu
;;                                           nu nu (- 1d0 nu)))
;;           (/ E (* (+ 1d0 nu) (- 1d0 (* 2d0 nu))))))
;;        (Ce (magicl:inv De3))
;;        (Ce-analytic (cl-mpm/fastmaths::fast-scale!
;;                      (cl-mpm/utils:matrix-from-list (list
;;                                                      1d0 (- nu) (- nu)
;;                                                      (- nu) 1d0 (- nu)
;;                                                      (- nu) (- nu) 1d0))
;;                      (/ 1d0 E)))
;;        )
;;   (pprint Ce)
;;   (pprint Ce-analytic))










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
       :energy-crit 1d-1
       :oobf-crit 1d-3
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

;; (let ((mp (aref (cl-mpm:sim-mps *sim*) 0))
;;       (total-strain 1d-2)
;;       (steps 1000)
;;       (strain 0d0)
;;       )
;;   (setf (cl-mpm/particle::mp-enable-plasticity mp) t)
;;   (defparameter *steps* (list))
;;   (defparameter *strains* (list))
;;   (defparameter *times* (list))
;;   (time
;;    (loop for i from 0 to steps
;;          while *run-sim*
;;          do
;;             (let ((t0 (get-internal-run-time)))
;;               (dotimes (j 100)
;;                 (let ((s (cl-mpm/utils:voigt-from-list (list strain 0d0 0d0
;;                                                              0d0 0d0 strain))))
;;                   (setf (cl-mpm/particle::mp-strain mp) s)
;;                   (cl-mpm/particle:constitutive-model mp s 1d-4)))
;;               (incf strain (/ total-strain steps))
;;               (let ((dt (- (get-internal-run-time) t0)))
;;                 (push i *steps*)
;;                 (push dt *times*)
;;                 (push strain *strains*)))))
;;   (vgplot:close-all-plots)
;;   (vgplot:plot *steps* *times*))

;; (let ((vec (cl-mpm/utils::vector-zeros))
;;       (iters 10000000000))
;;   (time
;;    (dotimes (i iters)
;;      (cl-mpm/utils::deep-copy vec)))
;;   (time
;;    (dotimes (i iters)
;;      (cl-mpm/utils::vector-zeros)))
;;   )

(defun test ()
  (setup :mps 2 :refine 1 :notch-length 2d0)
  (setf (cl-mpm::sim-mass-scale *sim*) 1d0)
  (setf (cl-mpm:sim-damping-factor *sim*)
        (* 0d0 (cl-mpm/setup::estimate-critical-damping *sim*)))
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk") *sim*)
  (let ((dt-scale 0.5d0)
        (steps 500))
    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (loop for i from 0 to 500
          while *run-sim*
          do
             (progn
               (format t "Step ~D/~D~%" i steps)
               (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" i)) *sim*)
               (cl-mpm/output:save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" i)) *sim*)
               (cl-mpm:update-sim *sim*)
               (plot)
               (setf (cl-mpm:sim-dt *sim*) (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
               (swank.live:update-swank)
               ))))



(defun test-static ()
  (let ((output-dir "./output/"))
    (uiop:ensure-all-directories-exist (list output-dir))
    (setup :mps 3 :refine 1)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 0.1d0
             (cl-mpm/setup::estimate-critical-damping *sim*)))
    (cl-mpm/dynamic-relaxation::save-conv-preamble output-dir)
    (cl-mpm/dynamic-relaxation:converge-quasi-static
     *sim*
     :dt-scale 0.5d0
     :damping-factor 0.1d0
     :kinetic-damping nil
     :conv-steps 1000
     :oobf-crit 1d-5
     :energy-crit 1d-5
     :post-iter-step
     (lambda (i energy oobf)
       (cl-mpm/dynamic-relaxation::save-conv *sim* output-dir i)
       (cl-mpm/output:save-vtk (merge-pathnames (format nil "sim_~5,'0d.vtk" i) output-dir) *sim*)
       (cl-mpm/output:save-vtk-nodes (merge-pathnames (format nil "sim_nodes_~5,'0d.vtk" i) output-dir) *sim*)
       (plot *sim*))))
  )

(defun test-quasi-time ()
  (setup :mps 2 :refine 2 :notch-length 1d0)
  ;; (setup-3d)
  (cl-mpm/setup::set-mass-filter *sim* 1.7d3 :proportion 1d-15)
  (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
  (let ((dt 50d0)
        (total-time 1000d0))
    (setf (cl-mpm::sim-damping-factor *sim*) 0d0)
    (setf (cl-mpm::sim-velocity-algorithm *sim*) :QUASI-STATIC)
    (cl-mpm/dynamic-relaxation::run-quasi-time
     *sim*
     :output-dir "./output/"
     :dt-scale 1d0
     ;; :substeps 10
     :plotter (lambda (sim) (plot sim))
     :conv-criteria 1d-3
     :substeps 20
     :dt dt
     :steps (round total-time dt)
     ;; :enable-plastic nil
     ;; :time total-time
     )))

(defun test-multi ()
  (setup :mps 2 :refine 4 :notch-length 1d0)
  ;; (setup-3d)
  ;; (cl-mpm/setup::set-mass-filter *sim* 1.7d3 :proportion 1d-15)
  (cl-mpm/setup::set-mass-filter *sim* 1d0 :proportion 1d-9)
  ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
  (let ((dt 50d0)
        (total-time 1000d0))
    (setf (cl-mpm::sim-damping-factor *sim*) 0d0)
    (setf (cl-mpm::sim-velocity-algorithm *sim*) :QUASI-STATIC)
    (cl-mpm/dynamic-relaxation::run-multi-stage
     *sim*
     :output-dir "./output/"
     :dt-scale 1d0
     ;; :substeps 10
     :plotter (lambda (sim) (plot sim))
     :conv-criteria 1d-3
     :substeps 20
     :dt dt
     :max-adaptive-steps 12
     :steps (round total-time dt)
     ;; :explicit-dt-scale 100d0
     ;; :explicit-dynamic-solver 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic
     :explicit-dt-scale 0.5d0
     :explicit-dynamic-solver 'cl-mpm/damage::mpm-sim-agg-damage
     ;; :enable-plastic nil
     ;; :time total-time
     )))


;; (defun test-real-time ()
;;   (setup :mps 2 :refine 1)
;;   (let ((dt 1d0)
;;         (total-time 1000d0)
;;         (output-dir "./output/")
;;         )
;;     (cl-mpm/output::save-simulation-parameters (merge-pathnames "./output/" "settings.json")
;;                                                *sim*
;;                                                (list :dt 1d0
;;                                                      :criteria-energy 1d-3
;;                                                      :criteria-oobf 1d-3
;;                                                      :criteria-hist 1d0
;;                                                      ))
;;     (setf (cl-mpm::sim-damping-factor *sim*) 0d0)
;;     (setf (cl-mpm::sim-velocity-algorithm *sim*) :BLEND)
;;     (change-class *sim* 'cl-mpm/damage::mpm-sim-agg-damage)
;;     (cl-mpm/dynamic-relaxation::save-timestep-preamble output-dir)
;;     (cl-mpm/dynamic-relaxation::save-conv-preamble output-dir)
;;     (cl-mpm/dynamic-relaxation::step-real-time *sim* 0
;;                     :output-dir "./output/"
;;                     :plotter (lambda (sim) (plot sim))
;;                     :dt-scale 0.9d0
;;                     :target-time (* 0.5d0 dt)
;;                     :enable-damage t
;;                     :enable-plastic t)))

(defun test-implicit-dynamic ()
  (setup :mps 2 :refine 1 :notch-length 4d0)
  (let* ((dt 10d0)
         (total-time 1000d0)
         (output-dir "./output/")
         (target-time 10d0)
         (substeps (floor target-time dt)))
    (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic)
    (uiop:ensure-all-directories-exist (list output-dir))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (cl-mpm/dynamic-relaxation::save-timestep-preamble output-dir)
    (cl-mpm/dynamic-relaxation::save-conv-preamble output-dir)
    (cl-mpm/dynamic-relaxation::elastic-static-solution *sim* :substeps 100)
    (setf (cl-mpm:sim-dt *sim*) dt)
    (setf (cl-mpm:sim-enable-damage *sim*) t)
    (loop for steps from 0 to 1000
     while *run-sim*
     do
     (let ((energy 0d0)
           (work 0d0)
           (oobf 0d0))
       (cl-mpm/dynamic-relaxation::save-timestep *sim* output-dir steps :DYNAMIC)
       (format t "Step ~d ~%" steps)
       (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" *sim-step*)) *sim*)
       (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
       (time
        (dotimes (i substeps)
          (cl-mpm::update-sim *sim*)))
       (incf *sim-step*)
       (plot *sim*)
       (vgplot:title (format nil "Time:~F"  (cl-mpm::sim-time *sim*)))
       (swank.live:update-swank)
       ))))


;; (defparameter *i* 7)
;; (let ((i *i*)
;;       (sim *sim*))
;;   (cl-mpm/output:save-vtk (merge-pathnames "./output/" (format nil "sim_step_~5,'0d.vtk" i)) sim)
;;   (cl-mpm/output:save-vtk-nodes (merge-pathnames "./output/" (format nil "sim_step_nodes_~5,'0d.vtk" i)) sim)
;;   (incf *i*))
