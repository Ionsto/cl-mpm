(defpackage :cl-mpm/examples/slump
  (:use :cl))
(in-package :cl-mpm/examples/slump)
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)
(declaim (optimize (debug 0) (safety 0) (speed 3)))
;(sb-ext:restrict-compiler-policy 'speed  0 0)
;(sb-ext:restrict-compiler-policy 'debug  3 3)
;(sb-ext:restrict-compiler-policy 'safety 3 3)
;; (setf *block-compile-default* t)


(defun cl-mpm/damage::length-localisation (local-length local-length-damaged damage)
  ;; (+ (* local-length (- 1d0 damage)) (* local-length-damaged damage))
  (* local-length (max (sqrt (- 1d0 damage)) 1d-10))
  ;; local-length
  )
;; (defmethod print-object ((object magicl:matrix) stream)
;;   (pprint object stream))

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-visco-elasto-plastic-damage) dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-damaged mesh mp dt fbar)
  (cl-mpm::update-stress-kirchoff mesh mp dt fbar)
  (cl-mpm::co-domain-corner-2d mesh mp dt)
  (cl-mpm::scale-domain-size mesh mp)
  )


(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-visco-elasto-plastic-damage) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
                     (strain cl-mpm/particle::mp-strain)
                     (strain-plastic cl-mpm/particle::mp-strain-plastic)
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
        ;; (setf damage-increment (cl-mpm/damage::tensile-energy-norm-pressure strain E nu de
        ;;                                                                     (* damage pressure 0d0)
        ;;                                                                     ;; pressure
        ;;                                                                     ;; 0d0
        ;;                                                                     ))
        ;; (setf damage-increment
        ;;       (* ;(- 1d0 damage)
        ;;          (cl-mpm/damage::tensile-energy-norm-pressure
        ;;           (cl-mpm/fastmaths:fast-.+ strain strain-plastic)
        ;;           E nu
        ;;           de
        ;;           (* damage pressure 0d0))))
        ;; (setf damage-increment (cl-mpm/damage::criterion-j2 stress))
        ;; (setf damage-increment (cl-mpm/damage::criterion-max-principal-stress stress))

        ;; (setf damage-increment
        ;;       (* 
        ;;        ;(- 1d0 damage)
        ;;        1d0
        ;;        (cl-mpm/damage::tensile-energy-norm strain E de)))
        ;; (setf damage-increment
        ;;       (max 0d0
        ;;            (cl-mpm/damage::jrucker-prager-criterion
        ;;             (magicl:scale stress (/ 1d0 (magicl:det def)))
        ;;             (* 50d0 (/ pi 180d0)))))
        ;; (setf damage-increment
        ;;       (max 0d0
        ;;            (cl-mpm/damage::criterion-dp
        ;;             (magicl:scale stress (/ 1d0 (magicl:det def)))
        ;;             (* angle (/ pi 180d0)))))

        (let ((total-stress (cl-mpm/constitutive:linear-elastic-mat
                             (cl-mpm/fastmaths:fast-.+
                              strain
                              (magicl:scale strain-plastic (- 1d0 damage))) de)))
          (setf damage-increment
                (max 0d0
                     (cl-mpm/damage::criterion-dp-pressure
                      ;; (magicl:scale total-stress (/ 1d0 (magicl:det def)))
                      stress
                      ;; total-stress
                      (* angle (/ pi 180d0))
                      0d0))))

        ;; (setf damage-increment
        ;;       (max 0d0
        ;;            (cl-mpm/damage::criterion-dp-pressure
        ;;             ;(magicl:scale stress (/ 1d0 (magicl:det def))) (* angle (/ pi 180d0))
        ;;             stress
        ;;             (* angle (/ pi 180d0))
        ;;             ;; (* (- pressure) damage)
        ;;             ;; (- pressure)
        ;;             0d0
        ;;             )))
        
        ;;Delocalisation switch
        (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
        (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
        ))))



(in-package :cl-mpm/examples/slump)
(defclass bc-water-damage (cl-mpm/buoyancy::bc-scalar)
  ((damage-rate
    :initform 1d0
    :initarg :damage-rate
    :accessor bc-water-damage-damage-rate)))

(defun make-bc-water-damage (sim datum rate)
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (with-accessors ((h cl-mpm/mesh:mesh-resolution))
        mesh
      (make-instance 'bc-water-damage
                     :index nil
                     :damage-rate rate
                     :sim sim))))
(defmethod cl-mpm/bc::apply-bc ((bc bc-water-damage) node mesh dt)
  (call-next-method)
  (when (cl-mpm::sim-enable-damage *sim*)
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
                    ;(cl-mpm/particle::mp-history-stress mp)
                    (cl-mpm/particle::mp-damage mp)
                    (abs (*
                          (bc-water-damage-damage-rate bc)
                          weathering dt)))
                   (setf
                    (cl-mpm/particle::mp-damage mp)
                    (min
                      (cl-mpm/particle::mp-damage mp)
                      1d0))))))))

(declaim (optimize (debug 3) (safety 3) (speed 0)))

;; (pushnew cl-mpm-pic *features*)
;; (pushnew :cl-mpm-fbar *features*)
;; (remove :cl-mpm-fbar *features*)
;; (setf *features* (delete :cl-mpm-pic *features*))
;; (asdf:compile-system :cl-mpm :force T)
;; (asdf:compile-system :cl-mpm)

(defun melt-sdf (sim sdf meltrate)
  (setf (cl-mpm:sim-mps sim)
        (lparallel:premove-if (lambda (mp)
                                (with-accessors ((pos cl-mpm/particle:mp-position)) mp
                                  (if t;(>= 0 (funcall sdf pos))
                                      (progn
                                        ;;Apply melt
                                        (with-accessors ((mass cl-mpm/particle:mp-mass)
                                                         (volume cl-mpm/particle:mp-volume)
                                                         (temp cl-mpm/particle::mp-boundary))
                                            mp
                                          (let ((density (/ mass volume)))
                                            ;;Apply some mass-meltrate
                                            (decf mass (* volume (abs (min 0 temp)) meltrate))
                                            ;;Keep density
                                            ;(setf volume (/ mass density))
                                            ;;If we overmelted, then remove particle
                                            (if (<= mass 1d2)
                                                t
                                                nil
                                                ))))
                                      nil)))
                              (cl-mpm:sim-mps sim))))

(defun apply-pullout (sim load-mps shear-rate)
  (with-accessors ((mps cl-mpm:sim-mps)

                   (mesh cl-mpm:sim-mesh))
      sim
    (loop for mp in load-mps
          do
             (cl-mpm::iterate-over-neighbours
              mesh
              mp
              (lambda (mesh mp node svp grad fsvp fgrad)
                (with-accessors (
                                 ;; (vel cl-mpm/particle:mp-velocity)
                                 )
                    mp
                  (with-accessors ((pos cl-mpm/mesh::node-position)
                                   (vel cl-mpm/mesh::node-velocity)
                                   (active cl-mpm/mesh::node-active)
                                   (acc cl-mpm/mesh::node-acceleration)
                                   )
                      node
                    (when active
                      (setf (magicl:tref vel 0 0)
                            shear-rate
                            ;; (magicl:tref vel 1 0)
                            ;; 0d0
                            )
                      )
                    )))))))

(defun max-v-sum (mp)
  (with-accessors ((vel cl-mpm/particle:mp-velocity))
      mp
    (magicl::sum (magicl:map #'abs vel))))

(defun find-max-cfl (sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps)
                   (dt cl-mpm:sim-dt))
      sim
    (let ((max-v (lparallel:preduce #'max
                                    (lparallel:pmapcar #'max-v-sum
                                                    mps))))
      (* dt (/ max-v (cl-mpm/mesh:mesh-resolution mesh))))))

(defun pescribe-velocity (sim load-mps vel)
  (let ((mps (cl-mpm:sim-mps sim)))
    (loop for mp in load-mps
          do
             (progn
               (loop for v in vel
                     for i from 0
                     do
                     (when v
                       (setf (magicl:tref (cl-mpm/particle:mp-velocity mp) i 0) v)
                       )
                     )))))
(defun increase-load (sim load-mps amount)
  (loop for mp in load-mps
        do (with-accessors ((pos cl-mpm/particle:mp-position)
                            (force cl-mpm/particle:mp-body-force)) mp
             ;(incf (magicl:tref force 0 0) amount)
             (magicl:.+ force amount force)
             )))

(defun length-from-def (sim mp dim)
  (let* ((mp-scale 2)
         (h-initial (magicl:tref (cl-mpm/particle::mp-domain-size mp) dim 0))
         ;(h-initial  (/ (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) mp-scale))
         )
    h-initial
    ;(* (magicl:tref (cl-mpm::mp-deformation-gradient mp) dim dim) h-initial)
    ))
(defun max-stress (mp)
  (declare (optimize (speed 0) (debug 3)))
  (multiple-value-bind (l v) (magicl:eig (cl-mpm::voight-to-matrix (cl-mpm/particle:mp-stress mp)))
    (cl-mpm/utils:get-stress (cl-mpm/particle:mp-stress mp) :xz)
    ;(apply #'max l)
    ;; (- (max 0d0 (apply #'max l))
    ;;    (max 0d0 (apply #'min l)))
    ;; (cl-mpm/fastmaths::voigt-tensor-reduce-simd (cl-mpm/particle::mp-velocity-rate mp))
    ;; (magicl:tref (cl-mpm/particle::mp-velocity-rate mp) 2 0)
    ;; (cl-mpm/particle::mp-damage-ybar mp)
    ;; (cl-mpm/constitutive::effective-strain-rate (cl-mpm/particle::mp-eng-strain-rate mp))
    ;; (cl-mpm/utils:get-stress (cl-mpm/particle:mp-stress mp) :xx)
    ;; (cl-mpm/particle::mp-time-averaged-visc mp)
    ;; (magicl:tref (cl-mpm/particle::mp-stress mp) 2 0)
    )
  )
(defun local-dist (sim mp)
  (with-accessors ((ll cl-mpm/particle::mp-true-local-length)) mp
    (with-accessors ((llm cl-mpm/particle::mp-true-local-length)) *dist-mp*
      ;(cl-mpm/damage::weight-func-mps mp *dist-mp* (* 0.5d0 (+ ll llm)))
      ;; (cl-mpm/damage::weight-func-mps mp *dist-mp* 500d0)
      ;; (cl-mpm/damage::weight-func-mps-damaged (cl-mpm:sim-mesh sim) mp *dist-mp* 500d0 )
      ;; (cl-mpm/damage::diff-damaged (cl-mpm:sim-mesh sim) mp *dist-mp*)
      (cl-mpm/damage::diff-squared mp *dist-mp*)
      ;; ll
      )))

(declaim (notinline plot))
(defun plot (sim)
  ;; (plot-disp-day)
  (plot-domain)
  )


(declaim (notinline plot-domain))
;; (defun plot (sim &optional (plot :damage)))
(defun plot-domain (sim &optional (plot :damage))
  (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms)))
    (vgplot:format-plot t "set object 1 rect from 0,0 to ~f,~f fc rgb 'blue' fs transparent solid 0.5 noborder behind" ms-x *water-height*)
    (vgplot:format-plot t "set style fill solid")
    (vgplot:format-plot t "set yrange [~f:~f]" *floor-datum* ms-y)
    (vgplot:format-plot t "set size ratio ~f" (/ (- ms-y *floor-datum*) ms-x))
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :colour-func #'cl-mpm/particle::mp-damage
     ;; :colour-func #'cl-mpm/particle::mp-true-visc
     ;; :colour-func (lambda (mp) (magicl:tref (cl-mpm/particle::mp-stress mp) 1 0))
     )
    ;; (vgplot:format-plot t "set yrange [~f:~f]" (* 2 h) ms-y)
    (let* ((normal (cl-mpm/penalty::bc-penalty-normal *floor-bc*))
           (datum (cl-mpm/penalty::bc-penalty-datum *floor-bc*))
           (dy (if (= (magicl:tref normal 1 0) 0)
                   0d0
                   (- (/ (magicl:tref normal 0 0) (magicl:tref normal 1 0)))))
           (offset datum )
           )
      ;; (pprint dy)
      (vgplot:format-plot t "replot (~f*x + ~f)~%" dy offset))
    ;; (vgplot:replot)
    ))


(defun remove-sdf (sim sdf)
      (setf (cl-mpm:sim-mps sim)
            (remove-if (lambda (mp)
                         (with-accessors ((pos cl-mpm/particle:mp-position)) mp
                           (>= 0 (funcall sdf pos))
                           ))
                       (cl-mpm:sim-mps sim))))

(defun damage-sdf (sim sdf &optional (d 1d0))
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
      (loop for mp across mps
            do (with-accessors ((pos cl-mpm/particle:mp-position)
                                 (damage cl-mpm/particle:mp-damage)) mp
                 ;; (setf damage (funcall sdf pos))
                  (when (>= 0 (funcall sdf pos))
                    (setf damage (coerce d 'double-float)))
                 ))))

(defun rectangle-sdf (position size)
  (lambda (pos)
      (let* ((position (magicl:from-list position '(2 1) :type 'double-float))
             (dist-vec (magicl:.- (magicl:map! #'abs (magicl:.- pos position))
                                  (magicl:from-list size '(2 1) :type 'double-float))))

        (+ (sqrt (magicl::sum
                  (magicl:map! (lambda (x) (* x x))
                               (magicl:map! (lambda (x) (max 0d0 x)) dist-vec))))
           (min (max (magicl:tref dist-vec 0 0)
                     (magicl:tref dist-vec 1 0)
                     ) 0d0)))))
(defun ellipse-sdf (position x-l y-l)
  (let ((aspect (/ x-l y-l)))
    (lambda (pos)
      (let* ((position (magicl:from-list position '(2 1) :type 'double-float))
             (dist-vec (magicl:.* (magicl:.- position pos) (magicl:from-list (list 1 aspect) '(2 1)
                                                                             :type 'double-float)))
             (distance (sqrt (magicl:tref (magicl:@ (magicl:transpose dist-vec)
                                                    dist-vec) 0 0))))
        (- distance x-l)))))
(defun line-sdf (position a b width)
  (let* ((start (magicl:from-list a '(2 1) :type 'double-float))
         (end (magicl:from-list b '(2 1) :type 'double-float))
         (pa (magicl:.- position start))
         (ba (magicl:.- end start))
         (h (min 1d0 (max 0d0 (/ (cl-mpm/fastmaths::dot pa ba)
                                 (cl-mpm/fastmaths::dot ba ba)
                                 ))))
         (v (magicl:.- pa (magicl:scale ba h)))
         )
    (- (sqrt (cl-mpm/fastmaths::dot v v)) width)))
(defun plane-sdf (position normal distance)
  (- distance (cl-mpm/fastmaths::dot position (cl-mpm/fastmaths::norm normal))))
(defun plane-point-sdf (position normal point)
  (let ((distance (cl-mpm/fastmaths::dot point (cl-mpm/fastmaths::norm normal))))
    (- distance (cl-mpm/fastmaths::dot position (cl-mpm/fastmaths::norm normal)))))

(declaim (notinline setup-test-column))
(defun setup-test-column (size block-size block-offset slope &optional (e-scale 1d0) (mp-scale 1d0)
                          &key (angle 0d0))
  (declare (optimize (speed 0)))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1 e-scale)
               (mapcar (lambda (s) (* s e-scale)) size)
               :sim-type 'cl-mpm/damage::mpm-sim-damage
               ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density *ice-density*)
         )
    (format t "~F ~F ~A mp count"
            e-scale
            mp-scale
            block-size
            (mapcar (lambda (e) (floor (* e e-scale mp-scale))) block-size)
            )

    ;; (incf (second block-offset) (* 2 h-x))

    (progn
      (let* ((length-scale (* 1 h))
             ;; (stress 0.3d6)
             (stress 10d3)
             (gf 5000d0)
             (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale stress 1d9))
             (ductility-ii (cl-mpm/damage::estimate-ductility-jirsek2004 (* 0.9d0 gf) length-scale stress 1d9))
             ;; (ductility 1d2)
             ;; (ductility 100d0)
             ;; (ductility 5d0)
             ;; (ductility-ii 5d0)
             ;(ductility 10d0)
             )
        (format t "~%Estimated ductility: ~F~%" ductility)
        (format t "~%Estimated ductility-II: ~F~%" ductility-ii)
        (when (< ductility 1d0)
          (error "Ductility can't be less than 1!"))
        (cl-mpm::add-mps
         sim
         (cl-mpm/setup::make-mps-from-list
          (cl-mpm/setup::make-block-mps-list
           ;; block-position
           block-offset
           block-size
           (mapcar (lambda (e) (floor (* e e-scale mp-scale))) block-size)
           density

           'cl-mpm/particle::particle-visco-elasto-plastic-damage
           :visc-factor 111d6
           :visc-power 3d0

           :E 1d9
           :nu 0.325d0

           :friction-angle 40.0d0

           :kt-res-ratio 1d-9
           :kc-res-ratio 1d0
           :g-res-ratio 1d-9

           :initiation-stress stress

           ;; :damage 1d0
           :delay-time 1d1
           :delay-exponent 1d0
           :ductility ductility
           :ductility-mode-2 ductility-ii
           :critical-damage 1d0;(- 1.0d0 1d-3)
           :damage-domain-rate 0.9d0;This slider changes how GIMP update turns to uGIMP under damage
           :local-length length-scale;(* 0.20d0 (sqrt 7))
           :local-length-damaged 10d-10

           :enable-plasticity nil
           :enable-damage nil
           :enable-viscosity nil
           :psi (* 00d0 (/ pi 180))
           :phi (* 40d0 (/ pi 180))
           :c 1000d3
           :gravity -9.8d0
           :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
           :index 0
           :angle angle
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
        (format t "Chalk plastic virgin angle: ~F~%"
                (* (/ 180 pi) angle-plastic))
        (format t "Chalk plastic residual angle: ~F~%"
                (* (/ 180 pi) angle-plastic-damaged)))
      (let ((mass-scale 1d4)
            (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim))))
        (setf
         (cl-mpm::sim-mass-scale sim) mass-scale
         (cl-mpm:sim-damping-factor sim)
              (* 0.1d0
                 (sqrt mass-scale)
                 (cl-mpm/setup::estimate-critical-damping sim))))
      (format t "Est dt: ~f~%" (cl-mpm/setup::estimate-elastic-dt sim))

      (setf (cl-mpm:sim-mass-filter sim) 1d-10)
      (setf (cl-mpm::sim-allow-mp-split sim) t)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-enable-fbar sim) nil)
      (setf (cl-mpm:sim-dt sim) 1d-8)
      (setf (cl-mpm:sim-bcs sim) (make-array 0))
      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-varfix
             (cl-mpm:sim-mesh sim)
             '(0 nil nil)
             '(0 nil nil)
             '(nil 0 nil)
             '(nil 0 nil)
             ;;Front back
             '(nil nil 0)
             '(nil nil 0)
             ;; '(0 nil 0)
             ;; '(0 nil 0)
             ))
      (format t "Bottom level ~F~%" h-y)

      (let* ((terminus-size (+ (second block-size) (* slope (first block-size))))
             (ocean-x 1000)
             ;; (ocean-y (+ h-y (* 0.0d0 terminus-size)))
             ;; (ocean-y (+ (second block-offset)
             ;;             (* terminus-size 0.9d0)))
             (ocean-y (+ (second block-offset)
                         (* terminus-size 1.0d0)
                         -100
                         ;; 50
                         ))
             ;; (ocean-y 0d0)
             (ocean-y (* (round ocean-y h-y) h-y))
             ;;          )
            ;(angle -1d0)
            )

        (pprint block-size)
        (format t "Ocean level ~a~%" ocean-y)
        (defparameter *water-height* ocean-y)

        (let ((floor-friction 0.5d0));0.8
          (defparameter *ocean-floor-bc*
            (cl-mpm/penalty::make-bc-penalty-point-normal
             sim
             (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
             (cl-mpm/utils:vector-from-list (list 0d0
                                                  (* 2 h-y)
                                                  0d0))
             (* 1d9 1d0)
             floor-friction))
          (let ((ang (* angle (/ pi 180))))
            (defparameter *floor-bc*
              (cl-mpm/penalty::make-bc-penalty-point-normal
               sim
               (cl-mpm/utils:vector-from-list (list (sin ang) (cos ang) 0d0))
               (cl-mpm/utils:vector-from-list (list (first block-offset) (second block-offset) 0d0))
               (* 1d9 1d-2)
               floor-friction))))
        (setf (cl-mpm/penalty::bc-penalty-damping *floor-bc*) 0.1d0)
        (defparameter *melt-bc*
          (make-bc-water-damage
           sim
           ocean-y
           1d0))

        ;; (loop for mp across (cl-mpm:sim-mps sim)
        ;;       do
        ;;          (let ((stress (cl-mpm/buoyancy::buoyancy-virtual-stress (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)
        ;;                                                                  ocean-y
        ;;                                                                  ;*water-density*
        ;;                                                                  900d0
        ;;                                                                  ))
        ;;                (scaler 1.0d0)
        ;;                )
        ;;            (cl-mpm/fastmaths:fast-scale!
        ;;             stress
        ;;             scaler)
        ;;            (setf (cl-mpm/particle::mp-strain mp)
        ;;                  (magicl:linear-solve (cl-mpm/particle::mp-elastic-matrix mp)
        ;;                                       stress)
        ;;                  (cl-mpm/particle:mp-stress mp) stress
        ;;                  )
        ;;            ))

        (setf (cl-mpm::sim-bcs-force-list sim)
              (list
               (cl-mpm/bc:make-bcs-from-list
                (list
                 (;; cl-mpm/buoyancy::make-bc-buoyancy-clip
                  cl-mpm/buoyancy::make-bc-buoyancy-body
                  sim
                  ocean-y
                  *water-density*
                  (lambda (pos)
                    (>= (magicl:tref pos 1 0) (* h-y 0))))))
               (cl-mpm/bc:make-bcs-from-list
                (list *floor-bc*)
                )
               ;; (cl-mpm/bc:make-bcs-from-list
               ;;  (list *ocean-floor-bc*)
               ;;  )
               ;; (cl-mpm/bc:make-bcs-from-list
               ;;  (list *melt-bc*))
               ))
        )
      (let ((normal (magicl:from-list (list (sin (- (* pi (/ angle 180d0))))
                                            (cos (+ (* pi (/ angle 180d0))))) '(2 1))))
        (defparameter *sliding-slope* (/ (- (magicl:tref normal 0 0))
                                         (magicl:tref normal 1 0)
                                         ))
        (defparameter *sliding-offset* (- h-y (* 0d0 *sliding-slope*))))
      sim)))

(defparameter *ice-density* 900d0)
(defparameter *water-density* 1000d0)
;; (defparameter *ice-density* 900)
;; (defparameter *water-density* 1000)
;Setup
(defun setup (&key (refine 1d0) (angle 0d0) (mps 2))
  (declare (optimize (speed 0)))
  (defparameter *run-sim* nil)
  (setf cl-mpm::*max-split-depth* 4)
  (let* ((mesh-size (/ 25 refine))
         (mps-per-cell mps)
         (slope 0d0)
         (shelf-height 200d0)
         (shelf-aspect 2.0)
         (runout-aspect 1.0)
         (shelf-length (* shelf-height shelf-aspect))
         (shelf-end-height (+ shelf-height (* (- slope) shelf-length )))
         (shelf-height-terminus shelf-height)
         (shelf-height shelf-end-height)
         (depth 200)
         (angle-offset (* shelf-length (sin (* angle (/ pi 180)))))
         (offset (list 0d0
                  ;(* 2d0 mesh-size)
                  (+
                   ;shelf-height
                   (* 2d0 mesh-size)
                   angle-offset)
                  0d0
                  )))
    (defparameter *floor-datum* (second offset))
    (defparameter *removal-point* (- (+ shelf-length (* runout-aspect shelf-height)) (* 2 mesh-size)))
    (defparameter *sim*
      (setup-test-column (list (+ shelf-length (* runout-aspect shelf-height))
                               (+ (* shelf-height 2)
                                  ;; (* shelf-height
                                  ;;    (cos (* angle (/ pi 180))))
                                  )
                               ;; depth
                               )
                         (list shelf-length
                               shelf-height
                               ;; depth
                               )
                         offset
                         slope
                         (/ 1 mesh-size)
                         mps-per-cell
                         :angle angle)
      )

    ;;Delete all the plotted frames
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
    ;; (loop for mp across (cl-mpm:sim-mps *sim*)
    ;;       do
    ;;          (setf (cl-mpm/particle:mp-damage mp) (random 0.2d0)))
    ;; (let ((size 50))
    ;;   (cl-mpm/setup:remove-sdf *sim*
    ;;                            (lambda (p)
    ;;                              (funcall
    ;;                               (cl-mpm/setup::rectangle-sdf
    ;;                                (list shelf-length
    ;;                                      (-
    ;;                                       (+ shelf-height-terminus (second offset))
    ;;                                       100))
    ;;                                (list
    ;;                                 size
    ;;                                 (/ size 2))) p))))
    ;; (cl-mpm/setup:remove-sdf *sim*
    ;;                          (lambda (p)
    ;;                            (funcall
    ;;                             (cl-mpm/setup::rectangle-sdf
    ;;                              (list shelf-length
    ;;                                    (+ shelf-height-terminus (second offset))
    ;;                                    )
    ;;                              (list
    ;;                               200d0
    ;;                               60d0)) p)))
    
    ;; (damage-sdf *sim* (lambda (p)
    ;;                     (line-sdf (magicl:from-list
    ;;                                (list
    ;;                                 (magicl:tref p 0 0)
    ;;                                 (magicl:tref p 1 0)
    ;;                                 )
    ;;                                '(2 1))
    ;;                               (list (- shelf-length shelf-height) shelf-height)
    ;;                               (list shelf-length 0d0)
    ;;                               mesh-size
    ;;                               )) 0.9d0)

    ;; (let ((width 50d0)
    ;;       (height (+ (second offset) shelf-height)))
    ;;   (cl-mpm/setup::apply-sdf *sim* (lambda (p)
    ;;                                    (cl-mpm/setup::line-sdf
    ;;                                     p
    ;;                                     (list 0d0 height 0d0)
    ;;                                     (list shelf-length height 0d0)
    ;;                                     width
    ;;                                     ))
    ;;                            (lambda (mp v)
    ;;                              (setf (cl-mpm/particle:mp-damage mp)
    ;;                                    (* 0.8d0 (exp (- (expt (/ (+ width v) width) 2))))
    ;;                                    )
    ;;                              )
    ;;                            ))

    ;; (let* ((normal (cl-mpm/fastmaths::norm
    ;;                 (magicl:from-list (list
    ;;                                    (- slope)
    ;;                                    1d0) '(2 1))))
    ;;        )
    ;;   (remove-sdf *sim*
    ;;               (lambda (p)
    ;;                 (plane-point-sdf
    ;;                  p
    ;;                  normal
    ;;                  (magicl:from-list (list shelf-length
    ;;                                          shelf-height-terminus)
    ;;                                    '(2 1) :type 'double-float)))))
    ;; (let* ((undercut-angle -70d0)
    ;;        (normal (magicl:from-list (list
    ;;                                   (cos (- (* pi (/ undercut-angle 180d0))))
    ;;                                   (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1)))
    ;;        (cutback (* shelf-height (/ 1d0 (cos (- (* pi (/ undercut-angle 180d0)))))))
    ;;        )
    ;;   (remove-sdf *sim* (lambda (p) (plane-point-sdf p
    ;;                                                  normal
    ;;                                                  (magicl:from-list (list (- shelf-length cutback)
    ;;                                                                          shelf-height)
    ;;                                                                    '(2 1) :type 'double-float)
    ;;                                                  ))))
    (let* ((undercut-angle 0d0)
           (normal (magicl:from-list (list
                                           (cos (- (* pi (/ undercut-angle 180d0))))
                                           (sin (- (* pi (/ undercut-angle 180d0))))
                                           0d0) '(3 1))))
      (cl-mpm/setup:remove-sdf *sim* (lambda (p) (plane-point-sdf
                                                  p
                                                  normal
                                                  (magicl:from-list (list shelf-length shelf-height 0)
                                                                    '(3 1) :type 'double-float)
                                                  ))))
    ;; (let* ((n (magicl:from-list (list shelf-height
    ;;                                   (* 3d0 shelf-length)) '(2 1) :type 'double-float))
    ;;       (d (plane-sdf (magicl:from-list (list 0d0 shelf-height) '(2 1) :type 'double-float) n 0d0))
    ;;       )
    ;;   (remove-sdf *sim* (lambda (p) (plane-sdf p
    ;;                                            n
    ;;                                            (- d)
    ;;                                            ))))
    )
  ;; (str:to-file #p"output/settings.json"
  ;;              (jonathan:to-json (list :ocean-height *water-height* :domain-size
  ;;                                      (cl-mpm/mesh::mesh-mesh-size (cl-mpm::sim-mesh *sim*))
  ;;                                      )))
  ;; (damage-sdf *sim* (lambda (p) (line-sdf p
  ;;                                         (list 0d0 0d0)
  ;;                                         (list 100d0 0d0)
  ;;                                         100d0
  ;;                                         )) 0.8))
  ;; (defparameter *sim* (setup-test-column '(1 1) '(1 1) '(0 0) 1 1))
  ;; (remove-sdf *sim* (ellipse-sdf (list 250 100) 20 40))
  ;; (remove-sdf *sim* (ellipse-sdf (list 2 50 100) 20 40))
  ;; (remove-sdf *sim* (ellipse-sdf (list 1.5 3) 0.25 0.5))
  (print (cl-mpm:sim-dt *sim*))
  (format t "Sim MPs: ~a~%" (length (cl-mpm:sim-mps *sim*)))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  (defparameter *x* 0d0)
  (defparameter *x-pos* '())
  (defparameter *cfl-max* '())
  (defparameter *sim-step* 0)
  ;(defparameter *run-sim* t)
  (defparameter *load-mps*
    (let* ((mps (cl-mpm:sim-mps *sim*))
           (least-pos
              (apply #'max (loop for mp across mps
                                 collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))))
      (loop for mp across mps when (>= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (- least-pos 0.001))
              collect mp)))
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (let ((x-min (loop for mp across mps
                       minimize (magicl:tref
                                 (cl-mpm/particle:mp-position mp)
                                 0 0))))
      (defparameter *far-field-mps*
        (loop for mp across mps
              when (= x-min (magicl:tref
                         (cl-mpm/particle:mp-position mp)
                         0 0))
                collect mp)))
    (let ((x-max (loop for mp across mps
                       maximize (magicl:tref
                                 (cl-mpm/particle:mp-position mp)
                                 0 0))))
      (defparameter *terminus-mps*
        (loop for mp across mps
              when (= x-max (magicl:tref
                             (cl-mpm/particle:mp-position mp)
                             0 0))
                collect mp)))
    ;; (increase-load *sim* *terminus-mps*
    ;;                (magicl:from-list (list (* 1d4) 0d0) '(2 1)))
    )
  (defparameter *dist-mp* (nth 0 *terminus-mps*))
  ;; (increase-load *sim* *load-mps* 1)
  ;; (increase-load *sim* *load-mps* 100)
  )
(defparameter *water-height* 0d0)


(defparameter *run-sim* nil)
(defun calculate-dt (courent c-value target-step)
  (let* (
         (cfl-dt (if (> courent 0d0) (/ c-value (/ courent (cl-mpm:sim-dt *sim*))) nil))
         (new-dt (/ c-value (if (> courent 0d0)
                                (/ courent (cl-mpm:sim-dt *sim*))
                                (/ c-value 1e-6))))
         (max-steps 1000)
         (sub-steps (max (min (floor (/ target-step new-dt)) max-steps) 1)))
    (when (> (floor (/ target-step new-dt)) max-steps)
        (format t "CFL requires more steps than max-steps~%"))
    (format t "C: ~f - steps: ~D - %dt: ~f~%" courent sub-steps new-dt)
    (format t "Cfl derived dt:~f~%" cfl-dt)
                                        ;(setf (cl-mpm:sim-dt *sim*) new-dt)
    (values new-dt sub-steps)))

(declaim (notinline run))
(defun run ()
  (declare (optimize (speed 0) (debug 3) (safety 3)))
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)
  (defparameter *run-sim* t)
    (vgplot:close-all-plots)
    (vgplot:figure)
  ;; (sleep 1)
  (with-open-file (stream (merge-pathnames "output/terminus_position.csv") :direction :output :if-exists :supersede)
    (format stream "Time (s),Terminus position~%")
    (loop for tim in (reverse *time*)
          for x in (reverse *x-pos*)
          do (format stream "~f, ~f ~%" tim x)))
  (let* ((target-time 1d1)
         (dt (cl-mpm:sim-dt *sim*))
         (dt-0 0d0)
         (dt-scale 1d0)
         (settle-steps 0)
         (damp-steps 0)
         (accelerate-target-time 1d2)
         (accelerate-mass-scale 1d6)
         (collapse-target-time 1d0)
         (collapse-mass-scale 1d4)
         (substeps (floor target-time dt))
         (sim-state :damp)
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
         (damping-0 (* 1d-3 (cl-mpm/setup::estimate-critical-damping *sim*)))
         (damage-0 0d0)
         (mass-0
           (lparallel:pmap-reduce #'cl-mpm/particle::mp-mass #'+ (cl-mpm:sim-mps *sim*))))

    (cl-mpm/output::save-simulation-parameters #p"output/settings.json"
                                               *sim*
                                               (list :dt target-time
                                                     :ocean-height *water-height*
                                                     ))

    ;; (setf (cl-mpm/penalty::bc-penalty-friction *floor-bc*) 0d0)
    (cl-mpm:iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (setf
        (cl-mpm/particle::mp-enable-viscosity mp) nil
        (cl-mpm/particle::mp-enable-plasticity mp) nil)))
    (setf (cl-mpm:sim-dt *sim*)
          (* dt-scale (cl-mpm/setup:estimate-elastic-dt *sim*)))
    (cl-mpm/dynamic-relaxation:converge-quasi-static
     *sim*
     :dt-scale dt-scale
     :energy-crit 1d-2
     :oobf-crit 1d-2
     :substeps 50
     :conv-steps 1000
     :dt-scale dt-scale
     :post-iter-step
     (lambda (i e o)
       (plot *sim*)
       (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_conv_~5,'0d.vtk" i)) *sim*)
       ))

    (cl-mpm::update-sim *sim*)
    (let* ((dt-e (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
           (substeps-e (floor target-time dt-e)))
      (format t "CFL dt estimate: ~f~%" dt-e)
      (format t "CFL step count estimate: ~D~%" substeps-e)
      (setf (cl-mpm:sim-dt *sim*) dt-e)
      (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
    (setf dt-0 (/ (cl-mpm:sim-dt *sim*) (sqrt (cl-mpm::sim-mass-scale *sim*))))

    (setf target-time accelerate-target-time
          (cl-mpm:sim-mass-scale *sim*) accelerate-mass-scale)
    (setf (cl-mpm:sim-dt *sim*)
          (* dt-scale (cl-mpm/setup:estimate-elastic-dt *sim*)))
    (setf substeps (floor target-time (cl-mpm:sim-dt *sim*)))

   (defparameter *oobf* 0)
   (defparameter *data-damage* 0)
   (defparameter *data-energy* 0)

   (defparameter *data-full-time* (list))
   (defparameter *data-full-disp* (list))
   (defparameter *data-full-damage* (list))
   (defparameter *data-full-oobf* (list))
   (defparameter *data-full-energy* (list))
   (with-open-file (stream (merge-pathnames "./output/timesteps.csv") :direction :output :if-exists :supersede)
		 (format stream "steps,time,damage,energy,oobf~%"))
   (let ((energy-crit 0d0)
         (power-crit 0d0))
     (time (loop for steps from 0 to 1000
                 while *run-sim*
                 do
                    (progn
                      (format t "Step ~d ~%" steps)
                      (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                      (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                      ;; (cl-mpm/output::save-vtk-cells (merge-pathnames (format nil "output/sim_cells_~5,'0d.vtk" *sim-step*)) *sim*)
                      (with-open-file (stream (merge-pathnames "output/timesteps.csv") :direction :output :if-exists :append)
                        (format stream "~D,~f,~f,~f,~f~%"
                                steps
                                *t*
                                *data-damage*
                                *data-energy*
                                *oobf*))
                      (let ((energy-estimate 0d0)
                            (work 0d0))
                        (time
                         (dotimes (i (1+ substeps))
                           (cl-mpm::update-sim *sim*)
                           ;; (incf work (cl-mpm/dynamic-relaxation::estimate-power-norm *sim*))
                           (incf work (cl-mpm/dynamic-relaxation::estimate-power-norm *sim*))
                           (incf *oobf* (cl-mpm/dynamic-relaxation::estimate-oobf *sim*))
                           (incf energy-estimate (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
                           ;; (when (> (abs work) 0d0)
                           ;;   (incf energy-estimate (/ (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*) work)))
                           (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))

                        (setf *oobf* (/ *oobf* substeps)
                              energy-estimate (/ energy-estimate substeps))
                        (setf work (/ (abs work) target-time))

                        (setf energy-estimate (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
                        (setf *oobf* (cl-mpm/dynamic-relaxation::estimate-oobf *sim*))
                        (setf energy-estimate (abs (/ energy-estimate work)))

                        ;; (setf energy-estimate (abs (/ energy-estimate work)))
                        ;; (setf energy-estimate (abs energy-estimate))

                        (setf *data-energy* energy-estimate)
                        ;; (setf *data-energy* energy-estimate)
                        (let ((damage-est
                                (-
                                 (lparallel:pmap-reduce (lambda (mp)
                                                          (*
                                                           1d0
					                                                 (cl-mpm/particle::mp-mass mp)
                                                           (cl-mpm/particle::mp-damage mp)))
                                                        #'+ (cl-mpm:sim-mps *sim*)
                                                        :initial-value 0d0)
                                 damage-0)))
                          (setf *data-damage* damage-est)
                          (format t "Total damage: ~E~%" damage-est))
                        (format t "Energy estimate: ~E~%" energy-estimate)
                        (format t "OOBF estimate: ~E~%" *oobf*)

                        (when (= steps damp-steps)
                          (setf sim-state :settle)
                          (setf (cl-mpm:sim-damping-factor *sim*)
                                damping-0))
                        (when (= steps settle-steps)
                          (setf (cl-mpm::sim-enable-damage *sim*) t)
                          (cl-mpm:iterate-over-mps
                           (cl-mpm:sim-mps *sim*)
                           (lambda (mp)
                             (setf
                              (cl-mpm/particle::mp-enable-viscosity mp) nil
                              (cl-mpm/particle::mp-enable-damage mp) t
                              (cl-mpm/particle::mp-enable-plasticity mp) t))))
                        (when (>= steps settle-steps)
                          (if (or
                               ;; (> energy-estimate 1d-3)
                               ;; (> *oobf* 1d-1)
                               nil
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
                        (print "calc adaptive")
                        (multiple-value-bind (dt-e substeps-e)
                            (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                          (format t "CFL dt estimate: ~f~%" dt-e)
                          (format t "CFL step count estimate: ~D~%" substeps-e)
                          (setf substeps substeps-e))
					              ;; (let* (;(dt-est (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
							          ;;        (dt-est (* dt-0 (sqrt (cl-mpm::sim-mass-scale *sim*))))
					 	            ;;        (substeps-est (floor target-time dt-est)))
					              ;;   (when (< substeps-est substeps)
					 	            ;;     (setf (cl-mpm:sim-dt *sim*) dt-est)
					 	            ;;     (setf substeps substeps-est)))
                        ;; (format t "CFL dt estimate: ~f~%" (cl-mpm:sim-dt *sim*))
                        ;; (format t "CFL step count estimate: ~D~%" substeps)
					              (cl-mpm/setup:remove-sdf
                         *sim*
                         (lambda (p)
                           (if (> (magicl:tref p 0 0) *removal-point*)
                               -1d0
                               1d0)))
                        (incf *sim-step*)
                        (push *t* *data-full-time*)
                        (push (lparallel:pmap-reduce
                               (lambda (mp) (cl-mpm/utils:varef (cl-mpm/particle::mp-displacement mp) 0))
                               #'max
                               (cl-mpm:sim-mps *sim*))
                              *data-full-disp*)
                        (push *data-energy* *data-full-energy*)
                        (push *oobf* *data-full-oobf*)
                        (plot *sim*)
                        (vgplot:title (format nil "Time:~F - Stage ~A - Energy ~E - OOBF ~E - Work ~E"  *t* sim-state energy-estimate *oobf* work))
                        (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                           :terminal "png size 1920,1080"
                                           )
                        )
                      
                      (format t "Mass loss: ~F~%"
                              (-
                               mass-0
                               (lparallel:pmap-reduce #'cl-mpm/particle::mp-mass #'+ (cl-mpm:sim-mps *sim*))))
                      (swank.live:update-swank)
                     )
                    ;; (progn
                    ;;   (let ((base-damping 1d-3))
                    ;;     (when (> steps settle-steps)
                    ;;       (setf (cl-mpm::sim-enable-damage *sim*) t))
                    ;;     (if (> steps damp-steps)
                    ;;         (if (> energy-crit 1d4)
                    ;;             (progn
                    ;;               (format t "Collapse ~%")
                    ;;               (setf
                    ;;                ;; (cl-mpm::sim-mass-scale *sim*) 1d3
                    ;;                target-time 1d0
                    ;;                (cl-mpm:sim-mass-scale *sim*) 1d0
                    ;;                ;; (cl-mpm::sim-damping-factor *sim*) 0d0;(* 1d4 0.001d0)
                    ;;                )
                    ;;               )
                    ;;             (progn
                    ;;               (format t "Accelerate ~%")
                    ;;               (setf 
                    ;;                ;; (cl-mpm::sim-mass-scale *sim*) 1d3
                    ;;                target-time 1d2
                    ;;                (cl-mpm:sim-mass-scale *sim*) 1d4
                    ;;                )
                    ;;               (setf
                    ;;                (cl-mpm::sim-damping-factor *sim*)
                    ;;                (*
                    ;;                 10d0
                    ;;                 damping-0))

                    ;;               ))
                    ;;         (progn
                    ;;           (setf (cl-mpm::sim-enable-damage *sim*) nil
                    ;;                 (cl-mpm:sim-mass-scale *sim*) 1d4)
                    ;;           (setf
                    ;;            (cl-mpm::sim-damping-factor *sim*)
                    ;;            (*
                    ;;             (cl-mpm:sim-mass-scale *sim*)
                    ;;             damping-0)))
                    ;;         ))
                    ;;   (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                    ;;     (format t "CFL dt estimate: ~f~%" dt-e)
                    ;;     (format t "CFL step count estimate: ~D~%" substeps-e)
                    ;;     (setf substeps substeps-e))

                    ;;   (setf (cl-mpm:sim-dt *sim*) (* dt-0 (sqrt (cl-mpm::sim-mass-scale *sim*))))
                    ;;   (setf substeps (floor target-time (cl-mpm:sim-dt *sim*)))
                    ;;   (format t "Actual dt estimate: ~f~%" (cl-mpm:sim-dt *sim*))
                    ;;   (format t "Actual step count estimate: ~D~%" substeps)
                    ;;   ;; (when (> steps 2)
                    ;;   ;;   (setf dt-scale 1d0)
                    ;;   ;;     )
                    ;;   (format t "Step ~d ~%" steps)
                    ;;   (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                    ;;   (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                    ;;   ;; (cl-mpm/output:save-csv (merge-pathnames (format nil "output/sim_~5,'0d.csv" *sim-step*)) *sim*)

                    ;;   (push *t* *time*)
                    ;;   ;; (let ((cfl (find-max-cfl *sim*)))
                    ;;   ;;   (push cfl *cfl-max*))
                    ;;   (setf *x*
                    ;;         (loop for mp across (cl-mpm:sim-mps *sim*)
                    ;;               maximize (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0)))
                    ;;   (push
                    ;;    *x*
                    ;;    *x-pos*)
                    ;;   (setf energy-crit 0d0
                    ;;         power-crit 0d0)
                    ;;   (let ((cfl 0))
                    ;;     (time (dotimes (i substeps)
                    ;;             ;; (increase-load *sim* *terminus-mps*
                    ;;             ;;                (magicl:from-list (list (* (cl-mpm:sim-dt *sim*) 1d4) 0d0) '(2 1)))
                    ;;             ;; (pescribe-velocity *sim* *terminus-mps* '(1d0 nil))

                    ;;             ;; (melt-sdf *sim* (rectangle-sdf (list 0 0) (list 5000 *water-level*))
                    ;;             ;;           (* (cl-mpm:sim-dt *sim*) 1d-3))
                    ;;             (cl-mpm::update-sim *sim*)
                    ;;             ;(incf energy-crit (* dt (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*)))
                    ;;             ;(incf energy-crit (* dt (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*)))
                    ;;             ;; (incf energy-crit (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*)))
                    ;;             (incf energy-crit (estimate-energy-crit *sim*))
                    ;;             (incf power-crit (estimate-power-crit *sim*))
                    ;;             ;; (setf cfl (max cfl (find-max-cfl *sim*)))
                    ;;             (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))
                    ;;     (setf energy-crit (/ energy-crit substeps))
                    ;;     ;; (loop for mp across (cl-mpm:sim-mps *sim*)
                    ;;     ;;       do
                    ;;     ;;          (with-accessors ((tv cl-mpm/particle::mp-time-averaged-visc))
                    ;;     ;;              mp
                    ;;     ;;            (setf tv (/ tv substeps))))
                    ;;     ;; (setf cfl (find-max-cfl *sim*))
                    ;;     (format t "CFL: ~f~%" cfl)
                    ;;     (push cfl *cfl-max*)
                    ;;     (format t "Energy norm: ~F~%" energy-crit)
                    ;;     ;; (let* ((dt-e (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
                    ;;     ;;        (substeps-e (floor target-time dt-e)))
                    ;;     ;;   (format t "CFL dt estimate: ~f~%" dt-e)
                    ;;     ;;   (format t "CFL step count estimate: ~D~%" substeps-e)
                    ;;     ;;   (setf (cl-mpm:sim-dt *sim*) dt-e)
                    ;;     ;;   (setf substeps substeps-e))
                    ;;     )

                    ;;   (cl-mpm/setup:remove-sdf
                    ;;    *sim*
                    ;;    (lambda (p)
                    ;;      (if (> (magicl:tref p 0 0) *removal-point*)
                    ;;          -1d0
                    ;;          1d0)))

                    ;;   (with-open-file (stream (merge-pathnames "output/terminus_position.csv") :direction :output :if-exists :append)
                    ;;     (format stream "~f, ~f ~%" *t* *x*))
                    ;;   (incf *sim-step*)
                    ;;   (plot *sim*)
                    ;;   (vgplot:title (format nil "Time:~F - Energy ~E - Power ~E"  *t* energy-crit power-crit))
                    ;;   (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                    ;;                      :terminal "png size 1920,1080"
                    ;;                      )
                    ;;   (swank.live:update-swank)
                    ;;   (sleep .1)
                    ;;   )
                 ))))
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
  ;; (plot-s-xx *far-field-mps* 125d0)
  ;(plot-disp)
  ;; (with-open-file (stream (merge-pathnames "output/far-field-stress.csv") :direction :output :if-exists :append)
  ;;   (format stream "~f, ~f ~%" *t* *x*)
  ;;   )
;  (cl-mpm/output:save-csv (merge-pathnames (format nil "output/sim_csv_~5,'0d.vtk" *sim-step*)) *sim*)
  ;; (vgplot:figure)
  ;; (vgplot:title "Terminus over time")
  ;; (vgplot:plot *time* *x-pos*)
  ;; (plot-cfl)

  (with-open-file (stream (merge-pathnames "output/terminus_position.csv") :direction :output :if-exists :supersede)
    (format stream "Time (s),Terminus position~%")
    (loop for tim in (reverse *time*)
          for x in (reverse *x-pos*)
          do (format stream "~f, ~f ~%" tim x)))

  (with-open-file (stream (merge-pathnames "output/far-field.csv") :direction :output :if-exists :supersede)
    (format stream "y,s_xx,s_an~%")
    (let* ((mps *far-field-mps*)
           (H 125)
           (rho-ice 900)
           ;(E (cl-mpm/particle::mp-e (first mps)))
           (nu (cl-mpm/particle::mp-nu (first mps)))
           (g (cl-mpm/particle::mp-gravity (first mps)))
           (s-xx (loop for mp in mps
                       collect
                       (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0)))
           (y (loop for mp in mps
                    collect
                    (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)))
           (s-an (mapcar
                  (lambda (y)
                    (* (/ nu (- 1 nu)) rho-ice g -1d0 (- y (/ H 2d0))))
                  y)))
      (loop for yi in y
            for si in s-xx
            for sai in s-an
            do
               (format stream "~f, ~f, ~f ~%" yi si sai))))
  )

;; (defmethod cl-mpm::update-node-forces ((sim cl-mpm::mpm-sim))
;;   (with-accessors ((damping cl-mpm::sim-damping-factor)
;;                    (mass-scale cl-mpm::sim-mass-scale)
;;                    (mesh cl-mpm::sim-mesh)
;;                    (dt cl-mpm::sim-dt))
;;       sim
;;     (cl-mpm::iterate-over-nodes
;;      mesh
;;      (lambda (node)
;;        ;; (cl-mpm::calculate-forces-cundall node damping dt mass-scale)
;;        ;; (cl-mpm::calculate-forces-cundall-conservative node damping dt mass-scale)
;;        ;; (cl-mpm::calculate-forces node damping dt mass-scale)
;;        ))))

(defun run-static ()
  (declare (optimize (speed 0) (debug 3) (safety 3)))
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk") *sim*)
  (defparameter *run-sim* t)
    (vgplot:close-all-plots)
    (vgplot:figure)

  (with-open-file (stream (merge-pathnames "output/terminus_position.csv") :direction :output :if-exists :supersede)
    (format stream "Time (s),Terminus position~%")
    (loop for tim in (reverse *time*)
          for x in (reverse *x-pos*)
          do (format stream "~f, ~f ~%" tim x)))

 (let* ((target-time 1d1)
        (dt (cl-mpm:sim-dt *sim*))
        (dt-0 0d0)
        (dt-scale 0.80d0)
        (substeps (floor target-time dt))
        (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
        (load-steps 20)
        )
   (setf (cl-mpm:sim-mass-scale *sim*) 1d0)
   (setf (cl-mpm:sim-damping-factor *sim*) 0.5d0)

   (cl-mpm/output::save-simulation-parameters #p"output/settings.json"
                                              *sim*
                                              (list :dt target-time
                                                    :ocean-height *water-height*
                                                    ))

    (cl-mpm::update-sim *sim*)
    (let* ((dt-e (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
           (substeps-e (floor target-time dt-e)))
      (format t "CFL dt estimate: ~f~%" dt-e)
      (format t "CFL step count estimate: ~D~%" substeps-e)
      (setf (cl-mpm:sim-dt *sim*) dt-e)
      (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
   (setf dt-0 (/ (cl-mpm:sim-dt *sim*) (sqrt (cl-mpm::sim-mass-scale *sim*))))
   (defparameter *oobf* 0)
   (defparameter *data-damage* 0)
   (defparameter *data-energy* 0)

   (defparameter *data-full-time* (list))
   (defparameter *data-full-damage* (list))
   (defparameter *data-full-oobf* (list))
   (defparameter *data-full-energy* (list))
   (with-open-file (stream (merge-pathnames "./output/timesteps.csv") :direction :output :if-exists :supersede)
		 (format stream "steps,time,damage,energy,oobf~%"))
   (let ((energy-crit 0d0)
         (power-crit 0d0))
     (time (loop for steps from 0 to load-steps
                 while *run-sim*
                 do
                    (progn
                      (let ((g0 (* -9.82d0
                                   0.50d0
                                   ;(/ 19 50)
                                   )))
                        (loop for mp across (cl-mpm:sim-mps *sim*)
                              do (setf (cl-mpm/particle:mp-gravity mp)
                                       (+
                                        g0
                                        (* (- -9.82d0 g0) (/ (1+ steps) load-steps))))))
                      (format t "Step ~d / ~D ~%" steps load-steps)
                      (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                      (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                      (cl-mpm/output::save-vtk-cells (merge-pathnames (format nil "output/sim_cells_~5,'0d.vtk" *sim-step*)) *sim*)
                      (setf (cl-mpm::sim-enable-damage *sim*) nil)
                      (cl-mpm/dynamic-relaxation::converge-quasi-static
                           *sim*
                           :energy-crit 1d-2
                           :oobf-crit 1d-2
                           :dt-scale dt-scale
                           :conv-steps 200
                           :substeps 20
                           :post-iter-step
                           (lambda (step ke fnorm)
                             (plot *sim*)
                             (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                             (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                             ))
                      (cl-mpm/damage::calculate-damage *sim*)

                      (swank.live:update-swank)
                      (incf *sim-step*)
                     )))))
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
  )

(defun plot-disp-day ()
  (let* ((df (lisp-stat:read-csv
              (uiop:read-file-string #P"./example_data/slump/stokes.csv")))
         (time-day (mapcar (lambda (x) (/ x (* 24 60 60))) *data-full-time*))
         ;; (time-month (mapcar (lambda (x) (/ x (* 32d0 24 60 60))) *time*))
         )
    ;; (vgplot:figure)
    (vgplot:title "S_{xx} over height")
    (vgplot:xlabel "Time (d)")
    (vgplot:ylabel "Displacment (m)")
    (vgplot:axis (list 0 (reduce #'max time-day)
                       0 (reduce #'max *data-full-disp*)))
    (vgplot:plot time-day *data-full-disp* "MPM"
                 (aops:each (lambda (x) (* x 32d0)) (lisp-stat:column df 'time)) (lisp-stat:column df 'disp) "stokes"
                 )))
(defun plot-disp ()
  (let* ((df (lisp-stat:read-csv
	            (uiop:read-file-string #P"./example_data/slump/stokes.csv")))
         (time-month (mapcar (lambda (x) (/ x (* 32d0 24 60 60))) *time*))
         ;; (time-month (mapcar (lambda (x) (/ x (* 32d0 24 60 60))) *time*))
         )
    ;; (vgplot:figure)
    (vgplot:title "S_{xx} over height")
    (vgplot:xlabel "Time (s)")
    (vgplot:ylabel "Displacment (m)")
    (vgplot:axis (list 0 (reduce #'max time-month)
                       0 (reduce #'max *x-pos*)))
    (vgplot:plot time-month *x-pos* "MPM"
                 (lisp-stat:column df 'time) (lisp-stat:column df 'disp) "stokes"
                 )))
(defun plot-s-xx (mps H)
  (let* (
         (rho-ice 900)
         (E (cl-mpm/particle::mp-e (first mps)))
         (nu (cl-mpm/particle::mp-nu (first mps)))
         (g (cl-mpm/particle::mp-gravity (first mps)))
         (s-xx (loop for mp in mps
                    collect
                    (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0)))
         (y (loop for mp in mps
                    collect
                    (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)))
         (s-an (mapcar
                (lambda (y)
                  (* (/ nu (- 1 nu)) rho-ice g -1d0 (- y (/ H 2d0))))
                y)))
    (vgplot:figure)
    (vgplot:title "S_{xx} over height")
    (vgplot:xlabel "Height (m)")
    (vgplot:ylabel "Longitudinal stress (MPa)")
    (vgplot:plot y (mapcar (lambda (x) (* 1d-6 x)) s-xx) "s_{xx} mpm"
                 y (mapcar (lambda (x) (* 1d-6 x)) s-an) "s_{xx} analytic")))
(defun plot-cfl ()
  (vgplot:figure)
  (vgplot:title "CFL over time")
  (vgplot:plot *time* *cfl-max*))


(defun profile ()
  (sb-profile:unprofile)
  (sb-profile:reset)
  (sb-profile:profile cl-mpm::update-sim
                      cl-mpm::reset-grid
                      cl-mpm::p2g
                      cl-mpm::filter-grid
                      cl-mpm::update-nodes
                      cl-mpm::apply-bcs
                      cl-mpm::g2p
                      cl-mpm::update-particle
                      cl-mpm::update-stress
                      cl-mpm::iterate-over-neighbours-shape
                      cl-mpm::iterate-over-neighbours-shape-linear
                      cl-mpm::p2g-mp
                      cl-mpm::g2p-mp
                      cl-mpm::p2g-mp-node
                      cl-mpm::g2p-mp-node
                      ;; cl-mpm::update-strain-kirchoff
                      cl-mpm/damage::calculate-damage
                      cl-mpm/damage::apply-damage
                      cl-mpm/damage::delocalise-damage
                      cl-mpm/damage::create-delocalisation-list
                      ;; cl-mpm/eigenerosion:update-fracture
                      ;; cl-mpm/eigenerosion::remove-material-damaged
                      ;; cl-mpm/eigenerosion::find-neighbours
                      )
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
(defun simple-time (&optional (k 4))
  (declare (optimize (speed 3) (debug 0) (safety 0)))
  (setf lparallel:*kernel* (lparallel:make-kernel k :name "custom-kernel"))
  (setup)
  (let ((mps (cl-mpm:sim-mps *sim*) )
        (a (magicl:random-hermitian 2))
        (b (magicl:zeros '(2 2)))
        (c (magicl:zeros '(2 2)))
        (iters 10000))
    (let ((mesh (cl-mpm::sim-mesh *sim*)))
      ;; (format t "Testing normal ~%")
      ;;  (time-form iters
      ;;             (cl-mpm::update-sim *sim*))
      (format t "Eig ~%")
       (time-form iters
                  (magicl::eig a))
      (format t "Eig real ~%")
      (time-form iters
                 (magicl:eig a :realtype))
      (format t "herm Eig ~%")
      (time-form iters
                 (magicl::hermitian-eig a)
                 )
      (format t "realpart Eig ~%")
      (time-form iters
                 (multiple-value-bind (l v) (magicl::eig a) (magicl:.realpart v)))
      ;; (time
      ;; (time
      ;;  (time-form iters
      ;;             (progn
      ;;               (lparallel:pdotimes
      ;;                ;; dotimes
      ;;                (i (length mps))
      ;;                (cl-mpm::iterate-over-neighbours-point-linear
      ;;                 mesh
      ;;                 (cl-mpm/particle:mp-position (aref mps i))
      ;;                 (lambda (m n s g)))))))
      ;; (format t "Testing simd ~%")
      ;; (time
      ;;  (time-form iters
      ;;             (progn
      ;;               (lparallel:pdotimes
      ;;                (i (length mps))
      ;;                (cl-mpm::iterate-over-neighbours-point-linear-simd
      ;;                 mesh
      ;;                 (cl-mpm/particle:mp-position (aref mps i))
      ;;                 (lambda (m n s g)))))))

      ;; (magicl.backends:with-backends (:LAPACK :BLAS :LISP)
      ;;   (time-form 100
      ;;              (cl-mpm::update-sim *sim*)))
      ;; (magicl.backends:with-backends (:SIMD :LAPACK :BLAS :LISP)
      ;;   (time-form 100
      ;;              (cl-mpm::update-sim *sim*)))

      ;; (format t "Normal ~%")
      ;; (time-form iters
      ;;            (magicl:.+ a b)
      ;;            )
      ;; (format t "Blas ~%")
      ;; (time-form iters
      ;;            (magicl.blas::.+-blas a b)
      ;;            )
      ;; (format t "Simd ~%")
      ;; (time-form iters
      ;;            (magicl.simd::.+-simd a b)
      ;;            )
      )))
(defun test-backend ()
  (let ((a (magicl:zeros '(2 1)))
        (b (magicl:zeros '(2 1)))
        (c (magicl:zeros '(2 1))))
    (magicl:.+ a b c)))
(defun time-diff ()
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      *sim*
    (let* ((mp-a (aref mps 0))
           (mp-b (aref mps (- (length mps) 1)))
           (len (cl-mpm/particle::mp-true-local-length mp-a))
           )
      (time
       (dotimes (i 1000)
         do
         (cl-mpm/damage::weight-func-mps-damaged mesh mp-a mp-b len)
         ;; (cl-mpm/damage::diff-squared mp-a mp-b)
         ;; (cl-mpm::iterate-over-neighbours-point-linear mesh (cl-mpm/particle:mp-position mp-a)
                                                       ;; (lambda (&rest a)))
         ;; (cl-mpm/damage::diff-damaged mesh mp-a mp-b)
         )))))
;; (let ((stress (magicl:zeros '(3 1)))
;;      (strain (magicl:zeros '(3 1)))
;;      (de (cl-mpm/constitutive::linear-elastic-matrix 1d0 0d0))
;;      )
;;  (format t "mat")
;;  (time
;;   (dotimes (i 100000)
;;       (cl-mpm/constitutive::maxwell-exp strain stress 1d0 0d0 1d0 1d0)))
;;  (format t "voight")
;;  (time
;;   (dotimes (i 100000)
;;     (cl-mpm/constitutive::maxwell-exp-v strain stress 1d0 0d0 1d0 1d0)))
;;  (format t "simd")
;;  (time
;;   (dotimes (i 100000)
;;     (cl-mpm/constitutive::maxwell-exp-v-simd strain stress 1d0 0d0 de 1d0  1d0)))
;;  )


(defun test ()
  (let ((dt-scale 0.5d0)
        (target-time 0.1d0)
        (refine 1d0)
        (substeps 1)
        )
    (setup :refine refine :mps 2)
    (setf (cl-mpm/penalty:bc-penalty-friction *floor-bc*) 0d0)
    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (setf
              (cl-mpm/particle::mp-initiation-stress mp) 1d10
              (cl-mpm/particle::mp-enable-plasticity mp) nil))
    (setf *run-sim* t)
    (cl-mpm/output::save-simulation-parameters #p"output/settings.json"
                                               *sim*
                                               (list :dt 0.1d0
                                                     :ocean-height *water-height*
                                                     ))
    (let ((mass-scale 1d0)
          (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
      (setf
       (cl-mpm::sim-mass-scale *sim*) mass-scale
       (cl-mpm:sim-damping-factor *sim*)
       (* 0.1d0
          ;; (sqrt mass-scale)
          (cl-mpm/setup::estimate-critical-damping *sim*)))
      (setf
       (cl-mpm:sim-dt *sim*)
       (cl-mpm/setup:estimate-elastic-dt *sim* :dt-scale dt-scale))
      (setf substeps (round target-time (cl-mpm:sim-dt *sim*)))
      )

    (setf (cl-mpm::sim-enable-damage *sim*) nil)
    (loop for i from 0 to 100
          while *run-sim*
          do (progn
               (format t "Step ~D/~D ~%" i 100)
               ;; (setf (cl-mpm:sim-dt *sim*) (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
               ;; (format t "CFL dt estimate: ~f~%" (cl-mpm:sim-dt *sim*))
               (loop for j from 0 to substeps
                     while *run-sim*
                     do (cl-mpm:update-sim *sim*))
               ;; (cl-mpm/damage::calculate-damage *sim*)
               ;; (let ((dy
               ;;         (lparallel:pmap-reduce (lambda (mp)
               ;;                                  (cl-mpm/particle::mp-damage-ybar mp))
               ;;                                #'max (cl-mpm:sim-mps *sim*))
               ;;         ))
               ;;   (format t "Max damage-ybar ~E~%" dy))
               (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" i)) *sim*)
               (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" i)) *sim*)
               (cl-mpm/output::save-vtk-cells (merge-pathnames (format nil "output/sim_cells_~5,'0d.vtk" i)) *sim*)
               (plot *sim*)
               (swank.live:update-swank)
               ))
    )
  )


(defun estimate-max-stress ()
  (let ((dt-scale 0.5d0)
        (target-time 1d2)
        (refine 1d0)
        )
    (setup :refine refine :mps 2)
    (setf (cl-mpm/penalty:bc-penalty-friction *floor-bc*) 0d0)
    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (setf
              (cl-mpm/particle::mp-initiation-stress mp) 1d10
              (cl-mpm/particle::mp-enable-plasticity mp) nil))
    (setf *run-sim* t)
    (cl-mpm/output::save-simulation-parameters #p"output/settings.json"
                                               *sim*
                                               (list :dt 0.1d0
                                                     :ocean-height *water-height*
                                                     ))
    (let ((mass-scale 1d2)
          (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
      (setf
       (cl-mpm::sim-mass-scale *sim*) mass-scale
       (cl-mpm:sim-damping-factor *sim*)
       (* 0.00d0
          ;; (sqrt mass-scale)
          (cl-mpm/setup::estimate-critical-damping *sim*))))

    (setf (cl-mpm::sim-enable-damage *sim*) nil)
    (cl-mpm/dynamic-relaxation:converge-quasi-static
     *sim*
     :dt-scale dt-scale
     :energy-crit 1d-3
     :oobf-crit 1d-3
     :substeps 50
     :conv-steps 1000
     :dt-scale dt-scale
     :post-iter-step
     (lambda (i e o)
       (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" i)) *sim*)
       (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" i)) *sim*)
       (cl-mpm/output::save-vtk-cells (merge-pathnames (format nil "output/sim_cells_~5,'0d.vtk" i)) *sim*)
       (cl-mpm/damage::calculate-damage *sim*)
       (plot *sim*)
       (let ((dy
               (lparallel:pmap-reduce (lambda (mp)
                                        (cl-mpm/particle::mp-damage-ybar mp))
                                      #'max (cl-mpm:sim-mps *sim*))
               ))
         (format t "Max damage-ybar ~E~%" dy))
       (swank.live:update-swank))
    ;; (loop for i from 0 to 100
    ;;       while *run-sim*
    ;;       do (progn
    ;;            (setf (cl-mpm:sim-dt *sim*) (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
    ;;            ;; (format t "CFL dt estimate: ~f~%" (cl-mpm:sim-dt *sim*))
    ;;            (loop for j from 0 to (* refine 5)
    ;;                  while *run-sim*
    ;;                  do (cl-mpm:update-sim *sim*))
    ;;            (cl-mpm/damage::calculate-damage *sim*)
    ;;            (let ((dy
    ;;                    (lparallel:pmap-reduce (lambda (mp)
    ;;                                             (cl-mpm/particle::mp-damage-ybar mp))
    ;;                                           #'max (cl-mpm:sim-mps *sim*))
    ;;                    ))
    ;;              (format t "Max damage-ybar ~E~%" dy))
    ;;            (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" i)) *sim*)
    ;;            (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" i)) *sim*)
    ;;            (cl-mpm/output::save-vtk-cells (merge-pathnames (format nil "output/sim_cells_~5,'0d.vtk" i)) *sim*)
    ;;            (plot *sim*)
    ;;            (swank.live:update-swank)
    ;;            ))
    )))




(defun test ()
  (setup)
  (sb-profile:profile "CL-MPM")
  (sb-profile:profile "CL-MPM/PARTICLE")
  (sb-profile:profile "CL-MPM/MESH")
  (sb-profile:profile "CL-MPM/SHAPE-FUNCTION")
  (sb-profile:profile "CL-MPM/DAMAGE")
  (sb-profile:reset)
  (time
   (dotimes (i 10)
     (cl-mpm::update-sim *sim*)))
  (sb-profile:report)
  )



;; (let* ((mp (aref (cl-mpm:sim-mps *sim*) 0))
;;        (stress (cl-mpm/particle::mp-undamaged-stress mp))
;;        (de (cl-mpm/particle::mp-elastic-matrix mp))
;;        (strain (cl-mpm/particle:mp-strain mp)))
;;   (pprint stress)
;;   (pprint strain)
;;   ;(pprint (magicl:linear-solve de stress))
;;   (pprint (magicl:linear-solve de stress))
;;   )

(defun cundall-test ()
  (let* ((step 0)
         (refine 1.0)
         (work 0d0)
         (dt-scale 0.8d0)
        )

    (defparameter data-cundall-step (list))
    (defparameter data-cundall-load (list))
    (defparameter data-cundall-energy (list))
    (defparameter data-visc-step (list))
    (defparameter data-visc-load (list))
    (defparameter data-visc-energy (list))
    (setup :refine (/ 1d0 refine))
    ;; (setf (cl-mpm:sim-mass-scale *sim*) 1d0)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 100d0
             (cl-mpm:sim-mass-scale *sim*)
             (expt (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)) 2)
             ))
    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (vgplot:close-all-plots)
    (vgplot:figure)
    (setf *run-sim* t)
    (time
     (loop for step from 0 to 200
           while *run-sim*
           do
              (progn
                (format t "Step ~D~%" step)
                (dotimes (i (floor (* refine 2)))
                  (setf cl-mpm/penalty::*debug-force* 0)
                  (cl-mpm:update-sim *sim*)
                  (incf work (cl-mpm/dynamic-relaxation:estimate-power-norm *sim*))
                  )
                (format t "Step ~D - Load ~E - OOBF ~E - KE ~E~%"
                        step
                        cl-mpm/penalty::*debug-force*
                        (cl-mpm/dynamic-relaxation::estimate-oobf *sim*)
                        (abs (/ (cl-mpm/dynamic-relaxation:estimate-energy-norm *sim*) work)))
                (setf (cl-mpm:sim-dt *sim*) (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
                (push (* (cl-mpm:sim-dt *sim*) step) data-visc-step)
                ;; (push cl-mpm/penalty::*debug-force* data-visc-load)

                (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" step)) *sim*)
                (push
                 (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*)
                 data-visc-load)
                (vgplot:plot
                 data-visc-step data-visc-load "Visc")
                (swank.live:update-swank))))
    )
  ;; (plot-cundall-test)
  )

(defun stop ()
  (setf *run-sim* nil)
  (setf cl-mpm/dynamic-relaxation::*run-convergance* nil))
