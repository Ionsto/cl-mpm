(defpackage :cl-mpm/examples/shear-box
  (:use :cl))
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)
;; (asdf:compile-system :cl-mpm/examples/shear-box :full t)
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)

(in-package :cl-mpm/examples/shear-box)
;(declaim (optimize (debug 3) (safety 3) (speed 0)))
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-mc) dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-damaged mesh mp dt fbar)
  (cl-mpm::update-stress-kirchoff mesh mp dt fbar)
  ;; (cl-mpm::update-stress-linear mesh mp dt fbar)
  )
(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-chalk-delayed) dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-damaged mesh mp dt fbar)
  (cl-mpm::update-stress-kirchoff mesh mp dt fbar)
  ;; (cl-mpm::update-stress-linear mesh mp dt fbar)
  )
(defun cl-mpm/damage::length-localisation (local-length local-length-damaged damage)
  ;; (+ (* local-length (- 1d0 damage)) (* local-length-damaged damage))
  ;; (* local-length (max (sqrt (- 1d0 damage)) 1d-10))
  local-length
  )
(defmethod cl-mpm/particle::constitutive-model ((mp cl-mpm/particle::particle-elastic-damage) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-slots ((stress cl-mpm/particle::stress)
               (stress-undamaged cl-mpm/particle::undamaged-stress)
               (de cl-mpm/particle::elastic-matrix))
      mp
    (cl-mpm/constitutive::linear-elastic-mat strain de stress-undamaged)
    (cl-mpm/utils::voigt-copy-into stress-undamaged stress)
    stress))



(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-chalk-delayed) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
                     (strain cl-mpm/particle::mp-strain)
                     (plastic-strain cl-mpm/particle::mp-strain-plastic)
                     (plastic-vm cl-mpm/particle::mp-strain-plastic-vm)
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
                     (g-r cl-mpm/particle::mp-shear-residual-ratio))
        mp
      (declare (double-float pressure damage))
      (progn
        ;; (setf damage-increment (cl-mpm/damage::criterion-max-principal-stress stress))
        (setf damage-increment (cl-mpm/damage::tensile-energy-norm strain E de))
        ;; (setf damage-increment
        ;;       (cl-mpm/damage::tensile-energy-norm (cl-mpm/fastmath:fast-.+
        ;;                                            strain
        ;;                                            (cl-mpm/fastmath:fast-scale-voigt
        ;;                                             plastic-strain
        ;;                                             1d0
        ;;                                             ;; (- 1d0 damage)
        ;;                                             ))
        ;;                                           E
        ;;                                           de))

        ;; (let ((es (cl-mpm/constitutive::linear-elastic-mat
        ;;            (cl-mpm/fastmath:fast-.+
        ;;             strain
        ;;             (magicl:scale plastic-strain
        ;;                           0d0
        ;;                           ;; (- 1d0 damage)
        ;;                           )
        ;;             )
        ;;                                                    de)))
        ;;   (setf damage-increment
        ;;         (max 0d0
        ;;              (cl-mpm/damage::criterion-dp
        ;;               es
        ;;               (* angle (/ pi 180d0))))))

        ;; (let ((es (cl-mpm/constitutive::linear-elastic-mat (cl-mpm/fastmath:fast-.+ strain plastic-strain)
        ;;                                                    de)))
        ;;   (setf damage-increment
        ;;         (max 0d0
        ;;              (cl-mpm/damage::criterion-dp;cl-mpm/damage::drucker-prager-criterion
        ;;               es
        ;;                                 ;(magicl:scale es (/ 1d0 (magicl:det def)))
        ;;               (* angle (/ pi 180d0))))))

        ;; (setf damage-increment
        ;;       (max 0d0
        ;;            (cl-mpm/damage::drucker-prager-criterion
        ;;             (magicl:scale stress (/ 1d0 (magicl:det def))) (* angle (/ pi 180d0)))))

        ;; (setf damage-increment
        ;;       (+ 
        ;;        ;; (max 0d0
        ;;        ;;      (cl-mpm/damage::criterion-dp
        ;;        ;;                          ;(magicl:scale stress (/ 1d0 (magicl:det def)))
        ;;        ;;       stress
        ;;        ;;       (* angle (/ pi 180d0))))
        ;;        (* E plastic-vm)
        ;;        ))

        ;; (setf damage-increment (cl-mpm/damage::tensile-energy-norm plastic-strain E de))
        ;; (setf damage-increment (cl-mpm/damage::tensile-energy-norm
        ;;                         (cl-mpm/fastmath:fast-.+
        ;;                          strain
        ;;                          plastic-strain) E de))
        ;; (setf damage-increment 0d0)
        ;; (setf damage-increment
        ;;       (* E (cl-mpm/utils:trace-voigt plastic-strain)))
        ;;Delocalisation switch
        (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
        (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
        ))))

(declaim (notinline plot))
(defun plot (sim)
  (plot-load-disp)
  ;; (plot-conv)
  ;; (plot-domain)
  )
(defun simple-plot-contact (sim &key (plot :point) (colour-func (lambda (mp) 0d0)) (contact-bcs nil))
  (declare (function colour-func))
  "A simple GIMP plot that display only the position and size of the MPs in a sim"
  (vgplot:format-plot t "set palette defined (0 'blue', 2 'red')")
  (vgplot:format-plot t "set ticslevel 0")
  (multiple-value-bind (x y lx ly c)
    (loop for mp across (cl-mpm:sim-mps sim)
          collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
          collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
          collect (magicl:tref (cl-mpm/particle::mp-domain-size mp) 0 0) into lx
          collect (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0) into ly
          collect (funcall colour-func mp) into c
          finally (return (values x y lx ly c)))
    (let* ((points
             (reduce #'append
                     (mapcar #'cl-mpm/penalty::bc-penalty-structure-contact-points contact-bcs)
                     )
             ;; (append
             ;;        ;; (cl-mpm/penalty::bc-penalty-structure-contact-points *shear-box-struct*)
             ;;        (cl-mpm/penalty::bc-penalty-structure-contact-points *shear-box-struct-left*)
             ;;        (cl-mpm/penalty::bc-penalty-structure-contact-points *shear-box-struct-right*)
             ;;        )
             )
           (c-x (loop for p in points collect (magicl:tref p 0 0)))
           (c-y (loop for p in points collect (magicl:tref p 1 0)))
           (c-v (loop for p in points collect 1d0)))
      (vgplot:format-plot t "set cbrange [~f:~f]" (reduce #'min c) (+ 1d-5 (reduce #'max c)))
      (cond
        ((eq plot :point)
         (vgplot:plot x y c ";;with points pt 7 lc palette"))
        ((eq plot :deformed)
         (if (and contact-bcs
                  points)
             (vgplot:plot x y lx ly c ";;with ellipses lc palette"
                          c-x c-y ";;with points pt 7")
             (vgplot:plot x y lx ly c ";;with ellipses lc palette"))))))

  ;; (setf (cl-mpm/penalty::bc-penalty-structure-contact-points *shear-box-struct*) nil)
  (loop for bcs in contact-bcs
        do (setf (cl-mpm/penalty::bc-penalty-structure-contact-points bcs) nil))
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms)))
    (vgplot:format-plot t "set xrange [~f:~f]" 0d0 ms-x)
    (vgplot:format-plot t "set yrange [~f:~f]" 0d0 ms-y)
    (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
  (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim))))
    (vgplot:format-plot t "set ytics ~f" h)
    (vgplot:format-plot t "set xtics ~f" h))
  (vgplot:replot))

(declaim (notinline plot-domain))
(defun plot-domain ()
  (let ((sim *sim*))
    ;; (vgplot:plot (mapcar (lambda (x) (* x 1d3)) *data-disp*) *data-v*)
    (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
    (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
           (ms-x (first ms))
           (ms-y (second ms))
           (offset (* 0d0 (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))))
      (vgplot:format-plot t "set object 1 rect from 0,~f to ~f,~f fc rgb 'black' fs transparent solid 0.8 noborder behind"
                          (+ (* 0.5d0 *box-size*) offset)
                          *box-size*
                          ms-y)
      (vgplot:format-plot t "set object 2 rect from 0,0 to ~f,~f fc rgb 'green' fs transparent solid 0.8 noborder behind"
                          (+ (* 1d0 *box-size*) *displacement-increment*)
                          (+ (* 0.5d0 *box-size*) offset)
                          )
      (vgplot:format-plot t "set object 3 rect from ~f,0 to ~f,~f fc rgb 'green' fs transparent solid 0.8 noborder behind"
                          (+ (* 2 *box-size*) *displacement-increment*)
                          ms-x
                          (+ (* 0.5d0 *box-size*) offset))
      (vgplot:format-plot t "set object 4 rect from ~f,~f to ~f,~f fc rgb 'black' fs transparent solid 0.8 noborder behind"
                          (* 2d0 *box-size*)
                          (+ (* 0.5d0 *box-size*) offset)
                          ms-x
                          ms-y))
    (vgplot:format-plot t "set style fill solid")
    ;; (cl-mpm/plotter:simple-plot
    ;;  *sim*
    ;;  :plot :deformed
    ;;  :colour-func
    ;;  (lambda (mp)
    ;;    (if (= 0 (cl-mpm/particle::mp-index mp))
    ;;        (cl-mpm/particle::mp-strain-plastic-vm mp)
    ;;        0d0))
    ;;  )
    (simple-plot-contact
     *sim*
     :plot :deformed
     :colour-func
     ;; #'cl-mpm/particle::mp-damage
     (lambda (mp ) (magicl:tref (cl-mpm/particle::mp-stress mp) 0 0))
     ;; (lambda (mp)
     ;;   (if (= 0 (cl-mpm/particle::mp-index mp))
     ;;       (cl-mpm/particle::mp-strain-plastic-vm mp)
     ;;       0d0))
     :contact-bcs
     (list
      ;; *shear-box-struct-left*
      ;; *shear-box-struct-right*
      ;; *shear-box-struct-left-static*
      ;; *shear-box-struct-right-dynamic*
      ;; *shear-box-struct-floor*
      *shear-box-struct*
      )
     )
    ;; (let* ((points (cl-mpm/penalty::bc-penalty-structure-contact-points *shear-box-struct-left*))
    ;;        (c-x (loop for p in points collect (magicl:tref p 0 0)))
    ;;        (c-y (loop for p in points collect (magicl:tref p 1 0)))
    ;;        (c-v (loop for p in points collect 1d0))
    ;;        )
    ;;   (when points
    ;;     (vgplot:plot
    ;;      c-x c-y c-v ";;with points pt 7"
    ;;      )))
    ))
(declaim (notinline plot-load-disp))
(defun plot-load-disp ()

  ;; (vgplot:figure)
  (let* ((width 0.06d0)
        (depth ;(cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))
          1d0
               )
        (max-load (reduce #'max (mapcar (lambda (x) (* depth (/ x width) 1d-3)) *data-v*)))
        )
    (vgplot:plot
     (mapcar (lambda (x) (* x 1d3)) *data-disp*)
     (mapcar (lambda (x) (* depth (/ x width) 1d-3)) *data-v*)
     "Load"
     (mapcar (lambda (x) (* x 1d3)) *data-disp*)
     (mapcar (lambda (x) (* x 5d-1)) *data-damage*)
     "Damage"
     (mapcar (lambda (x) (* x 1d3)) *data-disp*)
     (mapcar (lambda (x) (* x 50d0)) *data-plastic*)
     "Plastic"
     ))
  (vgplot:xlabel "Displacement (mm)")
  (vgplot:ylabel "Shear stress (kN/m^2)")
  )

(defparameter *enable-box-friction* t)
(defparameter *displacement-increment* 0d0)
(defparameter *box-size* 0d0)


(declaim (notinline setup-test-column))
(defun setup-test-column (size offset block-size &optional (e-scale 1) (mp-scale 1) &key (angle 0d0) (friction 0.1d0) (surcharge-load 72.5d3))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               ;; :sim-type 'cl-mpm::mpm-sim-usf
               :sim-type 'cl-mpm/damage::mpm-sim-damage
               ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         ;; (floor-offset (* h-y 2))
         (floor-offset 2d0)
         (density 1.7d3)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
         )
    (declare (double-float h density))
    (progn
      (let* ((angle-rad (* angle (/ pi 180)))
             ;; (init-stress 60d3)
             ;; (init-stress 100d3)
                                        ;(init-stress 100d3)
             (init-stress (* 1.0 131d3))
                                        ;(gf 48d0)
             ;; (gf 48d0)
             (gf 1d0)
             ;; (gf 100d0)
             (length-scale (* 1 h))
             ;; (length-scale 0.015d0)
             (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress 1d9))
             (ductility 5d0)
             ;; (ductility 1d8)
             ;; (ductility 1d1)
             ;; (ductility 1d4)
             )
        (format t "Estimated ductility ~E~%" ductility)
        (cl-mpm::add-mps
         sim
         (cl-mpm/setup::make-mps-from-list
          (cl-mpm/setup::make-block-mps-list
           offset
           block-size
           (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
           density

           ;; 'cl-mpm/particle::particle-vm
           ;; :E 1d9
           ;; :nu 0.24d0
           ;; :rho 100d3

           ;; 'cl-mpm/particle::particle-mc
           ;; ;; 'cl-mpm/particle::particle-dp
           ;; :E 1d9
           ;; :nu 0.24d0
           ;; ;; :psi 0d0
           ;; :psi (* 42d0 (/ pi 180))
           ;; :phi (* 42d0 (/ pi 180))
           ;; :c 131d3
           ;; :phi-r (* 30d0 (/ pi 180))
           ;; :c-r 0d0
           ;; :softening 0d0

           'cl-mpm/particle::particle-chalk-delayed;-grassl
           :E 1d9
           :nu 0.24d0
           :kt-res-ratio 1d-9
           :kc-res-ratio 1d0
           :g-res-ratio 1d-1
           ;; :g-res-ratio 1d-9
           ;; :damage 0.9d0
           :friction-angle 42d0
           :initiation-stress init-stress;18d3
           :delay-time 1d-3
           :delay-exponent 1d0
           :damage 0d0
           :ductility ductility
           :local-length length-scale
           :local-length-damaged 10d-10
           :enable-plasticity t
           :psi 0d0
           ;; :psi (* 42d0 (/ pi 180))
           :phi (* 60d0 (/ pi 180))
           :c 131d3

           ;; :phi (* 50d0 (/ pi 180))
           ;; :c (* 131d3 100d0)

           :index 0
           :gravity 0d0
           )))
        )
      (let* (;; (sur-height h-x)
             (sur-height (* 0.5 (second block-size)))
             (sur-size (list 0.06d0 sur-height))
                                        ;(load 72.5d3)
             (load surcharge-load)
             ;; (gravity 10d0)
             ;; (density (/ load (* gravity sur-height)))
             (gravity (if (> sur-height 0d0)
                          (/ load (* density sur-height))
                          0d0))
             (mp-surcharge t)
             )
        (format t "Gravity ~F~%" gravity)
        (if mp-surcharge
            (cl-mpm::add-mps
             sim
             (cl-mpm/setup::make-mps-from-list
              (cl-mpm/setup::make-block-mps-list
               (mapcar #'+ offset (list 0d0 (second block-size)))
               sur-size
               (mapcar (lambda (e) (* e e-scale 2)) sur-size)
               density
               'cl-mpm/particle::particle-elastic-damage
               :E 1d9
               :nu 0.24d0
               :initiation-stress 1d20
               :local-length 0d0
               :index 1
               :gravity (- gravity))))
            (progn
              (defparameter *pressure-bc*
                (cl-mpm/buoyancy::make-bc-pressure
                 sim
                 0d0
                 (- surcharge-load)
                 :clip-func
                 (lambda (pos)
                   (and
                    (> (cl-mpm/utils:varef pos 1)
                       (+ (second offset) (* 0.5d0 (second block-size))))
                    ;; (> (cl-mpm/utils:varef pos 0)
                    ;;    (first offset))
                    ;; (< (cl-mpm/utils:varef pos 0)
                    ;;    (+ (first offset) (first block-size)))
                    )

                   )))
              (cl-mpm:add-bcs-force-list
               sim
               *pressure-bc*)))
        )

      (let ((mp-0 (aref (cl-mpm:sim-mps sim) 0)))
        (when (typep mp-0 'cl-mpm/particle::particle-chalk-brittle)
          (let* (
                 (fc (cl-mpm/particle::mp-fc mp-0))
                 (ft (cl-mpm/particle::mp-ft mp-0))
                 (rc (cl-mpm/particle::mp-k-compressive-residual-ratio mp-0))
                 (rs (cl-mpm/particle::mp-shear-residual-ratio mp-0))
                 (angle-plastic (cl-mpm/particle::mp-phi mp-0))
                 (angle-plastic-damaged (atan (* (/ rs rc) (tan angle-plastic))))
                 )
            (format t "Chalk plastic virgin angle: ~F~%"
                    (* (/ 180 pi) angle-plastic))
            (format t "Chalk plastic residual angle: ~F~%"
                    (* (/ 180 pi) angle-plastic-damaged)))))


      (defparameter *mesh-resolution* h-x)
      (setf (cl-mpm:sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-enable-fbar sim) t)
      ;; (setf (cl-mpm::sim-mass-filter sim) 1d0)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (let ((ms 1d0))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        (setf (cl-mpm:sim-damping-factor sim)
              (* 0d0 (cl-mpm/setup::estimate-critical-damping sim))))

      (let ((dt-scale 1d0))
        (setf
         (cl-mpm:sim-dt sim)
         (* dt-scale h
            (sqrt (cl-mpm::sim-mass-scale sim))
            (sqrt (/ density (cl-mpm/particle::mp-p-modulus (aref (cl-mpm:sim-mps sim) 0)))))))

      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))
      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-varfix
             (cl-mpm:sim-mesh sim)
             '(0 nil nil)
             '(0 nil nil)
             '(nil 0 nil)
             '(nil 0 nil)
             '(nil nil 0)
             '(nil nil 0)))
      ;; (cl-mpm:iterate-over-mps
      ;;  (cl-mpm:sim-mps sim)
      ;;       do
      ;;       (let* ((pressure surcharge-load)
      ;;              (E 1d9)
      ;;              (nu 0.24d0)
      ;;              (strain (cl-mpm/utils:voigt-from-list (list 0d0 (- (/ pressure E))
      ;;                                                          0d0
      ;;                                                          0d0
      ;;                                                          0d0
      ;;                                                          0d0)))
      ;;              (stress (cl-mpm))
      ;;              )
      ;;         )
      ;;         (setf )
      ;;       )
      sim)))


(defparameter *sim* nil)
(defparameter *run-sim* t)
(defparameter *t* 0)
(defparameter *sim-step* 0)
(defparameter *refine* (/ 1d0 1d0))

(let ((refine (uiop:getenv "REFINE")))
  (when refine
    (setf *refine* (parse-integer (uiop:getenv "REFINE")))))

(defmethod cl-mpm::post-stress-step (mesh (mp cl-mpm/particle::particle-vm) dt)
  ;; (with-accessors ((ps cl-mpm/particle::mp-strain-plastic-vm)
  ;;                  (rho cl-mpm/particle::mp-rho))
  ;;     mp
  ;;   (let* ((rho-0 200d3)
  ;;          (rho-1 200d2)
  ;;          (soft 100d0)
  ;;          )
  ;;     (setf rho (+ rho-1 (* (- rho-0 rho-1) (exp (- (* soft ps))))))))
  )

(defmethod cl-mpm::post-stress-step (mesh (mp cl-mpm/particle::particle-mc) dt)
  ;; (with-accessors ((ps cl-mpm/particle::mp-strain-plastic-vm)
  ;;                  (c cl-mpm/particle::mp-c)
  ;;                  (phi cl-mpm/particle::mp-phi)
  ;;                  )
  ;;     mp
  ;;   (let ((phi_0 (* 42d0 (/ pi 180)))
  ;;         (phi_1 (* 30d0 (/ pi 180)))
  ;;         (c_0 131d3)
  ;;         (soft 1000d0;(* 100d0 *mesh-resolution*)
  ;;               ))
  ;;     (setf
  ;;      c (* c_0 (exp (- (* soft ps))))
  ;;      phi (+ phi_1 (* (- phi_0 phi_1) (exp (- (* soft ps))))))))
  )
(defmethod cl-mpm::post-stress-step (mesh (mp cl-mpm/particle::particle-chalk-delayed) dt))

(declaim (notinline make-penalty-box))
(defun make-penalty-box (sim left-x right-x height friction offset)
  (let* ((left-normal (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0)))
         (right-normal (cl-mpm/utils:vector-from-list (list -1d0 0d0 0d0)))
         (plane-normal (cl-mpm/utils:vector-from-list (list 0d0 -1d0 0d0)))
         (plane-normal-left (cl-mpm/utils:vector-from-list (list 0d0 1d0 0d0)))
         (epsilon (* 1d2 1d9))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (extra-height 0d0)
         (friction (tan (* 30d0 (/ pi 180))))
         ;; (friction 0.0d0)
         (damping 0d0))
    (defparameter *shear-box-left-static*
      (cl-mpm/penalty::make-bc-penalty-distance-point
       sim
       left-normal
       (cl-mpm/utils:vector-from-list (list left-x (+ (* 3d0 height) offset) 0d0))
       (* 2d0 height)
       epsilon
       0d0
       damping))
    (defparameter *shear-box-right-static*
      (cl-mpm/penalty::make-bc-penalty-distance-point
       sim
       right-normal
       (cl-mpm/utils:vector-from-list (list right-x (+ (* 3d0 height) offset) 0d0))
       (* 2d0 height)
       epsilon
       0d0
       damping))
    (defparameter *shear-box-left-slide*
      (cl-mpm/penalty::make-bc-penalty-distance-point
       sim
       plane-normal-left
       (cl-mpm/utils:vector-from-list (list (- left-x height) (+ height offset) 0d0))
       (* 1d0 height)
       epsilon
       0d0
       damping))
    (defparameter *shear-box-left-dynamic*
      (cl-mpm/penalty::make-bc-penalty-distance-point
       sim
       left-normal
       (cl-mpm/utils:vector-from-list (list left-x
                                            (+ (- (* 0.5d0 height)
                                                  extra-height)
                                               offset) 0d0))
       (+
        (* 0.5d0 height)
        (* 1 extra-height)
        )
       epsilon
       0d0
       damping))
    (defparameter *shear-box-right-dynamic*
      (cl-mpm/penalty::make-bc-penalty-distance-point
       sim
       right-normal
       (cl-mpm/utils:vector-from-list (list right-x
                                            (+
                                             (- (* 0.5d0 height)
                                                extra-height)
                                             offset) 0d0))
       (+
        (* 0.5d0 height)
        (* 1 extra-height)
        )
       epsilon
       0d0
       damping))
    (defparameter *shear-box-right-slide*
      (cl-mpm/penalty::make-bc-penalty-distance-point
       sim
       plane-normal
       (cl-mpm/utils:vector-from-list (list (+ right-x height) (+ height offset) 0d0))
       (* 1d0 height)
       epsilon
       0d0
       damping))
    (defparameter *shear-box-floor*
      (cl-mpm/penalty::make-bc-penalty-point-normal
       sim
       plane-normal-left
       (cl-mpm/utils:vector-from-list (list 0d0 offset 0d0))
       epsilon
       0d0
       damping))

    (defparameter *shear-box-struct-left*
      (cl-mpm/penalty::make-bc-penalty-structure
       sim
       epsilon
       0d0
       damping
       (list
        *shear-box-left-static*
        *shear-box-left-dynamic*
        *shear-box-left-slide*
        ;; *shear-box-floor*
        )))

    (defparameter *shear-box-struct-right*
      (cl-mpm/penalty::make-bc-penalty-structure
       sim
       epsilon
       0d0
       damping
       (list
        *shear-box-right-static*
        *shear-box-right-dynamic*
        *shear-box-right-slide*
        ;; *shear-box-floor*
        )))
    (defparameter *shear-box-struct-left-static*
      (cl-mpm/penalty::make-bc-penalty-structure
       sim
       epsilon
       0d0
       damping
       (list
        *shear-box-left-static*)))
    (defparameter *shear-box-struct-right-dynamic*
      (cl-mpm/penalty::make-bc-penalty-structure
       sim
       epsilon
       0d0
       damping
       (list
        *shear-box-right-dynamic*)))
    (defparameter *shear-box-struct-floor*
      (cl-mpm/penalty::make-bc-penalty-structure
       sim
       epsilon
       0d0
       damping
       (list
        *shear-box-floor*)))
    (defparameter *shear-box-struct*
      (cl-mpm/penalty::make-bc-penalty-structure
       sim
       epsilon
       0d0
       damping
       (list
        *shear-box-right-static*
        *shear-box-right-dynamic*
        *shear-box-right-slide*
        *shear-box-left-static*
        *shear-box-left-dynamic*
        *shear-box-left-slide*
        ;; *shear-box-floor*
        )))

    (defparameter *shear-box-controller*
      (cl-mpm/bc::make-bc-closure
       nil
       (lambda ()

         (setf (cl-mpm/penalty::bc-penalty-load *shear-box-left-dynamic*) 0d0)
         (setf (cl-mpm/penalty::bc-penalty-load *shear-box-right-dynamic*) 0d0)

         (setf (cl-mpm/penalty::bc-penalty-datum *shear-box-left-dynamic*)
               (+ left-x *displacement-increment*))

         (setf (magicl:tref (cl-mpm/penalty::bc-penalty-distance-center-point *shear-box-left-slide*) 0 0)
               (+ (- left-x height) *displacement-increment*))

         (setf (cl-mpm/penalty::bc-penalty-datum *shear-box-right-dynamic*)
               (- (+ right-x *displacement-increment*)))
         (let ((friction-bcs (list
                              *shear-box-right-static*
                              *shear-box-right-dynamic*
                              *shear-box-left-static*
                              *shear-box-left-dynamic*
                              )))
           (loop for bc in friction-bcs
                 do (setf (cl-mpm/penalty::bc-penalty-friction bc)
                          (if *enable-box-friction*
                              friction
                              0d0)))));)
       )
      )
    (cl-mpm:add-bcs-force-list
     *sim*
     (list
      (cl-mpm/bc:make-bcs-from-list
       (list
        *shear-box-controller*))
      (cl-mpm/bc:make-bcs-from-list
       (list
        ;; *shear-box-struct-left-static*
        ;; *shear-box-struct-left*
        ;; *shear-box-struct-right*
        ;; *shear-box-struct-right-dynamic*
        ;; *shear-box-struct-floor*
        *shear-box-struct*
        ))))))

(defun get-load ()
  (cl-mpm/penalty::bc-penalty-load *true-load-bc*)
  ;; (-
  ;;  (cl-mpm/penalty::bc-penalty-load *shear-box-left-dynamic*)
  ;;  (cl-mpm/penalty::bc-penalty-load *shear-box-right-dynamic*))
  )

(defun get-damage ()
  (lparallel:pmap-reduce
   (lambda (mp)
     (if (= (cl-mpm/particle::mp-index mp) 0)
         (*
          1d2
          (cl-mpm/particle::mp-mass mp)
          (cl-mpm/particle:mp-damage mp))
         0d0))
   #'+
   (cl-mpm:sim-mps *sim*)))
(defun get-plastic ()
  (lparallel:pmap-reduce
   (lambda (mp)
     (if (= (cl-mpm/particle::mp-index mp) 0)
         (cl-mpm/particle::mp-strain-plastic-vm mp)
         0d0))
   #'+
   (cl-mpm:sim-mps *sim*)))

(defun setup (&key (refine 1d0) (mps 2) (friction 0.0d0) (surcharge-load 72.5d3))
  (defparameter *displacement-increment* 0d0)
  (let* ((mps-per-dim mps)
         (mesh-size (/ 0.03d0 refine))
         (sunk-size 0.03d0)
         (box-size (* 2d0 sunk-size))
         (domain-size (* 3d0 box-size))
         (box-offset (* mesh-size 0d0))
         (offset (list box-size box-offset)))
    (setf *box-size* box-size)
    (defparameter *sim* (setup-test-column
                         (list domain-size (+ (* 2 box-size) box-offset))
                         offset
                         (list box-size box-size)
                         (/ 1d0 mesh-size)
                         mps-per-dim
                         :friction friction
                         :surcharge-load surcharge-load))
    (make-penalty-box *sim* box-size (* 2d0 box-size) sunk-size friction box-offset)
    ;; (cl-mpm/setup:remove-sdf
    ;;  *sim*
    ;;  (lambda (p)
    ;;    (if (> (- (cl-mpm/utils:varef p 1) box-offset) sunk-size)
    ;;        0d0
    ;;        1d0)))
    (defparameter *true-load-bc* *shear-box-left-dynamic*)
    ;; (setf (cl-mpm::sim-bcs-force-list *sim*)
    ;;       (list
    ;;        (cl-mpm/bc:make-bcs-from-list
    ;;         (list
    ;;          (cl-mpm/bc::make-bc-closure
    ;;           nil
    ;;           (lambda ()
    ;;             (apply-penalty-box box-size (* 2d0 box-size) sunk-size friction)))))))
    )
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (format t "Mesh-size: ~E~%" (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0)
  
  )

(defun stop ()
  (setf *run-sim* nil
        cl-mpm/dynamic-relaxation::*run-convergance* nil))

(declaim (notinline run))
(defun run (&optional (output-directory "./output/") &key (total-time 6d-2) (damping 1d0))
  (format t "Output dir ~A~%" output-directory)
  (ensure-directories-exist (merge-pathnames output-directory))
  (cl-mpm/output:save-vtk-mesh (merge-pathnames output-directory "mesh.vtk") *sim*)
  (cl-mpm/output::save-simulation-parameters
   (merge-pathnames output-directory "settings.json")
   *sim*)

  (defparameter *data-t* (list))
  (defparameter *data-disp* (list))
  (defparameter *data-v* (list))
  (defparameter *data-damage* (list))
  (defparameter *data-plastic* (list))
  (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :supersede)
    (format stream "disp,load~%"))
  (vgplot:close-all-plots)
  (let* ((displacment 1d-3)
         ;(total-time (* 50d0 displacment))
         (time-per-mm 100d0)
         (total-time (* time-per-mm displacment))
         (load-steps (round (* 100 (/ displacment 1d-3))))
         (target-time (/ total-time load-steps))
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 0.5d0)
         (enable-plasticity (cl-mpm/particle::mp-enable-plasticity (aref (cl-mpm:sim-mps *sim*) 0)))
         (disp-inc (/ displacment load-steps)))
    ;;Disp rate in test 4d-4mm/s -> 4d-7mm/s
    (format t "Loading rate: ~E~%" (/ displacment (* load-steps target-time)))
    (format t "Total time: ~E~%" total-time)
    (format t "Target time: ~E~%" target-time)
    (format t "Strain rate: ~E~%" (/ displacment total-time))

    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (when (typep mp 'cl-mpm/particle::particle-damage)
              (when (= (cl-mpm/particle::mp-index mp) 0)
                (setf (cl-mpm/particle::mp-delay-time mp) (* (* time-per-mm 1d-3) 1d-3)))))
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 0.05d0
             (sqrt (cl-mpm:sim-mass-scale *sim*))
             (cl-mpm/setup::estimate-critical-damping *sim*)))

    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (when (= (cl-mpm/particle::mp-index mp) 0)
               (setf (cl-mpm/particle::mp-enable-plasticity mp) nil)))
    ;;Turn off friction for settling
    (setf *enable-box-friction* nil)
    (setf (cl-mpm::sim-enable-damage *sim*) nil)

    (vgplot:figure)
    (cl-mpm/output:save-vtk (merge-pathnames output-directory (format nil "sim_conv_~5,'0d.vtk" 0)) *sim*)
    (cl-mpm/dynamic-relaxation:converge-quasi-static
     *sim*
     :energy-crit 1d-2
     :oobf-crit 1d-2
     :dt-scale 0.50d0
     :substeps 10
     :conv-steps 500
     :post-iter-step
     (lambda (i energy oobf)
       (plot-domain)
       (cl-mpm/output:save-vtk (merge-pathnames output-directory (format nil "sim_conv_~5,'0d.vtk" (+ 1 i))) *sim*)
       (cl-mpm/output:save-vtk-nodes (merge-pathnames output-directory (format nil "sim_nodes_conv_~5,'0d.vtk" (+ 1 i))) *sim*)
       ))
    (vgplot:figure)

    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (when (= (cl-mpm/particle::mp-index mp) 0)
               (setf (cl-mpm/particle::mp-enable-plasticity mp) enable-plasticity)))

    (setf *enable-box-friction* t)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (*
           1d-3
           ;; (sqrt (cl-mpm:sim-mass-scale *sim*))
           (cl-mpm/setup::estimate-critical-damping *sim*)))

    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))

    (setf substeps (round target-time (cl-mpm:sim-dt *sim*)))

    (when (slot-exists-p *sim* 'cl-mpm/damage::delocal-counter-max)
      (setf (cl-mpm/damage::sim-damage-delocal-counter-max *sim*) substeps))
    (defparameter *displacement-increment* 0d0)
    (format t "Substeps ~D~%" substeps)
    (vgplot:figure)
    (cl-mpm:update-sim *sim*)
    (let ((disp-av 0d0)
          (load-av 0d0)
          (d-av 0d0)
          (p-av 0d0))
      (push *t* *data-t*)
      (push disp-av *data-disp*)
      (push load-av *data-v*)
      (push d-av *data-damage*)
      (push p-av *data-plastic*)
      (setf load-av (get-load))
      (setf disp-av *displacement-increment*)
      (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :append)
        (format stream "~f,~f~%" disp-av load-av)))

    (setf (cl-mpm::sim-enable-damage *sim*) t)

    (setf cl-mpm/penalty::*debug-force* 0)
    (time (loop for steps from 0 below load-steps
                while *run-sim*
                do
                   (progn
                     (format t "Step ~d/~D~%" steps load-steps)
                     (cl-mpm/output:save-vtk (merge-pathnames output-directory (format nil "sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     ;; (cl-mpm/output::save-vtk-nodes (merge-pathnames output-directory (format nil "sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     (let ((load-av 0d0)
                           (disp-av 0d0)
                           (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
                       (time
                        (dotimes (i substeps)
                          (cl-mpm::update-sim *sim*)
                          (incf load-av (/ (get-load) substeps))
                          (incf disp-av (/ *displacement-increment* substeps))
                          (incf *displacement-increment* (/ disp-inc substeps))
                          (incf *t* (cl-mpm::sim-dt *sim*))))

                       ;; (setf load-av (get-load))
                       ;; (setf disp-av *displacement-increment*)

                       (push *t* *data-t*)
                       (push disp-av *data-disp*)
                       (push load-av *data-v*)
                       (push (get-damage) *data-damage*)
                       (push (get-plastic) *data-plastic*)
                       (format t "Disp ~E - Load ~E~%" disp-av load-av)
                       (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :append)
                         (format stream "~f,~f~%" disp-av load-av)))
                     (incf *sim-step*)
                     (plot *sim*)
                     (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                        :terminal "png size 1920,1080")
                     (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       (format t "CFL dt estimate: ~f~%" dt-e)
                       (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf substeps substeps-e))
                     (swank.live:update-swank))))

    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (progn
               (cl-mpm/fastmath::fast-zero (cl-mpm/particle::mp-acceleration mp))
               (cl-mpm/fastmath::fast-zero (cl-mpm/particle:mp-velocity mp))))
    (when nil
      (defparameter *true-load-bc* *shear-box-right-dynamic*)
      (time (loop for steps from 0 below load-steps
                  while *run-sim*
                  do
                     (progn
                       (format t "Step ~d ~%" steps)
                       (cl-mpm/output:save-vtk (merge-pathnames output-directory (format nil "sim_~5,'0d.vtk" *sim-step*)) *sim*)
                       (cl-mpm/output::save-vtk-nodes (merge-pathnames output-directory (format nil "sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                       (let ((load-av 0d0)
                             (disp-av 0d0)
                             (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
                         (time
                          (dotimes (i substeps)
                            (cl-mpm::update-sim *sim*)
                            (incf load-av (/ (get-load) substeps))
                            (incf disp-av (/ *displacement-increment* substeps))
                            (incf *displacement-increment* (/ (- disp-inc) substeps))
                            (incf *t* (cl-mpm::sim-dt *sim*))))

                         ;; (setf load-av cl-mpm/penalty::*debug-force*)
                         ;; (setf load-av (get-load))
                         ;; (setf disp-av *displacement-increment*)

                         (push *t* *data-t*)
                         (push disp-av *data-disp*)
                         (push load-av *data-v*)
                         (format t "Disp ~E - Load ~E~%" disp-av load-av)
                         (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :append)
                           (format stream "~f,~f~%" disp-av load-av)))
                       (incf *sim-step*)
                       (plot *sim*)
                       (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                          :terminal "png size 1920,1080")
                       (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                         (format t "CFL dt estimate: ~f~%" dt-e)
                         (format t "CFL step count estimate: ~D~%" substeps-e)
                         (setf substeps substeps-e))
                       (swank.live:update-swank)))))
    )
  (vgplot:figure)
  (plot-load-disp))


(defun run-static (&optional (output-directory "./output/"))
  (format t "Output dir ~A~%" output-directory)
  (ensure-directories-exist (merge-pathnames output-directory))
  (cl-mpm/output:save-vtk-mesh (merge-pathnames output-directory "mesh.vtk")
                          *sim*)
  (cl-mpm/output::save-simulation-parameters
   (merge-pathnames output-directory "settings.json")
   *sim*)

  (defparameter *data-t* (list))
  (defparameter *data-disp* (list))
  (defparameter *data-v* (list))
  (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :supersede)
    (format stream "disp,load~%"))
  (vgplot:close-all-plots)
  (let* ((load-steps 100)
         (dt (cl-mpm:sim-dt *sim*))
         (dt-scale 0.5d0)
         (displacment 1d-3)
         (enable-plasticity (cl-mpm/particle::mp-enable-plasticity (aref (cl-mpm:sim-mps *sim*) 0)))
         (disp-inc (/ displacment load-steps)))
    (loop for mp across (cl-mpm:sim-mps *sim*)
          when (typep mp 'cl-mpm/particle::particle-chalk-delayed)
          do (change-class mp 'cl-mpm/particle::particle-chalk-brittle))

    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 0.05d0
             ;; (sqrt (cl-mpm:sim-mass-scale *sim*))
             (cl-mpm/setup::estimate-critical-damping *sim*)))

    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (when (= (cl-mpm/particle::mp-index mp) 0)
               (setf (cl-mpm/particle::mp-enable-plasticity mp) enable-plasticity)))
    (vgplot:figure)
    (setf *enable-box-friction* nil)
    (cl-mpm/dynamic-relaxation:converge-quasi-static
     *sim*
     :energy-crit 1d-2
     :oobf-crit 1d-2
     :substeps 10
     :conv-steps 200
     :dt-scale dt-scale
     :post-iter-step
     (lambda (i energy oobf)
       (plot-domain)
       (cl-mpm/output:save-vtk (merge-pathnames output-directory (format nil "sim_conv_~5,'0d.vtk" (+ 1 i))) *sim*)
       )
     )
    (vgplot:figure)
    ;; (setf *enable-box-friction* t)
    (defparameter *displacement-increment* 0d0)
    (vgplot:figure)
    (setf cl-mpm/penalty::*debug-force* 0)

    (setf load-av (get-load))
    (setf disp-av *displacement-increment*)

    (push *t* *data-t*)
    (push disp-av *data-disp*)
    (push load-av *data-v*)
    (time (loop for steps from 0 below load-steps
                while *run-sim*
                do
                   (progn
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames output-directory (format nil "sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-nodes (merge-pathnames output-directory (format nil "sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     (let ((load-av 0d0)
                           (disp-av 0d0))
                       (incf *displacement-increment* disp-inc)
                       (let ((damage-0 0d0)
                             (damage-1 0d0))
                         (loop for i from 0 to 100
                               while (or (= i 0) (> damage-1 damage-0))
                               do
                                  (progn
                                    (setf damage-0 damage-1)
                                    (format t "Staggered solve: ~D~%" i)
                                    (cl-mpm/dynamic-relaxation:converge-quasi-static
                                     *sim*
                                     :energy-crit 1d-2
                                     :oobf-crit 1d-2
                                     :substeps 10
                                     :conv-steps 200
                                     :dt-scale dt-scale
                                     :post-iter-step
                                     (lambda (i energy oobf)))
                                    (cl-mpm/damage::calculate-damage *sim*)
                                    ;; (plot *sim*)
                                    (setf damage-1 (get-damage))
                                    (format t "Staggered damage diff: ~E~%" (- damage-1 damage-0))
                                    ))
                         (when (> damage-1 damage-0)
                           (break "Staggered error"))
                         )

                       (setf load-av (get-load))
                       (setf disp-av *displacement-increment*)

                       (push *t* *data-t*)
                       (push disp-av *data-disp*)
                       (push load-av *data-v*)
                       (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :append)
                         (format stream "~f,~f~%" disp-av load-av)))
                     (incf *sim-step*)
                     (plot *sim*)
                     (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                        :terminal "png size 1920,1080"
                                        )
                     (swank.live:update-swank))))
    )
  (vgplot:figure)
  (plot-load-disp))

(defun test-conv ()
  (setf *run-sim* t)
  (loop for refine in ;(list 1 2 4 8)
        (list 1 5 10 50)
        while *run-sim*
        do (progn
             (format t "Refine: ~D~%" refine)
             (setup
              :refine 4
              :mps 3
              ;:refine refine :mps 2
              ;; :refine 1 :mps (+ 1 refine)
                    )
             (run (format nil "output-~D/" refine)
                  :total-time (* 1d-2 refine)
                  )
             ;; (run-static (format nil "output-~D/" refine))
             ))
  (vgplot:close-all-plots)
  (plot-conv))

(defun test-time ()
  (setf *run-sim* t)
  (loop for refine in ;(list 1 2 4 8)
                   (list 1 5 10 50 100)
        while *run-sim*
        do (progn
             (format t "Refine: ~D~%" refine)
             (setup
              :refine 2
                                        ;:refine refine :mps 2
              ;; :refine 1 :mps (+ 1 refine)
              )
             (run (format nil "output-~D/" refine)
                  :total-time (* 1d-2 refine)
                  )
             ;; (run-static (format nil "output-~D/" refine))
             ))
  (vgplot:close-all-plots)
  (plot-conv))

(defun test-damping ()
  (setf *run-sim* t)
  (loop for refine in ;(list 1 2 4 8)
                   (list 1d0 0.1d0 0.01d0)
        while *run-sim*
        do (progn
             (format t "Refine: ~D~%" refine)
             (setup
              :refine 2
                                        ;:refine refine :mps 2
              ;; :refine 1 :mps (+ 1 refine)
              )
             (run (format nil "output-~D/" refine)
                  :damping refine
                  )
             ;; (run-static (format nil "output-~D/" refine))
             ))
  (vgplot:close-all-plots)
  (plot-conv))

(defun test-frictional ()
  (setf *run-sim* t)
  (loop for surcharge
          in (list 91d3 153d3 277d3)
        while *run-sim*
        do (progn
             (setup
              :refine 2
              :surcharge-load surcharge
                                        ;:refine refine :mps 2
              ;; :refine 1 :mps (+ 1 refine)
              )
             (run (format nil "output-~D/" surcharge))
             ;; (run-static (format nil "output-~D/" refine))
             ))
  (vgplot:close-all-plots)
  (plot-conv))


(defun plot-conv ()
  (let* ((folders (uiop:subdirectories "./"))
         (output-folders (loop for f in folders when (search "output-" (namestring f)) collect f))
         (name-list (list))
         (disp-list (list))
         (force-list (list)))
    (loop for f in output-folders
          do (progn
               (let ((df (lisp-stat:read-csv (uiop:merge-pathnames* f "disp.csv"))))
                 (push (first (last (pathname-directory f))) name-list)
                 (push (loop for x across (lisp-stat:column df 'load) collect x) force-list)
                 (push (loop for x across (lisp-stat:column df 'disp) collect x) disp-list)
                 )
               ))
    (apply #'vgplot:plot (reduce #'append (mapcar #'list
                                                    (mapcar (lambda (l) (mapcar (lambda (x) (* x 1d3)) l)) disp-list)
                                                    (mapcar (lambda (l) (mapcar (lambda (x) (* x (/ 1d-3 6d-2))) l)) force-list)
                                                    (mapcar (lambda (x) (format nil "~A" x)) name-list)
                                                    ))))
  (vgplot:xlabel "Displacement (mm)")
  (vgplot:ylabel "Shear stress (kN/m^2)")
  )


(defun profile ()
  (declare (optimize (speed 3) (debug 0) (safety 0)))
  (setup :refine 4)
  ;; (sb-profile:profile)
  ;; (sb-profile:reset)
  ;; (sb-profile:profile "CL-MPM")
  ;; (sb-profile:profile "CL-MPM/PARTICLE")
  ;; (sb-profile:profile "CL-MPM/CONSTITUTIVE")
  ;; (sb-profile:profile "CL-MPM/MESH")
  ;; (sb-profile:profile "CL-MPM/BC")
  ;; (sb-profile:profile "CL-MPM/CONSTITUTIVE")
  ;; (sb-profile:profile "CL-MPM/PENALTY")

  (setf (cl-mpm::sim-enable-damage *sim*) t)
  (setf (cl-mpm::sim-nonlocal-damage *sim*) t)
  (setf *enable-box-friction* t)
  (setf (cl-mpm:sim-damping-factor *sim*)
        (*
         1d-1
         (sqrt (cl-mpm:sim-mass-scale *sim*))
         (cl-mpm/setup::estimate-critical-damping *sim*))
        )
  ;; (when (slot-exists-p *sim* 'cl-mpm/damage::delocal-counter-max)
  ;;   (setf (cl-mpm/damage::sim-damage-delocal-counter-max *sim*) substeps))
  ;; (with-accessors ((mesh cl-mpm:sim-mesh)
  ;;                  (mps cl-mpm:sim-mps))
  ;;     *sim*
  ;;   (let* ((mp (aref mps 0)))
  ;;     (pprint
  ;;      (cl-mpm/particle::mp-local-length mp))
  ;;     (cl-mpm/damage::iterate-over-damage-bounds
  ;;      mesh
  ;;      mp
  ;;      (cl-mpm/particle::mp-local-length mp)
  ;;      (lambda (mesh mp node)))
  ;;     (time
  ;;      (loop repeat 1
  ;;            while *run-sim*
  ;;            do (progn
  ;;                 (
  ;;                  cl-mpm:iterate-over-mps
  ;;                  mps
  ;;                  (lambda (mp)
  ;;                    (cl-mpm/damage::iterate-over-damage-bounds
  ;;                     mesh
  ;;                     mp
  ;;                     (cl-mpm/particle::mp-local-length mp)
  ;;                     (lambda (mesh mp node)))
  ;;                    )))))
  ;;     ))

  (setf *displacement-increment* 0d0)
  (let ((disp-inc (/ 1d-3 10)))
    (time
     (loop repeat 100
           while *run-sim*
           do (progn
                (incf *displacement-increment* disp-inc)
                (cl-mpm::update-sim *sim*)
                (swank.live:update-swank)
                ))))
  ;; (sb-profile:report)
  )
(defun test-neighbour ()
  (declare (optimize (speed 3)))
  (setup :refine 64)
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (mps cl-mpm::sim-mps))
      *sim*
    (time
     (loop repeat 100
           while *run-sim*
           do
              (progn
                (swank.live:update-swank)
                (cl-mpm::iterate-over-mps
                 mps
                 (lambda (mp)
                   (cl-mpm::iterate-over-neighbours-shape-gimp-simd
                    ;; cl-mpm::iterate-over-neighbours-shape-linear
                    mesh
                    mp
                    (lambda (mesh mp node svp grads fsvp fgrads))))))))))

(defun run-conv ()
  (setf *run-sim* t)
  (loop for r in (list 1 2 4 8)
        while *run-sim*
        do
           (progn
             (setup :refine r :mps 4)
             (run)
             (vgplot:title (format nil "~D" r))
             (vgplot:print-plot (merge-pathnames (format nil "refine_~5,'0d.png" r)) :terminal "png size 1920,1080")
             (sleep 1d0)
             )))

(declaim (notinline tref-test))
;; (defun tref-test (x)
;;   (declare (cl-mpm/particle::particle x))
;;   (with-accessors ((pos cl-mpm/particle::mp-position))
;;       x
;;     (magicl:tref pos 0 0)))

(declaim (notinline tref-test-t))
(defun tref-test (x)
  (declare (optimize (speed 3)))
  (declare (cl-mpm/particle::particle x))
  (;with-accessors ((pos cl-mpm/particle::mp-position))
   with-slots ((pos cl-mpm/particle::position))
      x
    (declare (magicl::matrix/double-float pos))
    (magicl:tref pos 0 0)
    ))

(defun tref-test-t (x)
  (declare (optimize (speed 3)))
  (declare (cl-mpm/particle::particle x))
  (;with-accessors ((pos cl-mpm/particle::mp-position))
   with-slots ((pos cl-mpm/particle::position))
   x
   (declare (magicl::matrix/double-float pos))
   (cl-mpm/utils::mtref pos 0 0)
   ))



;; (let ((x (make-instance 'cl-mpm/particle::particle))
;;       (iters 100000000))
;;   (time
;;    (dotimes (i iters)
;;      (tref-test x)))
;;   (time
;;    (dotimes (i iters)
;;      (tref-test-t x))))


;; (defun row-major-mp (mp)
;;   (let* ((class (class-of mp))
;;          (all-slots (sb-mop:class-slots class)))
;;     (loop for slot in all-slots
;;           do
;;              (when (sb-mop:slot-boundp-using-class class mp slot)
;;                (let* ((value (sb-mop:slot-value-using-class class mp slot)))
;;                  (when value
;;                    (when (typep value 'magicl::matrix/double-float)
;;                      (pprint (magicl::matrix-layout value)))))))))

;; (loop for mp across (cl-mpm:sim-mps *sim*)
;;       do (row-major-mp mp))

(defun test-tref ()
  (let ((average-throughput 0d0)
        (iters 100)
        (runs 3)
        )
    (dotimes (i runs)
      (setup :refine 16 :mps 2)
      (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale 0.1d0))
      (setf (cl-mpm:sim-damping-factor *sim*)
            (* 0.1d0
               (sqrt (cl-mpm:sim-mass-scale *sim*))
               (cl-mpm/setup::estimate-critical-damping *sim*)))
      (setf (cl-mpm::sim-enable-damage *sim*) nil)
      (let ((dt (cl-mpm/utils::time-form
                 iters
                 (cl-mpm:update-sim *sim*))))
        (incf average-throughput (* dt iters))
        (format t "Throughput: ~F~%" dt)
        (format t "Throughput per mp: ~E~%" (/ dt (length (cl-mpm:sim-mps *sim*))))))
    (let ((time-per-mp (/ average-throughput (* iters runs (length (cl-mpm:sim-mps *sim*))))))
      (format t "Average time per mp: ~E~%" time-per-mp)
      (format t "Average throughput per mp: ~E~%" (/ 1d0 time-per-mp)))))

(defun test-staggered ()
  (setup :refine 4 :mps 2)
  (let* ((disp-inc 1d-3)
         (data-load (list))
         (data-disp (list))
         (data-damage (list))
         (data-iter (list))
         (load 0d0)
         (disp 0d0)
         (dt-scale 0.5d0)
         (enable-plasticity nil)
         (damage-0 0d0)
         (damage-new 0d0)
         )
    (vgplot:figure)
    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (change-class mp 'cl-mpm/particle::particle-chalk-brittle))

    (setf (cl-mpm::sim-enable-damage *sim*) nil)
    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 0.05d0
             ;; (sqrt (cl-mpm:sim-mass-scale *sim*))
             (cl-mpm/setup::estimate-critical-damping *sim*)))

    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (when (= (cl-mpm/particle::mp-index mp) 0)
               (setf (cl-mpm/particle::mp-enable-plasticity mp) enable-plasticity)))
    (cl-mpm/dynamic-relaxation:converge-quasi-static
     *sim*
     :energy-crit 1d-2
     :oobf-crit 1d-2
     :dt-scale 0.50d0
     :substeps 10
     :conv-steps 100
     :post-iter-step
     (lambda (i energy oobf)
       ))
    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (when (= (cl-mpm/particle::mp-index mp) 0)
               (setf (cl-mpm/particle::mp-enable-plasticity mp) t)))
    (incf *displacement-increment* disp-inc)
    (loop for i from 0 to 50
          while (and *run-sim*
                     (or (= i 0) (> damage-new damage-0)))
          do
             (progn
               (setf damage-0 (get-damage))
               (cl-mpm/dynamic-relaxation:converge-quasi-static
                *sim*
                :energy-crit 1d-2
                :oobf-crit 1d-2
                :dt-scale 0.50d0
                :substeps 10
                :conv-steps 100
                :post-iter-step
                (lambda (i energy oobf)))
               (cl-mpm/damage::calculate-damage *sim*)
               (setf damage-new (get-damage))
               (cl-mpm/output:save-vtk (merge-pathnames "./output/" (format nil "sim_~5,'0d.vtk" i)) *sim*)
               (plot *sim*)
               (setf load (get-load))
               (setf disp *displacement-increment*)
               (push load data-load)
               (push disp data-disp)
               (push (get-damage) data-damage)
               (push i data-iter)
               ;; (vgplot:plot data-iter data-load)
               (pprint load)
               (pprint disp)))
    (vgplot:figure)
    (vgplot:title "Load")
    (vgplot:plot data-iter data-load)
    (vgplot:figure)
    (vgplot:title "Disp")
    (vgplot:plot data-iter data-disp)
    (vgplot:figure)
    (vgplot:title "Damage")
    (vgplot:plot data-iter data-damage)
    )
  )

(defun test ()
  (setf *run-sim* t)
  (let ((refine 4))
    (loop for s in (list
                    10d4
                    20d4
                    30d4
                    )
          while *run-sim*
          do
             (progn
               (setup :refine refine :mps 2 :surcharge-load s)
               (run (format nil "../ham-shear-box/output-~D-~F/" refine s))))))

(defun test-refine ()
  (setf *run-sim* t)
  (loop for s in (list
                  2
                  4
                  8
                  ;; 2
                  ;; 16
                  )
        while *run-sim*
        do
           (progn
             (setup :refine s :mps 2 :surcharge-load 10d4)
             (run (format nil "../ham-shear-box/output-~D-10e4/" s)))))

(defun test-mc (exx eyy ezz eyz ezx exy)
  (let* ((eps (cl-mpm/constitutive::swizzle-coombs->voigt
               (cl-mpm/utils:voigt-from-list (list exx eyy ezz eyz ezx exy))))
         (E 1d0)
         (nu 0.1d0)
         (phi 1d0)
         (psi 0d0)
         (c 1d0)
         (de (cl-mpm/constitutive::linear-elastic-matrix E nu))
         (sig (magicl:@ de eps)))
    (multiple-value-bind (sig eps f) (cl-mpm/constitutive::mc-plastic sig de eps E nu phi psi c)
      ;; (pprint (cl-mpm/constitutive::swizzle-voigt->coombs sig))
      (pprint (cl-mpm/constitutive::swizzle-voigt->coombs eps))
      ;; (pprint (cl-mpm/utils:vector-from-list (multiple-value-list (cl-mpm/damage::principal-stresses-3d sig))))
      ;; (format t "Recalculated f ~E~%" (cl-mpm/constitutive::mc-yield-func (cl-mpm/utils:vector-from-list (multiple-value-list (cl-mpm/damage::principal-stresses-3d sig))) angle c))
      )
    ))
