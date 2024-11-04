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
      (when (cl-mpm/particle::mp-enable-damage mp)
        ;; (setf damage-increment (cl-mpm/damage::criterion-max-principal-stress stress))

        ;; (setf damage-increment (cl-mpm/damage::tensile-energy-norm strain E de))
        (setf damage-increment
              (max 0d0
                   (
                    ;cl-mpm/damage::criterion-dp-inscribe-coheasion
                    cl-mpm/damage::criterion-dp-coheasion
                    ;; cl-mpm/damage::criterion-dp-tensile
                    ;; stress
                    (magicl:scale stress (/ 1d0 (magicl:det def)))
                    (* angle (/ pi 180d0)))))

        ;; (setf damage-increment
        ;;       (max 0d0
        ;;            (;cl-mpm/damage::criterion-dp-tensile
        ;;             cl-mpm/damage::criterion-dp-coheasion
        ;;             ;; cl-mpm/damage::drucker-prager-criterion
        ;;             ;; stress
        ;;             (cl-mpm/fastmaths:fast-scale-voigt stress (/ 1d0 (magicl:det def)))
        ;;             (* angle (/ pi 180d0)))))

        ;; (setf damage-increment
        ;;       (cl-mpm/damage::tensile-energy-norm (cl-mpm/fastmaths:fast-.+
        ;;                                            strain
        ;;                                            (cl-mpm/fastmaths:fast-scale-voigt
        ;;                                             plastic-strain
        ;;                                             1d0
        ;;                                             ;; (- 1d0 damage)
        ;;                                             ))
        ;;                                           E
        ;;                                           de))

        ;; (let ((es (cl-mpm/constitutive::linear-elastic-mat
        ;;            (cl-mpm/fastmaths:fast-.+
        ;;             strai
        ;;             n
        ;;             (magicl:scale plastic-strain
        ;;                           0d0
        ;;                           ;; (- 1d0 damage)
        ;;                           )
        ;;             )
        ;;                                                    de)))
        ;; (setf damage-increment
        ;;       (max 0d0
        ;;            (sqrt
        ;;             (cl-mpm/constitutive::voigt-j2
        ;;              (cl-mpm/utils:deviatoric-voigt stress)))))

        ;; (let ((es (cl-mpm/constitutive::linear-elastic-mat (cl-mpm/fastmaths:fast-.+ strain plastic-strain)
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
        ;;                         (cl-mpm/fastmaths:fast-.+
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
  ;; (sleep 1)
  (plot-load-disp)
  ;; (plot-conv)
  ;; (plot-domain)
  ;; (vgplot:plot '(0) '(0))
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
     (mapcar (lambda (x) (* x 500d0)) *data-plastic*)
     "Plastic"
     (mapcar (lambda (x) (* x 1d3)) *data-disp*)
     (mapcar (lambda (x) (* x 1d3)) *data-energy*)
     ;; (mapcar (lambda (x) (* x 1d-1)) *data-energy*)
     "Energy"
     ;; (mapcar (lambda (x) (* x 1d3)) *data-disp*)
     ;; (mapcar (lambda (x) (* x 1d0)) *data-penalty-energy*)
     ;; "Penalty energy"
     ))
  (vgplot:xlabel "Displacement (mm)")
  (vgplot:ylabel "Shear stress (kN/m^2)")
  )

(defparameter *enable-box-friction* t)
(defparameter *displacement-increment* 0d0)
(defparameter *box-size* 0d0)


;;Lazy and unhyginic
(defmacro defmpgen (name &rest args)
  `(defmacro ,name ()
     `(cl-mpm::add-mps
       sim
       (cl-mpm/setup::make-mps-from-list
        (cl-mpm/setup::make-block-mps-list
         offset
         block-size
         (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
         density
         ,@',args
         :index 0
         :gravity 0d0
         )))))
;; (macrolet ((defmpgen (name &rest args)
;;              `(defmacro ,name ()
;;                 `(cl-mpm::add-mps
;;                   sim
;;                   (cl-mpm/setup::make-mps-from-list
;;                    (cl-mpm/setup::make-block-mps-list
;;                     offset
;;                     block-size
;;                     (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
;;                     density
;;                     ,@',args
;;                     :index 0
;;                     :gravity 0d0
;;                     )))))))

(defmpgen make-mps-elastic
  'cl-mpm/particle::particle-elastic
  :E 1d9
  :nu 0.24d0
  )
(defmpgen make-mps-vm

  'cl-mpm/particle::particle-vm
  :E 1d9
  :nu 0.24d0
  :rho 100d3)
(defmpgen make-mps-mc-residual
  'cl-mpm/particle::particle-mc
  :E 1d9
  :nu 0.24d0
  :psi (* 5d0 (/ pi 180))
  :phi (* 30d0 (/ pi 180))
  :c 0d0
  :phi-r (* 30d0 (/ pi 180))
  :c-r 0d0
  :softening 0d0)

(defmpgen make-mps-mc-peak
  'cl-mpm/particle::particle-mc
  :E 1d9
  :nu 0.24d0
  :psi (* 5d0 (/ pi 180))
  :phi (* 42d0 (/ pi 180))
  :c 131d3
  :phi-r (* 30d0 (/ pi 180))
  :c-r 0d0
  :softening 0d0)

(defmpgen make-mps-dp-peak
  'cl-mpm/particle::particle-dp
  :E 1d9
  :nu 0.24d0
  :psi (* 5d0 (/ pi 180))
  :phi (* 42d0 (/ pi 180))
  :c 131d3
  :phi-r (* 30d0 (/ pi 180))
  :c-r 0d0
  :softening 0d0)

(defmpgen make-mps-mc-softening
  'cl-mpm/particle::particle-mc
  :E 1d9
  :nu 0.24d0
  :psi (* 5d0 (/ pi 180))
  :phi (* 42d0 (/ pi 180))
  :c 131d3
  :phi-r (* 30d0 (/ pi 180))
  :c-r 0d0
  :softening 50d0)

(defmpgen make-mps-damage
  'cl-mpm/particle::particle-chalk-delayed
  :E 1d9
  :nu 0.24d0
  :kt-res-ratio 1d-9
  :kc-res-ratio 1d0
  :g-res-ratio 1d-9
  :friction-angle 42d0
  :initiation-stress init-stress;18d3
  :delay-time 1d-2
  :delay-exponent 1d0
  :damage 0.0d0
  :ductility ductility
  :local-length length-scale
  :local-length-damaged 10d-10
  :enable-damage t
  :enable-plasticity nil

  :psi (* 5d0 (/ pi 180))
  :phi (* 42d0 (/ pi 180))
  :c 131d3
  :phi-r (* 30d0 (/ pi 180))
  :c-r 0d0
  :softening 10d0)

(defmpgen make-mps-plastic-damage
  'cl-mpm/particle::particle-chalk-delayed
  :E 1d9
  :nu 0.24d0
  :kt-res-ratio 1d-9
  :kc-res-ratio 1d-2
  :g-res-ratio 1d-3
  :friction-angle 40d0
  :initiation-stress init-stress;18d3
  :delay-time 1d-2
  :delay-exponent 1d0
  :damage 0.0d0
  :ductility ductility
  :local-length length-scale
  :local-length-damaged 10d-10
  :enable-damage t
  :enable-plasticity t

  :psi (* 5d0 (/ pi 180))
  :phi (* 42d0 (/ pi 180))
  :c (* 131d3 10d0)
  :phi-r (* 30d0 (/ pi 180))
  :c-r 0d0
  :softening 0d0
  )

(declaim (notinline setup-test-column))
(defun setup-test-column (size offset block-size &optional (e-scale 1) (mp-scale 1)
                          &key
                            (angle 0d0)
                            (friction 0.1d0)
                            (surcharge-load 72.5d3)
                            (piston-scale 1d0)
                            (piston-mps 2)
                            )
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               ;; :sim-type 'cl-mpm::mpm-sim-usf
               :sim-type 'cl-mpm/damage::mpm-sim-damage
               ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 1.7d3)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size)))
    (declare (double-float h density))
    (progn
      (let* ((angle-rad (* angle (/ pi 180)))
             ;; (init-stress (* 1 131d3))
             ;(init-stress 131d3)
             (init-stress 131d3)
             ;; (init-stress 60d3)
             ;; (init-stress 100d3)
             ;; (init-stress 30d3)
             ;; (init-stress 60d3)
             ;; (init-stress 10d3)
             ;; (gf 5d0)
             (gf 5d0)
             ;; (length-scale (* 1 h))
             (length-scale (* 7.5d-3 1))
             (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress 1d9))
             ;; (ductility 30d0)
             ;; (ductility 1d60)
             ;; (ductility 10d0)
             ;; (ductility 5d0)
             )
        (format t "Estimated ductility ~E~%" ductility)
        ;; (make-mps-plastic-damage)
        ;; (make-mps-mc-residual)
        ;; (make-mps-mc-peak)
        (make-mps-mc-peak)
        ;; (make-mps-dp-peak)
        ;; (make-mps-damage)
        ;; (make-mps-plastic-damage)
        )
      (let* ((sur-height h-x)
             (sur-height (* 0.5 (second block-size)))
             (sur-size (list 0.06d0 sur-height))
                                        ;(load 72.5d3)
             (load surcharge-load)
             ;; (gravity 10d0)
             ;; (density (/ load (* gravity sur-height)))
             (gravity (if (> load 0d0)
                          (/ load (* density sur-height))
                          0d0))
             (mp-surcharge t))
        ;; (loop for mp across (cl-mpm:sim-mps sim)
        ;;       do
        ;;          (with-accessors ((stress cl-mpm/particle:mp-stress)
        ;;                           (strain cl-mpm/particle:mp-strain)
        ;;                           (E cl-mpm/particle::mp-E)
        ;;                           (nu cl-mpm/particle::mp-nu)
        ;;                           (de cl-mpm/particle::mp-elastic-matrix))
        ;;              mp
        ;;            (let* (
        ;;                   (strains (cl-mpm/utils:voigt-from-list (list
        ;;                                                           0d0
        ;;                                                           (- (/ (* surcharge-load (+ 1d0 nu) (- 1d0 (* nu 2)))
        ;;                                                                 (* E (- 1d0 nu))))
        ;;                                                           0d0
        ;;                                                           0d0
        ;;                                                           0d0
        ;;                                                           0d0)))
        ;;                   )
        ;;              (setf stress (magicl:@ de strains)
        ;;                    strain strains))))
        (format t "Gravity ~F~%" gravity)

        ;; (if mp-surcharge
        ;;     (cl-mpm::add-mps
        ;;      sim
        ;;      (cl-mpm/setup::make-mps-from-list
        ;;       (cl-mpm/setup::make-block-mps-list
        ;;        (mapcar #'+ offset (list 0d0 (second block-size)))
        ;;        sur-size
        ;;        (mapcar (lambda (e) (* e e-scale piston-mps)) sur-size)
        ;;        density
        ;;        'cl-mpm/particle::particle-elastic-damage
        ;;        :E (* 1d9 1d0)
        ;;        :nu 0.24d0
        ;;        :initiation-stress 1d20
        ;;        :local-length 0d0
        ;;        :index 1
        ;;        :gravity (- gravity))))
        ;;     (progn
        ;;       (defparameter *pressure-bc*
        ;;         (cl-mpm/buoyancy::make-bc-pressure
        ;;          sim
        ;;          0d0
        ;;          (- surcharge-load)
        ;;          :clip-func
        ;;          (lambda (pos)
        ;;            (and
        ;;             (> (cl-mpm/utils:varef pos 1)
        ;;                (+ (second offset) (* 0.5d0 (second block-size))))
        ;;             ;; (> (cl-mpm/utils:varef pos 0)
        ;;             ;;    (first offset))
        ;;             ;; (< (cl-mpm/utils:varef pos 0)
        ;;             ;;    (+ (first offset) (first block-size)))
        ;;             )

        ;;            )))
        ;;       (cl-mpm:add-bcs-force-list
        ;;        sim
        ;;        *pressure-bc*)))
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

;; (defmethod cl-mpm::post-stress-step (mesh (mp cl-mpm/particle::particle-vm) dt)
;;   ;; (with-accessors ((ps cl-mpm/particle::mp-strain-plastic-vm)
;;   ;;                  (rho cl-mpm/particle::mp-rho))
;;   ;;     mp
;;   ;;   (let* ((rho-0 200d3)
;;   ;;          (rho-1 200d2)
;;   ;;          (soft 100d0)
;;   ;;          )
;;   ;;     (setf rho (+ rho-1 (* (- rho-0 rho-1) (exp (- (* soft ps))))))))
;;   )

;; (defmethod cl-mpm::post-stress-step (mesh (mp cl-mpm/particle::particle-mc) dt)
;;   ;; (with-accessors ((ps cl-mpm/particle::mp-strain-plastic-vm)
;;   ;;                  (c cl-mpm/particle::mp-c)
;;   ;;                  (phi cl-mpm/particle::mp-phi)
;;   ;;                  )
;;   ;;     mp
;;   ;;   (let ((phi_0 (* 42d0 (/ pi 180)))
;;   ;;         (phi_1 (* 30d0 (/ pi 180)))
;;   ;;         (c_0 131d3)
;;   ;;         (soft 1000d0;(* 100d0 *mesh-resolution*)
;;   ;;               ))
;;   ;;     (setf
;;   ;;      c (* c_0 (exp (- (* soft ps))))
;;   ;;      phi (+ phi_1 (* (- phi_0 phi_1) (exp (- (* soft ps))))))))
;;   )
;; (defmethod cl-mpm::post-stress-step (mesh (mp cl-mpm/particle::particle-chalk-delayed) dt))


(defun make-penalty-object (bc)
  )
(defun save-json-penalty-box (filename sim)
  (with-accessors ((bcs cl-mpm:sim-bcs-force-list))
      sim
    (let ((bc-save-list (list)))
      (loop for bc-list in bcs
            do
               (loop for bc across bc-list
                     do
                        (typecase bc
                          (cl-mpm/penalty::bc-penalty-distance
                           (push bc bc-save-list)
                           )
                          (cl-mpm/penalty::bc-penalty-structure
                           (loop for sub-bc in (cl-mpm/penalty::bc-penalty-structure-sub-bcs bc)
                                 do
                                    (when (typep sub-bc 'cl-mpm/penalty::bc-penalty-distance)
                                      (push sub-bc bc-save-list)))) (t ))))
      (let ((json-object (list)))
        (loop for bc in bc-save-list
              do
                 (with-accessors ((center cl-mpm/penalty::bc-penalty-distance-center-point )
                                  (normal cl-mpm/penalty::bc-penalty-normal)
                                  (radius cl-mpm/penalty::bc-penalty-distance-radius))
                     bc
                   (push (list :position (list (cl-mpm/utils:varef center 0)
                                               (cl-mpm/utils:varef center 1)
                                               (cl-mpm/utils:varef center 2))
                               :normal (list (cl-mpm/utils:varef normal 0)
                                             (cl-mpm/utils:varef normal 1)
                                             (cl-mpm/utils:varef normal 2))
                               :radius radius) json-object))
              )
        (str:to-file
         filename
         (jonathan:to-json
          json-object
          ))))))

(defun save-vtk-penalty-box (filename sim)
  (with-accessors ((bcs cl-mpm:sim-bcs-force-list))
      sim
    (let ((bc-save-list (list)))
      (loop for bc-list in bcs
            do
               (loop for bc across bc-list
                     do
                        (typecase bc
                          (cl-mpm/penalty::bc-penalty-distance
                           (push bc bc-save-list)
                           )
                          (cl-mpm/penalty::bc-penalty-structure
                           (loop for sub-bc in (cl-mpm/penalty::bc-penalty-structure-sub-bcs bc)
                                 do
                                    (when (typep sub-bc 'cl-mpm/penalty::bc-penalty-distance)
                                      (push sub-bc bc-save-list)))) (t ))))
      (let ((json-object (list)))
        (loop for bc in bc-save-list
              do
                 (with-accessors ((center cl-mpm/penalty::bc-penalty-distance-center-point )
                                  (normal cl-mpm/penalty::bc-penalty-normal)
                                  (radius cl-mpm/penalty::bc-penalty-distance-radius))
                     bc
                   (push (list :position (list (cl-mpm/utils:varef center 0)
                                               (cl-mpm/utils:varef center 1)
                                               (cl-mpm/utils:varef center 2))
                               :normal (list (cl-mpm/utils:varef normal 0)
                                             (cl-mpm/utils:varef normal 1)
                                             (cl-mpm/utils:varef normal 2))
                               :radius radius) json-object))
              )
        (str:to-file
         filename
         (jonathan:to-json
          json-object
          ))))))

(defun save-vtk-penalty-box (filename sim)
  (with-accessors ((bcs cl-mpm:sim-bcs-force-list))
      sim
    (let ((bc-save-list (list)))
      (loop for bc-list in bcs
            do
               (loop for bc across bc-list
                     do
                        (typecase bc
                          (cl-mpm/penalty::bc-penalty-distance
                           (push bc bc-save-list)
                           )
                          (cl-mpm/penalty::bc-penalty-structure
                           (loop for sub-bc in (cl-mpm/penalty::bc-penalty-structure-sub-bcs bc)
                                 do
                                    (when (typep sub-bc 'cl-mpm/penalty::bc-penalty-distance)
                                      (push sub-bc bc-save-list))))
                          (t ))))
      (let ((json-object (list)))
        (loop for bc in bc-save-list
              do
                 (with-accessors ((center cl-mpm/penalty::bc-penalty-distance-center-point )
                                  (normal cl-mpm/penalty::bc-penalty-normal)
                                  (radius cl-mpm/penalty::bc-penalty-distance-radius))
                     bc
                   (let* ((line-dir (cl-mpm/penalty::2d-orthog normal))
                          (p1 (cl-mpm/fastmaths:fast-.+ center (cl-mpm/fastmaths:fast-scale-vector line-dir radius)))
                          (p2 (cl-mpm/fastmaths:fast-.- center (cl-mpm/fastmaths:fast-scale-vector line-dir radius))))
                     (push (list (list (cl-mpm/utils:varef p1 0)
                                       (cl-mpm/utils:varef p1 1)
                                       (cl-mpm/utils:varef p1 2))
                                 (list (cl-mpm/utils:varef p2 0)
                                       (cl-mpm/utils:varef p2 1)
                                       (cl-mpm/utils:varef p2 2))) json-object))))
        (with-open-file (fs filename :direction :output :if-exists :supersede)
          (format fs "# vtk DataFile Version 2.0~%")
          (format fs "Lisp generated vtk file, SJVS~%")
          (format fs "ASCII~%")
          (format fs "DATASET UNSTRUCTURED_GRID~%")

          (let* ((line-count (length json-object))
                 (point-count (* 2 line-count)))
            (format fs "POINTS ~d double~%" point-count)
            (loop for bc in json-object
                  do (loop for p in bc
                           do (format fs "~E ~E ~E ~%"
                                      (coerce (first p) 'single-float)
                                      (coerce (second p) 'single-float)
                                      0e0)))

            (format fs "~%")
            (format fs "CELLS ~D ~D~%" line-count (* 3 line-count))
            (loop for bc in json-object
                  for i from 0
                  do (format fs "~D ~D ~D~%" 2 (* i 2) (+ 1 (* i 2))))

            (format fs "CELL_TYPES ~D~%" line-count)
            (loop for bc in json-object
                  do (format fs "3~%"))
            ))))))

(defun 2d-orthog (vec)
  (cl-mpm/utils:vector-from-list (list (cl-mpm/utils:varef vec 1) (- (cl-mpm/utils:varef vec 0) 0))))

(defun make-beveled-right-angle (sim epsilon p-a corner p-b bevel-size)
  (let* ((diff-a (cl-mpm/fastmaths:fast-.- p-a corner))
         (diff-b (cl-mpm/fastmaths:fast-.- corner p-b))
         (normal-a (2d-orthog diff-a))
         (normal-b (2d-orthog diff-b))
         (length-a (- (cl-mpm/fastmaths:mag diff-a) bevel-size))
         (length-b (- (cl-mpm/fastmaths:mag diff-a) bevel-size))
         (corner-a (cl-mpm/fastmaths:fast-.+ p-a (cl-mpm/fastmaths:fast-scale! (cl-mpm/fastmaths:norm diff-a)
                                                                               length-a)))
         (corner-b (cl-mpm/fastmaths:fast-.+ p-b (cl-mpm/fastmaths:fast-scale! (cl-mpm/fastmaths:norm diff-b)
                                                                               length-b)))
         (la
           (cl-mpm/penalty::make-bc-penalty-line-segment sim p-a corner-a epsilon 0d0 0d0))
         (lb
           (cl-mpm/penalty::make-bc-penalty-line-segment sim corner-b p-b epsilon 0d0 0d0))
         (bevel
           (cl-mpm/penalty::make-bc-penalty-line-segment sim corner-a corner-b epsilon 0d0 0d0)))
    (values la bevel lb)
    ))


(declaim (notinline make-penalty-box))
(defun make-penalty-box (sim left-x right-x height friction-scale offset &key (epsilon-scale 1d2)
                                                                     (corner-size 1d0)
                                                                     )
  (let* ((left-normal (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0)))
         (right-normal (cl-mpm/utils:vector-from-list (list -1d0 0d0 0d0)))
         (plane-normal (cl-mpm/utils:vector-from-list (list 0d0 -1d0 0d0)))
         (plane-normal-left (cl-mpm/utils:vector-from-list (list 0d0 1d0 0d0)))
         (epsilon (* epsilon-scale 1d9))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (extra-height 0d0)
         (friction (* friction-scale (tan (* 30d0 (/ pi 180)))))
         ;; (friction 0.9d0)
         ;; (friction 1d0)
         ;; (friction 0.0d0)
         (gap-height (* 0 h))
         ;; (corner (* 0.5d0 h))
         (corner
           corner-size
           ;; (* 1d0 7.5d-3)
                 )
         (corner-x (* 1d0 corner))
         (corner-y (* 1d0 corner))
         (damping 0d0)
         (all-bcs (list))
         )
    (defparameter *box-left* left-x)
    (defparameter *box-size* (* 2d0 height))
    ;; (cl-mpm/penalty::make-bc-penalty-distance-point
    ;;  sim
    ;;  left-normal
    ;;  (cl-mpm/utils:vector-from-list (list left-x (+ (* 3d0 height)
    ;;                                                 offset
    ;;                                                 gap-height
    ;;                                                 ) 0d0))
    ;;  (* 2d0 height)
    ;;  epsilon
    ;;  0d0
    ;;  damping)
    (defparameter *shear-box-left-static*
      (cl-mpm/penalty::make-bc-penalty-line-segment
       sim
       (cl-mpm/utils:vector-from-list (list left-x (* 4 height) 0d0))
       (cl-mpm/utils:vector-from-list (list left-x (+ corner-y height) 0d0))
       epsilon
       0d0
       0d0))
    (defparameter *shear-box-left-static-bevel*
      (cl-mpm/penalty::make-bc-penalty-line-segment
       sim
       (cl-mpm/utils:vector-from-list (list left-x (+ corner-y height) 0d0))
       (cl-mpm/utils:vector-from-list (list (- left-x corner-x) height 0d0))
       epsilon
       0d0
       0d0))
    (defparameter *shear-box-right-static*
      (cl-mpm/penalty::make-bc-penalty-line-segment
       sim
       (cl-mpm/utils:vector-from-list (list right-x (+ corner-y height) 0d0))
       (cl-mpm/utils:vector-from-list (list right-x (* 4 height) 0d0))
       epsilon
       0d0
       0d0))
    (defparameter *shear-box-right-static-bevel*
      (cl-mpm/penalty::make-bc-penalty-line-segment
       sim
       (cl-mpm/utils:vector-from-list (list (+ right-x corner-x) height 0d0))
       (cl-mpm/utils:vector-from-list (list right-x (+ corner-y height) 0d0))
       epsilon
       0d0
       0d0))
    (defparameter *shear-box-right-static-slide*
      (cl-mpm/penalty::make-bc-penalty-line-segment
       sim
       (cl-mpm/utils:vector-from-list (list (+ right-x corner-x height) height 0d0))
       (cl-mpm/utils:vector-from-list (list (+ right-x corner-x) height 0d0))
       epsilon
       0d0
       0d0))

    (defparameter *shear-box-left-slide*
      (cl-mpm/penalty::make-bc-penalty-line-segment
       sim
       (cl-mpm/utils:vector-from-list (list (- left-x (+ corner-y height)) height 0d0))
       (cl-mpm/utils:vector-from-list (list (- left-x corner-x) height 0d0))
       epsilon
       0d0
       0d0))

    (defparameter *shear-box-left-bevel*
      (cl-mpm/penalty::make-bc-penalty-line-segment
       sim
       (cl-mpm/utils:vector-from-list (list (- left-x corner-x) height 0d0))
       (cl-mpm/utils:vector-from-list (list left-x (- height corner-y) 0d0))
       epsilon
       0d0
       0d0))
    (defparameter *shear-box-left-dynamic*
      (cl-mpm/penalty::make-bc-penalty-line-segment
       sim
       (cl-mpm/utils:vector-from-list (list left-x (- height corner-x) 0d0))
       (cl-mpm/utils:vector-from-list (list left-x 0d0 0d0))
       epsilon
       0d0
       0d0))
    (defparameter *shear-box-right-dynamic*
      (cl-mpm/penalty::make-bc-penalty-line-segment
       sim
       (cl-mpm/utils:vector-from-list (list right-x 0d0 0d0))
       (cl-mpm/utils:vector-from-list (list right-x (- height corner-y) 0d0))
       epsilon
       0d0
       0d0)
      )
    (defparameter *shear-box-right-dynamic-bevel*
      (cl-mpm/penalty::make-bc-penalty-line-segment
       sim
       (cl-mpm/utils:vector-from-list (list right-x (- height corner-y) 0d0))
       (cl-mpm/utils:vector-from-list (list (+ right-x corner-x) height 0d0))
       epsilon
       0d0
       0d0)
      )
    (defparameter *shear-box-floor*
      (cl-mpm/penalty::make-bc-penalty-distance-point
       sim
       plane-normal-left
       (cl-mpm/utils:vector-from-list (list (+ left-x height) offset 0d0))
       (* 2 height)
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
        ;; *shear-box-left-static*
        *shear-box-left-dynamic*
        *shear-box-left-slide*
        *shear-box-left-bevel*
        ;; *shear-box-floor*
        )))

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
        ;; *shear-box-right-static*
        ;; *shear-box-right-static-bevel*
        ;; *shear-box-right-static-slide*

        ;; *shear-box-right-dynamic*
        ;; *shear-box-right-dynamic-bevel*

        *shear-box-left-static*
        *shear-box-left-static-bevel*

        *shear-box-left-dynamic*
        *shear-box-left-slide*
        *shear-box-left-bevel*
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
        *shear-box-right-static-bevel*
        *shear-box-right-static-slide*

        *shear-box-right-dynamic*
        *shear-box-right-dynamic-bevel*

        ;; *shear-box-left-static*
        ;; *shear-box-left-static-bevel*

        ;; *shear-box-left-dynamic*
        ;; *shear-box-left-slide*
        ;; *shear-box-left-bevel*
        ;; *shear-box-floor*
        )))

    (dolist
        (bc
         (list
          *shear-box-right-static*
          *shear-box-right-static-bevel*
          *shear-box-right-static-slide*
          *shear-box-left-static*
          *shear-box-left-static-bevel*
          ))
      (cl-mpm/penalty::bc-increment-center bc (cl-mpm/utils:vector-from-list (list 0d0 gap-height 0d0))))
    (dolist
        (bc
         (list
          *shear-box-right-static*
          *shear-box-right-static-bevel*
          *shear-box-right-static-slide*
          *shear-box-right-dynamic*
          *shear-box-right-dynamic-bevel*
          *shear-box-left-static*
          *shear-box-left-static-bevel*
          *shear-box-left-dynamic*
          *shear-box-left-slide*
          *shear-box-left-bevel*
          ))
      (cl-mpm/penalty::bc-increment-center bc (cl-mpm/utils:vector-from-list (list 0d0 offset 0d0))))
    (defparameter *last-pos* 0d0)
    (defparameter *shear-box-controller*
      (cl-mpm/bc::make-bc-closure
       nil
       (lambda ()

         (setf (cl-mpm/penalty::bc-penalty-load *shear-box-left-dynamic*) 0d0)
         (setf (cl-mpm/penalty::bc-penalty-load *shear-box-right-dynamic*) 0d0)

         (let ((delta (- *displacement-increment* *last-pos*)))
           (dolist
               (bc
                (list
                 *shear-box-right-dynamic*
                 *shear-box-right-dynamic-bevel*
                 *shear-box-left-dynamic*
                 *shear-box-left-bevel*
                 *shear-box-left-slide*
                 ))
             (cl-mpm/penalty::bc-increment-center bc (cl-mpm/utils:vector-from-list (list delta 0d0 0d0))))
           (setf *last-pos* *displacement-increment*))
         (let ((friction-bcs (list
                              *shear-box-struct-left*
                              *shear-box-struct-right*
                              *shear-box-right-static*
                              *shear-box-right-static-bevel*
                              *shear-box-right-static-slide*
                              *shear-box-right-dynamic*
                              *shear-box-right-dynamic-bevel*
                              *shear-box-left-static*
                              *shear-box-left-static-bevel*
                              *shear-box-left-dynamic*
                              *shear-box-left-slide*
                              *shear-box-left-bevel*
                              )))
           ;; (break)
           ;; (pprint (if *enable-box-friction*
           ;;             friction
           ;;             0d0))
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
        ;; *shear-box-struct-right-static*
        ;; *shear-box-struct-left*
        ;; *shear-box-struct-right*
        ;; *shear-box-struct-right-dynamic*
        *shear-box-struct-floor*
        *shear-box-struct*
        *shear-box-struct-right*
        ))))))

(declaim (notinline get-load))
(defun get-load ()
  ;; (cl-mpm/penalty::bc-penalty-load *true-load-bc*)
  ;; (cl-mpm/penalty::bc-penalty-load *shear-box-left-dynamic*)
  (-
   (cl-mpm/penalty::bc-penalty-load *shear-box-left-dynamic*)
   (cl-mpm/penalty::bc-penalty-load *shear-box-right-dynamic*)))

(defun reset-load ()
  (setf 
   (cl-mpm/penalty::bc-penalty-load *true-load-bc*)
   0d0))

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
         (*
          (cl-mpm/particle::mp-mass mp)
          (cl-mpm/particle::mp-strain-plastic-vm mp))
         0d0))
   #'+
   (cl-mpm:sim-mps *sim*)))


(defgeneric get-piston-load (sim))
(defmethod get-piston-load ((sim cl-mpm:mpm-sim))
  (cl-mpm/penalty::bc-penalty-load *piston-penalty*))
(defmethod get-piston-load ((sim cl-mpm/mpi:mpm-sim-mpi))
  (cl-mpm/mpi:mpi-sum (cl-mpm/penalty::bc-penalty-load *piston-penalty*)))

(defun reset-piston-load ()
  (setf (cl-mpm/penalty::bc-penalty-load *piston-penalty*) 0d0))

(declaim (notinline make-piston))
(defun make-piston (box-size box-offset surcharge-load epsilon-scale piston-scale)
  (defparameter *piston-penalty*
    (cl-mpm/penalty::make-bc-penalty-distance-point 
     *sim*
     (cl-mpm/utils:vector-from-list (list 0d0 -1d0 0d0))
     (cl-mpm/utils:vector-from-list (list (* 1.5 box-size)
                                          (+ box-size box-offset) 0d0))
     (* 1.25d0 0.5d0 0.06d0)
     (* 1d9 epsilon-scale)
     0d0
     0d0))
  (defparameter *piston-confinement* 0d0)
  (let* ((control-stiffness
           ;; 1d9
           ;; 1d-12
           ;; 1d0
           ;; 1d-6
           (/ (* piston-scale 1d6) (cl-mpm/penalty::bc-penalty-epsilon *piston-penalty*))
           )
         (target-load surcharge-load)
         (window-size 1)
         (window-index 0)
         (window (make-array window-size :initial-element 0d0)))
    (defparameter *window* window)
    (defparameter *piston-controller*
      (cl-mpm/bc::make-bc-closure
       nil
       (lambda ()
         (let* ((inst-load (/ (get-piston-load *sim*) 0.06d0))
                (current-load inst-load))
           (setf (aref window (mod window-index window-size))
                 current-load)
           (incf window-index)
           (setf current-load (/ (loop for x across window sum x)
                                 window-size))
           ;; (format t "~A~%" window)
           ;; (format t "~F~%" current-load)
           (let* ((delta (- (- current-load target-load)))
                  (penalty-inc (* (cl-mpm:sim-dt *sim*) (* delta control-stiffness))))
             (declare (double-float penalty-inc))
             ;; (incf (the double-float (cl-mpm/penalty::bc-penalty-datum *piston-penalty*)) penalty-inc)
             (cl-mpm/penalty::bc-increment-center *piston-penalty* (cl-mpm/utils:vector-from-list  (list 0d0 (- penalty-inc) 0d0)))
             (incf *piston-confinement* current-load)
             (reset-piston-load)
             )
           )))))
  (when (> surcharge-load 0d0)
    (cl-mpm:add-bcs-force-list
     *sim*
     *piston-penalty*)
    (cl-mpm:add-bcs-force-list
     *sim*
     *piston-controller*)))
(declaim (notinline setup))
(defun setup (&key (refine 1d0) (mps 2) (friction 0.0d0) (surcharge-load 72.5d3) (epsilon-scale 1d2)
                (piston-scale 1d0)
                (piston-mps 2)
                )
  (defparameter *displacement-increment* 0d0)
  (let* ((mps-per-dim mps)
         (mesh-size (/ 0.03d0 refine))
         ;; (mesh-size (/ 0.03d0 refine)))
         (sunk-size 0.03d0)
         (box-size (* 2d0 sunk-size))
         (domain-size (* 3d0 box-size))
         (box-offset (* mesh-size 2d0))
         (offset (list box-size box-offset)))
    (setf *box-size* box-size)
    (defparameter *sim* (setup-test-column
                         (list domain-size (+ (* 2 box-size) box-offset))
                         offset
                         (list box-size box-size)
                         (/ 1d0 mesh-size)
                         mps-per-dim
                         :piston-scale piston-scale
                         :piston-mps piston-mps
                         :surcharge-load surcharge-load))
    (make-penalty-box *sim* box-size (* 2d0 box-size) sunk-size friction box-offset
                      :epsilon-scale epsilon-scale
                      :corner-size (* 0.125d0 mesh-size))
    (make-piston box-size box-offset surcharge-load epsilon-scale piston-scale)
    ;; (cl-mpm/setup:remove-sdf
    ;;  *sim*
    ;;  (lambda (p)
    ;;    (if (> (- (cl-mpm/utils:varef p 1) box-offset) sunk-size)
    ;;        0d0
    ;;        1d0)))

    ;; (cl-mpm/setup:remove-sdf
    ;;  *sim*
    ;;  (lambda (p)
    ;;    (if (<  (cl-mpm/utils:varef p 1) mesh-size)
    ;;        0d0
    ;;        1d0)))
    (defparameter *true-load-bc* *shear-box-left-dynamic*)
    ;; (cl-mpm/setup::damage-sdf *sim* (lambda (p) (cl-mpm/setup::line-sdf
    ;;                                              p
    ;;                                              (list box-size sunk-size 0d0)
    ;;                                              (list (* 2 box-size) sunk-size 0d0)
    ;;                                              mesh-size
    ;;                                              )) 1.0d0)

    (let ((mp (aref (cl-mpm:sim-mps *sim*) 0)))
      (when (typep mp 'cl-mpm/particle::particle-damage)
        (let* ((length-scale (cl-mpm/particle::mp-local-length mp))
               (width (* 1 length-scale)))
          (cl-mpm/setup::apply-sdf *sim*
                                   (lambda (p) (cl-mpm/setup::line-sdf
                                                p
                                                (list box-size (+ sunk-size box-offset) 0d0)
                                                (list (* 2 box-size) (+ sunk-size box-offset) 0d0)
                                                width
                                                ))
                                   (lambda (mp v)
                                     (when (< v 0d0)
                                       (setf (cl-mpm/particle:mp-damage mp)
                                             1d0
                                             ;; (exp (- (expt (/ (+ width v) width) 8)))
                                             ))
                                     )))))
    )
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (format t "Mesh-size: ~E~%" (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defun stop ()
  (setf *run-sim* nil
        cl-mpm/dynamic-relaxation::*run-convergance* nil))

(declaim (notinline run))
(defun run (&optional (output-directory "./output/")
            &key (total-time 6d-2)
              (damping 1d-2)
              (time-scale 1d0)
              (damage-time-scale 1d0)
              (sample-scale 1d0)
              (dt-scale 0.5d0)
              (displacment 0.5d-3)
              (skip-level nil)
              (enable-damage t)
              (enable-plasticity t)
              )

  ;; (when (not (uiop:file-exists-p output-directory)))
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
  (defparameter *data-energy* (list))
  (defparameter *data-penalty-energy* (list))
  (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :supersede)
    (format stream "disp,load,plastic,damage,energy~%"))
  ;; (vgplot:close-all-plots)
  (let* (;(displacment 1d-3)
         ;(total-time (* 50d0 displacment))
         (time-per-mm (* 100d0 time-scale))
         (total-time (* time-per-mm displacment))
         (load-steps (round (* sample-scale 500 (/ displacment 1d-3))))
         (target-time (/ total-time load-steps))
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale dt-scale)
         (load-0 0d0)
         ;; (enable-plasticity
         ;;   nil
         ;;   ;; t
         ;;   ;; (cl-mpm/particle::mp-enable-plasticity (aref (cl-mpm:sim-mps *sim*) 0))
         ;;   )
         ;; (enable-damage
         ;;   ;; nil
         ;;   t
         ;;   )
         (max-load 0d0)
         (average-surcharge 0d0)
         (total-steps 0)
         (disp-inc (/ displacment load-steps)))
    ;;Disp rate in test 4d-4mm/s -> 4d-7mm/s
    (format t "Loading rate: ~E~%" (/ displacment (* load-steps target-time)))
    (format t "Total time: ~E~%" total-time)
    (format t "Target time: ~E~%" target-time)
    (format t "Strain rate: ~E~%" (/ displacment total-time))

    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (when (typep mp 'cl-mpm/particle::particle-damage)
              (when (= (cl-mpm/particle::mp-index mp) 0)
                (setf (cl-mpm/particle::mp-delay-time mp) (* damage-time-scale (* time-per-mm 1d-3) 1d-4)))))

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

    ;; (vgplot:figure)
    (cl-mpm/output:save-vtk (merge-pathnames output-directory (format nil "sim_conv_~5,'0d.vtk" 0)) *sim*)
    (cl-mpm/dynamic-relaxation:converge-quasi-static
     *sim*
     :energy-crit 1d-2
     :oobf-crit 1d-2
     :dt-scale dt-scale
     :substeps 10
     :conv-steps 5000
     :post-iter-step
     (lambda (i energy oobf)
       ;; (plot-domain)
       (cl-mpm/output:save-vtk (merge-pathnames output-directory (format nil "sim_conv_~5,'0d.vtk" (+ 1 i))) *sim*)
       (cl-mpm/output:save-vtk-nodes (merge-pathnames output-directory (format nil "sim_nodes_conv_~5,'0d.vtk" (+ 1 i))) *sim*)
       (format t "Surcharge load ~E~%" (/ *piston-confinement* 10))
       (setf *piston-confinement* 0d0)
       ))
    ;; (vgplot:figure)

    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (when (= (cl-mpm/particle::mp-index mp) 0)
               (setf (cl-mpm/particle::mp-enable-plasticity mp) enable-plasticity)))
    ;; (cl-mpm:iterate-over-mps
    ;;  (cl-mpm:sim-mps *sim*)
    ;;  (lambda (mp)
    ;;    (when (= (cl-mpm/particle::mp-index mp) 0)
    ;;      (when (> (abs (- (cl-mpm/utils:varef (cl-mpm/particle::mp-position mp) 1) 0.03)) (* 0.5d0 0.03d0))
    ;;        (setf
    ;;         (cl-mpm/particle::mp-enable-damage mp) nil
    ;;         (cl-mpm/particle::mp-index mp) 3
    ;;         )))
    ;;    ))

    (setf *enable-box-friction* t)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (*
           damping
           ;; (sqrt (cl-mpm:sim-mass-scale *sim*))
           (cl-mpm/setup::estimate-critical-damping *sim*)))

    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))

    (setf substeps (round target-time (cl-mpm:sim-dt *sim*)))

    (when (slot-exists-p *sim* 'cl-mpm/damage::delocal-counter-max)
      (setf (cl-mpm/damage::sim-damage-delocal-counter-max *sim*) substeps))
    (defparameter *displacement-increment* 0d0)
    (format t "Substeps ~D~%" substeps)
    ;; (vgplot:figure)
    (setf cl-mpm/penalty::*debug-force* 0)
    (setf (cl-mpm::sim-enable-damage *sim*) enable-damage)

    (reset-load)
    (let ((disp-av 0d0)
          (load-av 0d0)
          (d-av 0d0)
          (p-av 0d0)
          (e-av 0d0))
      (format t "Running unloaded step ~D~%" substeps)
      (dotimes (i substeps)
        (cl-mpm::update-sim *sim*)
        (incf load-av (/ (get-load) substeps)))
      ;; (setf load-av (get-load))
      (setf load-0 load-av)
      (format t "Load-0 ~E~%" load-0)
      (push *t* *data-t*)
      (push disp-av *data-disp*)
      (push load-av *data-v*)
      (push d-av *data-damage*)
      (push p-av *data-plastic*)
      (push e-av *data-energy*)
      (setf load-av (get-load))
      (setf disp-av *displacement-increment*)
      (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :append)
        (format stream "~f,~f,~f,~f,~f~%" disp-av load-av p-av d-av e-av)))
    (reset-load)

    (time (loop for steps from 0 below load-steps
                while *run-sim*
                do
                   (progn
                     (format t "Step ~d/~D~%" steps load-steps)
                     (when (= (mod steps (ceiling sample-scale)) 0)
                       (cl-mpm/output:save-vtk (merge-pathnames output-directory (format nil "sim_~5,'0d.vtk" *sim-step*)) *sim*)
                       (save-json-penalty-box (merge-pathnames output-directory (format nil "sim_pb_~5,'0d.json" *sim-step*)) *sim*)
                       (save-vtk-penalty-box (merge-pathnames output-directory (format nil "sim_pb_~5,'0d.vtk" *sim-step*)) *sim*)
                       ;; (save-vtk-penalty-box (merge-pathnames output-directory (format nil "sim_box_~5,'0d.vtk" *sim-step*)) *sim*)
                       )
                     ;; (cl-mpm/output::save-vtk-nodes (merge-pathnames output-directory (format nil "sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     (let ((load-av 0d0)
                           (disp-av 0d0)
                           (high-resolution nil)
                           (p-av 0d0)
                           (d-av 0d0)
                           (e-av 0d0)
                           (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
                       (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :append)
                         (reset-load)
                         (time
                          (dotimes (i substeps)
                            (cl-mpm::update-sim *sim*)
                            ;; (cl-mpm/penalty::apply-penalty-nudge *sim* *shear-box-struct*)
                            (incf load-av (/ (get-load) substeps))
                            (incf disp-av (/ *displacement-increment* substeps))
                            (incf *displacement-increment* (/ disp-inc substeps))
                            (incf e-av (/ (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*) substeps))
                            (incf *t* (cl-mpm::sim-dt *sim*))

                            (when high-resolution
                              (push *t* *data-t*)
                              (push *displacement-increment* *data-disp*)
                              (push (get-load) *data-v*)
                              (push (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*) *data-energy*)
                              (push (get-damage) *data-damage*)
                              (push (get-plastic) *data-plastic*)
                              (format stream "~f,~f,~f,~f,~f~%" *displacement-increment* (get-load) (get-plastic) (get-damage) (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*)))
                            )))
                       (sb-ext:gc)

                       (setf load-av (get-load))
                       (setf disp-av *displacement-increment*)

                       (setf max-load (max max-load load-av))
                       (when skip-level
                         (when (and (< load-av (* max-load skip-level))
                                    (< 10 steps))
                           (skip)))

                       (format t "Surcharge load ~E~%" (/ *piston-confinement* substeps))
                       (incf average-surcharge (/ *piston-confinement* substeps))
                       (setf *piston-confinement* 0d0)
                       (unless high-resolution
                         (push *t* *data-t*)
                         (push disp-av *data-disp*)
                         (push load-av *data-v*)
                         (setf d-av (get-damage))
                         (setf p-av (get-plastic))
                         (push e-av *data-energy*)
                         (push d-av *data-damage*)
                         (push p-av *data-plastic*)
                         (format t "Disp ~E - Load ~E~%" disp-av load-av)
                         (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :append)
                           (format stream "~f,~f,~f,~f,~f~%" disp-av load-av p-av d-av e-av)))
                       )
                     (incf *sim-step*)
                     (plot *sim*)
                     (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                        :terminal "png size 1920,1080")
                     (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       (format t "CFL dt estimate: ~f~%" dt-e)
                       (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf substeps substeps-e))
                     (swank.live:update-swank))))

    (format t "Average surcharge ~F~%" (/ average-surcharge *sim-step*))
    ;; (loop for mp across (cl-mpm:sim-mps *sim*)
    ;;       do (progn
    ;;            (cl-mpm/fastmaths::fast-zero (cl-mpm/particle::mp-acceleration mp))
    ;;            (cl-mpm/fastmaths::fast-zero (cl-mpm/particle:mp-velocity mp))))
    ;; (when nil
    ;;   (defparameter *true-load-bc* *shear-box-right-dynamic*)
    ;;   (time (loop for steps from 0 below load-steps
    ;;               while *run-sim*
    ;;               do
    ;;                  (progn
    ;;                    (format t "Step ~d ~%" steps)
    ;;                    (cl-mpm/output:save-vtk (merge-pathnames output-directory (format nil "sim_~5,'0d.vtk" *sim-step*)) *sim*)
    ;;                    (cl-mpm/output::save-vtk-nodes (merge-pathnames output-directory (format nil "sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
    ;;                    (let ((load-av 0d0)
    ;;                          (disp-av 0d0)
    ;;                          (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
    ;;                      (time
    ;;                       (dotimes (i substeps)
    ;;                         (cl-mpm::update-sim *sim*)
    ;;                         (incf load-av (/ (get-load) substeps))
    ;;                         (incf disp-av (/ *displacement-increment* substeps))
    ;;                         (incf *displacement-increment* (/ (- disp-inc) substeps))
    ;;                         (incf *t* (cl-mpm::sim-dt *sim*))))

    ;;                      ;; (setf load-av cl-mpm/penalty::*debug-force*)
    ;;                      ;; (setf load-av (get-load))
    ;;                      ;; (setf disp-av *displacement-increment*)

    ;;                      (push *t* *data-t*)
    ;;                      (push disp-av *data-disp*)
    ;;                      (push load-av *data-v*)
    ;;                      (format t "Disp ~E - Load ~E~%" disp-av load-av)
    ;;                      (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :append)
    ;;                        (format stream "~f,~f~%" disp-av load-av)))
    ;;                    (incf *sim-step*)
    ;;                    (plot *sim*)
    ;;                    (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
    ;;                                       :terminal "png size 1920,1080")
    ;;                    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
    ;;                      (format t "CFL dt estimate: ~f~%" dt-e)
    ;;                      (format t "CFL step count estimate: ~D~%" substeps-e)
    ;;                      (setf substeps substeps-e))
    ;;                    (swank.live:update-swank)))))
    )
  ;; (vgplot:figure)
  (plot-load-disp))


(declaim (notinline run-static))
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
  (let* ((load-steps 5)
         (dt (cl-mpm:sim-dt *sim*))
         (dt-scale 0.5d0)
         (displacment 0.2d-3)
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
    (setf (cl-mpm::sim-enable-damage *sim*) nil)
    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (when (= (cl-mpm/particle::mp-index mp) 0)
               (setf (cl-mpm/particle::mp-enable-plasticity mp) enable-plasticity)))
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
                                     :energy-crit 1d-3
                                     :oobf-crit 1d-3
                                     :substeps 10
                                     :conv-steps 200
                                     :dt-scale dt-scale
                                     :post-iter-step
                                     (lambda (i energy oobf)))
                                    ;; (cl-mpm/damage::calculate-damage *sim*)
                                    ;; (plot *sim*)
                                    (setf damage-1 (get-damage))
                                    (format t "Staggered damage diff: ~E~%" (- damage-1 damage-0))
                                    ))
                         (when (> damage-1 damage-0)
                           (break "Staggered error")))

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
              :mps 4
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
  (setup :refine 2)
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
         1d-2
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
  (let ((s 20d4))
    (loop for r in (list
                    2
                    4
                    8
                    16
                    32
                    )
          while *run-sim*
          do
             (progn
               (setup :refine r :mps 2 :surcharge-load s)
               (run
                (format nil "../ham-shear-box/output-~D-~F/" r s)
                )
               (vgplot:title (format nil "~D" r))
               (vgplot:print-plot (merge-pathnames (format nil "refine_~5,'0d.png" r)) :terminal "png size 1920,1080")
               (sleep 1d0)
               ))))







(defun test-epsilon ()
  (setf *run-sim* t)
  (loop for eps in (list
                    ;; 1d2
                    1d3
                    ;; 1d4
                    )
        for dt-scale in (list
                         ;; 0.5d0
                         ;; 0.30d0
                         0.25d0
                         )
        while *run-sim*
        do
           (loop for r in (list 4)
                 while *run-sim*
                 do
                    (loop for s in (list
                                    10d4
                                    20d4
                                    ;; 30d4
                                    )
                          while *run-sim*
                          do
                             (progn
                               (setup :refine r :mps 3 :surcharge-load s
                                      :epsilon-scale eps
                                      ;:dt-scale (/ 0.5d0 (expt eps 1/3))
                                      )
                               (run (format nil "../ham-shear-box/output-~F_~F-~F/" eps r s)
                                    :time-scale 0.5d0
                                    :sample-scale 1d0
                                    :dt-scale dt-scale;0.5d0
                                    ))))))







(defparameter *ductility* 10d0)




(defun test-damage ()
  (setf *run-sim* t)
  (loop for d in (list
                  1.0d0
                  0.0d0
                  0.5d0
                  0.75d0
                  )
        do
           (loop for refine in (list 4)
                 do
                    (loop for s
                          ;; from 0d0 to 100d4 by 10d4
                            in
                            (list
                             10d4
                             20d4
                             ;; 30d4

                             ;; 15d4
                             ;; 25d4
                             ;; 35d4

                             )
                          while *run-sim*
                          do
                             (let ((scale 0.5d0))
                               (setup :refine refine :mps 4 :surcharge-load s)
                               (cl-mpm:iterate-over-mps
                                (cl-mpm:sim-mps *sim*)
                                (lambda (mp)
                                  (when (typep mp 'cl-mpm/particle::particle-chalk-delayed)
                                    (setf (cl-mpm/particle::mp-damage mp) d)
                                    ;; (setf (cl-mpm/particle::mp-shear-residual-ratio mp) d)
                                    ;; (setf (cl-mpm/particle::mp-k-compressive-residual-ratio mp) d)
                                    )
                                  ))
                               (run (format nil "../ham-shear-box/output-~F_~F-~F/" refine d s)
                                    :time-scale scale
                                    :sample-scale 1d0
                                    :damage-time-scale 1d0
                                    )
                               ;; (run-static (format nil "../ham-shear-box/output-~D-~F/" refine s))
                               )))))



(defparameter *skip* nil)

(defun test ()
  (setf *run-sim* t)
  (loop for refine in (list
                       ;; 2
                       ;; 4
                       8
                       ;; 16
                       ;; 32
                       ;; 4.5
                       ;; 8.5
                       ;; 16
                       ;; 32
                       )
        do
           (dolist (mps (list 3))
             (let (;(mps 2)
                   ;; (mps 2)
                   ;; (scale 0.5d0)
                   )
               (loop for s
                     ;; from 0d0 to 35d4 by 2.5d4
                       in
                       (list
                        ;; 5d4
                        ;; 10d4
                        20d4
                        30d4
                        )
                     while (and *run-sim*)
                     do
                        (let ((scale 1d0))
                          (setf *skip* nil)
                          (format t "Test ~D ~F" refine s)
                          (setup :refine refine :mps mps :surcharge-load s
                                 :epsilon-scale 1d2
                                 :piston-scale 0.1d0
                                 :piston-mps 2
                                 :friction 0d0)
                          (run (format nil "../ham-shear-box/output-~f_~D_~f-~F/" refine mps scale s)
                               :displacment 1d-3
                               :time-scale (* 1d0 scale)
                               :sample-scale (* 1d0 0.1d0)
                               :dt-scale 0.5d0
                               :damage-time-scale 1d0
                               ;; :skip-level 0.9d0
                               :enable-damage nil
                               :enable-plasticity t
                               )
                          (when *skip*
                            (setf *run-sim* t))))))))


(defun test-damage ()
  (setf *run-sim* t)
  (loop for refine in (list
                       ;; 2
                       ;; 4
                       ;; 8
                       ;; 16
                       ;; 32
                       )
        do
           (loop for d in (list 0d0
                                1d0
                                0.5d0)
                 do
                    (let (;(mps 2)
                          (mps 4)
                          (scale 0.5d0))
                      (loop for s
                            ;; from 0d0 to 100d4 by 10d4
                            ;; from 0d0 to 40d4 by 5d4
                              in
                              (list
                               10d4
                               20d4
                               30d4
                               )
                            while *run-sim*
                            do
                               (progn
                                 (format t "Test ~D ~F" refine s)
                                 (setup :refine refine :mps mps :surcharge-load s
                                        :epsilon-scale 1d2
                                        :friction 0d0)

                                 (cl-mpm:iterate-over-mps
                                  (cl-mpm:sim-mps *sim*)
                                  (lambda (mp)
                                    (when (typep mp 'cl-mpm/particle::particle-chalk-delayed)
                                      (setf (cl-mpm/particle::mp-damage mp) d))))

                                 (run (format nil "../ham-shear-box/output-~f_~D_~f_~f-~F/" refine mps scale d s)
                                      :displacment 0.5d-3
                                      :time-scale (* 1d0 scale)
                                      :sample-scale (* 1d0 1d0)
                                      :dt-scale 0.3d0
                                      :damage-time-scale 1d0)
                                 ;; (run-static (format nil "../ham-shear-box/output-~D-~F/" refine s))
                                 ))))))

(defun test-possion () (setf *run-sim* t)
  (loop for refine in (list
                       ;; 2
                       4
                       ;; 8
                       ;; 16
                       ;; 32
                       )
        do
           (loop for d in (list
                           0.1d0
                           0.2d0
                           0.30d0
                           0.25d0
                           )
                 do
                    (let (;(mps 2)
                          (mps 2)
                          (scale 2d0))
                      (loop for s
                            ;; from 0d0 to 100d4 by 10d4
                            ;; from 0d0 to 40d4 by 5d4
                              in
                              (list
                               10d4
                               20d4
                               30d4
                               )
                            while *run-sim*
                            do
                               (progn
                                 (setf *skip* nil)
                                 (format t "Test ~D ~F" refine s)
                                 (setup :refine refine :mps mps :surcharge-load s
                                        :epsilon-scale 1d3
                                        :friction 0d0)

                                 (cl-mpm:iterate-over-mps
                                  (cl-mpm:sim-mps *sim*)
                                  (lambda (mp)
                                    (when (typep mp 'cl-mpm/particle::particle-chalk-delayed)
                                      (setf (cl-mpm/particle::mp-nu mp) d))))

                                 (run (format nil "../ham-shear-box/output-~f_~D_~f_~f-~F/" refine mps scale d s)
                                      :displacment 0.10d-3
                                      :time-scale (* 1d0 scale)
                                      :sample-scale (* 1d0 5d0)
                                      :dt-scale 0.100d0
                                      :damage-time-scale 1d0)
                                 (when *skip*
                                   (setf *run-sim* t))
                                 ))))))

(defun skip ()
  (setf *skip* t
        *run-sim* nil))

(defun test-angle ()
  (setf *run-sim* t)
  (loop for refine in (list
                       ;; 2
                       4
                       ;; 8
                       ;; 16
                       ;; 4.5
                       ;; 8.5
                       ;; 16
                       ;; 32
                       )
        do
           (let (;(mps 2)
                 (mps 4)
                 ;; (scale 0.5d0)
                 )

             (loop for d in (list
                             ;; 10d0 20d0
                             ;; 10d0
                             40d0
                             ;; 40d0 50d0 60d0 70d0
                             )
                   do
                      (loop for s
                            ;; from 0d0 to 100d4 by 10d4
                            ;; from 0d0 to 40d4 by 5d4
                              in
                              (list
                               ;; 1d4
                               ;; 2d4
                               ;; 3d4
                               10d4
                               20d4
                               30d4

                               ;; 90d4
                               ;; 80d4
                               ;; 60d4
                               ;; 50d4
                               )
                            while (and *run-sim*)
                            do
                               (let ((scale 0.5d0))
                                 (setf *skip* nil)
                                 (format t "Test ~D ~F" refine s)
                                 (setup :refine refine :mps mps :surcharge-load s
                                        :epsilon-scale 1d4
                                        :piston-scale 1d0
                                        :piston-mps 2
                                        :friction 0d0)
                                 (cl-mpm:iterate-over-mps
                                  (cl-mpm:sim-mps *sim*)
                                  (lambda (mp)
                                    (when (typep mp 'cl-mpm/particle::particle-chalk-delayed)
                                      (setf (cl-mpm/particle::mp-friction-angle mp) d))))

                                 (run (format nil "../ham-shear-box/output-~f_~D_~f_MC-~F/" refine mps d s)
                                      :displacment 0.10d-3
                                      :time-scale (* 1d0 scale)
                                      :sample-scale (* 1d0 1d0)
                                      :dt-scale 0.100d0
                                      :damage-time-scale 1d0)
                                 (when *skip*
                                   (setf *run-sim* t))
                                 ))))))
