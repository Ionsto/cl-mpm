(defpackage :cl-mpm/examples/slump
  (:use :cl))
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)

(defun cl-mpm/damage::length-localisation (local-length local-length-damaged damage)
  ;; (+ (* local-length (- 1d0 damage)) (* local-length-damaged damage))
  (* local-length (max (sqrt (- 1d0 damage)) 1d-10))
  ;; local-length
  )
;; (defmethod print-object ((object magicl:matrix) stream)
;;   (pprint object stream))

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-chalk-delayed) dt fbar)
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
          ;; (setf damage-increment
          ;;       (cl-mpm/damage::tensile-energy-norm-pressure strain E nu de 0d0))
          ;;            (angle )
          ;(setf damage-increment (criterion-mc strain (* angle (/ pi 180d0)) E nu))
          ;; (setf damage-increment (criterion-dp strain (* angle (/ pi 180d0)) E nu))
          ;; (setf damage-increment (cl-mpm/damage::smooth-rankine-criterion (magicl:scale stress (/ 1d0 (magicl:det def)))))
          (setf damage-increment
                (max 0d0
                     (cl-mpm/damage::drucker-prager-criterion
                      (magicl:scale stress (/ 1d0 (magicl:det def))) (* angle (/ pi 180d0)))))
          )
        (when (>= damage 1d0)
          (setf damage-increment 0d0))
        ;;Delocalisation switch
        (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
        (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
        ))))
(in-package :cl-mpm/examples/slump)
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
                                                         (temp cl-mpm/particle::mp-boundary)
                                                         )
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
    ;; (cl-mpm/fastmath::voigt-tensor-reduce-simd (cl-mpm/particle::mp-velocity-rate mp))
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
;; (defun plot (sim &optional (plot :damage)))
(defun plot (sim &optional (plot :damage))
  (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:format-plot t "set object 1 rect from 0,0 to ~f,~f fc rgb 'blue' fs transparent solid 0.5 noborder behind" ms-x *water-height*))
  (vgplot:format-plot t "set style fill solid")
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :deformed
   :colour-func #'cl-mpm/particle::mp-damage
   )


  )
;; (defun plot (sim &optional (plot :damage))
;;   (declare (optimize (speed 0) (debug 3)))
;;   (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
;;   (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
;;          (ms-x (first ms))
;;          (ms-y (second ms))
;;          )
;;     (vgplot:format-plot t "set object 1 rect from 0,0 to ~f,~f fc rgb 'blue' fs transparent solid 0.5 noborder behind" ms-x *water-height*))

;;   (multiple-value-bind (x y c stress-y lx ly e density temp vx ybar dist)
;;     (loop for mp across (cl-mpm:sim-mps sim)
;;           collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
;;           collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
;;           collect (magicl:tref (cl-mpm::mp-velocity mp) 0 0) into vx
;;           collect (length-from-def sim mp 0) into lx
;;           collect (length-from-def sim mp 1) into ly
;;           collect (if (slot-exists-p mp 'cl-mpm/particle::damage) (cl-mpm/particle:mp-damage mp) 0) into c
;;           collect (if (slot-exists-p mp 'cl-mpm/particle::damage-ybar) (cl-mpm/particle::mp-damage-ybar mp) 0) into ybar
;;           collect (if (slot-exists-p mp 'cl-mpm/particle::temperature) (cl-mpm/particle:mp-temperature mp) 0) into temp
;;           collect (if (slot-exists-p mp 'cl-mpm/particle::strain-energy-density) (cl-mpm/particle::mp-strain-energy-density mp) 0) into e
;;           collect (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)) into density
;;           ;; collect (cl-mpm/particle:mp-volume mp) into density
;;           collect (max-stress mp) into stress-y
;;           collect 0d0 into dist
;;           finally (return (values x y c stress-y lx ly e density temp vx ybar dist)))

;;     (let* ((node-x '())
;;            (node-y '())
;;            (node-c '())
;;            (mesh (cl-mpm:sim-mesh *sim*))
;;            (nodes (cl-mpm/mesh:mesh-nodes mesh))
;;            )
;;       (dotimes (i (array-total-size nodes))
;;         (let ((n (row-major-aref nodes i)))
;;           (with-accessors ((index cl-mpm/mesh:node-index)
;;                            (boundary cl-mpm/mesh::node-boundary-node)
;;                            (boundary-s cl-mpm/mesh::node-boundary-scalar)
;;                            (active cl-mpm/mesh::node-active)
;;                            )

;;               n
;;             (let ((n-ratio (/ (cl-mpm/mesh::node-volume n) (cl-mpm/mesh::node-volume-true n))))
;;               (when active
;;                 (if boundary
;;                     ;; (print index)
;;                     (destructuring-bind (x y z) (cl-mpm/mesh:index-to-position mesh index)
;;                       ;; (push (nth 0 index) node-x)
;;                       ;; (push (nth 1 index) node-y)
;;                       (push x node-x)
;;                       (push y node-y)
;;                       (push 2 node-c)
;;                       ;; (push n-ratio node-c)
;;                       ;; (push boundary-s node-c)
;;                       )
;;                     (destructuring-bind (x y z) (cl-mpm/mesh:index-to-position mesh index)
;;                       (push x node-x)
;;                       (push y node-y)
;;                       ;; (push n-ratio node-c)
;;                       (push 0 node-c)
;;                       ;; (push boundary-s node-c)
;;                       ))
;;                 ;; (push (/ (cl-mpm/mesh::node-volume n) (cl-mpm/mesh::node-volume-true n)) node-c)
;;                 )))))
;;       (cond
;;         ((eq plot :point)
;;          ;; (vgplot:format-plot t "set cbrange [0:1]")
;;          ;; (vgplot:plot x y ";;with points pt 7")
;;                                         ;(vgplot:format-plot t "set cbrange [0:2]")
;;          (if node-x
;;              ;; (multiple-value-bind (xp yp) (max-spacial (map 'list #'identity (cl-mpm::sim-mps *sim*)))
;;              (progn
;;                (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min node-c) (+ 0.01 (apply #'max node-c)))
;;                (vgplot:plot
;;                 ;; x y ";;with points pt 7"
;;                 x y lx ly ";;with ellipses"
;;                 node-x node-y node-c ";;with points pt 7 lc palette"
;;                 ;; (list xp) (list yp) ";;with points"
;;                 ))
;;                ;)
;;              (vgplot:plot
;;               x y lx ly ";;with ellipses")
;;              )
;;          )
;;         ((eq plot :contact)
;;          (let* ((contact-points (cl-mpm/penalty::collect-contact-points-bc (cl-mpm:sim-mesh *sim*) (cl-mpm:sim-mps *sim*) *floor-bc*))
;;                 (c-x (mapcar (lambda (p) (magicl:tref p 0 0)) contact-points))
;;                 (c-y (mapcar (lambda (p) (magicl:tref p 1 0)) contact-points))
;;                 (c-c (mapcar (lambda (p) 1d0) contact-points))
;;            )
;;            ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 1e-20 (apply #'max stress-y)))
;;            (vgplot:format-plot t "set cbrange [~f:~f]" 0d0 (+ 1e-40 (apply #'max c)))
;;          ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min node-c) (+ 0.01 (apply #'max node-c)))
;;          (if c-x
;;              ;; (multiple-value-bind (xp yp) (max-spacial (map 'list #'identity (cl-mpm::sim-mps *sim*)))
;;              (vgplot:plot
;;               ;; x y ";;with points pt 7"
;;               x y c ";;with points pt 7 lc palette"
;;               c-x c-y ";;with points pt 7"
;;               ;; (list xp) (list yp) ";;with points"
;;               )
;;              (vgplot:plot
;;               ;; x y lx ly ";;with ellipses"
;;               x y c ";;with points pt 7 lc palette"
;;               )
;;              ))
;;          )
;;         ((eq plot :damage)
;;          (vgplot:format-plot t "set style fill solid")
;;          ;; (vgplot:format-plot t "set cbrange [0:1]")
;;          ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 1e-6 (apply #'max c)))
;;          (vgplot:format-plot t "set cbrange [~f:~f]" 0d0 (+ 1e-6 (apply #'max c)))
;;          ;; (vgplot:format-plot t "set cbrange [~f:~f]" (+ 0d0 (apply #'min c)) (+ 1e-6 (apply #'max c)))
;;          ;; (vgplot:plot x y c ";;with points pt 7 lc palette")
;;          (vgplot:plot x y lx ly c ";;with ellipses lc palette")
;;          )
;;         ((eq plot :dist)
;;          (vgplot:format-plot t "set cbrange [~f:~f]" 0d0 (+ 1e-6 (apply #'max dist)))
;;          (vgplot:plot x y dist ";;with points pt 7 lc palette")
;;          )
;;         ((eq plot :velocity)
;;          (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min vx) (+ 0.01 (apply #'max vx)))
;;          (vgplot:plot x y vx ";;with points pt 7 lc palette")
;;          )
;;         ((eq plot :temperature)
;;          (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min temp) (+ 0.01 (apply #'max temp)))
;;          ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 0.01 (apply #'max c)))
;;          (vgplot:plot x y temp ";;with points pt 7 lc palette")
;;          )
;;         ((eq plot :energy)
;;          (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min e) (+ 0.01 (apply #'max e)))
;;          (vgplot:plot x y e ";;with points pt 7 lc palette")
;;          )
;;         ((eq plot :stress)
;;          (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 1e-40 (apply #'max stress-y)))
;;          (vgplot:plot x y stress-y ";;with points pt 7 lc palette"))
;;         ((eq plot :damage-ybar)
;;          (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min ybar) (+ 1e-20 (apply #'max ybar)))
;;          (vgplot:plot x y ybar ";;with points pt 7 lc palette"))
;;         ((eq plot :deformed)
;;          ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 0.01 (apply #'max stress-y)))
;;          (vgplot:plot x y lx ly ";;with ellipses"))
;;         ((eq plot :density)
;;          (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min density) (+ 0.01 (apply #'max density)))
;;          (vgplot:plot x y e ";;with points pt 7 lc palette")
;;          )))
;;     )
;;   (vgplot:format-plot t "replot (~f*x + ~f)~%" *sliding-slope* *sliding-offset*)
;;   (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
;;          (ms-x (first ms))
;;          (ms-y (second ms))
;;          )
;;     (vgplot:axis (list 0 ms-x
;;                        0 ms-y))
;;     (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
;;     (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
;;       (vgplot:format-plot t "set ytics ~f" h)
;;       (vgplot:format-plot t "set xtics ~f" h))

;;   (vgplot:replot))

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
         (h (min 1d0 (max 0d0 (/ (cl-mpm/fastmath::dot pa ba)
                                 (cl-mpm/fastmath::dot ba ba)
                                 ))))
         (v (magicl:.- pa (magicl:scale ba h)))
         )
    (- (sqrt (cl-mpm/fastmath::dot v v)) width)))
(defun plane-sdf (position normal distance)
  (- distance (cl-mpm/fastmath::dot position (cl-mpm/fastmath::norm normal))))
(defun plane-point-sdf (position normal point)
  (let ((distance (cl-mpm/fastmath::dot point (cl-mpm/fastmath::norm normal))))
    (- distance (cl-mpm/fastmath::dot position (cl-mpm/fastmath::norm normal)))))

;; (defun test-iter ()
;;     (let ((iters 1000)
;;           (mps (cl-mpm:sim-mps *sim*))
;;           (mesh (cl-mpm:sim-mesh *sim*))
;;           )
;;       ;; (cl-mpm::update-sim *sim*)
;;       (format t "Direct iter")
;;       (time
;;        (loop repeat iters
;;              do
;;                 (loop for mp across mps
;;                       do
;;                          (cl-mpm::iterate-over-neighbours-shape-gimp
;;                           mesh
;;                           mp
;;                           (lambda (mesh mp node w g)
;;                             )))))
;;       (format t "Cached iter")
;;       (time
;;        (loop repeat iters
;;              do
;;                 (loop for mp across mps
;;                       do
;;                          (cl-mpm::iterate-over-neighbours
;;                           mesh
;;                           mp
;;                           (lambda (mesh mp node w g)
;;                             )))))
;;       ))

(declaim (notinline setup-test-column))
(defun setup-test-column (size block-size block-offset slope &optional (e-scale 1d0) (mp-scale 1d0)
                          &rest mp-args)
  (declare (optimize (speed 0)))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1 e-scale)
               (mapcar (lambda (s) (* s e-scale)) size)
               :sim-type 'cl-mpm/damage::mpm-sim-damage
               ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (angle 0d0)
         (density *ice-density*)
         )
    (progn
      (let ()
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-mps-from-list
               (cl-mpm/setup::make-block-mps-list
                ;; block-position
                block-offset
                block-size
                (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                density

                ;; 'cl-mpm/particle::particle-viscoplastic-damage
                ;; :E 1d9
                ;; :nu 0.3250d0
                ;; :visc-factor 111d6
                ;; :visc-power 3d0
                ;; :initiation-stress 0.3d6
                ;; :damage-rate 1d-5
                ;; :critical-damage 0.50d0
                ;; :local-length 100d0
                ;; :local-length-damaged 0.1d0
                ;; :damage 0.0d0
                
                'cl-mpm/particle::particle-visco-elasto-plastic-damage
                :visc-factor 111d6
                ;; :visc-factor 11d6
                :visc-power 3d0

                ;; 'cl-mpm/particle::particle-chalk-delayed
                :E 1d9
                :nu 0.33d0
                :enable-plasticity t

                :friction-angle 60.0d0

                :kt-res-ratio 1d-9
                :kc-res-ratio 1d-1
                :g-res-ratio 1d-5

                :fracture-energy 3000d0
                :initiation-stress 0.5d6;0.6d6
                :delay-time 1d4
                :ductility 10d0
                :critical-damage 1d0;(- 1.0d0 1d-3)
                :damage-domain-rate 0.9d0;This slider changes how GIMP update turns to uGIMP under damage
                :local-length 10d0;(* 0.20d0 (sqrt 7))
                :local-length-damaged 10d-10

                :psi (* 00d0 (/ pi 180))
                :phi (* 40d0 (/ pi 180))
                ;; :phi 0d0;(* 30d0 (/ pi 180))
                ;; :phi (* 70d0 (/ pi 180))
                :c 1000d3
                ;; :c 1000d3
                ;; :c 500d3

                :gravity -9.8d0
                ;; :gravity-axis (magicl:from-list '(0.5d0 0.5d0) '(2 1))
                :index 0
                ;:angle angle
                ;; :slope slope
                )))
        )
      (let ((mass-scale 1d6))
        (setf (cl-mpm::sim-mass-scale sim) mass-scale)
        (setf (cl-mpm:sim-damping-factor sim)
              ;; (* 1d-3 mass-scale)
              ;(* density 0.1d0)
              (* density 1d0)
              ;; 0d0
              ;; 0d0
              ;; 1d0
              ;; (* 0.00000001d0 mass-scale)
              ;; 0.1d0
              ;; 0.01d0
              ;; 0.1d0
              ;; 0.0d0
              ;; (*)1d0
              ;100d0
              )
        )
      (format t "Est dt: ~f~%"
              (* 1d0 h
                 (sqrt (cl-mpm::sim-mass-scale sim))
                 (sqrt (/ density (cl-mpm/particle::mp-p-modulus (aref (cl-mpm:sim-mps sim) 0))))))
      (setf (cl-mpm:sim-mass-filter sim) 1d-15)
      (setf (cl-mpm::sim-allow-mp-split sim) t)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) t)
      (setf (cl-mpm::sim-nonlocal-damage sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm:sim-dt sim) 1d-8)
      (setf (cl-mpm:sim-bcs sim) (make-array 0))
      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-var
             (cl-mpm:sim-mesh sim)
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0 nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 0 nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil nil 0)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil nil 0)))
             ;; (lambda (i) (cl-mpm/bc:make-bc-fixed (mapcar #'+ i '(0 0)) '(nil nil)))
             ;; (lambda (i) (cl-mpm/bc:make-bc-surface (mapcar #'+ i '(0 1))
             ;;                                        (magicl:from-list (list 0d0 1d0) '(2 1))))
             ;; (lambda (i) (cl-mpm/bc:make-bc-surface i (magicl:from-list (list 0d0 1d0) '(2 1))))
             ))
      (format t "Bottom level ~F~%" h-y)
      (let* ((terminus-size (+ (second block-size) (* slope (first block-size))))
             (ocean-x 1000)
            (ocean-y (+ h-y (* 0.90d0 0.5d0 terminus-size)))
             ;; (ocean-y (* (round (+ (* 0.0d0 0.5d0 terminus-size) h-y) h-y) h-y)
             ;;          )
            ;(angle -1d0)
            )

        ;; (loop for mp across (cl-mpm:sim-mps sim)
        ;;       do
        ;;          (with-accessors ((pos cl-mpm/particle:mp-position)
        ;;                           (stress cl-mpm/particle::mp-stress-kirchoff)
        ;;                           (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
        ;;                           (stress-cauchy cl-mpm/particle::mp-stress)
        ;;                           )
        ;;              mp
        ;;            (setf stress
        ;;                  (cl-mpm/utils:matrix-to-voight
        ;;                   (magicl:eye 2 :value (* 1d0 (cl-mpm/buoyancy::pressure-at-depth (magicl:tref pos 1 0) ocean-y *water-density*))))
        ;;                  stress-cauchy (magicl:scale stress 1d0)
        ;;                  undamaged-stress (magicl:scale stress 1d0)
        ;;                  )))

        (format t "Ocean level ~a~%" ocean-y)
        (defparameter *water-height* ocean-y)
        (defparameter *floor-bc*
          (cl-mpm/penalty::make-bc-penalty-point-normal
           sim
           (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
           (cl-mpm/utils:vector-from-list (list 00d0 (+ 0d0 h-y) 0d0))
           (* *ice-density* 1d5)
           0.9d0
           ;; 1d1
           ))
        (setf (cl-mpm::sim-bcs-force-list sim)
              (list
               (cl-mpm/bc:make-bcs-from-list
                (list
                 (cl-mpm/buoyancy::make-bc-buoyancy-clip
                  sim
                  ocean-y
                  *water-density*
                  (lambda (pos datum)
                    (>= (magicl:tref pos 1 0) h-y)
                    )
                  )
                 ;; (cl-mpm/bc::make-bc-closure
                 ;;  '(0 0)
                 ;;  (lambda ()
                 ;;    (setf (magicl:tref (cl-mpm/mesh::node-force
                 ;;                        (cl-mpm/mesh::get-node (cl-mpm:sim-mesh sim) '(24 1))
                 ;;                        ) 1 0)
                 ;;          1d3
                 ;;          )))
                 ))
               (cl-mpm/bc:make-bcs-from-list
                (list *floor-bc*)
                )
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
(defun setup ()
  (declare (optimize (speed 0)))
  (defparameter *run-sim* nil)
  (let* ((mesh-size 10)
         (mps-per-cell 2)
         (slope -0.02)
         (shelf-height 150)
         (shelf-aspect 1)
         (shelf-length (* shelf-height shelf-aspect))
         (shelf-end-height (+ shelf-height (* (- slope) shelf-length)))
         (shelf-height-terminus shelf-height)
         (shelf-height shelf-end-height)
         (offset (list 0
                       (* 1 mesh-size)
                       ))
         )
    (defparameter *sim*
      (setup-test-column (list (+ shelf-length (* 2 shelf-height))
                                                 (+ shelf-height 100)
                                                 )
                         (list shelf-length shelf-height
                               )
                         offset
                         slope
                         (/ 1 mesh-size) mps-per-cell))

    ;;Delete all the plotted frames
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
    ;; (remove-sdf *sim* (rectangle-sdf
    ;;                    (list shelf-length
    ;;                          shelf-height-terminus)
    ;;                    (list
    ;;                     200
    ;;                     (* 0.2 shelf-height-terminus))))
    ;; (damage-sdf *sim* (lambda (p) (line-sdf p
    ;;                                         (list (- shelf-length shelf-height) shelf-height)
    ;;                                         (list shelf-length 0d0)
    ;;                                         50d0
    ;;                                         )) 0.9d0)
    ;; (damage-sdf *sim* (lambda (p) (line-sdf p
    ;;                                         (list 0d0 shelf-height)
    ;;                                         (list shelf-length shelf-height)
    ;;                                         50d0
    ;;                                         )) 0.4d0)
    ;; (loop for mp across (cl-mpm:sim-mps *sim*)
    ;;       do
    ;;          (setf (cl-mpm/particle:mp-damage mp) (random 0.5d0)))
    ;; (let* ((normal (cl-mpm/fastmath::norm
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
    ;; (let* ((undercut-angle 10d0)
    ;;        (normal (magicl:from-list (list
    ;;                                        (cos (- (* pi (/ undercut-angle 180d0))))
    ;;                                        (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1))))
    ;;   (remove-sdf *sim* (lambda (p) (plane-point-sdf p
    ;;                                                  normal
    ;;                                                  (magicl:from-list (list shelf-length shelf-height)
    ;;                                                                    '(2 1) :type 'double-float)
    ;;                                                  ))))
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
  (sleep 1)
  (with-open-file (stream (merge-pathnames "output/terminus_position.csv") :direction :output :if-exists :supersede)
    (format stream "Time (s),Terminus position~%")
    (loop for tim in (reverse *time*)
          for x in (reverse *x-pos*)
          do (format stream "~f, ~f ~%" tim x)))
 (let* ((target-time 1d3)
         (dt (cl-mpm:sim-dt *sim*))
         (dt-scale 1d0)
         (substeps (floor target-time dt)))

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
   (let ((energy-crit 0d0))
     (time (loop for steps from 0 to 100
                 while *run-sim*
                 do
                    (progn
                      (let ((base-damping 1d-3))
                        (if (> steps 5)
                            (progn
                              (format t "Enabling damage~%")
                              (setf (cl-mpm::sim-enable-damage *sim*) t
                                    ;; (cl-mpm::sim-mass-scale *sim*) 1d3
                                    (cl-mpm::sim-damping-factor *sim*) (* 1d3 0.1)

                                    )
                              ;; (setf (cl-mpm::sim-damping-factor *sim*) base-damping
                                        ;dt-scale 1.0d0
                              );)
                            (progn
                              (setf (cl-mpm::sim-enable-damage *sim*) nil
                                    
                                    ;; (cl-mpm::sim-damping-factor *sim*)
                                    ;; (+ (* 1d-3 (cl-mpm::sim-mass-scale *sim*)
                                    ;;       ;; (exp (- steps))
                                    ;;       )
                                    ;;    base-damping)
                                    )
                              )
                            ))
                      (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                        (format t "CFL dt estimate: ~f~%" dt-e)
                        (format t "CFL step count estimate: ~D~%" substeps-e)
                        (setf substeps substeps-e))
                      ;; (when (> steps 2)
                      ;;   (setf dt-scale 1d0)
                      ;;     )
                      (format t "Step ~d ~%" steps)
                      (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                      ;; (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                      ;; (cl-mpm/output:save-csv (merge-pathnames (format nil "output/sim_~5,'0d.csv" *sim-step*)) *sim*)

                      (push *t* *time*)
                      ;; (let ((cfl (find-max-cfl *sim*)))
                      ;;   (push cfl *cfl-max*))
                      (setf *x*
                            (loop for mp across (cl-mpm:sim-mps *sim*)
                                  maximize (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0)))
                      (push
                       *x*
                       *x-pos*)
                      (setf energy-crit 0d0)
                      (let ((cfl 0))
                        (time (dotimes (i substeps)
                                ;; (increase-load *sim* *terminus-mps*
                                ;;                (magicl:from-list (list (* (cl-mpm:sim-dt *sim*) 1d4) 0d0) '(2 1)))
                                ;; (pescribe-velocity *sim* *terminus-mps* '(1d0 nil))

                                ;; (melt-sdf *sim* (rectangle-sdf (list 0 0) (list 5000 *water-level*))
                                ;;           (* (cl-mpm:sim-dt *sim*) 1d-3))
                                (cl-mpm::update-sim *sim*)
                                (incf energy-crit (* dt (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*)))
                                ;; (setf cfl (max cfl (find-max-cfl *sim*)))
                                (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))
                        ;; (loop for mp across (cl-mpm:sim-mps *sim*)
                        ;;       do
                        ;;          (with-accessors ((tv cl-mpm/particle::mp-time-averaged-visc))
                        ;;              mp
                        ;;            (setf tv (/ tv substeps))))
                        ;; (setf cfl (find-max-cfl *sim*))
                        (format t "CFL: ~f~%" cfl)
                        (push cfl *cfl-max*)
                        (format t "Energy norm: ~F~%" energy-crit)
                        (let* ((dt-e (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
                               (substeps-e (floor target-time dt-e)))
                          (format t "CFL dt estimate: ~f~%" dt-e)
                          (format t "CFL step count estimate: ~D~%" substeps-e)
                          (setf (cl-mpm:sim-dt *sim*) dt-e)
                          (setf substeps substeps-e))
                        )
                      (with-open-file (stream (merge-pathnames "output/terminus_position.csv") :direction :output :if-exists :append)
                        (format stream "~f, ~f ~%" *t* *x*))
                      (incf *sim-step*)
                      (plot *sim*)
                      (vgplot:title (format nil "Time:~F - Energy ~E"  *t* energy-crit))
                      (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                         :terminal "png size 1920,1080"
                                         )
                      (swank.live:update-swank)
                      (sleep .1)
                      )))))
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
(defun plot-disp-day ()
  (let* ((df (lisp-stat:read-csv
	            (uiop:read-file-string #P"stokes.csv")))
         (time-day (mapcar (lambda (x) (/ x (* 24 60 60))) *time*))
         ;; (time-month (mapcar (lambda (x) (/ x (* 32d0 24 60 60))) *time*))
         )
    ;; (vgplot:figure)
    (vgplot:title "S_{xx} over height")
    (vgplot:xlabel "Time (s)")
    (vgplot:ylabel "Displacment (m)")
    (vgplot:axis (list 0 (reduce #'max time-day)
                       0 (reduce #'max *x-pos*)))
    (vgplot:plot time-day *x-pos* "MPM"
                 (aops:each (lambda (x) (* x 32d0)) (lisp-stat:column df 'time)) (lisp-stat:column df 'disp) "stokes"
                 )))
(defun plot-disp ()
  (let* ((df (lisp-stat:read-csv
	            (uiop:read-file-string #P"stokes.csv")))
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



(let* ((mp (aref (cl-mpm:sim-mps *sim*) 0))
       (stress (cl-mpm/particle::mp-undamaged-stress mp))
       (de (cl-mpm/particle::mp-elastic-matrix mp))
       (strain (cl-mpm/particle:mp-strain mp)))
  (pprint stress)
  (pprint strain)
  ;(pprint (magicl:linear-solve de stress))
  (pprint (magicl:linear-solve de stress))
  )
