(defpackage :cl-mpm/examples/tpb
  (:use :cl))
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)
(in-package :cl-mpm/examples/tpb)

;; (ql:quickload :magicl)
;; (pushnew :cl-mpm-pic *features*)
;; (delete :cl-mpm-pic *features*)
;; (asdf:compile-system :cl-mpm :force T)
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(defun cl-mpm/damage::length-localisation (local-length local-length-damaged damage)
  ;; (+ (* local-length (- 1d0 damage)) (* local-length-damaged damage))
  (* local-length (max (sqrt (- 1d0 damage)) 1d-10))
  ;; local-length
  )

(defmethod cl-mpm::update-node-forces ((sim cl-mpm::mpm-sim))
  (with-accessors ((damping cl-mpm::sim-damping-factor)
                   (mass-scale cl-mpm::sim-mass-scale)
                   (mesh cl-mpm::sim-mesh)
                   (dt cl-mpm::sim-dt))
      sim
    ;; (break)
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (node)
       ;; (cl-mpm::calculate-forces-cundall-conservative node damping dt mass-scale)
       ;; (cl-mpm::calculate-forces-cundall node damping dt mass-scale)
       (cl-mpm::calculate-forces node damping dt mass-scale)
       ))))

(declaim (notinline plot))
(defun plot-domain (sim)
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :deformed
   ;; :colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :zz))
   ;; :colour-func #'cl-mpm/particle::mp-damage
   ;; :colour-func (lambda (mp)
   ;;                (* (cl-mpm/particle::mp-mass mp)
   ;;                   (cl-mpm/fastmaths::mag-squared (cl-mpm/particle::mp-velocity mp)))
   ;;                )
   :colour-func #'cl-mpm/particle::mp-index
   ;:colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-damage-ybar mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-strain-plastic-vm mp))
   )
  ;; (vgplot:format-plot t "replot (~f*x + ~f)~%" 0d0 (+ 0.2d0 *target-displacement*))
  )
(defun plot (sim)
  ;; (cl-mpm/plotter:simple-plot
  ;;  *sim*
  ;;  :plot :deformed
  ;;  ;; :colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :zz))
  ;;  :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp))
  ;;  ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-damage-ybar mp))
  ;;  ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-strain-plastic-vm mp))
  ;;  )
  ;; (vgplot:format-plot t "replot (~f*x + ~f)~%" 0d0 (+ 0.2d0 *target-displacement*))
  (plot-load-disp)
  ;; (plot-time-disp)
  )
(defun stop ()
  (setf *run-sim* nil)
  (setf cl-mpm/dynamic-relaxation::*run-convergance* nil))

(defun rectangle-sdf (position size)
  (lambda (pos)
    (let* (
           (pos (magicl:from-list (list
                                        (magicl:tref pos 0 0)
                                        (magicl:tref pos 1 0)
                                        ) '(2 1) :type 'double-float))
           (position (magicl:from-list position '(2 1) :type 'double-float))
           (dist-vec (magicl:.- (magicl:map! #'abs (magicl:.- pos position))
                                (magicl:from-list size '(2 1) :type 'double-float))))

      (+ (sqrt (magicl::sum
                (magicl:map! (lambda (x) (* x x))
                             (magicl:map! (lambda (x) (max 0d0 x)) dist-vec))))
         (min (max (magicl:tref dist-vec 0 0)
                   (magicl:tref dist-vec 1 0)
                   ) 0d0)))))

(defun apply-pullout (sim load-mps push-rate)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (loop for mp in load-mps
          do
             (cl-mpm::iterate-over-neighbours
              mesh
              mp
              (lambda (mesh mp node svp grad fsvp fgrad)
                (with-accessors ()
                    mp
                  (with-accessors ((pos cl-mpm/mesh::node-position)
                                   (vel cl-mpm/mesh::node-velocity)
                                   (active cl-mpm/mesh::node-active)
                                   (acc cl-mpm/mesh::node-acceleration)
                                   )
                      node
                    (when active
                      (setf (magicl:tref vel 1 0) push-rate)))))))))

(defparameter *current-load* 0d0)
(defun apply-force (sim load-mps push-rate)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (loop for mp in load-mps
          do
             (setf (magicl:tref (cl-mpm/particle:mp-body-force mp) 1 0)
                   (* *current-load* (/ 1d0 (* (cl-mpm/particle:mp-volume mp) (length load-mps))))))
    ))

(defun energy-norm (sim)
  (/ (loop for mp across (cl-mpm:sim-mps *sim*)
          sum (magicl:norm (cl-mpm/particle:mp-velocity mp))) (length (cl-mpm:sim-mps *sim*))))

(defun get-disp (load-mps)
  ;; (* *t* *tip-velocity*)
  (- (/ (loop for mp in load-mps
              sum (+
                   (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)
                   (* 0.5d0 (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0))
                   )) (length load-mps))
     *initial-surface*
     ))

(defun get-force-mps (sim load-mps)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (let ((force 0d0))
      (loop for mp in load-mps
            do
               (cl-mpm::iterate-over-neighbours
                mesh
                mp
                (lambda (mesh mp node svp &rest args)
                  (incf force
                        (* svp
                           (cl-mpm/fastmaths::mag (cl-mpm/mesh::node-force node))
                           )
                        ))))
      (/ force (length load-mps))
      )))

(defun get-reaction-force (load-nodes)
  ;; (cl-mpm/fastmaths::mag (cl-mpm/mesh::node-force (nth 0 load-nodes)))
  (loop for mp in load-nodes
        sum
        ;; (cl-mpm/fastmaths::mag (cl-mpm/mesh::node-force mp))
        (/ (- (magicl:tref (cl-mpm/mesh::node-force mp) 1 0)) 1
          ; mps
           )
        ))
(defparameter *target-displacement* 0d0)
(defun apply-disp-penalty (sim load-mps)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (let* ((penalty 1d5)
           (displacement *target-displacement*)
           (pos
             (get-disp load-mps))
           (force (* penalty (- displacement pos))))
      (incf *current-load* force)
      (loop for mp in load-mps
            do
               (setf (magicl:tref (cl-mpm/particle:mp-body-force mp) 1 0)
                     (* force (/ 1d0 (* (cl-mpm/particle:mp-volume mp) (length load-mps))))
                     )))))

                                        ;(defparameter *tip-velocity* -0.02d-3)
(defparameter *tip-velocity* -0.000d-3)

(declaim (notinline setup-test-column))
(defun setup-test-column (size block-size offset &optional (e-scale 1) (mp-scale 1))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               ;; 'cl-mpm::mpm-sim-usf
               ;; :sim-type 'cl-mpm/damage::mpm-sim-usl-damage
               :sim-type 'cl-mpm/damage::mpm-sim-damage
               ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 2.2d3)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
         )
    (declare (double-float h density))
    (progn
      (let* ()
        (let* (;(crack-scale 7d0)
               (crack-scale 1.0d0)
               ;; (length-scale 5.3d-3)
               ;; (length-scale 20d-3) ;;65.15
               ;; (length-scale 10d-3) ;;47.56
               (length-scale (* 5.4d-3 (sqrt 7))) ;;
               (gf 48d0)
               (init-stress 3.45d6)
               (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf
                                                                        (/ length-scale 1d0) init-stress 15.3d9))
               )
          (format t "Estimated ductility ~F~%" ductility)
          (format t "Actual local length ~F~%" (* crack-scale length-scale))
          (format t "Length/Mesh res ~F~%" (/ (* crack-scale length-scale) (* 2d0 h-x)))
          (cl-mpm/damage::estimate-ductility-jirsek2004 (* 48d0 0.5d0) (* length-scale (sqrt 7)) 3.45d6 15.3d9)
          (setf (cl-mpm:sim-mps sim)
                (cl-mpm/setup::make-mps-from-list
                 (append
                  (cl-mpm/setup::make-block-mps-list
                   offset
                   block-size
                   (mapcar (lambda (e) (round (* e mp-scale) h-x)) block-size)
                   density
                   ;; 'cl-mpm/particle::particle-concrete
                   'cl-mpm/particle::particle-limestone-delayed
                   :E 15.3d9
                   ;; :E 18d9
                   :nu 0.15d0
                   :fracture-energy (* 48d0 1d0)
                   :initiation-stress 3.45d6
                   ;;Material parameter
                   :local-length length-scale
                   :compression-ratio 8d0

                   ;:ductility 7.1d0
                   ;:ductility 10.16d0
                   :ductility ductility
                   ;; :ductility 3.0d0
                   ;:ductility 5.0d0
                   ;; :ductility 2.0d0
                   :delay-time 0.1d0

                   :critical-damage 1.000d0
                   :enable-plasticity nil
                   :c 0d0
                   :phi 0d0
                   :psi 0d0
                   ;; :local-length-damaged 0.01d0
                   :gravity 0.0d0
                   :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 0d0 0d0))
                   )
                  ;; impactors
                  )
                 ))))
      ;; (calculate-ductility-param 18d9  48d0 5.3d-3 3.4d6)
      ;calculate-ductility-param (E Gf l-c f-t)
      (setf (cl-mpm:sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-velocity-algorithm sim) :PIC)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (setf (cl-mpm::sim-mass-filter sim) 1d-10)
      (let ((ms 1d0))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        (setf (cl-mpm:sim-damping-factor sim)
              (* 1d-1
                 (cl-mpm/setup:estimate-critical-damping sim))))


      (dotimes (i 0)
        (dolist (dir (list :x :y))
          (cl-mpm::split-mps-criteria
           sim
           (lambda (mp h)
             (when
                 (and
                  (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                     (* (first block-size) 0.10)))
               dir)))))

      (let ((dt-scale 0.5d0))
        (setf
         (cl-mpm:sim-dt sim)
         (* dt-scale 
            (cl-mpm/setup:estimate-elastic-dt sim))))

      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))

      (let* ((crack-width 1d-2)
             ;; (crack-left (- (+ (first offset)  (* 0.5d0 (first block-size))) crack-width))
             ;; (crack-right (+ (+ (first offset) (* 0.5d0 (first block-size))) crack-width))
             (crack-left 0d0)
             (crack-right crack-width)
             (above-crack
               (loop for mp across (cl-mpm:sim-mps sim)
                     when
                     (and
                      (not (= (cl-mpm/particle::mp-index mp) 1))
                      (>= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) crack-left)
                      (<= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) crack-right))
                     collect mp)
               )
             (max-pos (loop for mp in above-crack
                            maximize (magicl:tref(cl-mpm/particle:mp-position mp) 1 0)))
             (min-x-pos (loop for mp in above-crack
                            minimize (magicl:tref(cl-mpm/particle:mp-position mp) 0 0)))
             )
        (defparameter *terminus-mps*
          (loop for mp in above-crack
                when (and
                      ;; (= min-x-pos (magicl:tref
                      ;;             (cl-mpm/particle:mp-position mp)
                      ;;             0 0))
                      (= max-pos (magicl:tref
                                      (cl-mpm/particle:mp-position mp)
                                      1 0)))
                  collect mp)))

      (setf *terminus-mps* (list (nth 0 *terminus-mps*)))
      ;; (setf *terminus-mps* (mapcar (lambda (i) (nth i *terminus-mps*))
      ;;                              (list (floor (- (length *terminus-mps*) 1) 2))))
      ;; (when (> (length *terminus-mps*) 2)
      ;;   (setf *terminus-mps* (mapcar (lambda (i) (nth i *terminus-mps*))
      ;;                                (list (floor (- (length *terminus-mps*) 1) 2)
      ;;                                      (+ (floor (- (length *terminus-mps*) 1) 2) 1)))))
      (loop for mp in *terminus-mps*
            do (setf (cl-mpm/particle::mp-index mp) 1))

      (let ((left-node-pos
              (list
               (round (first offset) h-x)
               (round (second offset) h-x)
               0))
            (right-node-pos
              (list
               (round (+ (first offset) (first block-size)) h-x)
               (round (second offset) h-x)
               0)
              ))

        (defparameter *fixed-nodes* (mapcar (lambda (id) (cl-mpm/mesh:get-node (cl-mpm:sim-mesh sim)
                                                                                     id))
                                                  (list
                                                   ;; left-node-pos
                                                   right-node-pos)))

        (format t "Fixed node ~A ~%" left-node-pos)
        (format t "Roller node ~A ~%" right-node-pos)

        (setf (cl-mpm:sim-bcs sim)
              (cl-mpm/bc::make-bcs-from-list
               (append
                (cl-mpm/bc::make-outside-bc-var-list
                 (cl-mpm:sim-mesh sim)
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil nil)))
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil nil)))
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil nil)))
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0))))

                (list
                 (cl-mpm/bc::make-bc-fixed
                  right-node-pos
                  '(nil 0 nil)))))))

      (defparameter *initial-surface*
        (loop for mp across (cl-mpm:sim-mps sim)
              when (not (= 1 (cl-mpm/particle::mp-index mp)))
              maximize (magicl:tref
                        (magicl:.+ (cl-mpm/particle:mp-position mp)
                                   (magicl:scale (cl-mpm/particle::mp-domain-size mp) 0.5d0)
                                   )
                        1 0
                        )))

      (defparameter *penalty*
        (cl-mpm/penalty::make-bc-penalty-distance-point
         sim
         (cl-mpm/utils:vector-from-list '(0d0 -1d0 0d0))
         (cl-mpm/utils:vector-from-list (list 0d0
                                              *initial-surface*
                                              0d0))
         h-x
         (* 15.3d9 1d1)
         0d0
         0d0
         )
        )
      (defparameter *current-inc* 0d0)
      (format t "~A~%" h-x)
      (setf (cl-mpm::sim-bcs-force-list sim)
            (list
             (cl-mpm/bc:make-bcs-from-list
              (list
               (cl-mpm/bc::make-bc-closure
                nil
                (lambda ()
                  (let ((delta (- *target-displacement* *current-inc*)))
                    (cl-mpm/penalty::bc-increment-center
                     *penalty*
                     (cl-mpm/utils:vector-from-list (list 0d0 delta 0d0)))
                    (setf *current-inc* *target-displacement* ))
                  ;; (setf
                  ;;  (cl-mpm/penalty::bc-penalty-datum *penalty*)
                  ;;  (- (+ *initial-surface* *target-displacement*)))
                  )
                )
               *penalty*
               ;; *floor-bc*
               ;; (cl-mpm/bc::make-bc-closure
               ;;  nil
               ;;  (lambda ()
               ;;    (with-accessors ((mesh cl-mpm:sim-mesh)
               ;;                     (dt cl-mpm::sim-dt))
               ;;        sim
               ;;      (let ((datum (* -1d0 (+ *initial-surface* *target-displacement*)))
               ;;            (normal (cl-mpm/utils:vector-from-list  '(0d0 -1d0 0d0))))
               ;;        ;; (format t "Datum: ~F~%" datum)
               ;;        (cl-mpm/penalty::apply-displacement-control-mps mesh (coerce *terminus-mps* 'vector )
               ;;                                         dt
               ;;                                         normal
               ;;                                         datum
               ;;                                         ;; (* density 1d6)
               ;;                                         (* 15.3d9 1)
               ;;                                         ;; (* 0d0 density 1d10)
               ;;                                         0d0))
               ;;      )))
               ))))

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


(defparameter *sim* nil)
(defparameter *run-sim* t)
(defparameter *t* 0)
(defparameter *sim-step* 0)
(defparameter *refine* (/ 1d0 1d0))
(let ((refine (uiop:getenv "REFINE")))
  (when refine
    (setf *refine* (parse-integer (uiop:getenv "REFINE")))
    ))

(declaim (notinline setup))
(defun setup (&key (undercut 0d0) (refine 1d0) (mps 2))
  (let* ((mesh-size (/ 0.0102 (* refine 1.0)))
         (mps-per-cell mps)
         (shelf-height 0.102d0)
         (shelf-length (* shelf-height 2.0))
         ;; (shelf-length 0.225d0)
         (domain-length (+ shelf-length
                           (* shelf-height 0.5d0)
                           ))
         (offset (list
                  0d0
                  (* shelf-height 0.5)
                       ))
         )
    (defparameter *sim*
      (setup-test-column (list domain-length
                               (* shelf-height 2.0))
                         (list shelf-length shelf-height)
                         offset
                         (/ 1d0 mesh-size) mps-per-cell))
    (let ((cut-depth (* 0.4d0 shelf-height))
          (cut-width
            ;; (/ 0.0102 1.0d0)
            mesh-size
                     ))
      (format t "Crack width:~F~%" (* 1d0 cut-width))
      (cl-mpm/setup::remove-sdf
       *sim*
       (rectangle-sdf
        (list
         0d0
         (+ (second offset) 0d0))
        (list
         cut-width
         ;; 0.0102
         ;; 10.0d-3
         ;2.5d-3
         ;; 1.33d-3
         ;; 10d-3
         ;; mesh-size
         cut-depth
         )))

      ;; (cl-mpm/setup::damage-sdf
      ;;  *sim*
      ;;  (rectangle-sdf
      ;;   (list
      ;;    0d0
      ;;                                   ;(+ (first offset) (* 0.5d0 shelf-length))
      ;;    (+ (second offset) 0d0)
      ;;    )
      ;;   (list
      ;;    ;; 3.0d-3
      ;;    10d-3
      ;;    1000d0
      ;;    ))
      ;;  0.9d0
      ;;  )
      )
    (format t "Total weight ~F~%"
            (loop for mp across (cl-mpm:sim-mps *sim*)
                  sum (* 9.8d0 (cl-mpm/particle:mp-mass mp))))

    (defparameter *current-load* 0d0)
    ;; (loop for mp across (cl-mpm:sim-mps *sim*)
    ;;       do
    ;;          (setf (cl-mpm/particle:mp-damage mp) (random 0.1d0)))
    ;; (cl-mpm/setup::damage-sdf
    ;;  *sim*
    ;;  (lambda (p)
    ;;    (cl-mpm/setup::line-sdf (magicl:from-list (list (magicl:tref p 0 0)
    ;;                                                    (magicl:tref p 1 0))
    ;;                                              '(2 1))
    ;;                            (list (- shelf-length shelf-height) shelf-height)
    ;;                            (list shelf-length soil-boundary)
    ;;                            10d0
    ;;                            )) 0.8d0)
    ;(let ((sdf
    ;        (lambda (p)
    ;          (cl-mpm/setup::line-sdf (magicl:from-list (list (magicl:tref p 0 0)
    ;                                                          (magicl:tref p 1 0))
    ;                                                    '(2 1))
    ;                                  (list (- shelf-length shelf-height) shelf-height)
    ;                                  (list shelf-length 0d0)
    ;                                  20d0
    ;                                  ))
    ;        ))
    ;  (loop for mp across (cl-mpm:sim-mps *sim*)
    ;        do (with-accessors ((pos cl-mpm/particle:mp-position)
    ;                            (damage cl-mpm/particle:mp-damage)) mp
    ;             (when (>= 0 (funcall sdf pos))
    ;               (setf damage (min 1d0 (max 0d0 (coerce (* (funcall sdf pos) -0.1d0) 'double-float)))))
    ;             )))


      (let* ((crack-width 1d-1)
             (crack-left 0d0)
             (crack-right crack-width)
             (crack-left-pos
               (loop for mp across (cl-mpm:sim-mps *sim*)
                     minimize (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
             (crack-height-pos
               (loop for mp across (cl-mpm:sim-mps *sim*)
                     when
                     (= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) crack-left-pos)
                     minimize (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)))
             (bottom-crack
               (nth 0
                    (loop for mp across (cl-mpm:sim-mps *sim*)
                          when
                          (and(= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) crack-left-pos)
                              (= (magicl:tref (cl-mpm/particle:mp-position mp) 1 0) crack-height-pos))
                          collect mp)))
             (slice-crack
               (loop for mp across (cl-mpm:sim-mps *sim*)
                     when
                     (= (magicl:tref (cl-mpm/particle:mp-position mp) 1 0) (magicl:tref (cl-mpm/particle:mp-position bottom-crack) 1 0))
                     collect mp))
             )
        (defparameter *mps-slice*
          slice-crack
          ))


    )
  (defparameter *data-full-time* '(0d0))
  (defparameter *data-full-load* '(0d0))
  (defparameter *data-full-reaction* '(0d0))
  (defparameter *target-displacement* 0d0)
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (format t "DOFs: ~D~%" (* 2d0 (array-total-size (cl-mpm/mesh:mesh-nodes (cl-mpm:sim-mesh *sim*)))))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defgeneric estimate-energy-crit (sim))
(defmethod estimate-energy-crit ((sim cl-mpm::mpm-sim))
  ;; (/)
  (loop for mp across (cl-mpm:sim-mps sim)
        summing (* (cl-mpm/particle::mp-mass mp)
                      ;; (cl-mpm/particle::mp-damage mp)
                      (cl-mpm/fastmaths::mag-squared (cl-mpm/particle::mp-velocity mp))))
  ;; (cl-mpm::sim-mass-scale sim)
  ;; (* (loop for mp across (cl-mpm:sim-mps *sim*)
  ;;          summing (cl-mpm/particle::mp-mass mp))
  ;;    ;; (cl-mpm::sim-mass-scale sim)
  ;;    )
  
  )

(defparameter *data-force* '())
(defparameter *data-displacement* '(0d0))
(defparameter *data-load* '(0d0))
(defparameter *data-time* '(0d0))
(defparameter *data-node-load* '(0d0))
(defparameter *target-displacement* 0d0)
(defparameter *data-averaged* t)

(defparameter *data-full-energy* '(0d0))
(defparameter *data-full-time* '(0d0))
(defparameter *data-full-load* '(0d0))

(defun run ()
  (vgplot:close-all-plots)
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)
  (defparameter *data-force* '())
  (defparameter *data-displacement* '(0d0))
  (defparameter *data-load* '(0d0))
  (defparameter *data-node-load* '(0d0))
  (defparameter *data-mp-load* '(0d0))
  (defparameter *data-time* '(0d0))
  (defparameter *target-displacement* 0d0)
  (defparameter *data-full-time* '(0d0))
  (defparameter *data-full-load* '(0d0))
  (defparameter *data-y* '(0d0))
  (defparameter *data-ybar* '(0d0))
  ;; (defmethod cl-mpm::update-node-forces ((sim cl-mpm::mpm-sim))
  ;;   (with-accessors ((damping cl-mpm::sim-damping-factor)
  ;;                    (mass-scale cl-mpm::sim-mass-scale)
  ;;                    (mesh cl-mpm::sim-mesh)
  ;;                    (dt cl-mpm::sim-dt))
  ;;       sim
  ;;     (cl-mpm::iterate-over-nodes
  ;;      mesh
  ;;      (lambda (node)
  ;;        (cl-mpm::calculate-forces node damping dt mass-scale)))))

  (let ((ms 1d6))
    (setf (cl-mpm::sim-mass-scale *sim*) ms)
    (* 1d-2
       (cl-mpm::sim-mass-scale *sim*)
       (cl-mpm/setup::estimate-critical-damping *sim*)))

  (with-open-file (stream (merge-pathnames "output/disp.csv") :direction :output :if-exists :supersede)
    (format stream "disp,load,load-mps~%"))

  (let* ((load-steps 100)
         (target-time (* 1d0 (/ 50 load-steps)))
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 1d0)
         (disp-step (/ -0.2d-3 load-steps))
         )

    (setf cl-mpm/penalty::*debug-force* 0d0)
    (setf cl-mpm/penalty::*debug-force-count* 0d0)
    (setf cl-mpm/damage::*enable-reflect-x* t)
    (time (cl-mpm::update-sim *sim*))
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                    (format t "CFL dt estimate: ~f~%" dt-e)
                    (format t "CFL step count estimate: ~D~%" substeps-e)
                    (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
    (when (slot-exists-p *sim* 'cl-mpm/damage::delocal-counter-max)
      (setf (cl-mpm/damage::sim-damage-delocal-counter-max *sim*) substeps))
    (setf (cl-mpm::sim-enable-damage *sim*) t)
    (time (loop for steps from 0 to load-steps
                while *run-sim*
                do
                   (progn
                     ;; (when (> *target-displacement* 0.05d-3)
                       ;; (setf (cl-mpm::sim-enable-damage *sim*) t)
                       ;; )
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)


                     (let ((average-force 0d0)
                           (average-disp 0d0)
                           (average-reaction 0d0))
                       (time
                        (dotimes (i substeps);)
                          (push
                           *t*
                           *data-full-time*)
                          (cl-mpm/penalty::reset-load *penalty*)
                          (cl-mpm::update-sim *sim*)

                          (incf average-disp (/ *target-displacement* substeps))
                          (incf average-force (/ (* (cl-mpm/penalty::resolve-load *penalty*)
                                                    2d0
                                                    ;; 1d-3
                                                    ;; (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))
                                                    )
                                                 substeps))
                          (incf *target-displacement* (/ disp-step substeps))
                          (setf *t* (+ *t* (cl-mpm::sim-dt *sim*))))
                        )
                       (setf average-reaction (* (get-reaction-force *fixed-nodes*) 2d0))
                       (setf average-reaction average-force)

                       (push
                        average-disp
                        *data-displacement*)
                       (push
                        average-force
                        *data-mp-load*)
                       (push
                        average-reaction
                        *data-load*)
                       (let ((ymax
                               (loop for mp across (cl-mpm:sim-mps *sim*)
                                     maximize (cl-mpm/particle::mp-damage-y-local mp)))
                             (ybarmax
                               (loop for mp across (cl-mpm:sim-mps *sim*)
                                     maximize (cl-mpm/particle::mp-damage-ybar mp)))
                             )
                         (format t "Ybar ratio: ~F~%" (/ ybarmax (max ymax 1d0)))
                         (push
                          ymax
                          *data-y*)
                         (push
                          ybarmax
                          *data-ybar*))
                       ;; (when *data-averaged*
                       ;;   (push
                       ;;    ;; (- *t*)
                       ;;    (get-disp *terminus-mps*)
                       ;;    *data-displacement*)

                       ;;   (push
                       ;;    (get-reaction-force *fixed-nodes*)
                       ;;    *data-load*)

                       ;;   (push
                       ;;    ;; (get-force-mps *sim* *terminus-mps*)
                       ;;    (/ cl-mpm/penalty::*debug-force*
                       ;;       (max 1 cl-mpm/penalty::*debug-force-count*))
                       ;;    *data-mp-load*))
                       ;; (incf *target-displacement* disp-step)


                       ;; (format t "Target load: ~f~%" (* *target-displacement* 20d9))
                       ;; (format t "Current load: ~f~%" (* (get-disp *terminus-mps*) 1d9))
                       ;; (format t "Node Load: ~f~%" (get-reaction-force *fixed-nodes*))
                       ;; (format t "Pen load: ~f~%" average-force)
                       
                       (with-open-file (stream (merge-pathnames "output/disp.csv") :direction :output :if-exists :append)
                         (format stream "~f,~f,~f~%"
                                 average-disp
                                 average-reaction
                                 average-force
                                 )))

                     ;; (setf cl-mpm/penalty::*debug-force*
                     ;;       (/ cl-mpm/penalty::*debug-force* substeps)
                     ;;       )

                     ;; (print cl-mpm/penalty::*debug-force*)

                     ;; (when (>= steps 5)
                     ;;   (setf (cl-mpm::sim-enable-damage *sim*) t)
                     ;;   (incf *current-load* -1d2)
                     ;;   (apply-force *sim* *terminus-mps* *current-load*))
                     ;; (when (reduce (lambda (a b) (and a b))
                     ;;               (mapcar (lambda (mp) (>= (cl-mpm/particle::mp-damage mp) 1d0)) *terminus-mps*))
                     ;;   (defparameter *run-sim* nil))
                     (format t "Target: ~f - Current: ~f Error: ~f - energy ~F~%"
                             *target-displacement*
                             (get-disp *terminus-mps*)
                             (* 100d0 (/ (- *target-displacement* (get-disp *terminus-mps*)) *target-displacement*))
                             (energy-norm *sim*))
                     ;; (format t "Debug load: ~f~%" cl-mpm/penalty::*debug-force*)
                     ;; (format t "Debug load count: ~f~%" cl-mpm/penalty::*debug-force-count*)

                     ;; (incf *target-displacement* -0.01d-3)

                     ;; (incf *current-load* (* -1d3 1d0))

                     ;; (loop for mp across (cl-mpm:sim-mps *sim*)
                     ;;       do
                     ;;          (when (>= (cl-mpm/particle:mp-damage mp) 1d0)
                     ;;            (let ((ms 1d2))
                     ;;              (setf (cl-mpm::sim-mass-scale *sim*) ms
                     ;;                    ;; target-time 1d0
                     ;;                    (cl-mpm:sim-damping-factor *sim*) (* 1d-8 ms)
                     ;;                    ))))

                     (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       (format t "CFL dt estimate: ~f~%" dt-e)
                       (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf substeps substeps-e))
                     ;; (setf (cl-mpm:sim-damping-factor *sim*)
                     ;;       (* (cl-mpm:sim-damping-factor *sim*) (expt 1d-3 1/40)))

                     (incf *sim-step*)
                     (plot *sim*)
                     (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                        :terminal "png size 1920,1080"
                                        )

                     (swank.live:update-swank)
                     ))))
  (format t "Estimated GF: ~F~%" (estimate-gf-no-plot))
  (vgplot:figure)
  (plot-load-disp)
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*))
(defun run-static ()
  (vgplot:close-all-plots)
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)
  (defparameter *data-force* '())
  (defparameter *data-displacement* '(0d0))
  (defparameter *data-load* '(0d0))
  (defparameter *data-node-load* '(0d0))
  (defparameter *data-mp-load* '(0d0))
  (defparameter *data-time* '(0d0))
  (defparameter *target-displacement* 0d0)
  (defparameter *data-full-time* '(0d0))
  (defparameter *data-full-load* '(0d0))
  (defparameter *data-y* '(0d0))
  (defparameter *data-ybar* '(0d0))

  (loop for mp across (cl-mpm:sim-mps *sim*)
        do (change-class mp 'cl-mpm/particle::particle-limestone))

  (with-open-file (stream (merge-pathnames "output/disp.csv") :direction :output :if-exists :supersede)
    (format stream "disp,load,load-mps~%"))

  (let* ((dt (cl-mpm:sim-dt *sim*))
         (load-steps 50)
         (dt-scale 0.5d0)
         (disp-step (/ -0.2d-3 load-steps))
         )


    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (setf (cl-mpm::sim-damping-factor *sim*) 
      (* 1d-2 (cl-mpm/setup::estimate-critical-damping *sim*)))
    (setf cl-mpm/penalty::*debug-force* 0d0)
    (setf cl-mpm/penalty::*debug-force-count* 0d0)
    (setf (cl-mpm/damage::sim-damage-delocal-counter-max *sim*) 0)
    (setf cl-mpm/damage::*enable-reflect-x* t)

    (time (loop for steps from 0 below load-steps
                while *run-sim*
                do
                   (progn
                       (setf (cl-mpm::sim-enable-damage *sim*) t)
                       ;; )
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)

                     (setf (cl-mpm::sim-enable-damage *sim*) nil)
                     (incf *target-displacement* disp-step)

                     (time
                      (progn
                        (handler-bind ((cl-mpm/dynamic-relaxation::non-convergence-error
                                         (lambda (c)
                                           (format t "System failed to converge - end sim")
                                           (setf *run-sim* nil)
                                           (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_fail_~5,'0d.vtk" *sim-step*)) *sim*)
                                           (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_fail_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                                           (error "Conv fail")
                                           ;(invoke-restart 'cl-mpm/dynamic-relaxation::continue)
                                           )))
                          (cl-mpm/dynamic-relaxation::converge-quasi-static
                           *sim*
                           :energy-crit 1d-2
                           :oobf-crit 1d-1
                           :dt-scale dt-scale
                           :conv-steps 200
                           :substeps 50
                           :post-iter-step
                           (lambda (step ke fnorm)
                             (let ((av (get-disp *terminus-mps*)))
                               ;; (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_conv_~5,'0d.vtk" step)) *sim*)
                               ;; (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_conv_nodes_~5,'0d.vtk" step)) *sim*)
                               (format t "Conv disp - Target: ~E - Current: ~E - Error: ~E~%"
                                       *target-displacement*
                                       av
                                       (/ (abs (- *target-displacement* av )) *target-displacement*)
                                       )))))
                        (cl-mpm/damage::calculate-damage *sim*)
                        )
                      )

                     (let ((average-force 0d0)
                           (average-disp 0d0)
                           (average-reaction 0d0));)
                        ;)
                       (setf average-disp (get-disp *terminus-mps*))
                       (setf average-force (* cl-mpm/penalty::*debug-force* 2.0d0))
                       (setf average-reaction (* (get-reaction-force *fixed-nodes*) 2d0))
                       ;; (incf *target-displacement* -0.01d-3)
                       (push
                        average-disp
                        *data-displacement*)
                       (push
                        average-force
                        *data-mp-load*)
                       (push
                        average-reaction
                        *data-load*)

                       (with-open-file (stream (merge-pathnames "output/disp.csv") :direction :output :if-exists :append)
                         (format stream "~f,~f~%"
                                 average-disp
                                 average-reaction
                                 average-force
                                 )))

                     (format t "Target: ~f - Current: ~f Error: ~f - energy ~F~%"
                             *target-displacement*
                             (get-disp *terminus-mps*)
                             (* 100d0 (/ (- *target-displacement* (get-disp *terminus-mps*)) *target-displacement*))
                             (energy-norm *sim*))
                     (let ((damage-est
                             (lparallel:pmap-reduce (lambda (mp)
                                                      (*
                                                       1d0
                                                       (cl-mpm/particle::mp-damage mp)))
                                                    #'+ (cl-mpm:sim-mps *sim*)
                                                    :initial-value 0d0)
                             ))
                       (format t "Total damage: ~E~%" damage-est))

                     ;; (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                     ;;   (format t "CFL dt estimate: ~f~%" dt-e)
                     ;;   (format t "CFL step count estimate: ~D~%" substeps-e)
                     ;;   (setf substeps substeps-e))

                     (incf *sim-step*)
                     (plot *sim*)
                     (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                        :terminal "png size 1920,1080"
                                        )

                     (swank.live:update-swank)
                     ))))
  (format t "Estimated GF: ~F~%" (estimate-gf-no-plot))
  (vgplot:figure)
  (plot-load-disp)
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*))

(declaim (notinline plot-time-disp))
(defun plot-time-disp ()
  ;; (vgplot:figure)
  (vgplot:xlabel "Displacement (mm)")
  (vgplot:ylabel "Load (N)")
  (vgplot:plot
   (mapcar (lambda (x) (* x -1d0)) *data-full-time*) (mapcar (lambda (x) (* x 0.1)) *data-full-load*) "mps"
   (mapcar (lambda (x) (* x -1d0)) *data-full-time*) (mapcar (lambda (x) (* x 0.1)) *data-full-reaction*) "nodes"
   (mapcar (lambda (x) (* x -1d0)) *data-full-time*) (mapcar (lambda (x) (* x 0.1)) (mapcar (lambda (x) (* x 1d6)) *data-full-energy*)) "energy"
   ))

;; (defun plot-load-disp ()
;;   ;; (vgplot:figure)
;;   (vgplot:xlabel "Displacement (mm)")
;;   (vgplot:ylabel "Load (N)")
;;   (vgplot:plot
;;    ;; (mapcar (lambda (x) (* x -1d3)) *data-displacement*) (mapcar (lambda (x) (* x 1.0)) *data-load*)
;;    (mapcar (lambda (x) (* x -1d0)) *data-time*) (mapcar (lambda (x) (* x 1.0)) *data-load*) "nodes"
;;    (mapcar (lambda (x) (* x -1d0)) *data-time*) (mapcar (lambda (x) (* x 1.0)) *data-mp-load*) "mps"
;;    ;; (mapcar (lambda (x) (* x -1d3)) *data-displacement*) (mapcar (lambda (x) (* x 1.0)) *data-mp-load*)
;;   ;; (vgplot:plot (mapcar (lambda (x) (* x 1d0)) *data-displacement*) *data-load*)
;;   ))
(defun plot-load-disp ()
  (let ((df (lisp-stat:read-csv
	           (uiop:read-file-string #P"example_data/tpb/load-disp.csv")))
        (fem (lisp-stat:read-csv
	           (uiop:read-file-string #P"example_data/tpb/load-disp-standard.csv")))
        )
    (vgplot:plot
     (lisp-stat:column df 'disp) (lisp-stat:column df 'load) "Experimental"
     (lisp-stat:column fem 'disp) (lisp-stat:column fem 'load) "FEM"
     ;; (mapcar (lambda (x) (* x -1d3)) *data-displacement*) *data-node-load* "node"
     (mapcar (lambda (x) (* x -1d3)) *data-displacement*) (mapcar (lambda (x) (* x 0.013)) *data-load*) "mpm-reaction"
     (mapcar (lambda (x) (* x -1d3)) *data-displacement*) (mapcar (lambda (x) (* x 0.013)) *data-mp-load*) "mpm-force"
     ;; (mapcar (lambda (x) (* x -1d3)) *data-displacement*) (mapcar (lambda (x) (* x 1d-4)) *data-y*) "mpm-y"
     ;; (mapcar (lambda (x) (* x -1d3)) *data-displacement*) (mapcar (lambda (x) (* x 1d-4)) *data-ybar*) "mpm-ybar"
     ;; (mapcar (lambda (x) (* x -1d3)) *data-displacement*) (mapcar (lambda (x) (* x -2d9)) *data-displacement*) "LE"
     )

    ;; (vgplot:format-plot t "set xrange [~f:~f]"
    ;;                     ;(* -1d3 (- 0.0000001d0 (reduce #'max *data-displacement*)))
    ;;                     0d0
    ;;                     (+ 1d-4 (* -1d3 (reduce #'min *data-displacement*)))
    ;;                     )
    ;; (vgplot:format-plot t "set yrange [~f:~f]"
    ;;                     (reduce #'min (mapcar #'min *data-load* *data-mp-load*))
    ;;                     (* 1.01 (reduce #'max (mapcar #'max *data-load* *data-mp-load*)))
    ;;                     )
    ;; (vgplot:axis (list nil nil nil nil))
    (vgplot:xlabel "Displacement (mm)")
    (vgplot:ylabel "Load (N)")
    )
  )

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

(defun simple-test ()
  (setup)
  ;; (sb-profile:profile "CL-MPM")
  ;; (sb-profile:profile "CL-MPM/PARTICLE")
  ;; (sb-profile:profile "CL-MPM/MESH")
  ;; (sb-profile:profile "CL-MPM/SHAPE-FUNCTION")
  ;; (sb-profile:reset)
  (time
   (dotimes (i 100)
         (cl-mpm::update-sim *sim*)))
  ;; (sb-profile:report)
  )

(defun simple-test ()
  (setup)
  (time
   (dotimes (i 100)
     (cl-mpm::update-stress (cl-mpm:sim-mesh *sim*) (cl-mpm:sim-mps *sim*) 1d0 nil)))
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
;; (lparallel:end-kernel)
;; (sb-ext::exit)
;; (uiop:quit)
;; (loop for mp across (cl-mpm:sim-mps *sim*)
;;       maximize (magicl:tref
;;                 (magicl:.+ (cl-mpm/particle:mp-position mp)
;;                                       (magicl:scale (cl-mpm/particle::mp-domain-size mp) 0.5d0)
;;                                       )
;;                 1 0
;;                 ))

(defun plot-stress-damage ()
  (vgplot:close-all-plots)
  (let* ((mp (aref (cl-mpm:sim-mps *sim*))))
    (with-accessors ((damage cl-mpm/particle:mp-damage)
                     (E cl-mpm/particle::mp-e)
                     (Gf cl-mpm/particle::mp-Gf)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (length cl-mpm/particle::mp-local-length)
                     ) mp
      (let* ((stress (loop for x from 1d4 to 10d6 by 1d4
                           collect x))
             (damage (mapcar (lambda (stress)
                               ;(cl-mpm/damage::brittle-concrete-d stress E Gf length init-stress)
                               ;; (cl-mpm/damage::brittle-concrete-linear-d stress E Gf length init-stress)
                               (cl-mpm/damage::damage-response-exponential stress E Gf length init-stress)
                               )
                             stress)))
        (vgplot:plot stress damage)))))

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

(defun vec-test ()
  (let ((iters 10000000))
    (time
     (dotimes (i iters)
       (let ((a (make-array 3 :element-type 'double-float :initial-element 0d0))
             (b (make-array 3 :element-type 'double-float :initial-element 0d0))
             )
         (loop for i fixnum from 0 to 2
               do (incf (aref a i) (aref b i))))
       ))
    (time
     (dotimes (i iters)
       (let ((a (static-vectors:make-static-vector 3 :element-type 'double-float :initial-element 0d0))
             (b (static-vectors:make-static-vector 3 :element-type 'double-float :initial-element 0d0))
             )
         (loop for i fixnum from 0 to 2
               do (incf (aref a i) (aref b i))))
       ))
    )
  )


(defun test ()
  (let ((iters 1000))
    (setf magicl:*default-allocator* #'magicl::c-allocator)
    (setup)
    (time
     (dotimes (i iters)
       (cl-mpm::update-sim *sim*)))
    (setf magicl:*default-allocator* #'magicl::lisp-allocator)
    (setup)
    (time
     (dotimes (i iters)
       (cl-mpm::update-sim *sim*))))
  )

(defun plot-crack-inter ()
  (vgplot:close-all-plots)
  (let* ((mesh (cl-mpm:sim-mesh *sim*))
         (mps-taken *mps-slice*)
         (x (loop for mp in mps-taken collect (magicl:tref (cl-mpm/particle::mp-position mp) 0 0)))
         (damage (loop for mp in mps-taken collect  (cl-mpm/particle::mp-damage mp)))
         )
    (vgplot:plot x (mapcar (lambda (x) (cl-mpm/damage::weight-func-mps mesh (nth 0 mps-taken) x (cl-mpm/particle::mp-local-length x))) mps-taken) ";;with points")
    ))
(defun plot-crack-damage ()
  (vgplot:close-all-plots)
  (let* ((mesh (cl-mpm:sim-mesh *sim*))
         (mps-taken *mps-slice*)
         (x (loop for mp in mps-taken collect (magicl:tref (cl-mpm/particle::mp-position mp) 0 0)))
         (damage (loop for mp in mps-taken collect  (cl-mpm/particle::mp-damage mp)))
         )
    (vgplot:plot x damage ";;with points")
    ))

(defun plot-interaction ()
  (vgplot:close-all-plots)
  (let* ((length 1d0)
         (scaler 2d0)
         (bounds (* scaler length))
         (x (loop for x from (- bounds) to bounds by 0.01d0 collect x)))
    (vgplot:plot x (mapcar (lambda (x) (cl-mpm/damage::weight-func (* x x) length)) x))
    ))



(defun estimate-gf-no-plot ()
  (let* (
         (strain (mapcar (lambda (x) (* x -1d0)) *data-displacement*))
         (stress (mapcar (lambda (x) (* x 0.013)) *data-load*))
         ;; (strain (mapcar (lambda (x) x) *data-disp*))
         ;; (stress (mapcar (lambda (x) x) *data-load*))
         (dstrain (mapcar #'- (butlast strain) (cdr strain)))
         (av-stress (mapcar #'+ (butlast stress) (cdr stress)))
         (energy (* 0.5d0 (reduce #'+ (mapcar #'* dstrain av-stress)))))
    (* energy (/ 1d0
                 (* (* 0.102d0 0.6d0) 0.013)
                 ;; (/ *bar-length* 2)
                 ;; (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))
                 ))))

(defun test-slot ()
  (declare (optimize (debug 0) (safety 0) (speed 3) (space 0)))
  (let ((mp (make-instance 'cl-mpm/particle::particle :nd 2))
        (iters 1000000000)
        (y 0d0))
    (time
     (dotimes (i iters)
       (with-accessors ((volume cl-mpm/particle::mp-volume))
           mp
           (setf y volume))
       ))
    (time
     (dotimes (i iters)
       (with-accessors ((volume cl-mpm/particle::mp-mass))
           mp
         (setf y volume))
       ))
    ;; (time
    ;;  (dotimes (i iters)
    ;;    (with-slots ((volume cl-mpm/particle::volume))
    ;;        mp
    ;;      (setf volume 0))
    ;;    ))
    ;; (time
    ;;  (dotimes (i iters)
    ;;    (setf (sb-mop:standard-instance-access mp 0) 0)
    ;;    ))
    ))
(defun estimate-damage-zone (sim)
  (loop for mp across (cl-mpm:sim-mps sim)
        when (>= (cl-mpm/particle::mp-damage mp) 0.90d0)
        maximize (magicl:tref (cl-mpm/particle::mp-position mp) 0 0)))




(defparameter *calibration-ductility* 
  (list
   2.0042060175857923
   2.436951134024035 
   2.723975956151442 
   3.15683146675204  
   4.096506576732223 
   4.964535875342916 
   6.6310461502795786
   8.516799231656632 
   11.491811513007194
   15.339048071138008
   19.840259647069864
   24.996329394101643
   30.224817711639417
   35.744415435311794
   40.75575843549393 
   44.53256351804118 
   50.27030010653037 
   55.57231093619771 
   58.9132798657607  
   ))


(defparameter *calibration-relative-disspiation-length*
  (list
   3.814488129867692
   3.8995725906014087
   3.996822488026833
   4.078867429490328
   4.206479401369255
   4.306739379621605
   4.422156476365529
   4.501088302450566
   4.589083489265088
   4.661836922063561
   4.728478198072885
   4.764691163131293
   4.806979486924731
   4.834055495145416
   4.855078223463543
   4.867044950663009
   4.888030880927016
   4.905999370753274
   4.917988176785211))

(defun interp (x-search x-list y-list)
  (when (< x-search (first x-list))
    (error "X-search is too small"))
  (when (> x-search (first (last x-list)))
    (error "X-search is too big"))
  (let ((np (position x-search x-list :test #'<=))
        (n (position x-search x-list :test #'> :from-end t)))
    (+ (nth n y-list) (* (- (nth np y-list) (nth n y-list)) (/ (- x-search (nth n x-list))
                                                               (- (nth np x-list) (nth n x-list))
                                                               )))))




(defun calculate-ductility-param (E Gf l-c f-t)
  (declare (optimize (debug 3) (speed 0)))
  (let* ((ductility-limit 5)
         (eta 1d0)
         (dissipation-ratio ductility-limit)
         (g-p (/ (* f-t f-t) (* 2 E))))
     (dotimes (i 10)
       (let* (
              (dissipation-length (* dissipation-ratio l-c))
              (local-disp (/ Gf dissipation-length)))
         (let ((new-ductility (/ local-disp g-p)))
           (format t "Iter: ~D - ductility ~F - disp-length ~F - local-disp: ~F~%" i eta dissipation-length local-disp)
           (format t "E_f ~F ~%" (/ (* f-t (+ eta 1)) (* 2 E)))
           (setf eta new-ductility)
           (setf dissipation-ratio (interp eta *calibration-ductility* *calibration-relative-disspiation-length*))
           )
         )))

  )

(defun estimate-gf (eta ft lc &optional (E 1d9))
  (let* ((gf (* eta (/ (expt ft 2) (* 2 E)))))
    (* gf lc)))


(defun test-energy (strain)
  (let* ((de (cl-mpm/constitutive::linear-elastic-matrix 1d0 0.1d0))
         (stress ;(magicl:@ de strain)
           (cl-mpm/constitutive::linear-elastic-mat strain de)
                 ))
    (format t "~A~%" stress)
    (multiple-value-bind (e1 e2 e3) (cl-mpm/damage::principal-stresses-3d strain)
      (multiple-value-bind (s1 s2 s3) (cl-mpm/damage::principal-stresses-3d stress)
        (format t "Energy ~F~%"
                (cl-mpm/fastmaths:dot strain stress)
                ;; (+
                                ;;  (* (max e1 0d0) s1)
                                ;;  (* (max e2 0d0) s2)
                                ;;  (* (max e3 0d0) s3)
                                ;;  )
                )
        ))
    (let ((strain+
            (multiple-value-bind (l v) (cl-mpm/utils::eig
                                        (cl-mpm/utils:voigt-to-matrix strain))
              ;; (loop for i from 0 to 2
              ;;       do
              ;;          (setf (nth i l) (max (nth i l) 0d0)))
              (cl-mpm/utils:matrix-to-voigt (magicl:@ v
                                          (magicl:from-diag l :type 'double-float)
                                          (magicl:transpose v))))))
      ;; (format t "Energy real ~A~%" (magicl:@ de strain+))
      (format t "Energy real ~F~%" (cl-mpm/fastmaths::dot strain+ (magicl:@ de strain+)))
      )
    )
  )

(defun cundall-test ()
  (let* ((load-steps 50)
         (disp-step (/ -0.2d-3 load-steps))
         (target 1d-7)
         (dt-scale 0.70d0)
         (step 0)
         (refine 2.0)
         (work 0d0)
        )

    (defparameter data-cundall-step (list))
    (defparameter data-cundall-load (list))
    (defparameter data-cundall-energy (list))
    (defparameter data-visc-step (list))
    (defparameter data-visc-load (list))
    (defparameter data-visc-energy (list))
    (setup :refine refine :mps 2)
    (setf (cl-mpm:sim-damping-factor *sim*)
          0.7d0)
    (setf *target-displacement* disp-step)
    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (vgplot:close-all-plots)
    (vgplot:figure)
    (time
     (loop for step from 0 to 500
                                        ;(floor (* refine 1000))
           while *run-sim*
           do
              (progn
                (format t "Step ~D~%" step)
                (dotimes (i (floor (* refine 10)))
                  (setf cl-mpm/penalty::*debug-force* 0)
                  (cl-mpm:update-sim *sim*)
                  (incf work (cl-mpm/dynamic-relaxation:estimate-power-norm *sim*))
                  (format t "Step ~D - Load ~E - OOBF ~E - KE ~E~%"
                          step
                          cl-mpm/penalty::*debug-force*
                          (cl-mpm/dynamic-relaxation::estimate-oobf *sim*)
                          (abs (/ (cl-mpm/dynamic-relaxation:estimate-energy-norm *sim*) work)))
                  (setf (cl-mpm:sim-dt *sim*) (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
                  ;; (format t "Target: ~E - Current: ~E Error: ~E~%"
                  ;;         *target-displacement*
                  ;;         (get-disp *terminus-mps*)
                  ;;         (* 100d0 (/ (- *target-displacement* (get-disp *terminus-mps*)) *target-displacement*)))
                  )
                (push (* (cl-mpm:sim-dt *sim*) step) data-visc-step)
                (push cl-mpm/penalty::*debug-force* data-visc-load)
                (vgplot:plot
                 data-visc-step data-visc-load "Visc")
                (swank.live:update-swank))))
    ;; (time
    ;;  (progn
    ;;    (cl-mpm/dynamic-relaxation::converge-quasi-static
    ;;     *sim*
    ;;     :energy-crit 1d-4
    ;;     :oobf-crit 1d-4
    ;;     :dt-scale dt-scale
    ;;     :conv-steps 500
    ;;     :substeps (* 10 1)
    ;;     :post-iter-step
    ;;     (lambda (step ke fnorm)
    ;;       (let ((av (get-disp *terminus-mps*))
    ;;             (load cl-mpm/penalty::*debug-force*)
    ;;             (energy (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
    ;;             )
    ;;         (format t "Conv disp - Target: ~E - Current: ~E - Error: ~E - LOAD ~E~%"
    ;;                 *target-displacement*
    ;;                 av
    ;;                 (abs (- *target-displacement* av ))
    ;;                 load)
    ;;         (push step data-cundall-step)
    ;;         (incf step)
    ;;         (push load data-cundall-load)
    ;;         (push energy data-cundall-energy)
    ;;         (vgplot:plot
    ;;          data-cundall-step data-cundall-load "Cundall")
    ;;         )))))

    ;; (setf step 0)
    ;; (setup)
    (incf *target-displacement* disp-step)
    ;; (change-class *sim* 'cl-mpm/damage::mpm-sim-usl-damage)
    )
  ;; (plot-cundall-test)
  )
(defun plot-cundall-test ()
  (vgplot:close-all-plots)
  (vgplot:figure)
  (vgplot:plot
   data-cundall-step data-cundall-load "Cundall"
   ;; data-visc-step data-visc-load "Visc"
   )
  ;; (vgplot:figure)
  ;; (vgplot:plot
  ;;  data-cundall-step data-cundall-energy "Cundall"
  ;;  ;; data-visc-step data-visc-energy "Visc"
  ;;  )
  )


(defun delete-nodes ()
  (setup)
  (defparameter *test-node* (sb-ext::make-weak-pointer (aref (cl-mpm/mesh:mesh-nodes (cl-mpm:sim-mesh *sim*)) 1 0 0)))
  (sb-ext:search-roots *test-node* :print :verbose)
  (room)
  (let* ((sim *sim*)
         (mesh (cl-mpm:sim-mesh *sim*)))
    (let ((bcs (cl-mpm:sim-bcs sim))
          (prune-count 0))
      (dotimes (i (array-total-size bcs))
        (let ((bc (row-major-aref bcs i)))
          (let ((index (cl-mpm/bc:bc-index bc)))
            ;;nil indexes indiciate global bcs
            (when index
              (let ((node (cl-mpm/mesh:get-node mesh index)))
                (when node
                  (when t
                    (incf prune-count)
                    (setf (row-major-aref bcs i) nil))))))))
      )
    ;;Trim out all nodes that we can get rid of
    (let ((nodes (cl-mpm/mesh:mesh-nodes mesh))
          (prune-count 0))
      (dotimes (i (array-total-size nodes))
        (let ((node (row-major-aref nodes i)))
          (when node
            (when t
              (incf prune-count)
              (setf (row-major-aref nodes i) nil)))))
      )
    ;;Cells
    (let ((cells (cl-mpm/mesh::mesh-cells mesh))
          (prune-count 0)
          )
      (dotimes (i (array-total-size cells))
        (let ((cell (row-major-aref cells i)))
          (when cell
            (when t
              (incf prune-count)
              (setf (row-major-aref cells i) nil)))))))
  (sb-ext:gc :full t)
  (room)
  (sb-ext:search-roots *test-node* :print :verbose)
  )
