(defpackage :cl-mpm/examples/lbar
  (:use :cl))
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)
(in-package :cl-mpm/examples/lbar)

;; (ql:quickload :magicl)
;; (pushnew :cl-mpm-pic *features*)
;; (delete :cl-mpm-pic *features*)
;; (asdf:compile-system :cl-mpm :force T)
(declaim (optimize (debug 0) (safety 0) (speed 3)))



(defmethod cl-mpm::update-node-forces ((sim cl-mpm::mpm-sim))
  (with-accessors ((damping cl-mpm::sim-damping-factor)
                   (mass-scale cl-mpm::sim-mass-scale)
                   (mesh cl-mpm::sim-mesh)
                   (dt cl-mpm::sim-dt))
      sim
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (node)
       (cl-mpm::calculate-forces node damping dt mass-scale)
       ;; (cl-mpm::calculate-forces-cundall-conservative node damping dt mass-scale)
       ;; (cl-mpm::calculate-forces-cundall node damping dt mass-scale)
       ))))

(defun cl-mpm/damage::length-localisation (local-length local-length-damaged damage)
  ;; (+ (* local-length (- 1d0 damage)) (* local-length-damaged damage))
  ;; (* local-length (max (sqrt (- 1d0 damage)) 1d-10))
  local-length
  )



(defun setup-visc-damping ()
  (defmethod cl-mpm::update-node-forces ((sim cl-mpm::mpm-sim))
    (with-accessors ((damping cl-mpm::sim-damping-factor)
                     (mass-scale cl-mpm::sim-mass-scale)
                     (mesh cl-mpm::sim-mesh)
                     (dt cl-mpm::sim-dt))
        sim
      (cl-mpm::iterate-over-nodes
       mesh
       (lambda (node)
         (cl-mpm::calculate-forces node damping dt mass-scale))))))

(defun plot-3d (sim)
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :point
   ;; :colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :zz))
   :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-damage-ybar mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-strain-plastic-vm mp))
   ;; :colour-func (lambda (mp) (magicl:tref (cl-mpm/particle::mp-displacement mp) 1 0))
   ))

(declaim (notinline plot))
(defun plot (sim)
  ;; (cl-mpm/plotter:simple-plot-3d
  ;;  *sim*
  ;;  :plot :point
  ;;  ;; :colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :zz))
  ;;  ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-index mp))
  ;;  :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp))
  ;;  ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-damage-ybar mp))
  ;;  ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-strain-plastic-vm mp))
  ;;  ;; :colour-func (lambda (mp) (magicl:tref (cl-mpm/particle::mp-displacement mp) 1 0))
  ;;  )
  ;; (vgplot:format-plot t "replot (~f*x + ~f)~%" 0d0 (+ 0.2d0 *target-displacement*))
  (plot-load-disp)
  ;; (plot-time-disp)
  )


(defun rectangle-sdf (position size)
  (lambda (pos)
    (let* ((pos (magicl:from-list (list
                                        (magicl:tref pos 0 0)
                                        (magicl:tref pos 1 0)) '(2 1) :type 'double-float))
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
          sum (magicl:norm (cl-mpm/particle:mp-velocity mp)))
     (length (cl-mpm:sim-mps *sim*))))

(defun get-disp (load-mps)
  ;; (* *t* *tip-velocity*)
  (-
   (/ (loop for mp in load-mps
            sum (-
                 (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)
                 (* 0.5d0 (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0))
                 )) (length load-mps))
     *initial-surface*))

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
  ;; (loop for mp in load-nodes
  ;;       sum
  ;;       ;; (cl-mpm/fastmaths::mag (cl-mpm/mesh::node-force mp))
  ;;       (- (magicl:tref (cl-mpm/mesh::node-force mp) 1 0))
  ;;       )
  0d0
  )

(defparameter *target-displacement* 0d0)

;(defparameter *tip-velocity* -0.02d-3)
(defparameter *tip-velocity* -0.000d-3)

(defclass mpm-sim-quasi-static-damage (cl-mpm::mpm-sim-quasi-static cl-mpm/damage::mpm-sim-damage) ())

(declaim (notinline setup-test-column))
(defun setup-test-column (size block-size offset &optional (e-scale 1) (mp-scale 1))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               ;:sim-type 'cl-mpm/damage::mpm-sim-usl-damage
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
      (let* ((scaler 1d0))
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-mps-from-list
               (append
                (cl-mpm/setup::make-block-mps-list
                 offset
                 block-size
                 (mapcar (lambda (e mp-s)
                           (round  (* e mp-s) h-x)
                           ) block-size
                             (list mp-scale
                                   mp-scale
                                   mp-scale
                                   ;; 2
                                   ))
                 density
                 'cl-mpm/particle::particle-limestone-delayed
                 :E 25.85d9
                 ;; :E 21.00d9
                 :nu 0.18d0

                 :fracture-energy 95d0
                 :initiation-stress  (* 2.7d6 1d0)
                 :critical-damage 1.0d0

                 :local-length (* 25d-3 scaler)
                 :local-length-damaged (* 25d-3 scaler)


                 ;; :ductility 10d0;6.8d0
                 ;; :ductility 6.0d0
                 ;; :ductility 12d0
                 :ductility 26.85d0
                 ;; :ductility 1.5d0
                 :compression-ratio 10d0
                 :gravity 0.0d0
                 :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 0d0 0d0))
                 :delay-time 0.01d0
                 :phi (* 35d0 (/ pi 180))
                 :psi 0d0
                 :c 3000d3
                 :enable-plasticity nil

                 )
                )
               )))
      ;; (cl-mpm/examples/tpb::calculate-ductility-param 21d9 95d0 (/ 25d-3 (sqrt 1)) 2.7d6)
      (setf (cl-mpm:sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) t)
      (setf (cl-mpm::sim-enable-fbar sim) nil)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (setf (cl-mpm::sim-mass-filter sim) 1d-15)

      (let ((ms 1d0))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        (setf (cl-mpm:sim-damping-factor sim) 0.1d0))

      (dotimes (i 0)
        (dolist (dir (list :x :y))
          (cl-mpm::split-mps-criteria
           sim
           (lambda (mp h)
             (when
                 (and
                  (> (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                     (+ (first offset) (* (first block-size) 0.45))
                     )
                  (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                     (+ (first offset) (* (first block-size) 0.55))
                     )
                  )
               dir
               )))))

      (let ((dt-scale 0.5d0))
        (setf
         (cl-mpm:sim-dt sim)
         (* dt-scale h
            (sqrt (cl-mpm::sim-mass-scale sim))
            (sqrt (/ density (cl-mpm/particle::mp-p-modulus (aref (cl-mpm:sim-mps sim) 0)))))))

      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))

      (let ((cut-size 0.25d0))
        (cl-mpm/setup::remove-sdf
         sim
         (rectangle-sdf
          (list
           (+ (first offset) (second block-size))
           (+ (second offset) 0d0))
          (list
           cut-size
           cut-size
           ))))

      ;;Right most mp
      (let* ((crack-width 0d0
                          ;(/ 0.05d0 4)
                          )
             (crack-pos
               (loop for mp across (cl-mpm:sim-mps sim)
                     when
                     (not (= (cl-mpm/particle::mp-index mp) 1))
                     maximize (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
             (above-crack
               (loop for mp across (cl-mpm:sim-mps sim)
                     when
                     (and
                      (not (= (cl-mpm/particle::mp-index mp) 1))
                      (>= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                         (- crack-pos crack-width)
                         ))
                     collect mp)
               )
             (min-pos (loop for mp in above-crack
                            minimize (magicl:tref(cl-mpm/particle:mp-position mp) 1 0)))
             )
        (defparameter *terminus-mps*
          (loop for mp in above-crack
                when (= min-pos (magicl:tref
                                 (cl-mpm/particle:mp-position mp)
                                 1 0))
                  collect mp))
        (format t "Loading mps ~D~%" (length *terminus-mps*))
        )
      ;;Correct mp
      ;; (let ((dist (loop for mp across (cl-mpm:sim-mps sim)
      ;;                   minimize (cl-mpm/fastmaths::mag-squared
      ;;                             (magicl:.*
      ;;                              (magicl:.- (cl-mpm/particle:mp-position mp)
      ;;                                         (cl-mpm/utils:vector-from-list (list 0.57d0 0.25d0 0d0))
      ;;                                         )
      ;;                              (cl-mpm/utils:vector-from-list (list 1d0 1d0 0d0)))))))
      ;;   (defparameter *terminus-mps*
      ;;     (loop for mp across (cl-mpm:sim-mps sim)
      ;;           when
      ;;           (<=
      ;;            (cl-mpm/fastmaths::mag-squared
      ;;             (magicl:.*
      ;;              (magicl:.- (cl-mpm/particle:mp-position mp)
      ;;                         (cl-mpm/utils:vector-from-list (list 0.57d0 0.25d0 0d0))
      ;;                         )
      ;;              (cl-mpm/utils:vector-from-list (list 1d0 1d0 0d0))))
      ;;            (+ dist 1d-10))
      ;;             collect mp))
      ;;   )

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
                                                  (list left-node-pos right-node-pos)
                                                  ))
        (format t "Fixed node ~A ~%" left-node-pos)
        (format t "Roller node ~A ~%" right-node-pos)
        (setf (cl-mpm:sim-bcs sim)
              (cl-mpm/bc::make-bcs-from-list
               (append
                (cl-mpm/bc::make-outside-bc-var-list
                 (cl-mpm:sim-mesh sim)
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil 0 nil)))
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 0 nil)))
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
                 )
                ;; (list
                ;;  (cl-mpm/bc::make-bc-fixed left-node-pos
                ;;                            '(0 0 nil))

                ;;  (cl-mpm/bc::make-bc-fixed right-node-pos
                ;;                            '(nil 0 nil)))
                ))))
      (defparameter *initial-surface*
        (loop for mp in *terminus-mps*
              when (= 1 (cl-mpm/particle::mp-index mp))
              minimizing (magicl:tref
                          (magicl:.- (cl-mpm/particle:mp-position mp)
                                   (magicl:scale (cl-mpm/particle::mp-domain-size mp) 0.5d0)
                                   )
                        1 0)))

      (format t "~A~%" h-x)
      (setf (cl-mpm::sim-bcs-force-list sim)
            (list
             (cl-mpm/bc:make-bcs-from-list
              (list
               (cl-mpm/bc::make-bc-closure
                nil
                (lambda ()
                  (with-accessors ((mesh cl-mpm:sim-mesh)
                                   (dt cl-mpm::sim-dt))
                      sim
                    (let ((datum (* 1d0 (+ *initial-surface* *target-displacement*)))
                          (normal (cl-mpm/utils:vector-from-list  '(0d0 1d0 0d0)))
                          (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim))))
                      (cl-mpm/penalty::apply-displacement-control-mps
                       mesh
                       (coerce *terminus-mps* 'vector)
                       dt
                       normal
                       datum
                       (* 25.85d9
                          1d2)
                       0d0)))))
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
(defparameter *refine* (/ 1d0 2d0))
;; (let ((refine (uiop:getenv "REFINE")))
;;   (when refine
;;     (setf *refine* (parse-integer (uiop:getenv "REFINE")))
;;     ))

(defun stop ()
  (setf *run-sim* nil)
  (setf cl-mpm/dynamic-relaxation::*run-convergance* nil))

(declaim (notinline setup))
(defun setup (&key (undercut 0d0) (refine 1d0) (mps 2))
  ;; (let ((mps-per-dim 4))
  ;;   (defparameter *sim* (setup-test-column '(16 16) '(8 8)  '(0 0) *refine* mps-per-dim)))
  ;; (defparameter *sim* (setup-test-column '(1 1 1) '(1 1 1) 1 1))
  (setf cl-mpm/dynamic-relaxation::*work* 0d0)

  (let* ((mesh-size (/ 0.025 (* refine 0.5d0)))
         (mps-per-cell mps)
         (shelf-height 0.500d0)
         (shelf-length 0.500d0)
         (domain-length (+ shelf-length
                           0.25))
         (offset (list
                  0.10
                  0d0
                  ;; 0d0
                  ))


         )
    (defparameter *sim*
      (setup-test-column (list domain-length
                               domain-length
                               ;; (* 3 mesh-size)
                               )
                         (list shelf-length
                               shelf-height
                               ;; mesh-size
                               )
                         offset
                         (/ 1d0 mesh-size) mps-per-cell))
    ;; (let ((cut-size 0.25d0))
    ;;   (cl-mpm/setup::remove-sdf
    ;;    *sim*
    ;;    (rectangle-sdf
    ;;     (list
    ;;      (+ (first offset) shelf-height)
    ;;      (+ (second offset) 0d0))
    ;;     (list
    ;;      cut-size
    ;;      cut-size
    ;;      ))))
    (format t "Total weight ~F~%"
            (loop for mp across (cl-mpm:sim-mps *sim*)
                  sum (* 9.8d0 (cl-mpm/particle:mp-mass mp))))

    (defparameter *current-load* 0d0)


    )
  (defparameter *target-displacement* 0d0)
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
  ;; (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))


(defparameter *data-force* '())
(defparameter *data-displacement* '(0d0))
(defparameter *data-load* '(0d0))
(defparameter *data-time* '(0d0))
(defparameter *data-node-load* '(0d0))
(defparameter *target-displacement* 0d0)
(defparameter *data-averaged* t)

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

  (setf *run-sim* t)
  (with-open-file (stream (merge-pathnames "output/disp.csv") :direction :output :if-exists :supersede)
    (format stream "disp,load~%"))
  ;; (let ((ms 1d4))
  ;;   (setf (cl-mpm::sim-mass-scale *sim*) ms)
  ;;   ;; (setf (cl-mpm::sim-damping-factor *sim*) (* ms 5d0))
  ;;   )
  (let ((mass-scale 1d4)
        (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
        (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh *sim*)))
        )
    (setf (cl-mpm::sim-mass-scale *sim*) mass-scale)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (*
           ;; 1d4
           0.1d0
           (sqrt mass-scale)
           (cl-mpm/setup::estimate-critical-damping *sim*)
           )))

  ;; (loop for mp across (cl-mpm:sim-mps *sim*)
  ;;       do (change-class mp 'cl-mpm/particle::particle-limestone))


  (let* ((target-time 0.10d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 0.5d0)
         (load-steps 100)
         (disp-total 0.8d-3)
         (disp-step (/ disp-total load-steps)))

    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup:estimate-elastic-dt *sim* :dt-scale dt-scale))
    (setf (cl-mpm/damage::sim-damage-delocal-counter-max *sim*) 100)
    (setf (cl-mpm::sim-enable-damage *sim*) nil)
    (cl-mpm::update-sim *sim*)
    (setf cl-mpm/penalty::*debug-force* 0d0)
    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (when (typep mp 'cl-mpm/particle::particle-damage)
               (when (= (cl-mpm/particle::mp-index mp) 0)
                 (setf (cl-mpm/particle::mp-delay-time mp) (* (/ target-time load-steps) 1d-1)))))

    ;; (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
    ;;                 (format t "CFL dt estimate: ~f~%" dt-e)
    ;;                 (format t "CFL step count estimate: ~D~%" substeps-e)
    ;;                 (setf substeps substeps-e))
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
      (format t "CFL dt estimate: ~f~%" dt-e)
      (format t "CFL step count estimate: ~D~%" substeps-e)
      (setf substeps substeps-e))


    (setf (cl-mpm::sim-enable-damage *sim*) t)
    (defparameter *data-energy* (list))
    (format t "Substeps ~D~%" substeps)
    (time (loop for steps from 0 to load-steps
                while *run-sim*
                do
                   (progn
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)


                     (let ((average-force 0d0)
                           (average-disp 0d0)
                           (average-reaction 0d0))
                       (time
                        (dotimes (i substeps) ;)

                          ;; (push
                          ;;  average-disp
                          ;;  *data-displacement*)
                          ;; (push
                          ;;  cl-mpm/penalty::*debug-force*
                          ;;  *data-load*)
                          (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm::sim-mesh *sim*)))
                                 (mesh-size (if (= nd 2)
                                                1d0
                                                (cl-mpm/mesh:mesh-resolution (cl-mpm::sim-mesh *sim*)))))
                            (push
                             (get-disp *terminus-mps*)
                             *data-displacement*)
                            (push
                             (/ cl-mpm/penalty::*debug-force* mesh-size)
                             *data-load*))
                          ;; (push
                          ;;  (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*)
                          ;;  *data-energy*)

                          (incf average-force (/ cl-mpm/penalty::*debug-force*
                                                 substeps))
                          (incf average-disp
                                (/ (get-disp *terminus-mps*) substeps)
                                )
                          (incf average-reaction (/ (get-reaction-force *fixed-nodes*) substeps))

                          (setf cl-mpm/penalty::*debug-force* 0d0)
                          (cl-mpm::update-sim *sim*)

                          (incf *target-displacement* (/ disp-step substeps))
                          (setf *t* (+ *t* (cl-mpm::sim-dt *sim*))))
                        )
                       ;;For quasi-
                       ;; (cl-mpm/damage::calculate-damage *sim*)

                       ;;Not averaged
                       (setf average-force (/ cl-mpm/penalty::*debug-force* 1d0))
                       (setf average-reaction (get-reaction-force *fixed-nodes*))
                       (setf average-disp (get-disp *terminus-mps*))

                       ;; (let* ((nd (cl-mpm/mesh:mesh-nd (cl-mpm::sim-mesh *sim*)))
                       ;;        (mesh-size (if (= nd 2)
                       ;;                       1d0
                       ;;                       (cl-mpm/mesh:mesh-resolution (cl-mpm::sim-mesh *sim*)))))
                       ;;   (push
                       ;;    average-disp
                       ;;    *data-displacement*)
                       ;;   (push
                       ;;    (/ average-force mesh-size)
                       ;;    *data-load*))
                       (push
                        (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*)
                        *data-energy*)
                       ;; (push
                       ;;  average-reaction
                       ;;  *data-load*)

                       (let ((mesh-size (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
                         (with-open-file (stream (merge-pathnames "output/disp.csv") :direction :output :if-exists :append)
                           (format stream "~f,~f~%"
                                   average-disp
                                   (/ average-force mesh-size)
                                   ))))

                     (format t "Target: ~E - Current: ~E Error: ~E - energy ~E~%"
                             *target-displacement*
                             (get-disp *terminus-mps*)
                             (* 100d0 (/
                                       (- *target-displacement* (get-disp *terminus-mps*))
                                       *target-displacement*))
                             (energy-norm *sim*))

                     (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       (format t "CFL dt estimate: ~f~%" dt-e)
                       (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf substeps substeps-e))

                     (incf *sim-step*)
                     ;; (plot *sim*)
                     (plot-load-disp)
                     (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                        :terminal "png size 1920,1080"
                                        )

                     (swank.live:update-swank)
                     ))))
  (vgplot:figure)
  ;; (plot-load-disp)
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*))
(defun run-static ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)
  (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes.vtk")) *sim*)
  (defparameter *data-force* '())
  (defparameter *data-displacement* '(0d0))
  (defparameter *data-load* '(0d0))
  (defparameter *data-node-load* '(0d0))
  (defparameter *data-mp-load* '(0d0))
  (defparameter *data-time* '(0d0))
  (defparameter *target-displacement* 0d0)
  (defparameter *data-full-time* '(0d0))
  (defparameter *data-full-load* '(0d0))
  (loop for mp across (cl-mpm:sim-mps *sim*)
        do (change-class mp 'cl-mpm/particle::particle-limestone))

  (with-open-file (stream (merge-pathnames "output/disp.csv") :direction :output :if-exists :supersede)
    (format stream "disp,load,reaction~%")
    (format stream "~f,~f~%" 0d0 0d0 0d0)
    )

  (let* ((target-time 0.5d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 0.5d0)
         (load-steps 25)
         (disp-step (/ 0.8d-3 load-steps))
         ;; (disp-step (/ 1.2d-3 load-steps))
         )
    ;; (setf (cl-mpm:sim-damping-factor *sim*) 0.7d0)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 1d-2
             (cl-mpm/setup::estimate-critical-damping *sim*)))

    (setf cl-mpm/penalty::*debug-force* 0d0)
    (setf cl-mpm/penalty::*debug-force-count* 0d0)
    (time (cl-mpm::update-sim *sim*))
    (setf *run-sim* t)

    (format t "Calculating dt~%")
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                    (format t "CFL dt estimate: ~f~%" dt-e)
                    (format t "CFL step count estimate: ~D~%" substeps-e)
                    (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
    (setf *target-displacement* 0d0)
    (time (loop for steps from 0 below load-steps
                while *run-sim*
                do
                   (progn
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     (let ((average-force 0d0)
                           (average-disp 0d0)
                           (average-reaction 0d0))

                       (incf *target-displacement* disp-step)
                       (setf (cl-mpm::sim-enable-damage *sim*) nil)
                       (time
                        (progn
                          ;; (dotimes (i 10)
                          ;;   (dotimes (j 20)
                          ;;     (setf cl-mpm/penalty::*debug-force* 0)
                          ;;     (cl-mpm:update-sim *sim*))
                          ;;   (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                          ;;     (format t "CFL dt estimate: ~E~%" dt-e)))
                          ;; (format t "Estimated KE ~E - OOBF ~E~%" (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*)(cl-mpm/dynamic-relaxation::estimate-oobf *sim*)
                          ;;         )
                          (handler-bind ((cl-mpm/dynamic-relaxation::non-convergence-error
                                           (lambda (c)
                                             (format t "System failed to converge")
                                             (setf *run-sim* nil)
                                             (cl-mpm/output:save-vtk (merge-pathnames "./output/" (format nil "sim_fail.vtk")) *sim*)
                                             (cl-mpm/output::save-vtk-nodes (merge-pathnames "./output/" (format nil "sim_fail_nodes.vtk")) *sim*)
                                             (invoke-restart 'cl-mpm/dynamic-relaxation::continue)
                                             )))
                            (cl-mpm/dynamic-relaxation::converge-quasi-static
                             *sim*
                             :energy-crit 1d-2
                             :oobf-crit 1d-2
                             :dt-scale dt-scale
                             :conv-steps 200
                             :substeps 50
                             :post-iter-step
                             (lambda (i fnorm oobf)
                               (let ((av (get-disp *terminus-mps*)))
                                 (format t "Conv disp - Target: ~E - Current: ~E - Error: ~E~%"
                                         *target-displacement*
                                         av
                                         (abs (- *target-displacement* av )))))))
                          (setf (cl-mpm::sim-enable-damage *sim*) t)
                          (cl-mpm/damage::calculate-damage *sim*)))
                       (incf average-disp (get-disp *terminus-mps*))
                       (incf average-force cl-mpm/penalty::*debug-force*)
                       (push
                        average-disp
                        *data-displacement*)
                       (push
                        average-force
                        *data-load*)

                       (plot-load-disp)
                       (let* ((mesh-size (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
                             (mesh-size 1d0)
                             )
                         (with-open-file (stream (merge-pathnames "output/disp.csv") :direction :output :if-exists :append)
                           (format stream "~f,~f,~f~%"
                                   average-disp
                                   (/ average-force mesh-size)
                                   (/ average-reaction mesh-size)
                                   ))))


                     (format t "Target: ~E - Current: ~E Error: ~E - energy ~E~%"
                             *target-displacement*
                             (get-disp *terminus-mps*)
                             (* 100d0 (/ (- *target-displacement* (get-disp *terminus-mps*)) *target-displacement*))
                             (energy-norm *sim*))
                     (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       (format t "CFL dt estimate: ~E~%" dt-e)
                       (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf substeps substeps-e))
                     (incf *sim-step*)
                     (swank.live:update-swank)
                     ))))
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*))

(defun plot-time-disp ()
  ;; (vgplot:figure)
  (vgplot:xlabel "Displacement (mm)")
  (vgplot:ylabel "Load (N)")
  (vgplot:plot
   (mapcar (lambda (x) (* x -1d0)) *data-full-time*) (mapcar (lambda (x) (* x 0.1)) *data-full-load*) "nodes"
   ))

(declaim (notinline plot-load-disp))
(defun plot-load-disp ()
  (let ((df (lisp-stat:read-csv
	           (uiop:read-file-string #P"example_data/lbar/load-disp.csv"))))
    (vgplot:xlabel "Displacement (mm)")
    (vgplot:ylabel "Load (N)")
    (let* ((x-model (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
           (x-experiment 0.1d0)
           (x-scale (/ x-experiment x-model))
           ;; (e-max (reduce #'max *data-energy*))
           (x-scale 0.1d0)
           )
      (vgplot:plot
       (lisp-stat:column df 'disp) (lisp-stat:column df 'load) "Data"
       ;; (mapcar (lambda (x) (* x -1d3)) *data-displacement*) *data-node-load* "node"
       (mapcar (lambda (x) (* x 1d3)) *data-displacement*) (mapcar (lambda (x) (* x x-scale)) *data-load*) "mpm-mps"
       ;; (mapcar (lambda (x) (* x 1d3)) *data-displacement*) (mapcar (lambda (x) (* 5000 (/ x e-max))) *data-energy*) "Energy"
       ;; (mapcar (lambda (x) (* x -1d3)) *data-displacement*) (mapcar (lambda (x) (* x -2d9)) *data-displacement*) "LE"
       ))

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
    )
  )

;; (setf lparallel:*kernel* (lparallel:make-kernel 2 :name "custom-kernel"))
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

(defun simple-time ()
  (setup)
  (cl-mpm:update-sim *sim*)
  (let* ((mps (cl-mpm:sim-mps *sim*))
         (mesh (cl-mpm:sim-mesh *sim*))
         (dt (cl-mpm:sim-dt *sim*)))
    (time
     (dotimes (i 100)
       (cl-mpm::p2g-force mesh mps))))
  ;; (time
  ;;  (dotimes (i 100)
  ;;    (cl-mpm::p2g-force *sim*)))
    )

(defun test ()
  (setup)
  (sb-profile:profile "CL-MPM")
  (sb-profile:profile "CL-MPM/PARTICLE")
  (sb-profile:profile "CL-MPM/MESH")
  (sb-profile:profile "CL-MPM/SHAPE-FUNCTION")
  (sb-profile:reset)

  ;; (let* ((mps (cl-mpm:sim-mps *sim*))
  ;;        (mesh (cl-mpm:sim-mesh *sim*)
  ;;        (dt (cl-mpm:sim-dt *sim*)))
  ;;   (time
  ;;    (dotimes (i 10)
  ;;      (cl-mpm::update-stress mesh mps dt))))

  (time
   (dotimes (i 10)
         (cl-mpm::update-sim *sim*)))
  (sb-profile:report)
  (sb-profile:unprofile)
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
                               (cl-mpm/damage::brittle-concrete-linear-d stress E Gf length init-stress)
                               )
                             stress)))
        (vgplot:plot stress damage)))))

(defun test-iter ()
  (setup)
  (let* ((mps (cl-mpm:sim-mps *sim*))
         (mesh (cl-mpm:sim-mesh *sim*))
         (dt (cl-mpm:sim-dt *sim*)))
    (time
     (dotimes (i 10)
       (cl-mpm::update-stress mesh mps dt)
       )
     )))

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

(defun est-eta ()
  (let ((E 21d9)
        (ft 2.7d6)
        (k 1.6d0)
        (R (/ 25d-3 (sqrt 1)))
        (Gf 95d0))
    (* (/ (* 2d0 E) ft)
       (+ (/ Gf (* k R ft))
          (/ ft (* 2 E))))))

;; (cl-mpm::iterate-over-nodes-serial
;;  (cl-mpm:sim-mesh *sim*)
;;  (lambda (n)
;;    (when (cl-mpm/mesh:node-active n)
;;      ;; (pprint (cl-mpm/mesh::node-velocity n))
;;      (pprint (cl-mpm/mesh::node-velocity n))
;;      (pprint (cl-mpm/fastmaths::mag (cl-mpm/mesh::node-velocity n)))
;;      ;; (format t "~a~%" )
;;      ))
;;  )



(defun cundall-test ()
  (defparameter *timesteps* (list))
  (defparameter *loads* (list))
  (setf *run-sim* t)
  (loop for refine in (list 1 2 3 4 5)
        do
           (let ((disp-step (/ 0.8d-3 50))
                 (target 1d2)
                 (dt-scale 0.8d0)
                 (step 0)
                 (time 0d0)
                 ;; (refine 1.0d0)
                 )

             (defparameter data-cundall-step (list))
             (defparameter data-cundall-load (list))
             (defparameter data-cundall-energy (list))
             (defparameter data-visc-step (list))
             (defparameter data-visc-load (list))
             (defparameter data-visc-energy (list))
             ;;2 302
             ;;3 286
             ;;4 285
             ;(setup :refine refine :mps 2)
             (setup :refine 1 :mps (+ refine 1))
             (setf (cl-mpm:sim-dt *sim*)
                   (cl-mpm/setup:estimate-elastic-dt *sim* :dt-scale dt-scale))

             (let ((mass-scale 1d4)
                   (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
                   (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh *sim*)))
                   )
               (setf (cl-mpm::sim-mass-scale *sim*) mass-scale)
               (setf (cl-mpm:sim-damping-factor *sim*)
                     (*
                      0d0
                      (sqrt mass-scale)
                      (cl-mpm/setup::estimate-critical-damping *sim*)
                      )))
             (incf *target-displacement* disp-step)
             ;; (change-class *sim* 'cl-mpm/damage::mpm-sim-usl-damage)
             ;; (setf (cl-mpm:sim-damping-factor *sim*) 0.7d0)
             ;; (vgplot:close-all-plots)
             ;; (vgplot:figure)
             (time
              (loop for step from 0 to 20
                                        ;(floor (* refine 1000))
                    while *run-sim*
                    do
                       (progn
                         (dotimes (i
                                   ;; 10
                                   ;; (round (* 1d0 refine))
                                   1
                                   ;; (floor 1d-5 (cl-mpm:sim-dt *sim*))
                                     )
                           (setf cl-mpm/penalty::*debug-force* 0)
                           (cl-mpm:update-sim *sim*)
                           (incf time (cl-mpm:sim-dt *sim*))
                           )

                         ;; (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" step)) *sim*)
                         ;; (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" step)) *sim*)

                         (format t "Step ~D - Load ~E - OOBF ~E - KE ~E - Disp ~E ~%"
                                 step
                                 cl-mpm/penalty::*debug-force*
                                 (cl-mpm/dynamic-relaxation::estimate-oobf *sim*)
                                 (abs (/ (cl-mpm/dynamic-relaxation:estimate-energy-norm *sim*) 1d0))
                                 (- (get-disp *terminus-mps*) *target-displacement*))
                         (push time data-visc-step)
                         (push cl-mpm/penalty::*debug-force* data-visc-load)
                         ;; (vgplot:plot
                         ;;  data-visc-step data-visc-load "Visc")
                         (ecase (length *timesteps*)
                               (0
                                (vgplot:plot
                                 data-visc-step data-visc-load "0"))
                               (1
                                (vgplot:plot
                                 (nth 0 *timesteps*) (nth 0 *loads*) "0"
                                 data-visc-step data-visc-load "1"))
                               (2
                                (vgplot:plot
                                 (nth 0 *timesteps*) (nth 0 *loads*) "0"
                                 (nth 1 *timesteps*) (nth 1 *loads*) "1"
                                 data-visc-step data-visc-load "2"))
                               (3
                                (vgplot:plot
                                 (nth 0 *timesteps*) (nth 0 *loads*) "0"
                                 (nth 1 *timesteps*) (nth 1 *loads*) "1"
                                 (nth 2 *timesteps*) (nth 2 *loads*) "2"
                                 data-visc-step data-visc-load "3"))
                               (4
                                (vgplot:plot
                                 (nth 0 *timesteps*) (nth 0 *loads*) "0"
                                 (nth 1 *timesteps*) (nth 1 *loads*) "1"
                                 (nth 2 *timesteps*) (nth 2 *loads*) "2"
                                 (nth 3 *timesteps*) (nth 3 *loads*) "3"
                                 data-visc-step data-visc-load "4"))
                               (t
                                (vgplot:plot
                                 ;; (nth 0 *timesteps*) (nth 0 *loads*) "0"
                                 ;; (nth 1 *timesteps*) (nth 1 *loads*) "1"
                                 ;; (nth 2 *timesteps*) (nth 2 *loads*) "2"
                                 ;; (nth 3 *timesteps*) (nth 3 *loads*) "3"
                                 ;; (nth 4 *timesteps*) (nth 4 *loads*) "4"
                                 data-visc-step data-visc-load "5"))
                               )
                         (swank.live:update-swank))))
             (push data-visc-step *timesteps*)
             (push data-visc-load *loads*)
             ))
  )
(defun quasi-static-test ()
  (defparameter *timesteps* (list))
  (defparameter *loads* (list))
  (setf *run-sim* t)
  (loop for refine in (list 1)
        do
           (let ((disp-step (/ 0.8d-3 50))
                 (target 1d2)
                 (dt-scale 0.5d0)
                 (step 0)
                 ;; (refine 1.0d0)
                 )

             (defparameter data-cundall-step (list))
             (defparameter data-cundall-load (list))
             (defparameter data-cundall-energy (list))
             (defparameter data-visc-step (list))
             (defparameter data-visc-load (list))
             (defparameter data-visc-energy (list))
             ;;2 302
             ;;3 286
             ;;4 285
                                        ;(setup :refine refine :mps 2)
             (setup :refine refine :mps 2)
             (setf (cl-mpm:sim-dt *sim*)
                   (cl-mpm/setup:estimate-elastic-dt *sim* :dt-scale dt-scale))

             (let ((mass-scale 1d0)
                   (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
               (setf (cl-mpm::sim-mass-scale *sim*) mass-scale)
               (setf (cl-mpm:sim-damping-factor *sim*)
                     0.01d0
                     ))
             (setup-visc-damping)
             (setf (cl-mpm:sim-damping-factor *sim*)
                   (* 1d-2
                      (cl-mpm/setup::estimate-critical-damping *sim*)))

             (incf *target-displacement* disp-step)
             ;; (change-class *sim* 'cl-mpm/damage::mpm-sim-usl-damage)
             ;; (setf (cl-mpm:sim-damping-factor *sim*) 0.7d0)
             (vgplot:close-all-plots)
             (vgplot:figure)
             (format t "Final load ~f~%" cl-mpm/penalty::*debug-force*)
             (progn
               (cl-mpm/dynamic-relaxation::converge-quasi-static
                *sim*
                ;; :energy-crit target
                :energy-crit 1d-2
                :oobf-crit 1d-2
                :dt-scale dt-scale
                :conv-steps 200
                :substeps 10            ;(* 50 refine)
                :post-iter-step
                (lambda (&rest args)
                  (let ((av (get-disp *terminus-mps*))
                        (load cl-mpm/penalty::*debug-force*)
                        (energy (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
                        )
                    (format t "Conv disp - Target: ~E - Current: ~E - Error: ~E - LOAD ~E~%"
                            *target-displacement*
                            av
                            (abs (- *target-displacement* av))
                            load
                            )
                    (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" step)) *sim*)
                    (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" step)) *sim*)
                    (push step data-visc-step)
                    (incf step)
                    (push load data-visc-load)
                    (push energy data-visc-energy)
                    (vgplot:plot
                     data-visc-step data-visc-load "Visc")
                    ))
                ))
             ))
  )
(defun plot-cundall-test ()
  (vgplot:close-all-plots)
  (vgplot:figure)
  (vgplot:plot
   ;; data-cundall-step data-cundall-load "Cundall"
   data-visc-step data-visc-load "Visc"
               )
  (vgplot:figure)
  (vgplot:plot
   ;; data-cundall-step data-cundall-energy "Cundall"
   data-visc-step data-visc-energy "Visc"
               )
  )

