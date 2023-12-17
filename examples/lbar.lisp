(defpackage :cl-mpm/examples/lbar
  (:use :cl))
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
;; (setf *block-compile-default* t)
(in-package :cl-mpm/examples/lbar)

(ql:quickload :magicl)
;; (pushnew :cl-mpm-pic *features*)
;; (delete :cl-mpm-pic *features*)
;; (asdf:compile-system :cl-mpm :force T)
(declaim (optimize (debug 3) (safety 3) (speed 0)))

(declaim (notinline plot))
(defun plot (sim)
  (cl-mpm/plotter:simple-plot-3d
   *sim*
   :plot :point
   ;; :colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :zz))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-index mp))
   :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-damage-ybar mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-strain-plastic-vm mp))
   ;; :colour-func (lambda (mp) (magicl:tref (cl-mpm/particle::mp-displacement mp) 1 0))
   )
  ;; (vgplot:format-plot t "replot (~f*x + ~f)~%" 0d0 (+ 0.2d0 *target-displacement*))
  ;; (plot-load-disp)
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
          sum (magicl:norm (cl-mpm/particle:mp-velocity mp))) (length (cl-mpm:sim-mps *sim*))))

(defun get-disp (load-mps)
  ;; (* *t* *tip-velocity*)
  (- (/ (loop for mp in load-mps
              sum (-
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
                           (cl-mpm/fastmath::mag (cl-mpm/mesh::node-force node))
                           )
                        ))))
      (/ force (length load-mps))
      )))

(defun get-reaction-force (load-nodes)
  ;; (cl-mpm/fastmath::mag (cl-mpm/mesh::node-force (nth 0 load-nodes)))
  (loop for mp in load-nodes
        sum
        ;; (cl-mpm/fastmath::mag (cl-mpm/mesh::node-force mp))
        (- (magicl:tref (cl-mpm/mesh::node-force mp) 1 0))
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
                     ))
      )))

;(defparameter *tip-velocity* -0.02d-3)
(defparameter *tip-velocity* -0.000d-3)

(declaim (notinline setup-test-column))
(defun setup-test-column (size block-size offset &optional (e-scale 1) (mp-scale 1))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               ;; 'cl-mpm::mpm-sim-usf
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
      (let* ((impactor-size (list 10d-3
                                  (* 0.99 h-x)))
            (impactors
              (cl-mpm/setup::make-block-mps-list
               (mapcar #'+ offset
                       (list
                        ;; 0d0
                        (- (* 0.5d0 (first block-size)) (* 0.5d0 (first impactor-size)))
                        (+ (second block-size) (* h-x 0.5))
                        ))
               impactor-size
               (mapcar (lambda (e)
                         (round  (* e mp-scale) h-x)) impactor-size)
               density
               'cl-mpm/particle::particle-elastic
               :E 20d9
               :nu 0.20d0
               :gravity -0.0d0
               :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
               :fixed-velocity (list 0 *tip-velocity*)
               :index 1
               )
              ))
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-mps-from-list
               (append
                (cl-mpm/setup::make-block-mps-list
                 offset
                 block-size
                 (mapcar (lambda (e)
                           (round  (* e mp-scale) h-x)
                           ) block-size)
                 density
                 'cl-mpm/particle::particle-limestone
                 :E 25.85d9
                 :nu 0.18d0
                 ;; :elastic-approxmation :plane-stress
                 :fracture-energy 95d0
                 :initiation-stress (* 2.7d6 1d0)
                 :critical-damage 1.000d0
                 :internal-length 25d-3
                 :local-length 25d-3
                 :local-length-damaged 25d-3
                 :compression-ratio 10d0
                 :gravity -0.0d0
                 :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
                 )
                ;; impactors
                )
               )))
      (setf (cl-mpm:sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) t)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (setf (cl-mpm::sim-mass-filter sim) 0d0)
      (let ((ms 1d6))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        (setf (cl-mpm:sim-damping-factor sim)
              (* 1d-3 density ms)
              ;; 1d0
              ))

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

      (let ((dt-scale 1d0))
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

      (let* ((crack-pos
               (loop for mp across (cl-mpm:sim-mps sim)
                     when
                     (not (= (cl-mpm/particle::mp-index mp) 1))
                     maximize (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
             (above-crack
               (loop for mp across (cl-mpm:sim-mps sim)
                     when
                     (and
                      (not (= (cl-mpm/particle::mp-index mp) 1))
                      (= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) crack-pos))
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
                  collect mp)))

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
                (list
                 (cl-mpm/bc::make-bc-fixed left-node-pos
                                           '(0 0 nil))

                 (cl-mpm/bc::make-bc-fixed right-node-pos
                                           '(nil 0 nil)))
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
               ;; *floor-bc*
               (cl-mpm/bc::make-bc-closure
                '(0 0 0)
                (lambda ()
                  (with-accessors ((mesh cl-mpm:sim-mesh)
                                   (dt cl-mpm::sim-dt))
                      sim
                    (let ((datum (* 1d0 (+ *initial-surface* *target-displacement*)))
                          (normal (cl-mpm/utils:vector-from-list  '(0d0 1d0 0d0))))
                      (cl-mpm/penalty::apply-displacement-control-mps mesh (coerce *terminus-mps* 'vector)
                                                       dt
                                                       normal
                                                       datum
                                                       (* density 1d5)
                                                       0d0)
                      )
                    )))
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
(let ((refine (uiop:getenv "REFINE")))
  (when refine
    (setf *refine* (parse-integer (uiop:getenv "REFINE")))
    ))

(declaim (notinline setup))
(defun setup (&key (undercut 0d0))
  ;; (let ((mps-per-dim 4))
  ;;   (defparameter *sim* (setup-test-column '(16 16) '(8 8)  '(0 0) *refine* mps-per-dim)))
  ;; (defparameter *sim* (setup-test-column '(1 1 1) '(1 1 1) 1 1))

  (let* ((mesh-size (/ 0.025 1.0d0))
         (mps-per-cell 2)
         (shelf-height 0.500d0)
         (shelf-length 0.500d0)
         ;; (shelf-length 0.225d0)
         (domain-length (+ shelf-length (* 8 mesh-size)))
         (offset (list
                  (* 2 mesh-size)
                  0d0
                  0d0))


         )
    (defparameter *sim*
      (setup-test-column (list domain-length
                               domain-length
                               (* 2 mesh-size))
                         (list shelf-length shelf-height
                               mesh-size)
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


    )
  (defparameter *target-displacement* 0d0)
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (defparameter *run-sim* t)
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

  (with-open-file (stream (merge-pathnames "output/disp.csv") :direction :output :if-exists :supersede)
    (format stream "disp,load,load-mps~%"))

  (let* ((target-time 0.5d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 1d0)
         (disp-step 0.008d-3)
         )

    (setf cl-mpm/penalty::*debug-force* 0d0)
    (setf cl-mpm/penalty::*debug-force-count* 0d0)
    (time (cl-mpm::update-sim *sim*))
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                    (format t "CFL dt estimate: ~f~%" dt-e)
                    (format t "CFL step count estimate: ~D~%" substeps-e)
                    (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
    ;; (incf *target-displacement* -0.000d-3)
    ;; (incf *target-displacement* disp-step)
    (time (loop for steps from 0 to 100
                while *run-sim*
                do
                   (progn
                     ;; (when (= steps 5)
                     ;;   (setf (cl-mpm::sim-enable-damage *sim*) t)
                     ;;   )
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     ;; (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)


                     (let ((average-force 0d0)
                           (average-disp 0d0)
                           (average-reaction 0d0))
                       (time
                        (dotimes (i substeps);)
                          (unless *data-averaged*
                            (push
                             ;; (- *t*)
                             (get-disp *terminus-mps*)
                             *data-displacement*)

                            (push
                             (get-reaction-force *fixed-nodes*)
                             *data-load*)

                            (push
                             ;; (get-force-mps *sim* *terminus-mps*)
                             (/ cl-mpm/penalty::*debug-force*
                                1d0
                                ;; (max 1 cl-mpm/penalty::*debug-force-count*)
                                )
                             *data-mp-load*)
                            (push
                             *t*
                             *data-time*))
                          (push
                           (get-reaction-force *fixed-nodes*)
                           *data-full-load*)
                          (push
                           *t*
                           *data-full-time*)

                          ;; (incf average-force (/ (get-force-mps *sim* *terminus-mps*) substeps))
                          ;; (incf average-force (/ cl-mpm/penalty::*debug-force* substeps))
                          (incf average-force (/
                                               (/ cl-mpm/penalty::*debug-force*
                                                  1d0
                                                  ;; (max 1 cl-mpm/penalty::*debug-force-count*)
                                                  )
                                               substeps
                                               ))
                          (incf average-reaction
                                (/ (get-reaction-force *fixed-nodes*) substeps)
                                )
                          (incf average-disp
                                (/ (get-disp *terminus-mps*) substeps)
                                )
                          (incf average-reaction (/ (get-reaction-force *fixed-nodes*) substeps))

                          (setf (cl-mpm::sim-enable-damage *sim*) nil)
                          (setf cl-mpm/damage::*delocal-counter-max* 0)
                          (when (= i (- substeps 1))
                            (setf (cl-mpm::sim-enable-damage *sim*) t))

                          (setf cl-mpm/penalty::*debug-force* 0d0)
                          (setf cl-mpm/penalty::*debug-force-count* 0d0)
                          (cl-mpm::update-sim *sim*)
                          ;; (incf *target-displacement* (* dt *tip-velocity*))
                          ;; (push
                          ;;  cl-mpm/penalty::*debug-force*
                          ;;  *data-load*)
                          ;; (incf *target-displacement* (/ -0.01d-3 substeps))
                          ;; (incf *target-displacement* (/ -0.001d-3 substeps))
                          (incf *target-displacement* (/ disp-step substeps))
                          (setf *t* (+ *t* (cl-mpm::sim-dt *sim*))))
                        )
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


                       (format t "Target load: ~f~%" (* *target-displacement* 20d9))
                       (format t "Current load: ~f~%" (* (get-disp *terminus-mps*) 1d9))
                       (format t "Node Load: ~f~%" (get-reaction-force *fixed-nodes*))
                       (format t "Pen load: ~f~%" average-force)
                       
                       (with-open-file (stream (merge-pathnames "output/disp.csv") :direction :output :if-exists :append)
                         (format stream "~f,~f~%"
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
  (vgplot:figure)
  (plot-load-disp)
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*))

(defun plot-time-disp ()
  ;; (vgplot:figure)
  (vgplot:xlabel "Displacement (mm)")
  (vgplot:ylabel "Load (N)")
  (vgplot:plot
   (mapcar (lambda (x) (* x -1d0)) *data-full-time*) (mapcar (lambda (x) (* x 0.1)) *data-full-load*) "nodes"
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
	           (uiop:read-file-string #P"example_data/lbar/load-disp.csv"))))
    (vgplot:xlabel "Displacement (mm)")
    (vgplot:ylabel "Load (N)")
    (let* ((x-model (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
           (x-experiment 0.1d0)
           (x-scale (/ x-experiment x-model)))
      (vgplot:plot
       (lisp-stat:column df 'disp) (lisp-stat:column df 'load) "Data"
       ;; (mapcar (lambda (x) (* x -1d3)) *data-displacement*) *data-node-load* "node"
       (mapcar (lambda (x) (* x 1d3)) *data-displacement*) (mapcar (lambda (x) (* x x-scale)) *data-mp-load*) "mpm-mps"
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
