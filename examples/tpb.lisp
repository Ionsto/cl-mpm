(defpackage :cl-mpm/examples/tpb
  (:use :cl))
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(in-package :cl-mpm/examples/tpb)

;; (pushnew :cl-mpm-pic *features*)
;; (delete :cl-mpm-pic *features*)
;; (asdf:compile-system :cl-mpm :force T)
(declaim (optimize (debug 3) (safety 3) (speed 0)))

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
  )


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
  (/ (loop for mp in load-mps
           sum (magicl:tref (cl-mpm/particle::mp-displacement mp) 1 0))
     (length load-mps))
  )

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
                        (* 1d0
                           (cl-mpm/fastmath::mag (cl-mpm/mesh::node-force node))
                           )
                        ))))
      force
      )))

(defun get-reaction-force (load-nodes)
  ;; (cl-mpm/fastmath::mag (cl-mpm/mesh::node-force (nth 0 load-nodes)))
  (/ (loop for mp in load-nodes
           sum
           ;; (cl-mpm/fastmath::mag (cl-mpm/mesh::node-force mp))
           (- (magicl:tref (cl-mpm/mesh::node-force mp) 1 0))
           )
     1d0
     ;(length load-nodes)
     )
  )

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

(defparameter *tip-velocity* -0.1d-3)
(defun setup-test-column (size block-size offset &optional (e-scale 1) (mp-scale 1))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               #'cl-mpm/shape-function:make-shape-function-bspline
               ;; 'cl-mpm::mpm-sim-usf
               'cl-mpm/damage::mpm-sim-damage
               ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 1.7d3)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
         )
    (declare (double-float h density))
    (progn
      (let* ((impactor-size (list h-x h-x))
            (impactors
              (cl-mpm/setup::make-block-mps-list
               (mapcar #'+ offset
                       (list
                        (- (* 0.5d0 (first block-size)) (* 0.5d0 (first impactor-size)))
                        (+ (second block-size) h-x)
                        ))
               impactor-size
               (mapcar (lambda (e) (* e e-scale mp-scale)) impactor-size)
               density
               'cl-mpm/particle::particle-elastic
               :E 20d9
               :nu 0.20d0
               :gravity -9.8d0
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
                 ;; (mapcar (lambda (e) (ceiling (* e e-scale mp-scale))) block-size)
                 (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                 density
                 'cl-mpm/particle::particle-concrete
                 :E 1d9
                 :nu 0.20d0
                 :fracture-energy 9d0
                 :initiation-stress (* 2.4d6 1d0)
                 :critical-damage 0.90d0
                 :local-length 5d-2
                 :local-length-damaged 5d-2
                 ;; :local-length-damaged 0.01d0
                 :gravity -9.8d0
                 :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
                 )
                ;; impactors
                )
               )))
      (setf (cl-mpm:sim-allow-mp-split sim) t)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (setf (cl-mpm::sim-mass-filter sim) 0d0)
      (let ((ms 1d5))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        ;; (setf (cl-mpm:si-damping-factor sim) (* 1d-2 density ms))
        ;; (setf (cl-mpm:sim-damping-factor sim) 10.0d0)
        ;; (setf (cl-mpm:sim-damping-factor sim) (* 1d-1 density ms))
        ;; (setf (cl-mpm:sim-damping-factor sim) (* 1d-2 density ms))
        (setf (cl-mpm:sim-damping-factor sim)
              ;; 1d4
              (* 0.01d0 density ms)))

      (dotimes (i 0)
        (dolist (dir (list :x :y))
          (cl-mpm::split-mps-criteria
           sim
           (lambda (mp h)
             (when
                 (and
                  (> (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                     120)
                  (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                     200 
                     )
                  (> (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                     10
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

      (let* ((crack-width h-x)
             (crack-left (- (+ (first offset)  (* 0.5d0 (first block-size))) crack-width))
             (crack-right (+ (+ (first offset) (* 0.5d0 (first block-size)) crack-width)))
             (above-crack
               (loop for mp across (cl-mpm:sim-mps sim)
                     when
                     (and
                      (not (= (cl-mpm/particle::mp-index mp) 1))
                      (> (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) crack-left)
                      (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) crack-right))
                     collect mp)
               )
             (max-pos (loop for mp in above-crack
                            maximize (magicl:tref(cl-mpm/particle:mp-position mp) 1 0)))
             )
        (defparameter *terminus-mps*
          (loop for mp in above-crack
                when (= max-pos (magicl:tref
                                 (cl-mpm/particle:mp-position mp)
                                 1 0))
                  collect mp)))
      ;; (setf *terminus-mps* (mapcar (lambda (i) (nth i *terminus-mps*))
      ;;                              (list (floor (- (length *terminus-mps*) 1) 2))))
      ;; (when (> (length *terminus-mps*) 2)
      ;;   (setf *terminus-mps* (mapcar (lambda (i) (nth i *terminus-mps*))
      ;;                                (list (floor (- (length *terminus-mps*) 1) 2)
      ;;                                      (+ (floor (- (length *terminus-mps*) 1) 2) 1)))))

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
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil 0 nil)))
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
                 )
                (list
                 (cl-mpm/bc::make-bc-fixed left-node-pos
                                           '(0 0 nil))

                 (cl-mpm/bc::make-bc-fixed right-node-pos
                                           '(nil 0 nil)))
                ))))

      (defparameter *floor-bc*
        (cl-mpm/penalty::make-bc-penalty-point-normal
         sim
         (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
         (cl-mpm/utils:vector-from-list (list 00d0 (second offset) 0d0))
         (* density 1d5)
         0.0d0
         ;; 1d1
         ))

      (setf (cl-mpm::sim-bcs-force-list sim)
            (list
             (cl-mpm/bc:make-bcs-from-list
              (list
               ;; *floor-bc*
               (cl-mpm/bc::make-bc-closure
                '(0 0 0)
                (lambda ()
                  (with-accessors ((mesh cl-mpm:sim-mesh)
                                   (dt cl-mpm::sim-dt)
                                   )
                      sim
                    (let ((datum (* -1d0 (+ 0.20d0 *target-displacement*)))
                          (normal (cl-mpm/utils:vector-from-list  '(0d0 -1d0 0d0))))
                      (cl-mpm/penalty::apply-force-mps mesh (coerce *terminus-mps* 'vector )
                                                       dt
                                                       normal
                                                       datum
                                                       (* density 2d5)
                                                       0d0))
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

(defun setup (&key (undercut 0d0))
  ;; (let ((mps-per-dim 4))
  ;;   (defparameter *sim* (setup-test-column '(16 16) '(8 8)  '(0 0) *refine* mps-per-dim)))
  ;; (defparameter *sim* (setup-test-column '(1 1 1) '(1 1 1) 1 1))

  (let* ((mesh-size 0.01)
         (mps-per-cell 2)
         (shelf-height 0.100d0)
         (shelf-length 0.55d0)
         (domain-length (+ shelf-length (* 4 mesh-size)))
         (offset (list (* 2 mesh-size) (* shelf-height 1)))
         )
    (defparameter *sim*
      (setup-test-column (list domain-length (+ mesh-size (* shelf-height 3)))
                         (list shelf-length shelf-height)
                         offset
                         (/ 1d0 mesh-size) mps-per-cell))
    ;; (let ((cut-depth (* 1 50d-3))
    ;;       (cut-width 2.5d-3))
    ;;   (cl-mpm/setup::remove-sdf
    ;;    *sim*
    ;;    (rectangle-sdf
    ;;     (list (+ (first offset) (* 0.5d0 shelf-length))
    ;;           (+ (second offset) 0d0)
    ;;           )
    ;;     (list
    ;;      2.5d-3
    ;;      ;; 5d-3
    ;;      cut-depth
    ;;      )))
    ;;   )
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

(defun run ()
  (vgplot:close-all-plots)
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)
  (defparameter *data-force* '())
  (defparameter *data-displacement* '(0d0))
  (defparameter *data-load* '(0d0))

  (let* ((target-time 1d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 1d0)
         )
    (cl-mpm::update-sim *sim*)
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                    (format t "CFL dt estimate: ~f~%" dt-e)
                    (format t "CFL step count estimate: ~D~%" substeps-e)
                    (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
    (defparameter *target-displacement* 0d0)
    (time (loop for steps from 0 to 500
                while *run-sim*
                do
                   (progn
                     ;; (when (= steps 5)
                     ;;   (setf (cl-mpm::sim-enable-damage *sim*) t)
                     ;;   )
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)

                     (setf cl-mpm/penalty::*debug-force* 0d0)
                     (setf cl-mpm/penalty::*debug-force-count* 0d0)
                     (let ((average-force 0d0))
                       (time
                        (dotimes (i substeps);)
                          ;; (push
                          ;;  ;; *target-displacement*
                          ;;  *t*
                          ;;  ;; (get-disp *terminus-mps*)
                          ;;  *data-displacement*)
                          ;; (push
                          ;;  ;; (get-reaction-force *fixed-nodes*)
                          ;;  (/ cl-mpm/penalty::*debug-force* substeps)
                          ;;  *data-load*)

                          ;; (incf average-force (/ (get-force-mps *sim* *terminus-mps*) substeps))
                          ;; (incf average-force (/ cl-mpm/penalty::*debug-force* substeps))
                          (incf average-force (/
                                               (/ cl-mpm/penalty::*debug-force*
                                                  (max 1 cl-mpm/penalty::*debug-force-count*))
                                               substeps
                                               ))
                          ;; (incf average-force (/ (get-reaction-force *fixed-nodes*) substeps))
                          (setf cl-mpm/penalty::*debug-force* 0d0)
                          (setf cl-mpm/penalty::*debug-force-count* 0d0)
                          (cl-mpm::update-sim *sim*)
                          ;; (push
                          ;;  cl-mpm/penalty::*debug-force*
                          ;;  *data-load*)
                          ;; (incf *target-displacement* (/ -0.01d-3 substeps))
                          ;; (incf *target-displacement* (/ -0.001d-3 substeps))
                          (setf *t* (+ *t* (cl-mpm::sim-dt *sim*))))
                        )
                       (push
                         (get-disp *terminus-mps*)
                        *data-displacement*)
                       (push
                        average-force
                        *data-load*)


                       (format t "Target load: ~f~%" (* *target-displacement* 1d9))
                       (format t "Current load: ~f~%" (* (get-disp *terminus-mps*) 1d9))
                       (format t "Node Load: ~f~%" (get-reaction-force *fixed-nodes*))
                       (format t "Pen load: ~f~%" average-force)
                       )
                     (incf *target-displacement* -0.01d-3)

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
  (plot-load-disp)
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*))

(defun plot-load-disp ()
  ;; (vgplot:figure)
  (vgplot:xlabel "Displacement (mm)")
  (vgplot:ylabel "Load (N)")
  (vgplot:plot (mapcar (lambda (x) (* x -1d3)) *data-displacement*) *data-load*)
  ;; (vgplot:plot (mapcar (lambda (x) (* x 1d0)) *data-displacement*) *data-load*)
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

(defun test ()
  (setup)
  ;; (sb-profile:profile "CL-MPM")
  ;; (sb-profile:profile "CL-MPM/PARTICLE")
  ;; (sb-profile:profile "CL-MPM/MESH")
  ;; (sb-profile:profile "CL-MPM/SHAPE-FUNCTION")
  ;; (sb-profile:reset)
  (time
   (dotimes (i 10)
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