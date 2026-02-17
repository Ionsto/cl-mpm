(defpackage :cl-mpm/examples/tpb
  (:use :cl))
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(setf cl-mpm/settings:*optimise-setting* cl-mpm/settings::*optimise-speed*)
;; (setf *block-compile-default* t)
;; (setf cl-mpm/settings:*optimise-setting* cl-mpm/settings::*optimise-debug*)
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)
(in-package :cl-mpm/examples/tpb)

;; (ql:quickload :magicl)
;; (pushnew :cl-mpm-pic *features*)
;; (delete :cl-mpm-pic *features*)
;; (asdf:compile-system :cl-mpm :force T)
(declaim (optimize (debug 3) (safety 3) (speed 0)))



(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-elastic-damage) dt fbar)
  (cl-mpm::update-stress-kirchoff-dynamic-relaxation mesh mp dt fbar))

(declaim (notinline plot))
(defun plot-domain (sim)
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :deformed
   :colour-func #'cl-mpm/particle::mp-index
   ;:colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp))
   )
  ;; (vgplot:format-plot t "replot (~f*x + ~f)~%" 0d0 (+ 0.2d0 *target-displacement*))
  )
(defun plot (sim)
  (plot-load-disp)
  ;; (plot-time-disp)
  )
(defun stop ()
  (setf *run-sim* nil)
  (setf cl-mpm/dynamic-relaxation::*run-convergance* nil))

(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-limestone) dt)
  (with-accessors (;(stress cl-mpm/particle::mp-undamaged-stress)
                   (y cl-mpm/particle::mp-damage-y-local)
                   (strain cl-mpm/particle::mp-strain)
                   (init-stress cl-mpm/particle::mp-initiation-stress)
                   (ybar cl-mpm/particle::mp-damage-ybar)
                   (de cl-mpm/particle::mp-elastic-matrix)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (E cl-mpm/particle::mp-e)
                   (nu cl-mpm/particle::mp-nu)
                   )
      mp
    (declare (double-float E nu))
    (progn
      (let* ((angle 30d0)
             (stress (cl-mpm/fastmaths:fast-scale!
                      (cl-mpm/constitutive:linear-elastic-mat strain de)
                      (/ 1d0 (magicl:det def)))))
        (setf y
              (cl-mpm/damage::tensile-energy-norm strain E de)
              ;; (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile stress (* angle (/ pi 180d0)))
              ;; (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile stress (* angle (/ pi 180d0)))
              )))))

(defparameter *current-load* 0d0)

(defparameter *target-displacement* 0d0)
(defparameter *tip-velocity* -0.000d-3)
(declaim (notinline setup-test-column))
(defun setup-test-column (size block-size offset &key (e-scale 1) (mp-scale 1)
                                                   (epsilon-scale 1d2)
                                                   )
  (let* ((sim (cl-mpm/setup::make-simple-sim
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               ;; :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
               :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul
               :args-list (list :enable-aggregate t
                                :enable-damage t)))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 2.2d3)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size)))
    (declare (double-float h density))
    (progn
      (let* ()
        (let* (;(crack-scale 7d0)
               (crack-scale 1.0d0)
               ;; (length-scale 5.3d-3)
               ;(length-scale (* 5.4d-3 (sqrt 7)))
               (length-scale (* 5.4d-3 1d0))
               (gf (* 1d0 48d0))
               (E 15.3d9)
               (init-stress 3.45d6)
               (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress E)))

          (format t "Estimated ductility ~F~%" ductility)
          (format t "Actual local length ~F~%" (* crack-scale length-scale))
          (format t "Mesh size ~F~%" (* h-x))
          (format t "Length/Mesh res ~F~%" (/ (* crack-scale length-scale) (* 2d0 h-x)))
          ;; (cl-mpm/damage::estimate-ductility-jirsek2004 (* 48d0 0.5d0) (* length-scale (sqrt 7)) 3.45d6 15.3d9)
          (cl-mpm:add-mps
           sim
           (cl-mpm/setup::make-block-mps
            offset
            block-size
            (mapcar (lambda (e) (round (* e mp-scale) h-x)) block-size)
            density
            ;; 'cl-mpm/particle::particle-concrete
            ;; 'cl-mpm/particle::particle-elastic-damage
            ;; :E 15.3d9
            ;; :nu 0.15d0
            ;; :initiation-stress 3.45d6
            ;; ;;Material parameter
            ;; :local-length length-scale
            ;; :ductility ductility
            'cl-mpm/particle::particle-limestone
            :E E
            :nu 0.15d0
            :initiation-stress init-stress
            ;;Material parameter
            :local-length length-scale
            :compression-ratio 8d0
            :ductility ductility
            ;; :gravity-axis (cl-mpm/utils:vector-zeros)
            ))))
      (setf (cl-mpm::sim-gravity sim) 0d0)
      ;; (calculate-ductility-param 18d9  48d0 5.3d-3 3.4d6)
      ;calculate-ductility-param (E Gf l-c f-t)
      (setf (cl-mpm:sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) t)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (cl-mpm/setup::set-mass-filter sim density :proportion 1d-15)
      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))
      (let* ((crack-width 1d-2)
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
             (max-pos (loop for mp in above-crack maximize (magicl:tref(cl-mpm/particle:mp-position mp) 1 0)))
             (min-x-pos (loop for mp in above-crack minimize (magicl:tref(cl-mpm/particle:mp-position mp) 0 0))))
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

        (let* ((hx h-x)
               (epsilon (* 15.3d9 epsilon-scale))
               (friction 0d0)
               (damping 0d0)
               (smoothness 1)
               (center (cl-mpm/utils:vector-from-list (list (- (first block-size) hx) (float (+ (second offset) 0d0) 0d0) 0d0)))
               (left (cl-mpm/fastmaths::fast-.+ center (cl-mpm/utils:vector-from-list (list (- hx) 0d0 0d0))))
               (left-down (cl-mpm/fastmaths::fast-.+ left (cl-mpm/utils:vector-from-list (list 0d0 (- hx) 0d0))))
               (right (cl-mpm/fastmaths::fast-.+ center (cl-mpm/utils:vector-from-list (list (+ hx) 0d0 0d0))))
               (right-down (cl-mpm/fastmaths::fast-.+ right (cl-mpm/utils:vector-from-list (list 0d0 (- hx) 0d0))))
               ;; (penlist
               ;;   (cl-mpm/penalty::make-bc-line-segments
               ;;    sim
               ;;    (list
               ;;     ;; (cl-mpm/utils:vector-from-list (list (- (first block-size) hx) (float (- (second offset) hx) 0d0) 0d0))
               ;;     ;; (cl-mpm/utils:vector-from-list (list (- (first block-size) hx) (float (+ (second offset) 0d0) 0d0) 0d0))
               ;;     ;; (cl-mpm/utils:vector-from-list (list (+ (first block-size) hx) (float (+ (second offset) 0d0) 0d0) 0d0))
               ;;     ;; (cl-mpm/utils:vector-from-list (list (+ (first block-size) hx) (float (- (second offset) hx) 0d0) 0d0))
               ;;     )
               ;;    epsilon friction damping))
               (left-smooth
                 (cl-mpm/penalty::make-bc-penalty-smooth-corner
                  sim
                  left-down
                  left
                  center
                  smoothness
                  epsilon
                  0d0
                  damping
                  ))
               (right-smooth
                 (cl-mpm/penalty::make-bc-penalty-smooth-corner
                  sim
                  center
                  right
                  right-down
                  smoothness
                  epsilon
                  0d0
                  damping
                  ))
               )
          (defparameter *penalty-point*
            (cl-mpm/penalty::make-bc-penalty-structure
             sim
             epsilon
             friction
             damping
             (list
              left-smooth
              right-smooth))
            
            ;; (cl-mpm/penalty::make-bc-penalty-distance-point
            ;;  sim
            ;;  (cl-mpm/utils:vector-from-list '(0d0 -1d0 0d0))
            ;;  (cl-mpm/utils:vector-from-list (list 0d0
            ;;                                       (+ (second block-size) (second offset))
            ;;                                       0d0))
            ;;  h-x
            ;;  (* 15.3d9 1d2)
            ;;  0d0
            ;;  0d0
            ;;  )
            ))
        (cl-mpm/penalty::bc-increment-center *penalty-point*
                                             (cl-mpm/utils:vector-from-list (list 0d0 1d-6 0d0)))

        (cl-mpm/setup::setup-bcs
         sim
         :left '(0 nil nil))

        (cl-mpm::add-bcs
         sim
         (cl-mpm/bc::make-bc-fixed
          right-node-pos
          '(nil 0 nil)))
        ;; (cl-mpm::add-bcs
        ;;  sim
        ;;  (cl-mpm/bc::make-bc-fixed
        ;;   (mapcar #'+ right-node-pos (list 0 -1 0))
        ;;   '(nil 0 nil)))
        )

      (defparameter *initial-surface*
        (loop for mp across (cl-mpm:sim-mps sim)
              when (not (= 1 (cl-mpm/particle::mp-index mp)))
              maximize (magicl:tref
                        (magicl:.+ (cl-mpm/particle:mp-position mp)
                                   (magicl:scale (cl-mpm/particle::mp-domain-size mp) 0.5d0)
                                   )
                        1 0
                        )))


      (defparameter *displacement* 0d0)
      

      (let* ((hx h-x)
             (ly 1d-2)
             (epsilon (* 15.3d9 1d2))
             (friction 0d0)
             (damping 0d0)
             (penlist
               (cl-mpm/penalty::make-bc-line-segments
                sim
                (list
                 (cl-mpm/utils:vector-from-list (list hx  (float (+ (second block-size) (second offset) ly) 0d0) 0d0))
                 (cl-mpm/utils:vector-from-list (list hx  (float (+ (second block-size) (second offset)) 0d0) 0d0))
                 (cl-mpm/utils:vector-from-list (list 0d0 (float (+ (second block-size) (second offset)) 0d0) 0d0))
                 )
                epsilon friction damping)))
        (defparameter *penalty*
          (cl-mpm/penalty::make-bc-penalty-structure
           sim
           epsilon
           friction
           damping
           penlist
           )
          ;; (cl-mpm/penalty::make-bc-penalty-distance-point
          ;;  sim
          ;;  (cl-mpm/utils:vector-from-list '(0d0 -1d0 0d0))
          ;;  (cl-mpm/utils:vector-from-list (list 0d0
          ;;                                       (+ (second block-size) (second offset))
          ;;                                       0d0))
          ;;  h-x
          ;;  (* 15.3d9 1d2)
          ;;  0d0
          ;;  0d0
          ;;  )
          ))
      (defparameter *displacement* 0d0)
      (defparameter *last-pos* 0d0)
      (defparameter *penalty-controller*
        (cl-mpm/bc::make-bc-closure
         nil
         (lambda ()
           (let ((delta (- *displacement* *last-pos*)))
             (cl-mpm/penalty::bc-increment-center *penalty* (cl-mpm/utils:vector-from-list (list 0d0 delta 0d0)))
             (setf *last-pos* *displacement*)))))
      (defparameter *current-inc* 0d0)
      (format t "~A~%" h-x)
      (cl-mpm:add-bcs-force-list
       sim
       (list
        ;; (cl-mpm/bc:make-bcs-from-list
        ;;  (list
        ;;   *penalty-controller*))
        (cl-mpm/bc:make-bcs-from-list
         (list
          *penalty-point*))
        (cl-mpm/bc:make-bcs-from-list
         (list
          *penalty*))))
      sim)))


(defparameter *sim* nil)
(defparameter *run-sim* t)
(defparameter *t* 0)
(defparameter *sim-step* 0)
(defparameter *refine* (/ 1d0 1d0))

(declaim (notinline setup))
(defun setup (&key
                (refine 1d0)
                (mps 2)
                (epsilon-scale 1d2))
  (let* ((mesh-size (/ 0.0102 (* refine 1.0)))
         (mps-per-cell mps)
         (shelf-height 0.102d0)
         (shelf-length (* shelf-height 2.0))
         ;; (shelf-length 0.225d0)
         (domain-length (+ (* shelf-length 2)
                           ;(* shelf-height 0.5d0)
                           ))
         (offset-y (* shelf-height 0.5))
         (offset (list 0d0 offset-y))
         )
    (defparameter *sim*
      (setup-test-column (list domain-length
                               (* shelf-height 2.0)
                               )
                         (list shelf-length shelf-height)
                         offset
                         :e-scale (/ 1d0 mesh-size)
                         :mp-scale mps-per-cell
                         :epsilon-scale epsilon-scale
                         ))
    (let ((cut-depth (* 0.4d0 shelf-height))
          (cut-width (* 1 mesh-size)))
      (format t "Crack width:~F~%" (* 1d0 cut-width))
      (cl-mpm/setup::remove-sdf
       *sim*
       (cl-mpm/setup::rectangle-sdf
        (list
         0d0
         (+ (second offset) 0d0)
         )
        (list
         cut-width
         ;; 0.0102
         cut-depth
         ))
       :refine 1))
    (format t "Total weight ~F~%"
            (loop for mp across (cl-mpm:sim-mps *sim*)
                  sum (* 9.8d0 (cl-mpm/particle:mp-mass mp))))

    (defparameter *current-load* 0d0)

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
     (mapcar (lambda (x) (* x -1d3)) *data-displacement*) (mapcar (lambda (x) (* x 0.013)) *data-load*) "load"
     ;; (mapcar (lambda (x) (* x -1d3)) *data-displacement*) (mapcar (lambda (x) (* x 0.013)) *data-mp-load*) "mpm-force"
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



(defun output-disp-header (output-dir)
  (with-open-file (stream (merge-pathnames "disp.csv" output-dir) :direction :output :if-exists :supersede)
    (format stream "disp,load~%")))

(defparameter *disp* 0d0)
(defun output-disp-data (output-dir)
  (with-open-file (stream (merge-pathnames "disp.csv" output-dir) :direction :output :if-exists :append)
    (format stream "~f,~f~%"
            *displacement*
            (* (cl-mpm/penalty::resolve-load *penalty*) 2d0))))
(defun test ()
  (defparameter *displacement* 0d0)
  (let* ((lstps 10)
         (total-disp -0.2d-3)
         (output-dir (format nil "./output-se/")))
    (setup :refine 1 :mps 3)
    (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) nil)
    (setf cl-mpm/damage::*enable-reflect-x* t)
    (setf (cl-mpm::sim-gravity *sim*) 0d4)
    (defparameter *data-displacement* '(0d0))
    (defparameter *data-load* '(0d0))
    (cl-mpm/setup::set-mass-filter *sim* 2.2d3 :proportion 1d-15)
    (setf (cl-mpm::sim-nonlocal-damage *sim*) t)
    (setf *disp* 0d0)
    (cl-mpm/dynamic-relaxation::run-adaptive-load-control
     *sim*
     :output-dir output-dir
     :plotter (lambda (sim)
                ;; (plot-domain sim)
                (plot-load-disp)
                (format t "Load ~E ~%" (cl-mpm/penalty::resolve-load *penalty*)))
     :loading-function (lambda (percent)
                         (setf *displacement* (* total-disp percent))
                         (let ((delta (- *displacement* *last-pos*)))
                           (cl-mpm/penalty::bc-increment-center *penalty* (cl-mpm/utils:vector-from-list (list 0d0 delta 0d0)))
                           (setf *last-pos* *displacement*))
                         )
     :pre-step (lambda ()
                 (output-disp-header output-dir)
                 (output-disp-data output-dir))
     :post-conv-step (lambda (sim)
                       ;;Save data
                       (push (* 2d0 (cl-mpm/penalty::resolve-load *penalty*)) *data-load*)
                       (push *displacement* *data-displacement*)
                       (output-disp-data output-dir))
     :load-steps lstps
     :enable-damage t
     :damping 1d0
     :substeps 20
     :criteria 1d-3
     :max-adaptive-steps 10
     :save-vtk-dr t
     :save-vtk-loadstep t
     :max-damage-inc 0.6d0
     :dt-scale 1d0)))
