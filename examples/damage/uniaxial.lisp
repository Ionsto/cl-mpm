(defpackage :cl-mpm/examples/damage/uniaxial
  (:use :cl
        :cl-mpm/example))
(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)
;; (sb-ext:restrict-compiler-policy 'speed  3 3)
;; (sb-ext:restrict-compiler-policy 'debug  0 0)
;; (sb-ext:restrict-compiler-policy 'safety 0 0)
;; (setf *block-compile-default* t)
(in-package :cl-mpm/examples/damage/uniaxial)

(declaim (optimize (debug 0) (safety 0) (speed 3)))

(declaim (notinline plot-domain))
(defun plot-domain (sim)
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :deformed
   :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp))))
(declaim (notinline plot))
(defun plot (sim)
  (vgplot:plot
   (reverse *data-displacement*)
   (reverse *data-load*)
   "MPM"
   )
  ;; (plot-domain sim)
  )
;; (defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-elastic-damage) dt)
;;   (with-accessors ((y cl-mpm/particle::mp-damage-y-local)
;;                    (strain cl-mpm/particle::mp-strain)
;;                    (init-stress cl-mpm/particle::mp-initiation-stress)
;;                    (ybar cl-mpm/particle::mp-damage-ybar)
;;                    (de cl-mpm/particle::mp-elastic-matrix)
;;                    (def cl-mpm/particle::mp-deformation-gradient)
;;                    (E cl-mpm/particle::mp-e)
;;                    (nu cl-mpm/particle::mp-nu)
;;                    (k cl-mpm/particle::mp-compression-ratio))
;;       mp
;;     (declare (double-float E nu))
;;     (progn
;;       (let* ((angle 60d0)
;;              (stress (cl-mpm/fastmaths:fast-scale!
;;                       (cl-mpm/constitutive:linear-elastic-mat strain de)
;;                       1d0
;;                       ;; (/ 1d0 (magicl:det def))
;;                       )))
;;         (setf y
;;               (cl-mpm/damage::tensile-energy-norm strain E de)
;;               ;; (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile stress (* angle (/ pi 180d0)))
;;               ;; (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile stress (* angle (/ pi 180d0)))
;;               ;; (cl-mpm/damage::criterion-modified-vm strain k E nu)
;;               )))))

(declaim (notinline setup))
(defun setup (&key (refine 1) (mps 2)
                (epsilon-scale 1d2)
                (gf 100d0)
                (gf-scale 1d0))
  (let* ((h (/ 1d0 refine))
         (L 10d0)
         (domain-size (list (* L 2)))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list L))
         (density 1d3)
         (E 1d9)
         (init-stress 1d5)
         (length-scale 1d0)
         (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 (* gf-scale gf) length-scale init-stress E))
         ;; (ductility 50d0)
         )
    (when (<= ductility 1d0)
      (error "Ductility must be greater than 1"))
    (format t "Ductility ~F~%" ductility)
    (setf *sim* (cl-mpm/setup::make-simple-sim
                 h
                 element-count
                 :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul
                 :args-list (list :enable-aggregate t
                                  :enable-fbar nil
                                  ;; :ghost-factor (* E 1d-3)
                                  )))
    (defparameter *h* h)
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0d0)
      block-size
      (mapcar (lambda (e) (* (/ e h) mps)) block-size)
      density
      'cl-mpm/particle::particle-elastic-damage
      :E E
      :nu 0d0
      ;; ;;Material parameter
      :initiation-stress init-stress
      :local-length length-scale
      :ductility ductility
      :residual-strength (- 1d0 1d-9)
      ;:residual-strength 1d0
      ;; :residual-strength 1d0;(- 1d0 1d-12)
      ))
    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (when (<= (abs (- (cl-mpm/utils:varef (cl-mpm/particle::mp-position mp) 0) (* 0.5d0 L)))
                 (* h 0.5d0))
         (setf (cl-mpm/particle::mp-initiation-stress mp)
               (* 0.9d0 (cl-mpm/particle::mp-initiation-stress mp))))))
    (setf (cl-mpm:sim-gravity *sim*) 0d0)
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)
    (cl-mpm/setup::setup-bcs
     *sim*
     ;:left '(0 nil nil)
     )
    (defparameter *penalty*
      (cl-mpm/penalty::make-bc-penalty-displacment
       *sim*
       (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0))
       (* E epsilon-scale)))
    (cl-mpm:add-bcs-force-list
     *sim*
     *penalty*)
    ))

(defun setup-plastic (&key (refine 1) (mps 2)
                        (epsilon-scale 1d2)
                        (gf-scale 1d0))
  (let* ((h (/ 1d0 refine))
         (L 10d0)
         (domain-size (list (* L 2)))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list L))
         (density 1d3)
         (E 1d9)
         (gf 100d0)
         (init-stress 1d5)
         (length-scale 1d0)
         (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress E gf-scale))
         ;; (ductility 10d0)
         )
    (when (<= ductility 1d0)
      (error "Ductility must be greater than 1"))
    (format t "Ductility ~F~%" ductility)
    (setf *sim* (cl-mpm/setup::make-simple-sim
                 h
                 element-count
                 :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul
                 :args-list (list :enable-aggregate t
                                  :enable-fbar t)))
    (defparameter *h* h)
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0d0)
      block-size
      (mapcar (lambda (e) (* (/ e h) mps)) block-size)
      density
      'cl-mpm/particle::particle-vm
      :E E
      :nu 0d0
      :rho (* (sqrt 2/3) init-stress)
      ;; :rho-r 1d0
      ;; :softening 1d0
      ))
    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (when (<= (abs (- (cl-mpm/utils:varef (cl-mpm/particle::mp-position mp) 0) (* 0.5d0 L)))
                 (* h 0.5d0))
         (setf (cl-mpm/particle::mp-rho mp)
               (* 0.9d0 (cl-mpm/particle::mp-rho mp))))))
    (setf (cl-mpm:sim-gravity *sim*) 0d0)
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)
    (cl-mpm/setup::setup-bcs
     *sim*
     ;:left '(0 nil nil)
     )
    (defparameter *penalty*
      (cl-mpm/penalty::make-bc-penalty-displacment
       *sim*
       (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0))
       (* E epsilon-scale)))
    (cl-mpm:add-bcs-force-list
     *sim*
     *penalty*)
    ))


(defun get-load ()
  (cl-mpm/penalty::resolve-load *penalty*))

(defun output-disp-header (output-dir)
  (with-open-file (stream (merge-pathnames "disp.csv" output-dir) :direction :output :if-exists :supersede)
    (format stream "disp,load~%")))

(defun output-disp-data (output-dir)
  (with-open-file (stream (merge-pathnames "disp.csv" output-dir) :direction :output :if-exists :append)
    (format stream "~f,~f~%"
            *displacement*
            (get-load))))

(defun run (&key (output-dir (format nil "./output/")))
  (let ((lstps 40)
        (total-disp 5d-3))
  (defparameter *displacement* 0d0)
    (defparameter *data-displacement* '(0d0))
    (defparameter *data-load* '(0d0))
    (vgplot:close-all-plots)
    (cl-mpm/dynamic-relaxation::run-adaptive-load-control
     *sim*
     :output-dir output-dir
     :plotter
     (lambda (sim)
       (plot *sim*)
       (format t "Load ~E ~%" (cl-mpm/penalty::resolve-load *penalty*)))
     :loading-function
     (lambda (percent)
       (setf *displacement* (* total-disp percent))
       (cl-mpm/penalty::bc-set-displacement
        *penalty*
        (cl-mpm/utils:vector-from-list (list *displacement* 0d0 0d0))))
     :pre-step
     (lambda ()
       (output-disp-header output-dir)
       (output-disp-data output-dir))
     :post-conv-step
     (lambda (sim)
       (push (get-load) *data-load*)
       (push *displacement* *data-displacement*)
       (output-disp-data output-dir))
     :load-steps lstps
     :enable-damage t
     :enable-plastic nil
     :damping (sqrt 2d0)
     :substeps 50
     :criteria 1d-6
     :max-adaptive-steps 0
     :save-vtk-dr t
     :save-vtk-loadstep t
     :max-damage-inc 1.1d0
     :true-stagger nil
     :dt-scale 0.4d0)))


(defun list-interp (p points)
  (let (;; (points (list 0d0 1d0 -1d0))
        )
    (flet ((interp (p0 p1 dx)
             (+ (* p0 (- 1 dx))
                (* p1 dx))))
      (let* ((point-count (- (length points) 1))
             (point-gap (/ 1d0 point-count))
             (p0 (floor p point-gap))
             (p1 (ceiling p point-gap)))
        (interp (nth p0 points) (nth p1 points) (/ (- p (* p0 point-gap)) point-gap))))))

(defun load-unload (p dist)
  (list-interp p (list 0d0 dist (* -0.5 dist))))

(defun multi-load-unload (p dist)
  (let ((factor 1.25))
    (list-interp p (list 0d0
                         dist
                         (* -0.5 dist)
                         (* factor dist)
                         (* -0.5 dist)
                         (* (expt factor 2) dist)
                         (* -0.5 dist)))))

(defun run-mcr (&key (output-dir (format nil "./output/")))
  (let ((lstps 200)
        (total-disp 2d-3))
    (defparameter *displacement* 0d0)
    (defparameter *data-displacement* '(0d0))
    (defparameter *data-load* '(0d0))
    (vgplot:close-all-plots)
    (cl-mpm/dynamic-relaxation::run-adaptive-load-control
     *sim*
     :output-dir output-dir
     :plotter
     (lambda (sim)
       (plot *sim*)
       (format t "Load ~E ~%" (cl-mpm/penalty::resolve-load *penalty*)))
     :loading-function
     (lambda (percent)
       (setf *displacement* (multi-load-unload percent total-disp))
       ;(setf *displacement* (* total-disp percent))
       (cl-mpm/penalty::bc-set-displacement
        *penalty*
        (cl-mpm/utils:vector-from-list (list *displacement* 0d0 0d0))))
     :pre-step
     (lambda ()
       (output-disp-header output-dir)
       (output-disp-data output-dir))
     :post-conv-step
     (lambda (sim)
       (push (get-load) *data-load*)
       (push *displacement* *data-displacement*)
       (output-disp-data output-dir))
     :load-steps lstps
     :enable-damage t
     :enable-plastic nil
     :damping (sqrt 2)
     :substeps 50
     :criteria 1d-6
     :max-adaptive-steps 0
     :save-vtk-dr nil
     :save-vtk-loadstep t
     :max-damage-inc 1.1d0
     :dt-scale 0.9d0)))

(defun test ()
  (cl-mpm/utils:set-workers 8)
  (let* ((output-dir (format nil "./output-3/")))
    (setup :refine 4 :mps 3 :gf 150d0)
    (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
    ;; (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
    (run :output-dir output-dir)
    ))


(defun test-ekl ()
  (let* ((output-dir (format nil "./output-ekl/")))
    (setup :refine 4 :mps 3)
    ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) nil)
    (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
    (run :output-dir output-dir)
    ))

(defun test-ll ()
  (cl-mpm/utils::set-workers 8)
  (let* ((refine 5)
         (gf 200d0)
         (mps 3))
    ;; (setup :refine refine :mps mps :gf gf)
    ;; (run :output-dir "./output-standard/")
    ;; (get-all-interactions :output-dir "./output-standard/")

    (setup :refine refine :mps mps :gf gf)
    (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
    (run :output-dir "./output-ll/")
    (get-all-interactions :output-dir "./output-ll/")

    ;; (setup :refine refine :mps mps :gf (* 4d0 gf))
    ;; (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
    ;; (ignore-errors
    ;;  (run :output-dir "./output-ekl/"))
    ;; (cl-mpm::reset-loadstep *sim*)
    ;; (cl-mpm:update-sim *sim*)
    ;; (get-all-interactions :output-dir "./output-ekl/")

    ;; (setup :refine refine :mps mps :gf (* 4d0 gf))
    ;; (setf (cl-mpm/damage::sim-enable-stress-based-length *sim*) t)
    ;; (ignore-errors
    ;;  (run :output-dir "./output-sb/"))
    ;; (cl-mpm::reset-loadstep *sim*)
    ;; (cl-mpm:update-sim *sim*)
    ;; (get-all-interactions :output-dir "./output-sb/")
    ))

(defun test-refine ()
  (dolist (r (list 1 2 4))
    (let* ((output-dir (format nil "./output-~D/" r)))
      (setup :refine r :mps 2)
      ;; (setf (cl-mpm::sim-nonlocal-damage *sim*) nil)
      (run :output-dir output-dir)
      )))


(defclass cl-mpm/particle::particle-elastic-damage-1d-tcs (cl-mpm/particle::particle-elastic-damage-tcs)
  ())

(defclass cl-mpm/particle::particle-elastic-damage-1d-spectral (cl-mpm/particle::particle-elastic-damage)
  ())

;; (defclass cl-mpm/particle::particle-elastic-damage-1d-spectral-strain (cl-mpm/particle::particle-elastic-damage)
;;   ())

(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-elastic-damage-1d-tcs) dt)
  (cl-mpm/damage::apply-tcs-1d mp))

(defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-elastic-damage-1d-spectral) dt)
  (cl-mpm/damage::apply-spectral-1d mp))

;; (defmethod cl-mpm/particle::post-damage-step ((mp cl-mpm/particle::particle-elastic-damage-spectral-strain) dt)
;;   (cl-mpm/damage::apply-tensile-strain-degredation mp))

(defun test-mcr ()
  (dolist (particle (list
                     ;; 'cl-mpm/particle::particle-elastic-damage
                     'cl-mpm/particle::particle-elastic-damage-1d-tcs
                     ;; 'cl-mpm/particle::particle-elastic-damage-1d-spectral
                     ;; 'cl-mpm/particle::particle-elastic-damage-spectral-strain
                     ))

    (let* ()
      (setup :refine 1 :mps 3
             :gf 100d0)
      (cl-mpm::iterate-over-mps
       (cl-mpm:sim-mps *sim*)
       (lambda (mp)
         (change-class mp particle)))
      ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
      ;; (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
      (run-mcr :output-dir (format nil "./output-~A/" particle)
               )
      )))
(defun plot-nonlocal-inter ()
  (let* ((find-pos (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0)))
         (mp (cl-mpm/setup::find-mp *sim* find-pos)))
    (multiple-value-bind (pos weights) (cl-mpm/damage::get-nonlocal-interactions *sim* mp)
      (let ((x (loop for p in pos collect (cl-mpm/utils:varef p 0)))
            (y (loop for p in pos collect (cl-mpm/utils:varef p 1))))
        (loop for p in pos
              do (format t "~A ~A~%" (cl-mpm/utils:varef p 0) (cl-mpm/utils:varef p 1)))
        (vgplot:format-plot t "set xrange [~f:~f]" -20d0 20d0)
        (vgplot:format-plot t "set yrange [~f:~f]" -20d0 20d0)
        (vgplot:3d-plot x y weights ";;with points lc palette")
        (vgplot:xlabel "x")
        (vgplot:ylabel "y")
        ))))

(defun get-all-interactions (&key (output-dir "./output/"))
  (cl-mpm::iterate-over-mps-serial
   (cl-mpm:sim-mps *sim*)
   (lambda (mp)
     (multiple-value-bind (pos weights) (cl-mpm/damage::get-nonlocal-interactions *sim* mp)
       (let ((x (loop for p in pos collect (cl-mpm/utils:varef p 0)))
             (y (loop for p in pos collect (cl-mpm/utils:varef p 1))))
         (with-open-file (stream (merge-pathnames (format nil "interaction_~D.csv" (cl-mpm/particle::mp-unique-index mp)) output-dir) :direction :output :if-exists :supersede)
           (format stream "x,w~%")
           (loop for ix in x
                 for w in weights
                 do (format stream "~F,~F~%"
                            (+ ix (cl-mpm/utils:varef
                                   (cl-mpm/particle::mp-position mp)
                                   0))
                            w)))
         ;; (loop for p in pos
         ;;       do (format t "~A ~A~%" (cl-mpm/utils:varef p 0) (cl-mpm/utils:varef p 1)))
         ;; (vgplot:format-plot t "set xrange [~f:~f]" -20d0 20d0)
         ;; (vgplot:format-plot t "set yrange [~f:~f]" -20d0 20d0)
         ;; (vgplot:3d-plot x y weights ";;with points lc palette")
         ;; (vgplot:xlabel "x")
         ;; (vgplot:ylabel "y")
         ))
     )
   )
  )

(defun plot-nonlocal-inter ()
  (loop for px from 0d0 to 0d0
        do
           (let* ((find-pos (cl-mpm/utils:vector-from-list (list px 0d0 0d0)))
                  (mp (cl-mpm/setup::find-mp *sim* find-pos)))
             (multiple-value-bind (pos weights) (cl-mpm/damage::get-nonlocal-interactions-ekl *sim* mp)
               (let ((x (loop for p in pos collect (cl-mpm/utils:varef p 0)))
                     (y (loop for p in pos collect (cl-mpm/utils:varef p 1)))
                     )
                 (loop for p in pos
                       for w in weights
                       do (format t "~A ~A ~A~%" (cl-mpm/utils:varef p 0) (cl-mpm/utils:varef p 1) w))
                 (vgplot:format-plot t "set xrange [~f:~f]" -1d0 1d0)
                 ;; (vgplot:format-plot t "set yrange [~f:~f]" -20d0 20d0)
                 (vgplot:plot x weights ";;with points")
                 ;; (vgplot:3d-plot x y weights ";;with points lc palette")
                 (vgplot:xlabel "x")
                 (vgplot:ylabel "y")
                 )))
           (sleep 1)
        ))

(defun run-adaptive (&key (output-dir (format nil "./output/"))
                       (lstps 40)
                       (adaptive-steps 10)
                       (max-damage 1.1d0)
                       )
  (let ((total-disp 5d-3))
  (defparameter *displacement* 0d0)
    (defparameter *data-displacement* '(0d0))
    (defparameter *data-load* '(0d0))
    (vgplot:close-all-plots)
    (cl-mpm/dynamic-relaxation::run-adaptive-load-control
     *sim*
     :output-dir output-dir
     :plotter
     (lambda (sim)
       (plot *sim*)
       (format t "Load ~E ~%" (cl-mpm/penalty::resolve-load *penalty*)))
     :loading-function
     (lambda (percent)
       (setf *displacement* (* total-disp percent))
       (cl-mpm/penalty::bc-set-displacement
        *penalty*
        (cl-mpm/utils:vector-from-list (list *displacement* 0d0 0d0))))
     :pre-step
     (lambda ()
       (output-disp-header output-dir)
       (output-disp-data output-dir))
     :post-conv-step
     (lambda (sim)
       (push (get-load) *data-load*)
       (push *displacement* *data-displacement*)
       (output-disp-data output-dir))
     :load-steps lstps
     :enable-damage t
     :enable-plastic nil
     :damping (sqrt 2d0)
     :substeps 50
     :criteria 1d-6
     :max-adaptive-steps adaptive-steps
     :save-vtk-dr t
     :save-vtk-loadstep t
     :max-damage-inc max-damage
     :true-stagger nil
     :dt-scale 0.9d0)))

(defun test-adaptive ()
  (cl-mpm/utils:set-workers 8)
  (dolist (lstps (list 8 80))
    (setup :refine 4 :mps 3 :gf 150d0)
    (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
    (run-adaptive :output-dir (format nil "./output-~D/" lstps)
                  :lstps lstps)
    )
  (dolist (max-damage (list 0.1d0 0.5d0 0.9d0))
    (let ((lstps 8))
      (setup :refine 4 :mps 3 :gf 150d0)
      (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
      (run-adaptive :output-dir (format nil "./output-adaptive-~F/" max-damage)
                    :lstps lstps
                    :adaptive-steps 10
                    :max-damage max-damage))))
