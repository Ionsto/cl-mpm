(defpackage :cl-mpm/examples/dr/quasi-time
  (:use :cl
   :cl-mpm/example
   :cl-mpm/utils)
  )
(in-package :cl-mpm/examples/dr/quasi-time)

;; (setf cl-mpm/settings:*optimise-setting* cl-mpm/settings::*optimise-speed*)
;; (sb-ext:restrict-compiler-policy 'speed  3 3)
;; (sb-ext:restrict-compiler-policy 'debug  0 0)
;; (sb-ext:restrict-compiler-policy 'safety 0 0)
(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)
(setf cl-mpm/settings:*optimise-setting* cl-mpm/settings::*optimise-debug*)


(defclass particle-voldeg-delayed (cl-mpm/particle::particle-elastic-damage-delayed)
  ())

(defmethod cl-mpm/particle::post-damage-step ((mp particle-voldeg-delayed) dt)
  ;; (apply-isotropic-degredation mp)
  (cl-mpm/damage::apply-vol-degredation mp)
  )
(defmethod cl-mpm::update-particle (mesh (mp particle-voldeg-delayed) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  (cl-mpm::update-domain-stretch mesh mp dt)
  ;; (cl-mpm::update-domain-polar mesh mp dt)
  ;; (cl-mpm::update-domain-corner mesh mp dt)
  )
(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-vm) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  ;; (cl-mpm::update-domain-polar mesh mp dt)
  ;; (cl-mpm::update-domain-polar mesh mp dt)
  ;; (cl-mpm::update-domain-corner mesh mp dt)
  )

(defun plot-domain ()
  (when *sim*
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :colour-func #'cl-mpm/particle::mp-damage
     )))

;;Can be :K-0 :K-UPDATED :P-ELASTIC P-PLASTIC-SCALAR P-PLASTIC-TANGENT


(defun setup (&key (refine 1) (mps 2))
  (defparameter *sim* nil)
  (let ((mps-per-dim mps)
        (size '(32 16))
        (block-size '(8 8))
        )
    (let* ((E 1d6)
           (density 1d3)
           (sim (cl-mpm/setup::make-simple-sim
                 (/ 1d0 refine)
                 (mapcar (lambda (x) (* x refine)) size)
                 ;; :sim-type 'mpm-sim-test
                 :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul
                 :args-list
                 (list
                  :enable-aggregate t
                  :enable-fbar t
                  :max-split-depth 6
                  ;; :mp-removal-size t
                  :enable-damage t
                  :enable-length-localisation t
                  :split-factor (/ 1.1d0 mps)
                  :enable-split nil
                  :gravity -10d0)))
           (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim))))
      (declare (double-float h density))
      (setf *sim* sim)
      (let* ((offset (* 2d0 h))
             (E 1d6)
             (length-scale 1d0)
             (init-stress 1d4)
             (gf 1000d0)
             (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress E)))
        (format t "Ductility ~E~%" ductility)
        (when (<= ductility 1d0)
          (error "Ductility too small"))
        (cl-mpm::add-mps
         sim
         (cl-mpm/setup::make-block-mps
          (list 0d0 offset 0d0)
          block-size
          (mapcar (lambda (e) (* (/ e h) mps)) block-size)
          density
          ;; 'cl-mpm/particle::particle-vm
          ;; :E E
          ;; :nu 0.3d0
          ;; :rho 20d3
          ;; :rho-r 2d3
          ;; :softening 0.1d0
          'particle-voldeg-delayed
          :E E
          :nu 0.3d0
          :initiation-stress init-stress
          :local-length length-scale
          :ductility ductility
          :delay-time 5d1
          :delay-exponent 1d0
          :residual-strength (- 1d0 1d-9)
          ))
        (cl-mpm::domain-sort-mps sim)
        (let ((dsize (first size))
              (eps-scale 1d1)
              (mu 0.05d0)
              )
          (defparameter *floor*
            (cl-mpm/penalty::make-bc-penalty-distance-point
             *sim*
             (cl-mpm/utils:vector-from-list (list 0d0 1d0 0d0))
             (cl-mpm/utils:vector-from-list (list (* dsize 0.5d0) offset 0d0))
             (* dsize 0.5d0)
             (* E eps-scale)
             mu
             0d0
             )))
        (when (> offset 0d0)
          (cl-mpm::add-bcs-force-list
           *sim*
           *floor*))
        (defparameter *density* density)
        (cl-mpm/setup::set-mass-filter sim density :proportion 1d-15))))
  (cl-mpm/setup::setup-bcs
   *sim*
   :left (list 0 nil nil)
   :bottom (list 0 0 nil)
   :top (list nil nil nil))
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defun run (&key (output-dir "./output-qs/")
              (dt-scale 1d0)
              (dt 0.5d0)
              (adaptive-steps 0)
              (min-adaptive-steps 0)
              (max-damage 1.1d0)
              )
  (let ((step 0))
    (cl-mpm/dynamic-relaxation::run-quasi-time
     *sim*
     :output-dir output-dir
     :plotter (lambda (sim) (plot-domain)
                (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step))
                                   :terminal "png size 1920,1080")
                (incf step))
     :damping (sqrt 2d0)
     :dt dt
     :total-time 90d0
     :conv-criteria 1d-6
     :max-adaptive-steps adaptive-steps
     :min-adaptive-steps min-adaptive-steps
     :max-damage-inc max-damage
     :substeps 20
     :save-vtk-dr t
     :save-vtk-loadstep t
     :dt-scale dt-scale)))

(defun run-explicit (&key (output-dir "./output/") (dt-scale 1d0)
                       (mass-scale 1d0)
                       (dt 0.5d0))
  (let ((step 0))
    (change-class *sim* 'cl-mpm/damage::mpm-sim-agg-damage)
    (cl-mpm/dynamic-relaxation::run-time
     *sim*
     :output-dir output-dir
     :plotter (lambda (sim) (plot-domain)
                (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step))
                                   :terminal "png size 1920,1080")
                (incf step)
                )
     :post-conv-step (lambda (sim)
                       (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
                       (setf (cl-mpm::sim-velocity-algorithm *sim*) :TBLEND)
                       (setf (cl-mpm::sim-mass-scale *sim*) mass-scale))
     :conv-criteria 1d-6
     :damping 0.01d0
     :dt dt
     :total-time 90d0
     :save-vtk-dr t
     :save-vtk-loadstep t
     :dt-scale 0.5d0)))
(defun run-adaptive-mass-explicit (&key (output-dir "./output/") (dt-scale 1d0)
                       (mass-scale 1d0)
                       (dt 0.5d0))
  (let ((step 0))
    (change-class *sim* 'cl-mpm/damage::mpm-sim-agg-damage)
    (cl-mpm/dynamic-relaxation::run-time
     *sim*
     :output-dir output-dir
     :plotter (lambda (sim) (plot-domain)
                (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step))
                                   :terminal "png size 1920,1080")
                (incf step)
                )
     :post-conv-step (lambda (sim)
                       (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
                       (setf (cl-mpm::sim-velocity-algorithm *sim*) :TBLEND)
                       (setf (cl-mpm::sim-mass-scale *sim*) mass-scale))
     :adaptive-mass-enabled t
     :adaptive-mass-scale 1d0
     :adaptive-hist 2d0
     :initial-lstps 2
     :conv-criteria 1d-6
     :damping 0.01d0
     :dt dt
     :total-time 90d0
     :save-vtk-dr t
     :save-vtk-loadstep t
     :dt-scale 0.5d0))
  )

;; (defun run-adaptive-mass-explicit (&key (output-dir "./output/") (dt-scale 1d0)
;;                        (mass-scale 1d0)
;;                        (dt 0.5d0))
;;   (let ((step 0))
;;     (cl-mpm/dynamic-relaxation::elastic-static-solution *sim* :crit 1d-6)
;;     (change-class *sim* 'cl-mpm/damage::mpm-sim-agg-damage)
;;     (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
;;     (setf (cl-mpm::sim-velocity-algorithm *sim*) :TBLEND)
;;     (setf (cl-mpm::sim-mass-scale *sim*) mass-scale)
;;     (uiop:ensure-all-directories-exist (list (merge-pathnames output-dir)))
;;     (cl-mpm/dynamic-relaxation::step-real-time
;;      *sim*
;;      0
;;      :output-dir output-dir
;;      :plotter (lambda (sim) (plot-domain)
;;                 (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step))
;;                                    :terminal "png size 1920,1080")
;;                 (incf step))
;;      ;; :post-conv-step (lambda (sim))
;;      :criteria 1d-3
;;      :damping 0.01d0
;;      :target-time dt
;;      :step-time 90d0
;;      :enable-mass-scaling t
;;      :mass-scaler 1d2
;;      ;; :save-vtk-dr t
;;      ;; :save-vtk-loadstep t
;;      :dt-scale 0.5d0)))

(defun test ()
  (cl-mpm::set-workers 8)
  (setup :refine 1 :mps 4)
  (run)
  )

(defun test-explicit ()
  (cl-mpm::set-workers 8)
  (setup :refine 1 :mps 4)
  (run-explicit))



(defun test-all ()
  (let ((r 1)
        (mps 6)
        (dt 0.5d0))
    (cl-mpm::set-workers 8)
    (setup :refine r :mps mps)
    (run :output-dir "./output-qs/" :dt dt)
    (setup :refine r :mps mps)
    (run-explicit :output-dir "./output-ms-1/" :mass-scale 1d0)
    (setup :refine r :mps mps)
    (run-explicit :output-dir "./output-ms-10/" :mass-scale 10d0)
    (setup :refine r :mps mps)
    (run-explicit :output-dir "./output-ms-100/" :mass-scale 100d0)
    (setup :refine r :mps mps)
    (run-adaptive-mass-explicit :output-dir "./output-ms-adaptive/" :mass-scale 1d0)
    ))

(defun test-adaptive ()
  (let ((r 1)
        (mps 6)
        (dt 0.5d0))
    (cl-mpm::set-workers 8)
    (setup :refine r :mps mps)
    (run :output-dir "./output-constant-10/" :max-damage 1.1d0 :dt 10d0)
    ;; (setup :refine r :mps mps)
    ;; (run :output-dir "./output-adaptive-0.5/"
    ;;      :adaptive-steps 4
    ;;      :min-adaptive-steps -4
    ;;      :max-damage 0.5d0)
    ;; (setup :refine r :mps mps)
    ;; (run :output-dir "./output-adaptive-0.1/"
    ;;      :adaptive-steps 4
    ;;      :min-adaptive-steps -4
    ;;      :max-damage 0.1d0)
    )
  )
