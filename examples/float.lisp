(defpackage :cl-mpm/examples/float
  (:use :cl))
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)
(setf *block-compile-default* t)
(in-package :cl-mpm/examples/float)
(declaim (optimize (debug 0) (safety 0) (speed 3)))

;; (pushnew :cl-mpm-fbar *features*)
;; (remove :cl-mpm-fbar *features*)
;; (setf *features* (delete :cl-mpm-pic *features*))
;; (asdf:compile-system :cl-mpm :force T)
;; (asdf:compile-system :cl-mpm)

(defmethod cl-mpm/particle::constitutive-model ((mp cl-mpm/particle::particle-elastic-damage) strain dt)
  "Strain intergrated elsewhere, just using elastic tensor"
  (with-accessors ((de cl-mpm/particle::mp-elastic-matrix)
                   (d cl-mpm/particle::mp-damage)
                   (p cl-mpm/particle::mp-pressure)
                   (def cl-mpm/particle::mp-deformation-gradient))
      mp
    ;; (cl-mpm/fastmath:fast-.+
    ;;  (cl-mpm/fastmath:fast-scale-voigt
    ;;   (cl-mpm/constitutive::linear-elastic-mat strain de) (- 1d0 d))
    ;;  ;; (cl-mpm/utils:voigt-zeros)
    ;;  (cl-mpm/utils::voigt-eye (* 1d0 p (cl-mpm/fastmath:det def))))

    (cl-mpm/constitutive::linear-elastic-mat strain de)
    )
  )




(defun length-from-def (sim mp dim)
  (let* ((mp-scale 2)
         (h-initial (magicl:tref (cl-mpm/particle::mp-domain-size mp) dim 0))
         ;(h-initial  (/ (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) mp-scale))
         )
    h-initial
    ;(* (magicl:tref (cl-mpm::mp-deformation-gradient mp) dim dim) h-initial)
    ))
(defun max-stress (mp)
  (multiple-value-bind (l v) (magicl:hermitian-eig (cl-mpm::voight-to-matrix (cl-mpm/particle:mp-stress mp)))
    (apply #'max l)
    (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0)))

(defparameter *water-height* 0d0)
(defun plot (sim &optional (plot :deformed))
  (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:format-plot t "set object 1 rect from 0,0 to ~f,~f fc rgb 'blue' fs transparent solid 0.5 noborder behind" ms-x *water-height*))
  (cl-mpm/plotter:simple-plot
   sim
   :plot :deformed
   ;; :colour-func (lambda (mp) (cl-mpm/utils:varef (cl-mpm/particle:mp-stress mp) 1))
   :colour-func #'cl-mpm/particle::mp-boundary
   ))



(defun create-new-mp (pos h-x h-y density gravity-angle)
  (cl-mpm::make-particle 2
                         'cl-mpm/particle::particle-viscoplastic-damage
                         ;; 'cl-mpm/particle::particle-glen
                         :E 1d8
                         :nu 0.3250d0

                         :visc-factor 11d6
                         :visc-power 3d0

                         :initiation-stress 0.33d6
                         :damage-rate 1d-18
                         :critical-damage 0.5d0
                         :local-length 50d0

                         :gravity -9.8d0
                         :gravity-axis (magicl:from-list (list (cos (* gravity-angle (/ pi 180)))
                                                               (sin (* gravity-angle (/ pi 180)))) '(2 1))
                         :volume (* h-x h-y)
                         :position pos
                         :mass (* h-x h-y density)
                         :size (magicl:from-list (list h-x h-y) '(2 1) :type 'double-float)
                         :size-0 (magicl:from-list (list h-x h-y) '(2 1) :type 'double-float)
                         ))



(declaim (notinline setup-test-column))
(defun setup-test-column (size block-size block-offset &optional (e-scale 1d0) (mp-scale 1d0)
                          &rest mp-args)
  (let* ((sim (cl-mpm/setup::make-block (/ 1 e-scale)
                                        (mapcar (lambda (s) (* s e-scale)) size)
                                        ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         ;(e-scale 1)
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 500d0)
         (gravity-angle -80)
         ;; (mass (/ (* 900 h-x h-y) (expt mp-scale 2)))
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size)))
    (progn
      (let ((block-position
              ;; (mapcar #'+ (list (* h-x (- (+ (/ 1 (* 2 mp-scale))) 0))
              ;;                   (* h-y (+ (/ 1d0 (* 2d0 mp-scale)))))
              ;;         block-offset)
              block-offset
              ))
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-block-mps
               block-position
               block-size
               (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
               density
               ;; ;; 'cl-mpm::make-particle
               ;; ;; 'cl-mpm/particle::particle-elastic-logspin
               ;; 'cl-mpm/particle::particle-elastic-damage
               ;; ;; 'cl-mpm/particle::particle-elastic
               ;; ;; 'cl-mpm/particle::particle-viscoplastic-damage
               ;; ;; 'cl-mpm/particle::particle-glen
               ;; :E 1d9
               ;; :nu 0.3250d0

               ;; :damage (- 1d0 0.999d0)
               ;; ;; :visc-factor 11d6
               ;; ;; :visc-power 3d0

               ;; ;; :initiation-stress 0.33d6
               ;; ;; :damage-rate 1d-25
               ;; ;; :critical-damage 0.5d0
               ;; ;; :local-length 50d0

                'cl-mpm/particle::particle-visco-elasto-plastic-damage
                :E 1d9
                :nu 0.325d0

                :friction-angle 30.0d0

                :kt-res-ratio 1d-10
                :kc-res-ratio 1d-2
                :g-res-ratio 1d-3

                ;; :fracture-energy 3000d0
                ;; :initiation-stress stress
                :damage 0d0

                ;; :delay-time 1d2
                ;; :delay-exponent 3d0
                ;; :ductility ductility
                ;; :ductility-mode-2 ductility-ii
                ;; :critical-damage 1d0;(- 1.0d0 1d-3)
                :damage-domain-rate 0.0d0;This slider changes how GIMP update turns to uGIMP under damage
                ;; :local-length length-scale;(* 0.20d0 (sqrt 7))
                ;; :local-length-damaged 10d-10

                :enable-plasticity nil
                :psi (* 00d0 (/ pi 180))
                :phi (* 40d0 (/ pi 180))
                :c 1000d3


               :gravity -9.8d0
               :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
               ;; :gravity-axis (magicl:from-list (list (cos (* gravity-angle (/ pi 180)))
               ;;                                       (sin (* gravity-angle (/ pi 180)))) '(2 1))
               :index 0
               )))

      (setf (cl-mpm:sim-damping-factor sim) (* 0d-2 (cl-mpm/setup:estimate-critical-damping sim)))
      (setf (cl-mpm:sim-mass-filter sim) 1d-15)
      (setf (cl-mpm::sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-mass-scale sim) 1d0)
      (setf (cl-mpm:sim-dt sim) 1d-2)
      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-var
             (cl-mpm:sim-mesh sim)
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0 nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0 nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil nil 0)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil nil 0)))
             ))
      (let ((ocean-y 200d0)
            (ocean-angle gravity-angle)
            (flow-rate 10))
        (defparameter *water-bc*
          (cl-mpm/buoyancy::make-bc-buoyancy-clip
           sim
           ocean-y
           1000d0
           (lambda (pos datum)
             ;; t;(>= (magicl:tref pos 1 0) h-y)
             t
             ;; (< (cl-mpm/utils:varef pos 1) datum)
             )))
        (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 1d-5)
        (setf *water-height* ocean-y)
        (cl-mpm:add-bcs-force-list
         sim
         *water-bc*))
      sim)))

;Setup
(defun setup (&key (refine 1) (mps 2))
  (defparameter *run-sim* nil)
  (let* ((mesh-size (/ 20 refine))
        (mps-per-cell mps)
        (block-size (list 20 20))
        (block-offset (list 0 200))
        (domain-size (list 20 500))
        )
    (defparameter *sim* (setup-test-column 
                         domain-size
                                           block-size
                                           block-offset (/ 1 mesh-size) mps-per-cell))
    (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
    ;; (defparameter *sim* (setup-pressure '(200 200) '(100 200) '(0 0) (/ 1 mesh-size) mps-per-cell))
    )
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
  (print (cl-mpm:sim-dt *sim*))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  (defparameter *x* 0d0)
  (defparameter *x-pos* '())
  (defparameter *s-xx* '())
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
  ;; (increase-load *sim* *load-mps* 1)
  ;; (increase-load *sim* *load-mps* 100)
  )


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

(defun stop ()
  (setf *run-sim* nil))
(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)
  (defparameter *run-sim* t)
    ;; (vgplot:close-all-plots)
  ;;   (vgplot:figure)
  ;; (sleep 1)
  ;; (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
  ;;        (ms-x (first ms))
  ;;        (ms-y (second ms))
  ;;        )
  ;;   (vgplot:axis (list 0 ms-x
  ;;                      0 ms-y))
  ;;   (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
  ;;   (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
  ;;     (vgplot:format-plot t "set ytics ~f" h)
  ;;     (vgplot:format-plot t "set xtics ~f" h))
  (vgplot:figure)
  (with-open-file (stream (merge-pathnames "output/terminus_position.csv") :direction :output :if-exists :supersede)
    (format stream "Time (s),Terminus position~%"))

  (let* ((target-time 0.5d0)
         (dt (cl-mpm:sim-dt *sim*))
         (dt-scale 1d0)
         (substeps (floor target-time dt)))
    (cl-mpm::update-sim *sim*)
    (let* ((dt-e (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
           (substeps-e (floor target-time dt-e)))
      (format t "CFL dt estimate: ~f~%" dt-e)
      (format t "CFL step count estimate: ~D~%" substeps-e)
      (setf (cl-mpm:sim-dt *sim*) dt-e)
      (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
    (time (loop for steps from 0 to 40
                while *run-sim*
                do
                   (progn
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     ;; (cl-mpm/output:save-csv (merge-pathnames (format nil "output/sim_csv_~5,'0d.csv" *sim-step*)) *sim*)

                     (push *t* *time*)
                     ;; (let ((cfl (find-max-cfl *sim*)))
                     ;;   (push cfl *cfl-max*))
                     (setf *x*
                           (loop for mp across (cl-mpm:sim-mps *sim*)
                                 maximize (magicl:tref (cl-mpm/particle::mp-displacement mp) 1 0)))
                     (push
                      *x*
                      *x-pos*)
                     (push
                      (/
                       (loop for mp across (cl-mpm:sim-mps *sim*)
                              sum (magicl:tref (cl-mpm/particle::mp-stress mp) 0 0))
                       (length (cl-mpm:sim-mps *sim*)))
                      *s-xx*)
                     (with-open-file (stream (merge-pathnames "output/terminus_position.csv") :direction :output :if-exists :append)
                       (format stream "~f, ~f ~%" *t* *x*))
                     (let ((cfl 0))
                       (time (dotimes (i substeps)
                               ;; (increase-load *sim* *terminus-mps*
                               ;;                (magicl:from-list (list (* (cl-mpm:sim-dt *sim*) 1d4) 0d0) '(2 1)))
                               ;; (pescribe-velocity *sim* *terminus-mps* '(1d0 nil))
                               (cl-mpm::update-sim *sim*)
                               ;; (setf cfl (max cfl (find-max-cfl *sim*)))
                               (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))
                       ;; (setf cfl (find-max-cfl *sim*))
                       (format t "CFL: ~f~%" cfl)
                       (push cfl *cfl-max*)
                         ;; (multiple-value-bind (dt-e substeps-e) (calculate-dt cfl 1d-4 target-time)
                         ;;   (format t "CFL dt estimate: ~f~%" dt-e)
                         ;;   (format t "CFL step count estimate: ~D~%" substeps-e)
                         ;;   (push cfl *cfl-max*)
                         ;;   (setf (cl-mpm:sim-dt *sim*) dt-e)
                         ;;   (setf substeps substeps-e)
                         ;;   )
                       ;; (let* ((dt-e (cl-mpm::calculate-min-dt *sim*))
                       ;;        (substeps-e (/ target-time dt-e)))
                       ;;   (setf substeps-e (min 1000 substeps-e))
                       ;;   (format t "CFL dt estimate: ~f~%" dt-e)
                       ;;   (format t "CFL step count estimate: ~D~%" substeps-e)
                       ;;   (setf (cl-mpm:sim-dt *sim*) dt-e)
                       ;;   (setf substeps substeps-e))
                         )
                     ;; (setf (cl-mpm:sim-damping-factor *sim*) (max 0.1d0 (/ (cl-mpm:sim-damping-factor *sim*) 1.1d0)))
                     ;; (setf (cl-mpm::sim-enable-damage *sim*) t)
                     (incf *sim-step*)
                     ;; (plot *sim*)
                     (vgplot:plot *time* *x-pos*)
                     (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                        :terminal "png size 1920,1080"
                                        )
                     (swank.live:update-swank)
                     (sleep .01)

                     ))))
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
  ;; (plot-s-xx *far-field-mps* 125d0)
  (plot-disp)
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
      (loop for yi in y
            for si in s-xx
            for sai in s-an
            do
               (format stream "~f, ~f, ~f ~%" yi si sai))))
  )
(defun plot-disp ()
  (let* ()
    ;; (vgplot:figure)
    (vgplot:title "S_{xx} time")
    (vgplot:xlabel "Time (s)")
    (vgplot:ylabel "S_{xx} (m)")
    (vgplot:plot *time* *x-pos*)
    ))
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

(defun run-conv ()
  (vgplot:close-all-plots)
  (setf *run-sim* t)
  (loop for r in (list 1 2 4 8)
        while *run-sim*
        do (progn
             (setup :refine r :mps 2)
             (run)
             (vgplot:title (format nil "~D" r))
             (vgplot:print-plot
              (merge-pathnames (format nil "refine_~5,'0d.png" r))
              :terminal "png size 1920,1080"))))
