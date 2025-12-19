(defpackage :cl-mpm/examples/bounce
  (:use :cl
   :cl-mpm/example
   :cl-mpm/utils))
(in-package :cl-mpm/examples/bounce)

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-elastic) dt fbar)
  ;(cl-mpm::update-stress-kirchoff-dynamic-relaxation mesh mp dt fbar)
  (cl-mpm::update-stress-linear mesh mp dt fbar)
  )
(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-elastic) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  ;; (cl-mpm::update-domain-det mesh mp dt)
  (cl-mpm::update-domain-stretch mesh mp dt)
  ;; (cl-mpm::update-domain-corner mesh mp dt)
  ;; (cl-mpm::co-domain-corner-2d mesh mp dt)
  ;; (cl-mpm::update-domain-polar-2d mesh mp dt)
  ;; (cl-mpm::scale-domain-size mesh mp)
  )

(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)
(declaim (optimize (debug 3) (safety 0) (speed 3)))
(defun plot (sim)
  (multiple-value-bind (x y)
    (loop for mp across (cl-mpm:sim-mps sim)
          collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
          collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
          finally (return (values x y)))
    (vgplot:plot x y ";;with points pt 7"))
  (vgplot:replot))
;; Setup
(defun setup (&key (mps 2) (refine 1))
  "Column setup"
  (let* ((L 10d0)
         (size (list 40 40))
         ;; (block-offset 20)
         (e-scale refine)
         (h (/ 1d0 e-scale))
         (block-size (list h L))
         (mp-scale mps)
         (sim (cl-mpm/setup::make-simple-sim
               (/ 1d0 e-scale)
               ;; (mapcar (lambda (x) (* x e-scale)) size)
               (list 1 (/ L h))
               ;; :sim-type 'cl-mpm/aggregate:mpm-sim-agg-usf
               :sim-type 'cl-mpm:mpm-sim-usf
               ;; :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic
               :args-list
               (append
                (list
                 ;; :enable-aggregate t
                 :vel-algo :FLIP
                 ))))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 1d4)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
         (offset (* h 0)))

    (declare (double-float h density))
    (format t "Element count ~D~%" (/ L h))
    (defparameter *sim* sim)
    (progn
      (let* ((E 1d7))
        (cl-mpm::add-mps
         sim
         (cl-mpm/setup::make-block-mps
          (mapcar #'+ (mapcar (lambda (x) 0) size) (list 0 0))
          block-size
          ;; (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
          (list 1 (* (second block-size) e-scale mp-scale))
          ;; (list h 10d0)
          density
          'cl-mpm/particle::particle-elastic
          :E E
          :nu 0d0)))
      (cl-mpm/setup::set-mass-filter sim density :proportion 1d-9)

      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))


      (setf (cl-mpm::sim-damping-factor sim) 0d0)
      (cl-mpm/setup::setup-bcs
       sim
       :top '(nil 0 nil)
       :bottom '(nil 0 nil))
      sim))

  (let ((v0 -0.1d0)
        (L 10d0)
        )

    (setf (cl-mpm::sim-gravity *sim*) 0d0)
    (cl-mpm:iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (setf (varef (cl-mpm/particle:mp-velocity mp) 1)
             (* v0 (sin (*
                         ;; 0.5
                         pi
                         (/ (varef (cl-mpm/particle:mp-position mp) 1)
                            L))))))))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  )
(defun setup-beam (&key (mps 2) (refine 1))
  (let* ((L 20d0)
         (size (list 40 40))
         (block-size (list L 5))
         (block-offset 20)
         (e-scale refine)
         (mp-scale mps)
         (sim (cl-mpm/setup::make-simple-sim
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               ;; :sim-type 'cl-mpm/aggregate:mpm-sim-agg-usf
               ;; :sim-type 'cl-mpm:mpm-sim-usf
               :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic
               :args-list
               (append
                (list
                 :enable-aggregate t
                 :vel-algo :FLIP
                 ))))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 1d2)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
         (offset (* h 0)))
    (declare (double-float h density))
    (defparameter *sim* sim)
    (progn
      (let* ((E 1d7))
        (cl-mpm::add-mps
         sim
         (cl-mpm/setup::make-block-mps
          (mapcar #'+ (mapcar (lambda (x) 0) size) (list 0 (+ offset block-offset) 0))
          block-size
          (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
          density
          'cl-mpm/particle::particle-elastic
          :E E
          :nu 0.30d0)))
      (cl-mpm/setup::set-mass-filter sim density :proportion 1d-15)
      ;; (setf (cl-mpm:sim-mass-filter *sim*) 0d0)

      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))

      (cl-mpm/setup::setup-bcs
       sim
       :left '(0 0 nil)
       :top '(nil nil nil)
       :bottom '(nil 0 nil))
      ;; (let ((v0 -0.1d0)
      ;;       )
      ;;   (cl-mpm:iterate-over-mps
      ;;    (cl-mpm:sim-mps *sim*)
      ;;    (lambda (mp)
      ;;      (setf (varef (cl-mpm/particle:mp-velocity mp) 1)
      ;;            (* v0 (sin (*
      ;;                        0.5
      ;;                        pi
      ;;                        (/ (varef (cl-mpm/particle:mp-position mp) 0)
      ;;                           L))))))))
      sim))

  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  )

(defun compute-energy (sim)
  ;; (cl-mpm::sim-stats-energy sim)
  ;; (compute-mesh-energy sim)
  (compute-mp-energy sim)
  )

(defun compute-mp-energy (sim)
  (cl-mpm::reduce-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (* 0.5d0
        (cl-mpm/particle::mp-mass mp)
        (cl-mpm/fastmaths::mag-squared (cl-mpm/particle::mp-velocity mp))))
   #'+))
(defun compute-mesh-energy (sim)
  (typecase
      sim
      (cl-mpm::mpm-sim
       (multiple-value-bind (energy oobf power) (cl-mpm/dynamic-relaxation::combi-stats sim)
         energy)
       )
    (cl-mpm/aggregate::mpm-sim-aggregated
     (multiple-value-bind (energy oobf power) (cl-mpm/dynamic-relaxation::combi-stats sim)
       energy)
     )
    )
  ;; (multiple-value-bind (energy oobf power) (cl-mpm/dynamic-relaxation::combi-stats-aggregated sim)
  ;;   energy)
  ;; (cl-mpm::reduce-over-mps
  ;;  (cl-mpm:sim-mps sim)
  ;;  (lambda (mp)
  ;;    (* 0.5d0
  ;;       (cl-mpm/particle::mp-mass mp)
  ;;       (cl-mpm/fastmaths::mag-squared (cl-mpm/particle::mp-velocity mp))))
  ;;  #'+)
  )

(defun compute-strain-energy (sim)
  (cl-mpm::reduce-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (* 0.5d0
        (cl-mpm/particle::mp-volume mp)
        (cl-mpm/fastmaths:dot
         ;; (cl-mpm/fastmaths:fast-scale-voigt (cl-mpm/particle::mp-stress mp)
         ;;                                    (* 1.0d0 (the double-float (cl-mpm/fastmaths:det-3x3 (cl-mpm/particle::mp-deformation-gradient mp))))
         ;;                                    )
         ;(cl-mpm/particle::mp-stress-kirchoff mp)
         (cl-mpm/particle::mp-stress mp)
         (cl-mpm/particle::mp-strain mp)
         ;; (cl-mpm/fastmaths:fast-.*
         ;;  ;; (cl-mpm/utils:voigt-from-list (list 1d0 1d0 1d0 2d0 2d0 2d0))
         ;;  (cl-mpm/particle::mp-strain mp))
         )))
   #'+))
(defun compute-gpe-diff (sim)
  (cl-mpm::reduce-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (* -1d0
        (cl-mpm::sim-gravity sim)
        (cl-mpm/particle::mp-mass mp)
        (varef (cl-mpm/particle::mp-displacement mp) 1)))
   #'+))

(defun save-timestep-preamble (output-dir)
  (with-open-file (stream (merge-pathnames output-dir "./timesteps.csv") :direction :output :if-exists :supersede)
    (format stream "steps,time,ke,se,gpe~%")))

(defun save-timestep (sim output-dir step)
  (unless (uiop:file-exists-p (merge-pathnames output-dir "timesteps.csv"))
    (save-timestep-preamble output-dir))
  (with-open-file (stream (merge-pathnames "timesteps.csv" output-dir) :direction :output :if-exists :append)
    (format stream "~D,~f,~f,~f,~f~%"
            step
            (cl-mpm::sim-time sim)
            ;; (cl-mpm::sim-stats-energy sim)
            (compute-energy sim)
            (compute-strain-energy sim)
            (compute-gpe-diff sim))))

(defun run (&key (output-dir "./output/")

              (dt-scale 1d0)
              (total-time 10d0)
              )
  (let* (;; (total-time 10d0)
         (target-time 1d0)
         ;(dt-scale 0.9d0)
         (substeps 0))
    (defparameter *run-sim* t)
    (defparameter *data-time* (list))
    (defparameter *data-ke* (list))
    (defparameter *data-se* (list))
    (defparameter *data-gpe* (list))

    (setf (cl-mpm::sim-dt-scale *sim*) dt-scale)
    (setf (cl-mpm:sim-dt *sim*)
          (* dt-scale (cl-mpm/setup::estimate-elastic-dt *sim*)))
    (setf substeps (round target-time (cl-mpm:sim-dt *sim*)))

    ;; (setf (cl-mpm::sim-damping-factor *sim*) 0d0)
    (uiop:ensure-all-directories-exist (list output-dir))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
    (cl-mpm/output::save-simulation-parameters (merge-pathnames "settings.json" output-dir)
                                               *sim*
                                               (list :dt target-time))
    (format t "Substeps ~D~%" substeps)
    (vgplot:close-all-plots)
    (let ((total-iter 0)
          (dt-accumulator 0d0))
      (time (loop for steps from 0 to (round total-time target-time)
                  while *run-sim*
                  do
                     (let ()
                       (format t "Step ~d substep est ~D~%" steps (floor target-time (cl-mpm::sim-dt *sim*)))
                       (cl-mpm/dynamic-relaxation::save-vtks *sim* output-dir steps)
                       (incf dt-accumulator target-time)
                       (time
                        (loop while (> dt-accumulator 0d0)
                              do
                                 (progn
                                   (save-timestep *sim* output-dir total-iter)
                                   (cl-mpm::update-sim *sim*)
                                   (push (cl-mpm::sim-time *sim*) *data-time*)
                                   (push (cl-mpm::sim-stats-energy *sim*) *data-ke*)
                                   ;; (push (compute-energy *sim*) *data-ke*)
                                   (push (compute-strain-energy *sim*) *data-se*)
                                   (push (compute-gpe-diff *sim*) *data-gpe*)
                                   (incf total-iter)
                                   (decf dt-accumulator (cl-mpm::sim-dt *sim*))
                                   ;; (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup:estimate-elastic-dt *sim* :dt-scale dt-scale))
                                   ;; (setf (cl-mpm:sim-dt *sim*) (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
                                   ;; (format t "Dt ~E~%" (cl-mpm:sim-dt *sim*))
                                   )
                          (when (> (cl-mpm::sim-time *sim*) 0.75d0)
                            (setf (cl-mpm::sim-gravity *sim*) 0d0))
                          ))
                       ;; (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       ;;   (format t "CFL dt estimate: ~f~%" dt-e)
                       ;;   (format t "CFL step count estimate: ~D~%" substeps-e)
                       ;;   (setf substeps substeps-e))
                       (plot-domain)
                       ;; (vgplot:plot
                       ;;  *data-time* *data-ke* "ke"
                       ;;  *data-time* *data-se* "se"
                       ;;  *data-time* *data-gpe* "gpe"
                       ;;  *data-time*
                       ;;  (mapcar #'+
                       ;;          *data-ke*
                       ;;          *data-se*
                       ;;          *data-gpe*
                       ;;          )
                       ;;  "error"
                       ;;  )
                       ;; (vgplot:title (format nil "Time:~F - KE ~E - OOBF ~E - Work ~E"  (cl-mpm::sim-time *sim*) energy oobf work))
                       (swank.live:update-swank)
                       ))))

    )
  )

(defun test ()
  (setup :mps 2)
  (run))


(defun test-range ()
  (let* ((classes (list
                   ;; 'cl-mpm/aggregate::mpm-sim-usf
                   ;; 'cl-mpm/aggregate::mpm-sim-agg-usf
                   ;; 'cl-mpm::mpm-sim-usl
                   ;; 'cl-mpm::mpm-sim-musl
                   'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic
                   ))
         (names (list
                 ;; "USF"
                 ;; "USF-AGG"
                 ;; "USL"
                 ;; "MUSL"
                 "IMPLICIT"
                 )))
    (loop for c in classes
          for n in names
          do
             (progn
               (let ((output-dir (format nil "output-~A/" n)))
                 (setup-beam :refine 0.5 :mps 3)
                 ;; (setup :refine 1 :mps 2)
                 ;; (setf (cl-mpm::sim-damping-factor *sim*) (* 0.01d0 (cl-mpm/setup:estimate-critical-damping *sim*)))
                 (change-class *sim* c)
                 (run :output-dir output-dir
                      ;; :dt-scale 0.5d0
                      :dt-scale 1d0
                      :total-time 100d0
                      )))))
  )

(defun test-dt-scale ()
  (let* ((dt-scales (list
                     ;; 1d0
                     0.5d0
                     0.1d0
                     0.01d0
                     ;; 2d0
                     ;; 10d0
                     )))
    (loop for d in dt-scales
          do
             (progn
               (let ((output-dir (format nil "output-~A/" d)))
                 (setup :mps 4)
                 ;; (setf (cl-mpm::sim-damping-factor *sim*) (* 0.01d0 (cl-mpm/setup:estimate-critical-damping *sim*)))
                 ;; (change-class *sim* c)
                 (change-class *sim* 'cl-mpm/aggregate:mpm-sim-agg-usf)
                 (run :dt-scale d :output-dir output-dir)))))
  )
(defun test-refine ()
  (let* ((dt-scales (list
                     0.5
                     1
                     2
                     )))
    (loop for d in dt-scales
          do
             (progn
               (let ((output-dir (format nil "output-~A/" d)))
                 (setup :refine d :mps 2)
                 (setf (cl-mpm::sim-damping-factor *sim*) (* 0d0 (cl-mpm/setup:estimate-critical-damping *sim*)))
                 ;; (change-class *sim* c)
                 (change-class *sim* 'cl-mpm/aggregate::mpm-sim-agg-usf)
                 (run :dt-scale 0.50d0
                      :total-time 100d0
                      :output-dir output-dir)))))
  )
