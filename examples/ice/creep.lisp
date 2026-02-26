(defpackage :cl-mpm/examples/ice/creep
  (:use :cl
   :cl-mpm/example))
(in-package :cl-mpm/examples/ice/creep)


(declaim (notinline plot-domain))
(defun plot-domain (&key (trial t))
  (when *sim*
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :trial trial
     ;; :colour-func (lambda (mp) (sqrt (cl-mpm/constitutive::voigt-j2 (cl-mpm/utils:deviatoric-voigt (cl-mpm/particle:mp-stress mp)))))
     ;; :colour-func #'cl-mpm/particle::mp-strain-plastic-vm
     ;; :colour-func (lambda (mp) (/ (cl-mpm/particle:mp-mass mp)
     ;;                              (cl-mpm/particle:mp-volume mp)))
     ;; :colour-func #'cl-mpm/particle::mp-damage
     :colour-func (lambda (mp) (cl-mpm/utils::varef (cl-mpm/particle::mp-stress mp ) 5))
     ))
  )

(defun setup (&key (refine 1) (mps 2)
                   (multigrid-refine 4)
                (ice-height 400d0)
                (aspect 1)
                )
  (let* ((density 918d0)
         (mesh-resolution (/ 10d0 refine))
         (ice-length (* aspect ice-height))
         (domain-size (list (+ ice-length (* 3 ice-height)) (* ice-height 2)))
         (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
         (block-size (list ice-length ice-height))
         (E 1d9))
    (defparameter *ice-length* ice-length)
    (setf *sim* (cl-mpm/setup::make-simple-sim mesh-resolution element-count
                                               :sim-type
                                               'cl-mpm/dynamic-relaxation::mpm-sim-dr-multigrid
                                               :args-list
                                               (list
                                                :enable-fbar t
                                                :enable-aggregate t
                                                :split-factor (* 1d0 (sqrt 2) (/ 1d0 mps))
                                                :max-split-depth 8
                                                :enable-split t
                                                :refinement multigrid-refine
                                                )))
    (setf mesh-resolution (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0 0)
      block-size
      (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
      density
      'cl-mpm/particle::particle-ice-linear
      :E 1d9
      :nu 0.24d0
      :viscosity 1.3d13
      :enable-viscosity t))
    (cl-mpm/setup:setup-bcs
     *sim*
     :left '(0 nil 0))
    :bottom '(nil 0 0)
    (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
    ))

(defun plot-time-disp ()
  (vgplot:figure)
  (vgplot:plot (reverse *data-time*) (reverse *data-disp*)
               ""
               ;; :xlabel "Time (months)"
               ;; :ylabel "Front position (m)"
               ;; :title "Ice Creep Test"
               )
  )
(defun save-csv (output-dir filename)
  (with-open-file (stream (merge-pathnames filename output-dir) :direction :output :if-exists :supersede)
    (format stream "time,disp~%")
    (loop for time in (reverse *data-time*)
          for disp in (reverse *data-disp*)
          do
             (format stream "~E, ~E~%" (float time 0e0) (float disp 0e0)))))

(defun test ()
  (setup
   :refine 0.25
   :ice-height 125d0
   :aspect 4d0
   :multigrid-refine 3)

  (let* ((day (* 24d0 60d0 60d0))
         (month (* 30d0 day))
         (dt (* 10d0 day))
         (total-time (* 6 month)))

    (defparameter *data-time* (list 0d0))
    (defparameter *data-disp* (list ))
    (flet ((get-disp ()
             (cl-mpm::reduce-over-mps
              (cl-mpm:sim-mps *sim*)
              (lambda (mp)
                (-
                 (cl-mpm/utils:varef
                  (cl-mpm/fastmaths:fast-.+
                   (cl-mpm/particle::mp-position mp)
                   (cl-mpm/particle::mp-domain-size mp)) 0)
                 500d0))
              #'max)))
      (push (get-disp) *data-disp*)
      (cl-mpm/dynamic-relaxation::run-quasi-time
       *sim*
       :output-dir "./output/"
       :dt dt
       :total-time total-time
       ;; :steps 1000
       :dt-scale 1d0
       :conv-criteria 1d-3
       :substeps 50
       :min-adaptive-steps -0
       :max-adaptive-steps 9
       :save-vtk-loadstep t
       :save-vtk-dr t

       :plotter (lambda (sim) (plot-domain)
                  (vgplot:title (format nil "T - ~E months~%" (/ (cl-mpm::sim-time sim) month))))

       :post-conv-step (lambda (sim) (plot-domain))
       :post-load-step (lambda (sim)
                         (push (/ (cl-mpm::sim-time sim)
                                  month) *data-time*)
                         (push 
                          (get-disp)
                               *data-disp*)
                         (save-csv "./analysis_scripts/ice/creep/" "data.csv")
                         )))
    (plot-time-disp)
    (save-csv "./analysis_scripts/ice/creep/" "data.csv")
    )
  )


(defun est-visco (A n)
  (* 3/2 Kn)
  )
