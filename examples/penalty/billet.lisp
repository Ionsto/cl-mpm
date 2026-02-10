(defpackage :cl-mpm/examples/billet
  (:use :cl
   :cl-mpm/example
   :cl-mpm/utils)
  )
(in-package :cl-mpm/examples/billet)
(declaim (notinline plot-domain))
(defun plot-domain ()
  (when *sim*
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :trial t
     :colour-func (lambda (mp) (cl-mpm/particle::mp-index mp)))))

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-elastic) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  (cl-mpm::update-domain-polar-2d mesh mp dt)
  ;; (cl-mpm::update-domain-max-corner-2d mesh mp dt)
  ;; (cl-mpm::scale-domain-size mesh mp)
  )

(defun setup (&key (refine 1) (mps 3))
  (let* ((L 10d0)
         (d 1d0)
         (domain-width 30d0)
         (h (/ 1d0 refine))
         (height 10d0)
         (domain-height (+ height (* 2 h)))
         (density 7750d3)
         (E 1d6)
         (domain-size (list domain-width domain-height))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list height height)))

    (setf *sim* (cl-mpm/setup::make-simple-sim
                 h
                 element-count
                 :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
                 :args-list (list :enable-aggregate t
                                  :enable-fbar t
                                  :enable-split t
                                  :split-factor (* 1.2d0 (sqrt 2) (/ 1d0 mps))
                                  )))
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0 0)
      block-size
      (mapcar (lambda (e) (* (/ e h) mps)) block-size)
      density
      'cl-mpm/particle::particle-vm
      :E E
      :nu 0.3d0
      :rho 20d3
      ;; 'cl-mpm/particle::particle-elastic
      ;; :E E
      ;; :nu 0.2d0
      ;; :index 0
      :gravity-axis (cl-mpm/utils:vector-zeros)))
    (cl-mpm/setup::set-mass-filter *sim* 1d0 :proportion 1d-15)
    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 nil nil))

    (let ((friction 0.5d0))
      (defparameter *penalty*
        (cl-mpm/penalty::make-bc-penalty-distance-point
         *sim*
         (cl-mpm/utils:vector-from-list '(0d0 -1d0 0d0))
         (cl-mpm/utils:vector-from-list (list (* 0.5d0 domain-width)
                                              height
                                              0d0))
         (/ domain-width 2)
         (* E 1d2)
         friction
         0d0)))
    (defparameter *current-inc* 0d0)
    (cl-mpm:add-bcs-force-list
     *sim*
     *penalty*)
    (setf (cl-mpm:sim-dt *sim*)
          (* 0.5d0 (cl-mpm/setup:estimate-elastic-dt *sim*)))
    (setf *run-sim* t))
  (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
  )

(defun run-time ()
  (change-class *sim* 'cl-mpm/aggregate::mpm-sim-agg-usf)
  ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic)
  (setf (cl-mpm::sim-velocity-algorithm *sim*) :BLEND)
  (let ((step 0))
    (cl-mpm/dynamic-relaxation::run-time
     *sim*
     :output-dir (format nil "./output/")
     :plotter
     (lambda (sim)
       (plot-domain)
       (vgplot:title (format nil "Step ~D" step))
       (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step)) :terminal "png size 1920,1080")
       (incf step))
     :damping 1d-3
     :dt 10d0
     :total-time 1000d0
     ;; :dt-scale 1d0
     :dt-scale (/ 0.5d0 (sqrt 1d2))
     :post-conv-step
     (lambda (sim)
       (let ((disp-rate -0.01d0))
         (defparameter *penalty-controller*
           (cl-mpm/bc::make-bc-closure
            nil
            (lambda ()
              (cl-mpm/penalty::bc-increment-center
               *penalty*
               (cl-mpm/utils:vector-from-list (list 0d0 (* disp-rate (cl-mpm:sim-dt *sim*)) 0d0))))))
         (cl-mpm::add-bcs-force-list
          *sim*
          *penalty-controller*)
         )))))
(defun run (&key (output-dir (format nil "./output/")))
  (let* ((lstps 20)
         (total-disp -5d0)
         (delta (/ total-disp lstps))
         (step 0))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (cl-mpm/dynamic-relaxation::run-load-control
     *sim*
     :output-dir output-dir
     :plotter (lambda (sim) (plot-domain))
     :loading-function (lambda (i)
                         (cl-mpm/penalty::bc-increment-center
                          *penalty*
                          (cl-mpm/utils:vector-from-list (list 0d0 delta 0d0))))
     :post-conv-step (lambda (sim)
                       (plot-domain)
                       (vgplot:title (format nil "Step ~D" step))
                       (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step)) :terminal "png size 1920,1080")
                       (incf step))
     :load-steps lstps
     :enable-plastic nil
     :kinetic-damping nil
     :damping 1d0
     :substeps 10
     :criteria 1d-3
     :save-vtk-dr t
     :save-vtk-loadstep t
     :dt-scale 1d0)))

(defun test ()
  (setup :mps 2 :refine 1)
  (run)
  ;; (dolist (r (list 1 2 3))
  ;;   (dolist (mps (list 2 4))
  ;;     (setup :mps mps :refine r)
  ;;     (run :output-dir (format nil "./output-~D-~D/" r mps))))
  ;; (run-time)
  )


;; (pprint (cl-mpm/fastmaths:mag
;;          (magicl:@ (cl-mpm/particle::mp-true-domain *mp*)
;;                    (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0)))))
