(defpackage :cl-mpm/examples/fbar/rigid-footing
  (:use :cl
   :cl-mpm/example
   :cl-mpm/utils)
  )
(in-package :cl-mpm/examples/fbar/rigid-footing)
(declaim (notinline plot-domain))
(defun plot-domain ()
  (when *sim*
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :trial t
     :colour-func (lambda (mp) (cl-mpm/particle::mp-index mp)))))

(defun get-load ()
  (* 2d0 (cl-mpm/penalty::resolve-load *penalty*)))

(declaim (notinline setup))
(defun setup (&key
                (refine 1)
                (mps 3)
                (enable-fbar t)
                (multigrid-refines 0)
                )
  (let* ((refine (* refine (expt 2 (- multigrid-refines))))
         (domain-width 5d0)
         (height domain-width)
         (width 0.5d0)
         (h (/ 0.5d0 refine))
         (domain-height (+ domain-width (* 2 h)))
         (density 1d3)
         (E 10d9)
         (nu 0.48d0)
         (domain-size (list domain-width domain-height))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list height height)))
    (setf
     *sim*
     (cl-mpm/setup::make-simple-sim
      h
      element-count
      :sim-type
       'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
      ;; 'cl-mpm/dynamic-relaxation::mpm-sim-dr-multigrid
      :args-list
      (list
       :enable-aggregate t
       :enable-split nil
       :enable-fbar enable-fbar
       :refinement multigrid-refines
       )))
    (setf h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0 0)
      block-size
      (mapcar (lambda (e) (* (/ e h) mps)) block-size)
      density
      'cl-mpm/particle::particle-vm
      :E E
      :nu nu
      :rho (* 1d6 (sqrt 3/2))
      ;; :rho 1d6
      ;; :rho 1d6
      ;; 'cl-mpm/particle::particle-mc
      ;; :E E
      ;; :nu 0.48d0
      ;; :phi 0d0
      ;; :psi 0d0
      ;; :c 1d6
      )
     )

    (setf (cl-mpm::sim-gravity *sim*) 0d0)
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)
    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 nil nil)
     :right '(0 nil nil)
     :bottom '(0 0 nil))

    (let* ((friction 0d0)
           (epsilon-scale 1d0)
           (epsilon (* (cl-mpm/particle::calculate-p-wave-modulus E nu) epsilon-scale))
           (width 0.5d0))
      (format t "Penalty parameter ~E~%" epsilon)
      (defparameter *penalty-down*
        (cl-mpm/penalty::make-bc-penalty-distance-point
         *sim*
         (cl-mpm/utils:vector-from-list '(0d0 -1d0 0d0))
         (cl-mpm/utils:vector-from-list (list (* 0.5d0 width)
                                              height
                                              0d0))
         (/ width 2)
         epsilon
         friction
         0d0))
      (defparameter *penalty-right*
        (cl-mpm/penalty::make-bc-penalty-distance-point
         *sim*
         (cl-mpm/utils:vector-from-list '(1d0 0d0 0d0))
         (cl-mpm/utils:vector-from-list (list
                                         width
                                         (+ height (/ width 2))
                                         0d0))
         (/ width 2)
         epsilon
         friction
         0d0))
      (defparameter *penalty*
        (cl-mpm/penalty::make-bc-penalty-structure
         *sim*
         epsilon
         friction
         0d0
         (list *penalty-down* *penalty-right*))))
    (cl-mpm:add-bcs-force-list
     *sim*
     *penalty*)
    )
  (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
  )

(defun plot-load-disp ()
  (vgplot:plot
   (mapcar (lambda (x) (* x -1d3)) *data-disp*)
   (mapcar (lambda (x) x) *data-load*)
   "Solution"
   )
  (vgplot:xlabel "Displacement (mm)")
  )

(defun save-csv (output-dir filename data-disp data-load)
  (with-open-file (stream (merge-pathnames filename output-dir) :direction :output :if-exists :supersede)
    (format stream "disp,load~%")
    (loop for disp in data-disp
          for load in data-load
          do (format stream "~E,~E~%" (float disp 0e0) (float load 0e0)))))
(declaim (notinline run))
(defun run (&key (output-dir (format nil "./output/"))
              (csv-dir (format nil "./output/"))
              (csv-filename (format nil "load-disp.csv")))
  (let* ((lstps 50)
         (total-disp -2d-3)
         (current-disp 0d0)
         (step 0))
    (defparameter *data-disp* (list 0d0))
    (defparameter *data-load* (list 0d0))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))

    (vgplot:close-all-plots)
    (time
     (cl-mpm/dynamic-relaxation::run-load-control
      *sim*
      :output-dir output-dir
      :plotter (lambda (sim) (plot-load-disp))
      :loading-function (lambda (i)
                          (setf current-disp (* i total-disp))
                          (cl-mpm/penalty::bc-set-displacement
                           *penalty*
                           (cl-mpm/utils:vector-from-list (list 0d0 current-disp 0d0))))
      :post-conv-step (lambda (sim)
                        (push current-disp *data-disp*)
                        (let ((load (get-load)))
                          (format t "Load ~E~%" load)
                          (push load *data-load*))
                        (plot-load-disp)
                        (save-csv csv-dir csv-filename *data-disp* *data-load*)
                        (incf step))
      :load-steps lstps
      :enable-plastic t
      :damping 1d0;(sqrt 2d0)
      :substeps 50
      :criteria 1d-6
      :save-vtk-dr nil
      :save-vtk-loadstep t
      :dt-scale 1d0))))

(defun test ()
  (setup :mps 3 :refine 1 :enable-fbar t :multigrid-refines 0)
  (run)
  ;; (save-csv "./examples/fbar/rigid-footing/" (format nil "data_fbaradjust_~A.csv" t) *data-disp* *data-load*)
  ;; (dolist (fbar (list t nil))
  ;;   (save-csv "./examples/fbar/rigid-footing/" (format nil "data_fbar_~A.csv" fbar) *data-disp* *data-load*)
  ;;   )
  )
