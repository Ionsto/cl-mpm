(defpackage :cl-mpm/examples/buoyancy-virtual-stress
  (:use :cl))
(in-package :cl-mpm/examples/buoyancy-virtual-stress)
(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)

(declaim (optimize (debug 3) (safety 3) (speed 2)))

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-elastic) dt fbar)
  (cl-mpm::update-stress-kirchoff-dynamic-relaxation mesh mp dt fbar)
  ;; (cl-mpm::update-stress-kirchoff mesh mp dt fbar)
  )

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-elastic) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  ;(cl-mpm::update-domain-polar-2d mesh mp dt)
  (cl-mpm::update-domain-deformation mesh mp dt)
  ;; (cl-mpm::update-domain-midpoint mesh mp dt)
  ;; (cl-mpm::scale-domain-size mesh mp)
  )

(defun plot (sim &optional (plot :deformed))
  ;; (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
  ;;        (ms-x (first ms))
  ;;        (ms-y (second ms))
  ;;        )
  ;;   (vgplot:axis (list 0 ms-x
  ;;                      0 ms-y))
  ;;   (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
  ;; (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim))))
  ;;   (vgplot:format-plot t "set ytics ~f" h)
  ;;   (vgplot:format-plot t "set xtics ~f" h))
  (cl-mpm/plotter::simple-plot
   sim
   :plot :deformed
   :colour-func (lambda (mp) (cl-mpm/utils:varef (cl-mpm/particle::mp-stress mp) 0))
   ;; :contact-bcs *penalty-bc*
   ))

(defparameter *original-configuration* (list))
(defparameter *original-size* (list))

(declaim (notinline setup-test-column))
(defun setup-test-column (size block-size &optional (e-scale 1) (mp-scale 1))
  (let ()
    (let* ((sim
             (cl-mpm/setup::make-simple-sim;
              (/ 1d0 e-scale)
              (mapcar (lambda (x) (* x e-scale)) size)
              :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
              ;; :sim-type 'cl-mpm/aggregate:mpm-sim-agg-usf
              ;; :sim-type 'cl-mpm:mpm-sim-usf
              :args-list (list
                          :enable-fbar nil
                          :enable-aggregate t
                          :mp-removal-size nil
                          :enable-split nil)))
           (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
           (h-x (/ h 1d0))
           (density 10d0))
      (progn
        (let ()
          (cl-mpm::add-mps
           sim
           (cl-mpm/setup::make-block-mps
            (list 0 0 0)
            block-size
            (mapcar (lambda (e mp-s) (round  (* e mp-s) h-x)) block-size
                    (list
                     mp-scale
                     mp-scale
                     1))
            density
            'cl-mpm/particle::particle-elastic
            :E 1d9
            :nu 0.49d0
            :gravity-axis (cl-mpm/utils:vector-from-list (list 0d0 0d0 0d0))
            )))
        (format t "MP count ~D~%" (length (cl-mpm:sim-mps sim)))
        (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
        (setf (cl-mpm::sim-gravity sim) -10d0)
        (setf (cl-mpm:sim-damping-factor sim)
              (* 0.1d0 (cl-mpm/setup::estimate-critical-damping sim)))
        (cl-mpm/setup::setup-bcs
         sim
         :top '(nil 0 nil)
         :bottom '(nil 0 nil))

        (defparameter *bc-pressure*
          (cl-mpm/buoyancy::make-bc-buoyancy-clip
           sim
           (second size)
           density
           (lambda (pos datum) t)))

        (cl-mpm:add-bcs-force-list
         sim
         *bc-pressure*)

        ;; (cl-mpm/setup::set-mass-filter sim density :proportion 1d-9)

        (defparameter *original-configuration*
          (loop for mp across (cl-mpm:sim-mps sim) collect (cl-mpm/utils:varef (cl-mpm/particle:mp-position mp) 1)))
        (defparameter *original-size*
          (loop for mp across (cl-mpm:sim-mps sim) collect (cl-mpm/utils:vector-copy (cl-mpm/particle::mp-domain-size mp))))

        sim))))
;Setup
(declaim (notinline setup))
(defun setup (&key (refine 0d0) (mps 2d0))
  ;; (defparameter *sim* (setup-test-column '(1 60) '(1 50) (/ 1 5) 2))
  (let* ((e (expt 2 refine))
         (L 1d0)
         (h (/ L e))
         (width
           (* 2 h)))
    (format t "H:~E~%" h)
    (defparameter
        *sim*
      (setup-test-column (list
                          (* 2 L)
                          L)
                         (list
                          L
                          L)
                         (/ 1d0 h)
                         mps
                         )))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defun compute-error (&optional sim)
  (let* ((mp-list (loop for mp across (cl-mpm:sim-mps *sim*) collect mp))
         ;(y-ref (loop for pos in *original-configuration* collect (float pos 1e0)))
         (y-ref (loop for mp in mp-list collect (cl-mpm/utils::varef (cl-mpm/particle::mp-position mp) 1)))
         (syy (loop for mp in mp-list collect (float (magicl:tref (cl-mpm/particle::mp-stress mp) 0 0) 1d0)))
         (vp-0-list (loop for size in *original-size* collect (float (* (cl-mpm/utils:varef size 0) (cl-mpm/utils:varef size 1)) 1e0)))
         (density 10d0)
         (gravity -10d0)
         (datum 1d0)
         (syy-ref (mapcar (lambda (x) (* 1d0 gravity density (- datum x))) y-ref))
         (err (/
               (* (magicl:norm
                   (magicl:.-
                    (magicl:from-list syy-ref (list (length syy-ref) 1))
                    (magicl:from-list syy (list (length syy) 1))))
                  (first vp-0-list))
               (abs (* (reduce #'+ vp-0-list) (* 1/2 (* datum datum) density gravity))))))
    ;; (pprint syy-ref)
    ;; (pprint syy)
    err))

(defun save-csv (output-file)
  (let* ((df (lisp-stat:make-df '(:error)
                                (list (make-array 1 :initial-element (compute-error *sim*))))))
    (lisp-stat:write-csv df output-file :add-first-row t))
  ;; (let* ((mp-list
  ;;          (loop for mp across (cl-mpm:sim-mps *sim*)
  ;;                collect mp))
  ;;        (y (loop for mp in mp-list collect (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)))
  ;;        (y-ref (loop for pos in *original-configuration*
  ;;                     collect (float pos 1e0)))
  ;;        (syy (loop for mp in mp-list collect (float (magicl:tref (cl-mpm/particle::mp-stress mp) 1 0) 1d0)))
  ;;        (rho 80d0)
  ;;        (E 1d5)
  ;;        (g 10d0)
  ;;        (vp-0-list (loop for size in *original-size*
  ;;                         collect (float (* (cl-mpm/utils:varef size 0) (cl-mpm/utils:varef size 1)) 1e0)))
  ;;        (pressure -1e4)
  ;;        (max-y 50)
  ;;        (syy-ref (mapcar (lambda (x) pressure) y-ref))
  ;;        (df (lisp-stat:make-df '(:y :syy :syy-ref :vp)
  ;;                               (list
  ;;                                (coerce y-ref '(vector single-float))
  ;;                                (coerce syy '(vector double-float))
  ;;                                (coerce syy-ref 'vector)
  ;;                                (coerce vp-0-list 'vector)))))
  ;;   (lisp-stat:write-csv df output-file :add-first-row t))
  )

(/ 0.664d0 (* 160 512))
(defparameter *run-sim* nil)
(defun run-conv ()
  ;; (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./analysis_scripts/column/data/")) do (uiop:delete-file-if-exists f))
  (setf *run-sim* t)
  (defparameter *data-refine* (list))
  (defparameter *data-error* (list))
  (loop for i in '(2 3 4 5 6 7 8 9)
        while *run-sim*
        do
           (let* (;(elements (expt 2 i))
                  (mps 2)
                  (refine i)
                  (elements (expt 2 refine))
                  (h (/ 1d0 elements))
                  )
             (let* ((e elements))
               (setup :refine i :mps mps)
               (format t "Running sim size ~a ~a ~%" refine elements)
               (format t "Sim dt: ~a ~%" (cl-mpm:sim-dt *sim*))
               (let* ((h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
                      (ms 1d0))
                 (setf (cl-mpm::sim-mass-scale *sim*) ms)
                 (cl-mpm/dynamic-relaxation::run-load-control
                  *sim*
                  :output-dir (merge-pathnames (format nil "./output-~A_~D/" i mps))
                  :load-steps 10
                  :substeps (* 10 refine)
                  :plotter #'plot
                  :damping (sqrt 2)
                  :save-vtk-dr nil
                  :save-vtk-loadstep t
                  :dt-scale 1d0
                  :criteria 1d-5
                  ;; :loading-function (lambda (f) (setf (nth 1 (cl-mpm/buoyancy::bc-pressure-pressures *bc-pressure*)) (* f -1d4)))
                  ))
               ;; (plot-sigma-yy)
               (push (compute-error *sim*) *data-error*)
               (push h *data-refine*))
             (vgplot:loglog (mapcar (lambda (x) (/ 1d0 x)) *data-refine*) *data-error*)
             (save-csv (merge-pathnames (format nil "./analysis_scripts/virtual_stress/buoyancy/data/data-~A_~D.csv" i mps)))
             )))

;; (defun run-)

(defun stop ()
  (setf (cl-mpm::sim-run-sim *sim*) nil)
  (setf *run-sim* nil))

(defun run (&key (plotter (lambda (sim) (plot sim))))
  (cl-mpm/output:save-vtk-mesh (asdf:system-relative-pathname "cl-mpm" "output/mesh.vtk") *sim*)
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
  (defparameter *run-sim* t)
    (vgplot:close-all-plots)
    (vgplot:figure)
  (sleep 1)
  (let* ((target-time 2d0)
         (dt (cl-mpm:sim-dt *sim*))
         (dt-scale 0.1d0)
         (substeps (floor target-time dt)))
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 1d0 (cl-mpm/setup:estimate-critical-damping *sim*)))
    (setf (cl-mpm:sim-dt *sim*) (* dt-scale (cl-mpm/setup:estimate-elastic-dt *sim*)))
    (format t "Substeps ~D~%" substeps)
    (time (loop for steps from 0 to 100
                while *run-sim*
                do
                   (progn
                     (format t "Step ~d ~%" steps)
                     (time
                      (dotimes (i substeps)
                        (cl-mpm::update-sim *sim*)
                        (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))
                     (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       (format t "CFL dt estimate: ~f~%" dt-e)
                       (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf substeps substeps-e))
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-cells (merge-pathnames (format nil "output/sim_cells_~5,'0d.vtk" *sim-step*)) *sim*)
                     (incf *sim-step*)
                     (funcall plotter *sim*)
                     (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" steps)))
                     (swank.live:update-swank)
                     (sleep .01)
                     ))))
  (vgplot:figure)
  (plot-sigma-yy)
  ;; (vgplot:title "Velocity over time")
  ;; (vgplot:plot *time* *velocity*)
  )
(defun plot-sigma-yy (&optional sim)
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (let* ((x-slice-pos (loop for mp across mps maximize (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
           (mp-list
             (loop for mp across mps
                     collect mp))
           (y (loop for mp in mp-list collect (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)))
           (y-ref (loop for pos in *original-configuration*
                        collect pos))

           (syy (loop for mp in mp-list collect (magicl:tref (cl-mpm/particle::mp-stress mp) 1 0)))
           (pressure -1d4)
           (max-y 50)
           (syy-ref (mapcar (lambda (x) pressure) y-ref))
           )
      (vgplot:plot syy y-ref "first;;with points pt 7"
                   syy-ref y-ref "Ref;;with points pt 7"
                   )
      (vgplot:legend)
      )))

;; (defun compute-error (&optional sim)
;;   (with-accessors ((mps cl-mpm:sim-mps))
;;       *sim*
;;     (let* ((x-slice-pos (loop for mp across mps maximize (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
;;            (mp-list
;;              (loop for mp across mps
;;                      collect mp))
;;            (y (loop for mp in mp-list collect (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)))
;;            (y-ref (loop for pos in *original-configuration*
;;                         collect pos))

;;            (vp-0-list (loop for size in *original-size*
;;                             collect (*  (cl-mpm/utils:varef size 0) (cl-mpm/utils:varef size 1))))
;;            (vl-0 (loop for vp-0 in vp-0-list sum vp-0))

;;            (syy (loop for mp in mp-list collect (magicl:tref (cl-mpm/particle::mp-stress mp) 1 0)))
;;            (pressure -1d4)
;;            (max-y 50)
;;            (syy-ref (mapcar (lambda (x) pressure) y-ref))
;;            )
;;       (loop for ref in syy-ref
;;             for val in syy
;;             for vp-0 in vp-0-list
;;             sum (/ (* (abs (- ref val)) vp-0) (* max-y vl-0))))))

(defun save-sigma-yy (&optional sim)
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (let* ((x-slice-pos (loop for mp across mps maximize (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
           (mp-list
             (loop for mp across mps
                   when (>= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (- x-slice-pos 0.001))
                     collect mp))
           (y (loop for mp in mp-list collect (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)))
           (syy (loop for mp in mp-list collect (magicl:tref (cl-mpm/particle::mp-stress mp) 1 0))))
      (with-open-file (stream (merge-pathnames "output/consolidation.csv") :direction :output :if-exists :supersede)
        (format stream "coord_y,sigma_yy~%")
        (loop for ymp in y
              for syymp in syy
              do (format stream "~f, ~f ~%" ymp syymp))))))

  ;; (defparameter *h*) 1d0
;; (defparameter *x* (loop for i from -1.5d0 to 1.5d0 by 0.1d0 collect i))
;; (defparameter *sf* (cl-mpm/shape-fuction:make-shape))
;; (vgplot:plot *x* (mapcar (lambda (x) (cl-mpm/shape-function::shape-bspline x *h*)) *x*))
;; (vgplot:plot *x* (mapcar (lambda (x) (cl-mpm/shape-function::shape-bspline-dsvp x *h*)) *x*))

(defun test-sf ()
  (let* ((sf (cl-mpm/shape-function:make-shape-function-bspline 2 1d0))
         (x (loop for i from -1.5d0 upto 1.5d0 by 0.01d0
                  collect i)))
    (with-accessors ((s cl-mpm/shape-function:svp)
                     (ds cl-mpm/shape-function:dsvp))
        sf
      (vgplot:figure)
      (vgplot:title "Shape function")
      (vgplot:plot x (mapcar (lambda (i) (funcall s i 0d0)) x))
      (vgplot:figure)
      (vgplot:title "Shape function grad")
      (vgplot:plot x (mapcar (lambda (i) (first (funcall ds i 0d0))) x))
      )))

(defun update-strain (f strain)
  (let ((df (magicl:.+ (magicl:eye 2) (cl-mpm::voight-to-matrix strain))))
    (magicl:@ df f)))

(defun test-bspline ()
  (vgplot:figure)
  (let ((x (loop for x from -3 upto 3 by 0.01 collect x))
        (nodes '(nil nil nil nil
                t
                t t t t))
        (knots ;'(0 0 0 1 2 3 4)
               ;; '(0 0 0 0 0.5d0 1.5d0 2.5d0 3.5d0) 
               '(-3.5d0 -2.5d0 -1.5d0 -0.5d0 0.5d0 1.5d0 2.5d0 3.5d0)
               )
        )
    (print (cl-mpm/shape-function::make-bspline-knots nodes 1d0))
    (print (cl-mpm/shape-function::make-bspline-knots '(t t t t t t t t t) 1d0))
    (vgplot:plot x (mapcar (lambda (x) (cl-mpm/shape-function::bspline knots x 2 2)) x) ""
                 x (mapcar (lambda (x) (cl-mpm/shape-function::bspline-dsvp knots x 3 2)) x) "")
    ;; (vgplot:plot x (mapcar (lambda (x) (cl-mpm/shape-function::nodal-bspline nodes x 0 1d0)) x) ""
    ;;              x (mapcar (lambda (x) (cl-mpm/shape-function::nodal-bspline-dsvp nodes x 1 1d0)) x)) ""
    ))



(defun test-dr ()
  (setup :mps 2 :refine 1)
  (cl-mpm/dynamic-relaxation::run-load-control
   *sim*
   :output-dir "./output/"
   :load-steps 10
   :plotter #'plot
   :damping 1d0
   :save-vtk-dr t
   :save-vtk-loadstep t
   :substeps 10
   :dt-scale 1d0
   :criteria 1d-5
   ;; :loading-function (lambda (f) (setf (nth 1 (cl-mpm/buoyancy::bc-pressure-pressures *bc-pressure*)) (* f -1d4)))
   ;:plotter #'plot
   ))





