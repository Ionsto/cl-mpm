(defpackage :cl-mpm/examples/column-virtual-stress
  (:use :cl))
(in-package :cl-mpm/examples/column-virtual-stress)
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

(defun get-disp (load-mps)
  (- (/
      (loop for mp in load-mps
            sum (+
                 (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)
                 (* 0.5d0 (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0))
                 )) (length load-mps))
     *initial-surface*))

(defun max-stress (mp)
  (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0))

(defun simple-plot-contact (sim &key (plot :point) (colour-func (lambda (mp) 0d0)) (contact-bcs nil))
  (declare (function colour-func))
  "A simple GIMP plot that display only the position and size of the MPs in a sim"
  (vgplot:format-plot t "set palette defined (0 'blue', 2 'red')")
  (vgplot:format-plot t "set ticslevel 0")
  (multiple-value-bind (x y lx ly c)
    (loop for mp across (cl-mpm:sim-mps sim)
          collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
          collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
          collect (magicl:tref (cl-mpm/particle::mp-domain-size mp) 0 0) into lx
          collect (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0) into ly
          collect (funcall colour-func mp) into c
          finally (return (values x y lx ly c)))
    (let* ((points (cl-mpm/penalty::bc-penalty-contact-points contact-bcs))
           (c-x (loop for p in points collect (magicl:tref p 0 0)))
           (c-y (loop for p in points collect (magicl:tref p 1 0)))
           (c-v (loop for p in points collect 1d0)))
      (vgplot:format-plot t "set cbrange [~f:~f]" (reduce #'min c) (+ 1d-5 (reduce #'max c)))
      (cond
        ((eq plot :point)
         (vgplot:plot x y c ";;with points pt 7 lc palette"))
        ((eq plot :deformed)
         (if (and contact-bcs
                  points)
             (vgplot:plot x y lx ly c ";;with ellipses lc palette"
                          c-x c-y ";;with points pt 7")
             (vgplot:plot x y lx ly c ";;with ellipses lc palette"))))))

  (setf (cl-mpm/penalty::bc-penalty-contact-points contact-bcs) nil)
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms)))
    (vgplot:format-plot t "set xrange [~f:~f]" 0d0 ms-x)
    (vgplot:format-plot t "set yrange [~f:~f]" 0d0 ms-y)
    (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
  (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim))))
    (vgplot:format-plot t "set ytics ~f" h)
    (vgplot:format-plot t "set xtics ~f" h))
  (vgplot:replot))


(defun plot (sim &optional (plot :deformed))
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:axis (list 0 ms-x
                       0 ms-y))
    (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
  (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim))))
    (vgplot:format-plot t "set ytics ~f" h)
    (vgplot:format-plot t "set xtics ~f" h))
  (cl-mpm/plotter::simple-plot
   sim
   :plot :deformed
   ;; :contact-bcs *penalty-bc*
   ))

(defparameter *original-configuration* (list))
(defparameter *original-size* (list))

(declaim (notinline setup-test-column))
(defun setup-test-column (size block-size &optional (e-scale 1) (mp-scale 1))
  (let ((nd (length block-size)))
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
           (h-y (/ h 1d0))
           (density 80)
           (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
           )
      (progn
        (let ((block-position (list (* h-x (+ 0 (/ 1d0 (* 2d0 mp-scale)) ));mp-scale->1
                                    (* h-y (+ 0 (/ 1d0 (* 2d0 mp-scale)) )))))
          (cl-mpm::add-mps
           sim
           (cl-mpm/setup::make-block-mps
            (list 0 0 0)
            block-size
            (mapcar (lambda (e mp-s)
                      (round  (* e mp-s) h-x)
                      ) block-size
                        (list
                         ;; mp-scale
                         ;1
                         mp-scale
                         mp-scale
                         1
                         ))
            density
            'cl-mpm/particle::particle-elastic
            :E 1d4
            :nu 0.0d0)))
        (format t "MP count ~D~%" (length (cl-mpm:sim-mps sim)))
        (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
        (setf (cl-mpm::sim-gravity sim) 0d0)
        (setf (cl-mpm:sim-damping-factor sim)
              (* 0.1d0 (cl-mpm/setup::estimate-critical-damping sim)))
        ;; (cl-mpm/setup::set-mass-filter sim density :proportion 0d-9)
        (setf (cl-mpm:sim-mass-filter sim) 0d0)
        (setf (cl-mpm:sim-dt sim) 1d-2)
        (cl-mpm/setup::setup-bcs
         sim
         :top '(nil nil nil)
         :bottom '(0 0 nil))
        (defparameter *bc-pressure*
          (cl-mpm/buoyancy::make-bc-pressure
           sim
           0d0
           -1d4))
        (cl-mpm:add-bcs-force-list
         sim
         *bc-pressure*)

        (defparameter *original-configuration*
          (loop for mp across (cl-mpm:sim-mps sim) collect (cl-mpm/utils:varef (cl-mpm/particle:mp-position mp) 1)))
        (defparameter *original-size*
          (loop for mp across (cl-mpm:sim-mps sim) collect (cl-mpm/utils:vector-copy (cl-mpm/particle::mp-domain-size mp))))

        sim))))
;; (setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))
(defun simple-time ()
  (setup)
  (time
   (dotimes (i 1000)
     (cl-mpm::update-sim *sim*))))
;Setup
(defun setup (&key (refine 0d0) (mps 2d0))
  ;; (defparameter *sim* (setup-test-column '(1 60) '(1 50) (/ 1 5) 2))
  (let* ((e (expt 2 (+ 4 refine)))
         (L 50d0)
         (h (/ L e))
         (width
           L
           ;; (* 2 h)
                )
         )
    (format t "H:~E~%" h)
    (defparameter
        *sim*
      (setup-test-column (list
                          width
                          (+ L (* 2 h))
                          ;; h
                          )
                         (list
                          width
                          L
                          ;; h
                          )
                         (/ 1d0 h)
                         mps
                         )))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defun compute-error (&optional sim)
  (let* ((mp-list (loop for mp across (cl-mpm:sim-mps *sim*) collect mp))
         (y-ref (loop for pos in *original-configuration* collect (float pos 1e0)))
         (syy (loop for mp in mp-list collect (float (magicl:tref (cl-mpm/particle::mp-stress mp) 1 0) 1d0)))
         (vp-0-list (loop for size in *original-size* collect (float (* (cl-mpm/utils:varef size 0) (cl-mpm/utils:varef size 1)) 1e0)))
         (pressure -1e4)
         (syy-ref (mapcar (lambda (x) pressure) y-ref))
         (err (/
               (* (magicl:norm
                   (magicl:.-
                    (magicl:from-list syy-ref (list (length syy-ref) 1))
                    (magicl:from-list syy (list (length syy) 1))))
                  (first vp-0-list))
                 (abs (* (reduce #'+ vp-0-list) pressure)))))
    err))

(defun save-csv (output-file)
  (let* ((df (lisp-stat:make-df '(:error)
                                (list (make-array 1 :initial-element (compute-error))))))
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

(defparameter *run-sim* nil)
(defun run-conv ()
  ;; (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./analysis_scripts/column/data/")) do (uiop:delete-file-if-exists f))
  (setf *run-sim* t)
  (defparameter *data-refine* (list))
  (defparameter *data-error* (list))
  (loop for i in '(2 3 4 5 6 7 8 9 10 12)
        while *run-sim*
        do
           (let* (;(elements (expt 2 i))
                  (refine i)
                  (elements (expt 2 refine))
                  (mps 2))
             (let* ((e elements)
                    (L 50d0)
                    (width L)
                    (h (/ L e)))
               (format t "H:~E~%" h)
               (defparameter
                   *sim*
                 (setup-test-column (list h (+ L (* 2 h)))
                                    (list h L)
                                    (/ 1d0 h)
                                    mps))
               (format t "Running sim size ~a ~a ~%" refine elements)
               (format t "Sim dt: ~a ~%" (cl-mpm:sim-dt *sim*))
               (let* ((h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
                      (ms 1d0))
                 (setf (cl-mpm::sim-mass-scale *sim*) ms)
                 (cl-mpm/dynamic-relaxation::run-load-control
                  *sim*
                  :output-dir (merge-pathnames (format nil "./output-~A_~D/" i mps))
                  :load-steps 10
                  :substeps (* 20 refine)
                  :plotter #'plot-sigma-yy
                  :damping (sqrt 2)
                  :save-vtk-dr t
                  :save-vtk-loadstep t
                  :dt-scale 1d0
                  :criteria 1d-9
                  :loading-function (lambda (f) (setf (nth 1 (cl-mpm/buoyancy::bc-pressure-pressures *bc-pressure*)) (* f -1d4)))
                  ))
               ;; (plot-sigma-yy)
               (push (compute-error *sim*) *data-error*)
               (push h *data-refine*))
             (vgplot:loglog (mapcar (lambda (x) (/ 1d0 x)) *data-refine*) *data-error*)
             (save-csv (merge-pathnames (format nil "./analysis_scripts/virtual_stress/column-virtual-stress/data/data-~A_~D.csv" i mps)))
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

(defun test ()
  (let ((a (make-array 2 :initial-contents '(0d0 0d0) :element-type 'double-float))
        (b (make-array 2 :initial-contents '(0d0 0d0) :element-type 'double-float))
        )
    ;(declare ((simple-array double-float 2) a b))
    (time
     (dotimes (i 1000000)
       (map '(vector double-float 2) #'+ a b)
       ))))


(defun cundall-test ()
  (let ((test-refines (list 0 1 2 3)))
    (defparameter *data-conv-energy* (list))
    (defparameter *data-conv-strain-energy* (list))
    (defparameter *data-conv-kinetic-energy* (list))
    (defparameter *data-conv-step* (list))
    (setf *run-sim* t)
    (loop for refine in test-refines
          do
             (let* ((step 0)
                    (work 0d0)
                    (dt-scale 0.80d0)
                    (time 0d0)
                    )
               (defparameter data-cundall-step (list))
               (defparameter data-cundall-load (list))
               (defparameter data-cundall-energy (list))
               (defparameter data-visc-step (list))
               (defparameter data-visc-load (list))
               (defparameter data-visc-energy (list))
               (defparameter data-visc-strain-energy (list))
               (defparameter data-visc-gravity-energy (list))
               (defparameter data-visc-kinetic-energy (list))
               (setup :refine refine :mps 2)
               (print (length (cl-mpm:sim-mps *sim*)))
               (setf (cl-mpm:sim-mass-scale *sim*) 1d0)
               (let* ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
                      (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh *sim*))))
                 (format t "H: ~E~%" h)
                 (format t "ND: ~D~%" nd)
                 (setf
                  (cl-mpm:sim-damping-factor *sim*)
                  (*
                   0.1d0
                   (cl-mpm/setup::estimate-critical-damping *sim*))))

               (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
               (format t "Critical dt ~E~%" (cl-mpm:sim-dt *sim*))
               (vgplot:close-all-plots)
               (vgplot:figure)
               (let*
                   ((top-y
                      (loop for mp across (cl-mpm:sim-mps *sim*)
                            maximize (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)))
                    (test-mp
                      (first
                       (loop for mp across (cl-mpm:sim-mps *sim*)
                             when (= top-y (magicl:tref (cl-mpm/particle::mp-position mp) 1 0))
                               collect mp)))
                    (total-disp -1d-6)
                    (load-steps 100)
                    (substeps (* 1 (expt 2 refine)))
                    (disp-inc (/ total-disp (* load-steps substeps)))
                    )
                 (setf *target-displacement* -1d0)

                 (format t "possible substeps?: ~E~%" (round 0.1d0 (cl-mpm:sim-dt *sim*)))
                 (time
                  (loop for step from 0 below load-steps
                        while *run-sim*
                        do
                           (progn
                             (format t "Step ~D~%" step)
                             (dotimes (i substeps
                                         ;(round 0.01d0 (cl-mpm:sim-dt *sim*))
                                       )
                               (setf cl-mpm/penalty::*debug-force* 0)
                               ;; (incf *target-displacement* disp-inc)
                               (cl-mpm:update-sim *sim*)
                               (incf work (cl-mpm/dynamic-relaxation:estimate-power-norm *sim*))
                               (incf time (cl-mpm:sim-dt *sim*))
                               )
                             (format t "Step ~D - Load ~E - OOBF ~E - KE ~E~%"
                                     step
                                     cl-mpm/penalty::*debug-force*
                                     (cl-mpm/dynamic-relaxation::estimate-oobf *sim*)
                                     (abs (/ (cl-mpm/dynamic-relaxation:estimate-energy-norm *sim*) work)))

                             (push time data-visc-step)
                             (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~2,'0d_~5,'0d.vtk" refine step)) *sim*)
                             (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~2,'0d_~5,'0d.vtk" refine step)) *sim*)
                             (push
                              ;; (/
                              ;;  (lparallel:pmap-reduce (lambda (mp)
                              ;;                           (* 0.5d0
                              ;;                              (cl-mpm/particle::mp-mass mp)
                              ;;                              (cl-mpm/fastmaths::mag-squared
                              ;;                               (cl-mpm/particle::mp-velocity mp))))
                              ;;                         #'+ (cl-mpm:sim-mps *sim*))
                              ;;  (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
                              ;; (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*)
                              (+
                               (magicl:tref (cl-mpm/particle::mp-position test-mp) 1 0)
                               (* 0.5d0
                                  (magicl:tref
                                   (cl-mpm/particle::mp-domain-size
                                    test-mp)
                                   1 0)))
                              data-visc-load)
                             (push
                              (lparallel:pmap-reduce (lambda (mp)
                                                       (* 0.5d0
                                                          (cl-mpm/particle::mp-mass mp)
                                                          (cl-mpm/fastmaths::mag-squared
                                                           (cl-mpm/particle::mp-velocity mp))))
                                                     #'+ (cl-mpm:sim-mps *sim*))
                              data-visc-kinetic-energy)
                             (push
                              (lparallel:pmap-reduce (lambda (mp)
                                                       (* 0.5d0
                                                          (cl-mpm/particle::mp-volume mp)
                                                          (cl-mpm/fastmaths:dot
                                                           (cl-mpm/particle::mp-stress mp)
                                                           (cl-mpm/particle::mp-strain mp))))
                                                     #'+ (cl-mpm:sim-mps *sim*))
                              data-visc-strain-energy)
                             (push
                              (lparallel:pmap-reduce (lambda (mp)
                                                       (*
                                                        -1d0
                                                        (cl-mpm/particle::mp-mass mp)
                                                        (cl-mpm/particle::mp-gravity mp)
                                                        (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)))
                                                     #'+ (cl-mpm:sim-mps *sim*))
                              data-visc-gravity-energy)
                             ;; (vgplot:plot
                             ;;  data-visc-step data-visc-load "Visc")
                             ;; (vgplot:plot
                             ;;  data-visc-step data-visc-kinetic-energy "KE"
                             ;;  data-visc-step data-visc-strain-energy "SE"
                             ;;  data-visc-step data-visc-gravity-energy "GPE"
                             ;;  data-visc-step (mapcar #'+
                             ;;                         data-visc-strain-energy
                             ;;                         data-visc-kinetic-energy
                             ;;                         data-visc-gravity-energy
                             ;;                         ) "Total"
                             ;;  )



                             ;; (apply #'vgplot:plot (reduce #'append (mapcar #'list *data-conv-step* *data-conv-energy* (list "1" "2" "3" "4"))))
                             ;; (plot *sim*)
                             (ecase (length *data-conv-step*)
                               (0
                                (vgplot:plot
                                 data-visc-step data-visc-load "0"))
                               (1
                                (vgplot:plot
                                 (nth 0 *data-conv-step*) (nth 0 *data-conv-energy*) "0"
                                 data-visc-step data-visc-load "1"))
                               (2
                                (vgplot:plot
                                 (nth 0 *data-conv-step*) (nth 0 *data-conv-energy*) "0"
                                 (nth 1 *data-conv-step*) (nth 1 *data-conv-energy*) "1"
                                 data-visc-step data-visc-load "2"))
                               (3
                                (vgplot:plot
                                 (nth 0 *data-conv-step*) (nth 0 *data-conv-energy*) "0"
                                 (nth 1 *data-conv-step*) (nth 1 *data-conv-energy*) "1"
                                 (nth 2 *data-conv-step*) (nth 2 *data-conv-energy*) "2"
                                 data-visc-step data-visc-load "3"))
                               (4
                                (vgplot:plot
                                 (nth 0 *data-conv-step*) (nth 0 *data-conv-energy*) "0"
                                 (nth 1 *data-conv-step*) (nth 1 *data-conv-energy*) "1"
                                 (nth 2 *data-conv-step*) (nth 2 *data-conv-energy*) "2"
                                 (nth 3 *data-conv-step*) (nth 3 *data-conv-energy*) "3"
                                 data-visc-step data-visc-load "4"))
                               (5
                                (vgplot:plot
                                 (nth 0 *data-conv-step*) (nth 0 *data-conv-energy*) "0"
                                 (nth 1 *data-conv-step*) (nth 1 *data-conv-energy*) "1"
                                 (nth 2 *data-conv-step*) (nth 2 *data-conv-energy*) "2"
                                 (nth 3 *data-conv-step*) (nth 3 *data-conv-energy*) "3"
                                 (nth 4 *data-conv-step*) (nth 4 *data-conv-energy*) "4"
                                 data-visc-step data-visc-load "5"))
                               (6
                                (vgplot:plot
                                 (nth 0 *data-conv-step*) (nth 0 *data-conv-energy*) "0"
                                 (nth 1 *data-conv-step*) (nth 1 *data-conv-energy*) "1"
                                 (nth 2 *data-conv-step*) (nth 2 *data-conv-energy*) "2"
                                 (nth 3 *data-conv-step*) (nth 3 *data-conv-energy*) "3"
                                 (nth 4 *data-conv-step*) (nth 4 *data-conv-energy*) "4"
                                 (nth 5 *data-conv-step*) (nth 5 *data-conv-energy*) "5"
                                 data-visc-step data-visc-load "6"))

                               (t
                                (vgplot:plot
                                 (nth 0 *data-conv-step*) (nth 0 *data-conv-energy*) "0"
                                 (nth 1 *data-conv-step*) (nth 1 *data-conv-energy*) "1"
                                 (nth 2 *data-conv-step*) (nth 2 *data-conv-energy*) "2"
                                 (nth 3 *data-conv-step*) (nth 3 *data-conv-energy*) "3"
                                 (nth 4 *data-conv-step*) (nth 4 *data-conv-energy*) "4"
                                 data-visc-step data-visc-load "5"))
                               )
                             ;)

                             (swank.live:update-swank)))))
               (push data-visc-load *data-conv-energy*)
               (push data-visc-step *data-conv-step*)
               )))
  ;; (vgplot:figure)
  ;; (vgplot:plot (nth 0 *data-conv-step*) (nth 0 *data-conv-energy*) "0"
  ;;              (nth 1 *data-conv-step*) (nth 1 *data-conv-energy*) "1"
  ;;              (nth 2 *data-conv-step*) (nth 2 *data-conv-energy*) "2"
  ;;              (nth 3 *data-conv-step*) (nth 3 *data-conv-energy*) "3"
  ;;              (nth 4 *data-conv-step*) (nth 4 *data-conv-energy*) "4"
  ;;              )
  )


;; (let ((h (* 0.3d0 (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
;;       (iter 1000000))
;;   (time (dotimes (i iter)
;;           (cl-mpm::iterate-over-neighbours-point-linear-3d (cl-mpm:sim-mesh *sim*) (cl-mpm/utils:vector-from-list (list h h 0d0)) (lambda (m n w g)))))
;;   (time (dotimes (i iter)
;;           (cl-mpm::iterate-over-neighbours-point-linear-2d (cl-mpm:sim-mesh *sim*) (cl-mpm/utils:vector-from-list (list h h 0d0)) (lambda (m n w g)))))
;;   (time (dotimes (i iter)
;;           (cl-mpm::iterate-over-neighbours-point-linear-simd (cl-mpm:sim-mesh *sim*) (cl-mpm/utils:vector-from-list (list h h 0d0)) (lambda (m n w g)))))
;;   )
(defun run-static (&key (output-dir "./output/")
                        (load-steps 10)
                     (dt-scale 0.5d0))

  (defparameter *run-sim* t)
  (uiop:ensure-all-directories-exist (list output-dir))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  (setf (cl-mpm:sim-damping-factor *sim*)
        (* 1d-1
           (cl-mpm/setup:estimate-critical-damping *sim*)))
  (with-open-file (stream (merge-pathnames output-dir "conv.csv") :direction :output :if-exists :supersede)
    (format stream "iter,step,damage,plastic,oobf,energy~%"))
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (let* ((load-0 0d0)
           ;; (load -5d0)
           (load (cl-mpm/particle::mp-gravity (aref (cl-mpm:sim-mps *sim*) 0)))
           (load-inc (/ (- load load-0) load-steps))
           (current-load load-0)
           (total-iter 0)
           )
      (loop for step from 0 below load-steps
            while *run-sim*
            do
               (incf current-load load-inc)
               (cl-mpm:iterate-over-mps
                mps
                (lambda (mp)
                  (setf (cl-mpm/particle:mp-gravity mp) current-load)))
               ;; (cl-mpm/penalty::bc-increment-center *bc-squish* (cl-mpm/utils:vector-from-list (list 0d0 load-inc 0d0)))
               (let ((crit (cl-mpm/setup:estimate-critical-damping *sim*)))
                 (setf (cl-mpm:sim-damping-factor *sim*)
                       (* 0.1d0 crit))
                 (defparameter *ke-last* 0d0)
                 (let ((conv-steps 0)
                       (substeps 1))
                   (time
                    (cl-mpm/dynamic-relaxation:converge-quasi-static
                     *sim*
                     :oobf-crit 1d-2
                     :energy-crit 1d-2
                     :kinetic-damping t
                     :damping-factor 1d-2
                     :dt-scale dt-scale
                     :substeps substeps
                     :conv-steps 1000
                     :post-iter-step
                     (lambda (i energy oobf)
                       (setf conv-steps (* substeps i))
                       (plot *sim*)
                       (vgplot:title (format nil "Step ~D - substep ~D - KE ~E - OOBF ~E"  step i energy oobf))
                       (format t "Substep ~D~%" i)
                       (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d_~5,'0d.vtk" step i)) *sim*)
                       (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d_~5,'0d.vtk" step i)) *sim*)
                       (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d_~5,'0d.vtk" step i)) *sim*)
                       (with-open-file (stream (merge-pathnames output-dir "conv.csv") :direction :output :if-exists :append)
                         (format stream "~D,~D,~f,~f,~f,~f~%" total-iter step (cl-mpm/dynamic-relaxation::get-plastic *sim*) (cl-mpm/dynamic-relaxation::get-damage *sim*)
                                 oobf energy))
                       (incf total-iter)
                       )))
                   (plot *sim*)
                   (vgplot:title (format nil "Step ~D - ~D" step conv-steps))
                   ))
               ;; (cl-mpm/penalty::bc-increment-center *bc-squish* (cl-mpm/utils:vector-from-list (list 0d0 load-inc 0d0)))
               (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~5,'0d.vtk" step)) *sim*)
               (cl-mpm/output::save-vtk-cells (merge-pathnames output-dir (format nil "sim_cells_~5,'0d.vtk" step)) *sim*)
               (cl-mpm::finalise-loadstep *sim*)
               (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step))
                                  :terminal "png size 1920,1080"
                                  )
               (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~5,'0d.vtk" step)) *sim*)
               (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim_p_~5,'0d.vtk" step)) *sim* )
               (sleep 0.1d0)
               (swank.live:update-swank)

            ))))



(defun save-test-vtks (&key (output-dir "./output/"))
  (cl-mpm/output:save-vtk (merge-pathnames "test.vtk" output-dir) *sim*)
  (cl-mpm/output:save-vtk-nodes (merge-pathnames "test_nodes_0.vtk" output-dir) *sim*)
  (cl-mpm/output:save-vtk-nodes (merge-pathnames "test_nodes_1.vtk" output-dir) *sim*)
  (cl-mpm/output:save-vtk-cells (merge-pathnames "test_cells.vtk" output-dir) *sim*)
  )


(defun test-dr ()
  (setup :mps 2 :refine 1)
  (cl-mpm/dynamic-relaxation::run-load-control
   *sim*
   :output-dir "./output/"
   :load-steps 10
   :plotter #'plot-sigma-yy
   :damping 1d0
   :kinetic-damping nil
   :adaptive-damping t
   :save-vtk-dr t
   :save-vtk-loadstep t
   :substeps 50
   :dt-scale 0.5d0
   :criteria 1d-5
   :loading-function (lambda (f) (setf (nth 1 (cl-mpm/buoyancy::bc-pressure-pressures *bc-pressure*)) (* f -1d4)))
   ;:plotter #'plot
   ))


(defmacro time-form (it form)
  `(progn
     (declaim (optimize speed))
     (let* ((iterations ,it)
            (start (get-internal-real-time)))
       (time
        (dotimes (i ,it)
          ,form))
       (let* ((end (get-internal-real-time))
              (units internal-time-units-per-second)
              (dt (/ (- end start) (* iterations units)))
              )
         (format t "Total time: ~f ~%" (/ (- end start) units)) (format t "Time per iteration: ~f~%" (/ (- end start) (* iterations units)))
         (format t "Throughput: ~f~%" (/ 1 dt))
         (format t "Time per MP: ~E~%" (/ dt (length (cl-mpm:sim-mps *sim*))))
         dt))))


(defun test-hotloop ()
  (format t "Setup~%")
  (time (setup :refine 10))
  (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (time (cl-mpm:update-sim *sim*))
  (setf (cl-mpm:sim-dt *sim*) (* 0.25d0 (cl-mpm::calculate-min-dt *sim*)))
  (let ((damping-factor 1d0))
    (time-form
     1000
     (progn
       (cl-mpm:update-sim *sim*)
       ;; (setf (cl-mpm:sim-damping-factor *sim*) (* damping-factor (cl-mpm/dynamic-relaxation::dr-estimate-damping *sim*)))
       ))))


(defun test-inv ()
  (let ((voigt (cl-mpm/utils:voigt-from-list (list 0d0 0d0 0d0 1d0 0d0 0d0)))
        (mat (cl-mpm/utils:matrix-from-list (list 1d0 2d0 3d0 4d0 5d0 6d0 7d0 8d0 0d0))))
    ;; (pprint (magicl:inv mat))
    ;; (pprint (cl-mpm/fastmaths::fast-inv-3x3 mat))
    (pprint (magicl:@ (cl-mpm/utils:voight-to-matrix voigt) (magicl:inv mat)))
    ;; (pprint (cl-mpm/utils::matrix-to-voight (magicl:linear-solve mat (cl-mpm/utils:voight-to-matrix voigt))))
    (pprint (cl-mpm/fastmaths::linear-solve-3x3-voigt mat voigt))
    ;; (time-form 1000000 (cl-mpm/fastmaths::fast-inv-3x3 mat))
    ;; (time-form 1000000 (magicl::inv mat))
    ))

(defun test-refine ()
  (loop for r in (list 2 3 4 5 6 7 8 9 10)
        do (progn
             (let* ((refine r)
                    (elements (expt 2 refine))
                    (mps 2)
                    (final-time 15))
             (let* ((e elements)
                    (L 50d0)
                    (h (/ L e)))
               (defparameter
                   *sim*
                 (setup-test-column (list h (+ L h))
                                    (list h L)
                                    (/ 1d0 h)
                                    mps))))
             (cl-mpm::iterate-over-mps
              (cl-mpm:sim-mps *sim*)
              (lambda (mp)
                (setf (cl-mpm/particle::mp-gravity mp) -1d0)))
             ;; (setf (cl-mpm::sim-mass-scale *sim*) 1d-2)
             (cl-mpm/dynamic-relaxation::run-load-control
              *sim*
              :output-dir (format nil "./output-~D_~F/" r (cl-mpm::sim-mass-scale *sim*))
              :plotter #'plot
              :load-steps 1
              :damping (* (sqrt (cl-mpm::sim-mass-scale *sim*)) 1d0)
              :substeps 50
              :criteria 1d-5
              :adaptive-damping t
              :kinetic-damping nil
              :save-vtk-dr t
              :dt-scale (/ 0.4d0 (sqrt 1d0))
              ))))



(defun test-inter ()
  (let ((data-x (list 0d0))
        (data-v (list 0d0))
        (time (list 0d0))
        (m 1d0)
        (x 0d0)
        (v 1d0)
        (stiffness 1d0)
        (final-time 20d0)
        (dt 0.01d0)
        (damping 1.1d0)
        (ct 0d0)
        )
    (setf *run-sim* t)
    (loop for i from 0 to (round final-time dt)
          while *run-sim*
          do (progn
               (dotimes (i 10)
                 (let* ((fi (* x stiffness -1))
                        (fd (* -1d0 v damping m))
                        (R (+ fi fd)))
                   ;; (setf v (+ v (* dt (/ R m))))
                   (let ((damping-f (exp (* -1d0 damping dt))))
                     (setf v (+ (* v damping-f) (* (/ fi m) (/ 1d0 damping) (- 1d0 damping-f)))))
                   (setf x (+ x (* dt v)))
                   (incf ct dt)
                   ))
               (push x data-x)
               (push v data-v)
               (push ct time)
               (vgplot:plot time data-x)
               (vgplot:format-plot t "set xrange [~f:~f]" 0d0 final-time)
               (swank.live:update-swank)))))
