(defpackage :cl-mpm/examples/column
  (:use :cl)
  (:import-from
    :cl-mpm/utils get-vector get-stress))
(in-package :cl-mpm/examples/column)
(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)
;; (sb-ext:restrict-compiler-policy 'speed  3 3)
;; (sb-ext:restrict-compiler-policy 'debug  0 0)
;; (sb-ext:restrict-compiler-policy 'safety 0 0)
;; (setf *block-compile-default* t)

(declaim (optimize (debug 3) (safety 3) (speed 2)))

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-elastic) dt fbar)
  (cl-mpm::update-stress-kirchoff mesh mp dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-dynamic-relaxation-incremental mesh mp dt fbar)
  )

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-elastic) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  ;(cl-mpm::update-domain-polar-2d mesh mp dt)
  (cl-mpm::update-domain-deformation mesh mp dt)
  ;; (cl-mpm::update-domain-midpoint mesh mp dt)
  ;; (cl-mpm::scale-domain-size mesh mp)
  )



(defun max-stress (mp)
  (get-stress (cl-mpm/particle:mp-stress mp) :xx))

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

;(defparameter *rho* 8000d0)
(defparameter *rho* 80d0)
(declaim (notinline setup-test-column))
(defun setup-test-column (size block-size &optional (e-scale 1) (mp-scale 1))
  (let ((nd (length block-size)))
    (let* ((sim
             (cl-mpm/setup::make-simple-sim
              (/ 1d0 e-scale)
              (mapcar (lambda (x) (* x e-scale)) size)
              :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
              ;; :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-usf
              :args-list (list
                          :enable-fbar nil
                          :enable-aggregate t
                          :mp-removal-size nil
                          :enable-split nil)))
           (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
           (density *rho*)
           (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
           )
      (progn
        (cl-mpm::add-mps
         sim
         (cl-mpm/setup::make-block-mps
          (list 0 0 0)
          block-size
          (mapcar (lambda (e mp-s) (round (* e mp-s) h)) block-size
                  (list mp-scale 1 1))
          density
          'cl-mpm/particle::particle-elastic
          :E 10d3
          :nu 0.0d0
          :gravity-axis (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0))
          ))
        (setf (cl-mpm:sim-gravity sim) -10d0)
        (format t "MP count ~D~%" (length (cl-mpm:sim-mps sim)))
        (setf (cl-mpm:sim-damping-factor sim)
              (* 0.1d0 (cl-mpm/setup::estimate-critical-damping sim)))
        (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
        (cl-mpm/setup::set-mass-filter sim density :proportion 1d-15)
        (setf (cl-mpm:sim-dt sim) (cl-mpm/setup:estimate-elastic-dt sim :dt-scale 0.5d0))
        (cl-mpm/setup::setup-bcs
         sim
         :left '(0 nil nil)
         :right '(nil nil nil)
         ;; :top '(nil nil nil)
         ;; :bottom '(0 0 nil)
         )
        (defparameter *original-configuration*
          (loop for mp across (cl-mpm:sim-mps sim) collect (cl-mpm/utils:get-vector (cl-mpm/particle:mp-position mp) :x)))
        (defparameter *original-size*
          (loop for mp across (cl-mpm:sim-mps sim) collect (cl-mpm/particle::mp-volume-0 mp)))
        sim))))

(defun simple-time ()
  (setup)
  (time
   (dotimes (i 1000)
     (cl-mpm::update-sim *sim*))))
;Setup
(defun setup (&key (refine 0d0) (mps 2d0))
  (let* ((e (expt 2 (+ 4 refine)))
         (L 50d0)
         (h (/ L e)))
    (format t "H:~E~%" h)
    (defparameter
        *sim*
      (setup-test-column (list
                          ;; h
                          (+ L (* 2 h))
                          ;; h
                          )
                         (list
                          ;; h
                          L
                          ;; h
                          )
                         (/ 1d0 h)
                         mps
                         )))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  (defparameter *sim-step* 0)
  )


(defun save-csv (output-file)
  (let* ((mp-list
           (loop for mp across (cl-mpm:sim-mps *sim*)
                 collect mp))
         (y (loop for mp in mp-list collect (cl-mpm/utils:get-vector (cl-mpm/particle::mp-position mp) :x)))
         (y-ref (loop for pos in *original-configuration*
                      collect (float pos 1e0)))

         (syy (loop for mp in mp-list collect (float (get-stress (cl-mpm/particle::mp-stress mp) :xx) 1e0)))
         (rho *rho*)
         (E 1d5)
         (g 10d0)
         (vp-0-list (loop for size in *original-size*
                          collect (float size 1e0)))
         (max-y 50)
         (syy-ref (mapcar (lambda (x) (float (* rho g (- x max-y)) 1e0)) y-ref))
         (df (lisp-stat:make-df '(:y :syy :syy-ref :vp)
                                (list
                                 (coerce y-ref '(vector single-float))
                                 (coerce syy 'vector)
                                 (coerce syy-ref 'vector)
                                 (coerce vp-0-list 'vector)))))
    (lisp-stat:write-csv df output-file :add-first-row t)))

(defparameter *run-sim* nil)
(defun run-conv ()
  (setf *run-sim* t)
  (defparameter *data-refine* (list))
  (defparameter *data-error* (list))
  (loop for i in
        '(4 5 6 7 8 9)
        while *run-sim*
        do
           (let* (;(elements (expt 2 i))
                  (refine i)
                  (elements (expt 2 refine))
                  (mps 2)
                  (final-time 15)
                  (step (list))
                  (res (list))
                  (total-step 0)
                  )
             (let* ((e elements)
                    (L 50d0)
                    (h (/ L e)))
               (format t "H:~E~%" h)
               (defparameter
                   *sim*
                 (setup-test-column (list (+ L (* 2 h))
                                          h)
                                    (list L
                                          h)
                                    (/ 1d0 h)
                                    mps))
               ;; (setf (cl-mpm::sim-gravity *sim*) (* -1 (cl-mpm::sim-gravity *sim*)))
               (format t "Running sim size ~a ~a ~%" refine elements)
               (format t "Sim dt: ~a ~%" (cl-mpm:sim-dt *sim*))
               (format t "Sim steps: ~a ~%" (/ final-time (cl-mpm:sim-dt *sim*)))
               (cl-mpm/dynamic-relaxation::run-load-control
                *sim*
                :output-dir (merge-pathnames (format nil "./output-~A_~D/" i mps))
                :load-steps 10
                :substeps (* 10 refine)
                :plotter (lambda (sim)
                           (vgplot:semilogy (reverse step) (reverse res))
                           ); #'plot-sigma-yy
                :post-iter-step (lambda (i o e)
                                  (push total-step step)
                                  (incf total-step)
                                  (push e res))
                :damping 1d0;(sqrt 2)
                :save-vtk-dr nil
                :save-vtk-loadstep nil
                :dt-scale 1d0
                :criteria 1d-5)
               ;; (plot-sigma-yy)
               (push (compute-error *sim*) *data-error*)
               (push h *data-refine*))
             (vgplot:loglog (mapcar (lambda (x) (/ 1d0 x)) *data-refine*) *data-error*)
             (save-csv (merge-pathnames (format nil "./analysis_scripts/column/data/data-~A_~D.csv" i mps)))
             )))

;; (defun run-)

(defun stop ()
  (setf (cl-mpm::sim-run-sim *sim*) nil)
  (setf *run-sim* nil))

(defun plot-sigma-yy (&optional sim)
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (let* ((mp-list
             (loop for mp across mps
                     collect mp))
           (y-ref (loop for pos in *original-configuration*
                        collect pos))
           (syy (loop for mp in mp-list collect (get-stress (cl-mpm/particle::mp-stress mp) :xx)))
           (rho *rho*)
           (E 1d5)
           (g 10d0)
           (max-y 50)
           (syy-ref (mapcar (lambda (x) (* rho g (- x max-y))) y-ref))
           )
      (vgplot:plot syy y-ref "first;;with points pt 7"
                   syy-ref y-ref "Ref;;with points pt 7"
                   )
      (vgplot:legend)
      )))

(defun compute-error (&optional sim)
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (let* ((mp-list
             (loop for mp across mps
                     collect mp))
           (y-ref (loop for pos in *original-configuration*
                        collect pos))
           (vp-0-list (loop for size in *original-size*
                            collect size))
           (vl-0 (loop for vp-0 in vp-0-list sum vp-0))

           (syy (loop for mp in mp-list collect (get-stress (cl-mpm/particle::mp-stress mp) :xx)))
           (rho *rho*)
           (E 1d5)
           (g 10d0)
           (max-y 50)
           (syy-ref (mapcar (lambda (x) (* rho g (- x max-y))) y-ref))
           )
      (loop for ref in syy-ref
            for val in syy
            for vp-0 in vp-0-list
            sum (/ (* (abs (- ref val)) vp-0) (* g rho max-y) vl-0)))))



(defun test-dr ()
  (setup :mps 2 :refine 1)
  ;; (change-class *sim* 'cl-mpm/aggregate::mpm-sim-agg-usf)
  ;; (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t)
  ;; (setf (cl-mpm::sim-velocity-algorithm *sim*) :FLIP)
  ;; (dotimes (i 100)
  ;;   (dotimes (j 10)
  ;;     (cl-mpm:update-sim *sim*))
  ;;   (cl-mpm/dynamic-relaxation::save-vtks *sim* "./output/" i))
  (cl-mpm/dynamic-relaxation::run-load-control
   *sim*
   :output-dir "./output/"
   :load-steps 40
   :plotter #'plot-sigma-yy
   :damping (sqrt 2)
   :save-vtk-dr t
   :save-vtk-loadstep t
   :substeps 10
   :dt-scale 1d0
   :criteria 1d-5
   ;:plotter #'plot
   )
  )


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


(defun test-refine ()
  (loop for r in (list 3 4 5 6 7 8 9 10)
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
             (cl-mpm/dynamic-relaxation::run-load-control
              *sim*
              :output-dir (format nil "./output-~D_~F/" r 2) 
              :plotter #'plot
              :load-steps 20
              :damping 1d0
              :substeps 50
              :criteria 1d-9
              :adaptive-damping t
              :kinetic-damping nil
              :save-vtk-dr nil
              :dt-scale 1d0
              ))))
