(defpackage :cl-mpm/examples/column
  (:use :cl))
(in-package :cl-mpm/examples/column)
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
;; (setf *block-compile-default* t)
(declaim (optimize (debug 3) (safety 3) (speed 2)))

(defun pescribe-velocity (sim load-mps vel)
  (let ((mps (cl-mpm:sim-mps sim)))
    (loop for mp in load-mps
          do
             (progn
               (setf (cl-mpm/particle:mp-velocity mp) vel)))))
(defun increase-load (sim load-mps amount)
  (loop for mp in load-mps
        do (with-accessors ((pos cl-mpm/particle:mp-position)
                            (force cl-mpm/particle:mp-body-force)) mp
             (incf (magicl:tref force 1 0) amount))))
(defun length-from-def (sim mp dim)
  (let* ((mp-scale 2)
         (h-initial (magicl:tref (cl-mpm/particle::mp-domain-size mp) dim 0))
         ;(h-initial  (/ (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) mp-scale))
         )
    h-initial
    ;(* (magicl:tref (cl-mpm::mp-deformation-gradient mp) dim dim) h-initial)
    ))

(defun get-disp (load-mps)
  (- (/
      (loop for mp in load-mps
            sum (+
                 (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)
                 (* 0.5d0 (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0))
                 )) (length load-mps))
     *initial-surface*))

(defun max-stress (mp)
  (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0)
  ;; (multiple-value-bind (l v) (magicl:hermitian-eig (cl-mpm::voight-to-matrix (cl-mpm/particle:mp-stress mp)))
  ;;   (apply #'max l))
  )

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
  (simple-plot-contact
   sim
   :plot :deformed
   :contact-bcs *penalty-bc*
   )
  ;; (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  ;; (multiple-value-bind (x y stress-y lx ly e)
  ;;   (loop for mp across (cl-mpm:sim-mps sim)
  ;;         collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
  ;;         collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
  ;;         collect (length-from-def sim mp 0) into lx
  ;;         collect (length-from-def sim mp 1) into ly
  ;;         ;; collect (/ (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0) 1e0) into stress-y
  ;;         collect (max-stress mp) into stress-y
  ;;         finally (return (values x y stress-y lx ly)))
  ;;   (cond
  ;;     (
  ;;      (eq plot :stress)
  ;;      (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 0.01 (apply #'max stress-y)))
  ;;      (vgplot:plot x y stress-y ";;with points pt 7 lc palette"))
  ;;     ((eq plot :deformed)
  ;;      ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 0.01 (apply #'max stress-y)))
  ;;      (vgplot:plot x y lx ly ";;with ellipses")))
  ;;   )

  ;; (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
  ;;        (ms-x (first ms))
  ;;        (ms-y (second ms))
  ;;        )
  ;;   (vgplot:axis (list 0 ms-x
  ;;                      0 ms-y))
  ;;   (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
  ;;   (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
  ;;     (vgplot:format-plot t "set ytics ~f" h)
  ;;     (vgplot:format-plot t "set xtics ~f" h))
  ;; (vgplot:replot)
  )

(declaim (notinline setup-test-column))
(defun setup-test-column (size block-size &optional (e-scale 1) (mp-scale 1))
  (let ((nd (length block-size)))
    (let* ((sim
             (cl-mpm/setup::make-simple-sim;
              (/ 1d0 e-scale)
              (mapcar (lambda (x) (* x e-scale)) size)
              :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
              ;; :sim-type 'cl-mpm::mpm-sim-usl
              ;; :args-list (list
              ;;             :enable-fbar t
              ;;             :enable-split nil)
              )
                )
           (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
           (h-x (/ h 1d0))
           (h-y (/ h 1d0))
           (density 80)
           (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
                                        ;(block-size (list (* (first size) h-x) (second block-size)))
           ;; (block-size
           ;;   (if (= nd 2)
           ;;       (list (* 1d0 h-x)
           ;;             (second block-size))
           ;;       (list (* 1d0 h-x)
           ;;             (second block-size)
           ;;             (* 1d0 h-x))))
           )
      (progn
        (let ((block-position (list (* h-x (+ 0 (/ 1d0 (* 2d0 mp-scale)) ));mp-scale->1
                                    (* h-y (+ 0 (/ 1d0 (* 2d0 mp-scale)) )))))
          ;; (print block-position)
          (cl-mpm::add-mps
           sim
           (cl-mpm/setup::make-block-mps
            ;; block-position
            (list 0 0 0)
            block-size
            (mapcar (lambda (e mp-s)
                      (round  (* e mp-s) h-x)
                      ) block-size
                        (list
                         mp-scale
                         ;; 1
                         mp-scale
                         1
                         ;; mp-scale
                                        ;mp-scale
                                        ;mp-scale
                         ))
            ;; (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                                        ;(list 1 (* e-scale mp-scale (second block-size)))
            density
            'cl-mpm/particle::particle-elastic
            :E 1d5
            :nu 0.0d0
            :gravity -10.0d0
            )))

        (setf (cl-mpm:sim-damping-factor sim)
              (* 
               0.1d0
               (cl-mpm/setup::estimate-critical-damping sim)))
        (setf (cl-mpm:sim-mass-filter sim) 1d-15)
        ;; (setf (cl-mpm::sim-ghost-factor sim) (* 1d5 1d-5))
        (setf (cl-mpm:sim-dt sim) 1d-2)
        (setf (cl-mpm:sim-bcs sim)
              (cl-mpm/bc::make-outside-bc-var (cl-mpm:sim-mesh sim)
                                              (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
                                              (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
                                              (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil nil)))
                                              (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 0 nil)))
                                              (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
                                              (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
                                              ))
        (let* ((crack-pos
               (loop for mp across (cl-mpm:sim-mps sim)
                     maximize (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)))
             (above-crack
               (loop for mp across (cl-mpm:sim-mps sim)
                     when
                     (and
                      (>= (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                          crack-pos))
                     collect mp)))
        (defparameter *terminus-mps* above-crack)
          (defparameter *initial-surface* (+ crack-pos
                                             (* 0.5d0
                                                (magicl:tref
                                                 (cl-mpm/particle::mp-domain-size
                                                  (first *terminus-mps*))
                                                 1 0))))
          (defparameter *target-displacement* 0d0))


        (defparameter *penalty-bc*
          (cl-mpm/penalty::make-bc-penalty-point-normal
           sim
           (cl-mpm/utils:vector-from-list  '(0d0 -1d0 0d0))
           (cl-mpm/utils:vector-from-list  (list 0d0 (- *initial-surface* 1d0) 0d0))
           (* 1d5 0.1)
           ;; (* 1d5 10)
           0d0))
        ;; (setf (cl-mpm::sim-bcs-force-list sim)
        ;;       (list
        ;;        (cl-mpm/bc:make-bcs-from-list
        ;;         (list
        ;;          *penalty-bc*
        ;;          ))))
        ;; (setf (cl-mpm::sim-bcs-force-list sim)
        ;;       (list
        ;;        (cl-mpm/bc:make-bcs-from-list
        ;;         (list
        ;;          (cl-mpm/bc::make-bc-closure
        ;;           nil
        ;;           (lambda ()
        ;;             (with-accessors ((mesh cl-mpm:sim-mesh)
        ;;                              (dt cl-mpm::sim-dt))
        ;;                 sim
        ;;               (let ((datum (- (+ *initial-surface* *target-displacement*)))
        ;;                     (normal (cl-mpm/utils:vector-from-list  '(0d0 -1d0 0d0))))
        ;;                 (cl-mpm/penalty::apply-displacement-control-mps
        ;;                  ;; cl-mpm/penalty::apply-force-mps
        ;;                  mesh
        ;;                  (coerce *terminus-mps* 'vector)
        ;;                  dt
        ;;                  normal
        ;;                  datum
        ;;                  (* 1d5 0.1)
        ;;                  0d0)
        ;;                 ))))))))

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
         (h (/ L e)))
    (format t "H:~E~%" h)
    (defparameter
        *sim*
      (setup-test-column (list
                          h
                          (+ L h)
                          ;; h
                          )
                         (list
                          h
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
  ;; (defparameter *run-sim* nil)
  ;; (defparameter *run-sim* t)

  ;; (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
  ;; (defparameter *load-mps-top*
  ;;   (let* ((mps (cl-mpm:sim-mps *sim*))
  ;;          (least-pos
  ;;             (apply #'max (loop for mp across mps
  ;;                                collect (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)))))
  ;;     (loop for id from 0 to (- (length mps) 1)
  ;;           when (>= (magicl:tref (cl-mpm/particle:mp-position (aref mps id)) 1 0) (- least-pos 0.001))
  ;;             collect (aref mps id))))
  ;; (defparameter *slice-mps*
  ;;   (let* ((mps (cl-mpm:sim-mps *sim*))
  ;;          (least-pos
  ;;            (apply #'max (loop for mp across mps
  ;;                               collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))))
  ;;     (loop for id from 0 to (- (length mps) 1)
  ;;           when (>= (magicl:tref (cl-mpm/particle:mp-position (aref mps id)) 0 0) (- least-pos 0.001))
  ;;             collect (aref mps id))))
  ;; (increase-load *sim* *load-mps* -100)
  ;; (increase-load *sim* *load-mps-top* -100)
  )


(defparameter *run-sim* nil)
(defun run-conv ()
  (setf *run-sim* t)
  (loop for i from 2 to 6;10
        while *run-sim*
        do
           (let ((elements (expt 2 i))
                 (final-time 15))
             (let* ((e elements)
                    (L 32d0)
                    (h (/ L e)))
               (format t "H:~E~%" h)
               (defparameter
                   *sim*
                 (setup-test-column (list
                                     h
                                     (+ L h)
                                     )
                                    (list
                                     h
                                     L
                                     )
                                    (/ 1d0 h)
                                    3
                                    )))
             ;; (defparameter *sim* (setup-test-column '(1 60) '(1 50) (/ elements 50) 2))
             (setf (cl-mpm:sim-dt *sim*)
                   (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale 0.5d0))
             (format t "Running sim size ~a ~%" elements)
             (format t "Sim dt: ~a ~%" (cl-mpm:sim-dt *sim*))
             (format t "Sim steps: ~a ~%" (/ final-time (cl-mpm:sim-dt *sim*)))
             (vgplot:figure)
             (let ((load-steps 10)
                   (substeps (round 1d0 (cl-mpm:sim-dt *sim*))))
               (time
                (loop for steps from 0 to load-steps;(round (/ final-time (cl-mpm:sim-dt *sim*)))
                      while *run-sim*
                      do
                         (progn
                           (dotimes (i substeps)
                             (cl-mpm::update-sim *sim*))
                           ;; (cl-mpm::finalise-loadstep *sim*)
                           (plot-sigma-yy)
                           (swank.live:update-swank)
                           ))))
             ;; (plot *sim*)
             (plot-sigma-yy)
             (sleep .01)
             ;; (cl-mpm/output:save-csv (merge-pathnames (format nil "conv_files/elements_~d.csv" elements)) *sim*)
             ;; (cl-mpm/output:save-vtk (merge-pathnames (format nil "conv_files/elements_~d.vtk" elements)) *sim*)
             ))
  )

(defun stop ()
  (setf *run-sim* nil))

(defun run ()
  (cl-mpm/output:save-vtk-mesh (asdf:system-relative-pathname "cl-mpm" "output/mesh.vtk")
                          *sim*)
  (defparameter *run-sim* t)
    (vgplot:close-all-plots)
    (vgplot:figure)
  (sleep 1)
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:axis (list 0 ms-x
                       0 ms-y))
    (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
    ;; (vgplot:axis (list 0 (nth 0 (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*))) 
    ;;                    0 (nth 1 (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))))
    (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
      (vgplot:format-plot t "set ytics ~f" h)
      (vgplot:format-plot t "set xtics ~f" h))
  (let* ((target-time 1d0)
         (dt (cl-mpm:sim-dt *sim*))
         (dt-scale 0.1d0)
         (substeps (floor target-time dt)))
    (cl-mpm::update-sim *sim*)
    (let* ((dt-e (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
           (substeps-e (floor target-time dt-e)))
      (format t "CFL dt estimate: ~f~%" dt-e)
      (format t "CFL step count estimate: ~D~%" substeps-e)
      (setf (cl-mpm:sim-dt *sim*) dt-e)
      (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
    (time (loop for steps from 0 to 100
                while *run-sim*
                do
                   (progn
                     (format t "Step ~d ~%" steps)
                     (time
                      (dotimes (i substeps)
                        ;; (pescribe-velocity *sim* *load-mps* (magicl:from-list '(0d0 -1d0) '(2 1) :type 'double-float))
                        ;; (pescribe-velocity *sim* *load-mps-top* (magicl:from-list '(0d0 -0.5d0) '(2 1)))
                        ;; (break)
                        ;; (increase-load *sim* *load-mps* (* -100 (cl-mpm:sim-dt *sim*)))
                        ;; (increase-load *sim* *load-mps-top* (* 100 (cl-mpm:sim-dt *sim*)))
                        (cl-mpm::update-sim *sim*)
                        ;; (cl-mpm/eigenerosion:update-fracture *sim*)
                        (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))
                        ;; (let ((h (/ (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)) 2)))
                        ;;   (setf *velocity* (cons (magicl:tref (cl-mpm/output::sample-point-velocity *sim* (list h (* h 2))) 1 0) *velocity*)))
                        ;; (setf *time*     (cons *t* *time*))
                        ))
                     (format t "Disp ~E ~%" (get-disp *terminus-mps*))

                     (let* ((dt-e (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
                            (substeps-e (floor target-time dt-e)))
                       (format t "CFL dt estimate: ~f~%" dt-e)
                       (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf (cl-mpm:sim-dt *sim*) dt-e)
                       (setf substeps substeps-e))
                     (cl-mpm/output:save-vtk (asdf:system-relative-pathname "cl-mpm" (format nil "output/sim_~5,'0d.vtk" *sim-step*))
                                             *sim*)

                     ;; (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     (incf *sim-step*)
                     (plot *sim*)
                     (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" steps)))
                     (swank.live:update-swank)
                     (sleep .01)
                     ))))
  (vgplot:figure)
  (plot-sigma-yy)
  ;; (vgplot:title "Velocity over time")
  ;; (vgplot:plot *time* *velocity*)
  )
(defun plot-sigma-yy ()
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (let* ((x-slice-pos (loop for mp across mps maximize (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
           (mp-list
             (loop for mp across mps
                   when (>= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (- x-slice-pos 0.001))
                     collect mp))
           (y (loop for mp in mp-list collect (magicl:tref (cl-mpm/particle::mp-position-trial mp) 1 0)))
           (y-ref (loop for mp in mp-list collect
                                          (magicl:tref
                                           (magicl:.-
                                            (cl-mpm/particle::mp-position mp)
                                            (cl-mpm/particle::mp-displacement mp)
                                            )
                                            1 0)))
           (syy (loop for mp in mp-list collect (magicl:tref (cl-mpm/particle::mp-stress mp) 1 0)))
           (rho 80d0)
           (E 1d5)
           (g 10d0)
           (max-y 32)
           ;; (max-y (+ (reduce #'max y-ref) (magicl:tref (cl-mpm/particle::mp-domain-size-0 (first mp-list)) 1 0)))
           (syy-ref (mapcar (lambda (x) (* rho g (- x max-y))) y-ref))
           )
      ;; (vgplot:figure)
      (vgplot:plot syy y-ref "first";"first;;with points pt 7"
                   syy-ref y-ref "reference";"Ref;;with points pt 7"
                   )
      (vgplot:legend)
      )))

(defun save-sigma-yy ()
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
