(defpackage :cl-mpm/examples/shear
  (:use :cl))
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* nil)

(in-package :cl-mpm/examples/shear)
(declaim (optimize (debug 2) (safety 2) (speed 2)))



(defun pescribe-velocity (sim load-mps vel)
  (let ((mps (cl-mpm:sim-mps sim)))
    (loop for mp in load-mps
          do
             (progn
               (setf (cl-mpm/particle:mp-velocity mp) vel)))))

(defun length-from-def (sim mp dim)
  (let* ((mp-scale 2)
         (h-initial (magicl:tref (cl-mpm/particle::mp-domain-size mp) dim 0))
         )
    h-initial))
(defun max-stress (mp)
  (multiple-value-bind (l v) (magicl:eig (cl-mpm::voight-to-matrix (cl-mpm/particle:mp-stress mp)))
    ;; (apply #'max l)
    (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0)
    ))
(defun plot-vel (sim)
  (let* ((node-x '())
         (node-y '())
         (node-dx '())
         (node-dy '())
         (mesh (cl-mpm:sim-mesh *sim*))
         (nodes (cl-mpm/mesh:mesh-nodes mesh))
         )
    (dotimes (i (array-total-size nodes))
      (let ((n (row-major-aref nodes i)))
        (with-accessors ((index cl-mpm/mesh:node-index)
                         (boundary cl-mpm/mesh::node-boundary-node)
                         (vel cl-mpm/mesh::node-velocity)
                         (active cl-mpm/mesh::node-active)
                         )

            n
          (when active
            (destructuring-bind (x y) (cl-mpm/mesh:index-to-position mesh index)
              (push x node-x)
              (push y node-y)
              (push (magicl:tref vel 0 0) node-dx)
              (push (magicl:tref vel 1 0) node-dy)
              )
            ))))
    (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
    (multiple-value-bind (x y c stress-y lx ly e density)
        (loop for mp across (cl-mpm:sim-mps sim)
              collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
              collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
              collect (length-from-def sim mp 0) into lx
              collect (length-from-def sim mp 1) into ly
              collect (if (slot-exists-p mp 'cl-mpm/particle::damage) (cl-mpm/particle:mp-damage mp) 0) into c
              collect (if (slot-exists-p mp 'cl-mpm/particle::strain-energy-density) (cl-mpm/particle::mp-strain-energy-density mp) 0) into e
              collect (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)) into density
              ;; collect (cl-mpm/particle:mp-volume mp) into density
              collect (max-stress mp) into stress-y
              finally (return (values x y c stress-y lx ly e density)))
      (vgplot:plot x y lx ly ";;with ellipses"
                   node-x node-y node-dx node-dy ";;with vectors")
      ))
  (vgplot:replot)
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:axis (list 0 ms-x
                       0 ms-y))
    (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
  (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
    (vgplot:format-plot t "set ytics ~f" h)
    (vgplot:format-plot t "set xtics ~f" h))
  (vgplot:replot)

  )
(defun plot (sim &optional (plot :deformed))
  (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  (multiple-value-bind (x y c stress-y lx ly e density)
    (loop for mp across (cl-mpm:sim-mps sim)
          collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
          collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
          collect (length-from-def sim mp 0) into lx
          collect (length-from-def sim mp 1) into ly
          collect (if (slot-exists-p mp 'cl-mpm/particle::damage) (cl-mpm/particle:mp-damage mp) 0) into c
          collect (if (slot-exists-p mp 'cl-mpm/particle::strain-energy-density) (cl-mpm/particle::mp-strain-energy-density mp) 0) into e
          collect (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)) into density
          ;; collect (cl-mpm/particle:mp-volume mp) into density
          collect (max-stress mp) into stress-y
          finally (return (values x y c stress-y lx ly e density)))
    (cond
      ((eq plot :damage)
       (vgplot:format-plot t "set cbrange [0:1]")
       ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 0.01 (apply #'max c)))
       (vgplot:plot x y c ";;with points pt 7 lc palette")
       )
      ((eq plot :energy)
       (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min e) (+ 0.01 (apply #'max e)))
       (vgplot:plot x y e ";;with points pt 7 lc palette")
       )
      (
       (eq plot :stress)
       (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 0.01 (apply #'max stress-y)))
       (vgplot:plot x y stress-y ";;with points pt 7 lc palette"))
      ((eq plot :deformed)
       ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 0.01 (apply #'max stress-y)))
       (vgplot:plot x y lx ly ";;with ellipses"))
      ((eq plot :density)
       (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min density) (+ 0.01 (apply #'max density)))
       (vgplot:plot x y e ";;with points pt 7 lc palette")
       ))
    )
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:axis (list 0 ms-x
                       0 ms-y))
    (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
    (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
      (vgplot:format-plot t "set ytics ~f" h)
      (vgplot:format-plot t "set xtics ~f" h))
  (vgplot:replot)
  )

(defun increase-load (sim load-mps amount)
  (loop for mp in load-mps
        do (with-accessors ((pos cl-mpm/particle:mp-position)
                            (force cl-mpm/particle:mp-body-force)) mp
             (incf (magicl:tref force 0 0) amount))))
(defun setup-test-column (size block-size block-offset &optional (e-scale 1d0) (mp-scale 1d0) (particle-type 'cl-mpm/particle::particle-elastic))
  (let* ((sim (cl-mpm/setup::make-block (/ 1 e-scale)
                                        (mapcar (lambda (s) (* s e-scale)) size)
                                        #'cl-mpm/shape-function:make-shape-function-linear))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         ;(e-scale 1)
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 1) ;(mass (/ (* 100 h-x h-y) (expt mp-scale 2)))
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size)))
    (progn
      (let ((block-position
              (mapcar #'+ (list (* h-x (- (+ (/ 1 (* 2 mp-scale))) 0))
                                (* h-y (+ (/ 1d0 (* 2d0 mp-scale)))))
                      block-offset)))
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-block-mps
               block-position
               block-size
               (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
               density
               'cl-mpm::make-particle
               particle-type
               :E 1d4
               :nu 0.3d0
               :gravity -0.0d0
               )))
      (setf (cl-mpm:sim-damping-factor sim) 0.01d0)
      ;; (setf (cl-mpm:sim-damping-factor sim) 1d-5)
      (setf (cl-mpm:sim-mass-filter sim) 1d-15)
      ;; (setf (cl-mpm:sim-mass-filter sim) 0d0)
      (setf (cl-mpm:sim-dt sim) 1d-2)
      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-var (cl-mpm:sim-mesh sim)
                                            (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil nil)))
                                            (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil nil)))
                                            (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil)))
                                            (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil nil)))
                                           ))
      (with-accessors ((mps cl-mpm:sim-mps))
          sim
        (macrolet ((loopovermps (type dim)
                     `(loop for mp across mps
                            ,type
                             (funcall #'magicl:tref
                                      (cl-mpm/particle:mp-position mp)
                                      ,dim 0))))
          (let* ((xmin (loopovermps minimize 0))
                 (xmax (loopovermps maximize 0))
                 (ymin (loopovermps minimize 1))
                 (ymax (loopovermps maximize 1))
                 )
            (defparameter *load-mps*
              (loop for mp across mps
                    when
                    (with-accessors ((pos cl-mpm/particle:mp-position))
                        mp
                      (or
                       t
                       (= (magicl:tref pos 0 0) xmin)
                       (= (magicl:tref pos 0 0) xmax)
                       (= (magicl:tref pos 1 0) ymin)
                       (= (magicl:tref pos 1 0) ymax)
                       ))
                    collect mp
                    )))
          ))
      ;; (increase-load sim *load-mps* 1d5)
      (defparameter *shear-rate* 0.1d0)
      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc:make-bcs-from-list
             (append
              (map 'list #'identity (cl-mpm:sim-bcs sim))
              (list
               (cl-mpm/bc::make-bc-closure
                '(0 0)
                (lambda ()
                  (apply-simple-shear
                   sim
                   *load-mps*
                   *shear-rate*)
                  )
                )))))
      sim)))

;Setup
(defun setup ()
  (declare (optimize (speed 0)))
  (let ((mesh-size 2.0)
        (mps-per-cell 6))
    (defparameter *sim* (setup-test-column (list (* 8 10) 8) '(8 8) '(0 0) (/ 1 mesh-size) mps-per-cell)))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  (defparameter *s-xx* '())
  (defparameter *s-yy* '())
  (defparameter *s-xy* '())
  (defparameter *sim-step* 0)
  (defparameter *run-sim* nil)
  (defparameter *run-sim* t)
  )

(defun mag-sq (a)
  (magicl::sum (magicl:.* a a)))

(defun calculate-energy-gravity (sim)
  (loop for mp across (cl-mpm:sim-mps sim)
        sum (* 9.8 (magicl:tref (cl-mpm/particle::mp-displacement mp) 1 0) (cl-mpm/particle:mp-mass mp))))

(defun calculate-energy-kinetic (sim)
  (loop for mp across (cl-mpm:sim-mps sim)
        sum (* 0.5 (cl-mpm/particle:mp-mass mp) (mag-sq (cl-mpm/particle:mp-velocity mp)))))

(defun calculate-energy-strain (sim)
  (loop for mp across (cl-mpm:sim-mps sim)
        sum
        (with-accessors ((volume cl-mpm/particle:mp-volume)
                         (stress cl-mpm/particle:mp-stress)
                         (strain cl-mpm/particle:mp-strain)
                         (def cl-mpm/particle::mp-deformation-gradient)
                         )
            mp
          (* 0.5 volume (magicl::sum (magicl:.* strain (magicl:scale stress (magicl:det def)))))
          ;; (* 0.5 volume (magicl::sum (magicl:.* strain (magicl:scale stress 1d0))))
          )))

(defun apply-simple-shear (sim load-mps shear-rate)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (loop for mp in load-mps
          do
             (cl-mpm::iterate-over-neighbours
              mesh
              mp
              (lambda (mesh mp node svp grad)
                (with-accessors (
                                 ;; (vel cl-mpm/particle:mp-velocity)
                                 )
                    mp
                  (with-accessors ((pos cl-mpm/mesh::node-position)
                                   (vel cl-mpm/mesh::node-velocity)
                                   (active cl-mpm/mesh::node-active)
                                   (acc cl-mpm/mesh::node-acceleration)
                                   )
                      node
                    (when active
                      (setf (magicl:tref vel 0 0)
                            (* shear-rate
                               (magicl:tref pos 1 0))
                            (magicl:tref vel 1 0)
                            0d0
                            )
                      ;; (setf (magicl:tref acc 0 0) 0d0
                      ;;       (magicl:tref acc 1 0) 0d0
                      ;;       )
                      )
                    ))))
             ;; (progn
             ;;   (setf (cl-mpm/particle:mp-velocity mp) vel))
          )))

(defparameter *run-sim* nil)

(defun run ()
  (declare (optimize (speed 1) (debug 2)))
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)
  (defparameter *run-sim* t)
    (vgplot:close-all-plots)
    (vgplot:figure)
  (format t "MP count ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (defparameter *energy-gpe* '())
  (defparameter *energy-ke* '())
  (defparameter *energy-se* '())
  (defparameter *s-xx* '())
  (defparameter *s-yy* '())
  (defparameter *s-xy* '())
  (sleep 1)
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:axis (list 0 ms-x
                       0 ms-y))
    (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
    (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
          (sub-steps (/ (/ 0.1 *shear-rate*) (cl-mpm:sim-dt *sim*))))
      ;; (vgplot:format-plot t "set ytics ~f" h)
      (vgplot:format-plot t "set xtics ~f" h)
      (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*))
                              *sim*)
      (incf *sim-step*)
      (format t "Substeps for ~a~%" sub-steps)
      (time (loop for steps from 0 to 50
                  while *run-sim*
                  do
                     (progn
                       (format t "Step ~d ~%" steps)
                       (let ((max-cfl 0))
                         (time (loop for i from 0 to (- sub-steps 1)
                                     while *run-sim*
                                     do
                                        (progn
                                          (cl-mpm::update-sim *sim*)
                                          (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))))
            (push *t* *time*)
                       (push (calculate-energy-strain *sim*) *energy-se*)
                       (push (calculate-energy-kinetic *sim*) *energy-ke*)
                       (push (calculate-energy-gravity *sim*) *energy-gpe*)
                       (with-accessors ((mps cl-mpm:sim-mps))
                           *sim*
                         (push
                          (/
                           (loop for mp across mps
                                 sum (magicl:tref (cl-mpm/particle::mp-stress mp) 0 0))
                           (length mps))
                          *s-xx*)
                         (push
                          (/
                           (loop for mp across mps
                                 sum (magicl:tref (cl-mpm/particle::mp-stress mp) 1 0))
                           (length mps))
                          *s-yy*)
                         (push
                          (/
                           (loop for mp across mps
                                 sum (magicl:tref (cl-mpm/particle::mp-stress mp) 2 0))
                           (length mps))
                          *s-xy*)
                         )

                       (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*))
                                               *sim*)
                       (incf *sim-step*)
                       (plot *sim*)
                       ;; (vgplot:print-plot (asdf:system-relative-pathname "cl-mpm" (format nil "output/frame_~5,'0d.png" steps)))
                       (swank.live:update-swank)
                       (sleep .1)

                       )))
      (cl-mpm/output:save-csv (merge-pathnames (format nil "output/final.csv" *sim-step*)) *sim*)
      (with-open-file (stream (merge-pathnames "output/energy.csv") :direction :output :if-exists :supersede)
        (format stream "Time (s),SE,KE,GPE~%")
        (loop for tim in (reverse *time*)
              for se in (reverse *energy-se*)
              for ke in (reverse *energy-ke*)
              for gpe in (reverse *energy-gpe*)
              do (format stream "~f, ~f, ~f, ~f ~%" tim se ke gpe)))
      ;; (vgplot:figure)
      ;; (vgplot:title "Velocity over time")
      ;; (vgplot:plot *time* *velocity*)
      (plot-energy)
      (plot-stress)
      ))
(defun plot-stress ()
  (vgplot:figure)
  (vgplot:title "Stress")
  (let* ((E 1d4)
         (shear (mapcar (lambda (x) (* x *shear-rate*)) *time*))
         (sxy-an (shear-stress-analytic E 0.3 shear)))
    (print (length shear))
    (print (length sxy-an))
    (vgplot:plot shear (mapcar (lambda (x) (/ x E)) *s-xx*) "s-xx"
                 shear (mapcar (lambda (x) (/ x E)) *s-yy*) "s-yy"
                 shear (mapcar (lambda (x) (/ x E)) *s-xy*) "s-xy"
                 shear (mapcar (lambda (x) (/ x E)) sxy-an) "s-xy-analytic"
                 )))
(defun plot-energy ()
  (vgplot:figure)
  (vgplot:title "Energy over time")
  (vgplot:plot *time* *energy-gpe* "GPE"
               *time* *energy-ke* "KE"
               *time* *energy-SE* "SE"
               *time* (mapcar #'+ *energy-gpe* *energy-ke* *energy-se*) "total energy"
               ))
(defun test-objective-rates ()
  (declare (optimize (speed 1) (debug 2)))
  (defparameter *run-sim* t)
    (vgplot:close-all-plots)
    (vgplot:figure)
  (format t "MP count ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (sleep 1)
  (setup)
  (defparameter *stress-xy-table* '())
  (defparameter *table-names* '())
  (defparameter *table-times* '())
  (loop for name in '("true" "inc" "logspin" "jaumann")
        for model in
        '(cl-mpm/particle::particle-elastic
          cl-mpm/particle::particle-elastic-inc
          cl-mpm/particle::particle-elastic-logspin
          cl-mpm/particle::particle-elastic-jaumann
          ;; cl-mpm/particle::particle-elastic-truesdale
          )
        do
           (progn
             (let ((mesh-size 1)
                   (mps-per-cell 2))
               (defparameter *sim* (setup-test-column (list (* 8 10) 8) '(8 8) '(0 0) (/ 1 mesh-size) mps-per-cell model)))
             (defparameter *t* 0)
             (defparameter *sim-step* 0)
             (defparameter *time* '())

             ;; (cl-mpm/output:save-vtk-mesh (merge-pathnames (format nil"output_~a/mesh.vtk" name))
             ;;                              *sim*)
             (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
                    (ms-x (first ms))
                    (ms-y (second ms))
                    )
               (vgplot:axis (list 0 ms-x
                                  0 ms-y))
               (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
             (format t "Running model ~a ~%" name)
             (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
                   (xy (list)))
               (vgplot:format-plot t "set xtics ~f" h)
               ;; (cl-mpm/output:save-vtk (merge-pathnames (format nil "output_~a/sim_~5,'0d.vtk" name *sim-step*)) *sim*)
               (incf *sim-step*)
               (time (loop for steps from 0 to 80
                           while *run-sim*
                           do
                              (progn
                                (format t "Step ~d ~%" steps)
                                (let ((max-cfl 0))
                                  (time (loop for i from 0 to 99
                                              while *run-sim*
                                              do
                                                 (progn
                                                   (cl-mpm::update-sim *sim*)
                                                   (setf *t* (+ *t* (cl-mpm::sim-dt *sim*))))
                                              ))
                                  )
                                (push *t* *time*)
                                (with-accessors ((mps cl-mpm:sim-mps))
                                    *sim*
                                  (push
                                   (/
                                    (loop for mp across mps
                                          sum (magicl:tref (cl-mpm/particle::mp-stress mp) 2 0))
                                    (length mps))
                                   xy))
                                ;; (cl-mpm/output:save-vtk (merge-pathnames (format nil "output_~a/sim_~5,'0d.vtk" name *sim-step*)) *sim*)
                                (incf *sim-step*)
                                (plot *sim*)
                                (swank.live:update-swank)
                                (sleep .01)

                                )))
               (push
                xy
                *stress-xy-table*
                )
               (push
                name
                *table-names*
                )
               ;; (cl-mpm/output:save-csv (merge-pathnames (format nil "output_~a/final.csv" name *sim-step*)) *sim*)
               )))
  (plot-stress-table))
(defun mat-log (m)
    (multiple-value-bind (l v) (magicl:eig m)
      (magicl:@
         v
         (magicl:from-diag (mapcar (lambda (x) (the double-float (log (the double-float x)))) l) :type 'double-float)
         (magicl:transpose v))))

(defun shear-stress-analytic (E nu shear)
  (loop for s in shear
        collect
        (let* ((F (magicl:from-list (list 1d0 0d0 s 1d0) '(2 2)))
               (b (mat-log (magicl:@ F (magicl:transpose F))))
               (s  (cl-mpm/constitutive:linear-elastic (cl-mpm/utils:matrix-to-voight b) E nu))
               )
          ;; b
          (magicl:tref s 2 0)
          )))
(defun plot-stress-table ()
  (vgplot:figure)
  ;; (vgplot:title "Stress")
  (vgplot:xlabel "Shear ratio")
  (vgplot:ylabel "Normalised shear stress")
  (let* ((E 1d4)
         (shear (mapcar (lambda (x) (* x *shear-rate*)) *time*))
         (sxy-an (shear-stress-analytic E 0.3 shear)))
    ;; (print shear)
    (print sxy-an)
  (vgplot:plot
   shear (mapcar (lambda (x) (/ x E)) (nth 0 *stress-xy-table*)) (nth 0 *table-names*)
   shear (mapcar (lambda (x) (/ x E)) (nth 1 *stress-xy-table*)) (nth 1 *table-names*)
   shear (mapcar (lambda (x) (/ x E)) (nth 2 *stress-xy-table*)) (nth 2 *table-names*)
   shear (mapcar (lambda (x) (/ x E)) (nth 3 *stress-xy-table*)) (nth 3 *table-names*)
   shear (mapcar (lambda (x) (/ x E)) sxy-an) "analytic"
               )))

(setf lparallel:*kernel* (lparallel:make-kernel 8 :name "custom-kernel"))
