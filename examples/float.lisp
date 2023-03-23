(defpackage :cl-mpm/examples/float
  (:use :cl))
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* nil)
(in-package :cl-mpm/examples/float)
(declaim (optimize (debug 3) (safety 2) (speed 2)))

(defun ellipse-sdf (position x-l y-l)
  (let ((aspect (/ x-l y-l)))
    (lambda (pos)
      (let* ((position (magicl:from-list position '(2 1) :type 'double-float))
             (dist-vec (magicl:.* (magicl:.- position pos) (magicl:from-list (list 1 aspect) '(2 1)
                                                                             :type 'double-float)))
             (distance (sqrt (magicl:tref (magicl:@ (magicl:transpose dist-vec)
                                                    dist-vec) 0 0))))
        (- distance x-l)))))

(defun increase-load (sim load-mps amount)
  (loop for mp in load-mps
        do (with-accessors ((pos cl-mpm/particle:mp-position)
                            (force cl-mpm/particle:mp-body-force)) mp
                                        ;(incf (magicl:tref force 0 0) amount)
             (magicl:.+ force amount force)
             )))

(defun remove-sdf (sim sdf)
  (setf (cl-mpm:sim-mps sim)
        (lparallel:premove-if (lambda (mp)
                                (with-accessors ((pos cl-mpm/particle:mp-position)) mp
                                  (>= 0 (funcall sdf pos))
                                  ))
                              (cl-mpm:sim-mps sim))))

(defun length-from-def (sim mp dim)
  (let* ((mp-scale 2)
         (h-initial (magicl:tref (cl-mpm/particle::mp-domain-size mp) dim 0))
         ;(h-initial  (/ (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) mp-scale))
         )
    h-initial
    ;(* (magicl:tref (cl-mpm::mp-deformation-gradient mp) dim dim) h-initial)
    ))
(defun max-stress (mp)
  (declare (optimize (speed 2) (debug 3)))
  (multiple-value-bind (l v) (magicl:eig (cl-mpm::voight-to-matrix (cl-mpm/particle:mp-stress mp)))
    ;; (/ (apply #'max l) 1d6)
    (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0)
    ))
(defun plot (sim &optional (plot :point))
  (declare (optimize (speed 0) (debug 3) (safety 3)))
  (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  (multiple-value-bind (x y c stress-y lx ly e density temp vx)
    (loop for mp across (cl-mpm:sim-mps sim)
          collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
          collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
          collect (magicl:tref (cl-mpm::mp-velocity mp) 0 0) into vx
          collect (length-from-def sim mp 0) into lx
          collect (length-from-def sim mp 1) into ly
          collect (if (slot-exists-p mp 'cl-mpm/particle::damage) (cl-mpm/particle:mp-damage mp) 0) into c
          collect (if (slot-exists-p mp 'cl-mpm/particle::temperature) (cl-mpm/particle:mp-temperature mp) 0) into temp
          collect (if (slot-exists-p mp 'cl-mpm/particle::strain-energy-density) (cl-mpm/particle::mp-strain-energy-density mp) 0) into e
          collect (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)) into density
          ;; collect (cl-mpm/particle:mp-volume mp) into density
          collect (max-stress mp) into stress-y
          finally (return (values x y c stress-y lx ly e density temp vx)))
    (let* ((node-x '())
           (node-y '())
           (node-c '())
           (mesh (cl-mpm:sim-mesh *sim*))
           (nodes (cl-mpm/mesh:mesh-nodes mesh))
           )
      (dotimes (i (array-total-size nodes))
        (let ((n (row-major-aref nodes i)))
          (with-accessors ((index cl-mpm/mesh:node-index)
                           (boundary cl-mpm/mesh::node-boundary-node))
              n
            (when boundary
              ;; (print index)
              (destructuring-bind (x y) (cl-mpm/mesh:index-to-position mesh index)
                ;; (push (nth 0 index) node-x)
                ;; (push (nth 1 index) node-y)
                (push x node-x)
                (push y node-y)
                (push 1 node-c)
                )))))
      (cond
        ((eq plot :point)
         ;; (vgplot:format-plot t "set cbrange [0:1]")
         (vgplot:plot x y ";;with points pt 7")
         (if node-x
             (vgplot:plot
              ;;x y ";;with points pt 7"
              x y lx ly ";;with ellipses"
              node-x node-y ";;with points pt 7")
             (vgplot:plot x y ";;with points pt 7"))
         )
        ((eq plot :damage)
         (vgplot:format-plot t "set cbrange [0:1]")
         (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 0.01 (apply #'max c)))
         ;; (vgplot:plot x y c ";;with points pt 7 lc palette")
         (if node-x
             (vgplot:plot
              x y c ";;with points pt 7 lc palette"
              node-x node-y ";;with points pt 7")
             (vgplot:plot x y c ";;with points pt 7 lc palette"))
         ;; (print (length node-x))
         ;; (vgplot:plot node-x node-y ";;with points pt 7")
         )
        ((eq plot :velocity)
         (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min vx) (+ 0.01 (apply #'max vx)))
         (vgplot:plot x y vx ";;with points pt 7 lc palette")
         )
        ((eq plot :temperature)
         (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min temp) (+ 0.01 (apply #'max temp)))
         ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 0.01 (apply #'max c)))
         (vgplot:plot x y temp ";;with points pt 7 lc palette")
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
         )))
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
  (vgplot:replot))



(defun setup-test-column (size block-size block-offset &optional (e-scale 1d0) (mp-scale 1d0)
                          &rest mp-args)
  (let* ((sim (cl-mpm/setup::make-block (/ 1 e-scale)
                                        (mapcar (lambda (s) (* s e-scale)) size)
                                        #'cl-mpm/shape-function:make-shape-function-bspline)) 
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         ;(e-scale 1)
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 500)
         ;; (mass (/ (* 900 h-x h-y) (expt mp-scale 2)))
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
               'cl-mpm/particle::particle-elastic
               :E 1d6
               ;:nu 0.3250d0
               :nu 0.0d0
               ;; :visc-factor 0.1d6
               ;; :visc-power 3d0
               ;; :critical-stress 1d6
               :gravity -0.0d0
               ;; :gravity 0d0
               ;; :body-force (magicl:from-list '(0d0 100d0) '(2 1))
               ;; :gravity-axis (magicl:from-list '(0.5d0 0.5d0) '(2 1))
               :index 0
               )))
      (setf (cl-mpm:sim-damping-factor sim) 0.6d0)
      (setf (cl-mpm:sim-mass-filter sim) 1d-15)
      (setf (cl-mpm::sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm:sim-dt sim) 1d-2)

      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-var
             (cl-mpm:sim-mesh sim)
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0)))
             )
            )
      ;; (setf (cl-mpm::sim-bcs-force sim)
      ;;       (cl-mpm/bc:make-bcs-from-list
      ;;       (list (cl-mpm/bc::make-bc-closure '(0 0)
      ;;                                         (lambda ()
      ;;                                           (cl-mpm/buoyancy::apply-bouyancy sim 300d0))))))
      (setf (cl-mpm::sim-bcs-force sim)
            (cl-mpm/bc:make-bcs-from-list
             (list (cl-mpm/buoyancy::make-bc-pressure
                    sim
                    0d0
                    -1d5
                    ))))

      ;; (let ((ocean-x 0)
      ;;       (ocean-y 300))
      ;;   (setf (cl-mpm::sim-bcs-force sim)
      ;;         (cl-mpm/bc:make-bcs-from-list
      ;;          (loop for x from (floor ocean-x h) to (floor (first size) h)
      ;;                append (loop for y from 0 to (floor ocean-y h)
      ;;                             collect (cl-mpm/bc::make-bc-buoyancy
      ;;                                      (list x y)
      ;;                                      (magicl:from-list (list 0d0 (* 9.8d0 (* 1.0d0 1000))) '(2 1))))))))
      sim)))

;Setup
(defun setup ()
  (let* ((shelf-length 1000)
         (shelf-height 200)
         (shelf-bottom 120)
         (notch-length 100)
         (notch-depth 10)
         (mesh-size 10)
         )
    (defparameter *csv-name* (merge-pathnames (format nil "output/surface_position_~D.csv" mesh-size)))
    ;; (defparameter *sim* (setup-test-column '(600 600)
    ;;                                        '(100 100)
    ;;                                        '(100 200) (/ 1d0 mesh-size) 2))

    (defparameter *sim* (setup-test-column
                         (list mesh-size 200)
                         (list mesh-size 100)
                         '(000 000) (/ 1d0 mesh-size) 2)))
    ;; (remove-sdf *sim* (rectangle-sdf (list shelf-length (+ shelf-height shelf-bottom)) (list notch-length notch-depth)))
    ;; (remove-sdf *sim* (cl-mpm/setup:rectangle-sdf '(500 300) '(100 50))

  ;;Hole in plate
  ;; (defparameter *sim* (setup-test-column '(500 500) '(500 500) '(000 0) (/ 1 10) 2))
  ;; (remove-sdf *sim* (lambda (p) (+ (funcall (ellipse-sdf (list 250 250) 100 100) p))))

  ;; (defparameter *sim* (setup-test-column '(500 500) '(200 200) '(0 0) (/ 1 20) 2))

  ;; (remove-sdf *sim* (ellipse-sdf (list 0 0) 100 100))
  ;; (remove-sdf *sim* (ellipse-sdf (list 000 000) 200 200))
  ;; (remove-sdf *sim* (cl-mpm/setup:rectangle-sdf (list 000 000) '(500 200)))
  ;; (remove-sdf *sim* (cl-mpm/setup:rectangle-sdf (list 000 000) '(200 500)))

  
  ;; (defparameter *sim* (setup-test-column '(500 400) '(300 100) '(000 250) (/ 1 25) 4))
  ;; (remove-sdf *sim* (rectangle-sdf (list 1000 225) '(100 50)))

  (format t "Simulation dt ~a~%" (cl-mpm:sim-dt *sim*))
  (format t "Simulation MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  (defparameter *x* 0d0)
  (defparameter *x-pos* '())
  (defparameter *s-min* '())
  (defparameter *s-max* '())
  (defparameter *s-average* '())
  (defparameter *sim-step* 0)
  (defparameter *run-sim* nil)
  (defparameter *load-mps*
    (let* ((mps (cl-mpm:sim-mps *sim*))
           (least-pos
              (apply #'max (loop for mp across mps
                                 collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))))
      (loop for mp across mps when (>= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (- least-pos 0.001))
              collect mp)))
  ;; (increase-load *sim* *load-mps* (magicl:from-list '(-1d0 0d0) '(2 1)))
  ;; (increase-load *sim* *load-mps* 100)
  )

(defparameter *run-sim* nil)

(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)
  (defparameter *run-sim* t)
  (vgplot:close-all-plots)
  ;; (sleep 1)
  (vgplot:figure)
  ;; (vgplot:plot '(10 10))
  (sleep 1)
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:axis (list 0 ms-x
                       0 ms-y))
    (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
    (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
      (vgplot:format-plot t "set ytics ~f" h)
      (vgplot:format-plot t "set xtics ~f" h))
  (with-open-file (stream *csv-name* :direction :output :if-exists :supersede)
    (format stream "Time (s),Surface position~%")
    (loop for tim in (reverse *time*)
          for x in (reverse *x-pos*)
          do (format stream "~f, ~f ~%" tim x)))
  (defparameter *notch-position* 0.1d0)
  (time (loop for steps from 0 to 100
                while *run-sim*
                do
                (progn
                  (format t "Step ~d ~%" steps)
                  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                  ;; (cl-mpm/output:save-csv (merge-pathnames (format nil "output/simcsv_~5,'0d.csv" *sim-step*)) *sim*)

                  (push *t* *time*)
                  (setf *x*
                        (loop for mp across (cl-mpm:sim-mps *sim*)
                              maximize (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)))
                  (push
                   *x*
                   *x-pos*)
                  (push
                   (loop for mp across (cl-mpm:sim-mps *sim*)
                         maximize (magicl:tref (cl-mpm/particle::mp-stress mp) 1 0))
                   *s-max*)
                  (push
                   (loop for mp across (cl-mpm:sim-mps *sim*)
                         minimize (magicl:tref (cl-mpm/particle::mp-stress mp) 1 0))
                   *s-min*)
                  (with-accessors ((mps cl-mpm:sim-mps))
                      *sim*
                    (push
                     (/
                      (loop for mp across mps
                            sum (magicl:tref (cl-mpm/particle::mp-stress mp) 1 0))
                      (length mps))
                     *s-average*))
                  (let ((max-cfl 0))
                    (time (dotimes (i 100)
                            ;; (increase-load *sim* *load-mps* (magicl:from-list (list (* (cl-mpm:sim-dt *sim*)
                                                                                       ;; 5d0) 0d0) '(2 1)))
                            ;; (pescribe-velocity *sim* *load-mps* '(1d0 nil))
                            (cl-mpm::update-sim *sim*)
                            ;; (remove-sdf *sim* (rectangle-sdf (list 1000 330) (list *notch-position* 30)))
                            ;; (incf *notch-position* (* (cl-mpm:sim-dt *sim*) 5d0))
                           (setf *t* (+ *t* (cl-mpm::sim-dt *sim*))))))
                  (incf *sim-step*)
                  (plot *sim*)
                  (swank.live:update-swank)
                  (sleep .01)
                  (with-open-file (stream  *csv-name* :direction :output :if-exists :append)
                          (format stream "~f, ~f ~%" *t* *x*))

                  )))
    (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*))
                                          *sim*)
  ;; (cl-mpm/output:save-csv (merge-pathnames (format nil "output/simcsv_~5,'0d.csv" *sim-step*)) *sim*)

  (plot-surface)
  (plot-stress)

  ;; (with-open-file (stream (merge-pathnames "output/terminus_position.csv") :direction :output :if-exists :supersede)
  ;;   (format stream "Time (s),Terminus position~%")
  ;;   (loop for tim in (reverse *time*)
  ;;         for x in (reverse *x-pos*)
  ;;         do (format stream "~f, ~f ~%" tim x)))
  )
(defun plot-surface ()
  (vgplot:figure)
  (vgplot:title "Top surface height")
  (format t "~a~%" (- (first *x-pos*) (first (last *x-pos*))))
  (let ((init-pos (first (last *x-pos*))))
    (vgplot:plot *time* (mapcar (lambda (x) (- x init-pos)) *x-pos*))))
(defun plot-stress ()
  (vgplot:figure)
  (format t "S_max ~a~%" (first *s-max*))
  (format t "S_min ~a~%" (first *s-min*))
  (format t "S_average ~a~%" (first *s-average*))
  (vgplot:title "S_{yy} min-max evolution")
  (vgplot:plot *time* *s-min* "s-min"
               *time* *s-max* "s-max"
               *time* *s-average* "s-average"
               ))

(defun run-conv ()
  (defparameter *elements* (list))
  (defparameter *conv-stress* (list))
    (loop for i from 2 to 6
          while *run-sim*
          do
             (let ((elements (expt 2 i))
                   (final-time 100))
               (let* ((size 100)
                      (mesh-size (/ size elements)))
                 (defparameter *sim* (setup-test-column (list mesh-size (+ size mesh-size))
                                                        (list mesh-size size)
                                                        '(000 000) (/ 1d0 mesh-size) 2)))
               ;; (defparameter *sim* (setup-test-column '(1 60) '(1 50) (/ elements 50) 2))
               (setf (cl-mpm:sim-dt *sim*) (* 1d-2 (/ 16 elements)))
               (format t "Running sim size ~a ~%" elements)
               (format t "Sim dt: ~a ~%" (cl-mpm:sim-dt *sim*))
               (format t "Sim steps: ~a ~%" (/ final-time (cl-mpm:sim-dt *sim*)))
               (time
                (loop for steps from 0 to (round (/ final-time (cl-mpm:sim-dt *sim*)))
                      do
                         (cl-mpm::update-sim *sim*)))
               (plot *sim*)
               (with-accessors ((mps cl-mpm:sim-mps))
                   *sim*
                 (push
                  (/
                   (loop for mp across mps
                         sum (magicl:tref (cl-mpm/particle::mp-stress mp) 1 0))
                   (length mps))
                  *conv-stress*))
               (push elements *elements*)

               (sleep .01)
               (cl-mpm/output:save-csv (merge-pathnames (format nil "conv_files/elements_~d.csv" elements)) *sim*)
               (cl-mpm/output:save-vtk (merge-pathnames (format nil "conv_files/elements_~d.vtk" elements)) *sim*)))
  (plot-conv)
  )
(defun plot-conv ()
  (vgplot:figure)
  (vgplot:title "Convergance stress")
  (vgplot:loglog *elements* (mapcar (lambda (x) (abs (- x -1d3))) *conv-stress*))
  )

(setf lparallel:*kernel* (lparallel:make-kernel 8 :name "custom-kernel"))

(defun profile ()
  (sb-profile:unprofile)
  (sb-profile:reset)
  (sb-profile:profile cl-mpm::update-sim
                      cl-mpm::reset-grid
                      cl-mpm::p2g
                      cl-mpm::filter-grid
                      cl-mpm::update-nodes
                      cl-mpm::apply-bcs
                      cl-mpm::g2p
                      cl-mpm::update-particle
                      cl-mpm::update-stress
                      cl-mpm::iterate-over-neighbours-shape
                      cl-mpm::iterate-over-neighbours-shape-linear
                      cl-mpm::p2g-mp
                      cl-mpm::g2p-mp
                      cl-mpm::p2g-mp-node
                      cl-mpm::g2p-mp-node
                      ;; cl-mpm::update-strain-kirchoff
                      cl-mpm/damage::calculate-damage
                      cl-mpm/damage::apply-damage
                      cl-mpm/damage::delocalise-damage
                      cl-mpm/damage::create-delocalisation-list
                      ;; cl-mpm/eigenerosion:update-fracture
                      ;; cl-mpm/eigenerosion::remove-material-damaged
                      ;; cl-mpm/eigenerosion::find-neighbours
                      )
  (loop repeat 100
        do (progn
             (cl-mpm::update-sim *sim*)
             ;; (cl-mpm/damage::calculate-damage (cl-mpm:sim-mesh *sim*)
             ;;                                  (cl-mpm:sim-mps *sim*)
             ;;                                  (cl-mpm:sim-dt *sim*)
             ;;                                  25d0)
             ;; (cl-mpm/eigenerosion:update-fracture *sim*)
             ))
  (sb-profile:report))
(defun simple-time ()
  (setf lparallel:*kernel* (lparallel:make-kernel 8 :name "custom-kernel"))
  (let ((mp (cl-mpm/particle:make-particle 2 'cl-mpm/particle:particle
                                           :size (magicl:from-list '(1d0 1d0) '(2 1)))))
    (format t "Old")
      ;; (time
      ;;  (let ((m (magicl:zeros '(2 2))))
      ;;  (lparallel:pdotimes (i 100000)
      ;;    (magicl::eighth m)))
      ;; )
  (time
   ;; (lparallel:pdotimes (i 1000000)
   (dotimes (i 1000)
     (with-accessors ((mesh cl-mpm:sim-mesh)
                      (mps cl-mpm:sim-mps)
                      (dt cl-mpm:sim-dt))
         *sim*
       (cl-mpm::update-sim *sim*)
       ;; (cl-mpm::update-stress-mp mesh (aref mps 0) dt)
         )
     ))
    (format t "new")
    ))
(defmacro time-form (it form)
  `(progn
     (declaim (optimize speed))
     (let* ((iterations ,it)
            (start (get-internal-real-time)))
       (dotimes (i ,it)
         ,form)
       (let* ((end (get-internal-real-time))
              (units internal-time-units-per-second)
              (dt (/ (- end start) (* iterations units)))
              )
         (format t "Total time: ~f ~%" (/ (- end start) units)) (format t "Time per iteration: ~f~%" (/ (- end start) (* iterations units)))
         (format t "Throughput: ~f~%" (/ 1 dt))
         dt))))
(defun time-parts ()
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (dt cl-mpm::dt)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               )
                *sim*
                (;time-form 1000
                 time
                           (dotimes (i 100)
                             (cl-mpm::reset-grid mesh)
                             (cl-mpm::p2g mesh mps)
                             (when (> mass-filter 0d0)
                               (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter *sim*)))
                             (cl-mpm::update-node-kinematics mesh dt)
                             (cl-mpm::apply-bcs mesh bcs dt)
                             (cl-mpm::update-stress mesh mps dt)
                             ;;   (when enable-damage
                             ;;     (cl-mpm/damage::calculate-damage mesh
                             ;;                                      mps
                             ;;                                      dt
                             ;;                                      50d0))

                             ;;   ;Map forces onto nodes
                               (cl-mpm::p2g-force mesh mps)
                               (cl-mpm::apply-bcs mesh bcs-force dt)
                               (cl-mpm::update-node-forces mesh (cl-mpm::sim-damping-factor *sim*) dt)
                             ;;   ;Reapply velocity BCs
                               (cl-mpm::apply-bcs mesh bcs dt)
                             ;;   ;;Also updates mps inline
                               (cl-mpm::g2p mesh mps dt)

                             ;;   (when remove-damage
                             ;;     (cl-mpm::remove-material-damaged *sim*))
                             ;;   (when split
                             ;;     (cl-mpm::split-mps *sim*))
                             (cl-mpm::check-mps mps)
                             )
                    )))
;(let ((stress (magicl:zeros '(3 1)))
;      (strain (magicl:zeros '(3 1)))
;      (de (cl-mpm/constitutive::linear-elastic-matrix 1d0 0d0))
;      )
;  (format t "mat")
;  (time
;   (dotimes (i 100000)
;       (cl-mpm/constitutive::maxwell-exp strain stress 1d0 0d0 1d0 1d0)))
;  (format t "voight")
;  (time
;   (dotimes (i 100000)
;     (cl-mpm/constitutive::maxwell-exp-v strain stress 1d0 0d0 1d0 1d0)))
;  (format t "simd")
;  (time
;   (dotimes (i 100000)
;     (cl-mpm/constitutive::maxwell-exp-v-simd strain stress 1d0 0d0 de 1d0  1d0)))
;  )

;; (cl-mpm/mesh::cell-iterate-over-neighbours *mesh* *cell*
;;                                            (lambda (m c n w g)
;;                                              (print n)
;;                                              (let ((n-pos (cl-mpm/mesh::node-position n))
;;                                                    (cent (cl-mpm/mesh::cell-centroid c)))

;;                                                ;; (print n)
;;                                                (print (magicl:.- cent n-pos))
;;                                                )
;;                                              )
;;                                            )
