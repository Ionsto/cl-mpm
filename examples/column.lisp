(defpackage :cl-mpm/examples/column
  (:use :cl))
(in-package :cl-mpm/examples/column)
(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)
(declaim (optimize (debug 3) (safety 3) (speed 2)))
(ql:quickload "vgplot")

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

(defun max-stress (mp)
  (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0)
  ;; (multiple-value-bind (l v) (magicl:eig (cl-mpm::voight-to-matrix (cl-mpm/particle:mp-stress mp)))
  ;;   (apply #'max l))
  )
(defun plot (sim &optional (plot :stress))
  (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  (multiple-value-bind (x y stress-y lx ly e)
    (loop for mp across (cl-mpm:sim-mps sim)
          collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
          collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
          collect (length-from-def sim mp 0) into lx
          collect (length-from-def sim mp 1) into ly
          ;; collect (/ (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0) 1e0) into stress-y
          collect (max-stress mp) into stress-y
          finally (return (values x y stress-y lx ly)))
    (cond
      (
       (eq plot :stress)
       (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 0.01 (apply #'max stress-y)))
       (vgplot:plot x y stress-y ";;with points pt 7 lc palette"))
      ((eq plot :deformed)
       ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 0.01 (apply #'max stress-y)))
       (vgplot:plot x y lx ly ";;with ellipses")))
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

(defun setup-test-column (size block-size &optional (e-scale 1) (mp-scale 1))
  (let* ((sim (cl-mpm/setup::make-block (/ 1d0 e-scale)
                                        (list (* 1 (first size)) (* e-scale (second size)))
                                        #'cl-mpm/shape-function::make-shape-function-linear)) 
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 100)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
        (block-size (list (* (first size) h-x) (second block-size)))
         )
    (progn
      (let ((block-position (list (* h-x (+ 0 (/ 1d0 (* 2d0 mp-scale)) ));mp-scale->1
                                  (* h-y (+ (/ 1d0 (* 2d0 mp-scale)) )))))
        (print block-position)
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-block-mps
               block-position
               block-size
               (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
               ;(list 1 (* e-scale mp-scale (second block-size)))
               density
               'cl-mpm::make-particle
               'cl-mpm/particle::particle-elastic
               :E 1d6 :nu 0.0d0)))

      (loop for mp across (cl-mpm::sim-mps sim)
            do (setf (cl-mpm/particle::mp-gravity mp) -10.0d0))
      (setf (cl-mpm:sim-damping-factor sim) 1d0)
      (setf (cl-mpm:sim-mass-filter sim) 1d-15)
      (setf (cl-mpm:sim-dt sim) 1d-3)
      ;; (setf (cl-mpm:sim-bcs sim) (cl-mpm/bc:make-outside-bc-nostick (cl-mpm/mesh:mesh-count (cl-mpm:sim-mesh sim))))
      ;; (setf (cl-mpm:sim-bcs sim) '())
      ;; (setf (cl-mpm:sim-bcs sim) (cl-mpm/bc:make-outside-bc (cl-mpm/mesh:mesh-count (cl-mpm:sim-mesh sim))))
      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-var (cl-mpm:sim-mesh sim)
                                           (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil)))
                                           (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil)))
                                           (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil)))
                                           (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil 0)))
                                           ))
      sim)))
(setf lparallel:*kernel* (lparallel:make-kernel 8 :name "custom-kernel"))
;Setup
(defun setup ()
  (defparameter *sim* (setup-test-column '(1 60) '(1 50) (/ 1 5) 4))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  (defparameter *sim-step* 0)
  (defparameter *run-sim* nil)
  (defparameter *run-sim* t)
  (defparameter *load-mps-top*
    (let* ((mps (cl-mpm:sim-mps *sim*))
           (least-pos
              (apply #'max (loop for mp across mps
                                 collect (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)))))
      (loop for id from 0 to (- (length mps) 1)
            when (>= (magicl:tref (cl-mpm/particle:mp-position (aref mps id)) 1 0) (- least-pos 0.001))
              collect (aref mps id))))
  (defparameter *slice-mps*
    (let* ((mps (cl-mpm:sim-mps *sim*))
           (least-pos
             (apply #'max (loop for mp across mps
                                collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))))
      (loop for id from 0 to (- (length mps) 1)
            when (>= (magicl:tref (cl-mpm/particle:mp-position (aref mps id)) 0 0) (- least-pos 0.001))
              collect (aref mps id))))
  ;; (increase-load *sim* *load-mps* -100)
  ;; (increase-load *sim* *load-mps-top* -100)
  )


(defparameter *run-sim* nil)

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
    (time (loop for steps from 0 to 100
                while *run-sim*
                do
                (progn
                  (format t "Step ~d ~%" steps)
                  (time
                   (dotimes (i 100)
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
                  (cl-mpm/output:save-vtk (asdf:system-relative-pathname "cl-mpm" (format nil "output/sim_~5,'0d.vtk" *sim-step*))
                                          *sim*)
                  (incf *sim-step*)
                  (plot *sim*)
                  (vgplot:print-plot (asdf:system-relative-pathname "cl-mpm" (format nil "output/frame_~5,'0d.png" steps)))
                  (swank.live:update-swank)
                  (sleep .01)
                  )))
    ;; (vgplot:figure)
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
           (y (loop for mp in mp-list collect (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)))
           (syy (loop for mp in mp-list collect (magicl:tref (cl-mpm/particle::mp-stress mp) 1 0))))
      (vgplot:figure)
      (vgplot:plot y syy ";;with points pt 7"))))
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
