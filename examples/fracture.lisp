(defpackage :cl-mpm/examples/fracture
  (:use :cl))
(in-package :cl-mpm/examples/fracture)
(declaim (optimize (debug 0) (safety 0) (speed 3)))
;; (declaim (optimize (debug 3) (safety 3) (speed 0)))
;; (ql:quickload "cl-mpm")
;; (ql:quickload "cl-mpm/setup")
;; (ql:quickload "cl-mpm/particle")
;; (ql:quickload "cl-mpm/bc")
(ql:quickload "vgplot")
;; (ql:quickload "swank.live")
;; (ql:quickload "cl-mpm/output")
;; (ql:quickload "magicl")
;; (ql:quickload "py4cl")
;; (setf py4cl:*python-command* "python3")
;; (py4cl:import-module "matplotlib.pyplot" :as "plt")
;; (ql:quickload "py4cl2")
;; (defun setup ()
;;   (py4cl2:defpymodule "matplotlib" nil :lisp-package "MPL")
;;   (mpl:use :backend "TKAgg")
;;   (py4cl2:defpymodule "matplotlib.pyplot" nil :lisp-package "PLT"))

(defun length-from-def (sim mp dim)
  (let* ((mp-scale 1)
         (h-initial  (/ (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) mp-scale)))
    (* (magicl:tref (cl-mpm::mp-deformation-gradient mp) dim dim) h-initial)))
(defun plot-lambda (sim &optional (c-val (lambda (x) (cl-mpm/particle:mp-damage x))))
  (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  (multiple-value-bind (x y c lx ly)
    (loop for mp across (cl-mpm:sim-mps sim)
          collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
          collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
          collect (length-from-def sim mp 0) into lx
          collect (length-from-def sim mp 1) into ly
          collect (funcall c-val mp) into c
          finally (return (values x y c lx ly)))
       (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 0.01 (apply #'max c)))
    ;; (vgplot:plot x y c ";;with points pt 7 lc palette")
    (break)
    (vgplot:plot x y lx ly ";;with ellipses")
    )
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:axis (list 0 ms-x
                       0 ms-y))
    (vgplot:format-plot t "set size ratio ~f" (/ ms-x ms-y)))
  (vgplot:replot))
(defun plot (sim &optional (plot :stress))
  (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  (multiple-value-bind (x y c stress-y lx ly)
    (loop for mp across (cl-mpm:sim-mps sim)
          collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
          collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
          collect (length-from-def sim mp 0) into lx
          collect (length-from-def sim mp 1) into ly
          collect (cl-mpm/particle:mp-damage mp) into c
          collect (/ (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0) 1e0) into stress-y
          finally (return (values x y c stress-y lx ly)))
    (cond
      ((eq plot :damage)
        (vgplot:format-plot t "set cbrange [0:1]")
        ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 0.01 (apply #'max c)))
        (vgplot:plot x y c ";;with points pt 7 lc palette")
       )
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
  ;; (vgplot:axis (list 0 (nth 0 (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim))) 
  ;;                    0 (nth 1 (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))))
  (vgplot:replot))
;; (defun plot-pyplot (sim)
;;   (multiple-value-bind (x y c)
;;     (loop for mp across (cl-mpm:sim-mps sim)
;;           collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
;;           collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
;;           collect (cl-mpm/particle:mp-damage mp) into c
;;           finally (return (values x y c)))
;;     (plt:scatter :x x :y y :c c))
;;   (plt:xlim 0 (first (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim))))
;;   (plt:ylim 0 (second (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim))))
;;   (plt:clim 0 1)
;;   ;; (plt:colorbar)
;;   (plt:show))

(defun remove-hole (sim position size)
      (setf (cl-mpm:sim-mps sim)
            (remove-if (lambda (mp)
                         (with-accessors ((pos cl-mpm/particle:mp-position)) mp
                           (let* ((position (magicl:from-list position '(2 1)))
                                  (dist-vec (magicl:.- position pos))
                                  (distance (sqrt (magicl:tref (magicl:@ (magicl:transpose dist-vec)
                                                                         dist-vec) 0 0))))
                             (>= size distance))))
                       (cl-mpm:sim-mps sim))))
(defun setup-test-column (size &optional (e-scale 1) (mp-scale 1))
  (let* ((sim (cl-mpm/setup::make-block (/ 1 e-scale)
                                        (mapcar (lambda (s) (* s e-scale)) size)
                                        #'cl-mpm::make-shape-function-bspline)) 
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         ;(e-scale 1)
         (h-x (/ h 1))
         (h-y (/ h 1))
         (mass (/ 1 (* e-scale mp-scale)))
         (block-size '(2 5))
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size)))
    (progn
      (setf (cl-mpm:sim-mps sim) 
            (cl-mpm/setup::make-block-mps
             (list (* h-x (- (+ (/ 1 (* 2 mp-scale)) e-scale) 0))
                   (* h-y (+ (/ 1 (* 2 mp-scale)) (* (- 8 (second block-size)) e-scale))))
             block-size
             (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
              'cl-mpm::make-particle-elastic-damage
              1e5 0.2d0
              :mass mass
              :critical-stress 1e2))
      ;; (remove-hole sim '(1d0 5.5d0) 0.5)
      (remove-hole sim '(2d0 5.5d0) 0.30)
      (setf (cl-mpm:sim-damping-factor sim) 0.01d0)
      (setf (cl-mpm:sim-mass-filter sim) 0.01d0)
      (let ((least-pos
              (apply #'min (loop for mp across (cl-mpm:sim-mps sim)
                                 collect (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)))))

        (loop for mp across (cl-mpm:sim-mps sim)
              do
                 (with-accessors ((pos cl-mpm/particle:mp-position)) mp
                   (when (<= (magicl:tref pos 1 0) (+ least-pos 0.001))
                     ;; (setf (cl-mpm/particle:mp-mass mp) (* 50 mass))
                     (setf (magicl:tref  (cl-mpm/particle:mp-body-force mp) 1 0) -100)
                     ))))
      (setf (cl-mpm:sim-dt sim) 1d-4)
      ;; ;(setf (cl-mpm:sim-bcs sim) (cl-mpm/bc:make-outside-bc-nostick (cl-mpm/mesh:mesh-count (cl-mpm:sim-mesh sim))))
      (setf (cl-mpm:sim-bcs sim) (cl-mpm/bc:make-outside-bc (cl-mpm/mesh:mesh-count (cl-mpm:sim-mesh sim))))
      sim)))
(setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))
;Setup
(defun setup ()
  (defparameter *sim* (setup-test-column '(4 8) 8 1))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  (defparameter *sim-step* 0)
  (defparameter *run-sim* nil)
  (defparameter *run-sim* t)
  )

(defparameter *run-sim* nil)

(defun run ()
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
                  (dotimes (i 100)
                    (cl-mpm::update-sim *sim*)
                    (cl-mpm/eigenerosion:update-fracture *sim*)
                    (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))
                    ;; (let ((h (/ (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)) 2)))
                    ;;   (setf *velocity* (cons (magicl:tref (cl-mpm/output::sample-point-velocity *sim* (list h (* h 2))) 1 0) *velocity*)))
                    ;; (setf *time*     (cons *t* *time*))
                    )
                  (plot *sim*)
                  (vgplot:print-plot (asdf:system-relative-pathname "cl-mpm" (format nil "output/frame_~5,'0d.png" steps)))
                  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*))
                                          *sim*)
                  (incf *sim-step*)
                  (swank.live:update-swank)
                  (sleep .01)
      (let ((least-pos
              (apply #'min (loop for mp across (cl-mpm:sim-mps *sim*)
                                 collect (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)))))

        (loop for mp across (cl-mpm:sim-mps *sim*)
              do
                 (with-accessors ((pos cl-mpm/particle:mp-position)) mp
                   (when (<= (magicl:tref pos 1 0) (+ least-pos 0.001))
                     ;; (setf (cl-mpm/particle:mp-mass mp) (* 50 mass))
                     (incf (magicl:tref  (cl-mpm/particle:mp-body-force mp) 1 0) -10)
                     ))))
                  )))
    ;; (vgplot:figure)
    ;; (vgplot:title "Velocity over time")
    ;; (vgplot:plot *time* *velocity*)
    )
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
                      )
  (loop repeat 10
        do (cl-mpm::update-sim *sim*))
  (sb-profile:report))

;; (defparameter *mat* (magicl:zeros '(3 1)))
;; (defparameter *vec* (3d-vectors:vec3 0 0 0))
;; (declaim (ftype (function (3d-vectors:vec3) (values)) zero-3d))
;; (defun zero-3d (x)
;;   (declare (optimize (speed 3)) (type 3d-vectors:vec3 x))
;;   (3d-vectors:nv* x 0)
;;   (values))
;; (declaim (ftype (function (magicl:matrix/double-float) (values)) zero-magicl))
;; (defun zero-magicl (x)
;;   (declare (optimize (speed 3)) (type magicl:matrix/double-float x))
;;   (setf x (magicl:scale x 0d0))
;;   (values))
;; (time-form (zero-3d *vec*) 100000)

;; (time-form (zero-magicl *mat*) 100000)
           

(defmacro time-form (form it)
  `(progn
    (declaim (optimize speed))
    (let* ((iterations ,it)
           (start (get-internal-real-time)))
      (dotimes (i iterations) ,form)
      (let* ((end (get-internal-real-time))
             (units internal-time-units-per-second)
             (dt (/ (- end start) (* iterations units)))
             )
        (format t "Total time: ~f ~%" (/ (- end start) units)) (format t "Time per iteration: ~f~%" (/ (- end start) (* iterations units)))
        (format t "Throughput: ~f~%" (/ 1 dt))
        dt))))
(with-slots ((mesh cl-mpm::mesh)
             (mps  cl-mpm::mps)
             (bcs  cl-mpm::bcs)
             (dt   cl-mpm::dt)
)
            *sim*
  (progn
    (let ((repeats 1000)
          (work (lambda (mesh mp node weight dsvp) (cl-mpm::p2g-mp-node mp node weight dsvp))))
      ;; (time (loop repeat repeats do
      ;;                       (cl-mpm::p2g-mp mesh (aref mps 0))))
      (time-form (cl-mpm::update-sim *sim*) 100)
      ;; (time (loop repeat repeats do
      ;;                       (cl-mpm::p2g-mp-new mesh (aref mps 0))))
      ;; (time (loop repeat repeats do
      ;;   (cl-mpm::iterate-over-neighbours mesh
      ;;                                    (aref mps 0)
      ;;                                    work)))
      ;; (time (loop repeat repeats do
      ;;   (cl-mpm::iterate-over-neighbours-shape-bspline mesh
      ;;                                                 (aref mps 0)
      ;;                                                 work))))
    ;; (time (loop repeat 10000 do
    ;;                       (cl-mpm::iterate-over-neighbours-shape-linear-macro mesh
    ;;                                                        (aref mps 0)
    ;;                                                        '(* weight weight))))
    )))


;; (time-form (cl-mpm::assemble-dsvp 1 '(1)) 100000)
