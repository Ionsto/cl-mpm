(declaim (optimize (debug 3) (safety 0) (speed 0)))
;(declaim (optimize (debug 3) (safety 3) (speed 0)))
(ql:quickload "cl-mpm")
(ql:quickload "cl-mpm/setup")
(ql:quickload "cl-mpm/particle")
(ql:quickload "cl-mpm/bc")
(ql:quickload "vgplot")
(ql:quickload "swank.live")
(ql:quickload "cl-mpm/output")
(ql:quickload "magicl")
;; (ql:quickload "py4cl")
;; (setf py4cl:*python-command* "python3")
;; (py4cl:import-module "matplotlib.pyplot" :as "plt")
(ql:quickload "py4cl2")
(py4cl2:defpymodule "matplotlib" nil :lisp-package "MPL")
(mpl:use :backend "TKAgg")
(py4cl2:defpymodule "matplotlib.pyplot" nil :lisp-package "PLT")
;; (plt:figure)
;; (plt:ion)
;; (plt:ioff)
;; (plt:plot '(1 2 3))
;; (plt:show)
(defun plot (sim)
  (multiple-value-bind (x y)
    (loop for mp across (cl-mpm:sim-mps sim)
          collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
          collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
          finally (return (values x y)))
    (vgplot:plot x y ";;with points pt 7"))
  (vgplot:replot))
(defun plot-pyplot (sim)
  (multiple-value-bind (x y c)
    (loop for mp across (cl-mpm:sim-mps sim)
          collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
          collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
          collect (cl-mpm/particle:mp-damage mp) into c
          finally (return (values x y c)))
    (plt:scatter :x x :y y :c c))
  (plt:xlim 0 (first (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim))))
  (plt:ylim 0 (second (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim))))
  (plt:clim 0 1)
  ;; (plt:colorbar)
  (plt:show))

(defun setup-test-column (size &optional (e-scale 1))
  (let* ((sim (cl-mpm/setup::make-block (/ 1 e-scale)
                                        (mapcar (lambda (s) (* s e-scale)) size)
                                        #'cl-mpm::make-shape-function-linear)) 
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         ;(e-scale 1)
         (h-x (/ h 1))
         (h-y (/ h 1))
         (block-size '(2 5))
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size)))
    (progn
      (setf (cl-mpm:sim-mps sim) 
            (cl-mpm/setup::make-block-mps
             (list (* h-x (+ 0.5 e-scale))
                   (* h-y (+ 0.5 (* (- 8 (second block-size)) e-scale))))
             block-size
             (mapcar (lambda (e) (* e e-scale)) block-size)
              'cl-mpm::make-particle-elastic-damage
              1e5 0d0 :mass 5))
      ;; (setf (cl-mpm:sim-mps sim) 
      ;;       (cl-mpm/setup::make-block-mps
      ;;        (list (* 1.5 h-x e-scale)
      ;;              (* 5   h-y e-scale))
      ;;         (mapcar (lambda (s) (* s e-scale)) (list 2 2)) 
      ;;         (mapcar (lambda (s) (* s 1)) '(2 2))
      ;;         'cl-mpm::make-particle-elastic-damage
      ;;         1e5 0d0 :mass 1))
      ;; (loop for mp across (cl-mpm:sim-mps sim) 
      ;;       do (progn
      ;;            (with-accessors ((pos cl-mpm/particle::mp-position)) mp
      ;;              (setf (magicl:tref pos 1 0) (+ 0.5d0 (magicl:tref pos 1 0))))))
      (setf (cl-mpm:sim-damping-factor sim) 0d0)
      (setf (cl-mpm:sim-mass-filter sim) 0.01d0)
      (loop for i from 1 to e-scale do
        (setf (cl-mpm/particle:mp-mass (aref (cl-mpm:sim-mps sim) (- (length (cl-mpm:sim-mps sim)) i))) 1000))
      (setf (cl-mpm:sim-dt sim) 1e-3)
      ;; ;(setf (cl-mpm:sim-bcs sim) (cl-mpm/bc:make-outside-bc-nostick (cl-mpm/mesh:mesh-count (cl-mpm:sim-mesh sim))))
      (setf (cl-mpm:sim-bcs sim) (cl-mpm/bc:make-outside-bc (cl-mpm/mesh:mesh-count (cl-mpm:sim-mesh sim))))
      sim)))
(setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))
;Setup
(progn
  (defparameter *sim* (setup-test-column '(4 8) 2))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  )


(defparameter *run-sim* nil)
(defparameter *run-sim* t)

(progn 
    (vgplot:close-all-plots)
    (vgplot:figure)
    (vgplot:axis (list 0 (nth 0 (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*))) 
                       0 (nth 1 (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))))
    (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
      (vgplot:format-plot t "set ytics ~f" h)
      (vgplot:format-plot t "set xtics ~f" h))
    (time (loop for steps from 0 to 10
                while *run-sim*
                do
                (progn
                  (format t "Step ~d ~%" steps)
                  (dotimes (i 100)
                    (cl-mpm::update-sim *sim*)
                    (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))
                    (let ((h (/ (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)) 2)))
                      (setf *velocity* (cons (magicl:tref (cl-mpm/output::sample-point-velocity *sim* (list h (* h 2))) 1 0) *velocity*)))
                    (setf *time*     (cons *t* *time*))
                    )
                  (plot *sim*)
                  (vgplot:print-plot (asdf:system-relative-pathname "cl-mpm" (format nil "output/frame_~5,'0d.png" steps)))
                  (swank.live:update-swank)
                  (sleep .02)
                  )))
    (vgplot:figure)
    (vgplot:title "Velocity over time")
    (vgplot:plot *time* *velocity*)
    )