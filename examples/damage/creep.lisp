(defpackage :cl-mpm/examples/damage/creep
  (:use :cl
        :cl-mpm/example
        ))
;; (sb-ext:restrict-compiler-policy 'speed  3 3)
;; (sb-ext:restrict-compiler-policy 'debug  0 0)
;; (sb-ext:restrict-compiler-policy 'safety 0 0)
;; (setf *block-compile-default* nil)
(in-package :cl-mpm/examples/damage/creep)
(declaim (optimize (debug 3) (safety 3) (speed 0)))

(defparameter *data-time* '())
(defparameter *data-disp* '())
(defparameter *data-load* '())
(defparameter *penalty* nil)


(declaim (notinline plot-domain))
(defun plot-domain (sim)
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :deformed
   :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp))))

(declaim (notinline plot))
(defun plot (sim)
  (when (> (length *data-disp*) 1)
    (vgplot:plot
     (reverse *data-time*)
     (reverse *data-disp*)
     "MPM")))

(declaim (notinline setup))
(defun setup (&key (refine 1) (mps 2)
                (epsilon-scale 1d2)
                (delay-time 1d0)
                (delay-exp 1d0)
                (gf 100d0)
                (gf-scale 1d0))
  (let* ((h (/ 1d0 refine))
         (L 10d0)
         (domain-size (list (* L 2)))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list L))
         (density 1d3)
         (E 1d9)
         (init-stress 1d5)
         (length-scale 1d0)
         (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 (* gf-scale gf) length-scale init-stress E))
         ;; (ductility 50d0)
         )
    (when (<= ductility 1d0)
      (error "Ductility must be greater than 1"))
    (format t "Ductility ~F~%" ductility)
    (setf *sim* (cl-mpm/setup::make-simple-sim
                 h
                 element-count
                 :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul
                 :args-list (list :enable-aggregate t
                                  :enable-fbar nil
                                  ;; :ghost-factor (* E 1d-3)
                                  )))
    (defparameter *h* h)
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0d0)
      block-size
      (mapcar (lambda (e) (* (/ e h) mps)) block-size)
      density
      'cl-mpm/particle::particle-elastic-damage-delayed
      :E E
      :nu 0d0
      ;; ;;Material parameter
      :initiation-stress init-stress
      :local-length length-scale
      :ductility ductility
      :residual-strength (- 1d0 1d-9)
      :delay-time delay-time
      :delay-exponent delay-exp))
    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (when (<= (abs (- (cl-mpm/utils:varef (cl-mpm/particle::mp-position mp) 0) (* 0.5d0 L)))
                 (* h 0.5d0))
         (setf (cl-mpm/particle::mp-initiation-stress mp)
               (* 0.9d0 (cl-mpm/particle::mp-initiation-stress mp))))))
    (setf (cl-mpm:sim-gravity *sim*) 0d0)
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)
    (cl-mpm/setup::setup-bcs
     *sim*
     ;:left '(0 nil nil)
     )
    (defparameter *penalty*
      (cl-mpm/buoyancy::make-bc-pressure
       *sim*
       2d5
       0d0))
    ;; (defparameter *penalty*
    ;;   (cl-mpm/penalty::make-bc-penalty-displacment
    ;;    *sim*
    ;;    (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0))
    ;;    (* E epsilon-scale)))
    (cl-mpm:add-bcs-force-list
     *sim*
     *penalty*)
    ))


(defun get-disp ()
  (- (cl-mpm::reduce-over-mps
      (cl-mpm:sim-mps *sim*)
      (lambda (mp)
        (cl-mpm/utils:varef (cl-mpm/fastmaths::fast-.+
                             (cl-mpm/particle::mp-position mp)
                             (cl-mpm/particle::mp-domain-size mp))
                            0))
      #'max)
     10d0))
(defun get-load ()
  (first (cl-mpm/buoyancy::bc-pressure-pressures *penalty*)))

(defun output-disp-header (output-dir)
  (with-open-file (stream (merge-pathnames "disp.csv" output-dir) :direction :output :if-exists :supersede)
    (format stream "time,disp,load~%")))

(defun output-disp-data (output-dir)
  (with-open-file (stream (merge-pathnames "disp.csv" output-dir) :direction :output :if-exists :append)
    (format stream "~f,~f~%"
            (cl-mpm:sim-time *sim*)
            (get-disp)
            (get-load)
            )))

(defun run (&key (output-dir (format nil "./output/")))
  (let ()
    (defparameter *data-time* '())
    (defparameter *data-disp* '())
    (defparameter *data-load* '())
    (vgplot:close-all-plots)
    (ensure-directories-exist (merge-pathnames output-dir))
    (cl-mpm/dynamic-relaxation::run-quasi-time
     *sim*
     :output-dir output-dir
     :dt 0.05d0
     :total-time 10d0
     :plotter
     (lambda (sim)
       (plot *sim*)
       (format t "Load ~E ~%" (get-load)))
     :post-conv-step
     (lambda (sim)
       (output-disp-header output-dir)
       (output-disp-data output-dir))
     :post-load-step
     (lambda (sim)
       (push (get-load) *data-load*)
       (push (get-disp) *data-disp*)
       (push (cl-mpm:sim-time *sim*) *data-time*)
       (output-disp-data output-dir))
     :substeps 50
     :max-adaptive-steps 5
     :min-adaptive-steps 0
     :save-vtk-dr t
     :save-vtk-loadstep t
     ;; :max-damage-inc 0.9d0
     ;; :true-stagger nil
     :dt-scale 1d0)))



(defun test ()
  (cl-mpm/utils:set-workers 8)
  (let* ((output-dir (format nil "./output/")))
    (setup :refine 4 :mps 4
           :gf 100d0
           :delay-time 1d0
           :delay-exp 1d0
           )
    ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
    ;; (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
    (run :output-dir output-dir)
    ))


