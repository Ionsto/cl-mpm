(defpackage :cl-mpm/examples/uniaxial
  (:use :cl
        :cl-mpm/example))
(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)
;; (sb-ext:restrict-compiler-policy 'speed  3 3)
;; (sb-ext:restrict-compiler-policy 'debug  0 0)
;; (sb-ext:restrict-compiler-policy 'safety 0 0)
;; (setf *block-compile-default* t)
(in-package :cl-mpm/examples/uniaxial)

(declaim (optimize (debug 0) (safety 0) (speed 3)))

(declaim (notinline plot-domain))
(defun plot-domain (sim)
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :deformed
   :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp))))
(declaim (notinline plot))
(defun plot (sim)
  (vgplot:plot
   (reverse *data-displacement*)
   (reverse *data-load*)
   "MPM"
   )
  ;; (plot-domain sim)
  )
(declaim (notinline setup))
(defun setup (&key (refine 1) (mps 2)
                (epsilon-scale 1d2)
                )
  (let* ((h (/ 1d0 refine))
         (L 10d0)
         (domain-size (list (* L 2) h))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list L h))
         (density 1d3)
         (E 1d9)
         (gf 100d0)
         (init-stress 1d5)
         (length-scale 1d0)
         (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress E))
         ;; (ductility 10d0)
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
      (list 0d0 0d0)
      block-size
      (mapcar (lambda (e) (* (/ e h) mps)) block-size)
      density
      'cl-mpm/particle::particle-elastic-damage
      :E E
      :nu 0d0
      ;; ;;Material parameter
      :initiation-stress init-stress
      :local-length length-scale
      :ductility ductility
      ))
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
      (cl-mpm/penalty::make-bc-penalty-displacment
       *sim*
       (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0))
       (* E epsilon-scale)))
    (cl-mpm:add-bcs-force-list
     *sim*
     *penalty*)
    ))


(defun get-load ()
  (/ (cl-mpm/penalty::resolve-load *penalty*) *h*)
  )

(defun output-disp-header (output-dir)
  (with-open-file (stream (merge-pathnames "disp.csv" output-dir) :direction :output :if-exists :supersede)
    (format stream "disp,load~%")))

(defun output-disp-data (output-dir)
  (with-open-file (stream (merge-pathnames "disp.csv" output-dir) :direction :output :if-exists :append)
    (format stream "~f,~f~%"
            *displacement*
            (get-load))))

(defun test ()
  (defparameter *displacement* 0d0)
  (let* ((lstps 10)
         (total-disp 5d-3)
         (output-dir (format nil "./output-ekl/")))
    (setup :refine 3 :mps 3)
    ;(setf (cl-mpm/damage::sim-enable-length-localisation *sim*) nil)
    (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
    (defparameter *data-displacement* '(0d0))
    (defparameter *data-load* '(0d0))
    (vgplot:close-all-plots)
    (cl-mpm/dynamic-relaxation::run-adaptive-load-control
     *sim*
     :output-dir output-dir
     :plotter (lambda (sim)
                ;; (plot-load-disp)
                (plot *sim*)
                (format t "Load ~E ~%" (cl-mpm/penalty::resolve-load *penalty*)))
     :loading-function (lambda (percent)
                         (setf *displacement* (* total-disp percent))
                         (cl-mpm/penalty::bc-set-displacement
                          *penalty*
                          (cl-mpm/utils:vector-from-list (list *displacement* 0d0 0d0))))
     :pre-step (lambda ()
                 (output-disp-header output-dir)
                 (output-disp-data output-dir))
     :post-conv-step (lambda (sim)
                       ;;Save data
                       (push (get-load) *data-load*)
                       (push *displacement* *data-displacement*)
                       (output-disp-data output-dir)
                       )
     :load-steps lstps
     :enable-damage t
     :damping 1d0;(sqrt 2)
     :substeps 10
     :criteria 1d-3
     :max-adaptive-steps 10
     :save-vtk-dr t
     :save-vtk-loadstep t
     :max-damage-inc 0.3d0
     :dt-scale 1d0)))
