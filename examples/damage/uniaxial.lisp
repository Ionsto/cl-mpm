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
(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-elastic-damage) dt)
  (with-accessors ((y cl-mpm/particle::mp-damage-y-local)
                   (strain cl-mpm/particle::mp-strain)
                   (init-stress cl-mpm/particle::mp-initiation-stress)
                   (ybar cl-mpm/particle::mp-damage-ybar)
                   (de cl-mpm/particle::mp-elastic-matrix)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (E cl-mpm/particle::mp-e)
                   (nu cl-mpm/particle::mp-nu)
                   (k cl-mpm/particle::mp-compression-ratio))
      mp
    (declare (double-float E nu))
    (progn
      (let* ((angle 30d0)
             (stress (cl-mpm/fastmaths:fast-scale!
                      (cl-mpm/constitutive:linear-elastic-mat strain de)
                      1d0
                      ;; (/ 1d0 (magicl:det def))
                      )))
        (setf y
              (cl-mpm/damage::tensile-energy-norm strain E de)
              ;; (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile stress (* angle (/ pi 180d0)))
              ;; (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile stress (* angle (/ pi 180d0)))
              ;; (cl-mpm/damage::criterion-modified-vm strain k E nu)
              )))))
(declaim (notinline setup))
(defun setup (&key (refine 1) (mps 2)
                (epsilon-scale 1d2)
                (gf-scale 1d0)
                )
  (let* ((h (/ 1d0 refine))
         (L 10d0)
         (domain-size (list (* L 2)))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list L))
         (density 1d3)
         (E 1d9)
         (gf 100d0)
         (init-stress 1d5)
         (length-scale 1d0)
         (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress E gf-scale))
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
      (list 0d0)
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

(defun setup-plastic (&key (refine 1) (mps 2)
                        (epsilon-scale 1d2)
                        (gf-scale 1d0))
  (let* ((h (/ 1d0 refine))
         (L 10d0)
         (domain-size (list (* L 2)))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list L))
         (density 1d3)
         (E 1d9)
         (gf 100d0)
         (init-stress 1d5)
         (length-scale 1d0)
         (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress E gf-scale))
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
                                  :enable-fbar t)))
    (defparameter *h* h)
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0d0)
      block-size
      (mapcar (lambda (e) (* (/ e h) mps)) block-size)
      density
      'cl-mpm/particle::particle-vm
      :E E
      :nu 0d0
      :rho (* (sqrt 2/3) init-stress)
      ;; :rho-r 1d0
      ;; :softening 1d0
      ))
    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (when (<= (abs (- (cl-mpm/utils:varef (cl-mpm/particle::mp-position mp) 0) (* 0.5d0 L)))
                 (* h 0.5d0))
         (setf (cl-mpm/particle::mp-rho mp)
               (* 0.9d0 (cl-mpm/particle::mp-rho mp))))))
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
  (cl-mpm/penalty::resolve-load *penalty*))

(defun output-disp-header (output-dir)
  (with-open-file (stream (merge-pathnames "disp.csv" output-dir) :direction :output :if-exists :supersede)
    (format stream "disp,load~%")))

(defun output-disp-data (output-dir)
  (with-open-file (stream (merge-pathnames "disp.csv" output-dir) :direction :output :if-exists :append)
    (format stream "~f,~f~%"
            *displacement*
            (get-load))))

(defun run (&key (output-dir (format nil "./output/")))
  (let ((lstps 50)
        (total-disp 5d-3))
  (defparameter *displacement* 0d0)
    (defparameter *data-displacement* '(0d0))
    (defparameter *data-load* '(0d0))
    (vgplot:close-all-plots)
    (cl-mpm/dynamic-relaxation::run-adaptive-load-control
     *sim*
     :output-dir output-dir
     :plotter
     (lambda (sim)
       (plot *sim*)
       (format t "Load ~E ~%" (cl-mpm/penalty::resolve-load *penalty*)))
     :loading-function
     (lambda (percent)
       (setf *displacement* (* total-disp percent))
       (cl-mpm/penalty::bc-set-displacement
        *penalty*
        (cl-mpm/utils:vector-from-list (list *displacement* 0d0 0d0))))
     :pre-step
     (lambda ()
       (output-disp-header output-dir)
       (output-disp-data output-dir))
     :post-conv-step
     (lambda (sim)
       (push (get-load) *data-load*)
       (push *displacement* *data-displacement*)
       (output-disp-data output-dir))
     :load-steps lstps
     :enable-damage nil
     :enable-plastic t
     :damping 1d0;(sqrt 2)
     :substeps 10
     :criteria 1d-3
     :max-adaptive-steps 10
     :save-vtk-dr t
     :save-vtk-loadstep t
     :max-damage-inc 0.1d0
     :dt-scale 1d0)))

(defun test ()
  (let* ((output-dir (format nil "./output-vm/")))
    (setup-plastic :refine 8 :mps 3)
    ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) nil)
    ;; (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
    (run :output-dir output-dir)
    ))

(defun test-ll ()
  (let* ((refine 4)
         (mps 3))
    (setup :refine refine :mps mps)
    (run :output-dir "./output-standard/")
    (setup :refine refine :mps mps)
    (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
    (run :output-dir "./output-ll/")

    (setup :refine refine :mps mps :gf-scale 1d0)
    (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
    (run :output-dir "./output-ekl/")
    ))
