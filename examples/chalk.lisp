(defpackage :cl-mpm/examples/chalk
  (:use :cl))
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(in-package :cl-mpm/examples/chalk)
(declaim (optimize (debug 3) (safety 3) (speed 0)))

(defun plot (sim)
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :point
   :colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :xy))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-damage-ybar mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-strain-plastic-vm mp))
   ))

(defun setup-test-column (size block-size offset &optional (e-scale 1) (mp-scale 1))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               #'cl-mpm/shape-function:make-shape-function-bspline
               'cl-mpm::mpm-sim-usf
               ;; 'cl-mpm/damage::mpm-sim-damage
               ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 1d3)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
         )
    (declare (double-float h density))
    (progn
      (let ()
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-mps-from-list
               (cl-mpm/setup::make-block-mps-list
                (mapcar (lambda (x) 0) size)
                ;; offset
                ;; '(0 0)
                block-size
                (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                density
                ;; 'cl-mpm/particle::particle-elastic-damage
                'cl-mpm/particle::particle-elastic
                ;; 'cl-mpm/particle::particle-vm
                :E 1d9
                :nu 0.3d0
                ;; :rho 20d3
                ;; :initiation-stress 1d3
                ;; :damage-rate 1d-6
                ;; :critical-damage 0.50d0
                ;; :local-length 2d0
                ;; :local-length-damaged 0.01d0
                ;; :damage 0.0d0
                :gravity -9.8d0
                :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
                ))))
      (setf (cl-mpm:sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (let ((ms 1d0))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        ;; (setf (cl-mpm:sim-damping-factor sim) (* 1d-2 density ms))
        ;; (setf (cl-mpm:sim-damping-factor sim) 10.0d0)
        (setf (cl-mpm:sim-damping-factor sim) 1.0d0)
        )

      (let ((dt-scale 1d0))
        (setf
         (cl-mpm:sim-dt sim)
         (* dt-scale h
            (sqrt (cl-mpm::sim-mass-scale sim))
            (sqrt (/ density (cl-mpm/particle::mp-p-modulus (aref (cl-mpm:sim-mps sim) 0)))))))

      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))
      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-var
             (cl-mpm:sim-mesh sim)
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil 0 nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil 0 nil)))
             ;; (lambda (i) nil)
             ;; (lambda (i) nil)
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
            ))

      ;; (let* ((terminus-size (second block-size))
      ;;        (ocean-y (* terminus-size 1.0d0)))
      ;;   (setf (cl-mpm::sim-bcs-force-list sim)
      ;;         (list
      ;;          (cl-mpm/bc:make-bcs-from-list
      ;;           (list
      ;;            (cl-mpm/buoyancy::make-bc-buoyancy-clip
      ;;             sim
      ;;             ocean-y
      ;;             1000d0
      ;;             (lambda (pos datum)
      ;;               t)
      ;;             ))))))

      sim)))


(defparameter *sim* nil)
(defparameter *run-sim* t)
(defparameter *t* 0)
(defparameter *sim-step* 0)
(defparameter *refine* (/ 1d0 2d0))
(let ((refine (uiop:getenv "REFINE")))
  (when refine
    (setf *refine* (parse-integer (uiop:getenv "REFINE")))  
    ))

(defun setup ()
  ;; (let ((mps-per-dim 4))
  ;;   (defparameter *sim* (setup-test-column '(16 16) '(8 8)  '(0 0) *refine* mps-per-dim)))
  ;; (defparameter *sim* (setup-test-column '(1 1 1) '(1 1 1) 1 1))

  (let* ((mesh-size 50)
         (mps-per-cell 2)
         (shelf-height 400)
         (shelf-aspect 2)
         (shelf-length (* shelf-height shelf-aspect))
         (shelf-height (+ shelf-height))
         (offset (list 0 (* 0 mesh-size)))
         )
    (defparameter *sim*
      (setup-test-column (list (+ shelf-length (* 2 shelf-height))
                               (+ shelf-height 100))
                         (list shelf-length shelf-height)
                         offset
                         (/ 1d0 mesh-size) mps-per-cell))
    (let* ((undercut-angle 10d0)
           (normal (magicl:from-list (list
                                           (cos (- (* pi (/ undercut-angle 180d0))))
                                           (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1))))
      (remove-sdf *sim* (lambda (p) (plane-point-sdf p
                                                     normal
                                                     (magicl:from-list (list shelf-length shelf-height)
                                                                       '(2 1) :type 'double-float)
                                                     ))))
    )
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)

  (let* ((target-time 1d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 1d0)
         )
    (cl-mpm::update-sim *sim*)
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                    (format t "CFL dt estimate: ~f~%" dt-e)
                    (format t "CFL step count estimate: ~D~%" substeps-e)
                    (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
    (time (loop for steps from 0 to 40
                while *run-sim*
                do
                   (progn
                     ;; (when (> steps 5)
                     ;;   (setf (cl-mpm::sim-enable-damage *sim*) t))
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (time
                      (dotimes (i substeps);)
                        (cl-mpm::update-sim *sim*)
                        (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))
                     ;; (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                     ;;   (format t "CFL dt estimate: ~f~%" dt-e)
                     ;;   (format t "CFL step count estimate: ~D~%" substeps-e)
                     ;;   (setf substeps substeps-e))
                     ;; (setf (cl-mpm:sim-damping-factor *sim*)
                     ;;       (* (cl-mpm:sim-damping-factor *sim*) (expt 1d-3 1/40)))

                     (incf *sim-step*)
                     (plot *sim*)
                     (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                        :terminal "png size 1920,1080"
                                        )

                     (swank.live:update-swank)
                     ))))
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*))


;; (setf lparallel:*kernel* (lparallel:make-kernel 8 :name "custom-kernel"))
;; (push (lambda ()
;;         (format t "Closing kernel~%")
;;         (lparallel:end-kernel))
;;       sb-ext:*exit-hooks*)
;; (setup)
;; (run)

;; (time
;;  (dotimes (i 1)
;;    (cl-mpm::update-stress (cl-mpm:sim-mesh *sim*)
;;                           (cl-mpm:sim-mps *sim*)
;;                           (cl-mpm:sim-dt *sim*))))

(defun test ()
  (setup)
  ;; (sb-profile:profile "CL-MPM")
  ;; (sb-profile:profile "CL-MPM/PARTICLE")
  ;; (sb-profile:profile "CL-MPM/MESH")
  ;; (sb-profile:profile "CL-MPM/SHAPE-FUNCTION")
  ;; (sb-profile:reset)
  (time
   (dotimes (i 10)
         (cl-mpm::update-sim *sim*)))
  (sb-profile:report)
  )
;; (lparallel:end-kernel)
;; (sb-ext::exit)
;; (uiop:quit)


;; (declaim (optimize (debug 0) (safety 0) (speed 3)))
;; (defun test ()
;;     (let ((iters 10000000))
;;       (let ((a (cl-mpm/utils:vector-zeros)))
;;         (time
;;          (lparallel:pdotimes (i iters)
;;            (magicl:.+ a (cl-mpm/utils:vector-zeros) a))))
;;       (let ((a (make-array 2 :element-type 'double-float)))
;;         (time
;;          (lparallel:pdotimes (i iters)
;;            (let ((b (make-array 2 :element-type 'double-float)))
;;              (loop for i fixnum from 0 to 1
;;                    do (incf (aref a i) (aref b i))))
;;            )))))
