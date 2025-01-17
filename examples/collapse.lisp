(defpackage :cl-mpm/examples/collapse
  (:use :cl))
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
;; (setf *block-compile-default* t)
;(sb-int:set-floating-point-modes :traps '(:overflow :invalid :inexact :divide-by-zero :underflow))
;; (sb-int:set-floating-point-modes :traps '(:overflow :divide-by-zero :underflow))

(in-package :cl-mpm/examples/collapse)
(declaim (optimize (debug 3) (safety 3) (speed 0)))

(defun plot-load-disp ()
  (vgplot:plot *data-steps* *data-energy*))
(defun plot (sim)
  ;; (plot-load-disp)
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :deformed
   :colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :xy))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-damage-ybar mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-strain-plastic-vm mp))
   )
  )

(defun stop ()
  (setf *run-sim* nil))

(defun setup-test-column (size block-size &optional (e-scale 1) (mp-scale 1))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               :sim-type 'cl-mpm::mpm-sim-usf
               ;; :sim-type 'cl-mpm::mpm-sim-usl
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
                block-size
                (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                density
                ;; 'cl-mpm/particle::particle-elastic-damage
                ;; 'cl-mpm/particle::particle-elastic
                'cl-mpm/particle::particle-vm
                :E 1d6
                :nu 0.3d0
                :rho 30d3
                ;; :initiation-stress 1d3
                ;; :damage-rate 1d-6
                ;; :critical-damage 0.50d0
                ;; :local-length 2d0
                ;; :local-length-damaged 0.01d0
                ;; :damage 0.0d0
                :gravity -10.0d0
                :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
                ))))
      (setf (cl-mpm:sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-enable-fbar sim) nil)
      (setf (cl-mpm::sim-mass-filter sim) 1d-10)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) t)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) t)
      (let ((ms 1d0))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        (setf (cl-mpm:sim-damping-factor sim)
              (* 1d-2
                 (cl-mpm/setup:estimate-critical-damping sim))))

      (let ((dt-scale 1d0))
        (setf
         (cl-mpm:sim-dt sim)
         (* dt-scale h
            (sqrt (cl-mpm::sim-mass-scale sim))
            (sqrt (/ density (cl-mpm/particle::mp-p-modulus (aref (cl-mpm:sim-mps sim) 0)))))))

      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))
      ;; (setf (cl-mpm:sim-bcs sim)
      ;;       (cl-mpm/bc::make-outside-bc-var
      ;;        (cl-mpm:sim-mesh sim)
      ;;        (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
      ;;        (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
      ;;        (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil 0 nil)))
      ;;        (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil 0 nil)))
      ;;        ;; (lambda (i) nil)
      ;;        ;; (lambda (i) nil)
      ;;        (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
      ;;        (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
      ;;       ))
      (setf
       (cl-mpm:sim-bcs sim)
       (cl-mpm/bc::make-outside-bc-varfix
        (cl-mpm:sim-mesh sim)
        '(0 nil 0)
        '(0 nil 0)
        '(nil 0 0)
        '(nil 0 0)
        '(nil nil 0)
        '(nil nil 0)))

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
(defparameter *refine* (/ 1d0 1d0))
(let ((refine (uiop:getenv "REFINE")))
  (when refine
    (setf *refine* (parse-integer (uiop:getenv "REFINE")))  
    ))

(defun setup (&key (refine 1) (mps 2))
  (let ((mps-per-dim mps))
    ;(defparameter *sim* (setup-test-column '(16 16 8) '(8 8 8) *refine* mps-per-dim))
    (defparameter *sim* (setup-test-column '(16 16) '(8 8) refine mps-per-dim))
    )
  ;; (defparameter *sim* (setup-test-column '(1 1 1) '(1 1 1) 1 1))
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)

  (defparameter *data-steps* (list))
  (defparameter *data-oobf* (list))
  (defparameter *data-energy* (list))
  (let* ((target-time 0.5d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 0.5d0)
         (dt-min (cl-mpm:sim-dt *sim*))
         )
    (cl-mpm::update-sim *sim*)
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                    (format t "CFL dt estimate: ~f~%" dt-e)
                    (format t "CFL step count estimate: ~D~%" substeps-e)
                    (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
    (time (loop for steps from 0 to 100
                while *run-sim*
                do
                   (progn
                     ;; (when t;(> steps 20)
                     ;;   (setf (cl-mpm::sim-damping-factor *sim*) 1d0)
                     ;;   ;; (setf (cl-mpm::sim-enable-damage *sim*) t)
                     ;;   )
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     (setf dt-min (cl-mpm::calculate-min-dt *sim*))
                     (time
                      (dotimes (i substeps)
                        (cl-mpm::update-sim *sim*)
                        (setf dt-min (min (cl-mpm::calculate-min-dt *sim*) dt-min))
                        (setf (cl-mpm:sim-dt *sim*) (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
                        (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))

                     (format t "Min dt est ~f~%" dt-min)
                     ;; (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                     ;;   (format t "CFL dt estimate: ~f~%" dt-e)
                     ;;   (format t "CFL step count estimate: ~D~%" substeps-e)
                     ;;   (setf substeps substeps-e))
                     (setf (cl-mpm:sim-dt *sim*) (* dt-min dt-scale))
                     (setf substeps (floor target-time (cl-mpm:sim-dt *sim*)))
                     (format t "Substeps ~D~%" substeps)
                     (let ((energy 0d0))
                       (setf energy (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
                       (push energy *data-energy*)
                       (push steps *data-steps*)
                       (format t "Energy ~F~%" energy))

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


(defmacro time-form (it form)
  `(progn
     (declaim (optimize speed))
     (let* ((iterations ,it)
            (start (get-internal-real-time)))
       (time
        (dotimes (i ,it)
              ,form))
       (let* ((end (get-internal-real-time))
              (units internal-time-units-per-second)
              (dt (/ (- end start) (* iterations units)))
              )
         (format t "Total time: ~f ~%" (/ (- end start) units)) (format t "Time per iteration: ~f~%" (/ (- end start) (* iterations units)))
         (format t "Throughput: ~f~%" (/ 1 dt))
         (format t "Time per MP: ~E~%" (/ dt (length (cl-mpm:sim-mps *sim*))))
         dt))))
(defun test-mult ()
  (time
   (cl-mpm:iterate-over-mps
    (cl-mpm:sim-mps *sim*)
    (lambda (mp)
      (with-accessors ((strain cl-mpm/particle:mp-strain)
                       (de cl-mpm/particle::mp-elastic-matrix)
                       (stress cl-mpm/particle::mp-stress-kirchoff))
          mp
        (cl-mpm/constitutive::linear-elastic-mat strain de stress))))))
(defun profile ()
  (setup :refine 16)
  ;; (sb-profile:unprofile)
  ;; (sb-profile:profile "CL-MPM")
  ;; ;; (sb-profile:profile "CL-MPM/PARTICLE")
  ;; ;; (sb-profile:profile "CL-MPM/MESH")
  ;; ;; (sb-profile:profile "CL-MPM/SHAPE-FUNCTION")
  ;; (sb-profile:reset)
  (time-form 100
             (progn
               (format t "~D~%" i)
                  (cl-mpm::update-sim *sim*)))
  ;; (time-form 1000000
  ;;            (progn
  ;;              ;; (cl-mpm:iterate-over-neighbours (cl-mpm:sim-mesh *sim*) (aref (cl-mpm:sim-mps *sim*) 0) (lambda (mesh mp &rest args) (setf (fill-pointer (cl-mpm/particle::mp-cached-nodes mp)) 0)))
  ;;              (cl-mpm::iterate-over-neighbours-shape-gimp-simd (cl-mpm:sim-mesh *sim*) (aref (cl-mpm:sim-mps *sim*) 0) (lambda (&rest args)))
  ;;              ))
  ;; (time
  ;;  (dotimes (i 100)
  ;;    (format t "~D~%" i)
  ;;    (cl-mpm::update-sim *sim*)))
  (format t "MPS ~D~%" (length (cl-mpm:sim-mps *sim*)))
  ;; (sb-profile:report)
  )

(defun profile ()
  (setup :refine 16)
  (time-form 100
             (progn
               (format t "~D~%" i)
               (cl-mpm::update-sim *sim*)))
  ;; (time-form
  ;;  100
  ;;  (progn
  ;;    (cl-mpm::check-mps *sim*)))
  ;; (time
  ;;  (dotimes (i 100)
  ;;    (format t "~D~%" i)
  ;;    (cl-mpm::update-sim *sim*)))
  (format t "MPS ~D~%" (length (cl-mpm:sim-mps *sim*)))
  ;; (sb-profile:report)
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


;; (cl-mpm::iterate-over-nodes-serial
;;  (cl-mpm:sim-mesh *sim*)
;;  (lambda (node)
;;    (loop for v across (magicl::storage (cl-mpm/mesh::node-velocity node))
;;          do
;;             (when (sb-ext::float-nan-p v)
;;               (format t "Found nan vel~%")
;;               (pprint node)
;;               ))))

(defun find-nans ()
  (cl-mpm::iterate-over-nodes
   (cl-mpm:sim-mesh *sim*)
   (lambda (mp)
     (loop for v across (magicl::storage (cl-mpm/mesh::node-velocity mp))
           do
              ;; (when (> v 0d0)
              ;;   (format t "~E~%" v))
              (when (or (sb-ext::float-nan-p v)
                        ;; (> v 1d10)
                        ;; (< v -1d10)
                        )
                (format t "Found nan vel~%")
                (pprint mp)
                )
           ))))



;;Old
;; Evaluation took:
;; 7.969 seconds of real time
;; 71.424494 seconds of total run time (53.687638 user, 17.736856 system)
;; [ Real times consist of 2.410 seconds GC time, and 5.559 seconds non-GC time. ]
;; [ Run times consist of 2.467 seconds GC time, and 68.958 seconds non-GC time. ]
;; 896.27% CPU
;; 33,467,579,605 processor cycles
;; 22,937,955,888 bytes consed

;; Total time: 7.969999 
;; Time per iteration: 0.07969999
;; Throughput: 12.547053
;; Time per MP: 4.8645017e-6
;; MPS 16384

;;New


(defun stop ()
  (setf *run-sim* nil)
  (setf cl-mpm/dynamic-relaxation::*run-convergance* nil))
(defun test-conv ()
  (defparameter *data-steps* (list))
  (defparameter *data-oobf* (list))
  (defparameter *data-energy* (list))
  (defparameter *data-name* (list))

  (vgplot:figure)
  (dolist (algo (list :FLIP :PIC :BLEND))
    (dolist (dt-scale (list 0.5d0))
      (let* ((steps (list))
             (energy (list))
             (oobf (list))
             ;; (dt-scale 0.5d0)
             (substeps (round 10 dt-scale))
             (path (merge-pathnames "./analysis_scripts/vel_algo/data/"))
             )
        (setup :refine 1 :mps 3)
        (setf (cl-mpm:sim-damping-factor *sim*) 
              (* 1d-1 (cl-mpm/setup::estimate-critical-damping *sim*)))
        (setf
         (cl-mpm::sim-velocity-algorithm *sim*)
         algo)
        (let ((name (format nil "~A_~E"
                            (cl-mpm::sim-velocity-algorithm *sim*)
                            dt-scale)))
          (vgplot:title name)
          (cl-mpm/dynamic-relaxation:converge-quasi-static
           *sim*
           :dt-scale dt-scale
           :energy-crit 1d-3
           :oobf-crit 1d-3
           :substeps substeps
           :conv-steps 1000
           :dt-scale dt-scale
           :pic-update nil
           :post-iter-step
           (lambda (i e o)
             ;; (plot *sim*)
             (push i steps)
             (push e energy)
             (push o oobf)
             ;; (vgplot:semilogy
             ;;  steps energy "energy"
             ;;  steps oobf "oobf")
             (apply #'vgplot:semilogy
                    (reduce #'append
                            (mapcar #'list
                                    (append *data-steps* (list steps))
                                    (append *data-oobf*  (list oobf))
                                    (append *data-name*  (list name))
                                    ))))
           )
          (with-open-file (stream (merge-pathnames path (format nil "~A.csv" name)) :direction :output :if-exists :supersede)
            (format stream "step,energy,oobf~%")
            (loop for s in steps
                  for e in energy
                  for o in oobf
                  do (format stream "~D,~F,~F~%" s e o)))
          (cl-mpm/output:save-vtk (merge-pathnames path (format nil "sim_~A.vtk" name)) *sim*)
          (cl-mpm/output:save-vtk-nodes (merge-pathnames path (format nil "sim_nodes_~A.vtk" name)) *sim*)

          (push steps *data-steps*)
          (push energy *data-energy*)
          (push oobf *data-oobf*)
          (push name *data-name*)
          )))))

(defun plot-all ()
  (apply #'vgplot:semilogy
         (reduce #'append
                 (mapcar #'list
                         *data-steps*
                         *data-oobf*
                         *data-name*
                         ))))

(defun plot-some ()
  (let* ((indexes (list 0 1 2
                        9 10 11
                        ))
         (steps (loop for i in indexes collect (nth i *data-steps*)))
         (energy (loop for i in indexes collect (nth i *data-energy*)))
         (oobf (loop for i in indexes collect (nth i *data-oobf*)))
         (name (loop for i in indexes collect (nth i *data-name*)))
        )
    (apply #'vgplot:semilogy
           (reduce #'append
                   (mapcar #'list
                           steps
                           oobf
                           name
                           )))))

;; (let ((mesh (cl-mpm:sim-mesh *sim*))
;;       (mp0 (aref (cl-mpm:sim-mps *sim*) 0)))
;;   (pprint (cl-mpm/particle:mp-position mp0))
;;   (cl-mpm::iterate-over-neighbours
;;    mesh mp0
;;    (lambda (mesh mp node svp grads fsvp fgrads)
;;      (with-accessors ((node-active cl-mpm/mesh:node-active)
;;                       (node-j-inc cl-mpm/mesh::node-jacobian-inc)
;;                       (node-volume cl-mpm/mesh::node-volume))
;;          node
;;        (pprint svp)))))
