(defpackage :cl-mpm/examples/rotate
  (:use :cl))
(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)
;; (sb-ext:restrict-compiler-policy 'speed  3 3)
;; (sb-ext:restrict-compiler-policy 'debug  0 0)
;; (sb-ext:restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* nil)
;(sb-int:set-floating-point-modes :traps '(:overflow :invalid :inexact :divide-by-zero :underflow))
;; (sb-int:set-floating-point-modes :traps '(:overflow :divide-by-zero :underflow))

(in-package :cl-mpm/examples/rotate)
(declaim (optimize (debug 3) (safety 3) (speed 0)))
(let* ((D (cl-mpm/utils:matrix-from-list
           (list 0d0 0.1d0 0d0
                 -0.1d0 0d0 0d0
                 0d0 0d0 0d0)))
       (dt 1d0)
       (domain (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0)))
       (omega (magicl:scale
                   (magicl:.-
                    D
                    (magicl:transpose D)) (* dt 0.5d0)))
           (R (magicl:.+
               (magicl:eye 3)
               omega))
           (dom (cl-mpm/utils:matrix-from-list
                 (list
                  (cl-mpm/utils:varef domain 0) 0d0 0d0
                  0d0 (cl-mpm/utils:varef domain 1) 0d0
                  0d0 0d0 (cl-mpm/utils:varef domain 2))))
           (DR (magicl:@
                R
                dom
                (magicl:transpose R)
                ))
           )
      (pprint DR)
      (setf (cl-mpm/utils:varef domain 0) (magicl:tref DR 0 0))
      (setf (cl-mpm/utils:varef domain 1) (magicl:tref DR 1 1))
      ;; (setf (cl-mpm/utils:varef domain 1) (magicl:tref dom 1 1))

      ;; (pprint D)
      ;; (pprint omega)
      ;; (setf domain
      ;;       (magicl:@
      ;;         domain))
      )

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-elastic) dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-noscale mesh mp dt fbar)
  (cl-mpm::update-stress-kirchoff mesh mp dt fbar)
  (cl-mpm::co-domain-corner-2d mesh mp dt)
  ;; (cl-mpm::update-domain-corner-2d mesh mp dt)
  ;; (with-accessors ((domain cl-mpm/particle::mp-domain-size)
  ;;                  (D cl-mpm/particle::mp-stretch-tensor))
  ;;     mp
  ;;   (let* ((omega (magicl:scale
  ;;                  (magicl:.-
  ;;                   D
  ;;                   (magicl:transpose D)) (* dt 0.5d0)))
  ;;          (R (magicl:.+
  ;;              (magicl:eye 3)
  ;;              omega))
  ;;          (dom (cl-mpm/utils:matrix-from-list
  ;;                (list
  ;;                 (cl-mpm/utils:varef domain 0) 0d0 0d0
  ;;                 0d0 (cl-mpm/utils:varef domain 1) 0d0
  ;;                 0d0 0d0 (cl-mpm/utils:varef domain 2))))
  ;;          (DR (magicl:@
  ;;               R
  ;;               dom
  ;;               (magicl:transpose R)
  ;;               ))
  ;;          ;; (DR (magicl:.+
  ;;          ;;      (magicl:@
  ;;          ;;       dom
  ;;          ;;       R
  ;;          ;;       )
  ;;          ;;      (magicl:@
  ;;          ;;       R
  ;;          ;;       dom)
  ;;          ;;      ))
  ;;          )
  ;;     ;; (pprint DR)
  ;;     (setf (cl-mpm/utils:varef domain 0) (magicl:tref DR 0 0))
  ;;     (setf (cl-mpm/utils:varef domain 1) (magicl:tref DR 1 1))
  ;;     ;; (setf (cl-mpm/utils:varef domain 1) (magicl:tref dom 1 1))

  ;;     ;; (pprint D)
  ;;     ;; (pprint omega)
  ;;     ;; (setf domain
  ;;     ;;       (magicl:@
  ;;     ;;         domain))
  ;;     ))
  ;; (cl-mpm::scale-domain-size mesh mp)
  ;; (cl-mpm::update-stress-linear mesh mp dt fbar)
  )

(declaim (notinline plot))
(defun plot (sim)
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

(declaim (notinline setup-test-column))
(defun setup-test-column (size block-size &optional (e-scale 1) (mp-scale 1))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               :sim-type 'cl-mpm::mpm-sim-usf
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
                (mapcar (lambda (x y) (- (* x 0.5)
                                       (* y 0.5)
                                       )) size block-size)
                block-size
                (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                density
                ;; 'cl-mpm/particle::particle-elastic-damage
                'cl-mpm/particle::particle-elastic
                ;; 'cl-mpm/particle::particle-vm
                :E 1d6
                :nu 0.3d0
                ;; :rho 20d3
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
      ;; (setf (cl-mpm::sim-mass-filter sim) 1d0)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) t)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) t)
      (let ((ms 1d0))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        (setf (cl-mpm:sim-damping-factor sim) (* 1d-2 density ms))
        )

      (let ((dt-scale 1d0))
        (setf
         (cl-mpm:sim-dt sim)
         (* dt-scale h
            (sqrt (cl-mpm::sim-mass-scale sim))
            (sqrt (/ density (cl-mpm/particle::mp-p-modulus (aref (cl-mpm:sim-mps sim) 0)))))))

      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))
      (let ((center-point
              (mapcar (lambda (x) (* x 0.5)) size)))
        (setf (cl-mpm:sim-bcs sim)
              (cl-mpm/bc::make-domain-bcs
               (cl-mpm:sim-mesh sim)
               (lambda (pos)
                 (let* ((dist (magicl:.-
                               (cl-mpm/utils:vector-from-list
                                (list (float (nth 0 pos) 0d0)
                                      (float (nth 1 pos) 0d0)
                                      (float (nth 2 pos) 0d0))
                                )
                               (cl-mpm/utils:vector-from-list
                                (list (float (nth 0 center-point) 0d0)
                                      (float (nth 1 center-point) 0d0)
                                      0d0))))
                        (vel (cl-mpm/penalty::2d-orthog dist)))
                   (cl-mpm/fastmaths:fast-scale! vel 0.1d0)
                   (cl-mpm/bc::make-bc-constant-velocity
                    pos
                    (list (cl-mpm/utils:varef vel 0)
                          (cl-mpm/utils:varef vel 1)
                          0d0)
                    ))))
              ;; (cl-mpm/bc::make-outside-bc-var
              ;;  (cl-mpm:sim-mesh sim)
              ;;  (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
              ;;  (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
              ;;  (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil 0 nil)))
              ;;  (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil 0 nil)))
              ;;  ;; (lambda (i) nil)
              ;;  ;; (lambda (i) nil)
              ;;  (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
              ;;  (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
              ;; )
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
(defparameter *refine* (/ 1d0 1d0))
(let ((refine (uiop:getenv "REFINE")))
  (when refine
    (setf *refine* (parse-integer (uiop:getenv "REFINE")))  
    ))

(defun setup (&key (refine 1) (mps 3))
  (let ((mps-per-dim mps)
        (block-size 5)
        )
    ;(defparameter *sim* (setup-test-column '(16 16 8) '(8 8 8) *refine* mps-per-dim))
    (defparameter *sim* (setup-test-column
                         '(10 10) (list block-size block-size)
                         refine mps-per-dim))
    )
  (cl-mpm::split-mps-criteria
   *sim*
   (lambda (mp h)
     :x))
  ;; (defparameter *sim* (setup-test-column '(1 1 1) '(1 1 1) 1 1))
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))


;; (defun set-velocity ()
;;   (cl-mpm:iterate-over-nodes

;;    )
;;   )

(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)

  (let* ((target-time 0.5d0)
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
    (time (loop for steps from 0 to 500
                while *run-sim*
                do
                   (progn
                     (when (> steps 5)
                       (setf (cl-mpm::sim-enable-damage *sim*) t))
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     (time
                      (dotimes (i substeps)
                        (cl-mpm::update-sim *sim*)
                        (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))

                     (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       (format t "CFL dt estimate: ~f~%" dt-e)
                       (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf substeps substeps-e))

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