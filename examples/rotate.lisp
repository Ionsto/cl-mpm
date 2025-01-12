(defpackage :cl-mpm/examples/rotate
  (:use :cl)
  (:import-from
   :cl-mpm/utils varef))
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


(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-elastic) dt fbar)
  ;; (cl-mpm::update-stress-kirchoff-noscale mesh mp dt fbar)
  (cl-mpm::update-stress-kirchoff-noscale mesh mp dt fbar)
  ;; (cl-mpm::co-domain-corner-2d mesh mp dt)
  ;; (cl-mpm::update-domain-deformation mesh mp dt)
  ;; (cl-mpm::update-domain-stretch mesh mp dt)
  (cl-mpm::update-domain-polar-2d mesh mp dt)

  ;; (cl-mpm::update-domain-corner-2d mesh mp dt)
  ;; (cl-mpm::update-domain-max-corner-2d mesh mp dt)
  ;; (cl-mpm::update-domain-midpoint mesh mp dt)
  ;; (cl-mpm::scale-domain-size mesh mp)
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
                ;; (mapcar (lambda (x) 0d0) size)
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
      (setf (cl-mpm::sim-mass-filter sim) 0d0)
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
        (pprint center-point)
        (setf (cl-mpm:sim-bcs sim)
              (cl-mpm/bc::make-domain-bcs
               (cl-mpm:sim-mesh sim)
               (lambda (pos)
                 (let* (
                        (loc
                          (cl-mpm/fastmaths:fast-scale!
                           (cl-mpm/utils:vector-from-list
                            (list (float (nth 0 pos) 0d0)
                                  (float (nth 1 pos) 0d0)
                                  0d0))
                           h-x))
                        (dist (magicl:.-
                               loc
                               (cl-mpm/utils:vector-from-list
                                (list (float (nth 0 center-point) 0d0)
                                      (float (nth 1 center-point) 0d0)
                                      0d0))
                               ))
                        (vel (cl-mpm/penalty::2d-orthog dist)))
                   (cl-mpm/fastmaths:fast-scale! vel 0.1d0)
                   (cl-mpm/bc::make-bc-constant-velocity
                    pos
                    (list (+ (cl-mpm/utils:varef vel 0)
                             ;; (* 0.10d0 (cl-mpm/utils:varef dist 0))
                             )
                          (cl-mpm/utils:varef vel 1)
                          0d0)
                    ))))))
      
      ;; (setup-simple-shear sim)

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
(defun setup-test-stretch (size block-size &optional (e-scale 1) (mp-scale 1))
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
                (mapcar (lambda (x y) (- (* x 0)
                                       (* y 0)
                                       )) size block-size)
                ;; (mapcar (lambda (x) 0d0) size)
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
      (setf (cl-mpm::sim-mass-filter sim) 0d0)
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
      (setup-simple-shear sim)
      sim)))


(defun setup-rotation (sim))

(defun setup-simple-shear (sim)
  (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim))))
    (setf (cl-mpm:sim-bcs sim)
          (cl-mpm/bc::make-domain-bcs
           (cl-mpm:sim-mesh sim)
           (lambda (pos)
             (let ((loc
                     (cl-mpm/fastmaths::fast-scale!
                      (cl-mpm/utils:vector-from-list
                       (list
                        (float (nth 0 pos) 0d0)
                        (float (nth 1 pos) 0d0)
                        (float (nth 2 pos) 0d0)
                        )
                       ) h)))
               (cl-mpm/bc::make-bc-constant-velocity
                pos
                (list (* 0.1d0 (cl-mpm/utils:varef loc 1))
                      0d0               ;(cl-mpm/utils:varef loc 1)
                      0d0)
                )))))))


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
        (domain-size (list 25 25))
        (offset (list 0 0))
        )
    ;(defparameter *sim* (setup-test-column '(16 16 8) '(8 8 8) *refine* mps-per-dim))
    (defparameter *sim* (setup-test-column
                         domain-size
                         (list block-size block-size)
                         refine mps-per-dim))
    )
  ;; (cl-mpm::split-mps-vector
  ;;  *sim*
  ;;  (lambda (mp)
  ;;    (cl-mpm/fastmaths:norm (cl-mpm/utils:vector-from-list (list 1d0 1d0 0d0)))))
  (dotimes (i 1)
    (cl-mpm::split-mps-criteria
     *sim*
     (lambda (mp h)
       :x)))
  ;; (defparameter *sim* (setup-test-column '(1 1 1) '(1 1 1) 1 1))
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))



(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)

  (defparameter *data-step* (list))
  (defparameter *data-volume* (list))
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
    (time (loop for steps from 0 to 100
                while *run-sim*
                do
                   (progn
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     (time
                      (dotimes (i substeps)
                        (cl-mpm::update-sim *sim*)
                        ;; (cl-mpm::split-mps-eigenvalue *sim*)
                        (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))

                     (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       (format t "CFL dt estimate: ~f~%" dt-e)
                       (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf substeps substeps-e))

                     (let ((ar (calculate-raster-area *sim* :sub-res 10)))
                       (format t "Area ratio ~A" ar)
                       (push ar *data-volume*)
                       (push *sim-step* *data-step*))
                     (incf *sim-step*)
                     (plot *sim*)
                     ;; (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                     ;;                    :terminal "png size 1920,1080"
                     ;;                    )

                     (swank.live:update-swank)
                     ))))
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*))



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

(defparameter *raster-array* (make-array 0 :initial-element nil))
(defun calculate-raster-area (sim &key (sub-res 10))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps))
      sim
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (mc (cl-mpm/mesh::mesh-count mesh))
           (raster-grid (if (= (array-total-size *raster-array*)
                               (* (* (nth 0 mc) sub-res) (* (nth 1 mc) sub-res)))
                            (progn (loop for v across (make-array (array-total-size *raster-array*) :displaced-to *raster-array*) do (setf v 0)) *raster-array*)
                            (setf *raster-array*
                                  (make-array (list (* (nth 0 mc) sub-res)
                                                    (* (nth 1 mc) sub-res))
                                              :initial-element nil))))
           (raster-h (/ h sub-res))
          )
      (cl-mpm:iterate-over-mps
       mps
       (lambda (mp)
         (with-accessors ((pos cl-mpm/particle:mp-position)
                          (length cl-mpm/particle::mp-domain-size))
             mp
           (let ((len (cl-mpm/fastmaths:fast-scale-vector length 0.5d0)))
             (loop for x from (round (- (varef pos 0) (varef len 0)) raster-h)
                     upto
                     (round (+ (varef pos 0) (varef len 0)) raster-h)
                   do
                      (loop for y from (round (- (varef pos 1) (varef len 1)) raster-h)
                              upto
                              (round (+ (varef pos 1) (varef len 1)) raster-h)
                            do (setf (aref raster-grid x y) t))
                   )))))
      (let ((total-domain-vol (lparallel:pmap-reduce
                               (lambda (mp)
                                 (* (varef (cl-mpm/particle::mp-domain-size mp) 0)
                                    (varef (cl-mpm/particle::mp-domain-size mp) 1)))
                               #'+
                               mps))
            (raster-domain-vol (* (expt raster-h 2)
                                  (lparallel:pcount-if #'identity (make-array (array-total-size raster-grid) :displaced-to raster-grid)))))
        ;; (pprint total-domain-vol)
        ;; (pprint raster-domain-vol)
        (/ raster-domain-vol total-domain-vol)))))

