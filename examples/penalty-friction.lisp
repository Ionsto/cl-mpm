(defpackage :cl-mpm/examples/penalty-friction
  (:use :cl))
(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)
;; (sb-ext:restrict-compiler-policy 'speed  3 3)
;; (sb-ext:restrict-compiler-policy 'debug  0 0)
;; (sb-ext:restrict-compiler-policy 'safety 0 0)
;(sb-int:set-floating-point-modes :traps '(:overflow :invalid :inexact :divide-by-zero :underflow))
;; (sb-int:set-floating-point-modes :traps '(:overflow :divide-by-zero :underflow))

;; (setf *block-compile-default* t)
(in-package :cl-mpm/examples/penalty-friction)
(declaim (optimize (debug 3) (safety 3) (speed 0)))

(defun plot (sim)
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :point
   ;:colour-func (lambda (mp) (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0))
   :colour-func (lambda (mp) (magicl:tref (cl-mpm/particle::mp-penalty-frictional-force mp) 0 0))
   ))

(defun setup-test-column (size block-size &optional (e-scale 1) (mp-scale 1) &key (angle 0d0) (friction 0.1d0))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               :sim-type 'cl-mpm::mpm-sim-usl
               ;; 'cl-mpm/damage::mpm-sim-damage
               ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         ;; (floor-offset (* h-y 2))
         (floor-offset 2d0)
         (density 1d3)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
         )
    (declare (double-float h density))
    (progn
      (let ((angle-rad (* angle (/ pi 180))))
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-mps-from-list
               (cl-mpm/setup::make-block-mps-list
                ;; (mapcar (lambda (x) 0) size)
                (list
                 1d0
                 floor-offset
                      0)
                block-size
                (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                density
                'cl-mpm/particle::particle-elastic
                :E 1d9
                :nu 0.30d0
                :gravity -10.0d0
                :gravity-axis (cl-mpm/utils:vector-from-list (list
                                                              (- (sin angle-rad))
                                                              (cos angle-rad)
                                                              0d0))
                ))))
      (setf (cl-mpm:sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      ;; (setf (cl-mpm::sim-mass-filter sim) 1d0)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) t)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) t)
      (let ((ms 1d0))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        (setf (cl-mpm:sim-damping-factor sim)
              (* 0d0 (cl-mpm/setup::estimate-critical-damping sim))))

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
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
            ))
      (defparameter *floor-bc*
        (cl-mpm/penalty::make-bc-penalty-point-normal
         sim
         (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
         (cl-mpm/utils:vector-from-list (list 0d0 floor-offset 0d0))
         (* 1d9 1d0)
         friction)
        )
      (setf (cl-mpm::sim-bcs-force-list sim)
            (list
             (cl-mpm/bc:make-bcs-from-list
              (list
               *floor-bc*))))
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

(defun setup (&key (refine 1d0) (mps 2) (angle 0d0) (friction 0.3d0))
  (let ((mps-per-dim mps))
    (defparameter *sim* (setup-test-column
                         '(16 8)
                         '(2 2)
                         (* 1d0 refine)
                         mps-per-dim
                         :angle angle
                         :friction friction
                         )))
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (format t "Mesh-size: ~E~%" (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
  ;; (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)

  (defparameter *data-t* (list))
  (defparameter *data-v* (list))
  (let* ((target-time 0.1d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 1.0d0)
         )
    (cl-mpm::update-sim *sim*)
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                    (format t "CFL dt estimate: ~f~%" dt-e)
                    (format t "CFL step count estimate: ~D~%" substeps-e)
                    (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
    (vgplot:figure)
    (time (loop for steps from 0 to 40
                while *run-sim*
                do
                   (progn
                     (when (> steps 5)
                       (setf (cl-mpm::sim-enable-damage *sim*) t))
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     (let ((v-av 0d0))
                       (time
                        (dotimes (i substeps)
                          (cl-mpm::update-sim *sim*)
                          (incf v-av
                                (/
                                 (lparallel:pmap-reduce
                                  (lambda (mp)
                                    (magicl:tref (cl-mpm/particle::mp-velocity mp) 0 0))
                                  #'+
                                  (cl-mpm:sim-mps *sim*))
                                 (* (length (cl-mpm:sim-mps *sim*)) substeps)))
                          (incf *t* (cl-mpm::sim-dt *sim*))))
                       (format t "Average V: ~F~%" v-av)
                       (push *t* *data-t*)
                       (push v-av *data-v*))
                     (incf *sim-step*)
                     ;; (vgplot:plot *data-t* *data-v*)
                     (plot *sim*)
                     (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                        :terminal "png size 1920,1080"
                                        )

                     (swank.live:update-swank)
                     ))))
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*))

(defun test-convergance ()
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
  (defparameter *conv-data-t* (list))
  (defparameter *conv-data-v* (list))
  (defparameter *conv-data-refine* (list))
  (vgplot:figure)
  (loop for refine in (list 1 2 3 4 5)
        do
           (progn
             (setup :angle 10d0 :refine refine :mps 2)
             ;; (setup :angle 10d0 :mps (+ refine 1))
             (defparameter *data-t* (list))
             (defparameter *data-v* (list))
             (defparameter *sim-step* 0)
             (let* ((target-time 0.1d0)
                    (dt-scale 0.5d0)
                    (dt (* dt-scale (cl-mpm/setup::estimate-elastic-dt *sim*)))
                    (substeps (floor target-time dt)))
               (format t "CFL dt estimate: ~f~%" (cl-mpm:sim-dt *sim*))
               (format t "Substeps ~D~%" substeps)
               (time (loop for steps from 0 to 40
                           while *run-sim*
                           do
                              (progn
                                (format t "Step ~d ~%" steps)
                                (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~2,'0d_~5,'0d.vtk" refine *sim-step*)) *sim*)
                                (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_~2,'0d_nodes_~5,'0d.vtk" refine *sim-step*)) *sim*)
                                (incf *sim-step*)
                                (let ((v-av 0d0)
                                      (disp-av 0d0))
                                  (time
                                   (dotimes (i substeps)
                                     (cl-mpm::update-sim *sim*)
                                     (incf v-av
                                           (/
                                            (lparallel:pmap-reduce
                                             (lambda (mp)
                                               (magicl:tref (cl-mpm/particle::mp-velocity mp) 0 0))
                                             #'+
                                             (cl-mpm:sim-mps *sim*))
                                            (* (length (cl-mpm:sim-mps *sim*)) substeps)))
                                     (incf disp-av
                                           (/
                                            (lparallel:pmap-reduce
                                             (lambda (mp)
                                               (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0))
                                             #'+
                                             (cl-mpm:sim-mps *sim*))
                                            (* (length (cl-mpm:sim-mps *sim*)) substeps)))
                                     (incf *t* (cl-mpm::sim-dt *sim*))))
                                  (push *t* *data-t*)
                                  (push v-av *data-v*)
                                  ;; (push disp-av *data-v*)
                                  )
                                (apply #'vgplot:plot (reduce #'append (mapcar #'list
                                                                              (append *conv-data-t* (list *data-t*))
                                                                              (append *conv-data-v* (list *data-v*))
                                                                              (mapcar (lambda (x) (format nil "~A" x)) (append *conv-data-refine* (list refine)))
                                                                              )))
                                (swank.live:update-swank)
                                )))))

             (push *data-t* *conv-data-t*)
             (push *data-v* *conv-data-v*)
             (push refine *conv-data-refine*)
             (apply #'vgplot:plot (reduce #'append (mapcar #'list
                                                           *conv-data-t*
                                                           *conv-data-v*
                                                           (mapcar (lambda (x) (format nil "~A" x)) *conv-data-refine*)
                                                           )))
             ))

(defun test-varying-friction ()
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
  (defparameter *conv-data-t* (list))
  (defparameter *conv-data-v* (list))
  (defparameter *conv-data-refine* (list))
  (vgplot:figure)
  (loop for refine from 0d0 upto 1d0 by 0.1d0
        do
           (progn
             (setup :angle 10d0 :refine 1)
             (setf (cl-mpm/penalty::bc-penalty-friction *floor-bc*) refine)
             ;; (setup :angle 10d0 :mps (+ refine 1))
             (defparameter *data-t* (list))
             (defparameter *data-v* (list))
             (defparameter *sim-step* 0)
             (let* ((target-time 0.1d0)
                    (dt (cl-mpm/setup::estimate-elastic-dt *sim*))
                    (substeps (floor target-time dt)))
               (format t "CFL dt estimate: ~f~%" (cl-mpm:sim-dt *sim*))
               (format t "Substeps ~D~%" substeps)
               (time (loop for steps from 0 to 10
                           while *run-sim*
                           do
                              (progn
                                (format t "Step ~d ~%" steps)
                                (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~2,'0d_~5,'0d.vtk" refine *sim-step*)) *sim*)
                                (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_~2,'0d_nodes_~5,'0d.vtk" refine *sim-step*)) *sim*)
                                (incf *sim-step*)
                                (let ((v-av 0d0)
                                      (disp-av 0d0))
                                  (time
                                   (dotimes (i substeps)
                                     (cl-mpm::update-sim *sim*)
                                     (incf v-av
                                           (/
                                            (lparallel:pmap-reduce
                                             (lambda (mp)
                                               (magicl:tref (cl-mpm/particle::mp-velocity mp) 0 0))
                                             #'+
                                             (cl-mpm:sim-mps *sim*))
                                            (* (length (cl-mpm:sim-mps *sim*)) substeps)))
                                     (incf *t* (cl-mpm::sim-dt *sim*))))
                                  (setf disp-av
                                        (/
                                         (lparallel:pmap-reduce
                                          (lambda (mp)
                                            (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0))
                                          #'+
                                          (cl-mpm:sim-mps *sim*))
                                         (length (cl-mpm:sim-mps *sim*))))

                                  (push *t* *data-t*)
                                  (push v-av *data-v*))
                                (apply #'vgplot:plot (reduce #'append (mapcar #'list
                                                                              (append *conv-data-t* (list *data-t*))
                                                                              (append *conv-data-v* (list *data-v*))
                                                                              (mapcar (lambda (x) (format nil "~A" x)) (append *conv-data-refine* (list refine)))
                                                                              )))
                                (swank.live:update-swank)
                                )))))

             (push *data-t* *conv-data-t*)
             (push *data-v* *conv-data-v*)
             (push refine *conv-data-refine*)
             (apply #'vgplot:plot (reduce #'append (mapcar #'list
                                                           *conv-data-t*
                                                           *conv-data-v*
                                                           (mapcar (lambda (x) (format nil "~A" x)) *conv-data-refine*)
                                                           )))
             ))
(defun test-varying-angle ()
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
  (defparameter *conv-data-t* (list))
  (defparameter *conv-data-v* (list))
  (defparameter *conv-data-refine* (list))
  (vgplot:figure)
  (loop for refine in (list 5d0 10d0 20d0 30d0 40d0)
        do
           (progn
             (setup :angle refine :refine 1)
             (setf (cl-mpm/penalty::bc-penalty-friction *floor-bc*) 0.2d0)
             ;; (setup :angle 10d0 :mps (+ refine 1))
             (defparameter *data-t* (list))
             (defparameter *data-v* (list))
             (defparameter *sim-step* 0)
             (let* ((target-time 0.1d0)
                    (dt (cl-mpm/setup::estimate-elastic-dt *sim*))
                    (substeps (floor target-time dt)))
               (format t "CFL dt estimate: ~f~%" (cl-mpm:sim-dt *sim*))
               (format t "Substeps ~D~%" substeps)
               (time (loop for steps from 0 to 10
                           while *run-sim*
                           do
                              (progn
                                (format t "Step ~d ~%" steps)
                                (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~2,'0d_~5,'0d.vtk" refine *sim-step*)) *sim*)
                                (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_~2,'0d_nodes_~5,'0d.vtk" refine *sim-step*)) *sim*)
                                (incf *sim-step*)
                                (let ((v-av 0d0)
                                      (disp-av 0d0))
                                  (time
                                   (dotimes (i substeps)
                                     (cl-mpm::update-sim *sim*)
                                     (incf v-av
                                           (/
                                            (lparallel:pmap-reduce
                                             (lambda (mp)
                                               (magicl:tref (cl-mpm/particle::mp-velocity mp) 0 0))
                                             #'+
                                             (cl-mpm:sim-mps *sim*))
                                            (* (length (cl-mpm:sim-mps *sim*)) substeps)))
                                     (incf *t* (cl-mpm::sim-dt *sim*))))
                                  (setf disp-av
                                        (/
                                         (lparallel:pmap-reduce
                                          (lambda (mp)
                                            (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0))
                                          #'+
                                          (cl-mpm:sim-mps *sim*))
                                         (length (cl-mpm:sim-mps *sim*))))

                                  (push *t* *data-t*)
                                  (push v-av *data-v*))
                                (apply #'vgplot:plot (reduce #'append (mapcar #'list
                                                                              (append *conv-data-t* (list *data-t*))
                                                                              (append *conv-data-v* (list *data-v*))
                                                                              (mapcar (lambda (x) (format nil "~A" x)) (append *conv-data-refine* (list refine)))
                                                                              )))
                                (swank.live:update-swank)
                                )))))

             (push *data-t* *conv-data-t*)
             (push *data-v* *conv-data-v*)
             (push refine *conv-data-refine*)
             (apply #'vgplot:plot (reduce #'append (mapcar #'list
                                                           *conv-data-t*
                                                           *conv-data-v*
                                                           (mapcar (lambda (x) (format nil "~A" x)) *conv-data-refine*)
                                                           )))
             ))
