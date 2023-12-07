(defpackage :cl-mpm/examples/chalk
  (:use :cl))
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)
(in-package :cl-mpm/examples/chalk)
(declaim (optimize (debug 3) (safety 3) (speed 0)))

(defun plot (sim)
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :deformed
   ;; :colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :xx))
   :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-damage-ybar mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-yield-func mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-strain-plastic-vm mp))
   ;; :colour-func #'cl-mpm/particle::mp-strain-plastic-vm
   ))

(defun setup-test-column (size block-size offset &optional (e-scale 1) (mp-scale 1))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               #'cl-mpm/shape-function:make-shape-function-bspline
               ;; 'cl-mpm::mpm-sim-usf
               'cl-mpm/damage::mpm-sim-damage-nd-2
               ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 1.7d3)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
         )
    (declare (double-float h density))
    (progn
      (let ()
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-mps-from-list
               (cl-mpm/setup::make-block-mps-list
                ;; (mapcar (lambda (x) 0) size)
                offset
                ;; '(0 0)
                block-size
                (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                density
                ;; 'cl-mpm/particle::particle-elastic
                ;; 'cl-mpm/particle::particle-chalk-brittle
                'cl-mpm/particle::particle-chalk-delayed
                ;; 'cl-mpm/particle::particle-limestone
                ;; 'cl-mpm/particle::particle-mc
                ;; 'cl-mpm/particle::particle-vm
                :E 1d9
                :nu 0.35d0
                ;; :psi (* 40d0 (/ pi 180))
                ;; :phi (* 00d0 (/ pi 180))
                ;; :c 0.4d6

                :rho 1d6
                :enable-plasticity nil

                :coheasion 1d5
                :friction-angle 40d0

                :fracture-energy 1000d0
                :initiation-stress 100d3
                :delay-time 1d4
                ;; :compression-ratio 8d0


                :critical-damage 1d0;0.999d0
                :local-length 10d0
                ;; :local-length-damaged 1d0
                :local-length-damaged 1d0

                :gravity -9.8d0
                :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
                ))))
      (setf (cl-mpm:sim-allow-mp-split sim) t)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-enable-fbar sim) nil)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) t)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (setf (cl-mpm::sim-mass-filter sim) 1d-15)
      (let ((ms 1d5))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        ;; (setf (cl-mpm:sim-damping-factor sim) (* 1d-2 density ms))
        ;; (setf (cl-mpm:sim-damping-factor sim) 10.0d0)
        (setf (cl-mpm:sim-damping-factor sim) (* 1d-1 ms))
        )

      (dotimes (i 0)
        (dolist (dir (list :x :y))
          (cl-mpm::split-mps-criteria
           sim
           (lambda (mp h)
             (when
                 (and
                  (> (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                     120)
                  (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                     200
                     )
                  (> (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                     20
                     )
                  )
               dir
               )))))

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

      (defparameter *floor-bc*
        (cl-mpm/penalty::make-bc-penalty-point-normal
         sim
         (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
         (cl-mpm/utils:vector-from-list (list 00d0 (+ h-y) 0d0))
         (* density 1d5)
         0.9d0
         ;; 1d1
         ))

      ;; (setf (cl-mpm::sim-bcs-force-list sim)
      ;;       (list
      ;;        (cl-mpm/bc:make-bcs-from-list
      ;;         (list
      ;;          *floor-bc*
      ;;          ))))

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

(defmethod cl-mpm::post-stress-step (mesh (mp cl-mpm/particle::particle-vm) dt)
  (with-accessors ((ps cl-mpm/particle::mp-strain-plastic-vm)
                   (rho cl-mpm/particle::mp-rho))
      mp
    ;; (let ((rho_0 1d6)
    ;;       (rho_1 1d1)
    ;;       (soft 5d1))
    ;;   (setf rho (max rho_1
    ;;                  (* rho_0 (exp (- (* soft ps)))))))
    ))

(defparameter *sim* nil)
(defparameter *run-sim* t)
(defparameter *t* 0)
(defparameter *sim-step* 0)
(defparameter *refine* (/ 1d0 2d0))
(let ((refine (uiop:getenv "REFINE")))
  (when refine
    (setf *refine* (parse-integer (uiop:getenv "REFINE")))
    ))

(defun setup (&key (undercut 0d0))
  ;; (let ((mps-per-dim 4))
  ;;   (defparameter *sim* (setup-test-column '(16 16) '(8 8)  '(0 0) *refine* mps-per-dim)))
  ;; (defparameter *sim* (setup-test-column '(1 1 1) '(1 1 1) 1 1))

  (let* ((mesh-size 10)
         (mps-per-cell 2)
         (shelf-height 100)
         (soil-boundary 50)
         (shelf-aspect 2)
         (runout-aspect 2)
         (shelf-length (* shelf-height shelf-aspect))
         (domain-length (+ shelf-length (* runout-aspect shelf-height)))
         (shelf-height (+ shelf-height soil-boundary))
         (offset (list 0 (* 0 mesh-size)))
         )
    (defparameter *sim*
      (setup-test-column (list domain-length
                               (+ shelf-height 100))
                         (list domain-length shelf-height)
                         offset
                         (/ 1d0 mesh-size) mps-per-cell))

    (let* ((undercut-angle undercut)
           (normal (magicl:from-list (list
                                      (cos (- (* pi (/ undercut-angle 180d0))))
                                      (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1))))
      (cl-mpm/setup:remove-sdf *sim*
                               (lambda (p)
                                 (if (> (magicl:tref p 1 0) soil-boundary)
                                     (cl-mpm/setup::plane-point-sdf
                                      (magicl:from-list (list (magicl:tref p 0 0)
                                                              (magicl:tref p 1 0)) '(2 1))
                                      normal
                                      (magicl:from-list (list shelf-length soil-boundary)
                                                        '(2 1) :type 'double-float))
                                     1d0)
                                 )))
    ;; (cl-mpm/setup::damage-sdf *sim* (lambda (p) (cl-mpm/setup::line-sdf
    ;;                                              (magicl:from-list (list (magicl:tref p 0 0)
    ;;                                                                      (magicl:tref p 1 0)) '(2 1))
    ;;                                              (list (- shelf-length shelf-height) shelf-height)
    ;;                                              (list shelf-length soil-boundary)
    ;;                                              10d0
    ;;                                              )) 0.9d0)
    )
  ;; (loop for mp across (cl-mpm:sim-mps *sim*)
  ;;       do (setf (cl-mpm/particle::mp-damage mp) (random 0.1d0)))
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)

  (let* ((target-time 1d3)
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
                     (when (= steps 1)
                       (setf (cl-mpm::sim-enable-damage *sim*) t)
                       ;; (setf target-time 1d2)
                       ;; (setf (cl-mpm::sim-mass-scale *sim*) 1d2
                       ;;       target-time 1d1)
                       ;; (setf target-time 1d1
                       ;;       dt-scale 0.5d0)
                       ;; (let ((ms (cl-mpm::sim-mass-scale *sim*)))
                       ;;   (setf (cl-mpm:sim-damping-factor *sim*) (* 0d-8 ms))
                       ;;   )
                       (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                         (format t "CFL dt estimate: ~f~%" dt-e)
                         (format t "CFL step count estimate: ~D~%" substeps-e)
                         (setf substeps substeps-e))
                       )

                     (let ((new-target-time 1d2))
                       (when (not (= target-time new-target-time))
                         (loop for mp across (cl-mpm:sim-mps *sim*)
                               do
                                  (when (>= (cl-mpm/particle:mp-damage mp) 0.99d0)
                                    (setf target-time new-target-time
                                          dt-scale 1.0d0
                                          (cl-mpm::sim-mass-scale *sim*) 1d4
                                          )

                                    (let ((ms (cl-mpm::sim-mass-scale *sim*)))
                                      (setf (cl-mpm:sim-damping-factor *sim*) (* 1d-4 ms)))
                                    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                                      (format t "CFL dt estimate: ~f~%" dt-e)
                                      (format t "CFL step count estimate: ~D~%" substeps-e)
                                      (setf substeps substeps-e))
                                    (loop-finish)

                                    ))))
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     ;; (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     (time
                      (let ((current-damage (cl-mpm::sim-enable-damage *sim*)))
                        ;; (setf (cl-mpm::sim-enable-damage *sim*) nil)
                        (dotimes (i substeps)
                          ;; (if (= i (- substeps 1))
                          ;;     (setf (cl-mpm::sim-enable-damage *sim*) current-damage)
                          ;;     (setf cl-mpm/damage::*delocal-counter-max* 0)
                          ;;     )
                          (cl-mpm::update-sim *sim*)
                          (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)));)
                        ;; (setf (cl-mpm::sim-enable-damage *sim*) current-damage)
                        ;; (when current-damage
                        ;;   (setf (cl-mpm::sim-enable-damage *sim*) t)
                        ;;   (cl-mpm::update-sim *sim*)
                          )

                        ))

                     ;; (loop for mp across (cl-mpm:sim-mps *sim*)
                     ;;       do
                     ;;          (when (>= (cl-mpm/particle:mp-damage mp) 1d0)
                     ;;            (let ((ms 1d2))
                     ;;              (setf (cl-mpm::sim-mass-scale *sim*) ms
                     ;;                    ;; target-time 1d0
                     ;;                    (cl-mpm:sim-damping-factor *sim*) (* 1d-8 ms)
                     ;;                    ))))

                     (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       (format t "CFL dt estimate: ~f~%" dt-e)
                       (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf substeps substeps-e))
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
  (sb-profile:profile "CL-MPM")
  (sb-profile:profile "CL-MPM/PARTICLE")
  (sb-profile:profile "CL-MPM/MESH")
  (sb-profile:profile "CL-MPM/SHAPE-FUNCTION")
  (sb-profile:reset)
  (time
   (dotimes (i 100)
         (cl-mpm::update-sim *sim*)))
  (sb-profile:report)
  )
(defun test-undercut ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output_chalk/mesh.vtk")
                               *sim*)
  (loop for c in (list 0d0 10d0 20d0 30d0 40d0 50d0)
        while *run-sim*
        do
           (progn
             (setup :undercut (- c))
             (run)
             (cl-mpm/output:save-vtk (merge-pathnames (format nil "output_chalk/chalk_~5,'0d.vtk" c)) *sim*)
             )))

(defun profile ()
  (sb-profile:unprofile)
  (sb-profile:reset)
  (sb-profile:profile "CL-MPM")
  (sb-profile:profile "CL-MPM/PARTICLE")
  (sb-profile:profile "CL-MPM/MESH")
  (loop repeat 100
        do (progn
             (cl-mpm::update-sim *sim*)
             ;; (cl-mpm/damage::calculate-damage (cl-mpm:sim-mesh *sim*)
             ;;                                  (cl-mpm:sim-mps *sim*)
             ;;                                  (cl-mpm:sim-dt *sim*)
             ;;                                  25d0)
             ;; (cl-mpm/eigenerosion:update-fracture *sim*)
             ))
  (sb-profile:report))

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


(defmacro time-form (it form)
  `(progn
     (declaim (optimize speed))
     (let* ((iterations ,it)
            (start (get-internal-real-time)))
       (dotimes (i ,it)
         ,form)
       (let* ((end (get-internal-real-time))
              (units internal-time-units-per-second)
              (dt (/ (- end start) (* iterations units)))
              )
         (format t "Total time: ~f ~%" (/ (- end start) units)) (format t "Time per iteration: ~f~%" (/ (- end start) (* iterations units)))
         (format t "Throughput: ~f~%" (/ 1 dt))
         dt))))

(defun test ()
  (let ((data-cores '())
        (data-dt '()))
    (loop for k from 1 to 8
          do
             (progn
               (setf lparallel:*kernel* (lparallel:make-kernel k :name "custom-kernel"))
               (format t "Kernel: ~D~%" k)
               (let* ((iterations 100000)
                      (start (get-internal-real-time)))
                 (let ((a (cl-mpm/utils:vector-zeros)))
                   (time
                    (lparallel:pdotimes (i iterations)
                      (cl-mpm::update-strain-kirchoff (cl-mpm:sim-mesh *sim*) (aref (cl-mpm:sim-mps *sim*) 0) 0d0 nil)
                      )))
                 (let* ((end (get-internal-real-time))
                        (units internal-time-units-per-second)
                        (dt (/ (- end start) (* iterations units)))
                        )
                   (format t "Total time: ~f ~%" (/ (- end start) units)) (format t "Time per iteration: ~f~%" (/ (- end start) (* iterations units)))
                   (format t "Throughput: ~f~%" (/ 1 dt))
                   (push (/ 1d0 dt) data-dt)
                   (push k data-cores)
                   ))

               ;; (let ((iters 100000))
               ;;   (let ((a (cl-mpm/utils:vector-zeros)))
               ;;     (time
               ;;      (lparallel:pdotimes (i iters)
               ;;        (magicl:@ (magicl:eye 3) (cl-mpm/utils:vector-zeros))))))
               (lparallel:end-kernel)
               ))
    (vgplot:close-all-plots)
    (vgplot:plot data-cores data-dt))
  )


(defun test-mc (exx eyy ezz eyz ezx exy)
  (let* ((eps (cl-mpm/utils:voigt-from-list (list exx eyy ezz eyz ezx exy)))
         (E 1d0)
         (nu 0.0d0)
         (angle 0.1d0)
         (de (cl-mpm/constitutive::linear-elastic-matrix E nu))
         (sig (magicl:@ de eps)))
    (cl-mpm/constitutive::mc-plastic sig de eps E nu angle angle 0d0)))




(defun plot-stress-damage ()
  (vgplot:close-all-plots)
  (let* ((mp (aref (cl-mpm:sim-mps *sim*))))
    (with-accessors ((damage cl-mpm/particle:mp-damage)
                     (E cl-mpm/particle::mp-e)
                     (Gf cl-mpm/particle::mp-Gf)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (length cl-mpm/particle::mp-local-length)
                     ) mp
      (let* ((stress (loop for x from 1d5 to 1d6 by 1d4
                           collect x))
             (damage (mapcar (lambda (stress)
                               ;(cl-mpm/damage::brittle-concrete-linear-d stress E Gf length init-stress)
                               (cl-mpm/damage::test-damage stress 2d5 1d5)
                               ;; (cl-mpm/damage::damage-response-exponential stress E Gf length init-stress)
                               )
                             stress)))
        (vgplot:plot stress damage)))))


(defun plot-interaction ()
  (vgplot:close-all-plots)
  (let* ((length 1d0)
         (x (loop for x from -5 to 5 by 0.1d0 collect x)))
    (vgplot:plot x (mapcar (lambda (x) (cl-mpm/damage::weight-func x length)) x))
    ))

(let* ((s (cl-mpm/utils::stress-from-list (list 500d3 0d0 0d0 0d0 0d0 0d0)))
       (j2 (cl-mpm/constitutive::voigt-j2
            (cl-mpm/utils::deviatoric-voigt s)))
       (p (cl-mpm/constitutive::voight-trace s))
       (ft 300d3)
       (fc 500d3)
       (angle (atan (* 3 (/ (- fc ft) (+ fc ft)))))
       (k (* (/ 2d0 (sqrt 3)) (/ (* ft fc) (+ ft fc))))
       (s_1 (-
             (* (/ 3d0 (+ 3 (tan angle)))
                (+ (sqrt (* 3 j2)) (* 1/3 (tan angle) p)))
             k
             )))
  (format t "~A~%" s_1))

(defun test (s1 s2 s3)
  (let ((damage-inc-mat (cl-mpm/utils:matrix-zeros))
        (damage-tensor (cl-mpm/utils:voight-to-matrix (cl-mpm/utils:voigt-from-list (list 0.1d0 0.5d0 0d0 0d0 0d0 0d0))))
        (ybar-tensor (cl-mpm/utils:matrix-zeros))
        (E 1d9)
        (length 10d0)
        (Gf 1000d0)
        (init-stress 100d3)
        (cauchy-undamaged (cl-mpm/utils:voight-to-matrix (cl-mpm/utils:voigt-from-list (list s1 s2 s3 0d0 0d0 0d0)))))
    (multiple-value-bind (ls v) (cl-mpm/utils:eig cauchy-undamaged)
      (loop for i from 0 to 2
            do
               (let* ((sii (nth i ls))
                      (vii (magicl::column v i))
                      (vsi (magicl:@ vii (magicl:transpose vii)))
                      (dii (magicl::trace (magicl:@ damage-tensor vsi)))
                      ;; (dii 0d0)
                      (new-damage (cl-mpm/damage::damage-response-exponential sii E Gf length init-stress))
                      (damage-increment (- (max dii new-damage) dii))
                      )
                 ;; (when (> damage-increment 0d0)
                 ;;   (break))
                 (format t "Sii ~F : new-damage ~F ~%" sii new-damage)
                 (format t "Damage  dim ~D tensor ~A~%" i (magicl:@ damage-tensor vsi))
                 (format t "Damage  dim ~D vsi ~A~%" i vsi)
                 ;; (magicl:.+ ybar-tensor
                 ;;            (magicl:scale! vsi sii)
                 ;;            ybar-tensor)
                 (magicl:.+ damage-inc-mat
                            (magicl:scale! vsi damage-increment)
                            ;; (magicl:scale! vsi (* (/ dt tau) (- 1d0 (exp (- (* 1d0 (abs (- new-damage dii))))))))
                            damage-inc-mat))))
    (format t "Damage inc ~A~%" damage-inc-mat)
    (magicl:.+ damage-tensor
               damage-inc-mat
               damage-tensor)

    (format t "Damage tensor ~A~%" damage-tensor)
    )
  )
