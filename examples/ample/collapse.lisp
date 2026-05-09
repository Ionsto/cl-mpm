(defpackage :cl-mpm/examples/collapse
  (:use :cl))

(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)
(setf cl-mpm/settings:*optimise-setting* cl-mpm/settings::*optimise-debug*)
;; (setf cl-mpm/settings:*optimise-setting* cl-mpm/settings::*optimise-speed*)
;; (sb-ext:restrict-compiler-policy 'speed  3 3)
;; (sb-ext:restrict-compiler-policy 'debug  0 0)
;; (sb-ext:restrict-compiler-policy 'safety 0 0)
;; (setf *block-compile-default* t)
;(sb-int:set-floating-point-modes :traps '(:overflow :invalid :inexact :divide-by-zero :underflow))
;; (sb-int:set-floating-point-modes :traps '(:overflow :divide-by-zero :underflow))

(in-package :cl-mpm/examples/collapse)
;; (declaim (optimize (debug 0) (safety 0) (speed 3)))
(declaim (optimize (debug 3) (safety 3) (speed 0)))

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-elastic) dt fbar)
  (cl-mpm::update-stress-kirchoff mesh mp dt fbar))

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-elastic) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  (cl-mpm::update-domain-stretch mesh mp dt)
  ;; (cl-mpm::update-domain-polar mesh mp dt)
  )

(defun plot-load-disp ()
  (vgplot:semilogy *data-steps* *data-energy*))

(declaim (notinline plot))
(defun plot (sim)
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
         (h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
         (ms-x (first ms))
         (ms-y (second ms))))
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :deformed
   :trial nil
   ;; :colour-func #'cl-mpm/particle::mp-index
   :colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :xy))
  ))

(defparameter *eta* 1d0)

(defun stop ()
  (setf *run-sim* nil))

(declaim (notinline setup-test-column))
(defun setup-test-column (size block-size sim-type &optional (e-scale 1) (mp-scale 1) (multigrid-refinement 0))
  (let* ((multigrid-enabled (> multigrid-refinement 0))
         (E 1d6)
         (density 1d3)
         (sim (cl-mpm/setup::make-simple-sim
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               ;; :sim-type 'cl-mpm/aggregate::mpm-sim-agg-usf
               :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static
               ;; :sim-type (if multigrid-enabled
               ;;               'cl-mpm/dynamic-relaxation::mpm-sim-dr-multigrid
               ;;               'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul)
               :args-list
               (append
                (list
                 :split-factor (* (sqrt 2) (/ 1d0 mp-scale))
                 :enable-fbar nil
                 :enable-aggregate nil
                 :ghost-factor nil;(* E 1d-0)
                 ;; :ghost-factor (* density 1d-6)
                 :mass-update-count 1
                 :damping-update-count 1
                 :max-split-depth 6
                 :enable-split nil
                 :gravity -10d0
                 )
                (when multigrid-enabled
                  (list :refinement multigrid-refinement))
                )))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
         (offset (* h 0)))
    (declare (double-float h density))
    (progn
      (let* ((angle 40d0)
             )
        (cl-mpm::add-mps
         sim
         (cl-mpm/setup::make-block-mps
          (mapcar #'+
                  (mapcar (lambda (x) 0) size)
                  (list 0d0 offset 0))
          block-size
          (mapcar (lambda (e) (* (/ e h) mp-scale)) block-size)
          density
          'cl-mpm/particle::particle-elastic
          :E E
          :nu 0.3d0
          ;; 'cl-mpm/particle::particle-vm
          ;; :E E
          ;; :nu 0.3d0
          ;; :rho 20d3
          ;; 'cl-mpm/particle::particle-mc
          ;; :E 1d6
          ;; :nu 0.3d0
          ;; :phi (* 30d0 (/ 3.14d0 180d0))
          ;; :psi (* 00d0 (/ 3.14d0 180d0))
          ;; :c 10d3
          ;; :enable-plasticity t
          ;; 'cl-mpm/particle::particle-elastic-damage-delayed
          ;; :E 1d9
          ;; :nu 0.3d0
          ;; :initiation-stress 1d4
          ;; :local-length (* 2 h)
          ;; :ductility 50d0
          ;; :delay-time 100d0
          ;; :delay-exponent 2d0

          ;; 'cl-mpm/particle::particle-chalk-erodable
          ;; :E E
          ;; :nu 0.24d0
          ;; :enable-damage t
          ;; :enable-plasticity t
          ;; :friction-angle angle

          ;; :kt-res-ratio 1d0
          ;; :kc-res-ratio 0d0
          ;; :g-res-ratio 0.51d0

          ;; :initiation-stress 1d3
          ;; :delay-time 1d0
          ;; :delay-exponent 1d0
          ;; :ductility ductility
          ;; :local-length h
          ;; :psi (* 0d0 (/ pi 180))
          ;; :phi (* angle (/ pi 180))
          ;; :c (* init-c oversize)
          ;; :softening 0d0
          ;; :gravity -10d0
          ;; :gravity-axis (cl-mpm/utils:vector-from-list '(-1d0 1d0 0d0))
          )))
      ;; (format t "Charictoristic time ~E~%" (/ ))
      ;; (setf (cl-mpm:sim-allow-mp-split sim) t)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      ;; (setf (cl-mpm::sim-enable-fbar sim) t)
      (defparameter *density* density)
      (cl-mpm/setup::set-mass-filter sim density :proportion 1d-15)

      ;; (let ((dt-scale 1d0))
      ;;   (setf
      ;;    (cl-mpm:sim-dt sim)
      ;;    (* dt-scale h
      ;;       (sqrt (cl-mpm::sim-mass-scale sim))
      ;;       (sqrt (/ density (cl-mpm/particle::mp-p-modulus (aref (cl-mpm:sim-mps sim) 0)))))))
      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))

      (cl-mpm/setup::setup-bcs
       sim
       :top '(nil nil nil)
       :bottom '(nil 0 nil)
       ;; :left '(0 nil nil)
       ;; :front '(nil nil 0)
       ))
    sim)
  )


(defparameter *sim* nil)
(defparameter *run-sim* t)
(defparameter *t* 0)
(defparameter *sim-step* 0)
(defparameter *refine* (/ 1d0 1d0))
(let ((refine (uiop:getenv "REFINE")))
  (when refine
    (setf *refine* (parse-integer (uiop:getenv "REFINE")))
    ))

(declaim (notinline setup))
(defun setup (&key (refine 1) (mps 2)
                (sim-type 'cl-mpm:mpm-sim-usf)
                (multigrid-refine 0))
  (defparameter *sim* nil)
  (let ((mps-per-dim mps))
    (setf *sim* (setup-test-column '(16 16) '(8 8) sim-type refine mps-per-dim multigrid-refine))
    ;; (setf *sim* (setup-test-column '(8 8) '(8 8) sim-type refine mps-per-dim multigrid-refine))
    ;; (setf *sim* (setup-test-column '(32 32 16) '(8 8 8) sim-type refine mps-per-dim multigrid-refine))
    ;; (setf *sim* (setup-test-column '(64 32) '(16 16) sim-type refine mps-per-dim multigrid-refine))
    )
  ;; (cl-mpm/setup::initialise-stress-self-weight
  ;;  *sim*
  ;;  15d0
  ;;  )
  (cl-mpm/setup::setup-bcs
   *sim*
   :left (list 0 nil nil)
   :bottom (list nil 0 nil)
   :top (list nil nil nil)
   )
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defun setup-scaling (&key (refine 1) (mps 2)
                           (refine-0 2)
                (sim-type 'cl-mpm:mpm-sim-usf)
                (multigrid-refine 0))
  (defparameter *sim* nil)
  (let ((mps-per-dim mps))
    (setf *sim* (setup-test-column (list (* refine 8) 9) (list (* refine 8) 8) sim-type refine-0 mps-per-dim multigrid-refine))
    ;; (setf *sim* (setup-test-column '(8 8) '(8 8) sim-type refine mps-per-dim multigrid-refine))
    ;; (setf *sim* (setup-test-column '(32 32 16) '(8 8 8) sim-type refine mps-per-dim multigrid-refine))
    ;; (setf *sim* (setup-test-column '(64 32) '(16 16) sim-type refine mps-per-dim multigrid-refine))
    )
  ;; (cl-mpm/setup::initialise-stress-self-weight
  ;;  *sim*
  ;;  15d0
  ;;  )
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defun save-csv (output-file)
  (let* ((mp-list
           (loop for mp across (cl-mpm:sim-mps *sim*)
                 collect mp))
         (y (loop for mp in mp-list collect (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)))
         (y-ref (loop for pos in *original-configuration*
                      collect (float pos 1e0)))
         (syy (loop for mp in mp-list collect (float (magicl:tref (cl-mpm/particle::mp-stress mp) 1 0) 1e0)))
         (rho 80d0)
         (E 1d5)
         (g 10d0)
         (vp-0-list (loop for size in *original-size*
                          collect (float (* (cl-mpm/utils:varef size 0) (cl-mpm/utils:varef size 1)) 1e0)))
         (pressure -1e4)
         (max-y 50)
         (syy-ref (mapcar (lambda (x) pressure) y-ref))
         (df (lisp-stat:make-df '(:y :syy :syy-ref :vp)
                                (list
                                 (coerce y-ref '(vector single-float))
                                 (coerce syy 'vector)
                                 (coerce syy-ref 'vector)
                                 (coerce vp-0-list 'vector)))))
    (lisp-stat:write-csv df output-file :add-first-row t)))



(defun save-data (name)
  (with-open-file (stream (merge-pathnames name) :direction :output :if-exists :supersede)
    (format stream "time,oobf,energy,damage~%")
    (loop for time in (reverse *data-t*)
          for oobf in (reverse *data-oobf*)
          for energy in (reverse *data-energy*)
          for damage in (reverse *data-damage*)
          do (format stream "~f,~f,~f,~f~%" time oobf energy damage))))


(defmacro time-form (it form)
  `(progn
     (declaim (optimize speed))
     (let* ((iterations ,it)
            (start (get-internal-real-time))
            (gc-start sb-ext:*gc-real-time*)
            )
       (time
        (dotimes (i ,it)
              ,form))
       (let* ((end (get-internal-real-time))
              (gc-end sb-ext:*gc-real-time*)
              (units internal-time-units-per-second)
              (dt (/ (- end start) (* iterations units))))
         (format t "Total time: ~f ~%" (/ (- end start) units))
         (format t "Time per iteration: ~f~%" (/ (- end start) (* iterations units)))
         (format t "Total gc time: ~f ~%" (/ (- gc-end gc-start) units))
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
  (format t "MPS ~D~%" (length (cl-mpm:sim-mps *sim*)))
  ;; (sb-profile:report)
  )

(defun test-real-time ()
  (setup :mps 2 :refine 1)
  (let ((output-dir (format nil "./output-TFLIP/")))
    (vgplot:close-all-plots)
    ;; (change-class *sim* 'cl-mpm/aggregate::mpm-sim-agg-usf)
    (change-class *sim* 'cl-mpm/aggregate::mpm-sim-agg-usf)
    ;; (change-class *sim* 'cl-mpm/aggregate::mpm-sim-musl)
    ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic)
    (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t)
    ;(setf (cl-mpm::sim-velocity-algorithm *sim*) :TFLIP)
    ;; (setf (cl-mpm::sim-velocity-algorithm *sim*) :TBLEND)
    ;; (setf (cl-mpm::sim-ghost-factor *sim*) (* 1d3 1d-6))
    ;; (setf (cl-mpm::sim-velocity-algorithm *sim*) :PIC)
    (setf (cl-mpm::sim-velocity-algorithm *sim*) :BLEND)
    (cl-mpm/setup::set-mass-filter *sim* *density* :proportion 1d-9)
    (cl-mpm:update-sim *sim*)
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (let ((step 0))
      (cl-mpm/dynamic-relaxation::run-time
       *sim*
       :output-dir output-dir
       :plotter (lambda (sim)
                  (plot *sim*)
                  (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step)) :terminal "png size 1920,1080")
                  (incf step))
       :total-time 100d0
       :damping 1d0
       :dt 0.5d0
       :initial-quasi-static nil
       :dt-scale 0.9d0))))

(defun test-dt ()
  (dolist (r (list 1.5d0 1d0 0.5d0 0.25d0))
    (setup :refine 1d0 :mps 3)
    (let ((output-dir (format nil "./output-agg-wrong-~F/" r)))
      (format t "Running with dt scale ~F~%" r)
      (vgplot:close-all-plots)
      ;; (change-class *sim* 'cl-mpm/aggregate::mpm-sim-agg-usf)
      (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t)
      (setf (cl-mpm::sim-velocity-algorithm *sim*) :FLIP)
      (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
      (let ((step 0))
        (cl-mpm/dynamic-relaxation::run-time
         *sim*
         :output-dir output-dir
         :plotter (lambda (sim)
                    (plot *sim*)
                    (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step)) :terminal "png size 1920,1080")
                    (incf step))
         :total-time 100d0
         :damping 1d-1
         :dt 1d0
         :initial-quasi-static nil
         :dt-scale r)))))

(defun timing-test ()
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (mps cl-mpm::sim-mps))
      *sim*
    (cl-mpm::kill-workers)
    (cl-mpm::make-workers)
    ;; (cl-mpm::omp
    ;;  8
    ;;  (lambda (i)
    ;;    (pprint (lparallel:kernel-worker-index)))
    ;;  )
    (let* ((iters 1000)
           (dt 1d0))
      ;; (time
      ;;  (cl-mpm::update-stress-omp mesh mps dt nil))
      ;; (time
      ;;  (dotimes (i iters)
      ;;    (cl-mpm::iterate-over-nodes
      ;;     mesh
      ;;     (lambda (n)
      ;;       (cl-mpm/mesh::reset-node-force n)))
      ;;    ))
      ;; (let* ((mut (sb-thread:make-mutex))
      ;;        (nodes (cl-mpm/mesh::mesh-active-nodes mesh))
      ;;        (inner (length nodes))
      ;;        (iters 100))
      ;;   (time
      ;;    (dotimes (j iters)
      ;;      (lparallel:pdotimes (i inner)
      ;;        (cl-mpm/mesh::reset-node-force (aref nodes i)))))
      ;;   (time
      ;;    (dotimes (j iters)
      ;;      (cl-mpm::omp
      ;;       inner
      ;;       (lambda (i)
      ;;         (cl-mpm/mesh::reset-node-force (aref nodes i))
      ;;         )
      ;;       ))))
      ;; (time-form
      ;;  iters
      ;;  (cl-mpm::iterate-over-nodes
      ;;   mesh
      ;;   (lambda (n)
      ;;     (cl-mpm/mesh::reset-node-force n))))
      ;; (time-form
      ;;  iters
      ;;  (cl-mpm::iterate-over-nodes-omp
      ;;   mesh
      ;;   (lambda (n)
      ;;     (cl-mpm/mesh::reset-node-force n))))
      )
    )
  )


(defun test-scaling ()
  (let ((data-perf (list))
        (data-threads (list))
        )
    (vgplot:close-all-plots)
    (dolist (threads (list 8))
      (lparallel:end-kernel)
      (setf lparallel:*kernel* (lparallel:make-kernel threads))
      (sb-ext:gc :full t)
      (let ((t0 (get-internal-real-time)))
        (let ((agg nil)
              (r (* 1 threads)))
          (setup-scaling :mps 2
                         :refine r
                         :multigrid-refine 0
                         :refine-0 1)
          (cl-mpm::domain-sort-mps *sim*)
          ;; (cl-mpm::compact-mesh-active *sim*)
          (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul)
          ;; (change-class *sim* 'cl-mpm/implicit::mpm-sim-implicit)
          (cl-mpm/setup::set-mass-filter *sim* *density* :proportion 1d-15)
          (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
                (cl-mpm::sim-ghost-factor *sim*) nil;(* 1d6 1d-3)
                (cl-mpm::sim-enable-fbar *sim*) nil)
          ;; (if agg
          ;;     (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
          ;;           (cl-mpm::sim-ghost-factor *sim*) nil
          ;;           (cl-mpm::sim-enable-fbar *sim*) nil)
          ;;     (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) nil
          ;;           (cl-mpm::sim-ghost-factor *sim*) (* 1d6 1d-2)
          ;;           (cl-mpm::sim-enable-fbar *sim*) nil))
          (format t "threads ~d~%" threads)
          (setf (cl-mpm::sim-gravity *sim*) 0d0)
          (cl-mpm:update-sim *sim*)
          (format t "start test ~d~%" threads)
          (let (;; (testdata (make-array (* 1000000 r) :element-type 'double-float :initial-element 0d0))
                (mesh (cl-mpm:sim-mesh *sim*))
                (mps (cl-mpm:sim-mps *sim*))
                (mp (aref (cl-mpm:sim-mps *sim*) 0)))
            ;; (declare ((simple-array double-float *) testdata))
            (let ((sim *sim*)
                  (bcs-force (cl-mpm::sim-bcs-force *sim*))
                  (bcs-force-list (cl-mpm::sim-bcs-force-list *sim*))
                  (bcs (cl-mpm::sim-bcs *sim*))
                  (dt 1d0)
                  (dt-loadstep 1d0))
              ;; (sb-sprof:with-profiling (:max-samples 100000 :sample-interval 0.00001 :report :graph :mode :alloc))
              (let ((dt  (time-form
                          10
                          (progn
                            ;; (lparallel:pdotimes (i (length testdata))
                            ;;   (setf (aref testdata i) 0d0))
                            (cl-mpm:update-sim *sim*)
                            ;; (cl-mpm::reset-nodes-force sim)
                            ;; (cl-mpm::update-stress mesh mps dt-loadstep nil)
                            ;; (cl-mpm::p2g-force-fs sim)
                            ;; (cl-mpm::apply-bcs mesh bcs-force dt)
                            ;; (loop for bcs-f in bcs-force-list
                            ;;       do (cl-mpm::apply-bcs mesh bcs-f dt-loadstep))

                            ;; (when ghost-factor
                            ;;   (cl-mpm/ghost::apply-ghost-cached sim)
                            ;;   (cl-mpm::apply-bcs mesh bcs dt))

                            ;; (incf solve-count)
                            ;; (cl-mpm/dynamic-relaxation::update-node-fictious-mass sim)
                            ;; ;;update our nodes after force mapping
                            ;; (cl-mpm::update-node-forces sim)
                            ;; (cl-mpm::apply-bcs mesh bcs dt)
                            ;; (cl-mpm::update-nodes sim)
                            ;; (cl-mpm::update-filtered-cells sim)
                            ;; (cl-mpm::g2p mesh mps dt 0d0 :trial)
                            ;; (cl-mpm::update-stress mesh mps 1d0 nil)
                            ))))
                (push (float (/ 1d0 dt) 0d0) data-perf)
                (push threads data-threads)
                (vgplot:semilogx
                 ;; data-threads data-threads "ideal"
                 data-threads (mapcar (lambda (x) (/ x (first (last data-perf)))) data-perf) "scaling")
                )
              ;; )
          ;;   )
          (pprint data-perf)
          )
        (let ((t1 (get-internal-real-time)))
          (format t "threads ~d: time ~e~%" threads (/ (- t1 t0) internal-time-units-per-second))))))
    (pprint data-threads)
    (pprint data-perf)
    )))

(defun test ()
  (let ((agg t)
        (r 1))
    (setup :mps 4 :refine r :multigrid-refine 0)
    ;; (change-class *sim* 'cl-mpm/implicit::mpm-sim-implicit)
    (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul)
    (cl-mpm/setup::set-mass-filter *sim* *density* :proportion 1d-15)
    (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) nil
          (cl-mpm::sim-ghost-factor *sim*) nil
          (cl-mpm::sim-enable-fbar *sim*) nil)
    ;; (if agg
    ;;     (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
    ;;           (cl-mpm::sim-ghost-factor *sim*) nil
    ;;           (cl-mpm::sim-enable-fbar *sim*) nil)
    ;;     (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) nil
    ;;           (cl-mpm::sim-ghost-factor *sim*) (* 1d6 1d-2)
    ;;           (cl-mpm::sim-enable-fbar *sim*) nil))
    
    ;; (cl-mpm:update-sim *sim*)
    ;; (let ((mesh (cl-mpm:sim-mesh *sim*))
    ;;       (mp (aref (cl-mpm:sim-mps *sim*) 0)))
    ;;   ;; (time-form 1000
    ;;   ;;            ;; (cl-mpm::update-)
    ;;   ;;            ;; (cl-mpm::iterate-over-neighbours-cached
    ;;   ;;            ;;  mesh
    ;;   ;;            ;;  mp
    ;;   ;;            ;;  (lambda (node svp grads fsvp fgrads)
    ;;   ;;            ;;    )
    ;;   ;;            ;;  )
    ;;   ;;            )
    ;;   ;; (time
    ;;   ;;   (dotimes (i 1000000000)
    ;;   ;;     (cl-mpm::iterate-over-neighbours-cached
    ;;   ;;      mesh
    ;;   ;;      mp
    ;;   ;;      (lambda (node svp grads fsvp fgrads))
    ;;   ;;      )
    ;;   ;;     ;; (cl-mpm::calculate-strain-rate-disp (cl-mpm::sim-mesh *sim*)
    ;;   ;;     ;;                                     (aref (cl-mpm:sim-mps *sim*) 0)
    ;;   ;;     ;;                                     1d0)
    ;;   ;;     )
    ;;   ;;   ;; (cl-mpm:update-sim *sim*)
    ;;   ;;   )
    ;;   )
    (setf (cl-mpm::sim-gravity *sim*) -20d0)
    (time
     (cl-mpm/dynamic-relaxation::run-load-control
      *sim*
      :output-dir (format nil "./output-agg-~a-~d/" agg r)
      :plotter #'plot
      :load-steps 20
      :damping 1d0;(sqrt 2)
      :substeps 50
      :criteria 1d-9
      :save-vtk-dr nil
      :save-vtk-loadstep nil
      :dt-scale 1d0))
    ))

(defun test-load-control-agg ()
  (dolist (agg (list nil))
    (dolist (r (list 2 3 4))
      (setup :mps 2 :refine r :multigrid-refine 0)
      (cl-mpm/setup::set-mass-filter *sim* *density* :proportion 1d-15)
      (when agg
        (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
              (cl-mpm::sim-ghost-factor *sim*) nil))
      (setf (cl-mpm::sim-gravity *sim*) -1d0)
      (cl-mpm/dynamic-relaxation::run-load-control
       *sim*
       :output-dir (format nil "./output-agg-~a-~d/" agg r)
       :plotter #'plot
       :load-steps 3
       :damping 1d0;(sqrt 2)
       :substeps 10
       :criteria 1d-3
       :save-vtk-dr t
       :save-vtk-loadstep t
       :dt-scale 2d0))))


;; (let ((v (cl-mpm/aggregate::assemble-global-vec *sim* #'cl-mpm/mesh:node-force 0))
;;       (iters 1000000)
;;       )
;;   (time-form
;;    iters
;;    (let ((bcs (cl-mpm/aggregate::assemble-global-bcs *sim* 0)))
;;      (cl-mpm/fastmaths::fast-.* v bcs v))
;;    )
;;   (time-form
;;    iters
;;    (cl-mpm/aggregate::apply-global-bcs *sim* v 0))
;;   )


(defun timing-test ()
  (let* ((mesh (cl-mpm:sim-mesh *sim*))
         (nodes (cl-mpm/mesh::mesh-active-nodes mesh)))
    (pprint (length nodes))
    ;; (time-form
    ;;  100
    ;;  (cl-mpm::omp
    ;;   (length nodes)
    ;;   (lambda (n)
    ;;     ;; (dotimes (i 1000)
    ;;     ;;   (sqrt i))
    ;;     (cl-mpm/mesh::reset-node-force (aref nodes n))
    ;;     )))
    (let ((iters 1000))
      (time-form
       iters
       (lparallel:pdotimes (i (length nodes))
         (cl-mpm/mesh::reset-node-force (aref nodes i))))
      (time-form
       iters
       (cl-mpm::omp
        (length nodes)
        (lambda (n)
          (cl-mpm/mesh::reset-node-force (aref nodes n))))))
    ))

(defun test-load-control ()
  (dolist (agg (list nil))
    (dolist (r (list 1))
      (setup :mps 4 :refine r :multigrid-refine 0)
      (cl-mpm::domain-sort-mps *sim*)
      ;; (change-class *sim* 'cl-mpm/implicit::mpm-sim-implicit)
      (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static)
      ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-dr-usf)
      ;; (setf (cl-mpm::sim-velocity-algorithm *sim*) :tblend)
      (if agg
          (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
                (cl-mpm::sim-ghost-factor *sim*) nil)
          (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) nil
                (cl-mpm::sim-ghost-factor *sim*) nil))
      (cl-mpm/setup::set-mass-filter *sim* *density* :proportion 0d-9)
      (setf (cl-mpm::sim-gravity *sim*) -10d0)
      (let ((step (list))
            (res (list))
            (total-step 0))
        (vgplot:close-all-plots)
        (cl-mpm:update-sim *sim*)
        (time-form
         100
         (cl-mpm:update-sim *sim*))
        ;; (time
        ;;  (cl-mpm/dynamic-relaxation::run-load-control
        ;;   *sim*
        ;;   :output-dir (format nil "./output-dr-~a-~d/" (if agg "agg" "noagg") r)
        ;;   ;; :plotter (lambda (sim) (vgplot:semilogy (reverse step) (reverse res)))
        ;;   :plotter #'plot
        ;;   :load-steps 10
        ;;   :damping (sqrt 2d0)
        ;;   :substeps 50;(round (* 20 r))
        ;;   :criteria 1d-9
        ;;   :conv-steps 5000
        ;;   ;; :sub-conv-steps 1000
        ;;   :post-iter-step (lambda (i o e)
        ;;                     (push total-step step)
        ;;                     (incf total-step)
        ;;                     (push e res))
        ;;   :save-vtk-dr nil
        ;;   :save-vtk-loadstep t
        ;;   :dt-scale 1d0))
        ))))


(defun printstiff () 
  (setup :mps 2 :refine 0.25 :multigrid-refine 0)
  (change-class *sim* 'cl-mpm/implicit::mpm-sim-implicit)
  (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
        (cl-mpm::sim-ghost-factor *sim*) nil)
  (cl-mpm::update-sim *sim*)
  (cl-mpm::reset-node-displacement *sim*)
  (format t "printing stiffness ~%")
  (loop for vals in (cl-mpm/utils::sparse-to-coordinates
                     ;; (cl-mpm/implicit::assemble-stiffness-forward-difference *sim*)
                     ;; (cl-mpm/implicit::sim-global-k *sim*)
                     ;; (cl-mpm/implicit::assemble-stiffness-sparse *sim*)
                     ;; (cl-mpm/implicit::assemble-sparse-e *sim*)
                     (cl-mpm/utils::mat-to-sparse (cl-mpm/implicit::assemble-implicit-e *sim*))
                     ;; (cl-mpm/aggregate::sim-global-sparse-e *sim*)
                     )
        do (format t "~a ~a ~a~%" (nth 0 vals)
                   (nth 1 vals)
                   (nth 2 vals))))


;; (cl-mpm::iterate-over-nodes-serial
;;  (cl-mpm:sim-mesh *sim*)
;;  (lambda (n)
;;    (when (cl-mpm/mesh::node-active n)
;;      (format t "~a ~a~%" (cl-mpm/mesh::node-stiffness-fd n) (cl-mpm/mesh::node-agg-fd n)))))

;; (let ((bcs (cl-mpm/implicit::assemble-global-bcs *sim*)))
;;   (format t "bcs ~%")
;;   (loop for n across (cl-mpm/implicit::sim-nodes-fd *sim*)
;;         for i from 0
;;         do (format t "~a ~a ~a~%" (cl-mpm/mesh:node-index n)
;;                    (= 0d0 (cl-mpm/utils:varef bcs (+ 0 (* i 2))))
;;                    (= 0d0 (cl-mpm/utils:varef bcs (+ 1 (* i 2)))))))

;; (let ((bcs (cl-mpm/implicit::assemble-internal-bcs *sim*)))
;;   (format t "bcs ~%")
;;   (loop for n across (cl-mpm/aggregate::sim-agg-nodes-fdc *sim*)
;;         for i from 0
;;         do (format t "~a ~a ~a~%" (cl-mpm/mesh:node-index n)
;;                    (= 0d0 (cl-mpm/utils:varef bcs (+ 0 (* i 2))))
;;                    (= 0d0 (cl-mpm/utils:varef bcs (+ 1 (* i 2)))))))

(defun test-load-control-implicit ()
  (dolist (agg (list t))
    (dolist (r (list 1))
      (setup :mps 3 :refine r :multigrid-refine 0)
      ;; (cl-mpm::domain-sort-mps *sim*)
      (change-class *sim* 'cl-mpm/implicit::mpm-sim-implicit)
      ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static)
      ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-dr-usf)
      ;; (setf (cl-mpm::sim-velocity-algorithm *sim*) :tblend)
      (if agg
          (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
                (cl-mpm::sim-ghost-factor *sim*) nil)
          (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) nil
                (cl-mpm::sim-ghost-factor *sim*) nil))
      (cl-mpm/setup::set-mass-filter *sim* *density* :proportion 1d-15)
      (setf (cl-mpm::sim-gravity *sim*) -10d0)
      (let ((step (list))
            (res (list))
            (total-step 0))
        (vgplot:close-all-plots)
        (time
         (cl-mpm/dynamic-relaxation::run-load-control
          *sim*
          :output-dir (format nil "./output-nr-~a-~d/" (if agg "agg" "noagg") r)
          ;; :plotter (lambda (sim) (vgplot:semilogy (reverse step) (reverse res)))
          :plotter #'plot
          :load-steps 30
          :damping (sqrt 2)
          :substeps 1;(round (* 20 r))
          :criteria 1d-9
          :conv-steps 5000
          ;; :sub-conv-steps 1000
          :post-iter-step (lambda (i o e)
                            (push total-step step)
                            (incf total-step)
                            (push e res))
          :save-vtk-dr t
          :save-vtk-loadstep t
          :dt-scale 1d0))
        ))))


(defun test-3d ()
  (setup :mps 3 :refine 1)
  (cl-mpm/setup::set-mass-filter *sim* *density* :proportion 1d-9)

  (setf (cl-mpm::sim-ghost-factor *sim*) (* 1d6 1d-4))
  (setf (cl-mpm::sim-gravity *sim*) -1000d0)
  (vgplot:close-all-plots)
  (let ((step (list))
        (res (list))
        (total-step 0))
    (cl-mpm/dynamic-relaxation::run-load-control
     *sim*
     :output-dir (format nil "./output/")
     :load-steps 10
     :damping (sqrt 2)
     :substeps 10
     :criteria 1d-9
     :save-vtk-dr t
     :save-vtk-loadstep nil
     :plotter (lambda (sim) (vgplot:semilogy (reverse step) (reverse res)))
     :post-iter-step (lambda (i o e)
                       (push total-step step)
                       (incf total-step)
                       (push e res))

     :dt-scale 0.5d0)))

(defun test-3d-real-time ()
  (setup :mps 2 :refine 1)
  (let ((output-dir (format nil "./output-tflip/")))
    (vgplot:close-all-plots)
    (change-class *sim* 'cl-mpm/aggregate::mpm-sim-agg-usf)
    (setf (cl-mpm::sim-velocity-algorithm *sim*) :tflip)
    (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) nil)
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (let ((step 0))
      (cl-mpm/dynamic-relaxation::run-time
       *sim*
       :output-dir output-dir
       :plotter (lambda (sim)
                  (plot *sim*)
                  (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step)) :terminal "png size 1920,1080")
                  (incf step))
       :total-time 100d0
       :damping 1d-2
       :dt 0.1d0
       :initial-quasi-static nil
       :dt-scale 0.9d0))))






(defun get-damage ()
  (lparallel:pmap-reduce
   (lambda (mp)
     (if (typep mp 'cl-mpm/particle::particle-damage)
         (*
          (cl-mpm/particle::mp-mass mp)
          (cl-mpm/particle:mp-damage mp))
         0d0))
   #'+
   (cl-mpm:sim-mps *sim*)))
(defun get-plastic ()
  (lparallel:pmap-reduce
   (lambda (mp)
     (if (typep mp 'cl-mpm/particle::particle-plastic)
         (*
          (cl-mpm/particle::mp-mass mp)
          (cl-mpm/particle::mp-strain-plastic-vm mp))
         0d0))
   #'+
   (cl-mpm:sim-mps *sim*)))


(defun save-test-vtks (&key (output-dir "./output/"))
  (cl-mpm/output:save-vtk (merge-pathnames "test.vtk" output-dir) *sim*)
  (cl-mpm/output:save-vtk-nodes (merge-pathnames "test_nodes_0.vtk" output-dir) *sim*)
  (cl-mpm/output:save-vtk-nodes (merge-pathnames "test_nodes_1.vtk" output-dir) *sim*)
  (cl-mpm/output:save-vtk-cells (merge-pathnames "test_cells.vtk" output-dir) *sim*))

(require 'sb-sprof)

(defun spy (mat)
  (let ((pos-x (list))
        (pos-y (list)))
    (vgplot:close-all-plots)
    (vgplot:figure)
    (destructuring-bind (lx ly) (magicl:shape mat)
      (loop for x from 0 below lx
            do (loop for y from 0 below ly
                     do (when (> (abs (cl-mpm/utils::mtref mat x y)) 0d0)
                          (push x pos-x)
                          (push y pos-y)
                          ;; (push (- ly (+ 1 y)) pos-y)
                          )))
      (let ((ad 0.25d0))
        (vgplot:format-plot t "set xrange [~f:~f]" (- ad) (- (+ lx ad) 1d0))
        (vgplot:format-plot t "set yrange [~f:~f]" (- ad) (- (+ ly ad) 1d0)))
      (vgplot:plot pos-x pos-y ";;with points"))))


(defun setup-ghost (&key (refine 1) (mps 2)
                (sim-type 'cl-mpm:mpm-sim-usf)
                (multigrid-refine 0))
  (defparameter *sim* nil)
  (let ((mps-per-dim mps))
    (setf *sim* (setup-test-column '(4 1) '(2 1) sim-type refine mps-per-dim multigrid-refine)))
  ;; (cl-mpm/setup::initialise-stress-self-weight
  ;;  *sim*
  ;;  15d0
  ;;  )
  (format t "mps: ~d~%" (length (cl-mpm:sim-mps *sim*)))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))



(defun ghost-test ()
  (setup-ghost :refine 1)
  (setf (cl-mpm::sim-gravity *sim*) 0d0)
  (cl-mpm::iterate-over-mps
   (cl-mpm:sim-mps *sim*)
   (lambda (mp)
     (incf (cl-mpm/utils:varef (cl-mpm/particle::mp-position mp) 0) 1d-5)))
  (setf (cl-mpm::sim-ghost-factor *sim*) nil)
  (cl-mpm:update-sim *sim*)
  (dolist (x (list 2))
    (dolist (y (list 0 1))
      (let ((node (cl-mpm/mesh::get-node (cl-mpm:sim-mesh *sim*) (list x y 0))))
        (setf (cl-mpm/utils:varef (cl-mpm/mesh::node-displacment node) 1) 1d0))))
  (cl-mpm/ghost:apply-ghost *sim* 1d0)
  (cl-mpm::iterate-over-nodes-serial
   (cl-mpm:sim-mesh *sim*)
   (lambda (node)
     (when (cl-mpm/mesh::node-active node)
       (format t "node ~a ~e ~e~%" (cl-mpm/mesh:node-index node)
               (cl-mpm/utils:varef (cl-mpm/mesh::node-ghost-force node) 0)
               (cl-mpm/utils:varef (cl-mpm/mesh::node-ghost-force node) 1))))))

(defun mat-to-sparse (mat)
  (let ((v (make-array 0 :adjustable t :fill-pointer 0 :element-type 'double-float))
        (r (make-array 0 :adjustable t :fill-pointer 0 :element-type 'fixnum))
        (c (make-array 0 :adjustable t :fill-pointer 0 :element-type 'fixnum)))
    (loop for x below (magicl:nrows mat) do
      (loop for y below (magicl:ncols mat)
            when (> (abs (cl-mpm/utils::mtref mat x y)) 0d0)
              do (progn
                   (vector-push-extend (cl-mpm/utils::mtref mat x y) v)
                   (vector-push-extend x r)
                   (vector-push-extend y c))))
    (values v r c)))

(defun test-sparse-mat ()
  (let* ((m (cl-mpm/utils:matrix-from-list (list 0d0 0d0 2d0
                                                 0d0 0d0 0d0
                                                 2d0 0d0 1d0)))
         (bcs (cl-mpm/utils::vector-from-list (list 1d0 0d0 1d0)))
         )
    (multiple-value-bind (v r c) (mat-to-sparse m)
      ;; (pprint v)
      ;; (pprint r)
      ;; (pprint c)
      (let ((sm (cl-mpm/utils::build-sparse-matrix v r c 3 3))
            (x (cl-mpm/utils:vector-from-list (list 5d0 1d0 2d0))))
        (pprint m)
        (pprint sm)
        (pprint (cl-mpm/utils::sparse-to-mat sm))
        ;; (pprint (magicl:@ sm (magicl:vector-from-list (list 1d0 1d0 1d0))))
        ;; (pprint (magicl:@ m x))
        ;; (pprint (cl-mpm/fastmaths::fast-@-sparse-mat-dense-vec-masked sm x bcs bcs))
        (pprint (cl-mpm/implicit::linear-solve-with-bcs m x bcs))
        (pprint (cl-mpm/implicit::sparse-linear-solve-with-bcs sm x bcs))
        ))))

(defun test-summing ()
  (let ((v #(1d0 2d0 1d0 1d0 1d0))
        (r #(3 3 3 0 0))
        (c #(0 0 0 3 3)))
    (let ((ms (cl-mpm/utils::build-sparse-matrix v r c 4 4)))
      (pprint (cl-mpm/utils::sparse-to-mat ms))))
  )
;; (let ((a (cl-mpm/utils:vector-zeros))
;;       (b (cl-mpm/utils:vector-zeros)))
;;   (time
;;    (let ((as (cl-mpm/utils::fast-storage a))
;;          (bs (cl-mpm/utils::fast-storage b)))
;;      (dotimes (i 100000)
;;        (loop for i from 0 to 2
;;              do (setf (aref as i) (aref bs i)))
;;        ;; (cl-mpm/utils::copy-into a b)
;;        ))

;;    ))


;; (let ((a (cl-mpm/utils::matrix-eye 1d0))
;;       (out (cl-mpm/utils::matrix-zeros)))
;;   (time-form 100000000 (cl-mpm/fastmaths::fast-inv-3x3 a out)))

;; (let ((a (cl-mpm/utils::matrix-from-list (list 1d0 2d0 3d0 2d0 5d0 6d0 7d0 8d0 9d0)))
;;       (out (cl-mpm/utils::matrix-zeros))
;;       )
;;   (pprint (magicl:inv a))
;;   (pprint (cl-mpm/fastmaths::fast-inv-3x3 a))
;;   (pprint (cl-mpm/fastmaths::fast-inv-3x3 a out)))

;; (dotimes (i 100000000)
;;   (let ((v (cl-mpm/utils:voigt-from-list (loop repeat 6 collect (* (random 1d0) (- (random 2d0) 1d0))))))
;;     (pprint v)
;;     (cl-mpm/fastmaths::eigenvalues-3x3 v))
;;   )


;; (loop for i from 0 below (expt 4 2)
;;       do (pprint (cl-mpm::morton-to-index i 2)))



;; (defun pwork (function &optional arguments)
;;   "distribute a task to each lparallel worker and wait until all finish.
;; arguments may be a vector of size equal or smaller than (kernel-worker-count).
;; all results are collected in a vector of the same size."
;;   (let* ((nr-workers (lparallel:kernel-worker-count))
;;          (arguments (or arguments
;;                         (make-array nr-workers :initial-element nil)))
;;          ;; in case fewer arguments than threads have been provided
;;          (n (min (length arguments) nr-workers))
;;          (channel (make-channel))
;;          (from-workers (make-queue))
;;          (to-workers (make-queue))
;;          (results (make-array n :initial-element nil)))
;;     (when (< n (length arguments))
;;       (error "more work than workers.  systematic error?"))
;;     (loop repeat n
;;           do (submit-task channel (lambda ()
;;                                     (push-queue t from-workers)
;;                                     (pop-queue to-workers)
;;                                     (cons *worker-id*
;;                                           (apply function
;;                                                  (aref arguments *worker-id*))))))
;;     (loop repeat n do (pop-queue from-workers))
;;     ;; at this point every worker should wait for work
;;     ;; and the next line sets them doing the work
;;     (loop repeat n do (push-queue t to-workers))
;;     ;; and collect the results
;;     (loop repeat n
;;           for (id . result) = (receive-result channel) do
;;             (setf (aref results id) result))
;;     results))

;; (progn
;;   (setup :refine 8)
;;   (let* ((value (cl-mpm/mesh::mesh-nodes (cl-mpm:sim-mesh *sim*)))
;;          (len (array-total-size value))
;;          (v (make-array (array-total-size value) :element-type 'cl-mpm/mesh::node :displaced-to value))
;;          (vs (make-array len :element-type 'cl-mpm/mesh::node))
;;          (iters 1000000))
;;     (declare ((vector cl-mpm/mesh::node *) v)
;;              ((simple-array cl-mpm/mesh::node) vs))
;;     (loop for va across v
;;           for i from 0
;;           do (setf (aref vs i) va))
;;     (time-form
;;      iters
;;      (dotimes (i len) (aref v i)))
;;     (time-form
;;      iters
;;      (dotimes (i len) (aref vs i)))
;;     (pprint (* iters len))
;;     ))



  ;; (let* ((e 1d0)
  ;;        (nu 0.2d0)
  ;;        (rho 1d0)
  ;;        (strain (cl-mpm/utils:voigt-from-list (list 0d0 0d0 0d0
  ;;                                                    0d0 0d0 2d0)))
  ;;        (de (cl-mpm/constitutive:linear-elastic-matrix e nu))
  ;;        (stress (magicl:@ de strain))
  ;;        )
  ;;   (cl-mpm/constitutive::plastic-vm-tangent stress de strain rho e nu)
  ;;   ;; (multiple-value-bind ())
  ;;   )

(defun make-stats-csv (output-dir filename)
  (with-open-file (stream (merge-pathnames filename output-dir) :direction :output :if-exists :supersede)
    (format stream "refine,lstps,time,error~%")))
(defun save-stats-csv (output-dir filename refine lstps time error)
  (with-open-file (stream (merge-pathnames filename output-dir) :direction :output :if-exists :append)
    (format stream "~d,~d,~e,~e~%" refine lstps (float time 0e0) (float error 0e0))))


(defun compare-conv-lstp ()
  (let ((agg t)
        (refine 1)
        (mps 4)
        (lstp-list '(1 4 8 16 32 48 64 96 128)))
    ;; (run-conv-lstp-implicit  :refine refine :agg agg :lstp-list lstp-list :mps mps)
    (run-conv-lstp :refine refine :agg agg :lstp-list lstp-list :mps mps)
    ))

(defun run-conv-lstp (&key (refine 4) (agg nil)
                           (mps 3) (lstp-list '(1 2 4 8 16 32 64)))
  (let ((data-output-dir "./examples/ample/collapse/")
        (filename "stats.csv"))
    (setf *run-sim* t)
    (make-stats-csv data-output-dir filename)
    (loop for i in
          lstp-list
          while *run-sim*
          do
             (let* ((elements (expt 2 refine)))

               (setup :refine refine :mps mps)
               (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) agg)
               (setf lparallel:*debug-tasks-p* nil)
               ;; (ignore-errors)
               (let ((start (get-internal-real-time)))
                 (cl-mpm/dynamic-relaxation::run-load-control
                  *sim*
                  :output-dir (merge-pathnames (format nil "./output-dr-~a_~f_~d/" i refine mps))
                  :load-steps i
                  :substeps (round (* 10 (expt 2 refine)))
                  :plotter (lambda (sim))
                  :damping (sqrt 2d0)
                  :save-vtk-dr nil
                  :save-vtk-loadstep t
                  :conv-steps 1000
                  :dt-scale 1d0
                  :criteria 1d-9)
                 (let* ((end (get-internal-real-time))
                        (units internal-time-units-per-second)
                        (dt (/ (- end start) units)))
                   (save-stats-csv data-output-dir filename refine i dt (compute-error *sim*))))))))

(defun compute-error (sim)
  0d0)

(defun run-conv-lstp-implicit (&key (refine 4) (agg nil)
                                 (mps 3)
                                 (lstp-list '(1 2 4 8 16 32 64)))
  (let ((data-output-dir "./examples/ample/collapse/")
        (filename "stats_nr.csv"))
    (setf *run-sim* t)
    (make-stats-csv data-output-dir filename)
    (loop for i in
          lstp-list
          while *run-sim*
          do
             (let* ((elements (expt 2 refine)))
               (setup :refine refine :mps mps)
               (change-class *sim* 'cl-mpm/implicit::mpm-sim-implicit)
               (cl-mpm:iterate-over-mps
                (cl-mpm:sim-mps *sim*)
                (lambda (mp)
                  (change-class mp 'cl-mpm/particle::particle-vm-implicit)))
               (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) agg)
               (setf lparallel:*debug-tasks-p* nil)
               (ignore-errors
                (let ((start (get-internal-real-time)))
                  (cl-mpm/dynamic-relaxation::run-load-control
                   *sim*
                   :output-dir (merge-pathnames (format nil "./output-nr-~a_~d/" i mps))
                   :load-steps i
                   :substeps 1
                   :plotter (lambda (sim))
                   :conv-steps 500
                   :damping (sqrt 2)
                   :save-vtk-dr nil
                   :save-vtk-loadstep t
                   :criteria 1d-9
                   :dt-scale 1d0)
                  (let* ((end (get-internal-real-time))
                         (units internal-time-units-per-second)
                         (dt (/ (- end start) units)))
                    (save-stats-csv data-output-dir filename refine i dt (compute-error *sim*)))))))))


(defun compute-max-extent (sim)
  (cl-mpm::reduce-over-mps
   (cl-mpm:sim-mps sim)
   (lambda (mp)
     (cl-mpm/utils:varef (cl-mpm/fastmaths:fast-.+
                          (cl-mpm/particle::mp-domain-size mp)
                          (cl-mpm/particle::mp-position mp)) 0))
   #'max))



(defun run-conv ()
  (let ((agg t)
        (mps 4)
        (refine-list (list 0.5 1 2 4 8 16)))
    (run-conv-disp :agg agg :refine-list refine-list :mps mps)
    ;; (run-conv-disp-implict :agg agg :refine-list refine-list :mps mps)
    ))

(defun run-conv-disp (&key (refine-list (list 1)) (agg nil)
                           (mps 3))
  (let ((data-output-dir "./examples/ample/collapse/")
        (filename "convergence.csv")
        (i 30)
        )
    (setf *run-sim* t)
    (make-stats-csv data-output-dir filename)
    (loop for refine in refine-list
          while *run-sim*
          do
             (let* ()
               (setup :refine refine :mps mps)
               (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) agg)
               (cl-mpm/setup::set-mass-filter *sim* *density* :proportion 1d-3)
               (setf (cl-mpm::sim-gravity *sim*) -10d0)
               (setf lparallel:*debug-tasks-p* nil)
               ;; (ignore-errors)
               (let ((start (get-internal-real-time)))
                 (cl-mpm/dynamic-relaxation::run-load-control
                  *sim*
                  :output-dir (merge-pathnames (format nil "./output-dr-~a_~f_~d/" i refine mps))
                  :load-steps i
                  :substeps (round (* 10 (expt 2 refine)))
                  :plotter (lambda (sim))
                  :damping (sqrt 2d0)
                  :save-vtk-dr t
                  :save-vtk-loadstep t
                  :conv-steps 1000
                  :dt-scale 2d0
                  :criteria 1d-9)
                 (let* ((end (get-internal-real-time))
                        (units internal-time-units-per-second)
                        (dt (/ (- end start) units)))
                   (save-stats-csv data-output-dir
                                   filename
                                   refine
                                   i
                                   dt
                                   (compute-max-extent *sim*)
                                   )))))))

(defun run-conv-disp-implict (&key (refine-list (list 1)) (agg nil)
                           (mps 3))
  (let ((data-output-dir "./examples/ample/collapse/")
        (filename "convergence.csv")
        (i 30)
        )
    (setf *run-sim* t)
    (make-stats-csv data-output-dir filename)
    (loop for refine in refine-list
          while *run-sim*
          do
             (let* ()
               (setup :refine refine :mps mps)
               (change-class *sim* 'cl-mpm/implicit::mpm-sim-implicit)
               (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) agg)
               (cl-mpm/setup::set-mass-filter *sim* *density* :proportion 1d-3)
               (setf (cl-mpm::sim-gravity *sim*) -10d0)
               (setf lparallel:*debug-tasks-p* nil)
               ;; (ignore-errors)
               (let ((start (get-internal-real-time)))
                 (cl-mpm/dynamic-relaxation::run-load-control
                  *sim*
                  :output-dir (merge-pathnames (format nil "./output-nr-~a_~f_~d/" i refine mps))
                  :load-steps i
                  :substeps 1
                  :plotter (lambda (sim))
                  :damping (sqrt 2d0)
                  :save-vtk-dr t
                  :save-vtk-loadstep t
                  :conv-steps 1000
                  :dt-scale 1d0
                  :criteria 1d-9)
                 (let* ((end (get-internal-real-time))
                        (units internal-time-units-per-second)
                        (dt (/ (- end start) units)))
                   (save-stats-csv data-output-dir
                                   filename
                                   refine
                                   i
                                   dt
                                   (compute-max-extent *sim*)
                                   )))))))
