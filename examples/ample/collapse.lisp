(defpackage :cl-mpm/examples/collapse
  (:use :cl))

(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)
;; (sb-ext:restrict-compiler-policy 'speed  3 3)
;; (sb-ext:restrict-compiler-policy 'debug  0 0)
;; (sb-ext:restrict-compiler-policy 'safety 0 0)
;; (setf *block-compile-default* nil)
;(sb-int:set-floating-point-modes :traps '(:overflow :invalid :inexact :divide-by-zero :underflow))
;; (sb-int:set-floating-point-modes :traps '(:overflow :divide-by-zero :underflow))

(in-package :cl-mpm/examples/collapse)
;; (declaim (optimize (debug 0) (safety 0) (speed 3)))
(declaim (optimize (debug 3) (safety 3) (speed 0)))

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-elastic) dt fbar)
  (cl-mpm::update-stress-kirchoff-dynamic-relaxation mesh mp dt fbar))

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-elastic) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  (cl-mpm::update-domain-polar-2d mesh mp dt))

(defun plot-load-disp ()
  (vgplot:semilogy *data-steps* *data-energy*))
(declaim (notinline plot))
(defun plot (sim)
  ;; (cl-mpm::g2p (cl-mpm:sim-mesh *sim*)
  ;;              (cl-mpm:sim-mps *sim*)
  ;;              (cl-mpm:sim-dt *sim*)
  ;;              0d0
  ;;              :TRIAL)
  ;; (plot-load-disp)
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
         (h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
         (ms-x (first ms))
         (ms-y (second ms))
        )
    ;; (vgplot:format-plot t "set object 1 rect from ~f,~f to ~f,~f fc rgb 'black' fs transparent solid 0.5 noborder behind" 0 ms-y ms-x (- datum))
    )
  ;; (format t "~E~%" (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :deformed
   :trial nil
   :colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :xy))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-damage-ybar mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-strain-plastic-vm mp))
   )
  )

(defparameter *eta* 1d0)

(defun stop ()
  (setf *run-sim* nil))

(declaim (notinline setup-test-column))
(defun setup-test-column (size block-size sim-type &optional (e-scale 1) (mp-scale 1) (multigrid-refinement 0))
  (let* (;; (multigrid-refinement 3)
         (multigrid-enabled (> multigrid-refinement 0))
         ;; (e-scale (if multigrid-enabled (/ e-scale (expt 2 (- multigrid-refinement 1)))
         ;;              e-scale))
         ;; (mp-scale (if multigrid-enabled (* mp-scale (expt 2 (- multigrid-refinement 1)))
         ;;               mp-scale))
         (E 1d6)
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
                 :split-factor (* 2d0 (/ 1d0 mp-scale))
                 :enable-fbar t
                 :enable-aggregate nil
                 :ghost-factor (* E 1d-4)
                 :max-split-depth 6
                 :enable-split t
                 :gravity -10d0
                 )
                (when multigrid-enabled
                  (list :refinement multigrid-refinement))
                )))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 1d3)
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
          ;; 'cl-mpm/particle::particle-elastic
          ;; :E 1d6
          ;; :nu 0.30d0
          'cl-mpm/particle::particle-vm
          :E E
          :nu 0.3d0
          :rho 20d3
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

      (let ((dt-scale 1d0))
        (setf
         (cl-mpm:sim-dt sim)
         (* dt-scale h
            (sqrt (cl-mpm::sim-mass-scale sim))
            (sqrt (/ density (cl-mpm/particle::mp-p-modulus (aref (cl-mpm:sim-mps sim) 0)))))))

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
    (setf *sim* (setup-test-column '(32 16) '(8 8) sim-type refine mps-per-dim multigrid-refine)))
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

(defun test-real-time ()
  (setup :mps 2)
  (let ((output-dir (format nil "./output-0.25/")))
    (vgplot:close-all-plots)
    ;; (change-class *sim* 'cl-mpm/aggregate::mpm-sim-agg-usf)
    (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic)
    (setf (cl-mpm::sim-velocity-algorithm *sim*) :BLEND)
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (let ((step 0))
      (cl-mpm/dynamic-relaxation::run-time
       *sim*
       :output-dir output-dir
       :plotter (lambda (sim)
                  (plot *sim*)
                  (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step)) :terminal "png size 1920,1080")
                  (incf step)
                  )
       :total-time 100d0
       :damping 1d-3
       :dt 1d0
       :initial-quasi-static nil
       :dt-scale 5d0))))

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

(defun test-load-control ()
  (setup :mps 3 :refine 1 :multigrid-refine 0)
  (cl-mpm/setup::set-mass-filter *sim* *density* :proportion 1d-15)
  ;; (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t)
  (setf (cl-mpm::sim-gravity *sim*) -10d0)
  ;; (setf (cl-mpm::sim-gravity *sim*) -1000d0)
  (cl-mpm/dynamic-relaxation::run-load-control
   *sim*
   :output-dir (format nil "./output/")
   :plotter #'plot
   :load-steps 50
   :damping 1d0;(sqrt 2)
   :substeps 10
   :criteria 1d-3
   :save-vtk-dr t
   :save-vtk-loadstep t
   :dt-scale 1d0))


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
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
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
  ;; (dolist (x (list 0 1))
  ;;   (dolist (y (list 0 1))
  ;;     (let ((node (cl-mpm/mesh::get-node (cl-mpm:sim-mesh *sim*) (list x y 0))))
  ;;       (setf (cl-mpm/utils:varef (cl-mpm/mesh::node-displacment node) 0) 1d0)
  ;;       )))
  (dolist (x (list 2))
    (dolist (y (list 0 1))
      (let ((node (cl-mpm/mesh::get-node (cl-mpm:sim-mesh *sim*) (list x y 0))))
        (setf (cl-mpm/utils:varef (cl-mpm/mesh::node-displacment node) 1) 1d0))))
  (cl-mpm/ghost:apply-ghost *sim* 1d0)
  (cl-mpm::iterate-over-nodes-serial
   (cl-mpm:sim-mesh *sim*)
   (lambda (node)
     (when (cl-mpm/mesh::node-active node)
       (format t "Node ~A ~E ~E~%" (cl-mpm/mesh:node-index node)
               (cl-mpm/utils:varef (cl-mpm/mesh::node-ghost-force node) 0)
               (cl-mpm/utils:varef (cl-mpm/mesh::node-ghost-force node) 1))))
   )
  )

;; (let* ((nx 1d0)
;;        (ny 0d0)
;;        (nz 0d0)
;;        (dsvp-adjuster (magicl:transpose! (cl-mpm/utils::arb-matrix-from-list
;;                                           (list
;;                                            nx  0d0 0d0
;;                                            0d0 ny  0d0
;;                                            0d0 0d0 nz
;;                                            0d0 nz  ny
;;                                            nz  0d0 nx
;;                                            ny  nx  0d0)
;;                                           3
;;                                           6)))
;;        (gradient-terms (cl-mpm/utils:vector-zeros))
;;        (grads (list 1d0 2d0 3d0))
;;        (dsvp
;;          (cl-mpm/fastmaths::fast-scale!
;;           (magicl:transpose (cl-mpm/shape-function::assemble-dsvp-3d grads))
;;           1d0)))
;;   (pprint dsvp)
;;   (pprint dsvp-adjuster)
;;   (pprint (magicl:@ dsvp dsvp-adjuster))
;;   )
