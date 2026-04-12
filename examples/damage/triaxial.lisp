(defpackage :cl-mpm/examples/damage/triaxial
  (:use :cl
   :cl-mpm/example
        :cl-mpm/utils))
(in-package :cl-mpm/examples/damage/triaxial)

(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)
(setf cl-mpm/settings::*optimise-setting* cl-mpm/settings::*optimise-speed*)

(declaim (notinline plot-domain))


(defun plot-domain ()
  (when *sim*
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :trial t
     :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp)))))

(defun get-load ()
  (* 2d0 (cl-mpm/penalty::resolve-load *penalty*)))

(declaim (notinline setup))
(defun setup (&key
                (refine 1)
                (mps 3)
                (enable-fbar t)
                (multigrid-refines 0)
                (gf 10d0)
                (kt 1d0)
                (kc 0d0)
                (angle 40d0)
                (angle-r 0d0)
                (oversize-factor (- 1d0 1d-2))
                )
  (let* ((width 80d-3)
         (height 170d-3)
         (h (/ 0.01d0 refine))
         (offset width)
         (domain-height (* height 1.25d0))
         (density 1d3)
         (E 1.3d7)
         (nu 0.30d0)
         (domain-size (list (* 3 width)
                            domain-height
                            (* 3 width)))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list width height width)))
    (format t "Mesh size ~E~%" h)
    (setf
     *sim*
     (cl-mpm/setup::make-simple-sim
      h
      element-count
      :sim-type
       'cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul
      ;; 'cl-mpm/dynamic-relaxation::mpm-sim-dr-multigrid
      :args-list
      (list
       ;; :enable-aggregate t
       ;; :ghost-factor nil
       ;; :mass-update-count 1
       ;; :enable-aggregate nil
       :ghost-factor (* E 1d-2)
       :enable-aggregate nil
       ;; :ghost-factor nil
       :enable-split t
       :max-split-depth 2
       :enable-fbar enable-fbar
       ;; :refinement multigrid-refines
       )))

    (setf h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
    (let* ((init-stress (cl-mpm/damage::mohr-coloumb-coheasion-to-tensile 40d3 angle))
           ;; (L 10d-3)
           (L h)
           (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf L init-stress E)))
      (format t "Ductility ~E~%" ductility)
      (assert (> ductility 1d0))
      (pprint block-size)
      (pprint (mapcar (lambda (e) (* (/ e h) mps)) block-size))
      (cl-mpm:add-mps
       *sim*
       (cl-mpm/setup:make-block-mps
        (list offset 0d0 offset)
        block-size
        (mapcar (lambda (e) (* (/ e h) mps)) block-size)
        density
        ;; 'cl-mpm/particle::particle-damage-frictional
        ;; :E E
        ;; :nu nu
        ;; :local-length L
        ;; :ductility ductility
        ;; :friction-angle (cl-mpm/utils:deg-to-rad angle)
        ;; :initiation-stress init-stress
        ;; :friction-model :MC
        ;; :kt-res-ratio 1d0
        ;; :kc-res-ratio 0d0
        ;; :g-res-ratio 1d0

        'cl-mpm/particle::particle-plastic-damage-frictional
        :E E
        :nu nu
        :local-length L
        :ductility ductility
        :friction-angle (cl-mpm/utils:deg-to-rad angle)
        :residual-friction (cl-mpm/utils:deg-to-rad angle-r)
        :initiation-stress init-stress
        :friction-model :DP
        :oversize oversize-factor
        :kt-res-ratio kt
        :kc-res-ratio kc
        :psi 0d0
        :plastic-damage-evolution t


        ;; 'cl-mpm/particle::particle-vm
        ;; :E E
        ;; :nu nu
        ;; :rho (* 1d6 (sqrt 3/2))
        ;; :rho-r (* 1d6 (sqrt 3/2))
        ;; :softening 10d0
        ;; :rho 1d6
        ;; :rho 1d6
        ;; 'cl-mpm/particle::particle-mc
        ;; :E E
        ;; :nu nu
        ;; :phi (cl-mpm/utils:deg-to-rad 30d0)
        ;; :psi (cl-mpm/utils:deg-to-rad 5d0)
        ;; :c 1d5
        )
       ))
    (let ((angle 45d0)
          (scale 4d0)
          )
      (cl-mpm/setup::remove-sdf
       *sim*
       (lambda (p)
         (-
          (funcall
           (cl-mpm/setup::circle-sdf
            (list (* width 1.5d0)
                  0d0
                  (* width 1.5d0)
                  )
            (* 0.5d0 width))
           (cl-mpm/fastmaths:fast-.*
            p
            (cl-mpm/utils:vector-from-list (list 1d0 0d0 1d0))))))
       :refine 0)

      ;; (cl-mpm/setup::apply-sdf
      ;;  *sim*
      ;;  (cl-mpm/setup::-sdf
      ;;   (list (* width 1.5d0) (* height 0.5d0) (* width 1.5d0))
      ;;   10d-3
      ;;   1d-3)
      ;;  (lambda (mp sdf)
      ;;    ;; (setf (cl-mpm/particle::mp-damage mp) sdf)
      ;;    (when (<= sdf 0d0)
      ;;      ;; (pprint sdf)
      ;;      ;; (setf (cl-mpm/particle::mp-damage mp) sdf)
      ;;      ;(cl-mpm/damage::set-mp-damage mp sdf)
      ;;      (cl-mpm/damage::set-mp-damage mp 0.99d0)
      ;;      )
      ;;    ))
      (cl-mpm/setup::apply-sdf
       *sim*
       (cl-mpm/setup::circle-sdf
        (list (* width 1.5d0) (* height 0.5d0) (* width 1.5d0))
        (* 4d-3 scale))
       (lambda (mp sdf)
         (when (<= sdf 0d0)
           (cl-mpm/damage::set-mp-damage mp 0.99d0))))
      ;; (cl-mpm/setup::apply-sdf
      ;;  *sim*
      ;;  (cl-mpm/setup::ellipse-sdf
      ;;   (list (* width 1.5d0) (* height 0.5d0) 0d0)
      ;;   (* 4d-3 scale)
      ;;   (* 4d-3 scale))
      ;;   (lambda (mp sdf)
      ;;     (when (<= sdf 0d0)
      ;;       (pprint sdf)
      ;;       (cl-mpm/damage::set-mp-damage mp 0.99d0))))
      )


    (setf (cl-mpm::sim-gravity *sim*) 0d0)
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)
    (cl-mpm/setup::setup-bcs
     *sim*
     ;; :left '(0 nil nil)
     ;; :right '(0 nil nil)
     :bottom '(nil 0 nil))
    (cl-mpm::add-bcs
     *sim*
     (cl-mpm/bc::make-bc-fixed
      (list (round offset h) 0 (round (+ (* 0.5d0 width) offset) h))
      '(0 0 0)))

    (let* ((friction 0d0)
           (epsilon-scale 1d1)
           (normal (cl-mpm/utils:vector-from-list '(0d0 -1d0 0d0)))
           (epsilon (* (cl-mpm/particle::calculate-p-wave-modulus E nu) epsilon-scale))
           )
      ;; (defparameter *penalty*
      ;;   (cl-mpm/penalty::make-bc-penalty-distance-point
      ;;    *sim*
      ;;    normal
      ;;    (cl-mpm/utils:vector-from-list (list (+ offset (* 0.5 width)) height 0d0))
      ;;    width
      ;;    epsilon
      ;;    0d0
      ;;    0d0
      ;;    ;; :clip-function (lambda (p) (cl-mpm/penalty::clip-radial
      ;;    ;;                             p
      ;;    ;;                             normal
      ;;    ;;                             (cl-mpm/utils:vector-from-list (list offset height 0d0))
      ;;    ;;                             penalty-width)
      ;;    ))
      (defparameter *penalty*
        (cl-mpm/penalty::make-bc-penalty-displacment
         *sim*
         (cl-mpm/fastmaths:fast-scale! normal -1d0)
         epsilon
         ;; :clip-function (lambda (p) (cl-mpm/penalty::clip-radial
         ;;                             p
         ;;                             normal
         ;;                             (cl-mpm/utils:vector-from-list (list offset height 0d0))
         ;;                             penalty-width)
         ))
      )
    (cl-mpm:add-bcs-force-list
     *sim*
     *penalty*)
    )
  (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
  )

(defun output-disp-header (output-dir)
  (with-open-file (stream (merge-pathnames "disp.csv" output-dir) :direction :output :if-exists :supersede)
    (format stream "disp,load~%")))

(defun output-disp-data (output-dir)
  (with-open-file (stream (merge-pathnames "disp.csv" output-dir) :direction :output :if-exists :append)
    (format stream "~f,~f~%"
            *displacement*
            (get-load))))

(defun plot-load-disp ()
  (vgplot:plot
   (mapcar (lambda (x) (* x -1d3)) *data-disp*)
   (mapcar (lambda (x) x) *data-load*)
   "Solution"
   )
  (vgplot:xlabel "Displacement (mm)")
  )

(defun save-csv (output-dir filename data-disp data-load)
  (with-open-file (stream (merge-pathnames filename output-dir) :direction :output :if-exists :supersede)
    (format stream "disp,load~%")
    (loop for disp in data-disp
          for load in data-load
          do (format stream "~E,~E~%" (float disp 0e0) (float load 0e0)))))
(declaim (notinline run))
(defun run (&key (output-dir (format nil "./output/"))
              (refine 1)
              (tensile nil)
              (enable-plastic nil)
              (csv-dir nil)
              (csv-filename (format nil "load-disp.csv")))
  (unless csv-dir
    (setf csv-dir output-dir))
  (let* ((lstps 20)
         (total-disp -10d-3)
         (current-disp 0d0)
         (step 0))
    (when tensile
      (setf total-disp (abs total-disp)))
    (defparameter *data-disp* (list 0d0))
    (defparameter *data-load* (list 0d0))
    (defparameter *displacement* 0d0)
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))

    (vgplot:close-all-plots)
    (time
     (cl-mpm/dynamic-relaxation::run-adaptive-load-control
      *sim*
      :output-dir output-dir
      :plotter (lambda (sim)
                 (plot-domain))
      :loading-function (lambda (i)
                          (setf current-disp (* i total-disp))
                          (cl-mpm/penalty::bc-set-displacement
                           *penalty*
                           (cl-mpm/utils:vector-from-list (list 0d0 current-disp 0d0))))
      :post-conv-step (lambda (sim)
                        (push current-disp *data-disp*)
                        (let ((load (get-load)))
                          (format t "Load ~E~%" load)
                          (push load *data-load*))
                        ;; (plot-load-disp)
                        (save-csv csv-dir csv-filename *data-disp* *data-load*)
                        ;; (output-disp-data output-dir)
                        (incf step))
      :load-steps lstps
      :enable-plastic enable-plastic
      :enable-damage t
      :damping (sqrt 1d0)
      :min-adaptive-steps 0
      :max-adaptive-steps 3
      :adaption-constant 4
      :max-damage-inc 0.50d0
      :min-damage-inc 0.05d0
      :substeps (round (* refine 100))
      :criteria 1d-3
      :save-vtk-dr t
      :save-vtk-loadstep t
      :dt-scale 1d0))))

(defun test ()
  (dolist (refine (list 1))
    (setup :mps 3 :refine refine :enable-fbar t :multigrid-refines 0
           :angle 30d0
           :angle-r 0d0
           :gf 60d0
           )
    (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
    ;; (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
    (run :output-dir (format nil "./output-agg-~D/" refine)
         :enable-plastic nil
         :refine refine))
  ;; (save-csv "./examples/fbar/rigid-footing/" (format nil "data_fbaradjust_~A.csv" t) *data-disp* *data-load*)
  ;; (dolist (fbar (list t nil))
  ;;   (save-csv "./examples/fbar/rigid-footing/" (format nil "data_fbar_~A.csv" fbar) *data-disp* *data-load*)
  ;;   )
  )

(defun test-oversize ()
  (let ((refine 1)
        (mps 3))
    (dolist (oversize (list
                       (- 1d0 1d-1)
                       ;; (- 1d0 1d-2)
                       ;; (- 1d0 1d-3)
                       ))
      (setup :mps mps
             :refine refine
             :enable-fbar t
             :multigrid-refines 0
             :gf 100d0
             :angle 30d0
             :angle-r 10d0
             :oversize-factor oversize)
      ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
      (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
      (run :output-dir (format nil "./output-pdinc-tensile-~F/" oversize)
           :refine refine
           :tensile t))))

(defun test-tension-compression ()
  (let ((refine 1)
        (mps 2))
    (dolist (tension (list
                      nil
                      t
                       ))
      (setup :mps mps
             :refine refine
             :enable-fbar t
             :multigrid-refines 0
             :gf 100d0
             :angle 15d0
             :angle-r 0d0
             :kt (- 1d0 1d-6)
             )
      ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
      (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
      (run :output-dir (format nil "./output-isotropic-~D-~A/" refine (if tension "tension" "compression"))
           :refine refine
           :tensile tension
           :enable-plastic nil
           )))
  )
;; (defun profile ()
;;   (setup :refine 0.5 :mps 2)
;;   (cl-mpm::update-sim *sim*)
;;   (sb-profile:profile "CL-MPM")
;;   (sb-profile:profile "CL-MPM/PARTICLE")
;;   (sb-profile:profile "CL-MPM/AGGREGATE")
;;   (sb-profile:profile "CL-MPM/MESH")
;;   (sb-profile:profile "CL-MPM/SHAPE-FUNCTION")
;;   (sb-profile:profile "CL-MPM/DAMAGE")
;;   (sb-profile:reset)
;;   (setf (cl-mpm::sim-enable-damage *sim*) nil)
;;   (time
;;    (dotimes (i 10)
;;      (cl-mpm::update-sim *sim*)))
;;   (sb-profile:report)
;;   )
