(defpackage :cl-mpm/examples/damage/biaxial
  (:use :cl
   :cl-mpm/example
        :cl-mpm/utils))
(in-package :cl-mpm/examples/damage/biaxial)

;; (sb-ext:restrict-compiler-policy 'speed  3 3)
;; (sb-ext:restrict-compiler-policy 'debug  0 0)
;; (sb-ext:restrict-compiler-policy 'safety 0 0)
(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)
;; (defmethod cl-mpm/dynamic-relaxation::damage-increment-criteria ((sim cl-mpm/dynamic-relaxation::mpm-sim-dr-ul))
;;   (cl-mpm/dynamic-relaxation::damage-increment-criteria sim))

(declaim (notinline plot-domain))


(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-plastic-damage-frictional) dt)
  (with-accessors ((strain cl-mpm/particle::mp-strain)
                   (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                   (E cl-mpm/particle::mp-e)
                   (de cl-mpm/particle::mp-elastic-matrix)
                   (y cl-mpm/particle::mp-damage-y-local)
                   (angle cl-mpm/particle::mp-friction-angle)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (model cl-mpm/particle::mp-friction-model)
                   (pd-inc cl-mpm/particle::mp-plastic-damage-evolution)
                   (ps-vm cl-mpm/particle::mp-strain-plastic-vm)
                   )
      mp
    (let ((stress
            undamaged-stress
            ;; (cl-mpm/fastmaths:fast-scale-voigt undamaged-stress (/ 1d0 (cl-mpm/fastmaths:det-3x3 def)))
                  )
          (ps-y (sqrt (* E (expt ps-vm 2))))
          )
      (setf
       y
       (+
        (if pd-inc ps-y 0d0)
        (ecase model
          (:SE (cl-mpm/damage::tensile-energy-norm strain e de))
          (:RANKINE (cl-mpm/damage::criterion-rankine stress))
          (:RANKINE-SMOOTH (cl-mpm/damage::criterion-rankine-smooth stress))
          (:PRINC-STRAIN (cl-mpm/damage::criterion-max-principal-strain strain e))
          (:PRINC-STRESS (cl-mpm/damage::criterion-max-principal-stress stress))
          (:DP (cl-mpm/damage::drucker-prager-criterion stress angle))
          (:MC (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile stress angle))
          (:MCR (cl-mpm/damage::criterion-mohr-coloumb-rankine-stress-tensile stress angle))
          ))))))


(defun plot-domain ()
  (when *sim*
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :trial t
     :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp)))))

(defun get-load ()
  (* 1d0 (abs (cl-mpm/penalty::resolve-load *penalty*))))

(declaim (notinline setup))
(defun setup (&key
                (refine 1)
                (mps 3)
                (enable-fbar t)
                (multigrid-refines 0)
                (gf 10d0)
                (kt (- 1d0 1d-6))
                (kc 0d0)
                (angle 40d0)
                (angle-r 0d0)
                (oversize-factor (- 1d0 1d-2))
                (local-length 0.01d0)
                (epsilon-scale 1d2)
                (model :MC)
                )
  (let* ((width   80d-3)
         (height 170d-3)
         (h (/ 0.01d0 refine))
         (offset width)
         (domain-height (* height 1.25))
         (density 1d3)
         (E 1.3d7)
         (nu 0.30d0)
         (domain-size (list (* 3 width) domain-height))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list width height)))
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
       :enable-aggregate nil
       :ghost-factor (* (cl-mpm/utils::calculate-p-wave-modulus E nu) 1d-6)
       :mass-update-count 1
       :enable-split nil
       :max-split-depth 8
       :split-factor (* 1d0 (sqrt 2) (/ 1d0 mps))
       :enable-fbar enable-fbar
       ;; :refinement multigrid-refines
       )))

    (setf h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
    (let* (;; (init-stress (cl-mpm/damage::mohr-coloumb-coheasion-to-tensile 40d3 angle))
           (init-stress 46d3)
           ;; (L 10d-3)
           ;; (L h)
           (L local-length)
           (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf L init-stress E)))
      (format t "Ductility ~E~%" ductility)
      (format t "Init stress ~E~%" init-stress)
      (assert (> ductility 1d0))
      (pprint block-size)
      (pprint (mapcar (lambda (e) (* (/ e h) mps)) block-size))
      (cl-mpm:add-mps
       *sim*
       (cl-mpm/setup:make-block-mps
        (list width 0d0)
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

        ;; 'cl-mpm/particle::particle-plastic-damage-frictional
        ;; 'cl-mpm/particle::particle-fpd-isotropic
        'cl-mpm/particle::particle-fpd-tcs
        :E E
        :nu nu
        :local-length L
        :ductility ductility
        :friction-angle (cl-mpm/utils:deg-to-rad angle)
        :residual-friction (cl-mpm/utils:deg-to-rad angle-r)
        :initiation-stress init-stress
        :friction-model model
        :oversize oversize-factor
        :kt-res-ratio kt
        :kc-res-ratio kc
        :psi (cl-mpm/utils:deg-to-rad angle-r)
        :plastic-damage-evolution nil
        :residual-strength (- 1d0 1d-3)

        ;; :phi (cl-mpm/utils:deg-to-rad angle)
        ;; :phi-r (cl-mpm/utils:deg-to-rad angle-r)
        ;; ;; :c-r 0d0
        ;; :softening 1d0


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
          (scale 4d0))
      ;; (cl-mpm/setup::refine-sdf
      ;;  *sim*
      ;;  (cl-mpm/setup::circle-sdf
      ;;   (list (* width 1.5d0) (* height 0.5d0) 0d0)
      ;;   4d-3)
      ;;  :refine 1
      ;;  )
      (cl-mpm/setup::apply-sdf
       *sim*
        (cl-mpm/setup::ellipse-sdf
         (list (* width 1.5d0) (* height 0.5d0) 0d0)
         (* 2d-3 scale)
         (* 2d-3 scale)
         :transform-matrix (cl-mpm/utils::rotation-matrix angle))
       ;; (cl-mpm/setup::circle-sdf
       ;;  (list (* width 1.5d0) (* height 0.5d0) 0d0)
       ;;  ;; (* scale 2d-3)
       ;;  )
       (lambda (mp sdf)
         ;; (setf (cl-mpm/particle::mp-damage mp) sdf)
         (when (<= sdf 0d0)
           ;; (pprint sdf)
           (cl-mpm/damage::set-mp-damage mp 0.99d0)
           ;; (cl-mpm/damage::set-mp-damage mp (/ (abs sdf)))
           )
         )
       )
      ;; (cl-mpm/setup::remove-sdf
      ;;  *sim*
      ;;  (cl-mpm/setup::ellipse-sdf
      ;;   (list (* width 1.5d0) (* height 0.5d0) 0d0)
      ;;   (* 4d-3 scale)
      ;;   (* 1d-3 scale)
      ;;   :transform-matrix (cl-mpm/utils::rotation-matrix angle))
      ;;  :refine 1)
      ;; (cl-mpm/setup::apply-sdf
      ;;  *sim*
      ;;  (cl-mpm/setup::ellipse-sdf
      ;;   (list (* width 1.5d0) (* height 0.5d0) 0d0)
      ;;   (* 2d-3 scale)
      ;;   (* 2d-3 scale)
      ;;   :transform-matrix (cl-mpm/utils::rotation-matrix angle))
      ;;  (lambda (mp sdf)
      ;;    (when (<= sdf 0d0)
      ;;      ;; (pprint sdf)
      ;;      (cl-mpm/damage::set-mp-damage mp 0.99d0))))
      )


    (setf (cl-mpm::sim-gravity *sim*) 0d0)
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-9)
    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 nil nil)
     :right '(0 nil nil)
     :bottom '(nil 0 nil))
    (cl-mpm::add-bcs
     *sim*
     (cl-mpm/bc::make-bc-fixed
      (list (round offset h) 0 0)
      '(0 0 nil)))

    (let* ((friction 0d0)
           (normal (cl-mpm/utils:vector-from-list '(0d0 -1d0 0d0)))
           (epsilon (* (cl-mpm/particle::calculate-p-wave-modulus E nu) epsilon-scale)))
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
         epsilon))
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
              (enable-damage t)
              (csv-dir nil)
              (lstps 5)
              (total-disp -5d-3)
              (csv-filename (format nil "load-disp.csv")))
  (unless csv-dir
    (setf csv-dir output-dir))
  (let* ((current-disp 0d0)
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
      :enable-damage enable-damage
      :damping (sqrt 2d0)
      :min-adaptive-steps 0
      :max-adaptive-steps 2
      :adaption-constant 4
      :max-damage-inc 0.90d0
      :min-damage-inc 0.1d0
      :substeps (round (* refine 50))
      ;; :sub-conv-steps 50
      :sub-conv-steps 1000
      :criteria 1d-3
      :true-stagger nil
      :save-vtk-dr t
      :save-vtk-loadstep t
      :dt-scale 1d0))))

(defun test ()
  (dolist (refine (list 1))
    (setup :mps 3
           :refine refine
           :enable-fbar t
           :angle 16d0
           :angle-r 10d0
           :gf 40d0
           :kt 1d0
           :model :DP
           :local-length 2d-2
           :oversize-factor (- 1d0 1d-2))

    (cl-mpm::iterate-over-mps
     (cl-mpm::sim-mps *sim*)
     (lambda (mp)
       (cl-mpm/fastmaths:fast-.+ (cl-mpm/particle::mp-position mp)
                                 (cl-mpm/utils::vector-from-list (list 1d-6 0d0 0d0))
                                 (cl-mpm/particle::mp-position mp)
                                 ))
     )
    ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
    ;; (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
    (run :output-dir (format nil "./output-ghost-~D/" refine)
         :enable-plastic nil
         :enable-damage t
         :refine refine)))


;; (defun test-plastic ()
;;   (dolist (refine (list 2))
;;     (setup :mps 3
;;            :refine refine
;;            :enable-fbar t
;;            :angle 16d0
;;            :angle-r 10d0
;;            :gf 40d0
;;            :kt 1d0
;;            :local-length 1d-2
;;            :oversize-factor 0d0
;;            )
;;     (cl-mpm::iterate-over-mps
;;      (cl-mpm:sim-mps *sim*)
;;      (lambda (mp)
;;        (change-class mp 'cl-mpm/particle::particle-mc)
;;        (when (> (cl-mpm/particle::mp-damage mp) 0d0)
;;          (setf (cl-mpm/particle::mp-c mp) (* (cl-mpm/particle::mp-c mp) 0.99d0)))))
;;     ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
;;     ;; (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
;;     (run :output-dir (format nil "./output-agg-~D/" refine)
;;          :enable-damage nil
;;          :enable-plastic t
;;          :refine refine)))

(defun test-oversize ()
  (let ((refine 0.5)
        (mps 4))
    (let ((gf 1d-2)
          (E 1.3d7)
          (nu 0.30d0))
      (dolist (oversize (list
                         0d0
                         ;; (- 1d0 1d-1)
                         ;; (- 1d0 1d-2)
                         ;; (- 1d0 1d-3)
                         ))
        (setup :mps mps
               :refine refine
               :enable-fbar t
               :gf 100d0
               :angle 30d0
               :angle-r 00d0
               :oversize-factor oversize)
        (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) nil
              (cl-mpm::sim-ghost-factor *sim*) (* gf (cl-mpm/utils::calculate-p-wave-modulus E nu)))
        ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
        ;; (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
        (run :output-dir (format nil "./output-~F/" oversize)
             :refine refine
             :enable-damage nil
             :enable-plastic t
             :lstps 20
             :total-disp -20d-3
             ;; :tensile t
             )))))


;;   (let ((refine 0.5)
;;         (mps 3))
;;     (dolist (gf (list 1d0 1d-2 1d-4 1d-6))
;;       (let ((E 1.3d7)
;;             (nu 0.30d0))
;;         (dolist (oversize (list 0d0))
;;           (setup :mps mps
;;                  :refine refine
;;                  :enable-fbar t
;;                  :gf 100d0
;;                  :angle 30d0
;;                  :angle-r 00d0
;;                  :oversize-factor oversize)
;;           ;; (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) nil
;;           ;;       (cl-mpm::sim-ghost-factor *sim*) (* gf (cl-mpm/utils::calculate-p-wave-modulus E nu)))
;;           (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
;;                 (cl-mpm::sim-ghost-factor *sim*) nil)
;;           ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
;;           ;; (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
;;           (run :output-dir (format nil "./output-~F/" gf)
;;                :refine refine
;;                :enable-damage nil
;;                :enable-plastic nil
;;                :lstps 3
;;                :total-disp -2d-3
;;                ;; :tensile t
;;                ))))))

(defun test-tension-compression ()
  (let ((refine 2)
        (mps 4))
    (dolist (tension (list
                      nil
                      ;; t
                      ))
      (setup :mps mps
             :refine refine
             :enable-fbar t
             :multigrid-refines 0
             :gf 60d0
             :angle 15d0
             :angle-r 10d0
             :kt (- 1d0 1d-9)
             :kc 0d0
             )
      (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
      ;; (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
      (run :output-dir (format nil "./output-~D-~A/" refine (if tension "tension" "compression"))
           :refine refine
           :tensile tension
           :enable-plastic nil
           :total-disp -10d-3
           )))
  )

(defun test-degredation ()
  (let ((refine 2)
        (mps 3))
    (dolist (particle (list
                       ;; 'cl-mpm/particle::particle-fpd-isotropic
                       ;; 'cl-mpm/particle::particle-fpd-tcs
                       ;; 'cl-mpm/particle::particle-fpd-isotropic
                       'cl-mpm/particle::particle-fpd-tcs-strain
                       ;; 'cl-mpm/particle::particle-fpd-spectral
                       ;; 'cl-mpm/particle::particle-fpd-spectral-strain
                       ))
      (setup :mps mps
             :refine refine
             :enable-fbar t
             :multigrid-refines 0
             :gf 10d0
             :angle 15d0
             :angle-r 10d0
             :oversize-factor (- 1d0 1d-3)
             ;; :kt (- 1d0 1d-3)
             )
      (cl-mpm::iterate-over-mps
       (cl-mpm:sim-mps *sim*)
       (lambda (mp)
         (change-class mp particle)))
      (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
      ;; (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
      (run :output-dir (format nil "./output-~D-~A/" refine particle)
           :refine refine
           :enable-plastic t))))


(defun test-model ()
  (let ((refine 1)
        (mps 3))
    (dolist (model (list
                    :PRINC-STRAIN
                    :RANKINE
                    :RANKINE-SMOOTH
                    :DP
                    :MC
                    :SE
                    ))
      (setup :mps mps
             :refine refine
             :gf 40d0
             :angle 30d0
             :angle-r 00d0
             :model model)
      ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
      ;; (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
      (run :output-dir (format nil "./output-~A-~D/" model refine)
           :refine refine
           :enable-plastic nil
           :enable-damage t
           :tensile nil))))

(defun test-agg ()
  (let ((refine 0.5))
    (dolist (agg (list nil))
      (setup :mps 3
             :refine refine
             :enable-fbar nil
             :angle 30d0
             :angle-r 0d0
             :gf 10d0
             :kt 1d0;(- 1d0 1d-6)
             )
      (when agg
        (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t
              (cl-mpm::sim-ghost-factor *sim*) nil))
      ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
      ;; (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
      (run :output-dir (format nil "./output-~D-~A/" refine agg)
           ;; :enable-plastic t
           :refine refine))))

(defun test-models ()
  (dolist (refine (list 1))
    (dolist (type (list :NORMAL :EKL :LL))
      (setup :mps 3
             :refine refine
             :enable-fbar t
             :angle 16d0
             :angle-r 0d0
             :gf 40d0
             :kt 1d0
             :local-length 1d-2
             )
      (case type
        (:NORMAL )
        (:EKL (setf (cl-mpm/damage::sim-enable-ekl *sim*) t))
        (:LL (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t))
        )
      ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
      ;; (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
      (run :output-dir (format nil "./output-~A-~D/" type refine)
           :enable-plastic nil
           :refine refine))))

(defun test-angle ()
   (dolist (refine (list 1))
     (dolist (model (list :MC
                          :DP
                          ))
       (dolist (angle (list
                       ;; 5d0
                       15d0
                       30d0
                       45d0
                       60d0
                       ))
         (setup :mps 3
                :refine refine
                :enable-fbar t
                :angle angle
                :angle-r 1d0
                :gf 40d0
                :model model
                :epsilon-scale 1d2
                )
         ;; (cl-mpm::iterate-over-mps
         ;;  (cl-mpm::sim-mps *sim*)
         ;;  (lambda (mp)
         ;;    (cl-mpm/fastmaths:fast-.+ (cl-mpm/particle::mp-position mp)
         ;;                              (cl-mpm/utils::vector-from-list (list 1d-6 0d0 0d0))
         ;;                              (cl-mpm/particle::mp-position mp)
         ;;                              ))
         ;;  )
         ;; (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
         ;; (setf (cl-mpm/damage::sim-enable-ekl *sim*) t)
         (let ((output-dir (format nil "./output-~A-~F-~D/" model angle refine)))
           (format t "Testing ~A~%" output-dir)
           (time
            (run :output-dir output-dir
                 :lstps 10
                 :total-disp -5d-3
                 :enable-damage t
                 :refine refine)))))))


;; (let ((h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
;;       (nd (cl-mpm/mesh::mesh-nd (cl-mpm:sim-mesh *sim*)))
;;       )
;;   (pprint (/ (cl-mpm::reduce-over-mps
;;               (cl-mpm:sim-mps *sim*)
;;               (lambda (mp)
;;                 (multiple-value-bind (l v) (cl-mpm/utils:eig (cl-mpm/particle::mp-true-domain mp))
;;                   (let* ((abs-l (mapcar #'abs l))
;;                          (max-l (reduce #'max abs-l)))
;;                     max-l)))
;;               #'max)
;;              h)))





(pprint
 (cl-mpm/implicit::tensor-2nd-partial-deriv
  (cl-mpm/utils:matrix-from-list
   (list 1d0 0d0 0d0
         0d0 2d0 0d0
         0d0 0d0 3d0)
   ) #'log (lambda (x) (/ 1d0 x))))
