(defpackage :cl-mpm/examples/beam
  (:use :cl
   :cl-mpm/example
   :cl-mpm/utils
   ))
(in-package :cl-mpm/examples/beam)
(declaim (notinline plot-domain))

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-elastic) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  (cl-mpm::update-domain-deformation mesh mp dt)
  ;; (cl-mpm::update-domain-stretch mesh mp dt)
  ;; (cl-mpm::update-domain-polar mesh mp dt)
  )

(defun plot-domain ()
  (when *sim*
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :trial t
     :colour-func (lambda (mp) (cl-mpm/particle::mp-index mp)))))

(defun setup (&key (refine 1) (mps 3))
  (let* ((L 10d0)
         (d 1d0)
         (dsize 11d0)
         (density 7750d3)
         (mesh-resolution (/ 0.25d0 refine))
         (domain-size (list dsize dsize))
         (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
         (depth 1d0)
         (offset-y (- 11 (* 2 d)))
         (center-y (+ offset-y (/ d 2)))
         (block-size (list L d)))
    (setf *sim* (cl-mpm/setup::make-simple-sim
                 mesh-resolution
                 element-count
                 ;; :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-multi grid
                 :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
                 ;; :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static
                 ;; :sim-type 'cl-mpm/implicit::mpm-sim-implicit
                 :args-list (list
                             :mass-update-count 1
                             :mp-removal-size nil
                             ;; :refinement 2
                             ;; :ghost-factor (* 12d6 1d-4)
                             :enable-aggregate t
                             :ghost-factor nil
                             :enable-split nil
                             )))
    (setf mesh-resolution (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
    (format t "Mesh resolution ~E~%" mesh-resolution)
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0 offset-y)
      block-size
      (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
      density
      'cl-mpm/particle::particle-elastic
      :E 12d6
      :nu 0.2d0
      :index 0
      :gravity-axis (cl-mpm/utils:vector-zeros)
      )
     )
    (cl-mpm::domain-sort-mps *sim*)
    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (setf (cl-mpm/particle::mp-index mp) 0)))
    (let ((mass 0d0))
      (cl-mpm:iterate-over-mps
       (cl-mpm:sim-mps *sim*)
       (lambda (mp)
         (when
             (and
              (<= (abs (- (+ (varef (cl-mpm/particle:mp-position mp) 0)
                             (* 0.5d0 (varef (cl-mpm/particle::mp-domain-size mp) 0)))
                          L)) 1d-3)
              (< (abs (-  (varef (cl-mpm/particle:mp-position mp) 1) center-y)) (+ 1d-3 (* 1d0 (varef (cl-mpm/particle::mp-domain-size mp) 1)))))
           (setf mass (cl-mpm/particle::mp-mass mp))
           (setf
            ;; (cl-mpm/particle::mp-gravity mp) (/ -100d3 (* 2d0 (cl-mpm/particle::mp-mass mp)))
            (cl-mpm/particle::mp-index mp) 1)
           (setf (cl-mpm/particle::mp-gravity-axis mp) (cl-mpm/utils:vector-from-list (list 0d0 1d0 0d0)))
           )))
      (setf (cl-mpm:sim-gravity *sim*) (/ -100d3 (* 2d0 mass))))

    (cl-mpm/setup::set-mass-filter *sim* density :proportion 0d-9)
    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 nil nil))
    (cl-mpm::add-bcs
     *sim*
     (cl-mpm/bc::make-bc-fixed
      (list 0 (round center-y mesh-resolution) 0)
      (list 0 0 nil)))
    (setf *run-sim* t))
  (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
  )


(defun get-pos ()
  (let ((pos (cl-mpm/utils:vector-zeros)))
    (cl-mpm:iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (when (= 1 (cl-mpm/particle::mp-index mp))
         (setf pos (cl-mpm/fastmaths:fast-.+
                    pos
                    (cl-mpm/fastmaths:fast-.+
                     (cl-mpm/particle::mp-displacement mp)
                     (cl-mpm/particle::mp-displacement-increment mp)
                     ))))))
    (cl-mpm/fastmaths::fast-scale! pos -0.5d0)
    pos))


(defparameter *analytic-x*
  (list
   0d0
   0.083d0
   0.162d0
   0.235d0
   0.302d0
   0.494d0
   0.603d0
   0.670d0
   0.714d0
   0.744d0
   0.767d0
   0.785d0
   0.799d0
   0.811d0
   ;; 1.000d0
   ))

(defparameter *analytic-y*
  (list
   0d0
   0.004d0
   0.016d0
   0.034d0
   0.056d0
   0.160d0
   0.255d0
   0.329d0
   0.388d0
   0.434d0
   0.472d0
   0.504d0
   0.531d0
   0.555d0
   ;; 1.000d0
   ))


(defmethod cl-mpm::update-sim :after ((sim cl-mpm::mpm-sim))
  ;; (with-accessors ((mesh cl-mpm::sim-mesh)
  ;;                  (mps cl-mpm::sim-mps)
  ;;                  (dt cl-mpm::sim-dt)
  ;;                  (damping cl-mpm::sim-damping))
  ;;     sim
  ;;   (cl-mpm::g2p mesh mps 0d0 0d0 :TRIAL))
  ;; (let ((pos (get-pos)))
  ;;   (format t "Disp ~E ~E ~%" (cl-mpm/utils::varef pos 0) (cl-mpm/utils::varef pos 1))
  ;;   (push (cl-mpm/utils::varef pos 0) *data-y*)
  ;;   (push (cl-mpm/utils::varef pos 1) *data-x*))
  )

(defparameter *data-x* (list))
(defparameter *data-y* (list))
(defparameter *data-type* (list))
(defun run (&key (output-dir "./output/")
              (substeps 1)
              (load-steps 10)
              (dt-scale 0.9d0)
              )
  (defparameter *data-x* (list 0d0))
  (defparameter *data-y* (list 0d0))
  (defparameter *data-type* (list nil))
  (let ((i 1))
    (ensure-directories-exist output-dir)
    (cl-mpm/output::save-vtk (merge-pathnames  output-dir (format nil "sim_~5,'0d.vtk" 0)) *sim*)
    (time
     (cl-mpm/dynamic-relaxation::run-load-control
      *sim*
      :output-dir output-dir
      :plotter (lambda (sim)
                 ;; (plot-domain)
                 ;;(plot-analytic)
                 )
      :load-steps load-steps
      :damping (sqrt 2d0)
      :dt-scale dt-scale
      :substeps substeps
      :conv-steps 1000
      :criteria 1d-9
      ;; :sub-conv-steps 1000
      :save-vtk-dr t
      :save-vtk-loadstep t
      :post-iter-step
      (lambda (sim o e)
        ;; (with-accessors ((mesh cl-mpm::sim-mesh)
        ;;                  (mps cl-mpm::sim-mps)
        ;;                  (dt cl-mpm::sim-dt)
        ;;                  (damping cl-mpm::sim-damping))
        ;;     *sim*
        ;;   (cl-mpm::g2p mesh mps 0d0 0d0 :TRIAL))
        ;; (let ((pos (get-pos)))
        ;;   (format t "Disp ~E ~E ~%" (cl-mpm/utils::varef pos 0) (cl-mpm/utils::varef pos 1))
        ;;   (push (cl-mpm/utils::varef pos 0) *data-y*)
        ;;   (push (cl-mpm/utils::varef pos 1) *data-x*)
        ;;   (push nil *data-type*)
        ;;   ;; (cl-mpm/output::save-vtk (merge-pathnames  output-dir (format nil "sim_~5,'0d.vtk" i)) *sim*)
        ;;   (incf i))
        )
      :post-conv-step
      (lambda (sim)
        ;; (let ((pos (get-pos)))
        ;;   (format t "Disp ~E ~E ~%" (cl-mpm/utils::varef pos 0) (cl-mpm/utils::varef pos 1))
        ;;   (push (cl-mpm/utils::varef pos 0) *data-y*)
        ;;   (push (cl-mpm/utils::varef pos 1) *data-x*)
        ;;   (push t *data-type*)
        ;;   ;; (cl-mpm/output::save-vtk (merge-pathnames  output-dir (format nil "sim_~5,'0d.vtk" i)) *sim*)
        ;;   (incf i))
        ))))
  (output-disp output-dir)
  (plot-analytic))

(defun plot-analytic ()
  (vgplot:plot (mapcar (lambda (x) (* x 10d0)) *analytic-x*)
               (mapcar (lambda (x) (* x 10d0)) *analytic-y*)
               "Analytic"
               *data-x*
               *data-y*
               ;; "Newton-Raphson;;with points"
               "Newton-Raphson"
               ))

(defun output-disp (output-dir)
  (with-open-file (stream (merge-pathnames "disp.csv" output-dir) :direction :output :if-exists :supersede)
    (format stream "iter,x,y,type~%")
    (loop
      for i from 0
      for x in *data-x*
      for y in *data-y*
      for type in *data-type*
      do (format stream "~D,~f,~f,~D~%" i x y (if type 1 0)))))

(defun test-both-agg ()
  (cl-mpm/utils:set-workers 8)
  (dolist (agg (list t nil))
    (setup :mps 3 :refine 1)
    (setf (cl-mpm::sim-ghost-factor *sim*) nil
          (cl-mpm/aggregate::sim-enable-aggregate *sim*) agg)
    (cl-mpm/setup::set-mass-filter *sim* 7750d3 :proportion 1d-15)
    (run :output-dir (format nil "./output-~A/" agg) :substeps 400 :dt-scale 0.9d0 :load-steps 10)))

(defun test-agg ()
  (cl-mpm/utils:set-workers 8)
  (setup :mps 3 :refine 1)
  (setf (cl-mpm::sim-ghost-factor *sim*) nil
        (cl-mpm/aggregate::sim-enable-aggregate *sim*) t)
  (cl-mpm/setup::set-mass-filter *sim* 7750d3 :proportion 1d-9)
  (run :output-dir "./output-agg/" :substeps 400 :dt-scale 0.5d0 :load-steps 10))

(defun test-noagg ()
  (cl-mpm/utils:set-workers 8)
  (setup :mps 3 :refine 1)
  (setf (cl-mpm::sim-ghost-factor *sim*) nil
        (cl-mpm/aggregate::sim-enable-aggregate *sim*) nil)
  (cl-mpm/setup::set-mass-filter *sim* 7750d3 :proportion 1d-9)
  (run :output-dir "./output-noagg/" :substeps 400 :dt-scale 0.5d0 :load-steps 10))

;; (defun test-cuts ()
;;   (cl-mpm/utils:set-workers 16)
;;   (setup :mps 3 :refine 1)
;;   (setf (cl-mpm::sim-ghost-factor *sim*) nil
;;         (cl-mpm/aggregate::sim-enable-aggregate *sim*) t)
;;   (run :output-dir "./output-agg/" :substeps 100))


;; (defun test-cuts-noagg ()
;;   (cl-mpm/utils:set-workers 16)
;;   (setup :mps 3 :refine 1)
;;   (setf (cl-mpm::sim-ghost-factor *sim*) nil
;;         (cl-mpm/aggregate::sim-enable-aggregate *sim*) nil)
;;   (run :output-dir "./output-noagg/" :substeps 100))

(defun test ()
  (cl-mpm/utils:set-workers 16)
  (setup :mps 3 :refine 1)
  (change-class *sim* 'cl-mpm/implicit::mpm-sim-implicit)
  (run :output-dir "./output-nr/" :substeps 1)
  ;; (setup :mps 3 :refine 1)
  ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static)
  ;; (run :output-dir "./output-dr/" :substeps 400)
  )




;; (cl-mpm::iterate-over-mps-serial
;;  (cl-mpm:sim-mps *sim*)
;;  (lambda (mp)
;;    (pprint (cl-mpm/dynamic-relaxation::estimate-ul-enhancement mp 2))))
