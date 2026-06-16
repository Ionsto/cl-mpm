(defpackage :cl-mpm/examples/penalty/billet
  (:use :cl
   :cl-mpm/example
        :cl-mpm/utils)
  )
(in-package :cl-mpm/examples/penalty/billet)
(declaim (notinline plot-domain))
(defun plot-domain ()
  (when *sim*
    (cl-mpm/plotter:simple-plot
     *sim*
     :plot :deformed
     :trial t
     :colour-func (lambda (mp) (cl-mpm/particle::mp-index mp)))))

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-elastic) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  ;; (cl-mpm::update-domain-polar-2d mesh mp dt)
  ;; (cl-mpm::update-domain-polar mesh mp dt)
  ;; (cl-mpm::update-domain-max-corner-2d mesh mp dt)
  ;; (cl-mpm::update-domain-midpoint mesh mp dt)
  (cl-mpm::update-domain-stretch mesh mp dt)
  ;; (cl-mpm::scale-domain-size mesh mp)
  )

(defun setup (&key (refine 1) (mps 3)
                (friction 1d0)
                )
  (let* ((L 10d0)
         (d 1d0)
         (domain-width 30d0)
         (h (/ 1d0 refine))
         (height 10d0)
         (domain-height (+ height (* 2 h)))
         (density 1000d3)
         (E 1d6)
         (domain-size (list domain-width domain-height))
         (element-count (mapcar (lambda (x) (round x h)) domain-size))
         (block-size (list height height)))

    (setf *sim* (cl-mpm/setup::make-simple-sim
                 h
                 element-count
                 :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
                 :args-list (list :enable-aggregate t
                                  :enable-fbar t
                                  :enable-split nil
                                  :max-split-depth 6
                                  :mp-removal-size nil
                                  :split-factor (* 1.1d0 (sqrt 1d0) (/ 1d0 mps))
                                  )))
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list 0 0)
      block-size
      (mapcar (lambda (e m) (* (/ e h) m)) block-size
              (list (* 1 mps)
                    mps))
      density
      'cl-mpm/particle::particle-vm
      :E E
      :nu 0.3d0
      :rho 20d3
      ;; 'cl-mpm/particle::particle-elastic
      ;; :E E
      ;; :nu 0.2d0
      ;; :index 0
      ;; :gravity-axis (cl-mpm/utils:vector-zeros)
      ))
    (setf (cl-mpm::sim-gravity *sim*) 0d0)
    (cl-mpm/setup::set-mass-filter *sim* 1d3 :proportion 1d-15)
    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(0 nil nil)
     :bottom '(nil 0 nil)
     )

    (let ((epsilon-scale 1d2)
          )
      (defparameter *penalty*
        (cl-mpm/penalty::make-bc-penalty-distance-point
         *sim*
         (cl-mpm/utils:vector-from-list '(0d0 -1d0 0d0))
         (cl-mpm/utils:vector-from-list (list (* 0.5d0 domain-width)
                                              height
                                              0d0))
         (/ domain-width 2)
         (* E epsilon-scale)
         friction
         0d0)))
    (defparameter *current-inc* 0d0)
    (cl-mpm:add-bcs-force-list
     *sim*
     *penalty*)
    (setf (cl-mpm:sim-dt *sim*)
          (* 0.5d0 (cl-mpm/setup:estimate-elastic-dt *sim*)))
    (setf *run-sim* t))
  (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
  )

(defparameter *penalty-controller* nil)
(defun run-time ()
  (change-class *sim* 'cl-mpm/aggregate::mpm-sim-agg-usf)
  ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic)
  (setf (cl-mpm::sim-velocity-algorithm *sim*) :BLEND)
  (let ((step 0))
    (cl-mpm/dynamic-relaxation::run-time
     *sim*
     :output-dir (format nil "./output/")
     :plotter
     (lambda (sim)
       (plot-domain)
       (vgplot:title (format nil "Step ~D" step))
       (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step)) :terminal "png size 1920,1080")
       (incf step))
     :damping 1d-3
     :dt 10d0
     :total-time 1000d0
     ;; :dt-scale 1d0
     :dt-scale (/ 0.5d0 (sqrt 1d2))
     :post-conv-step
     (lambda (sim)
       (let ((disp-rate -0.01d0))
         (defparameter *penalty-controller*
           (cl-mpm/bc::make-bc-closure
            nil
            (lambda ()
              (cl-mpm/penalty::bc-increment-center
               *penalty*
               (cl-mpm/utils:vector-from-list (list 0d0 (* disp-rate (cl-mpm:sim-dt *sim*)) 0d0))))))
         (cl-mpm::add-bcs-force-list
          *sim*
          *penalty-controller*)
         )))))
(defun run (&key (output-dir (format nil "./output/")))
  (let* ((lstps 10)
         (total-disp -5d0)
         (step 0))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (defparameter *disp* 0d0)
    (cl-mpm/dynamic-relaxation::run-load-control
     *sim*
     :pre-step
     (lambda ()
       (with-open-file (stream (merge-pathnames "./load.csv" output-dir) :direction :output :if-exists :supersede)
         (format stream "step,load~%")))
     :output-dir output-dir
     :plotter (lambda (sim) (plot-domain))
     :loading-function
     (lambda (i)
       (setf *disp* (* i total-disp))
       (cl-mpm/penalty::bc-set-displacement
        *penalty*
        (cl-mpm/utils:vector-from-list (list 0d0 (* i total-disp) 0d0))))
     :post-conv-step
     (lambda (sim)
       (plot-domain)
       (vgplot:title (format nil "Step ~D" step))
       (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" step)) :terminal "png size 1920,1080")
       (with-open-file (stream (merge-pathnames output-dir "./load.csv") :direction :output :if-exists :append)
         (format stream "~F,~F~%" *disp* (cl-mpm/penalty::bc-penalty-load *penalty*)))
       (incf step))
     :load-steps lstps
     :enable-plastic t
     :damping (sqrt 2d0)
     :substeps 100
     :conv-steps 100
     :criteria 1d-9
     :save-vtk-dr nil
     :save-vtk-loadstep t
     :dt-scale 0.9d0)))

(defun test ()
  (cl-mpm/utils:set-workers 8)
  (setup :mps 6 :refine 1)
  (run)
  ;; (dolist (r (list 1 2 3))
  ;;   (dolist (mps (list 2 4))
  ;;     (setup :mps mps :refine r)
  ;;     (run :output-dir (format nil "./output-~D-~D/" r mps))))
  ;; (run-time)
  )

(defparameter *split-cartesian* nil)
(defun cl-mpm::split-mps (sim)
  "Split mps that match the split-criteria"
  (dotimes (d (cl-mpm/mesh::mesh-nd (cl-mpm::sim-mesh sim)))
    (if *split-cartesian*
        (cl-mpm::split-mps-cartesian sim)
        (cl-mpm::split-mps-eigenvalue sim))))

(defun test-mps ()
  (cl-mpm/utils:set-workers 8)
  (dolist (mps (list 2 3 4 5 6))
    (setup :mps mps :refine 1)
    (setf (cl-mpm::sim-allow-mp-split *sim*) nil)
    (setf (cl-mpm::sim-mp-removal-size *sim*) nil)
    (ignore-errors
     (run :output-dir (format nil "./output-~D/" mps)))))

(defun test-splitting ()
  (cl-mpm/utils:set-workers 8)
  (dolist (split (list nil t))
    (dolist (mps (list 3))
      (setup :mps mps :refine 1)
      (setf (cl-mpm::sim-allow-mp-split *sim*) t)
      (setf *split-cartesian* split)
      (setf (cl-mpm::sim-mp-removal-size *sim*) nil)
      (run :output-dir (format nil "./output-split-~A/" (if *split-cartesian* "cartesian" "eigenvalue"))))))

(defun test-friction ()
  (let ((mps 6))
    (cl-mpm/utils:set-workers 8)
    (setup :mps mps :refine 1 :friction 1d0)
    (run :output-dir "./output-noslip/")
    (setup :mps mps :refine 1 :friction 0d0)
    (run :output-dir "./output-slip/"))
  )


;; (pprint (cl-mpm/fastmaths:mag
;;          (magicl:@ (cl-mpm/particle::mp-true-domain *mp*)
;;                    (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0)))))


;; (let* ((angle-r (cl-mpm/utils::deg-to-rad 170d0))
;;        (S (cl-mpm::make-scaling-matrix (cl-mpm/utils::vector-from-list (list 1d0 1d0 0d0)) 0.5d0))
;;        (M (cl-mpm/utils::matrix-from-list (list 1d0 0.5d0 0d0
;;                                                 0d0 1d0 0d0
;;                                                 0d0 0d0 1d0))))
;;   (pprint (cl-mpm/utils::diagonal (magicl:@ S (cl-mpm/utils::matrix-from-diag-vec (cl-mpm/utils:vector-from-list (list 1d0 1d0 0d0))))))
;;   (let ((vec (cl-mpm/utils::vector-from-list (list 1d0 0d0 0d0))))
;;     (pprint (cl-mpm/fastmaths::dot vec (magicl:@ S M vec))))
;;   (let ((vec (cl-mpm/utils::vector-from-list (list 0d0 1d0 0d0))))
;;     (pprint (cl-mpm/fastmaths::dot vec (magicl:@ S M vec))))
;;   ;; (multiple-value-bind (l v) (magicl::eig M)
;;   ;;   (setf v (magicl:.realpart v))
;;   ;;   (pprint l)
;;   ;;   (pprint v)
;;   ;;   ;; (pprint (magicl:@ S M))
;;   ;;   (let* ((S (cl-mpm::make-scaling-matrix (cl-mpm/utils::matrix-column v 0) 0.5d0))
;;   ;;          (M (magicl:@ S M))
;;   ;;          ;; (r 0d0)
;;   ;;          (x (list))
;;   ;;          (y (list)))
;;   ;;     (loop for r from 0d0 to (* 2 pi) by 0.01d0
;;   ;;           do (let* ((vec (cl-mpm/utils::vector-from-list (list (cos r) (sin r) 0d0)))
;;   ;;                     (distance (cl-mpm/fastmaths::dot vec (magicl:@ M vec)))
;;   ;;                     (result (cl-mpm/fastmaths::fast-scale-vector vec (/ 1d0 (sqrt distance))))
;;   ;;                     )
;;   ;;                (push (cl-mpm/utils::varef result 0) x)
;;   ;;                (push (cl-mpm/utils::varef result 1) y)))
;;   ;;     (vgplot:plot x y ";;with points")
;;   ;;     (vgplot:format-plot t "set xrange [~f:~f]" -2d0 2d0)
;;   ;;     (vgplot:format-plot t "set yrange [~f:~f]" -2d0 2d0)
;;   ;;     ))
;;   )

;; (let ((dF (cl-mpm/utils::matrix-from-list (list 1d0 0d0 0d0
;;                                                 0d0 2d0 0d0
;;                                                 0d0 0d0 1d0
;;                                                 ))))
;;   (multiple-value-bind (u s vt) (magicl:svd dF)
;;     (let* (;; (vt (magicl:transpose vt))
;;            (R (magicl:@ u vt))
;;            (U (magicl:@ (magicl:transpose vt) s vt))
;;            )
;;       ;; (pprint R)
;;       ;; (pprint U)
;;       (let ((test (cl-mpm/utils::matrix-from-list (list 1d0 0d0 0d0
;;                                                         0d0 1d0 0d0
;;                                                         0d0 0d0 1d0))))
;;         (pprint (magicl:@ dF test (magicl:transpose dF)))
;;         (pprint (magicl:@ R U test (magicl:transpose R)))
;;         ;; (pprint (magicl:@ R U test))
;;         )
;;       ;; (setf true-domain (magicl:@ R (magicl:@ true-domain U) (magicl:transpose R)))
;;       ;; (setf true-domain (magicl:@ R true-domain (magicl:transpose R) U))
;;       )))
