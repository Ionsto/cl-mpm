(defpackage :cl-mpm/examples/penalty/sliding
  (:use :cl
   :cl-mpm/example
        :cl-mpm/utils))
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(in-package :cl-mpm/examples/penalty/sliding)
(declaim (optimize (debug 3) (safety 3) (speed 0)))


(declaim (notinline setup))
(defun setup (&key (refine 1) (mps 3)
                (epsilon-scale 1d1)
                (mu 0.25d0))
  (let* ((d 8d0)
         ;; (h (/ d 2))
         (h 1d0)
         (dsize 50d0)
         (density 2d3)
         (mesh-resolution (/ 0.125d0 refine))
         (offset (* 2d0 mesh-resolution))
         (domain-size (list dsize (+ offset (* 2 d))))
         (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
         (block-size (list d h))
         (E 1d9)
         (nu 0.3d0)
         )
    (setf *sim* (cl-mpm/setup::make-simple-sim
                 mesh-resolution
                 element-count
                 :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static
                 :args-list
                 (list
                  :enable-aggregate t
                  :enable-fbar t)))
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup:make-block-mps
      (list offset offset)
      block-size
      (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
      density
      'cl-mpm/particle::particle-elastic
      :E E
      :nu nu
      :index 0))
    (setf (cl-mpm::sim-gravity *sim*) -10d0)
    (defparameter *normal-force* (* 10d0 density (* d h)))
    (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-15)

    (cl-mpm/setup::setup-bcs
     *sim*
     :left '(nil 0 nil))

    (let* ((epsilon (* (cl-mpm/particle::calculate-p-wave-modulus E nu) epsilon-scale))
           ;; (epsilon (* E 1d2))
           )
      ;; (defparameter *loader*
      ;;   (cl-mpm/penalty::make-bc-penalty-displacment
      ;;    *sim*
      ;;    (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0))
      ;;    (* E epsilon-scale)))
      (defparameter *loader*
        (cl-mpm/penalty::make-bc-penalty-distance-point
         *sim*
         (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0))
         (cl-mpm/utils:vector-from-list (list offset (+ offset (* d 0.5d0)) 0d0))
         (* d 0.5d0)
         epsilon
         0d0
         0d0
         ))
      (defparameter *follower*
        (cl-mpm/penalty::make-bc-penalty-distance-point
         *sim*
         (cl-mpm/utils:vector-from-list (list -1d0 0d0 0d0))
         (cl-mpm/utils:vector-from-list (list (+ d offset) (+ offset (* d 0.5d0)) 0d0))
         (* d 0.5d0)
         epsilon
         0d0
         0d0
         ))
      (defparameter *floor*
        (cl-mpm/penalty::make-bc-penalty-distance-point 
         *sim*
         (cl-mpm/utils:vector-from-list (list 0d0 1d0 0d0))
         (cl-mpm/utils:vector-from-list (list (* dsize 0.5 ) offset 0d0))
         (* dsize 0.5d0)
         epsilon
         mu
         0d0
         )))
    (cl-mpm::add-bcs-force-list
     *sim*
     *floor*)
    (cl-mpm::add-bcs-force-list
     *sim*
     *loader*)
    (cl-mpm::add-bcs-force-list
     *sim*
     *follower*)
    (setf *run-sim* t))
  (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))

  )

(defun plot-vel ()
  (let ((times (reverse *data-time*))
        (vel (reverse *data-vel*))
        )
  (when (> (length *data-time*) 1)
    (let ((dt (- (second times) (first times))))
      (let ((acc (mapcar (lambda (x) (/ x dt)) (mapcar #'- (rest vel) (butlast vel)))))
        (vgplot:plot
         ;; times
         ;; vel
         (rest times)
         acc
         "Velocity"))))))

(defun run-time ()
  ;; (change-class *sim* 'cl-mpm/aggregate::mpm-sim-agg-usf)
  (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-implict-dynamic)
  (setf (cl-mpm::sim-velocity-algorithm *sim*) :FLIP)
  (defparameter *data-time* (list 0d0))
  (defparameter *data-vel* (list 0d0))
  (let ((grav-angle (cl-mpm/particle::mp-gravity-axis (aref (cl-mpm::sim-mps *sim*) 0))))
    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (setf (cl-mpm/particle::mp-gravity-axis mp) (cl-mpm/utils:vector-from-list (list 0d0 1d0 0d0)))))
    (cl-mpm/dynamic-relaxation::run-time
     *sim*
     :output-dir (format nil "./output/")
     :plotter (lambda (sim) (plot-vel))
     :damping 0d-4
     :dt 0.01d0
     :total-time 10d0
     :dt-scale 1d0
     :post-iter-step
     (lambda (sim)
       (push (cl-mpm::sim-time *sim*) *data-time*)
       (push (cl-mpm::reduce-over-mps
              (cl-mpm:sim-mps *sim*)
              (lambda (mp)
                (*
                 (cl-mpm/particle::mp-volume mp)
                 (cl-mpm/utils::varef (cl-mpm/particle::mp-velocity mp) 0)))
              #'+)
             *data-vel*))
     :post-conv-step
     (lambda (sim)
       (cl-mpm::iterate-over-mps
        (cl-mpm:sim-mps *sim*)
        (lambda (mp)
          (setf (cl-mpm/particle::mp-gravity-axis mp) grav-angle)))))))



(declaim (notinline get-load))
(defun get-load ()
  (/
   (-
    (cl-mpm/penalty::resolve-load *loader*)
    (cl-mpm/penalty::resolve-load *follower*))
   *normal-force*))

(defun output-disp-header (output-dir)
  (with-open-file (stream (merge-pathnames "disp.csv" output-dir) :direction :output :if-exists :supersede)
    (format stream "disp,load~%")))

(defparameter *displacement* 0d0)
(defun output-disp-data (output-dir)
  (with-open-file (stream (merge-pathnames "disp.csv" output-dir) :direction :output :if-exists :append)
    (format stream "~f,~f~%"
            *displacement*
            (get-load))))
(defun plot ()
  (vgplot:plot
   *data-disp*
   *data-load*
   "Load")
  )

(declaim (notinline run))
(defun run (&key (output-dir "./output/") (results-dir nil))
  (unless results-dir
    (setf results-dir output-dir))
  (let ((total-disp 1d-1)
        (current-disp 0d0))
    (defparameter *data-disp* (list 0d0))
    (defparameter *data-load* (list 0d0))
    (defparameter *displacement* 0d0)
    (vgplot:close-all-plots)
    (cl-mpm/dynamic-relaxation::run-load-control
     *sim*
     :output-dir output-dir
     :plotter
     (lambda (sim)
       (plot))
     :load-steps 50
     :damping 1d0
     :substeps 10
     :criteria 1d-3
     :save-vtk-dr nil
     :save-vtk-loadstep t
     :loading-function
     (lambda (percent)
       (setf current-disp (* total-disp percent))
       (setf *displacement* current-disp)
       (cl-mpm/penalty::bc-set-displacement
        *loader*
        (cl-mpm/utils:vector-from-list
         (list
          current-disp
          0d0
          0d0)))
       (cl-mpm/penalty::bc-set-displacement
        *follower*
        (cl-mpm/utils:vector-from-list
         (list
          current-disp
          0d0
          0d0))))
     :pre-step
     (lambda ()
       (output-disp-header output-dir)
       (output-disp-data output-dir))
     :post-conv-step
     (lambda (sim)
       (push current-disp *data-disp*)
       (push (get-load) *data-load*)
       (output-disp-data output-dir))
     :dt-scale 1d0))
  )

(defun test ()
  (let* ((mu 0.5d0))
    (setup
     :refine 1
     :mps 2
     :mu mu)
    (format t "Maximum normalised force ~E~%" mu)
    (run)))

(defun angle-calc (angle mu)
  (let* ((psi (* (- 90d0 (abs angle)) (/ pi 180)))
         (gravity 10d0))
    (format t "~%Estimated acceleration ~F~%" (* gravity (max 0d0 (- (cos psi) (* mu (sin psi)))) ))))


(defun test-eps ()
  (let* ((mu 0.5d0))
    (dolist (eps '(1d-2 1d-1 1d0))
      (setup :refine 1
             :mps 2
             :mu mu
             :epsilon-scale 1d3)
      (setf (cl-mpm::sim-gravity *sim*) -10d0)
      (setf cl-mpm/penalty::*friction-epsilon-scale* eps)
      (run :output-dir (format nil "./output-eps-~E/" eps)))))


(defun test-refine ()
  (let* ((mu 0.5d0))
    (dolist (eps '(2d0 4d0))
      (setup :refine eps
             :mps 2
             :mu mu
             :epsilon-scale 1d2)
      (setf cl-mpm/penalty::*friction-epsilon-scale* 1d0)
      (run :output-dir (format nil "./output-refine-~E/" eps)))))


;; (let ((mesh (cl-mpm:sim-mesh *sim*)))
;;   (cl-mpm/aggregate::iterate-over-cell-shape-local
;;    mesh
;;    (cl-mpm/mesh::get-cell mesh (list 0 0 0))
;;    (cl-mpm/utils:vector-from-list (list 0.1d0 0.1d0 0.1d0))
;;    (lambda (n w g)
;;      (pprint w))))
