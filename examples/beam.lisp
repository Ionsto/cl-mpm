(defpackage :cl-mpm/examples/beam
  (:use :cl))
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)

(in-package :cl-mpm/examples/beam)
(declaim (optimize (debug 0) (safety 0) (speed 3)))


(defun find-max-cfl (sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps)
                   (dt cl-mpm:sim-dt))
      sim
    (let ((max-v (loop for mp across mps
                       maximize (with-accessors ((vel cl-mpm/particle:mp-velocity)) mp
                                  (magicl::sum (magicl:map #'abs vel))))))
      (* dt (/ max-v (cl-mpm/mesh:mesh-resolution mesh))))))

(defun pescribe-velocity (sim load-mps vel)
  (let ((mps (cl-mpm:sim-mps sim)))
    (loop for mp in load-mps
          do
             (progn
               (setf (cl-mpm/particle:mp-velocity mp) vel)))))
(defun increase-load (sim load-mps amount)
  (loop for mp in load-mps
        do (with-accessors ((pos cl-mpm/particle:mp-position)
                            (force cl-mpm/particle:mp-body-force)) mp
             (incf (magicl:tref force 0 0) amount))))

(defun length-from-def (sim mp dim)
  (let* ((mp-scale 2)
         (h-initial (magicl:tref (cl-mpm/particle::mp-domain-size mp) dim 0))
         ;(h-initial  (/ (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) mp-scale))
         )
    (* (magicl:tref (cl-mpm::mp-deformation-gradient mp) dim dim) h-initial)))
(defun max-stress (mp)
  (multiple-value-bind (l v) (magicl:eig (cl-mpm::voight-to-matrix (cl-mpm/particle:mp-stress mp)))
    (apply #'max l)
    ;(magicl:tref (cl-mpm/particle:mp-stress mp) 2 0)
    ))
(defun plot (sim &optional (plot :damage))
  (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  (multiple-value-bind (x y c stress-y lx ly e density)
    (loop for mp across (cl-mpm:sim-mps sim)
          collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
          collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
          collect (length-from-def sim mp 0) into lx
          collect (length-from-def sim mp 1) into ly
          collect (if (slot-exists-p mp 'cl-mpm/particle::damage) (cl-mpm/particle:mp-damage mp) 0) into c
          collect (if (slot-exists-p mp 'cl-mpm/particle::strain-energy-density) (cl-mpm/particle::mp-strain-energy-density mp) 0) into e
          collect (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)) into density
          ;; collect (cl-mpm/particle:mp-volume mp) into density
          collect (max-stress mp) into stress-y
          finally (return (values x y c stress-y lx ly e density)))
    (cond
      ((eq plot :damage)
       (vgplot:format-plot t "set cbrange [0:1]")
       ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 0.01 (apply #'max c)))
       (vgplot:plot x y c ";;with points pt 7 lc palette")
       )
      ((eq plot :energy)
       (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min e) (+ 0.01 (apply #'max e)))
       (vgplot:plot x y e ";;with points pt 7 lc palette")
       )
      (
       (eq plot :stress)
       (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 0.01 (apply #'max stress-y)))
       (vgplot:plot x y stress-y ";;with points pt 7 lc palette"))
      ((eq plot :deformed)
       ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 0.01 (apply #'max stress-y)))
       (vgplot:plot x y lx ly ";;with ellipses"))
      ((eq plot :density)
       (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min density) (+ 0.01 (apply #'max density)))
       (vgplot:plot x y e ";;with points pt 7 lc palette")
       ))
    )
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:axis (list 0 ms-x
                       0 ms-y))
    (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
    (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
      (vgplot:format-plot t "set ytics ~f" h)
      (vgplot:format-plot t "set xtics ~f" h))
  (vgplot:replot))

(defun remove-sdf (sim sdf)
      (setf (cl-mpm:sim-mps sim)
            (remove-if (lambda (mp)
                         (with-accessors ((pos cl-mpm/particle:mp-position)) mp
                           (>= 0 (funcall sdf pos))
                           ))
                       (cl-mpm:sim-mps sim))))
(defun ellipse-sdf (position x-l y-l)
  (let ((aspect (/ x-l y-l)))
    (lambda (pos)
      (let* ((position (magicl:from-list position '(2 1) :type 'double-float))
             (dist-vec (magicl:.* (magicl:.- position pos) (magicl:from-list (list 1 aspect) '(2 1)
                                                                             :type 'double-float)))
             (distance (sqrt (magicl:tref (magicl:@ (magicl:transpose dist-vec)
                                                    dist-vec) 0 0))))
        (- distance x-l)))))

(defun setup-test-column (size block-size block-offset &optional (e-scale 1d0) (mp-scale 1d0))
  (let* ((sim (cl-mpm/setup::make-block (/ 1 e-scale)
                                        (mapcar (lambda (s) (* s e-scale)) size)
                                        #'cl-mpm/shape-function:make-shape-function-linear)) 
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         ;(e-scale 1)
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (mass (/ (* 100 h-x h-y) (expt mp-scale 2)))
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size)))
    (progn
      (let ((block-position
              (mapcar #'+ (list (* h-x (- (+ (/ 1 (* 2 mp-scale))) 0))
                                (* h-y (+ (/ 1d0 (* 2d0 mp-scale)))))
                      block-offset)))
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-block-mps
               block-position
               block-size
               (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
               'cl-mpm::make-particle
               ;; 'cl-mpm/particle::particle-viscoelastic-fracture

               'cl-mpm/particle::particle-elastic-damage
               :E 1d7
               :nu 0.325d0

               ;; :E 1e6 :nu 1d8

               ;; 'cl-mpm/particle::particle-viscoplastic
               ;; :E 1e7
               ;; :nu 0.325
               ;; :visc-factor 1e6
               ;; :visc-power 3

               ;; Fluid
               ;; 'cl-mpm/particle::particle-fluid
               ;; :stiffness 1e5
               ;; :adiabatic-index 6
               ;; :rest-density 900d0
               ;; :viscosity 1.02e-3

               ;; :E 1e6 :nu 0.33
               :mass mass
               :critical-stress 1d6
               ;; :fracture-toughness 5d0
               :gravity -9.8d0
               )))
      (let ((prev-mps (cl-mpm:sim-mps sim)))
        (setf (cl-mpm:sim-mps sim) (make-array 1000 :fill-pointer 0 :adjustable t))
        (loop for mp across prev-mps
              do (vector-push-extend mp (cl-mpm:sim-mps sim))))
      (setf (cl-mpm:sim-damping-factor sim) 1d-3)
      (setf (cl-mpm:sim-mass-filter sim) 1d-12)
      (setf (cl-mpm:sim-dt sim) 1d-2)

      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-var (cl-mpm:sim-mesh sim)
                                            (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 0)))
                                            (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil)))
                                            (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0)))
                                            (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0)))
                                            ;; (lambda (i) (cl-mpm/bc:make-bc-friction i
                                            ;;                                         (magicl:from-list '(0d0 1d0)
                                            ;;                                                           '(2 1))
                                            ;;                                         0.25d0))
                                            ))

      ;; (let ((step-x 200)
      ;;       (step-y 400)
      ;;       )

      ;;   (loop for x from 0 to (round step-x h)
      ;;         do (loop for y from (- (round step-y h) 1) to (round step-y h)
      ;;                  do (push
      ;;                      (cl-mpm/bc:make-bc-fixed (list x y) '(nil 0))
      ;;                           (cl-mpm:sim-bcs sim))))
      ;;   )
      sim)))

;Setup
(defun setup ()
  (defparameter *sim* (setup-test-column '(700 600) '(500 100) '(0 400) (/ 1 50) 2))
  ;; (defparameter *sim* (setup-test-column '(1 1) '(1 1) '(0 0) 1 1))
  ;;(remove-sdf *sim* (ellipse-sdf (list 400 100) 10 40))
  ;; (remove-sdf *sim* (ellipse-sdf (list 1.5 3) 0.25 0.5))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  (defparameter *sim-step* 0)
  (defparameter *run-sim* nil)
  (defparameter *run-sim* t)
  (defparameter *load-mps*
    (let* ((mps (cl-mpm:sim-mps *sim*))
           (least-pos
              (apply #'max (loop for mp across mps
                                 collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))))
      (loop for id from 0 to (- (length mps) 1) when (>= (magicl:tref (cl-mpm/particle:mp-position (aref mps id)) 0 0) (- least-pos 0.001))
              collect (aref mps id))))
  ;; (increase-load *sim* *load-mps* 1)
  ;; (increase-load *sim* *load-mps* 100)
  )

(defun remove-material-damaged (sim)
  "Remove material points that have strain energy density exceeding fracture toughness"
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (delete-if (lambda (mp)
                 (with-accessors ((damage cl-mpm/particle:mp-damage))
                     mp
                   (>= damage 1d0))) mps)))


(defparameter *run-sim* nil)

(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)
  (defparameter *run-sim* t)
    (vgplot:close-all-plots)
    (vgplot:figure)
  (sleep 1)
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:axis (list 0 ms-x
                       0 ms-y))
    (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
    ;; (vgplot:axis (list 0 (nth 0 (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*))) 
    ;;                    0 (nth 1 (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))))
    (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
      ;; (vgplot:format-plot t "set ytics ~f" h)
      (vgplot:format-plot t "set xtics ~f" h)
      (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*))
                              *sim*)
    (time (loop for steps from 0 to 100
                while *run-sim*
                do
                (progn
                  (format t "Step ~d ~%" steps)
                  ;;; Adaptive timestepping
                  ;; (let* ((target-dt 1e-1)
                  ;;        (courent (find-max-cfl *sim*))
                  ;;        (c-value 0.001d0)
                  ;;        (cfl-dt (if (> courent 0d0) (/ c-value (/ courent (cl-mpm:sim-dt *sim*))) nil))
                  ;;        (new-dt (/ c-value (if (> courent 0d0)
                  ;;                               (/ courent (cl-mpm:sim-dt *sim*))
                  ;;                               (/ c-value 1e-6))))
                  ;;        (sub-steps (max (min (floor (/ target-dt new-dt)) 100) 1)))
                  ;;   (format t "C: ~f - steps: ~D - %dt: ~f~%" courent sub-steps new-dt)
                  ;;   (format t "Cfl derived dt:~f~%" cfl-dt)
                  ;;   (setf (cl-mpm:sim-dt *sim*) new-dt))
                  ;; (break)
                  (let ((max-cfl 0))
                    (time (dotimes (i 100)
                           ;; (pescribe-velocity *sim* *load-mps* (magicl:from-list '(0.5d0 0d0) '(2 1)))
                           (cl-mpm::update-sim *sim*)
                            ;; (cl-mpm/damage::calculate-damage (cl-mpm:sim-mesh *sim*)
                            ;;                               (cl-mpm:sim-mps *sim*)
                            ;;                               (cl-mpm:sim-dt *sim*)
                            ;;                               25d0
                            ;;                               )
                           (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))
                           ))
                    ;; (format t "Max cfl: ~f~%" max-cfl)
                    )
                  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*))
                                          *sim*)
                  (incf *sim-step*)
                  (plot *sim*)
                  ;; (vgplot:print-plot (asdf:system-relative-pathname "cl-mpm" (format nil "output/frame_~5,'0d.png" steps)))
                  (swank.live:update-swank)
                  (sleep .01)

                  )))
    ;; (vgplot:figure)
    ;; (vgplot:title "Velocity over time")
    ;; (vgplot:plot *time* *velocity*)
    ))

(setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))

(defun profile ()
  (sb-profile:unprofile)
  (sb-profile:reset)
  (sb-profile:profile cl-mpm::update-sim
                      cl-mpm::reset-grid
                      cl-mpm::p2g
                      cl-mpm::filter-grid
                      cl-mpm::update-nodes
                      cl-mpm::apply-bcs
                      cl-mpm::g2p
                      cl-mpm::update-particle
                      cl-mpm::update-stress
                      cl-mpm::iterate-over-neighbours-shape
                      cl-mpm::iterate-over-neighbours-shape-linear
                      cl-mpm::p2g-mp
                      cl-mpm::g2p-mp
                      cl-mpm::p2g-mp-node
                      cl-mpm::g2p-mp-node
                      cl-mpm/damage::calculate-damage
                      cl-mpm/damage::apply-damage
                      cl-mpm/damage::delocalise-damage
                      cl-mpm/damage::create-delocalisation-list
                      ;; cl-mpm/eigenerosion:update-fracture
                      ;; cl-mpm/eigenerosion::remove-material-damaged
                      ;; cl-mpm/eigenerosion::find-neighbours
                      )
  (loop repeat 100
        do (progn
             (cl-mpm::update-sim *sim*)
             (cl-mpm/damage::calculate-damage (cl-mpm:sim-mesh *sim*)
                                              (cl-mpm:sim-mps *sim*)
                                              (cl-mpm:sim-dt *sim*)
                                              25d0)
             ;; (cl-mpm/eigenerosion:update-fracture *sim*)
             ))
  (sb-profile:report))


(require 'sb-simd)
(declaim
 (inline simd-accumulate)
 (ftype (function ((simple-array double-float)
                   (simple-array double-float)
                   (simple-array double-float)
                   )
                  (values)) simd-accumulate))
;; (defun simd-mult (a b)
;;   (declare (type sb-simd:f64vec a b))
;;   (setf (sb-simd-avx:f64.2-aref a 0)
;;         (sb-simd-avx:f64.2+
;;          (sb-simd-avx:f64.2-aref a 0)
;;          (sb-simd-avx:f64.2-aref b 0)
;;          ))
;;   a
;;   )
(declaim
 (inline mult-simd)
 (ftype (function (magicl:matrix/double-float
                   magicl:matrix/double-float
                   magicl:matrix/double-float)
                  (values)) mult-simd))
(defun mult-simd (a b res)
  (declare (type magicl:matrix/double-float a b res))
  (let ((a-s (magicl::storage a))
        (b-s (magicl::storage b))
        (res-s (magicl::storage res))
        )
    (declare (type sb-simd:f64vec a-s b-s res-s))
    ;; (declare (type (simple-array double-float) a-s b-s res-s))
    (loop for i from 0 to 2
          do
             (setf (aref res-s i)
                   (sb-simd-avx:f64.2-horizontal+
                    (sb-simd-avx:f64.2*
                     (sb-simd-avx:f64.2-aref b-s 0)
                     (sb-simd-avx:f64.2-aref a-s (* 2 i)))
                    ))
             ;; (loop for j from 0 to 1
             ;;       do (incf (aref res-s i) (* (aref a-s (+ (* 2 i) j))
             ;;                                  (aref b-s j)
             ;;                                  )))
          )))
;; (declaim
;;  (inline mult)
;;  (ftype (function (magicl:matrix/double-float
;;                    magicl:matrix/double-float
;;                    magicl:matrix/double-float)
;;                   (values)) mult))
;; (defun mult (a b res)
;;   (declare (type magicl:matrix/double-float a b res))
;;   (let ((a-s (magicl::storage a))
;;         (b-s (magicl::storage b))
;;         (res-s (magicl::storage res))
;;         )
;;     (declare (type (simple-array double-float) a-s b-s res-s))
;;     (loop for i from 0 to 2
;;           do (loop for j from 0 to 1
;;                    do (incf (aref res-s i) (* (aref a-s (+ i (* 3 j)))
;; ;;                                               (aref b-s j)))))))
(let ((mp (cl-mpm/particle:make-particle 2 'cl-mpm/particle:particle :mass 1d0 :volume 1d0 :gravity 0d0 :body-force (magicl:zeros '(2 1))))
      (node (cl-mpm/mesh::make-node '(0 0) '(0 0)))
      (dsvp (cl-mpm/shape-function:assemble-dsvp-2d '(0d0 2d0)))
      (a (magicl:zeros '(3 2)))
      (b (magicl:zeros '(2 1)))
      (res (magicl:zeros '(2 1)))
      )
  (time
   (dotimes (i 1000000)
     (setf (cl-mpm/particle:mp-stress mp) (magicl:from-list '(1d0 1d0 1d0) '(3 1)))
     ;(cl-mpm::p2g-mp-node mp node 1d0 dsvp)
     ;; (cl-mpm::p2g-mp-node mp node 1d0 dsvp)
     ;; (cl-mpm::det-int-force mp dsvp res)
     (cl-mpm::mult-force dsvp (cl-mpm/particle:mp-stress mp) 1d0 res)
     ;; (cl-mpm::det-int-force mp dsvp res)
     )
   ))
;;  (let ((a (magicl:from-list '(1 1) '(2 1) :type 'double-float))
;;        (b (magicl:from-list '(1 0 0 -1 0.5 0.5) '(3 2) :type 'double-float))
;;        (res (magicl:zeros '(3 1) :type 'double-float)))
;;    (format t "Magicl ~%")
;;    (format t "~A ~%" (magicl:@ b a))
;;    (format t "lisp ~%")
;;    (mult b a res)
;;    (format t "~A ~%" res)
;;    (magicl:scale! res 0d0)
;;    (format t "lisp-simd ~%")
;;    (mult-simd (magicl:transpose b) a res)
;;    (format t "~A ~%" res)
;;    (let ((iter 10000000))
;;      (format t "Magicl ~%")
;;      (time
;;       (dotimes (i iter)
;;         (magicl:@ b a)))
;;      (format t "lisp ~%")
;;      (time
;;       (dotimes (i iter)
;;         (mult b a res)))
;;      ;; (format t "lisp-simd ~%")
;;      ;; (time
;;      ;;  (dotimes (i iter)
;;      ;;    (mult-simd (magicl:transpose b) a res)))
;;      )
;;    )