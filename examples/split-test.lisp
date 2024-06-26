(defpackage :cl-mpm/examples/split-test
  (:use :cl))
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)
(in-package :cl-mpm/examples/split-test)
;; (pushnew :cl-mpm-pic *features*)
;; (asdf:compile-system :cl-mpm :force T)
;; (asdf:compile-system :cl-mpm :force T)

(defun get-disp (load-mps)
  ;; (* *t* *tip-velocity*)
  (- (/ (loop for mp in load-mps
              sum (+
                   (magicl:tref (cl-mpm/particle::mp-position mp) 0 0)
                   (* 0.5d0 (magicl:tref (cl-mpm/particle::mp-domain-size mp) 0 0))
                   )) (length load-mps))
     *initial-surface*
     ))

(defun max-v-sum (mp)
  (with-accessors ((vel cl-mpm/particle:mp-velocity))
      mp
    (magicl::sum (magicl:map #'abs vel))))

(defun find-max-cfl (sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps)
                   (dt cl-mpm:sim-dt))
      sim
    (let ((max-v (lparallel:preduce #'max
                                    (lparallel:pmapcar #'max-v-sum
                                                    mps))))
      (* dt (/ max-v (cl-mpm/mesh:mesh-resolution mesh))))))

(defun pescribe-velocity (sim load-mps vel)
  (let ((mps (cl-mpm:sim-mps sim)))
    (loop for mp in load-mps
          do
             (progn
               (loop for v in vel
                     for i from 0
                     do
                     (when v
                       (setf (magicl:tref (cl-mpm/particle:mp-velocity mp) i 0) v)
                       )
                     )))))

(defun increase-load (sim load-mps amount)
  (loop for mp in load-mps
        do (with-accessors ((pos cl-mpm/particle:mp-position)
                            (vol cl-mpm/particle::mp-volume)
                            (force cl-mpm/particle:mp-body-force)) mp
             ;(incf (magicl:tref force 0 0) amount)
             (magicl:.+ force (magicl:scale amount (/ 1d0 vol)) force)
             )))

(defun length-from-def (sim mp dim)
  (let* ((mp-scale 2)
         (h-initial (magicl:tref (cl-mpm/particle::mp-domain-size mp) dim 0))
         ;(h-initial  (/ (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) mp-scale))
         )
    h-initial
    ;(* (magicl:tref (cl-mpm::mp-deformation-gradient mp) dim dim) h-initial)
    ))
(defun max-stress (mp)
  (multiple-value-bind (l v) (magicl:hermitian-eig (cl-mpm::voight-to-matrix (cl-mpm/particle:mp-stress mp)))
    (apply #'max l)
    (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0)))
(defun max-spacial (mps)
  (let* ((x (loop for mp in mps collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
         (y (loop for mp in mps collect (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)))
         (s (loop for mp in mps collect (magicl:tref
                                         (cl-mpm/utils::deviatoric-voigt
                                          (cl-mpm/particle:mp-stress mp)) 0 0)))
         (max-s (reduce #'max s))
         (ipos (position max-s s))
         (x-pos (nth ipos x))
         (y-pos (nth ipos y))
         )
    (values x-pos y-pos))
  )
(declaim (notinline plot))
(defun plot (sim)
  (cl-mpm/plotter:simple-plot *sim* :plot :deformed
                                    :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp)
                                                   ;; (if
                                                                  ;; 1d0
                                                                  ;; 0d0)
                                                         ))
  ;; (vgplot:close-all-plots)
  ;; (plot-load-energy)
  )
;; (defun plot (sim &optional (plot :damage))
;;   (declare (optimize (speed 0) (debug 3) (safety 3)))
;;   (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
;;   (multiple-value-bind (x y c stress-y lx ly e density temp vx)
;;     (loop for mp across (cl-mpm:sim-mps sim)
;;           collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
;;           collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
;;           collect (magicl:tref (cl-mpm::mp-velocity mp) 0 0) into vx
;;           collect (length-from-def sim mp 0) into lx
;;           collect (length-from-def sim mp 1) into ly
;;           collect (if (slot-exists-p mp 'cl-mpm/particle::damage) (cl-mpm/particle:mp-damage mp) 0) into c
;;           collect (if (slot-exists-p mp 'cl-mpm/particle::temperature) (cl-mpm/particle:mp-temperature mp) 0) into temp
;;           collect (if (slot-exists-p mp 'cl-mpm/particle::strain-energy-density) (cl-mpm/particle::mp-strain-energy-density mp) 0) into e
;;           collect (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)) into density
;;           ;; collect (cl-mpm/particle:mp-volume mp) into density
;;           collect (max-stress mp) into stress-y
;;           finally (return (values x y c stress-y lx ly e density temp vx)))
;;     (let* ((node-x '())
;;            (node-y '())
;;            (node-c '())
;;            (mesh (cl-mpm:sim-mesh *sim*))
;;            (nodes (cl-mpm/mesh:mesh-nodes mesh))
;;            )
;;       (dotimes (i (array-total-size nodes))
;;         (let ((n (row-major-aref nodes i)))
;;           (with-accessors ((index cl-mpm/mesh:node-index)
;;                            (boundary cl-mpm/mesh::node-boundary-node)
;;                            (boundary-s cl-mpm/mesh::node-boundary-scalar)
;;                            (active cl-mpm/mesh::node-active)
;;                            )

;;               n
;;             (let ((n-ratio (/ (cl-mpm/mesh::node-volume n) (cl-mpm/mesh::node-volume-true n))))
;;               (when active
;;                 (if boundary
;;                     ;; (print index)
;;                     (destructuring-bind (x y) (cl-mpm/mesh:index-to-position mesh index)
;;                       ;; (push (nth 0 index) node-x)
;;                       ;; (push (nth 1 index) node-y)
;;                       (push x node-x)
;;                       (push y node-y)
;;                       (push 2 node-c)
;;                       ;; (push n-ratio node-c)
;;                       ;; (push boundary-s node-c)
;;                       )
;;                     (destructuring-bind (x y z) (cl-mpm/mesh:index-to-position mesh index)
;;                       (push x node-x)
;;                       (push y node-y)
;;                       ;; (push n-ratio node-c)
;;                       (push 0 node-c)
;;                       ;; (push boundary-s node-c)
;;                       ))
;;                 ;; (push (/ (cl-mpm/mesh::node-volume n) (cl-mpm/mesh::node-volume-true n)) node-c)
;;                 )))))
;;      (cond
;;         ((eq plot :point)
;;          ;; (vgplot:format-plot t "set cbrange [0:1]")
;;          ;; (vgplot:plot x y ";;with points pt 7")
;;          ;(vgplot:format-plot t "set cbrange [0:2]")
;;          (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min node-c) (+ 0.01 (apply #'max node-c)))
;;          (if node-x
;;              (multiple-value-bind (xp yp) (max-spacial (map 'list #'identity (cl-mpm::sim-mps *sim*)))
;;                (vgplot:plot
;;                 ;; x y ";;with points pt 7"
;;                 x y lx ly ";;with ellipses"
;;                 node-x node-y node-c ";;with points pt 7 lc palette"
;;                 (list xp) (list yp) ";;with points"
;;                 ))
;;              (vgplot:plot x y ";;with points pt 7"))
;;          )
;;         ((eq plot :damage)
;;          (vgplot:format-plot t "set cbrange [0:1]")
;;          (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 0.01 (apply #'max c)))
;;          (vgplot:plot x y c ";;with points pt 7 lc palette")
;;          ;; (if node-x
;;          ;;     (vgplot:plot
;;          ;;      x y c ";;with points pt 7 lc palette"
;;          ;;      node-x node-y ";;with points pt 7")
;;          ;;     (vgplot:plot x y c ";;with points pt 7 lc palette"))
;;          )
;;         ((eq plot :velocity)
;;          (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min vx) (+ 0.01 (apply #'max vx)))
;;          (vgplot:plot x y vx ";;with points pt 7 lc palette")
;;          )
;;         ((eq plot :temperature)
;;          (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min temp) (+ 0.01 (apply #'max temp)))
;;          ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 0.01 (apply #'max c)))
;;          (vgplot:plot x y temp ";;with points pt 7 lc palette")
;;          )
;;         ((eq plot :energy)
;;          (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min e) (+ 0.01 (apply #'max e)))
;;          (vgplot:plot x y e ";;with points pt 7 lc palette")
;;          )
;;         (
;;          (eq plot :stress)
;;          (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 0.01 (apply #'max stress-y)))
;;          (vgplot:plot x y stress-y ";;with points pt 7 lc palette"))
;;         ((eq plot :deformed)
;;          ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 0.01 (apply #'max stress-y)))
;;          (vgplot:plot x y lx ly ";;with ellipses"))
;;         ((eq plot :density)
;;          (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min density) (+ 0.01 (apply #'max density)))
;;          (vgplot:plot x y e ";;with points pt 7 lc palette")
;;          )))
;;     )
;;   (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
;;          (ms-x (first ms))
;;          (ms-y (second ms))
;;          )
;;     (vgplot:axis (list 0 ms-x
;;                        0 ms-y))
;;     (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
;;     (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
;;       (vgplot:format-plot t "set ytics ~f" h)
;;       (vgplot:format-plot t "set xtics ~f" h))
;;   (vgplot:replot))

(defun remove-sdf (sim sdf)
      (setf (cl-mpm:sim-mps sim)
            (remove-if (lambda (mp)
                         (with-accessors ((pos cl-mpm/particle:mp-position)) mp
                           (>= 0 (funcall sdf pos))
                           ))
                       (cl-mpm:sim-mps sim))))

(defun damage-sdf (sim sdf &optional (new-damage 1d0))
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
      (loop for mp across mps
            do (with-accessors ((pos cl-mpm/particle:mp-position)
                                 (damage cl-mpm/particle:mp-damage)) mp
                  (when (>= 0 (funcall sdf pos))
                    (setf damage new-damage))))))

(defun rectangle-sdf (position size)
  (lambda (pos)
      (let* ((position (magicl:from-list position '(2 1) :type 'double-float))
             (dist-vec (magicl:.- (magicl:map! #'abs (magicl:.- pos position))
                                  (magicl:from-list size '(2 1) :type 'double-float))))

        (+ (sqrt (magicl::sum
                  (magicl:map! (lambda (x) (* x x))
                               (magicl:map! (lambda (x) (max 0d0 x)) dist-vec))))
           (min (max (magicl:tref dist-vec 0 0)
                     (magicl:tref dist-vec 1 0)
                     ) 0d0)))))
(defun ellipse-sdf (position x-l y-l)
  (let ((aspect (/ x-l y-l)))
    (lambda (pos)
      (let* ((position (magicl:from-list position '(2 1) :type 'double-float))
             (dist-vec (magicl:.* (magicl:.- position pos) (magicl:from-list (list 1 aspect) '(2 1)
                                                                             :type 'double-float)))
             (distance (sqrt (magicl:tref (magicl:@ (magicl:transpose dist-vec)
                                                    dist-vec) 0 0))))
        (- distance x-l)))))


(defun setup-test-column (size block-size block-offset &optional (e-scale 1d0) (mp-scale 1d0)
                          &rest mp-args)
  (let* ((sim (cl-mpm/setup::make-block (/ 1 e-scale)
                                        (mapcar (lambda (s) (* s e-scale)) size)
                                        #'cl-mpm/shape-function:make-shape-function-bspline
                                        'cl-mpm::mpm-sim-usf
                                        ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         ;(e-scale 1)
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 2d3)
         ;; (mass (/ (* 900 h-x h-y) (expt mp-scale 2)))
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size)))
    (progn
      (let ((block-position
              block-offset))
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-block-mps
               block-position
               block-size
               (mapcar (lambda (e) (round (* e mp-scale) h-x)) block-size)
               density
                'cl-mpm/particle::particle-elastic
                :E 1d6
                :nu 0.30d0
                 :index 0
               )))
      (let ((mass-scale 1d0))
        (setf (cl-mpm::sim-mass-scale sim) mass-scale)
        ;; (setf (cl-mpm:sim-damping-factor sim)
        ;;       ;; (* 0.01d0 mass-scale)
        ;;       )
        )
      (setf (cl-mpm:sim-damping-factor sim) (* density 1d-3))
      (setf (cl-mpm:sim-mass-filter sim) 1d-15)
      (setf (cl-mpm::sim-allow-mp-split sim) t)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-nonlocal-damage sim) nil)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) t)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) t)
      (setf (cl-mpm:sim-dt sim) 1d-4)
      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-var
             (cl-mpm:sim-mesh sim)
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil 0)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil 0)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0 0)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0 0)))
             )
            )
      (defparameter *target-displacement* 0d0)
      (defparameter *initial-surface*
        (loop for mp across (cl-mpm:sim-mps sim)
              when (not (= 1 (cl-mpm/particle::mp-index mp)))
              maximize (magicl:tref
                        (magicl:.+ (cl-mpm/particle:mp-position mp)
                                   (magicl:scale (cl-mpm/particle::mp-domain-size mp) 0.5d0)
                                   )
                        0 0
                        )))
      (let ((mps (cl-mpm:sim-mps sim)))
        (let ((x-max (loop for mp across mps
                           maximize (magicl:tref
                                     (cl-mpm/particle:mp-position mp)
                                     0 0))))
          (defparameter *terminus-mps*
            (loop for mp across mps
                  when (= x-max (magicl:tref
                                 (cl-mpm/particle:mp-position mp)
                                 0 0))
                    collect mp))))
      ;; (setf (cl-mpm::sim-bcs-force-list sim)
      ;;       (list
      ;;        (cl-mpm/bc:make-bcs-from-list
      ;;         (list
      ;;          ;; *floor-bc*
      ;;          (cl-mpm/bc::make-bc-closure
      ;;           '(0 0 0)
      ;;           (lambda ()
      ;;             (with-accessors ((mesh cl-mpm:sim-mesh)
      ;;                              (dt cl-mpm::sim-dt)
      ;;                              )
      ;;                 sim
      ;;               (let ((datum (* -1d0 (+ *initial-surface* *target-displacement*)))
      ;;                     (normal (cl-mpm/utils:vector-from-list  '(-1d0 0d0 0d0))))
      ;;                 ;; (format t "Datum: ~F~%" datum)
      ;;                 (cl-mpm/penalty::apply-displacement-control-mps mesh (coerce *terminus-mps* 'vector )
      ;;                                                  dt
      ;;                                                  normal
      ;;                                                  datum
      ;;                                                  (* density 1d2)
      ;;                                                  0d0))
      ;;               )))
      ;;          ))))
      (defparameter *load-bc*
        (cl-mpm/buoyancy::make-bc-pressure
         sim
         0d0
         0d0
         ))
      (setf (cl-mpm::sim-bcs-force-list sim)
            (list
             (cl-mpm/bc:make-bcs-from-list
              (list
               *load-bc*
               ))))

      ;; (defparameter *shear-rate* .1d0)
      ;; (setf (cl-mpm:sim-bcs sim)
      ;;       (cl-mpm/bc:make-bcs-from-list
      ;;        (append
      ;;         (map 'list #'identity (cl-mpm:sim-bcs sim))
      ;;         (list
      ;;          (cl-mpm/bc::make-bc-closure
      ;;           '(0 0)
      ;;           (lambda ()
      ;;             ;; (apply-pullout
      ;;             ;;  sim
      ;;             ;;  *terminus-mps*
      ;;             ;;  *shear-rate*)
      ;;             )
      ;;           )))))
      (defparameter *pressure-inc-rate* 0d4)
      (defparameter *fatigue-load* 0d5)
      (defparameter *fatigue-period* 8d0)
      ;; (defparameter *load-bc*
      ;;   (cl-mpm/buoyancy::make-bc-pressure
      ;;    sim
      ;;    0.93d6
      ;;    ;; 0.20d6
      ;;    0d0
      ;;    0d0
      ;;    ))
      ;; (setf (cl-mpm::sim-bcs-force sim)
      ;;       (cl-mpm/bc:make-bcs-from-list
      ;;        (list *load-bc*)))

      sim)))

(defparameter *bar-length* 10d0)
;Setup
(defun setup ()
  (defparameter *run-sim* nil)
  (let ((mesh-size 1.0000)
        ;;At 2x2 we get 5630J
        ;;At 4x4 we get ??
        (mps-per-cell 2)
        (bar-length *bar-length*))
    (defparameter *sim* (setup-test-column (list 30 mesh-size)
                                           (list bar-length mesh-size) '(000 0) (/ 1 mesh-size) mps-per-cell))
    )
  ;; (remove-sdf *sim* (ellipse-sdf (list 1.5 3) 0.25 0.5))
  (print (cl-mpm:sim-dt *sim*))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  (defparameter *x* 0d0)
  (defparameter *x-pos* '())

  (defparameter *data-load* '())
  (defparameter *data-disp* '())
  (defparameter *max-stress* '())
  (defparameter *max-damage* '())
  (defparameter *max-x* '())
  (defparameter *energy-dissipation* '())

  (defparameter *cfl-max* '())
  (defparameter *sim-step* 0)
  ;(defparameter *run-sim* t)
  (defparameter *load-mps*
    (let* ((mps (cl-mpm:sim-mps *sim*))
           (least-pos
              (apply #'max (loop for mp across mps
                                 collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))))
      (loop for mp across mps when (>= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (- least-pos 0.001))
              collect mp)))
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (let ((x-min (loop for mp across mps
                       minimize (magicl:tref
                                 (cl-mpm/particle:mp-position mp)
                                 0 0))))
      (defparameter *far-field-mps*
        (loop for mp across mps
              when (= x-min (magicl:tref
                         (cl-mpm/particle:mp-position mp)
                         0 0))
                collect mp)))
    (let ((y-min (loop for mp across mps
                       minimize (magicl:tref
                                 (cl-mpm/particle:mp-position mp)
                                 1 0))))
      (defparameter *bottom-mps*
        (loop for mp across mps
              when (= y-min (magicl:tref
                             (cl-mpm/particle:mp-position mp)
                             1 0))
                collect mp)))
    (loop for mp in *terminus-mps*
          do (setf (cl-mpm/particle::mp-split-depth mp) 10));)
    ;; (increase-load *sim* *terminus-mps*
    ;;                (magicl:from-list (list (* 1d5) 0d0) '(2 1)))
    )
  ;; (increase-load *sim* *load-mps* 1)
  ;; (increase-load *sim* *load-mps* 100)
  )


(defparameter *run-sim* nil)
(defun calculate-dt (courent c-value target-step)
  (let* (
         (cfl-dt (if (> courent 0d0) (/ c-value (/ courent (cl-mpm:sim-dt *sim*))) nil))
         (new-dt (/ c-value (if (> courent 0d0)
                                (/ courent (cl-mpm:sim-dt *sim*))
                                (/ c-value 1e-6))))
         (max-steps 1000)
         (sub-steps (max (min (floor (/ target-step new-dt)) max-steps) 1)))
    (when (> (floor (/ target-step new-dt)) max-steps)
        (format t "CFL requires more steps than max-steps~%"))
    (format t "C: ~f - steps: ~D - %dt: ~f~%" courent sub-steps new-dt)
    (format t "Cfl derived dt:~f~%" cfl-dt)
                                        ;(setf (cl-mpm:sim-dt *sim*) new-dt)
    (values new-dt sub-steps)))


(defun apply-pullout (sim load-mps shear-rate)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (loop for mp in load-mps
          do
             (cl-mpm::iterate-over-neighbours
              mesh
              mp
              (lambda (mesh mp node svp grad fsvp fgrad)
                (with-accessors (
                                 ;; (vel cl-mpm/particle:mp-velocity)
                                 )
                    mp
                  (with-accessors ((pos cl-mpm/mesh::node-position)
                                   (vel cl-mpm/mesh::node-velocity)
                                   (active cl-mpm/mesh::node-active)
                                   (acc cl-mpm/mesh::node-acceleration)
                                   )
                      node
                    (when active
                      (setf (magicl:tref vel 0 0)
                            shear-rate
                            (magicl:tref vel 1 0)
                            0d0
                            )
                      )
                    )))))))

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
    (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
      (vgplot:format-plot t "set ytics ~f" h)
      (vgplot:format-plot t "set xtics ~f" h))
  ;; (with-open-file (stream (merge-pathnames "output/terminus_position.csv") :direction :output :if-exists :supersede)
  ;;   (format stream "Time (s),Terminus position~%")
  ;;   (loop for tim in (reverse *time*)
  ;;         for x in (reverse *x-pos*)
  ;;         do (format stream "~f, ~f ~%" tim x)))
  ;; (dotimes (i 1000)
  ;;   (cl-mpm::update-sim *sim*))

  (let* ((target-time 1d0)
         (dt (cl-mpm:sim-dt *sim*))
         (dt-scale 0.5d0)
         (substeps (floor target-time dt))
         (disp-step 0.1d0)
         (pressure-step 1d4)
         )
    (cl-mpm::update-sim *sim*)

    (let* ((dt-e (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
           (substeps-e (floor target-time dt-e)))
      (format t "CFL dt estimate: ~f~%" dt-e)
      (format t "CFL step count estimate: ~D~%" substeps-e)
      (setf (cl-mpm:sim-dt *sim*) dt-e)
      (setf substeps substeps-e))

    (format t "Substeps ~D~%" substeps)
    (time (loop for steps from 0 to 100
                while *run-sim*
                do
                   (progn
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     ;(cl-mpm/output:save-csv (merge-pathnames (format nil "output/sim_csv_~5,'0d.vtk" *sim-step*)) *sim*)

                     (push *t* *time*)
                     ;; (let ((cfl (find-max-cfl *sim*)))
                     ;;   (push cfl *cfl-max*))
                     (with-accessors ((mps cl-mpm:sim-mps))
                         *sim*
                       (setf *x*
                             (loop for mp across mps
                                   maximize (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0)))
                       (push
                        *x*
                        *x-pos*)

                       (push
                        (loop for mp across mps
                              maximize (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0))
                        *max-x*)
                       (push
                        (/ (loop for mp in *far-field-mps*
                                 sum
                                 (magicl:tref (cl-mpm/particle::mp-stress mp) 0 0))
                           (length *far-field-mps*))
                        *max-stress*)
                       (push
                        (/ (loop for mp across mps
                                 sum (cl-mpm/particle::mp-damage mp))
                           (length mps))
                        *max-damage*)
                       (let ((m-d (loop for mp across mps maximize (cl-mpm/particle::mp-damage mp))))
                         (when (>= m-d 1d0)
                           (setf *run-sim* nil)))
                       )

                     (format t "Core loop~%")
                     (let ((average-force 0d0)
                           (average-disp 0d0)
                           (intergral-energy 0d0)
                           )
                       (time (loop for i from 1 to substeps
                                   while *run-sim*
                                   do
                                      (progn
                                        ;; (incf (first (cl-mpm/buoyancy::bc-pressure-pressures *load-bc*))
                                        ;;       (* (cl-mpm:sim-dt *sim*) *pressure-inc-rate*))
                                        ;; (setf (first (cl-mpm/buoyancy::bc-pressure-pressures *load-bc*))
                                        ;;       (* *fatigue-load* (sin (/ (* *t* 2 3.14) *fatigue-period*))))
                                        (setf cl-mpm/penalty::*debug-force* 0d0)
                                        (cl-mpm::update-sim *sim*)
                                        (incf average-force (/
                                                             (/ cl-mpm/penalty::*debug-force*
                                                                -1.0d0
                                                                ;; (max 1 cl-mpm/penalty::*debug-force-count*)
                                                                )
                                                             substeps
                                                             ))
                                        (incf average-disp
                                              (/
                                               (get-disp *terminus-mps*)
                                               substeps
                                               )
                                              )

                                        (incf (nth 0 (cl-mpm/buoyancy::bc-pressure-pressures *load-bc*)) (/ pressure-step substeps))
                                        ;; (incf *target-displacement* (/ disp-step substeps))
                                        ;; (with-accessors ((mps cl-mpm:sim-mps))
                                        ;;     *sim*
                                        ;;   (loop for mp across mps
                                        ;;         when (>= (cl-mpm/particle::mp-damage mp) 1d0)
                                        ;;           do (setf *run-sim* nil)))
                                        ;; (setf cfl (max cfl (find-max-cfl *sim*)))
                                        (setf *t* (+ *t* (cl-mpm::sim-dt *sim*))))))
                       (when *run-sim*
                         (push average-force *data-load*)
                         (push average-disp *data-disp*)
                         ;; 
                         )
                       )

                     (format t "Target: ~f - Current: ~f Error: ~f~%"
                             *target-displacement*
                             (get-disp *terminus-mps*)
                             (* 100d0 (/ (- *target-displacement* (get-disp *terminus-mps*)) *target-displacement*)))
                     (format t "Total energy dissipated: ~F~%" (reduce #'+ *energy-dissipation*))

                       (let* ((dt-e (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
                              (substeps-e (floor target-time dt-e)))
                         (format t "CFL dt estimate: ~f~%" dt-e)
                         (format t "CFL step count estimate: ~D~%" substeps-e)
                         (setf (cl-mpm:sim-dt *sim*) dt-e)
                         (setf substeps substeps-e))
                     (incf *sim-step*)
                     (plot *sim*)
                     
                     (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*)))
                     (swank.live:update-swank)
                     (sleep .01)
                     ))))
  (format t "Finished~%")
  (format t "Total energy dissipated: ~F~%" (reduce #'+ *energy-dissipation*))
  ;; (plot-load-energy)
  (plot-load-disp)
  (with-open-file (stream (merge-pathnames "output/load-disp.csv") :direction :output :if-exists :supersede)
    (format stream "disp,load~%")
    (loop for x in (reverse *data-disp*)
          for l in (reverse *data-load*)
          do (format stream "~f, ~f ~%" x l)))
  ;; (vgplot:figure)
  ;; (vgplot:plot *time* *energy-dissipation* "Energy dissipation")
  ;; (plot-stress-damage-time)
  ;; (plot-creep-damage)
  )
(defun run-conv ()
  (defparameter *elements* '())
  (defparameter *ttf* '())
  (defparameter *energy* '())

  (defparameter *terminus-mps* '())
  (loop for i from 2 to 8
        do
           (let ((elements (expt 2 i))
                 (final-time 10)
                 (energy 0d0)
                 (failure-time 0d0)
                 (time (list))
                 (s-x (list))
                 (damage-x (list))
                 (bar-length 5)
                 )
             (let ((mesh-size (/ 5d0 elements))
                   (mps-per-cell 4))
               (defparameter *sim* (setup-test-column (list 10 mesh-size)
                                                      (list bar-length mesh-size) '(0 0) (/ 1 mesh-size) mps-per-cell))
               (damage-sdf *sim* (rectangle-sdf (list (/ bar-length 2) 0)
                                                (list
                                                 0.2
                                                 mesh-size)) 0.10d0)
               (with-accessors ((mps cl-mpm:sim-mps))
                   *sim*
                 (let ((x-max (loop for mp across mps
                                    maximize (magicl:tref
                                              (cl-mpm/particle:mp-position mp)
                                              0 0))))
                   (defparameter *terminus-mps*
                     (loop for mp across mps
                           when (= x-max (magicl:tref
                                          (cl-mpm/particle:mp-position mp)
                                          0 0))
                             collect mp))))
               )
             (defparameter *run-sim* t)
             ;(defparameter *sim* (setup-test-column '(1 60) '(1 50) (/ elements 50) 2))
             (setf (cl-mpm:sim-dt *sim*) (* 1d-3 (/ 16 elements)))
             (format t "Running sim size ~a ~%" elements)
             (format t "Sim dt: ~a ~%" (cl-mpm:sim-dt *sim*))
             (format t "Sim steps: ~a ~%" (/ final-time (cl-mpm:sim-dt *sim*)))
             (time
              (loop for steps from 0 to (round (/ final-time (cl-mpm:sim-dt *sim*)))
                    when *run-sim*
                    do
                       (progn
                         (cl-mpm::update-sim *sim*)
                         (with-accessors ((mps cl-mpm:sim-mps ))
                             *sim*
                           (incf
                            energy
                            (loop for mp across mps
                                       sum
                                       (with-accessors ((inc cl-mpm/particle::mp-damage-increment)
                                                        (vol cl-mpm/particle::mp-volume)) mp
                                         (* inc vol)
                                         )))
                           (loop for mp across mps
                                 when (>= (cl-mpm/particle::mp-damage mp) 1d0)
                                   do (setf *run-sim* nil)))
                         (incf failure-time (cl-mpm:sim-dt *sim*))
                         )))
             (plot *sim*)
             (sleep .01)
             (push failure-time
                   *ttf*)
             (push energy
                   *energy*)
             (push elements
                   *elements*)
             (format t "Failure at ~a with energy dissipation ~a ~%" failure-time energy)
             (cl-mpm/output:save-csv (merge-pathnames (format nil "conv_files/elements_~d.csv" elements)) *sim*)
             (cl-mpm/output:save-vtk (merge-pathnames (format nil "conv_files/elements_~d.vtk" elements)) *sim*)
             )
           (vgplot:figure)
           (vgplot:plot *elements* *ttf* "Time to failure")
           (vgplot:figure)
           (vgplot:plot *elements* *energy* "Energy to failure")
        ))

(defun estimate-gf ()
  (let* ((strain (mapcar (lambda (x) (/ x *bar-length*)) *data-disp*))
        (stress (mapcar (lambda (x) (/ x (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))) *data-load*))
         (dstrain (mapcar #'- (butlast strain) (cdr strain)))
         (av-stress (mapcar #'+ (butlast stress) (cdr stress)))
         (energy (* 0.5d0 (reduce #'+ (mapcar #'* dstrain av-stress))))
        )
    (vgplot:figure)
    (vgplot:title "Stress-strain")
    (vgplot:xlabel "Strain")
    (vgplot:ylabel "Stress (Pa)")
    (vgplot:plot strain stress)
    (format t "Estimated fracture energy G_f:~F~%" energy)
  ))
(defun plot-velocity ()
  (let ((x (loop for mp in *bottom-mps*
                 collect (magicl:tref
                          (cl-mpm/particle:mp-position mp)
                          0 0)))
        (v (loop for mp in *bottom-mps*
                 collect (magicl:tref
                          (cl-mpm/particle:mp-velocity mp)
                          0 0)))

        )
    (vgplot:figure)
    (vgplot:axis '(t t 0 t))
    (vgplot:plot
     x v "velocity"
     )))
(defun plot-damage ()
    (let ((x (loop for mp in *bottom-mps*
                   collect (magicl:tref
                            (cl-mpm/particle:mp-position mp)
                            0 0)))
          (d (loop for mp in *bottom-mps*
                   collect (cl-mpm/particle:mp-damage mp)))
          )
      (vgplot:figure)
      (vgplot:axis '(t t 0 1))
      (vgplot:plot
       x d "damage"
       )))
(defun plot-load-disp ()
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (let* ()
      (vgplot:figure)
      (vgplot:xlabel "Displacement")
      (vgplot:ylabel "Load")
      (vgplot:plot
       (mapcar (lambda (x) (* 1d0 x)) *data-disp*) *data-load*
       ;(mapcar (lambda (x) (/ x 1d0)) *time*) *max-damage* "MPM"
       ;; (mapcar (lambda (x) (/ x 1d0)) *time*) *max-damage* "MPM"
       )
      )
    ))
(defun plot-load-energy ()
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (let* ()
      (vgplot:figure)
      (vgplot:xlabel "Displacement")
      (vgplot:ylabel "Normalised load-energy")
      (vgplot:plot
       (mapcar (lambda (x) (* 1d0 x)) *data-disp*) (mapcar (lambda (x) (/ x (reduce #'max *data-load*))) *data-load*) "load"
       (mapcar (lambda (x) (* 1d0 x)) *data-disp*) (mapcar (lambda (x) (/ x (reduce #'max *energy-dissipation*))) *energy-dissipation*) "energy"
                                        ;(mapcar (lambda (x) (/ x 1d0)) *time*) *max-damage* "MPM"
       ;; (mapcar (lambda (x) (/ x 1d0)) *time*) *max-damage* "MPM"
       )
      (vgplot:figure)
      (vgplot:title "Load-disp")
      (vgplot:xlabel "Displacement")
      (vgplot:ylabel "Load")
      (vgplot:plot (mapcar (lambda (x) (* 1d0 x)) *data-disp*) *data-load* "load")
      (vgplot:figure)
      (vgplot:title "Load-energy")
      (vgplot:xlabel "Displacement")
      (vgplot:ylabel "Energy")
      (vgplot:plot (mapcar (lambda (x) (* 1d0 x)) *data-disp*) *energy-dissipation* "energy")
      (vgplot:figure)
      (vgplot:title "Stress-strain")
      (vgplot:xlabel "Strain")
      (vgplot:ylabel "Stress Pa")
      (vgplot:plot (mapcar (lambda (x) (/ x 5d0)) *data-disp*)
                   (mapcar (lambda (x) (/ x (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))) *data-load*) "Load")
      )
    ))

(defun plot-creep-damage ()
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (let* ((s-y 0.2d6)
          (s-x (first *max-x*))
          (len 500d0)
          (max-time (reduce #'max *time*))
          (df (lisp-stat:read-csv
	             (uiop:read-file-string #P"creep.csv")))
          (max-time-data (reduce #'max (lisp-stat:column df 'time)))
          )
      (vgplot:figure)
      (vgplot:xlabel "Normalised time to failure")
      (vgplot:ylabel "Damage")
      (vgplot:plot
       (mapcar (lambda (x) (/ x max-time)) *time*) *max-damage* "MPM"
       (aops:each (lambda (x) (/ x max-time-data)) (lisp-stat:column df 'time)) (lisp-stat:column df 'damage) "data - 0.93"
       )
      )
    ))
(defun plot-stress-damage-time ()
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (let ((s-y 0.2d6)
          (s-x (first *max-x*))
          (len 500d0)
          (df (lisp-stat:read-csv
	             (uiop:read-file-string #P"creep.csv"))))
      (vgplot:figure)
      (vgplot:plot
       *time* (mapcar (lambda (x) (/ x s-x)) *max-x*) "Far field extension"
       *time* (mapcar (lambda (x) (/ x s-y)) *max-stress*) "max stress"
       *time* *max-damage* "max damage"
       )
      (vgplot:figure)
      (vgplot:plot
       *time* *max-damage* "max damage"
       ;; (aops:each (lambda (x) (* x (/ 150 111))) (lisp-stat:column df 'time) ) (lisp-stat:column df 'damage) "0.93"
       )
      (vgplot:figure)
      (vgplot:plot
       (mapcar (lambda (x) (/ x len)) *max-x*) (mapcar (lambda (x) (/ x s-y)) *max-stress*) "max stress"
       (mapcar (lambda (x) (/ x len)) *max-x*) *max-damage* "max damage"
       )
      )
    ))
(defun plot-s-xx (mps H)
  (let* (
         (rho-ice 900)
         (E (cl-mpm/particle::mp-e (first mps)))
         (nu (cl-mpm/particle::mp-nu (first mps)))
         (g (cl-mpm/particle::mp-gravity (first mps)))
         (s-xx (loop for mp in mps
                    collect
                    (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0)))
         (y (loop for mp in mps
                    collect
                    (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)))
         (s-an (mapcar
                (lambda (y)
                  (* (/ nu (- 1 nu)) rho-ice g -1d0 (- y (/ H 2d0))))
                y)))
    (vgplot:figure)
    (vgplot:title "S_{xx} over height")
    (vgplot:xlabel "Height (m)")
    (vgplot:ylabel "Longitudinal stress (MPa)")
    (vgplot:plot y (mapcar (lambda (x) (* 1d-6 x)) s-xx) "s_{xx} mpm"
                 y (mapcar (lambda (x) (* 1d-6 x)) s-an) "s_{xx} analytic")))
(defun plot-cfl ()
  (vgplot:figure)
  (vgplot:title "CFL over time")
  (vgplot:plot *time* *cfl-max*))

(setf lparallel:*kernel* (lparallel:make-kernel 8 :name "custom-kernel"))

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
                      ;; cl-mpm::update-strain-kirchoff
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
             ;; (cl-mpm/damage::calculate-damage (cl-mpm:sim-mesh *sim*)
             ;;                                  (cl-mpm:sim-mps *sim*)
             ;;                                  (cl-mpm:sim-dt *sim*)
             ;;                                  25d0)
             ;; (cl-mpm/eigenerosion:update-fracture *sim*)
             ))
  (sb-profile:report))

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
(defun simple-time ()
  (setf lparallel:*kernel* (lparallel:make-kernel 8 :name "custom-kernel"))
  (let ((mp (cl-mpm/particle:make-particle 2 'cl-mpm/particle:particle
                                           :size (magicl:from-list '(1d0 1d0) '(2 1))))
        (mesh (cl-mpm::sim-mesh *sim*)))
  (time-form 1000
     (with-accessors ((mesh cl-mpm:sim-mesh)
                      (mps cl-mpm:sim-mps)
                      (dt cl-mpm:sim-dt))
         *sim*
       (cl-mpm::update-sim *sim*)
         ))))
;(let ((stress (magicl:zeros '(3 1)))
;      (strain (magicl:zeros '(3 1)))
;      (de (cl-mpm/constitutive::linear-elastic-matrix 1d0 0d0))
;      )
;  (format t "mat")
;  (time
;   (dotimes (i 100000)
;       (cl-mpm/constitutive::maxwell-exp strain stress 1d0 0d0 1d0 1d0)))
;  (format t "voight")
;  (time
;   (dotimes (i 100000)
;     (cl-mpm/constitutive::maxwell-exp-v strain stress 1d0 0d0 1d0 1d0)))
;  (format t "simd")
;  (time
;   (dotimes (i 100000)
;     (cl-mpm/constitutive::maxwell-exp-v-simd strain stress 1d0 0d0 de 1d0  1d0)))
;  )

;; (defun dot (x)
;;   (sqrt (magicl::sum (magicl:.* x x))))
;; (defun mp-sdf (mp x)
;;   (with-accessors ((pos cl-mpm/particle:mp-position)
;;                    (size cl-mpm/particle::mp-domain-size))
;;       mp
;;     (let ((r (max (magicl:tref size 0 0) (magicl:tref size 1 0))))
;;       (- (dot (magicl:.- pos x)) r))))

;; (defun draw-state (file)
;;   (declare (optimize (safety 3) (debug 3)))
;;   (let* ((resolution 10)
;;          (size (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
;;          (width (round (first size) resolution))
;;          (height (round (second size) resolution))
;;          (png (make-instance 'zpng:png
;;                               :color-type :truecolor
;;                               :width width
;;                               :height height))
;;          (image (zpng:data-array png))
;;          (max 255))
;;     (lparallel:pdotimes (y height)
;;       (dotimes (x width)
;;         ;; (print (type-of image))
;;         (setf (aref image y x 0) 255
;;               )
;;         (let ((d 1e10)
;;               (ipos (magicl:scale (magicl:from-list (list x y) '(2 1) :type 'double-float) resolution)))
;;           (loop for mp across (cl-mpm:sim-mps *sim*)
;;                 do
;;                    (setf d (min d (mp-sdf mp ipos)))
;;                 )
;;           (when (< d 0d0)
;;             (setf (aref image y x 0) 255
;;                   (aref image y x 1) 0
;;                   (aref image y x 2) 0)
;;             )))
;;       (zpng:write-png png file)
;;       )))
