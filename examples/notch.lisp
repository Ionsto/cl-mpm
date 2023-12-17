(defpackage :cl-mpm/examples/notch
  (:use :cl))
;; (sb-ext:restrict-compiler-policy 'speed  3 3)
;; (sb-ext:restrict-compiler-policy 'debug  0 0)
;; (sb-ext:restrict-compiler-policy 'safety 0 0)
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)
;; (setf *block-compile-default* nil)
(in-package :cl-mpm/examples/notch)
(declaim (optimize (debug 3) (safety 3) (speed 0)))


(defun find-max-cfl (sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps)
                   (dt cl-mpm:sim-dt))
      sim
    (let ((max-v (loop for mp across mps
                       maximize (with-accessors ((vel cl-mpm/particle:mp-velocity)) mp
                                  (magicl::sum (magicl:map #'abs vel))))))
      (* dt (/ max-v (cl-mpm/mesh:mesh-resolution mesh))))))

(defun find-average-v (sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps)
                   (dt cl-mpm:sim-dt))
      sim
    (let ((sum-v (loop for mp across mps
                       sum (with-accessors ((vel cl-mpm/particle:mp-velocity)) mp
                             (the double-float (sqrt (the double-float (magicl::sum (magicl:.* vel vel)))))))))
      (/ sum-v (length mps)))))

(defun pescribe-velocity (sim load-mps vel)
  (let ((mps (cl-mpm:sim-mps sim)))
    (loop for mp in load-mps
          do
             (progn
               (loop for v in vel
                     for i from 0
                     do
                     (when v
                       (setf (magicl:tref (cl-mpm/particle:mp-velocity mp) i 0) v)))))))
(defun increase-load (sim load-mps amount)
  (loop for mp in load-mps
        do (with-accessors ((pos cl-mpm/particle:mp-position)
                            (force cl-mpm/particle:mp-body-force)) mp
             ;(incf (magicl:tref force 0 0) amount)
             (magicl:.+ force amount force)
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
  ;; (declare (opt (speed 2) (debug 3)))
  (multiple-value-bind (l v) (magicl:hermitian-eig (cl-mpm::voight-to-matrix (cl-mpm/particle:mp-stress mp)))
    ;; (/ (apply #'max l) 1d6)
    (/ (- (apply #'max l) (cl-mpm/particle::mp-pressure mp)) 1d6)
    ;; (cl-mpm/particle::mp-damage-ybar mp)
    ;; (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0)
    )
  ;; (let ((s (cl-mpm/utils::deviatoric-voigt (cl-mpm/particle:mp-stress mp))))
  ;;   (magicl:tref s 0 0))
  )
(defparameter *water-height* 0d0)
(defun plot (sim &optional (plot :damage))
  (declare (optimize (speed 0) (debug 3) (safety 3)))
  (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  (vgplot:title (format nil "Time ~,5f h" (/ *t* (* 60 60))))
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:format-plot t "set object 1 rect from 0,0 to ~f,~f fc rgb 'blue' fs transparent solid 0.5 noborder behind" ms-x *water-height*))
  (multiple-value-bind (x y c stress-y lx ly e density temp vx)
    (loop for mp across (cl-mpm:sim-mps sim)
          collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
          collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
          collect (magicl:tref (cl-mpm::mp-velocity mp) 0 0) into vx
          collect (length-from-def sim mp 0) into lx
          collect (length-from-def sim mp 1) into ly
          collect (if (slot-exists-p mp 'cl-mpm/particle::damage) (cl-mpm/particle:mp-damage mp) 0) into c
          collect (if (slot-exists-p mp 'cl-mpm/particle::temperature) (cl-mpm/particle:mp-temperature mp) 0) into temp
          collect (if (slot-exists-p mp 'cl-mpm/particle::strain-energy-density) (cl-mpm/particle::mp-strain-energy-density mp) 0) into e
          collect (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)) into density
          ;; collect (cl-mpm/particle:mp-volume mp) into density
          collect (max-stress mp) into stress-y
          finally (return (values x y c stress-y lx ly e density temp vx)))
    (let* ((node-x '())
           (node-y '())
           (node-c '())
           (mesh (cl-mpm:sim-mesh *sim*))
           (nodes (cl-mpm/mesh:mesh-nodes mesh))
           )
      (dotimes (i (array-total-size nodes))
        (let ((n (row-major-aref nodes i)))
          (with-accessors ((index cl-mpm/mesh:node-index)
                           (boundary cl-mpm/mesh::node-boundary-node)
                           (boundary-s cl-mpm/mesh::node-boundary-scalar)
                           (active cl-mpm/mesh::node-active)
                           )

              n
            (let ((n-ratio (/ (cl-mpm/mesh::node-volume n) (cl-mpm/mesh::node-volume-true n))))
              (when active
                (if boundary
                    ;; (print index)
                    (destructuring-bind (x y z) (cl-mpm/mesh:index-to-position mesh index)
                      ;; (push (nth 0 index) node-x)
                      ;; (push (nth 1 index) node-y)
                      (push x node-x)
                      (push y node-y)
                      (push 2 node-c)
                      ;; (push n-ratio node-c)
                      ;; (push boundary-s node-c)
                      )
                    (destructuring-bind (x y z) (cl-mpm/mesh:index-to-position mesh index)
                      (push x node-x)
                      (push y node-y)
                      ;; (push n-ratio node-c)
                      (push 0 node-c)
                      ;; (push boundary-s node-c)
                      ))
                ;; (push (/ (cl-mpm/mesh::node-volume n) (cl-mpm/mesh::node-volume-true n)) node-c)
                )))))
     (cond
        ((eq plot :point)
         ;; (vgplot:format-plot t "set cbrange [0:1]")
         ;; (vgplot:plot x y ";;with points pt 7")
         ;(vgplot:format-plot t "set cbrange [0:2]")
         (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min node-c) (+ 0.01 (apply #'max node-c)))
         (if node-x
             (multiple-value-bind (xp yp) (max-spacial (map 'list #'identity (cl-mpm::sim-mps *sim*)))
               (vgplot:plot
                ;; x y ";;with points pt 7"
                x y lx ly ";;with ellipses"
                node-x node-y node-c ";;with points pt 7 lc palette"
                (list xp) (list yp) ";;with points"
                ))
             (vgplot:plot x y ";;with points pt 7"))
         )
        ((eq plot :damage)
         (vgplot:format-plot t "set cbrange [0:1]")
         (vgplot:format-plot t "set style fill solid")
         (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 1d-10 (apply #'max c)))
         (vgplot:plot x y lx ly c ";;with ellipses lc palette")
         ;; (if node-x
         ;;     (vgplot:plot
         ;;      x y c ";;with points pt 7 lc palette"
         ;;      node-x node-y ";;with points pt 7")
         ;;     (vgplot:plot x y c ";;with points pt 7 lc palette"))
         )
        ((eq plot :velocity)
         (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min vx) (+ 0.01 (apply #'max vx)))
         (vgplot:plot x y vx ";;with points pt 7 lc palette")
         )
        ((eq plot :temperature)
         (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min temp) (+ 0.01 (apply #'max temp)))
         ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 0.01 (apply #'max c)))
         (vgplot:plot x y temp ";;with points pt 7 lc palette")
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
         )))
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
            (lparallel:premove-if (lambda (mp)
                         (with-accessors ((pos cl-mpm/particle:mp-position)) mp
                           (>= 0 (funcall sdf pos))
                           ))
                       (cl-mpm:sim-mps sim))))

(defun melt-sdf (sim sdf meltrate)
  (setf (cl-mpm:sim-mps sim)
        (lparallel:premove-if (lambda (mp)
                                (with-accessors ((pos cl-mpm/particle:mp-position)) mp
                                  (if t;(>= 0 (funcall sdf pos))
                                      (progn
                                        ;;Apply melt
                                        (with-accessors ((mass cl-mpm/particle:mp-mass)
                                                         (volume cl-mpm/particle:mp-volume)
                                                         (temp cl-mpm/particle::mp-boundary)
                                                         )
                                            mp
                                          (let ((density (/ mass volume)))
                                            ;;Apply some mass-meltrate
                                            (decf mass (* volume (abs (min 0 temp)) meltrate))
                                            ;;Keep density
                                            ;(setf volume (/ mass density))
                                            ;;If we overmelted, then remove particle
                                            (if (<= mass 1d2)
                                                t
                                                nil
                                                ))))
                                      nil)))
                              (cl-mpm:sim-mps sim))))

(defun damage-sdf (sim sdf)
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
      (loop for mp across mps
            do (with-accessors ((pos cl-mpm/particle:mp-position)
                                 (damage cl-mpm/particle:mp-damage)) mp
                  (when (>= 0 (funcall sdf pos))
                    (setf damage 1d0))))))

(defun rectangle-sdf (position size)
  (lambda (pos)
      (let* ((position (cl-mpm/utils:vector-from-list position))
             (dist-vec (magicl:.- (magicl:map! #'abs (magicl:.- pos position))
                                  (cl-mpm/utils:vector-from-list size))))
        (setf (magicl:tref dist-vec 2 0) 0d0)
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


(defun ice-init-stress (nu rho-ice rho-water ice-datum water-datum y g ice-height)
  (let ((h (- y (- ice-datum ice-height)))
        (hw (- ice-height (- ice-datum water-datum)))
        )
    ;; (format t "HW: ~f~%" hw)
    ;; (format t "h ~f~%" h)
    (cl-mpm/utils::stress-from-list (list
                                     (-
                                      (* -1d0 rho-water g (max 0d0 (- water-datum y)))
                                      0d0
                                      ;; (* (/ nu (- 1d0 nu))
                                      ;;   rho-ice g (- h (/ ice-height 2d0)))
                                      ;; 0d0
                                      ;; (* 0.5d0 rho-water g (/ (expt hw 2d0) ice-height))
                                      )
                                     (* -1d0 rho-ice g (- ice-height h))
                                     0d0
                                     ))))
(declaim (notinline setup-test-column))
(defun setup-test-column (size block-size block-offset &optional (e-scale 1d0) (mp-scale 1d0)
                          &rest mp-args)
  (let* ((sim (cl-mpm/setup::make-block (/ 1 e-scale)
                                        (mapcar (lambda (s) (* s e-scale)) size)
                                        :sim-type 'cl-mpm/damage::mpm-sim-damage
                                        ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         ;(e-scale 1)
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density *ice-density*)
         ;; (mass (/ (* 900 h-x h-y) (expt mp-scale 2)))
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size)))
    (progn
      (let ((block-position block-offset))
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-block-mps
               block-position
               block-size
 (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
               density
               ;; 'cl-mpm/particle::particle-viscoelastic
               ;; 'cl-mpm/particle::particle-elastic-damage
               ;; 'cl-mpm/particle::particle-elastic
               ;; 'cl-mpm/particle::particle-elastic-logspin
                'cl-mpm/particle::particle-viscoplastic-damage
               ;; 'cl-mpm/particle::particle-elastic
               :E 1d9
               :nu 0.3250d0
               :visc-factor 111d6
               :visc-power 3d0


               :initiation-stress 0.10d6
               :damage-rate 1d-6
               ;:critical-damage 0.56d0
               :critical-damage 0.5d0
               :local-length 100d0
               :local-length-damaged 0.01d0
               ;; :local-length-damaged 50d0
               :damage 0.0d0

               :gravity 9.8d0
               :gravity-axis (cl-mpm/utils:vector-from-list (list 0d0 -1d0 0d0))
               :index 0
               )))

      ;; (cl-mpm::split-mps-criteria
      ;;  sim
      ;;  (lambda (mp h)
      ;;    (when
      ;;        (or
      ;;         (> (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
      ;;              (- *ice-length* (* 1d0 h)))
      ;;         (< (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
      ;;            (+ (second block-offset) (* 1d0 h)))
      ;;         )
      ;;      :xy
      ;;      )))

      (let ((ms 1d7))
        (setf (cl-mpm:sim-damping-factor sim)
              (* 1d-5 ms)
              ;; 0d0
              ;; 0.5d0
              )
        (setf (cl-mpm::sim-mass-scale sim) ms))
      (setf (cl-mpm:sim-mass-filter sim) 1d-15);1d-15
      (setf (cl-mpm::sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) t)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) t)
      (setf (cl-mpm::sim-enable-damage sim) t)
      (setf (cl-mpm:sim-dt sim) 1d-2)
      (setf (cl-mpm:sim-bcs sim)
            (append
             (cl-mpm/bc::make-outside-bc-var
              (cl-mpm:sim-mesh sim)
              (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil nil)))
              (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil nil)))
              (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0 nil)))
              (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0 nil)))
              (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil nil 0)))
              (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil nil 0)))
              )))
      (let ((water-line 300))
        (loop for mp across (cl-mpm:sim-mps sim)
              do
                 (with-accessors ((pos cl-mpm/particle:mp-position)
                                  (damage cl-mpm/particle::mp-damage)
                                  (stress cl-mpm/particle::mp-stress-kirchoff)
                                  (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
                                  (stress-cauchy cl-mpm/particle::mp-stress)
                                  )
                     mp
                   (setf stress
                         ;(ice-init-stress nu rho-ice rho-water ice-datum water-datum y g ice-size)
                         ;; (ice-init-stress 0.325d0
                         ;;                  *ice-density*
                         ;;                  *water-density*
                         ;;                  (+ (second block-offset) (second block-size))
                         ;;                  water-line
                         ;;                  (magicl:tref pos 1 0)
                         ;;                  9.8d0
                         ;;                  (second block-size)
                         ;;                  )
                         (cl-mpm/utils:matrix-to-voight
                          (magicl:eye 3 :value (* 1d0 (cl-mpm/buoyancy::pressure-at-depth (magicl:tref pos 1 0) water-line *water-density*))))
                         stress-cauchy (magicl:scale stress 1d0)
                         undamaged-stress (magicl:scale stress 1d0)
                         damage (random 0.1d0)
                         )))

        (let ((ocean-x 1000d0)
              (ocean-y 300d0))
          (setf *water-height* ocean-y)
          (setf (cl-mpm::sim-bcs-force-list sim)
                (list
                 (cl-mpm/bc:make-bcs-from-list
                  (list
                   (cl-mpm/buoyancy::make-bc-buoyancy
                    sim
                    300d0
                    *water-density*
                    )))))))
      sim)))

(defparameter *ice-length* 1000)
(defparameter *ice-density* 850d0)
(defparameter *water-density* 1028d0)
(defparameter *average-window* 5)
;Setup
(defun setup (&optional (notch-length 100))
  (let* ((shelf-length *ice-length*)
         (shelf-height 200)
         (density-ratio (/ *ice-density* *water-density*))
         (shelf-bottom (- 300 (* density-ratio shelf-height)));;120
         (notch-length notch-length)
         (notch-depth 30);0
         (mesh-size 25)
         (mps-per-cell 2)
         (offset 00)
         )
    (format t "Shelf-bottom ~f~%" shelf-bottom)
    ;; (format t "Shelf-bottom ~a~%" shelf-bottom)
    (defparameter *sim* (setup-test-column (list (+ shelf-length 500) 500
                                                 )
                                           (list shelf-length shelf-height
                                                 )
                                           (list offset shelf-bottom
                                                 ) (/ 1 mesh-size) mps-per-cell))
    ;; (defparameter *sim* (setup-test-column (list (+ shelf-length 500) 500)
    ;;                                        (list shelf-length shelf-height)
    ;;                                        (list offset 50) (/ 1 mesh-size) mps-per-cell))
    ;;Bench calving
    (remove-sdf *sim* (rectangle-sdf (mapcar (lambda (x) (coerce x 'double-float))
                                             (list (+ shelf-length offset)
                                                   (+ shelf-height shelf-bottom)
                                                   0d0))
                                     (mapcar (lambda (x) (coerce x 'double-float))
                                             (list notch-length notch-depth 0d0)))))

  (format t "Simulation dt ~a~%" (cl-mpm:sim-dt *sim*))
  (format t "Simulation MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
  ;;Remove png snapshots 
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))

  (defparameter *velocity* '())
  (defparameter *t* 0)
  (defparameter *time* '())
  (defparameter *x* 0d0)
  (defparameter *x-pos* '())
  (defparameter *sim-step* 0)
  ;; (defparameter *run-sim* nil)
  (defparameter *load-mps*
    (let* ((mps (cl-mpm:sim-mps *sim*))
           (least-pos
              (apply #'max (loop for mp across mps
                                 collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))))
      (loop for mp across mps when (>= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (- least-pos 0.001))
              collect mp)))
  (defparameter *bottom-mps*
    (let* ((mps (cl-mpm:sim-mps *sim*))
           (least-pos
             (apply #'min (loop for mp across mps
                                collect (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)))))
      (loop for mp across mps when (<= (magicl:tref (cl-mpm/particle:mp-position mp) 1 0) (+ least-pos 0))
            collect mp)))

  )
(defun average-window (mps window)
  (let* ((x (loop for mp in mps collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
         (s (loop for mp in mps collect (magicl:tref
                                         (cl-mpm/utils::deviatoric-voigt
                                          (cl-mpm/particle:mp-stress mp)) 0 0)))
         (slen (length s))
         (sfirst (first s))
         (slast (first (last s)))
         )
    (loop repeat window
          do
             (progn
               (push sfirst s)
               (push slast (cdr (last s)))))
    (loop for i from 0 below slen
          collect
          (loop for iw from (- window) to window
                sum (/(nth (- (+ i window) iw) s) (+ 1 (* 2 window)))
                ))))
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

(defun plot-bottom-sxx (mps)
  ;; (vgplot:figure)
  (let* ((x (loop for mp in mps collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
         (s (loop for mp in mps collect (magicl:tref
                                                 (cl-mpm/utils::deviatoric-voigt
                                                  (cl-mpm/particle:mp-stress mp)) 0 0)))
         (x (mapcar (lambda (z) (- z *ice-length*)) x))
         (av (average-window mps *average-window*))
         ;; (av (loop for s1 in (append (list (first s)) s)
         ;;           for s2 in s
         ;;           for s3 in (append s (last s))
         ;;           collect (/ (+ s1 s2 s3) 3d0)))
         (xav x)
         (max-s (reduce #'max s))
         (ipos (position max-s s))
         (x-pos (nth ipos x))

         (max-s-av (reduce #'max av))
         (ipos-av (position max-s-av av))
         (x-pos-av (nth ipos-av xav))
        )
    (format t "Average max pos: ~A~%" (- x-pos-av))
    (format t "Real max pos: ~A~%" (- x-pos))
    (vgplot:plot x s "Stress"
                 xav av "Average"
                 (list x-pos) (list max-s) "Max;;with points pt 7"
                 (list x-pos-av) (list max-s-av) "Max;;with points pt 7")
  ))

(defparameter *run-sim* nil)

(defun sim-mass (sim)
  (reduce
   #'+
   (lparallel:pmapcar
    (lambda (mp) (with-accessors ((mass cl-mpm/particle:mp-mass)) mp
                   mass))
    (cl-mpm:sim-mps sim))))

(defun test-notch-length ()
  (defparameter *run-sim* t)
  (defparameter *start-mass* (sim-mass *sim*))
  (defparameter *notch-length* '())
  (defparameter *mass-loss* '())
  (vgplot:close-all-plots)
  (vgplot:figure)
  ;; (vgplot:plot '(10 10))
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
  (with-open-file (stream (merge-pathnames "output/notch_mass.csv") :direction :output :if-exists :supersede)
    (format stream "Notch Length,Mass Loss~%"))
  (defparameter *notch-position* 1d0)
  (loop for notch-length in '(10 20 30 40 50 100)
        do (progn
             (format t "Notch length ~A~%" notch-length)
             (setup notch-length)
             ;; (setup 100)
             ;; (defparameter *run-sim* t)
             (time (loop for steps from 0 to 40
                         while *run-sim*
                         do
                            (progn
                              (format t "Step ~d ~%" steps)
                              (time (dotimes (i 1000)
                                      (cl-mpm::update-sim *sim*)
                                      (remove-sdf *sim* (rectangle-sdf (list 1100 0 0) (list 300 1000 0)))
                                      (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))
                              (incf *sim-step*)
                              (plot *sim*)
                              (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                                 :terminal "png size 1920,1080"
                                                 )

                              (swank.live:update-swank)
                              (sleep .01)
                              )))
             ;; (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*))
                                     ;; *sim*)
             ;; (cl-mpm/output:save-csv (merge-pathnames (format nil "output/simcsv_~5,'0d.csv" *sim-step*)) *sim*)
             ;; (vgplot:figure)
             ;; (vgplot:title "Terminus over time")
             ;; (vgplot:plot *time* *x-pos*)

             (defparameter *end-mass* (sim-mass *sim*)))
           (let ((mass-loss (- *end-mass* *start-mass*)))
             (format t "Mass lost ~a~%" mass-loss)
             (push notch-length *notch-length*)
             (push mass-loss *mass-loss*)
             (with-open-file (stream (merge-pathnames "output/notch_mass.csv") :direction :output :if-exists :append)
               (format stream "~f, ~f ~%" notch-length mass-loss)))
        )

  (vgplot:figure)
  (vgplot:title "Mass loss vs notch length")
  (vgplot:plot *notch-length* *mass-loss*)
  )
(defun test-notch-length-max-stress ()
  (defparameter *run-sim* t)
  (defparameter *start-mass* (sim-mass *sim*))
  (defparameter *notch-length* '())
  (defparameter *mass-loss* '())
  (with-open-file (stream (merge-pathnames "output_notch/notch_txx.csv") :direction :output :if-exists :supersede)
    (format stream "length,txx,x,ydisp~%"))
  (vgplot:figure)
  (defparameter *notch-position* 0.1d0)
  (ensure-directories-exist "./output_notch/")
  (defparameter *run-sim* t)
  (loop for notch-length from 0 to 100 by 50
        do (progn
             (format t "Notch length ~A~%" notch-length)
             (setup notch-length)
             (cl-mpm/output:save-vtk-mesh (merge-pathnames "output_notch/mesh.vtk") *sim*)
             (setf (cl-mpm::sim-allow-mp-damage-removal *sim*) nil)
             (setf (cl-mpm::sim-enable-damage *sim*) nil)
             (let* ((target-time 1d0)
                    (substeps (floor target-time (cl-mpm::sim-dt *sim*))))
               (cl-mpm::update-sim *sim*)
               (let* ((dt-e (* 1d0 (cl-mpm::calculate-min-dt *sim*)))
                      (substeps-e (floor target-time dt-e)))
                 (format t "CFL dt estimate: ~f~%" dt-e)
                 (format t "CFL step count estimate: ~D~%" substeps-e)
                 (setf (cl-mpm:sim-dt *sim*) dt-e)
                 (setf substeps substeps-e))
               (format t "Substeps ~D~%" substeps)
               ;; (setf *run-sim* t)
               (setf *sim-step* 0)
               (time (loop for steps from 0 to 10
                           while *run-sim*
                           do
                           (progn
                             (format t "Step ~d ~%" steps)
                             (time (dotimes (i substeps)
                                     (cl-mpm::update-sim *sim*)
                                     (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))
                             (swank.live:update-swank)
                            (let ((av-v (find-average-v *sim*)))
                            (format t "~A~%" av-v)
                              )
                (let* ((slicer (floor (length *bottom-mps*) 2))
                       (test-mps (subseq *bottom-mps* slicer))
                       (x (loop for mp in test-mps collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
                       (s (average-window test-mps *average-window*))
                       ;; (s (loop for mp in test-mps collect (magicl:tref (cl-mpm/utils::deviatoric-voigt
                       ;;                                       (cl-mpm/particle:mp-stress mp)) 0 0)))
                       (ydisp (reduce #'max (loop for mp across (cl-mpm::sim-mps *sim*)
                                                  collect (magicl:tref (cl-mpm/particle::mp-displacement mp) 1 0))))
                       (max-s (reduce #'max s))
                       (ipos (position max-s s))
                       (x-pos (nth ipos x)))
                  (plot-bottom-sxx test-mps)
                  (format t "X conv: ~f~%" (- *ice-length* x-pos))
                  (format t "Y conv: ~f~%" ydisp)
                  )
                (incf *sim-step*)))))

                ;;Take the 2nd half of the domain - taking off the end-effect
                (let* ((slicer (floor (length *bottom-mps*) 2))
                        (test-mps (subseq *bottom-mps* slicer))
                        (x (loop for mp in test-mps collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
                        (s (loop for mp in test-mps collect (magicl:tref (cl-mpm/utils::deviatoric-voigt
                                                             (cl-mpm/particle:mp-stress mp)) 0 0)))

                       (ydisp (reduce #'max (loop for mp across (cl-mpm::sim-mps *sim*)
                                                  collect (magicl:tref (cl-mpm/particle::mp-displacement mp) 1 0))))
                        (max-s (reduce #'max s))
                        (ipos (position max-s s))
                        (x-pos (nth ipos x)))
                  (plot-bottom-sxx test-mps)
                  (format t "X pos final: ~f~%" (- *ice-length* x-pos))
                  (format t "Y pos final: ~f~%" ydisp)
				                 (with-open-file (stream (merge-pathnames "output_notch/notch_txx.csv") :direction :output :if-exists :append)
                 (format stream "~f, ~f, ~f, ~f ~%" notch-length max-s x-pos ydisp)))

             (cl-mpm/output:save-vtk (merge-pathnames (format nil "output_notch/sim_~a.vtk" notch-length)) *sim*)
             (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
             ;; (cl-mpm/output:save-csv (merge-pathnames (format nil "output_notch/final_~a.csv" notch-length)) *sim*)

             (let* ((test-mps *bottom-mps*)
                    (x (loop for mp in test-mps collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
                    (s (loop for mp in test-mps collect (magicl:tref (cl-mpm/utils::deviatoric-voigt
                                                                      (cl-mpm/particle:mp-stress mp)) 0 0))))
				       (with-open-file (stream (merge-pathnames (format nil "output_notch/bottom_txx_~a.csv" notch-length)) :direction :output :if-exists :supersede)
                 (format stream "x,t_xx~%")
                 (loop for sval in s
                       for xval in x
                       do (format stream "~f, ~f~%" xval sval))))

        )))
(defun plot-conv ()
  (vgplot:figure)
  (vgplot:plot
   *time* *conv-max-point* "Max point")
  (vgplot:figure)
  (vgplot:plot
   *time* *conv-vel* "Max vel")
  (vgplot:figure)
  (vgplot:plot
   *conv-max-point* *conv-vel* "vel-conv")
  )
(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)
  (defparameter *run-sim* t)
  (defparameter *start-mass* (sim-mass *sim*))
  (defparameter *conv-vel* (list))
  (defparameter *conv-max-point* (list))
  (vgplot:close-all-plots)
  ;; (sleep 1)
  ;; (vgplot:figure)
  ;; ;; (vgplot:plot '(10 10))
  ;; (sleep 1)
  ;; (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))
  ;;        (ms-x (first ms))
  ;;        (ms-y (second ms))
  ;;        )
  ;;   (vgplot:axis (list 0 ms-x
  ;;                      0 ms-y))
  ;;   (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
  ;;   (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
  ;;     (vgplot:format-plot t "set ytics ~f" h)
  ;;     (vgplot:format-plot t "set xtics ~f" h))
  (vgplot:figure)
  (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
    (vgplot:format-plot t "set xtics ~f" h))
  (with-open-file (stream (merge-pathnames "output/terminus_position.csv") :direction :output :if-exists :supersede)
    (format stream "Time (s),Terminus position~%")
    (loop for tim in (reverse *time*)
          for x in (reverse *x-pos*)
          do (format stream "~f, ~f ~%" tim x)))
  (defparameter *notch-position* 1d0)
  (let* ((target-time 1d4)
         (substeps (floor (/ target-time (cl-mpm::sim-dt *sim*))))
         (dt-scale 1d0))
    (cl-mpm::update-sim *sim*)
    (let* ((dt-e (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
           (substeps-e (floor target-time dt-e))
           )
      ;; (setf substeps-e (min 10000 substeps-e))
      (format t "CFL dt estimate: ~f~%" dt-e)
      (format t "CFL step count estimate: ~D~%" substeps-e)
      (setf (cl-mpm:sim-dt *sim*) dt-e)
      (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
    (time (loop for steps from 0 to 100
                while *run-sim*
                do
                   (progn
                     (if (> steps 1)
                         (setf (cl-mpm::sim-enable-damage *sim*) t)
                         (setf (cl-mpm::sim-enable-damage *sim*) nil))
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     ;; (cl-mpm/output:save-csv (merge-pathnames (format nil "output/sim_~5,'0d.csv" *sim-step*)) *sim*)
                     ;; (cl-mpm/output:save-csv (merge-pathnames (forma

                     (push *t* *time*)
                     ;; (setf *x*
                     ;;       (loop for mp across (cl-mpm:sim-mps *sim*)
                     ;;             maximize (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0)))
                     (setf *x* (sim-mass *sim*))
                     (push
                      *x*
                      *x-pos*)
                     (let ((max-cfl 0))
                       ;; (remove-sdf *sim* (rectangle-sdf (list 1000 330) (list *notch-position* 40)))
                       (time (loop repeat substeps
                                   while *run-sim*
                                   do
                                      (progn
                                        (cl-mpm::update-sim *sim*)
                                        ;(melt-sdf *sim* (rectangle-sdf (list 0 0) (list 1500 300)) (* (cl-mpm:sim-dt *sim*) 1e0))
                                        (incf *notch-position* (* (cl-mpm:sim-dt *sim*) 5d0))
                                        (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))))
                     (let* ((dt-e (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
                            (substeps-e (floor target-time dt-e))
                            )
                       (format t "CFL dt estimate: ~f~%" dt-e)
                       (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf (cl-mpm:sim-dt *sim*) dt-e)
                       (setf substeps substeps-e))
                     (let ((av-v (find-average-v *sim*)))
                       (format t "~A~%" av-v)
                       ;; (push av-v *conv-vel*)
                       )

                     ;; (let* ((dt-e (cl-mpm::calculate-min-dt *sim*))
                     ;;        (substeps-e (/ target-time dt-e))
                     ;;        )
                     ;;   (setf substeps-e (min 1000 substeps-e))
                     ;;   (format t "CFL dt estimate: ~f~%" dt-e)
                     ;;   (format t "CFL step count estimate: ~D~%" substeps-e)
                     ;;   (setf (cl-mpm:sim-dt *sim*) dt-e)
                     ;;   (setf substeps substeps-e)
                     ;;   )
                     (incf *sim-step*)
                     (plot *sim*)
                     (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                        :terminal "png size 1920,1080"
                                        )

                     ;; (plot-bottom-sxx *bottom-mps*)
                     ;; (vgplot:print-plot (merge-pathnames (format nil "output/bottom_dist_~5,'0d.png" *sim-step*)))
                     ;; (let* ((slicer (floor (length *bottom-mps*) 2))
                     ;;        (test-mps (subseq *bottom-mps* slicer))
                     ;;        (test-mps *bottom-mps*)
                     ;;        (x (loop for mp in test-mps collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
                     ;;        (s (loop for mp in test-mps collect (magicl:tref (cl-mpm/utils::deviatoric-voigt
                     ;;                                                          (cl-mpm/particle:mp-stress mp)) 0 0)))
                     ;;        (ydisp (reduce #'max (loop for mp across (cl-mpm::sim-mps *sim*)
                     ;;                                   collect (magicl:tref (cl-mpm/particle::mp-displacement mp) 1 0))))
                     ;;        (max-s (reduce #'max s))
                     ;;        (ipos (position max-s s))
                     ;;        (x-pos (nth ipos x))
                     ;;        )
                     ;;   (format t "Max at: ~f ~%" (- *ice-length* x-pos))
                     ;;   (format t "Max disp: ~f ~%" ydisp)
                     ;;   (push (- *ice-length* x-pos) *conv-max-point*))
                     (vgplot:replot)
                     (swank.live:update-swank)
                     (sleep .01)
                     ;; (with-open-file (stream (merge-pathnames "output/terminus_position.csv") :direction :output :if-exists :append)
                     ;;   (format stream "~f, ~f ~%" *t* *x*)
                     ;;   )

                     ))))
    (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*))
                                          *sim*)
  ;; (cl-mpm/output:save-csv (merge-pathnames (format nil "output/simcsv_~5,'0d.csv" *sim-step*)) *sim*)
    ;; (vgplot:figure)
    ;; (vgplot:title "Terminus over time")
    ;; (vgplot:plot *time* *x-pos*)

  (defparameter *end-mass* (sim-mass *sim*))
  (format t "Mass lost ~a~%" (- *end-mass* *start-mass*))
  (with-open-file (stream (merge-pathnames "output/terminus_position.csv") :direction :output :if-exists :supersede)
    (format stream "Time (s),Terminus position~%")
    (loop for tim in (reverse *time*)
          for x in (reverse *x-pos*)
          do (format stream "~f, ~f ~%" tim x)))
  )

(setf lparallel:*kernel* (lparallel:make-kernel 8 :name "custom-kernel"))

(defun profile ()
  (sb-profile:unprofile)
  (sb-profile:reset)
  (sb-profile:profile cl-mpm::update-sim
                      cl-mpm::reset-grid
                      cl-mpm::p2g
                      cl-mpm::filter-grid
                      cl-mpm::update-node-kinematics
                      cl-mpm::apply-bcs
                      cl-mpm::update-stress
                      cl-mpm::p2g-force
                      cl-mpm::update-forces
                      cl-mpm::g2p
                      cl-mpm::update-particle
                      cl-mpm/buoyancy::apply-non-conforming-nuemann
                      cl-mpm/buoyancy::locate-mps-cells
                      cl-mpm/buoyancy::apply-force-mps
                      cl-mpm/buoyancy::apply-force-cells
                      cl-mpm/damage::calculate-damage
                      cl-mpm/damage::apply-damage
                      cl-mpm/damage::delocalise-damage
                      cl-mpm/damage::create-delocalisation-list
                      )
  (loop repeat 1
        do (progn
             (cl-mpm::update-sim *sim*)
             ;; (cl-mpm/damage::calculate-damage (cl-mpm:sim-mesh *sim*)
             ;;                                  (cl-mpm:sim-mps *sim*)
             ;;                                  (cl-mpm:sim-dt *sim*)
             ;;                                  25d0)
             ;; (cl-mpm/eigenerosion:update-fracture *sim*)
             ))
  (sb-profile:report))

;; (let ((iters 1000000)
;;       (a (magicl:zeros '(2 1))))
(declaim (inline make-stress-vector))
(defun make-stress-vector ()
  (magicl::make-matrix/double-float 3 1 3 :column-major (make-array 3 :element-type 'double-float)))

(defun test-stress ()
  (with-accessors ((mps cl-mpm::sim-mps)
                   (mesh cl-mpm::sim-mesh)
                   (dt cl-mpm::sim-dt))
      *sim*
    (let ((m (magicl:zeros '(2 2))))
      (time
       (loop repeat 100
             do
                ;; (magicl::storage m)
                ;; (magicl::matrix/double-float-storage m)
                (cl-mpm::update-stress mesh mps dt)
             )))))
(defun test-mult ()
  (let ((iters 100000)
        (a (cl-mpm/utils:matrix-from-list '(1d0 0d0 0d0 2d0)))
        (b (cl-mpm/utils:matrix-from-list '(1d0 0d0 0d0 2d0)))
        (c (magicl:zeros '(2 1)))
        )
    (time
     (loop repeat iters
           do
              (progn (setf c (magicl:@ a b)))))
    (time
     (loop repeat iters
           do
              (progn (magicl:mult a b :target c))))
    (print (magicl:@ a b))
    (magicl:mult a b :target c)
    (print c)
    )
  )

;;   )
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
  (setup)
  (time
   (dotimes (i 100)
     (with-accessors ((mesh cl-mpm:sim-mesh)
                      (mps cl-mpm:sim-mps)
                      (dt cl-mpm:sim-dt))
         *sim*
       (cl-mpm::update-sim *sim*))))
  ;; (time-form 1000
  ;;            (progn
  ;;              (cl-mpm::update-sim *sim*)))
  )

(defun time-parts ()
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (dt cl-mpm::dt)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               )
                *sim*
                (time
                           (dotimes (i 100)
                             (cl-mpm::reset-grid mesh)
                             (cl-mpm::p2g mesh mps)
                             ;; (when (> mass-filter 0d0)
                             ;;   (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter *sim*)))
                             ;; (cl-mpm::update-node-kinematics mesh dt)
                             ;; (cl-mpm::apply-bcs mesh bcs dt)
                             ;; (cl-mpm::update-stress mesh mps dt)
                             ;;   (when enable-damage
                             ;;     (cl-mpm/damage::calculate-damage mesh
                             ;;                                      mps
                             ;;                                      dt
                             ;;                                      50d0))

                             ;;   ;Map forces onto nodes
                               (cl-mpm::p2g-force mesh mps)
                             ;;   (cl-mpm::apply-bcs mesh bcs-force dt)
                             ;;   (cl-mpm::update-node-forces mesh (cl-mpm::sim-damping-factor *sim*) dt)
                             ;; ;;   ;Reapply velocity BCs
                             ;;   (cl-mpm::apply-bcs mesh bcs dt)
                             ;; ;;   ;;Also updates mps inline
                             ;;   (cl-mpm::g2p mesh mps dt)

                             ;; ;;   (when remove-damage
                             ;; ;;     (cl-mpm::remove-material-damaged *sim*))
                             ;; ;;   (when split
                             ;; ;;     (cl-mpm::split-mps *sim*))
                             ;; (cl-mpm::check-mps mps)
                             )
                    )))
