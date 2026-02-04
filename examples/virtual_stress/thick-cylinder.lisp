(defpackage :cl-mpm/examples/virtual-stress/thick-cylinder
  (:use :cl))
(in-package :cl-mpm/examples/virtual-stress/thick-cylinder)
(sb-ext:restrict-compiler-policy 'speed  3 3)
(sb-ext:restrict-compiler-policy 'debug  0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
;; (sb-ext:restrict-compiler-policy 'speed  0 0)
;; (sb-ext:restrict-compiler-policy 'debug  3 3)
;; (sb-ext:restrict-compiler-policy 'safety 3 3)

(declaim (optimize (debug 3) (safety 3) (speed 2)))



;; (defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-elastic) dt fbar)
;;   (cl-mpm::update-stress-linear mesh mp dt fbar))

(declaim (notinline plot))
(defun plot (sim &optional (plot :deformed))
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:axis (list 0 ms-x
                       0 ms-y))
    (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
  (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim))))
    (vgplot:format-plot t "set ytics ~f" h)
    (vgplot:format-plot t "set xtics ~f" h))
  (cl-mpm/plotter::simple-plot
   sim
   :plot :deformed
   ;:colour-func (lambda (mp) (compute-radial-stress mp))
   :colour-func (lambda (mp) (abs (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :xy)))
   ;; :contact-bcs *penalty-bc*
   ))

(declaim (notinline setup-test-column))
(defun setup-test-column (&key (e-scale 1) (mp-scale 1)
                            (sym nil)
                            )
  (let* ((nd 2)
         (a 1d0)
         (b 5d0)
         (offset-scalar 10d0)
         (size (list (* offset-scalar 2) (* offset-scalar 2)))
         (block-size (list (* b 2) (* b 2)))
         (offset (list (- offset-scalar b)
                       (- offset-scalar b)))
         ;; (offset-scalar 0d0)
         ;; (size (list (+ b 1) (+ b 1)))
         ;; (block-size (list b b))
         ;; (offset (list offset-scalar
         ;;               offset-scalar))
         )
    (when sym
      (setf
       offset-scalar 0d0)
      (setf
       size (list (+ b 1) (+ b 1))
       block-size (list b b)
       offset (list offset-scalar
                     offset-scalar)))

    (defparameter *origin* (cl-mpm/utils:vector-from-list (list offset-scalar offset-scalar 0d0)))
    (let* ((sim
             (cl-mpm/setup::make-simple-sim;
              (/ 1d0 e-scale)
              (mapcar (lambda (x) (* x e-scale)) size)
              :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
              :args-list (list
                          :enable-fbar nil
                          :enable-aggregate t
                          :mp-removal-size nil
                          :enable-split nil)))
           (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
           (h-x (/ h 1d0))
           (h-y (/ h 1d0))
           (density 1d3)
           (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size)))
      (progn
        (let ((block-position (list (* h-x (+ 0 (/ 1d0 (* 2d0 mp-scale)) ));mp-scale->1
                                    (* h-y (+ 0 (/ 1d0 (* 2d0 mp-scale)) )))))
          (cl-mpm::add-mps
           sim
           (cl-mpm/setup::make-block-mps
            offset
            block-size
            (mapcar (lambda (e mp-s)
                      (round  (* e mp-s) h-x)) block-size
                        (list
                         mp-scale
                         mp-scale
                         1
                         ))
            density
            'cl-mpm/particle::particle-linear-elastic
            :E 1d9
            :nu 0.2d0)))
        (let ((refine 0))
          (cl-mpm/setup:remove-sdf
           sim
           (lambda (p)
             (funcall
              (cl-mpm/setup::circle-sdf
               (list offset-scalar offset-scalar 0d0)
               a)
              p))
           :refine refine)
          (cl-mpm/setup:remove-sdf
           sim
           (lambda (p)
             (- (funcall
                 (cl-mpm/setup::circle-sdf
                  (list offset-scalar
                        offset-scalar
                        0d0)
                  b)
                 p)))
           :refine refine))
        (format t "MP count ~D~%" (length (cl-mpm:sim-mps sim)))
        (setf (cl-mpm::sim-gravity sim) 0d0)
        (cl-mpm/setup::setup-bcs
         sim
         :bottom '(nil 0 nil)
         :left '(0 nil nil))
        (defparameter *bc-pressure*
          (cl-mpm/buoyancy::make-bc-pressure
           sim
           100d3
           100d3
           :clip-func (lambda (pos)
                        (<=
                         (funcall
                          (cl-mpm/setup::circle-sdf
                           (list offset-scalar offset-scalar 0d0)
                           (* 0.5d0 (+ a b)))
                          pos)
                         0d0))))
        (cl-mpm:add-bcs-force-list
         sim
         *bc-pressure*)
        sim))))
(defun setup (&key (refine 0d0) (mps 2d0) (sym nil))
  ;; (defparameter *sim* (setup-test-column '(1 60) '(1 50) (/ 1 5) 2))
  (let* ((e (expt 2 refine))
         (h (/ 1d0 e)))
    (format t "H:~E~%" h)
    (defparameter
        *sim*
      (setup-test-column
      :e-scale (/ 1d0 h)
                         :mp-scale mps
                         :sym sym
                         )))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defun compute-hoop-stress (mp)
  (let ((pos (cl-mpm/particle::mp-position mp))
        (stress (cl-mpm/particle::mp-stress mp)))
    (let* ((x (magicl:tref pos 0 0))
           (y (magicl:tref pos 1 0))
           (r (sqrt (+ (* x x) (* y y))))
           (sx (cl-mpm/utils:get-stress stress :xx))
           (sy (cl-mpm/utils:get-stress stress :yy))
           (sxy (cl-mpm/utils:get-stress stress :xy))
           (theta (atan y x))
           )
      (+ (* sx (expt (cos theta) 2))
         (* sy (expt (sin theta) 2))
         (* 2 sxy (sin theta) (cos theta))))))

(defun compute-radial-stress (mp)
  (let ((pos (cl-mpm/particle::mp-position mp))
        (stress (cl-mpm/particle::mp-stress mp)))
    (let* ((x (magicl:tref pos 0 0))
           (y (magicl:tref pos 1 0))
           (r (sqrt (+ (* x x) (* y y))))
           (sx (cl-mpm/utils:get-stress stress :xx))
           (sy (cl-mpm/utils:get-stress stress :yy))
           (sxy (cl-mpm/utils:get-stress stress :xy))
           (theta (atan y x))
           )
      (+ (* sx (expt (sin theta) 2))
         (* sy (expt (cos theta) 2))
         (* -2 sxy (sin theta) (cos theta))))))

(defun compute-radial-stress-normal (mp)
  (let ((pos (cl-mpm/particle::mp-position mp))
        (stress (cl-mpm/particle::mp-stress mp)))
    (let* ((n (cl-mpm/fastmaths:norm
               (cl-mpm/fastmaths:fast-.-
                pos
                *origin*))))
      (magicl:@ (magicl:transpose n) (cl-mpm/utils:voigt-to-matrix stress) n))))
(defun compute-radial-stress-analytic (mp)
  (let ((radius (cl-mpm/fastmaths:mag
                 (cl-mpm/fastmaths:fast-.-
                  (cl-mpm/particle:mp-position mp)
                  *origin*))))
    (let ((b 5d0)
          (a 1d0)
          (pressure 100e3)
          )
      (* (/ pressure (- (expt (/ b a) 2) 1))
         (- 1 (expt (/ b radius) 2))))))

(defun compute-hoop-stress-normal (mp)
  (let ((pos (cl-mpm/particle::mp-position mp))
        (stress (cl-mpm/particle::mp-stress mp)))
    (let* ((n (cl-mpm/setup::2d-orthog (cl-mpm/fastmaths:norm pos))))
      (pprint n)
      (magicl:@ (magicl:transpose n) (cl-mpm/utils:voigt-to-matrix stress) n))))


(defun compute-error (&optional sim)
  (let* ((mp-list (loop for mp across (cl-mpm:sim-mps *sim*) collect mp))
         (vp-list (loop for mp in mp-list collect (cl-mpm/particle::mp-volume mp)))
         (pressure 100e3)
         (syy (mapcar #'compute-radial-stress mp-list))
         (syy-ref (mapcar #'compute-radial-stress-analytic mp-list))
         (err (/
               (* (magicl:norm
                   (cl-mpm/fastmaths:fast-scale!
                    (magicl:.*
                     (magicl:from-list vp-list (list (length syy) 1))
                     (magicl:.-
                      (magicl:from-list syy-ref (list (length syy-ref) 1))
                      (magicl:from-list syy (list (length syy) 1)))
                     )
                    (/ 1d0 pressure))))
                 (abs (* (reduce #'+ vp-list))))))
    err))

(defun save-csv (output-file)
  (let* ((df (lisp-stat:make-df '(:error)
                                (list (make-array 1 :initial-element (compute-error))))))
    (lisp-stat:write-csv df output-file :add-first-row t))
  )


(defun loading-function (f)
  (let ((load (* f 100d3)))
    (setf (cl-mpm/buoyancy::bc-pressure-pressures *bc-pressure*) (list load load 0d0)))
  )

(defparameter *run-sim* nil)
(defun run-conv ()
  (setf *run-sim* t)
  (defparameter *data-refine* (list))
  (defparameter *data-error* (list))
  (loop for i in '(1 2 3 4 5)
        while *run-sim*
        do
           (let* ((refine i)
                  (mps 2))
             (let* ()
               (setup
                  :refine i
                  :mp-scale mps
                  :sym t)
               (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) nil)
               (cl-mpm/setup::set-mass-filter *sim* 1d3 :proportion 1d-9)
               (format t "Running sim size ~a ~a ~%" refine mps)
               (cl-mpm/dynamic-relaxation::run-load-control
                *sim*
                :output-dir (merge-pathnames (format nil "./output-~A_~D/" i mps))
                :load-steps 10
                :substeps (* 20 refine)
                :plotter #'plot
                :damping (sqrt 2)
                :save-vtk-dr t
                :save-vtk-loadstep t
                :dt-scale 1d0
                :criteria 1d-9
                :loading-function #'loading-function)
               ;; (plot-sigma-yy)
               ;; (push (compute-error *sim*) *data-error*)
               ;; (push h *data-refine*)
               )
             ;; (vgplot:loglog (mapcar (lambda (x) (/ 1d0 x)) *data-refine*) *data-error*)
             (save-csv (merge-pathnames (format nil "./analysis_scripts/virtual_stress/thick-cylinder/data/data-~A_~D.csv" i mps)))
             )))

;; (defun run-)

(defun stop ()
  (setf (cl-mpm::sim-run-sim *sim*) nil)
  (setf *run-sim* nil))

(defun plot-sigma-yy (&optional sim)
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (let* ((x-slice-pos (loop for mp across mps maximize (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
           (mp-list
             (loop for mp across mps
                     collect mp))
           (y (loop for mp in mp-list collect (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)))
           (y-ref (loop for pos in *original-configuration*
                        collect pos))

           (syy (loop for mp in mp-list collect (magicl:tref (cl-mpm/particle::mp-stress mp) 1 0)))
           (pressure -1d4)
           (max-y 50)
           (syy-ref (mapcar (lambda (x) pressure) y-ref))
           )
      (vgplot:plot syy y-ref "first;;with points pt 7"
                   syy-ref y-ref "Ref;;with points pt 7"
                   )
      (vgplot:legend)
      )))

;; (defun compute-error (&optional sim)
;;   (with-accessors ((mps cl-mpm:sim-mps))
;;       *sim*
;;     (let* ((x-slice-pos (loop for mp across mps maximize (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
;;            (mp-list
;;              (loop for mp across mps
;;                      collect mp))
;;            (y (loop for mp in mp-list collect (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)))
;;            (y-ref (loop for pos in *original-configuration*
;;                         collect pos))

;;            (vp-0-list (loop for size in *original-size*
;;                             collect (*  (cl-mpm/utils:varef size 0) (cl-mpm/utils:varef size 1))))
;;            (vl-0 (loop for vp-0 in vp-0-list sum vp-0))

;;            (syy (loop for mp in mp-list collect (magicl:tref (cl-mpm/particle::mp-stress mp) 1 0)))
;;            (pressure -1d4)
;;            (max-y 50)
;;            (syy-ref (mapcar (lambda (x) pressure) y-ref))
;;            )
;;       (loop for ref in syy-ref
;;             for val in syy
;;             for vp-0 in vp-0-list
;;             sum (/ (* (abs (- ref val)) vp-0) (* max-y vl-0))))))

(defun save-sigma-yy (&optional sim)
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (let* ((x-slice-pos (loop for mp across mps maximize (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
           (mp-list
             (loop for mp across mps
                   when (>= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (- x-slice-pos 0.001))
                     collect mp))
           (y (loop for mp in mp-list collect (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)))
           (syy (loop for mp in mp-list collect (magicl:tref (cl-mpm/particle::mp-stress mp) 1 0))))
      (with-open-file (stream (merge-pathnames "output/consolidation.csv") :direction :output :if-exists :supersede)
        (format stream "coord_y,sigma_yy~%")
        (loop for ymp in y
              for syymp in syy
              do (format stream "~f, ~f ~%" ymp syymp))))))


(defun test-dr ()
  (setup :mps 4 :refine 2 :sym t)
  (setf (cl-mpm/aggregate::sim-enable-aggregate *sim*) t)
  (cl-mpm/setup::set-mass-filter *sim* 1d3 :proportion 1d-9)
  ;; (cl-mpm:update-sim *sim*)
  ;; (plot *sim*)
  (cl-mpm/dynamic-relaxation::run-load-control
   *sim*
   :output-dir "./output/"
   :load-steps 1
   :plotter #'plot
   :damping 1d0
   :kinetic-damping nil
   :save-vtk-dr t
   :save-vtk-loadstep t
   :substeps 10
   :dt-scale 1d0
   :criteria 1d-9
   :loading-function #'loading-function))
