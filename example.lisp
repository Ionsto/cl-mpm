(ql:quickload "cl-mpm")
(ql:quickload "vgplot")
(ql:quickload "swank.live")
(ql:quickload "magicl")


(defparameter *nD* 2)

(defun dostep (sim n)
  (progn
    (dotimes (i n)
      (cl-mpm::update-sim sim)
      )
      (when  (cl-mpm:sim-mps sim)
      (let* ( 
             (pos (loop for mp in (cl-mpm:sim-mps sim) collect (slot-value mp 'cl-mpm::position)))
             (sig (loop for mp in (cl-mpm:sim-mps sim) collect (slot-value mp 'cl-mpm::stress)))
             (x (loop for p in pos collect (magicl:tref p 0 0)))
             (sigy (loop for s in sig collect (magicl:tref s 1 0)))
             (y (loop for p in pos collect (magicl:tref p 1 0))))
        (vgplot:plot x y ";;with points pt 7")
        ))
      (vgplot:replot)
      (swank.live:update-swank)
      (sleep .05)
      ))
(defun make-column (size)
  (let* ((nD 2)
         (mp-spacing 1d0)
         (sim (cl-mpm:make-mpm-sim (list 1 size) 1 1e-3 
                                                  (cl-mpm::make-shape-function-linear nD))))
    (progn (setf (cl-mpm:sim-mps sim) 
                 (loop for i from 0 to (- size 1) collect 
                       (cl-mpm::make-particle-elastic nD
                                                      1e2
                                                      0
                                                      :pos (list 0.5 (+ 0.5d0 (* mp-spacing i)))
                                                      :volume mp-spacing)))
          (setf (cl-mpm:sim-bcs sim) (cl-mpm/bc:make-outside-bc (cl-mpm:mesh-mesh-size (cl-mpm:sim-mesh sim)))) 
           sim)))
(progn 
  (defparameter *mp-spacing* 1d0)
  (defparameter *sim* (cl-mpm:make-mpm-sim '(1 32) 1 1e-3 
                                           (cl-mpm::make-shape-function-linear *nD*)))
  (setf (cl-mpm:sim-mps *sim*) 
        (loop for i from 0 to 31 collect 
              (cl-mpm::make-particle-elastic *nD*
                                             1e2
                                             0
                                             :pos (list 0.5 (+ 0.5d0 (* *mp-spacing* i)))
                                             :volume *mp-spacing*)))
  (setf (cl-mpm:sim-bcs *sim*) (cl-mpm/bc:make-outside-bc (cl-mpm:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))))

(defmacro time-form-noprint (form it)
  `(progn
    (declaim (optimize speed))
    (let* ((iterations ,it)
           (start (get-internal-real-time)))
      (dotimes (i iterations) ,form)
      (let* ((end (get-internal-real-time))
             (units internal-time-units-per-second)
             )
        (/ (- end start) (* iterations units))))))
(defmacro time-form (form it)
  `(progn
    (declaim (optimize speed))
    (let* ((iterations ,it)
           (start (get-internal-real-time)))
      (dotimes (i iterations) ,form)
      (let* ((end (get-internal-real-time))
             (units internal-time-units-per-second)
             )
        (format t "Total time: ~f ~%" (/ (- end start) units)) (format t "Time per iteration: ~f~%" (/ (- end start) (* iterations units)))
        (format t "Throughput: ~f~%" (/ 1 (/ (- end start) (* iterations units))))))))


(let ((mesh (cl-mpm:sim-mesh *sim*))
      (mps (cl-mpm:sim-mps *sim*))
      )
  (time-form (cl-mpm::update-sim *sim*) 1000))

(vgplot:figure)
(let* ((sizes (loop for i from 1 to 8 collect (expt 2 i)))
       (times (loop for size in sizes collect 
                    (let ((sim (make-column size)))
                      (time-form-noprint (cl-mpm::update-sim sim) 100))))
       (throughput (loop for dt in times collect (/ 1 dt)))
       )
  (vgplot:xlabel "Size")
  (vgplot:ylabel "Throughput")
  (vgplot:plot sizes throughput)
  (vgplot:figure)
  (vgplot:loglog sizes times)
  (defparameter *perf-sizes* sizes)
  (defparameter *perf-times* times)
  (defparameter *perf-throughput* throughput)
  )
(let ((s  (mapcar (lambda (x) (log x)) *perf-sizes*))
      (tm (mapcar (lambda (x) (log x)) *perf-times*))
      (a 5)
      (b 7))
  (/ (- (nth a tm) (nth b tm)) (-(nth a s) (nth b s))))

(/ 1d0 (time-form-noprint (cl-mpm::update-sim *sim*) 100))

(let ((mesh (cl-mpm:sim-mesh *sim*))
      (mps (cl-mpm:sim-mps *sim*))
      (bcs (cl-mpm:sim-bcs *sim*))
      (dt (cl-mpm:sim-dt *sim*))
      )
  (format t "Reset-grid ~%")
  (time-form (cl-mpm::reset-grid mesh) 10000) 
  (format t "P2G ~%")
  (time-form (cl-mpm::p2g mesh mps) 1000) 
  (format t "Filter ~%")
  (time-form (cl-mpm::filter-grid mesh 1e-3) 10000) 
  (format t "Update-nodes ~%")
  (time-form (cl-mpm::update-nodes mesh dt) 1000) 
  (format t "Apply-bcs ~%")
  (time-form (cl-mpm::apply-bcs mesh bcs) 1000) 
  (format t "G2P ~%")
  (time-form (cl-mpm::p2g mesh mps) 1000) 
  (format t "Update-particles ~%")
  (time-form (cl-mpm::update-particle mps dt) 1000)
  (format t "Update-stress ~%")
  (time-form (cl-mpm::update-stress mesh mps dt) 1000)
  )

(let ((mesh (cl-mpm:sim-mesh *sim*))
      (mps (cl-mpm:sim-mps *sim*))
      (bcs (cl-mpm:sim-bcs *sim*))
      (dt (cl-mpm:sim-dt *sim*))
      )
  (format t "Update-nodes ~%")
  (time-form (cl-mpm::update-nodes mesh dt) 10000) 
  (format t "Update-nodes-serial ~%")
  (time-form (cl-mpm::update-nodes-serial mesh dt) 10000) 
  (format t "Update-stress ~%")
  (time-form (cl-mpm::update-stress mesh mps dt) 1000)
  (format t "Update-stress-serial ~%")
  (time-form (cl-mpm::update-stress-serial mesh mps dt) 1000)
  )

(setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))
(defparameter *sim* (make-column 4))
(time-form (cl-mpm::update-sim *sim*) 100)

(time (dotimes (i 100) (cl-mpm::update-sim *sim*)))
(/ 1 (/ 12.850 10000))


(vgplot:figure)
(vgplot:axis (list 0 (nth 0 (cl-mpm:mesh-mesh-size (cl-mpm:sim-mesh *sim*))) 
                   0 (nth 1 (cl-mpm:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))))

(loop for mp in (cl-mpm:sim-mps *sim*) do
      (setf (cl-mpm::mp-E mp) 1e1))
(dotimes (k 100)
  (format t "Frame ~d ~%" k)
  (dostep *sim* 10))



(defun const-newt (density mcst)
  (destructuring-bind (rest_density eos_stiffness eos_power) mcst
    (* eos_stiffness (expt (- (/ density rest_density) 1) eos_power))))
(defun const-ideal (density mcst)
  (destructuring-bind (R temp) mcst
    (* density R temp)))

(defparameter *density* (loop for i from 1 to 1.4 by 0.01 collect i))
(defparameter *mcst-newt* '(1.22 9.5288e04 1.4))
(defparameter *mcst-ideal* '(287.05 300))
(defparameter *pressure-newt* (mapcar (lambda (density) (const-newt density *mcst-newt*)) *density*))
(defparameter *pressure-ideal* (mapcar (lambda (density) (const-ideal density *mcst-ideal*)) *density*))
*pressure-ideal*
*pressure-newt*
(vgplot:figure)
(vgplot:plot *density* *pressure-newt*
             *density* *pressure-ideal*) 
