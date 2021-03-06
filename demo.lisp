(ql:quickload "cl-mpm")
(ql:quickload "cl-mpm/setup")
(ql:quickload "vgplot")
(ql:quickload "swank.live")

(ql:quickload "magicl")

(declaim (optimize (speed 3)))

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
;(defun make-column (size)
;  (let* ((nD 2)
;         (mp-spacing 1d0)
;         (sim (cl-mpm:make-mpm-sim (list 1 size) 1 1e-3 
;                                                  (cl-mpm::make-shape-function-linear nD))))
;    (progn (setf (cl-mpm:sim-mps sim) 
;                 (loop for i from 0 to (- size 1) collect 
;                       (cl-mpm::make-particle-elastic nD
;                                                      1e2
;                                                      0
;                                                      :pos (list 0.5 (+ 0.5d0 (* mp-spacing i)))
;                                                      :volume mp-spacing)))
;          (setf (cl-mpm:sim-bcs sim) (cl-mpm/bc:make-outside-bc (cl-mpm:mesh-size (cl-mpm:sim-mesh sim)))) 
;           sim)))


;(progn 
;  (defparameter *nD* 2)
;  (defparameter *mp-spacing* 1d0)
;  (defparameter *sim* (cl-mpm:make-mpm-sim '(8 8) 1 1e-3 
;                                           (cl-mpm::make-shape-function-linear *nD*)))
;  (setf (cl-mpm:sim-mps *sim*) 
;        (loop for x from 0 to 4 append 
;              (loop for y from 0 to 4 collect 
;                     (cl-mpm::make-particle-elastic *nD*
;                                                    1e2
;                                                    0
;                                                    :pos (list (+ 0.5d0 (* *mp-spacing* x))
;                                                               (+ 0.5d0 (* *mp-spacing* y)))
;                                                    :volume *mp-spacing*))))
;  (setf (cl-mpm:sim-bcs *sim*) (cl-mpm/bc:make-outside-bc (cl-mpm:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))))
;
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
             (dt (/ (- end start) (* iterations units)))
             )
        (format t "Total time: ~f ~%" (/ (- end start) units)) (format t "Time per iteration: ~f~%" (/ (- end start) (* iterations units)))
        (format t "Throughput: ~f~%" (/ 1 dt))
        dt))))


(defun test-throughput ()
  (progn
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
      ;(defparameter *perf-sizes* sizes)
      ;(defparameter *perf-times* times)
      ;(defparameter *perf-throughput* throughput)
      )))
;(let ((s  (mapcar (lambda (x) (log x)) *perf-sizes*))
;      (tm (mapcar (lambda (x) (log x)) *perf-times*))
;      (a 5)
;      (b 7))
;  (/ (- (nth a tm) (nth b tm)) (-(nth a s) (nth b s))))

(defun test-sub-steps (sim)
  (let ((mesh (cl-mpm:sim-mesh sim))
        (mps (cl-mpm:sim-mps sim))
        (bcs (cl-mpm:sim-bcs sim))
        (dt (cl-mpm:sim-dt sim))
        (reset-grid-time 0)
        (p2g-time 0)
        (filter-grid-time 0)
        (update-nodes-time 0)
        (apply-bcs-time 0)
        (g2p-time 0)
        (update-particle-time 0)
        (update-stress-time 0)
        )
    (format t "Reset-grid ~%")
    (setf reset-grid-time (time-form (cl-mpm::reset-grid mesh) 10000)) 
    (format t "P2G ~%")
    (setf p2g-time (time-form (cl-mpm::p2g mesh mps) 1000)) 
    (format t "Filter ~%")
    (setf filter-grid-time  (time-form (cl-mpm::filter-grid mesh 1e-3) 10000)) 
    (format t "Update-nodes ~%")
    (setf update-nodes-time (time-form (cl-mpm::update-nodes mesh dt) 1000)) 
    (format t "Apply-bcs ~%")
    (setf apply-bcs-time (time-form (cl-mpm::apply-bcs mesh bcs) 1000)) 
    (format t "G2P ~%")
    (setf g2p-time (time-form (cl-mpm::g2p mesh mps) 1000)) 
    (format t "Update-particles ~%")
    (setf update-particle-time (time-form (cl-mpm::update-particle mps dt) 1000))
    (format t "Update-stress ~%")
    (setf update-stress-time (time-form (cl-mpm::update-stress mesh mps dt) 1000))
    (let* ((time-list (list reset-grid-time p2g-time filter-grid-time update-nodes-time apply-bcs-time g2p-time update-particle-time update-stress-time))
           (total-time (/ (reduce #'+ time-list) 100)))
      (format t "Percentile usage:~%")
      (format t "Reset-grid ~f ~%" (/ reset-grid-time total-time))
      (format t "p2g ~f ~%" (/ p2g-time total-time))
      (format t "filter-grid ~f ~%" (/ filter-grid-time total-time))
      (format t "update-nodes ~f ~%" (/ update-nodes-time total-time))
      (format t "apply-bcs ~f ~%" (/ apply-bcs-time total-time))
      (format t "g2p ~f ~%" (/ g2p-time total-time))
      (format t "update-particle ~f ~%" (/ update-particle-time total-time))
      (format t "update-stress ~f ~%" (/ update-stress-time total-time))
      )))


(setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))
(progn
    (defparameter *sim* (cl-mpm/setup::make-column 1 64))
    (setf (cl-mpm:sim-mps *sim*) 
           (cl-mpm/setup::make-column-mps-elastic
             63 
             (list 1 63)
             1e4 0d0)))
;(test-sub-steps *sim*)
(let* ( (mesh (cl-mpm:sim-mesh *sim*))
        (mps (cl-mpm:sim-mps *sim*)))
  (time-form (cl-mpm::p2g mesh mps) 5000)
  (time-form (cl-mpm::p2g-array mesh mps) 5000)
  )


;(let ((mesh (cl-mpm:sim-mesh *sim*))
;      (mps (cl-mpm:sim-mps *sim*))
;      (bcs (cl-mpm:sim-bcs *sim*))
;      (dt (cl-mpm:sim-dt *sim*))
;      )
;  (format t "p2g ~%")
;  (time-form (cl-mpm::p2g mesh mps) 10000) 
;  (format t "p2g-serial ~%")
;  (time-form (cl-mpm::p2g-serial mesh mps) 10000) 
;  ;(format t "Update-nodes ~%")
;  ;(time-form (cl-mpm::update-nodes mesh dt) 10000) 
;  ;(format t "Update-nodes-serial ~%")
;  ;(time-form (cl-mpm::update-nodes-serial mesh dt) 10000) 
;  ;(format t "Update-stress ~%")
;  ;(time-form (cl-mpm::update-stress mesh mps dt) 1000)
;  ;(format t "Update-stress-serial ~%")
;  ;(time-form (cl-mpm::update-stress-serial mesh mps dt) 1000)
;  )


;(setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))
;(defparameter *sim* (make-column 16))
;(time-form (cl-mpm::update-sim *sim*) 100)
;(time (dotimes (i 100) (cl-mpm::update-sim *sim*)))
;
;
;(vgplot:figure)
;(vgplot:axis (list 0 (nth 0 (cl-mpm:mesh-mesh-size (cl-mpm:sim-mesh *sim*))) 
;                   0 (nth 1 (cl-mpm:mesh-mesh-size (cl-mpm:sim-mesh *sim*)))))
;
;(loop for mp in (cl-mpm:sim-mps *sim*) do
;      (setf (cl-mpm::mp-E mp) 1e1))
;(dotimes (k 50)
;  (format t "Frame ~d ~%" k)
;  (dostep *sim* 100))
;
