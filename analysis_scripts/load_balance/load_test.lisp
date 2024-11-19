#|
(progn
(ql:quickload :magicl)
(ql:quickload :vgplot))
(ql:quickload :lparallel)
(ql:quickload :swank.live)
|#


(defun vector-from-list (p)
  (magicl:from-list p '(2 1)))

(defun vector-zero ()
  (vector-from-list (list 0d0 0d0)))
(defun varef (v i)
  (aref (magicl::matrix/double-float-storage v) i))

(defclass particle ()
  ((position
    :accessor particle-position
    :initarg :pos
    :initform (vector-zero))
   (index
    :accessor particle-index
    :initform 0
    :initarg :index
    )))

(defclass domain ()
  ((index
    :accessor domain-index
    :initarg :index
    :initform -1)
   (pos
    :accessor domain-position
    :initarg :position
    :initform (list 0 0))
   (particles
    :accessor domain-particles
    :initarg :particles
    :initform (list))
   (particle-count
    :accessor domain-particle-count
    :initform 0
    )
   (bounds
    :accessor domain-bounds
    :initarg :bounds
    :initform (list (list 0d0 1d0)
                    (list 0d0 1d0)))))

(defclass sim ()
  ((domain-count
    :accessor sim-domain-count
    :initform (list 0d0 0d0)
    )
   (domains
    :accessor sim-domains
    :initform (list)
    )
   
   (particles
    :accessor sim-particles
    :initform (list)
    )))

;; (defparameter *sim* nil)
(defun setup ()
  (defparameter *sim* (make-instance 'sim))
  (setf (sim-particles *sim*)
        (loop repeat 1000
              collect (make-instance 'particle
                                     :pos (vector-from-list
                                           (list (* (random 1d0) (random 1d0))
                                                 (sin (random 1d0))
                                                 )))))
  (let* ((dc (list 4 4))
         (size (mapcar (lambda (x) (/ 1d0 x)) dc))
         (index 0)
         )
    (setf (sim-domain-count *sim*) dc)
    (setf (sim-domains *sim*)
          (loop for x from 0 below (first dc)
                append
                (loop for y from 0 below (second dc)
                      collect (prog1
                                  (make-instance 'domain
                                                 :index index
                                                 :position (list x y)
                                                 :bounds (list (list (* x (first size)) (* (1+ x) (first size)))
                                                               (list (* y (second size)) (* (1+ y) (second size)))))
                                (incf index))))))
  )

(defun get-domain (sim x y)
  (let ((out nil))
    (loop for d in (sim-domains sim)
          while (not out)
          when (and
                (= x (first (domain-position d)))
                (= y (second (domain-position d))))
            do (setf out d)
          )
    out))


(defun plot ()

  (vgplot:format-plot t "unset object")
  (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  (multiple-value-bind (x y c)
      (loop for p in (sim-particles *sim*)
            collect (varef (particle-position p) 0) into x
            collect (varef (particle-position p) 1) into y
            collect 0d0 into c
            finally (return (values x y c)))
    (loop for d in (sim-domains *sim*)
          do
             (destructuring-bind ((xl xu) (yl yu)) (domain-bounds d)
               (vgplot:format-plot t "set object ~D rect from ~f,~f to ~f,~f fc transparent fs rgb 'black' solid 0.8 noborder behind"
                                   (1+ (domain-index d))
                                   xl yl
                                   xu yu
                                   ))
          )
    (vgplot:plot x y c ";;lc palette"))
  (vgplot:format-plot t "set xrange [~f:~f]" 0d0 1d0)
  (vgplot:format-plot t "set yrange [~f:~f]" 0d0 1d0)
  (vgplot:format-plot t "set size ratio ~f" 1d0)
  )


(defun in-bounds (domain pos)
  (with-accessors ((bounds domain-bounds))
      domain
    (let ((in t))
      (loop for d from 0 to 1
            while in
            do
               (destructuring-bind (bl bu) (nth d bounds)
                 (setf in (and in
                               (<= bl (nth d pos))
                               (> bu (nth d pos))))))
      in)))
(defun in-bounds-particle (domain particle)
  (with-accessors ((bounds domain-bounds))
      domain
    (with-accessors ((pos particle-position))
        particle
      (in-bounds domain (list (varef pos 0) (varef pos 1))))))

;; (defun iterate-over-particles (sim func)
;;   (with-accessors ((particles sim-particles))
;;       sim
;;     (lparallel:pmapcar func particles)))
(defun count-particles (sim domain)
  (with-accessors ((particles sim-particles))
      sim
    (let ((iter 0))
      (loop for p in particles
            do (when (in-bounds-particle domain p)
                 (incf iter)))
      iter)))

(defun update-counts ()
  (lparallel:pmapcar
   (lambda (d)
     (setf (domain-particle-count d) (count-particles *sim* d)))
   (sim-domains *sim*)))

(defun load-balance ()
  (let ((step-size 1d-1))
    (update-counts)
    (let ((dc (sim-domain-count *sim*)))
      (loop for d in (list 0 1)
            do
               (let* ((arr (make-array (nth d dc)))
                      (size (make-array (nth d dc)))
                      )
                 (loop for dom in (sim-domains *sim*)
                       do
                          (progn
                            (incf (aref arr (nth d (domain-position dom)))
                                  (domain-particle-count dom))
                            (setf (aref size (nth d (domain-position dom)))
                                  (- (apply #'- (nth d (domain-bounds dom))))))
                       )
                 (let ((incr (make-array (- (nth d dc) 1))))
                   (loop for i from 0 below (- (nth d dc) 1)
                         do (setf (aref incr i)
                                  (* step-size
                                     (min (aref size (1+ i)) (aref size i))
                                     (if (or (> (aref arr (1+ i)) 0)
                                             (> (aref arr i) 0))
                                         (/
                                          (abs (- (aref arr (1+ i)) (aref arr i)))
                                          (max (aref arr (1+ i)) (aref arr i)))
                                         0d0
                                         )
                                     (if (> (aref arr (1+ i)) (aref arr i))
                                         -1
                                         1
                                         ))))
                   (format t "Mps ~A~%" arr)
                   (pprint size)
                   (pprint incr)
                   (loop for dom in (sim-domains *sim*)
                         do
                            (let* ((v (nth d (domain-position dom))))
                              (when (>= (- v 1) 0)
                                (incf (first (nth d (domain-bounds dom)))
                                      (- (aref incr (- v 1)))))
                              (when (< v (- (nth d dc) 1))
                                (incf (second (nth d (domain-bounds dom)))
                                      (- (aref incr v))))
                              )
                         )
                   )
                 )))))


(defun test ()
  (setup)
  (dotimes (i 50)
    (plot)
    (dotimes (k 10)
      (load-balance))
    (update-counts)
    (let ((pmin (loop for d in (sim-domains *sim*) minimize (domain-particle-count d)))
          (pmax (loop for d in (sim-domains *sim*) maximize (domain-particle-count d)))
          )
      (format t "Min: ~D~%" pmin)
      (format t "Max: ~D~%" pmax)
      (format t "Occ: ~F~%" (* (if (> pmin 0) (/ pmax pmin) 0d0) 100d0)))
    (sleep 0.1d0)
    (swank.live:update-swank)))
