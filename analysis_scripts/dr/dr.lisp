(ql:quickload :vgplot)
(ql:quickload :magicl)
(ql:quickload :lisp-stat)
(ql:quickload :swank.live)
(ql:quickload :lparallel)

(defun ddx (vec index h)
  (cond
    ((= index 0) 0d0)
    ((= (1+ index) (length vec)) 0d0)
    (t (/ (+ (aref vec (- index 1))
             (* -2 (aref vec index))
             (aref vec (+ index 1))
             )
          (expt h 2)))))
(defun stop ()
  (defparameter *run-sim* nil))

(defun iterate-over-vec (vec func)
  (lparallel:pdotimes (i (- (length vec) 1))
    (funcall func i)))

(defun solution (x length)
  (let* ((r (make-array (length x) :initial-element 0d0))
         ;; (max-x length)
         (du (/ 1d0 (- (length x) 1))))
    (pprint du)
    (loop
      for i from 0 below (length x)
      do (setf (aref r i) (- 1d0 (* i du))))
    r))

(defun solve-parabolic (u h c)
  (let ((dt (* 0.5d0 (/ (expt h 2) c)))
        (ddu (make-array (length u) :initial-element 0d0)))
    (dotimes (i 10)
      (setf (aref u 0) 1d0)
      (setf (aref u (- (length u) 1)) 0d0)
      (iterate-over-vec
       u
       (lambda (i) (setf (aref ddu i) (ddx u i h))))
      (iterate-over-vec
       u
       (lambda (i)
         (incf (aref u i) (* dt c (aref ddu i))))))))
(defun solve-hyperbolic (u v h c damp)
  (let ((dt (* 0.5d0 (/ (expt h 1) (sqrt c))))
        (ddu (make-array (length u) :initial-element 0d0)))
    (dotimes (i 10)
      (setf (aref u 0) 1d0)
      (setf (aref u (- (length u) 1)) 0d0)
      (setf (aref v 0) 0d0)
      (setf (aref v (- (length u) 1)) 0d0)
      (iterate-over-vec
       u
       (lambda (i) (incf (aref v i)
                         (* dt
                            (-
                             (* c (ddx u i h))
                             (* damp (aref v i)))))))
      (iterate-over-vec
       u
       (lambda (i)
         (incf (aref u i) (* dt (aref v i))))))))

(defun residual (u r)
  (/
   (reduce #'+ (aops:map-array (aops:each #'- u r) (lambda (x) (expt x 2)) ))
   (reduce #'+ (aops:each #'* r r))))
(defun test ()
  (let* ((length 10d0)
         (resolution 100)
         (h (/ length (- resolution 1)))
         (c 1d0)
         (damp 0.25d0)
         (steps 1000)
         ;; (dt (* 0.95d0 (/ (expt h 1) (sqrt c))))
         ;; (dt )
         (x (make-array resolution :initial-contents
                        (loop
                          repeat resolution
                          for x from 0d0 by h collect x)))
         (bcs (list 1 0))
         (u (make-array resolution :initial-element 0d0))
         (u-para (make-array resolution :initial-element 0d0))
         (v (make-array resolution :initial-element 0d0))
         (r (make-array steps :initial-element 0d0 :fill-pointer 0))
         (r-para (make-array steps :initial-element 0d0 :fill-pointer 0))
         (data-step (make-array steps :initial-element 0d0  :fill-pointer 0))
         (res (solution x length))
         )

    (defparameter *u* u)
    (defparameter *x* x)
    (defparameter *run-sim* t)
    ;; (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (loop for i from 0 below steps
          while *run-sim*
          do
             (progn
               (vector-push i data-step)
               (solve-parabolic u-para h c)
               (solve-hyperbolic u v h c damp)
               (vector-push (residual u res) r)
               (vector-push (residual u-para res) r-para)
               ;; (vgplot:plot
               ;;  (coerce x 'list) (coerce u 'list) "Parabolic"
               ;;  (coerce x 'list) (coerce u-para 'list) "Hyperbolic"
               ;;              )
               (vgplot:xlabel "Spacial position (x)")
               (vgplot:ylabel "Primary variable (u)")
               (vgplot:semilogy
                (coerce data-step 'list) (coerce r-para 'list) "Parabolic"
                (coerce data-step 'list) (coerce r 'list) "Hyperbolic")
               (vgplot:print-plot (merge-pathnames (format nil "outframes-conv/frame_~5,'0d.png" i))
                                  :terminal "png size 1920,1080")

               (format t "~E~%" (aref r i))
               (sleep 0.5)
               (swank.live:update-swank)))))
