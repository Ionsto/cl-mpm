(ql:quickload :vgplot)
(ql:quickload :swank.live)
(ql:quickload :lparallel)

(defclass node ()
  ((eta
    :accessor node-eta
    :initform 0d0
    )
   (eta-dash-prev
    :accessor node-eta-dash-prev
    :initform 0d0)
   (eta-dash
    :accessor node-eta-dash
    :initform 0d0)
   (pos
    :accessor node-pos
    :initarg :pos)
   (bc
    :accessor node-bc
    :initform nil)
   (alpha
    :accessor node-alpha
    :initform 1d0
    :initarg :alpha)
   (eta-dt
    :accessor node-eta-dt
    :initform 0d0)
   (eta-dt-dt
    :accessor node-eta-dt-dt
    :initform 0d0)
   ))
(defclass sim ()
  ((nodes
    :accessor sim-nodes
    :initarg :nodes
    )
   (dt
    :accessor sim-dt
    :initarg :dt)
   (h
    :accessor sim-h
    :initarg :h)
   ))
(defun get-eta-dash (nodes)
  (loop for node in nodes collect (node-eta-dash node)))
(defun plot-data (nodes)
  (let ((x (loop for node across nodes collect (node-pos node)))
        (eta-dash (loop for node across nodes collect (node-eta-dash node)))
        )
    (vgplot:plot x eta-dash)))
(defun iterate-over-nodes (nodes func)
  (declare (function func))
  (lparallel:pdotimes (i (length nodes))
    (funcall func (aref nodes i) i)))
(defun iterate-over-domain-nodes (nodes func)
  (declare (function func))
  (lparallel:pdotimes (i (- (length nodes) 2))
    (funcall func (aref nodes (1+ i)) (1+ i))))
(defparameter *run-sim* t)
(defun stop ()
  (setf *run-sim* nil))
(defun update (nodes h dt)
  (iterate-over-domain-nodes
   nodes
   (lambda (node x)
     (let ((node-right (aref nodes (+ x 1)))
           (node-left (aref nodes (- x 1))))
       (setf (node-eta-dt node)
             (-
              (- (node-eta-dash node) (node-eta node))
              (* (node-alpha node)
                 (/ (- (+ (node-eta-dash node-left)
                          (node-eta-dash node-right))
                       (* 2 (node-eta-dash node))) (expt h 2))))))

     )))

(defun integrate (nodes h dt)
  (iterate-over-nodes
   nodes
   (lambda (node x)
     (let ((bc (node-bc node)))
       (when bc
         (setf (node-eta-dash node) bc
               (node-eta-dash-prev node) bc
               (node-eta-dt node) 0d0
               (node-eta node) bc)))))
  (iterate-over-domain-nodes
   nodes
   (lambda (node i)
     (setf (node-eta-dash-prev node)
           (node-eta-dash node))
     (incf (node-eta-dash node)
           (* -1d0 dt (node-eta-dt node))))))


(defun criteria (sim)
  (let ((f-sum 0d0)
        (eta-sum 0d0))
    (loop for node across (sim-nodes sim)
          do (progn
               (incf f-sum (abs (- (node-eta-dash node)
                                   (node-eta-dash-prev node))))
               (incf eta-sum (node-eta node))
               ))
    ;; f-sum
    (abs (/ f-sum eta-sum))
    ))
(defun test ()
  (defparameter *run-sim* t)
  (let* ((node-count 1000)
         (domain-size 100d0)
         (steps 100)
         (h (/ domain-size (- node-count 1)))
         (l 50d0)
         (alpha l;; (/ 1 (expt l 2))
                )
         ;(dt (/ (expt h 2) domain-size))
         (dt (* 0.4d0 (/ (expt h 2) alpha)))
         (f nil)
         (nodes (make-array node-count :initial-contents (loop for x from 0 below node-count collect (make-instance 'node
                                                                                                                    :pos (* x h)
                                                                                                                    :alpha alpha)))))
    (defparameter *sim* (make-instance 'sim
                                       :nodes nodes
                                       :dt dt
                                       :h h))
    (setf (node-bc (aref nodes 0)) 1d0)
    ;; (setf (node-bc (aref nodes 40)) 0.8d0)
    ;; (setf (node-bc (aref nodes 40)) 0.8d0)
    (loop for x from 0 to steps
          while (and *run-sim*
                     (or
                      ;; (< x 2)
                      (if f (> f 1d-5)
                          t)))
          do
             (progn
               (format t "Step ~D~%" x)
               (plot-data nodes)
               (dotimes (i 100)
                 (update nodes h dt)
                 (integrate nodes h dt))
               (setf f (criteria *sim*))
               (format t "Criteria ~E~%" f)
               (sleep 0.1)
               (swank.live:update-swank)
               ))))
