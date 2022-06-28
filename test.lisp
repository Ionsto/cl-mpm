(ql:quickload "magicl")
(ql:quickload "lparallel")

(defstruct data
  a b)
(defun make-data-init ()
  (let ((d (make-data)))
    (setf (data-a d) (magicl:zeros '(2 2)))
    (setf (data-b d) (magicl:zeros '(2 2)))
    d))

(defun do-work (item)
  (magicl:@ (data-a item) (data-b item)))

(defmacro time-form-noprint (form it)
  `(progn
    ;(declaim (optimize speed))
    (let* ((iterations ,it)
           (start (get-internal-real-time)))
      (dotimes (i iterations) ,form)
      (let* ((end (get-internal-real-time))
             (units internal-time-units-per-second)
             )
        (/ (- end start) (* iterations units))))))

(defun loop-array (items)
  (lparallel:pdotimes (i (length items)) 
    (do-work (aref items i))))

(defun loop-list (items)
  (lparallel:pdotimes (i (length items)) 
    (do-work (nth i items))))

(progn 
  (defparameter *item-count* 10000)
  (defparameter *item-list* (loop for i from 1 to *item-count* collect (make-data-init)))
  (defparameter *item-array* (make-array *item-count* :initial-contents *item-list*)))
(setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))
(time (loop-array *item-array*))
(time (loop-list *item-list*))
