(defpackage :cl-mpm/shape-function
  (:use :cl)
  (:export
   #:order
   #:svp
   #:dsvp
   #:shape-function
   #:shape-function-linear
   #:shape-function-bspline
   #:shape-function-gimp
   #:make-shape-function-linear
   #:make-shape-function-bspline
   #:assemble-dsvp
   #:assemble-dsvp-2d
   #:shape-bspline
   #:shape-bspline-dsvp
   ))
(in-package :cl-mpm/shape-function)
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(defmacro shape-linear-form (x)
  `(quote (- 1d0 (abs ,x))))
;; (defun shape-linear (x h)
;;   (- 1d0 (abs (/ x h))))

(declaim (inline shape-linear)
         (ftype (function (double-float double-float) double-float) shape-linear))
(defun shape-linear (x h)
  (declare (type double-float x h))
  (the double-float (- 1d0 (abs (/ x h)))))

(declaim (inline shape-linear-dsvp)
         (ftype (function (double-float double-float) double-float) shape-linear-dsvp))
(defun shape-linear-dsvp (x h)
  (declare (type double-float x h))
  (the double-float (/ (if (> x 0d0)
                               -1d0
                               1d0) h)
  ;; (the double-float (/ (signum x) h)
       ))

(declaim (inline shape-gimp)
         (ftype (function (double-float double-float double-float) double-float) shape-gimp))
(defun shape-gimp (x l h)
  (declare (type double-float x l h))
  (cond
    ((and (< (- (+ h l)) x) (<= x (- l h)))
     (/ (expt (+ x h l) 2d0) (* 4d0 h l)))
    ((and (< (- l h) x) (<= x (- l)))
     (+ 1d0 (/ x h)))
    ((and (< (- l) x) (<= x l))
     (- 1d0 (/ (+ (* x x) (* l l)) (* 2d0 h l))))
    ((and (< l x) (<= x (- h l)))
     (- 1d0 (/ x h)))
    ((and (< (- h l) x) (<= x (+ h l)))
     (/ (expt (- (+ h l) x) 2d0) (* 4d0 h l)))
    (t
     0d0)))
(declaim (inline shape-gimp-dsvp)
         (ftype (function (double-float double-float double-float) double-float) shape-gimp-dsvp))
(defun shape-gimp-dsvp (x l h)
  (declare (type double-float x l h))
  (cond
    ((and (< (- (+ h l)) x) (<= x (- l h)))
     (/ (+ x h l) (* 2d0 h l))
     )
    ((and (< (- l h) x) (<= x (- l)))
     (/ 1d0 h))
    ((and (< (- l) x) (<= x l))
     (- (/ x (* h l)))
     )
    ((and (< l x) (<= x (- h l)))
     (- (/ 1d0 h)))
    ((and (< (- h l) x) (<= x (+ h l)))
     (- (/ (- (+ h l) x) (* 2d0 h l))))
    (t
     0d0)))




(declaim (inline shape-bspline)
         (ftype (function (double-float double-float) double-float) shape-bspline))
(defun shape-bspline (x h)
  (if (< (abs x) (/ h 2))
          (- (/ 3 4) (expt (/ (abs x) h) 2))
          (* (/ 1 8) (expt (- 3 (/ (* 2 (abs x)) h)) 2))))
(declaim (inline shape-bspline-dsvp)
         (ftype (function (double-float double-float) double-float) shape-bspline-dsvp))
(defun shape-bspline-dsvp (X H)
  (IF (< (ABS X) (/ H 2d0))
    (-
     (* (* 2d0 (/ (ABS X) H))
        (/
         (*
          (IF (> X 0d0)
              1d0
              -1d0)
          H)
         (EXPT H 2d0))))
    (* 1/8
       (* (* 2d0 (- 3d0 (/ (* 2d0 (ABS X)) H)))
          (-
           (/
            (*
             (* 2d0
                (IF (> X 0d0)
                    1d0
                    -1d0))
             H)
            (EXPT H 2d0)))))))

(defun bspline (knots eta index poly)
  (if (= poly 0)
      (if (and (>= eta (nth index knots)) (< eta (nth (+ index 1) knots)))
          1d0
          0d0)
      (+
       (if (/= (nth index knots) (nth (+ poly index) knots))
           (* (/ (- eta (nth index knots))
                 (- (nth (+ poly index) knots) (nth index knots)))
              (bspline knots eta index (- poly 1)))
           0)
       (if (/= (nth (+ 1 index) knots) (nth (+ 1 poly index) knots))
           (* (/ (- (nth (+ index poly 1) knots) eta)
                 (- (nth (+ poly index 1) knots) (nth (+ 1 index) knots)))
              (bspline knots eta (+ 1 index) (- poly 1)))
           0)
       )))
(defun bspline-dsvp (knots eta index poly)
  (let ((dsvp-knots knots
                   ; (cdr (butlast knots))
                    )
        (index (- index 1)))
    (-
     (let ((kpi (nth (+ poly index) knots))
           (ki (nth index knots)))
       (if (/= kpi ki)
           (* (/ poly (- kpi ki))
              (bspline dsvp-knots eta index (- poly 1)))
           0))
     (let ((kpi (nth (+ poly index 1) knots))
           (ki (nth (+ 1 index) knots)))
       (if (/= kpi ki)
           (* (/ poly (- kpi ki))
              (bspline dsvp-knots eta (+ 1 index) (- poly 1)))
           0))
     )))
;; (defun bspline-dsvp (knots eta index poly)
;;   (let ((dsvp-knots (cdr (butlast knots)))
;;         (index-t (- index 1)))
;;      (let ((kpi (nth (+ poly index 1) knots))
;;            (ki (nth (+ 1 index) knots)))
;;        (if (/= kpi ki)
;;            (* (/ poly (- kpi ki))
;;               (bspline dsvp-knots eta (+ 0 index) (- poly 1)))
;;            0)
;;      )))
(defun make-half-knots (nodes inc)
  (loop for n in nodes
        for ni in (cdr nodes)
        for pos = 0 then (incf pos inc)
        with i = 0
        collect (progn
                  (cond
                    ((and n ni)
                     (setf i (/ (+ pos pos inc) 2))
                     )
                    ((and (not n) (not ni))
                     i
                     )
                    ((and n (not ni))
                     (setf i pos)
                     )
                    ;; ((and (not n) ni)
                    ;;  (setf i (+ pos inc))
                    ;;  )
                    (t
                     i)
                    )
                  ;; (if (and n ni)
                  ;;         (setf i (/ (+ pos pos 1) 2))
                  ;;         (if ni
                  ;;             i
                  ;;             (setf i (+ pos 0))))
                  )))
(defun make-bspline-knots (nodes h)
  (let* ((mid-node-id (round (- (length nodes) 1) 2)))
    (nreverse (append
                     (nreverse (make-half-knots (subseq nodes mid-node-id) h))
                     (make-half-knots (nreverse (subseq nodes 0 (+ 1 mid-node-id))) (- h))
                     ))))
(defun nodal-bspline (nodes eta node h)
  ;; (setf (nth 4 nodes) (and (nth 3 nodes) (nth 5 nodes)))
  ;; (unless (nth 3 nodes)
  ;;   (incf node 0))
  ;; (unless (nth 5 nodes)
  ;;   (incf node 0))
  ;; (cond
  ;;   ((reduce #'eq (mapcar #'eq '(nil nil nil t
  ;;                                     t
  ;;                                     nil nil nil nil) nodes))
  ;;     ;; (print "Left linear")
  ;;     (shape-linear eta h)
  ;;     )
  ;;   ((reduce #'eq (mapcar #'eq '(nil nil nil nil
  ;;                                t
  ;;                                t nil nil nil) nodes))
  ;;    ;; (print "Right linear")
  ;;    (shape-linear eta h)
  ;;    )
  ;;   (t
  ;;    (setf (nth 4 nodes) (and (nth 3 nodes) (nth 5 nodes)))
  ;;    (bspline (make-bspline-knots nodes h) (+ eta) (+ node 2) 2)
  ;;    ))
  (setf (nth 4 nodes) (and (nth 3 nodes) (nth 5 nodes)))
  (bspline (make-bspline-knots nodes h) (+ eta) (+ node 2) 2)
  )
(defun nodal-bspline-dsvp (nodes eta node h)
  ;; (setf (nth 4 nodes) (and (nth 3 nodes) (nth 5 nodes)))
  ;; (unless (nth 3 nodes)
  ;;   (incf node -1))
  ;; (unless (nth 5 nodes)
  ;;   (incf node 1))
  (bspline-dsvp (make-bspline-knots nodes h) (+ eta) (+ node 2) 2))

;; (let ((x (loop for x from -2d0 upto 2d0 by 0.01d0 collect x))
;;       (h 1))
;;   (vgplot:figure)
;;   (vgplot:plot
;;    x
;;    (mapcar (lambda (y) (shape-gimp (+ y h) 0.01d0 1d0)) x) "-1"
;;    x
;;    (mapcar (lambda (y) (shape-gimp (+ y 0) 0.01d0 1d0)) x) "0"
;;    x
;;    (mapcar (lambda (y) (shape-gimp (- y h) 0.01d0 1d0)) x) "1"
;;    ))


(defmacro create-svp (arg form)
  `(lambda (,arg) ,form))

(defmacro create-dsvp (arg form)
  (let ((dx (symbolic-derivation:derive arg form)))
    `(lambda (,arg) ,dx)))


(defun svp-1d (svp dsvp)
  svp)
(defun svp-2d (svp dsvp)
  (lambda (x y) (* (funcall svp x) (funcall svp y))))
(defun svp-3d (svp dsvp)
  (lambda (x y z) (* (funcall svp x) (funcall svp y) (funcall svp y))))

(defun dsvp-1d (svp dsvp)
  (lambda (x) (list (funcall dsvp x))))
(defmacro execute-dsvp-2d (svp dsvp x y)
  `(let (
          (wx (,svp ,x))
          (wy (,svp ,y))
          (dx (,dsvp ,x))
          (dy (,dsvp ,y))
          )
                                        ;(d/dx d/dy)
      (list (* wy dx) (* wx dy))))
(defun dsvp-2d (svp dsvp)
  (lambda (x y) (let (
                      (wx (funcall svp x))
                      (wy (funcall svp y))
                      (dx (funcall dsvp x))
                      (dy (funcall dsvp y))
                      )
                  ;(d/dx d/dy)
                  (list (* wy dx) (* wx dy)))))
(defun dsvp-3d (svp dsvp)
  (lambda (x y z) 
    (let (
          (wx (funcall svp x))
          (wy (funcall svp y))
          (wz (funcall svp z))
          (dx (funcall dsvp x))
          (dy (funcall dsvp y))
          (dz (funcall dsvp z))
          )
      ;(d/dx d/dy d/dz)
      (list (* wy wz dx) (* wx wz dy) (* wx wy dz)))))

(defun nd-svp (nD svp dsvp)
  (case  nD
    (1 (svp-1d svp dsvp))
    (2 (svp-2d svp dsvp))
    (3 (svp-3d svp dsvp))))

(defun nd-dsvp (nD svp dsvp)
  (case  nD
    (1 (dsvp-1d svp dsvp))
    (2 (dsvp-2d svp dsvp))
    (3 (dsvp-3d svp dsvp))))

(defclass shape-function ()
  ((nD
     :accessor nD
     :initarg :nD)
   (order
     :accessor order
     :initarg :order)
   (svp
     :accessor svp
     :initarg :svp)
   (dsvp
     :accessor dsvp
     :initarg :dsvp)
   ))

(defclass shape-function-linear (shape-function)
  ((order :initform 1)))

(defclass shape-function-gimp (shape-function)
  ((order :initform 1)))

(defclass shape-function-bspline (shape-function)
  ((order :initform 2)))

(defclass shape-function-bspline-c2 (shape-function)
  ((order :initform 3)))

(defmacro make-shape-function (arg shape-form nD order &optional (shape-class 'shape-function))
  `(let ((svp (create-svp ,arg ,shape-form))
         (dsvp (create-dsvp ,arg ,shape-form)))
     (make-instance ,shape-class
                    :nD ,nD
                    :order ,order
                    :svp (nd-svp ,nD svp dsvp)
                    :dsvp (nd-dsvp ,nD svp dsvp))))


(defun make-shape-function-linear (nD h)
  (make-shape-function x (- 1d0 (abs (/ x h))) nD 1 'shape-function-linear))

(defun make-shape-function-bspline (nD h)
  (make-shape-function x
      (if (< (abs x) (/ h 2))
          (- (/ 3 4) (expt (/ (abs x) h) 2))
          (* (/ 1 8) (expt (- 3 (/ (* 2 (abs x)) h)) 2)))
                       nD 2 'shape-function-bspline))

(defun make-shape-function-gimp (nD h)
  (make-instance
   'shape-function-gimp
   :nD nD))

;; (symbolic-derivation:derive 'x
;;       '(if (< (abs x) (/ h 2)) 
;;           (- (/ 3 4) (expt (/ (abs x) h) 2))
;;           (* (/ 1 8) (expt (- 3 (/ (* 2 (abs x)) h)) 2))))
;(defun assemble-dsvp (dsvp)
;  "Assemble d/di to the strain-displacement matrix"
;  (let ((nD 2)
;        (dx (aref dsvp 0))
;        (dy (aref dsvp 1)))
;    (magicl:from-list (list dx 0d0 0d0 dy dy 0d0) '(3 2) :type 'double-float)))
(defun assemble-dsvp-1d (dsvp)
  "Assemble d/di to the strain-displacement matrix"
  (let ((dx (nth 0 dsvp)))
    (magicl:from-list (list dx) '(1 1) :type 'double-float)))

(declaim
 (inline assemble-dsvp-2d)
 (ftype (function (list) magicl:matrix/double-float) assemble-dsvp-2d))
(defun assemble-dsvp-2d (dsvp)
  "Assemble d/di to the strain-displacement matrix"
  (let* ((dx (nth 0 dsvp))
         (dy (nth 1 dsvp)))
    (magicl:from-array (make-array 6 :initial-contents (list dx 0d0
                                                             0d0 dy
                                                             dy dx))
                       '(3 2) :type 'double-float)))

(declaim
 (inline assemble-dstretch-2d)
 (ftype (function (list) magicl:matrix/double-float) assemble-dstretch-2d))
(defun assemble-dstretch-2d (dsvp)
  "Assemble d/di to the strain-displacement matrix"
  (let* ((dx (nth 0 dsvp))
         (dy (nth 1 dsvp)))
    (magicl:from-array (make-array 8 :initial-contents (list dx 0d0
                                                             0d0 dy
                                                             dy 0d0
                                                             0d0 dx
                                                             ))
                       '(4 2) :type 'double-float)))

(declaim
 (inline assemble-dstretch-2d-prealloc)
 (ftype (function (list magicl:matrix/double-float) magicl:matrix/double-float) assemble-dstretch-2d-prealloc))
(defun assemble-dstretch-2d-prealloc (dsvp result)
  "Assemble d/di to the strain-displacement matrix"
  (let* ((dx (nth 0 dsvp))
         (dy (nth 1 dsvp))
         (s (magicl::storage result)))
    (declare ;((simple-array double-float 8) s)
     ((simple-array double-float *) s)
     (double-float dx dy))
    (setf (aref s 0) dx
          ;(aref s 1) 0d0
          ;(aref s 2) 0d0
          (aref s 3) dy
          (aref s 4) dy
          ;(aref s 5) 0
          ;(aref s 6) 0
          (aref s 7) dx
          ))
    ;(magicl:from-array (make-array 8 :initial-contents (list dx 0d0
    ;                                                         0d0 dy
    ;                                                         dy 0d0
    ;                                                         0d0 dx
    ;                                                         ))
                       ;; '(4 2) :type 'double-float))
  result)

;; (time
;;  (let ((a (magicl:zeros '(1000 1)))
;;        (b (magicl:zeros '(1000 1))))
;;        (dotimes (i 1000000)
;;          (setf a (magicl:.+ a b))
;;          )))
;; (time
;;  (let ((a (magicl:zeros '(1000 1)))
;;        (b (magicl:zeros '(1000 1))))
;;    (dotimes (i 1000000)
;;      (magicl:.+ a b a)
;;      )))

(declaim
 (inline assemble-vorticity-2d)
 (ftype (function (list) magicl:matrix/double-float) assemble-vorticity-2d))
(defun assemble-vorticity-2d (dsvp)
  "Assemble d/di to the strain-displacement matrix"
  (let ((dx (nth 0 dsvp))
        (dy (nth 1 dsvp)))
    (magicl:from-array (make-array 6 :initial-contents (list 0d0 0d0
                                                             0d0 0d0
                                                             dy (- dx))) '(3 2) :type 'double-float)))

(defun assemble-vorticity-3d (dsvp)
  "Assemble d/di to the strain-displacement matrix"
  (let ((dx (nth 0 dsvp))
        (dy (nth 1 dsvp))
        (dz (nth 2 dsvp)))
    (magicl:from-list (list dx  0d0 0d0;xx
                            0d0 dy  0d0;yy
                            0d0 0d0  dz;zz
                            dy  (- dx)  0d0;yx
                            0d0 dz  (- dy) ;yz
                            dz  0d0 (- dx) ;xy
                            ) '(6 3) :type 'double-float)))

(defun assemble-dsvp-3d (dsvp)
  "Assemble d/di to the strain-displacement matrix"
  (let ((dx (nth 0 dsvp))
        (dy (nth 1 dsvp))
        (dz (nth 2 dsvp)))
    (magicl:from-list (list dx  0d0 0d0;xx
                            0d0 dy  0d0;yy
                            0d0 0d0  dz;zz
                            dy  dx  0d0;yx
                            0d0 dz  dy ;yz
                            dz  0d0 dx ;xy
                            ) '(6 3) :type 'double-float)))

(defun assemble-dsvp (nD dsvp)
  (case nD
    (1 (assemble-dsvp-1d dsvp))
    (2 (assemble-dsvp-2d dsvp))
    (3 (assemble-dsvp-3d dsvp))))


