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
;; (declaim (optimize (debug 3) (safety 3) (speed 0)))

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

(declaim (inline shape-gimp-fast)
         (ftype (function (double-float double-float double-float) double-float) shape-gimp-fast))
(defun shape-gimp-fast (x l h)
  (declare (type double-float x l h))
  ;; (when (= 0d0 l)
  ;;   (error "Zero length domain!"))
  (let ((ax (abs x)))
    (declare (type double-float ax))
    (cond
      ((and (< (- h l) ax) (<= ax (+ h l)))
       (/ (expt (- (+ h l) ax) 2d0) (* 4d0 h l)))
      ((and (< l ax) (<= ax (- h l)))
       (- 1d0 (/ ax h)))
      ((<= ax l)
       (- 1d0 (/ (+ (* x x) (* l l)) (* 2d0 h l))))
      (t
       0d0))))

(declaim (inline shape-gimp-branchless)
         (ftype (function (double-float double-float double-float) double-float) shape-gimp-fast))
(defun shape-gimp-branchless (x l h)
  (declare (type double-float x l h))
  (let ((ax (abs x))
        (res 0d0))
    (declare (type double-float ax res))
    (when (and (< (- h l) ax) (<= ax (+ h l)))
      (incf res (/ (expt (- (+ h l) ax) 2d0) (* 4d0 h l))))
    (when (and (< l ax) (<= ax (- h l)))
      (incf res(- 1d0 (/ ax h))))
    (when (<= ax l)
      (incf res (- 1d0 (/ (+ (* x x) (* l l)) (* 2d0 h l)))))
    res))

;; (let* ((h 1.0d0)
;;        (x (loop for x from -2d0 upto 2d0 by 0.01d0 collect (* h x)))
;;        (node 0)
;;        (l 1.2d0)
;;        )
;;   (vgplot:figure)
;;   (vgplot:plot
;;    x
;;    (mapcar (lambda (y) (shape-gimp-fast y (* l 0.5d0) h)) x) "0"
;;    x
;;    ;; (mapcar (lambda (y) (shape-gimp-dsvp y (* l 0.5d0) h)) x) "0"
;;    ;; x
;;    (mapcar (lambda (y) (shape-gimp-fbar y (* l 0.5d0) h)) x) "fbar"
;;    ))

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

(declaim (inline shape-gimp-fbar)
         (ftype (function (double-float double-float double-float) double-float) shape-gimp-fbar))
(defun shape-gimp-fbar (x l h)
  (declare (type double-float x l h))
  (cond
    ((and (< (- (+ h l)) x) (<= x (- l h)))
     (/ (+ x h l) (* 4d0 l)))
    ((and (< (- l h) x) (<= x (- h l)))
     0.5d0)
    ((and (< (- h l) x) (<= x (+ h l)))
     (/ (- (+ h l) x) (* 4d0 l)))
    (t
     0d0)))
(declaim (inline shape-gimp-fbar-dsvp)
         (ftype (function (double-float double-float double-float) double-float) shape-gimp-fbar-dsvp))
(defun shape-gimp-fbar-dsvp (x l h)
  (declare (type double-float x l h))
  (cond
    ((and (< (- (+ h l)) x) (<= x (- l h)))
     (/ 1d0 (* 4d0 l)))
    ((and (< (- l h) x) (<= x (- h l)))
     0d0)
    ((and (< (- h l) x) (<= x (+ h l)))
     (- (/ 1d0 (* 4d0 l))))
    (t
     0d0)))


;; (let* ((x (loop for x from -4d0 to 4d0 by 0.01d0 collect x))
;;        (h 0.10d0)
;;        (offset 0.0d0)
;;        (fbar (mapcar (lambda (x) (shape-gimp-fbar (+ x offset) h 1d0)) x))
;;        (gimp (mapcar (lambda (x) (shape-gimp-fast (+ x offset) h 1d0)) x))
;;        )
;;   (vgplot:close-all-plots)
;;   (vgplot:figure)
;;   (vgplot:plot x fbar ""
;;                x gimp ""))



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
                    (t
                     i)
                    ))))
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
  ;; (print nodes)
  ;; (let ((nc (count t nodes)))
  ;;   (cond
  ;;     ((= nc 1)
  ;;      (if (nth 3 nodes)
  ;;          (bspline (make-bspline-knots nodes h) (+ eta) (+ node 2) 1)
  ;;          (bspline (make-bspline-knots nodes h) (+ eta) (+ node 3) 1)))
  ;;     (t
  ;;      (bspline (make-bspline-knots nodes h) (+ eta) (+ node 2) 2))
  ;;     )
  ;;   )
  (bspline (make-bspline-knots nodes h) (+ eta) (+ node 2) 2)
  )
(defun nodal-bspline-dsvp (nodes eta node h)
  (setf (nth 4 nodes) (and (nth 3 nodes) (nth 5 nodes)))
  ;; (let ((nc (count t nodes)))
  ;;   (cond
  ;;     ((= nc 1)
  ;;      (if (nth 3 nodes)
  ;;          (bspline-dsvp (make-bspline-knots nodes h) (+ eta) (+ node 3) 1)
  ;;          (bspline-dsvp (make-bspline-knots nodes h) (+ eta) (+ node 4) 1)))
  ;;     (t
  ;;      (bspline-dsvp (make-bspline-knots nodes h) (+ eta) (+ node 3) 2))
  ;;     )
  ;;   )
  (bspline-dsvp (make-bspline-knots nodes h) (+ eta) (+ node 3) 2)
  ;(bspline-dsvp (make-bspline-knots nodes h) (+ eta) (+ node 3) 2)
  )

;; (let* ((h 1.0d0)
;;       (x (loop for x from -2d0 upto 2d0 by 0.01d0 collect (* h x)))
;;       (node 0)
;;       ;; (node-list '(nil nil nil nil T t nil nil nil))
;;        (node-list '(nil nil nil t T t t t t))
;;        )
;;   (vgplot:figure)
;;   (vgplot:plot
;;    x
;;    (mapcar (lambda (y)
;;              (nodal-bspline '(nil nil nil nil T t t t t) (+ y h) 0 h)) x) "-1a"
;;    x
;;    (mapcar (lambda (y)
;;              (nodal-bspline node-list (+ y) -1 h)) x) "-1"
;;     ;; x
;;     ;; (mapcar (lambda (y)
;;     ;;           (nodal-bspline-dsvp node-list (+ y) -1 h)
;;     ;;           ) x) "d-1"
;;    x
;;    (mapcar (lambda (y)
;;              (nodal-bspline node-list (+ y) 0 h)
;;              ) x) "0"
;;    ;; x
;;    ;; (mapcar (lambda (y)
;;    ;;           (nodal-bspline-dsvp node-list (+ y) 0 h)
;;    ;;           ) x) "d0"
;;    x
;;    (mapcar (lambda (y)
;;              (nodal-bspline node-list (+ y) 1 h)
;;              ) x) "1"
;;    ;; x
;;    ;; (mapcar (lambda (y)
;;    ;;           (nodal-bspline-dsvp node-list (+ y) 1 h)
;;    ;;           ) x) "d1"
;;   ;; x
;;   ;; (mapcar (lambda (y) (nodal-bspline-dsvp
;;   ;;                      node-list
;;   ;;                      (+ y h) node 1d0)) x) "-1"
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
                    :svp nil;(nd-svp ,nD svp dsvp)
                    :dsvp nil;(nd-dsvp ,nD svp dsvp)
                    )))


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
 (inline assemble-dsvp-2d-prealloc)
 (ftype (function (list magicl:matrix/double-float) magicl:matrix/double-float) assemble-dsvp-2d-prealloc))
(defun assemble-dsvp-2d-prealloc (dsvp mat)
  "Assemble d/di to the strain-displacement matrix"
  (let* ((dx (nth 0 dsvp))
         (dy (nth 1 dsvp)))
    (setf
     (magicl:tref mat 0 0) dx
     (magicl:tref mat 1 1) dy
     (magicl:tref mat 2 0) dy
     (magicl:tref mat 2 1) dx))
  mat)

(declaim
 (inline assemble-dstretch-2d)
 (ftype (function (list) magicl:matrix/double-float) assemble-dstretch-2d))
(defun assemble-dstretch-2d (dsvp)
  "Assemble d/di to the strain-displacement matrix"
  (let* ((dx (nth 0 dsvp))
         (dy (nth 1 dsvp)))
    (cl-mpm/utils:stretch-dsvp-3d-zeros)
    ))

(declaim
 (inline assemble-dstretch-2d-prealloc)
 (ftype (function (list magicl:matrix/double-float) magicl:matrix/double-float) assemble-dstretch-2d-prealloc))
(defun assemble-dstretch-2d-prealloc (dsvp result)
  (declare (list dsvp)
           (magicl:matrix/double-float result)
           ;; (optimize (speed 3) (safety 0) (debug 0))
           )
  "Assemble d/di to the strain-displacement matrix"
  (destructuring-bind (dx dy dz) dsvp
    (let* ((s (magicl::matrix/double-float-storage result)))
      (declare (double-float dx dy dz))
      (setf
                                        ;dx/dx
       (aref s (+ 0 (* 9 0))) dx
                                        ;dy/dy
       (aref s (+ 1 (* 9 1))) dy
                                        ;dz/dz
       ;; (aref s (+ 2 (* 9 2))) dz
                                        ;Dy/dx
       (aref s (+ 3 (* 9 0))) dy
                                        ;Dx/dy
       (aref s (+ 4 (* 9 1))) dx
                                        ;; dz/dx
       ;; (aref s (+ 5 (* 9 0))) dz 

       ;; (aref s (+ 6 (* 9 2))) dx
       ;;                                  ;Dz/dy
       ;; (aref s (+ 7 (* 9 1))) dz
                                        ;; Dy/dz
       ;; (aref s (+ 8 (* 9 2))) dy
       )
       ;);)
      result))
  
  result)

(defun assemble-dstretch-3d (dsvp)
  (assemble-dstretch-3d-prealloc dsvp (cl-mpm/utils::stretch-dsvp-3d-zeros)))

(declaim (inline assemble-dstretch-3d-prealloc)
         (ftype (function (list magicl:matrix/double-float) magicl:matrix/double-float)
                assemble-dstretch-3d-prealloc))
(defun assemble-dstretch-3d-prealloc (dsvp result)
  (declare (list dsvp)
           (magicl:matrix/double-float result)
           ;; (optimize (speed 3) (safety 0) (debug 0))
           )
  "Assemble d/di to the strain-displacement matrix"
  (destructuring-bind (dx dy dz) dsvp
    (let* ((s (magicl::matrix/double-float-storage result)))
      (declare (double-float dx dy dz))
      (setf
                                        ;dx/dx
       (aref s (+ 0 (* 9 0))) dx
                                        ;dy/dy
       (aref s (+ 1 (* 9 1))) dy
                                        ;dz/dz
       (aref s (+ 2 (* 9 2))) dz
                                        ;Dy/dx
       (aref s (+ 3 (* 9 0))) dy
                                        ;Dx/dy
       (aref s (+ 4 (* 9 1))) dx
                                        ;dz/dx
       (aref s (+ 5 (* 9 0))) dz 
                                        ;dx/dz
       (aref s (+ 6 (* 9 2))) dx
                                        ;Dz/dy
       (aref s (+ 7 (* 9 1))) dz
                                        ;Dy/dz
       (aref s (+ 8 (* 9 2))) dy
       )
       ;);)
      result)))

;; (time
;;  (let ((a (magicl:zeros '(1000 1)))
;;        (b (magicl:zeros '(1000 1))))
;;        (dotimes (i 1000000)
;;          (setf a (cl-mpm/fastmath::fast-.+ a b))
;;          )))
;; (time
;;  (let ((a (magicl:zeros '(1000 1)))
;;        (b (magicl:zeros '(1000 1))))
;;    (dotimes (i 1000000)
;;      (cl-mpm/fastmath::fast-.+ a b a)
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
  (assemble-dsvp-3d-prealloc dsvp (cl-mpm/utils::dsvp-3d-zeros))
  )


(defun assemble-dsvp-3d-prealloc (dsvp mat)
  "Assemble d/di to the strain-displacement matrix"
  (destructuring-bind (dx dy dz) dsvp
      (setf
       (magicl:tref mat 0 0) dx
       (magicl:tref mat 1 1) dy
       (magicl:tref mat 2 2) dz

       (magicl:tref mat 3 1) dz
       (magicl:tref mat 3 2) dy

       (magicl:tref mat 4 0) dz
       (magicl:tref mat 4 2) dx

       (magicl:tref mat 5 0) dy
       (magicl:tref mat 5 1) dx
       )
    mat)
  )
;; (defun @-assemble-dsvp-3d-prealloc (dsvp voigt result)
;;   "Assemble d/di to the strain-displacement matrix"
;;   (let ((res (cl-mpm/utils::fast-storage result))
;;         (v (cl-mpm/utils::fast-storage voigt)))
;;     (declare ((simple-array double-float (3)) res)
;;              ((simple-array double-float (6)) v)
;;              )
;;     (destructuring-bind (dx dy dz) dsvp
;;       (declare (double-float dx dy dz))
;;       (setf
;;        (aref res 0) (+
;;                      (* dx (aref v 0))
;;                      (* dz (aref v 4))
;;                      (* dy (aref v 5)))
;;        (aref res 1) (+
;;                      (* dy (aref v 1))
;;                      (* dz (aref v 3))
;;                      (* dx (aref v 5)))
;;        (aref res 2) (+
;;                      (* dz (aref v 2))
;;                      (* dy (aref v 3))
;;                      (* dx (aref v 4)))

;;        )
;;       ;; (setf
;;       ;;  (magicl:tref mat 0 0) dx
;;       ;;  (magicl:tref mat 1 1) dy
;;       ;;  (magicl:tref mat 2 2) dz

;;       ;;  (magicl:tref mat 3 1) dz
;;       ;;  (magicl:tref mat 3 2) dy

;;       ;;  (magicl:tref mat 4 0) dz
;;       ;;  (magicl:tref mat 4 2) dx

;;       ;;  (magicl:tref mat 5 0) dy
;;       ;;  (magicl:tref mat 5 1) dx
;;       ;;  )
;;       ;; mat
;;       ))
;;   )
 

(declaim (ftype (function (list magicl:matrix/double-float magicl:matrix/double-float) (values)) @-combi-assemble-dstretch-3d))
(defun @-combi-assemble-dstretch-3d (grads vel stretch)
  "Assemble d/di to the strain-displacement matrix"
  (let ((res (cl-mpm/utils::fast-storage stretch))
        (v (cl-mpm/utils::fast-storage vel)))
    (declare ((simple-array double-float (9)) res)
             ((simple-array double-float (3)) v)
             )
    (destructuring-bind (dx dy dz) grads
      (declare (double-float dx dy dz))
      (macrolet ((component (x y comp-pairs)
                   (declare (fixnum x y))
                   `(progn
                      (setf (aref res ,(the fixnum (+ y (* x 3))))
                            (+
                             (aref res ,(the fixnum (+ y (* x 3))))
                             ,(loop for comp in comp-pairs
                                    append
                                    (destructuring-bind (grad index) comp
                                      `(* ,grad (aref v ,index)))
                                    )))
                      ))

                 )
        (component 0 0 ((dx 0)))
        (component 1 1 ((dy 1)))
        (component 2 2 ((dz 2)))
        (component 0 1 ((dx 1)))
        (component 0 2 ((dx 2)))
        (component 1 0 ((dy 0)))
        (component 2 0 ((dz 0)))
        (component 1 2 ((dy 2)))
        (component 2 1 ((dz 1)))
        )))
  )



(defun assemble-dsvp (nD dsvp)
  (case nD
    (1 (assemble-dsvp-1d dsvp))
    (2 (assemble-dsvp-2d dsvp))
    (3 (assemble-dsvp-3d dsvp))))

(defun grads-3d (weights linear-grads)
  "Take weights and gradients of weights and assemble them into"
  (mapcar #'* linear-grads
          (list
           (* (nth 1 weights) (nth 2 weights))
           (* (nth 0 weights) (nth 2 weights))
           (* (nth 0 weights) (nth 1 weights)))))
(defun grads-2d (weights linear-grads)
  (mapcar #'* linear-grads (nreverse weights))
  )

