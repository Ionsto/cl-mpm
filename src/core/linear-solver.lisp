(defpackage :cl-mpm/linear-solver
  (:use :cl
        :cl-mpm/fastmaths)
  (:import-from
   :magicl tref .+ .-
   )
  (:import-from
   :cl-mpm/utils varef)
  (:export
   #:solve-richardson
   #:solve-conjugant-gradients))

;; (declaim (optimize (debug 0) (safety 0) (speed 3)))
(declaim (optimize (debug 3) (safety 3) (speed 0)))
(in-package :cl-mpm/linear-solver)

(defun estimate-max-eigenvector (a-operator b-size &key (tol 1d-9)(max-iters 10000))
  (declare (function a-operator))
  (let* ((crit tol)
         (err crit)
         (bn (magicl:rand (list b-size 1)))
         (bn1 (cl-mpm/utils::arb-matrix  b-size 1))
         )
    (loop for i from 0 below max-iters
          while (>= err crit)
          do (progn
               (setf bn1 (funcall a-operator bn))
               (let ((norm (cl-mpm/fastmaths:mag bn1)))
                 (cl-mpm/fastmaths::fast-scale! bn1 (/ 1d0 norm))

                 (setf err
                       (/ (cl-mpm/fastmaths:mag (cl-mpm/fastmaths:fast-.- bn bn1))
                        norm))
                 (setf bn bn1))))
    bn))
(defun estimate-max-eigenvalue (a-operator b-size &key (max-iters 100))
  (declare (function a-operator))
  (let ((vn (estimate-max-eigenvector a-operator b-size :max-iters max-iters)))
    (/
     (cl-mpm/fastmaths:dot (funcall a-operator vn) vn)
     (cl-mpm/fastmaths:dot vn vn))))

(defun solve-richardson (A-operator b
                         &key (tol 1d-9) (max-iters 1000)
                                        )
  "Solve problems of the form Ax = b, where the computation of A(x) is supplied as a lambda function, and b is some vector"
  (let ((b-norm (cl-mpm/fastmaths::mag b))
        (vector-size (magicl:nrows b))
        )
    (if (= b-norm 0d0)
        ;;Trivial case of 0 being the answer
        (cl-mpm/utils::arb-matrix vector-size 1)
        ;;Nontrivial case
        (let* ((x (magicl:rand (list vector-size 1)))
               (xn (cl-mpm/utils::deep-copy x))
               (r (cl-mpm/utils::deep-copy x))
               (crit tol)
               (err crit)
               (eigen-max (estimate-max-eigenvalue A-operator vector-size)))

          (setf r (cl-mpm/fastmaths::fast-.-
                   (funcall a-operator x)
                   b))
          (loop for i from 0 to max-iters
                while (>= err crit)
                do
                   (progn
                     (cl-mpm/fastmaths::fast-.+
                      x
                      (cl-mpm/fastmaths::fast-scale!
                       r
                       (/ -1d0 eigen-max))
                      x)
                     (cl-mpm/utils::copy-into x xn)
                     (setf r (cl-mpm/fastmaths::fast-.-
                              (funcall a-operator x)
                              b))
                     (setf err
                           (cl-mpm/fastmaths:mag r)
                           ;; (/
                           ;;  (cl-mpm/fastmaths:mag
                           ;;   (cl-mpm/fastmaths::fast-.-
                           ;;    (funcall a-operator x)
                           ;;    b))
                           ;;  b-norm)
                           )
                     ;; (format t "error ~E~%" err)
                     ))

          x))))




;; (defun project-x (sim x)
;;   (cl-mpm::iterate-over-nodes
;;    (cl-mpm:sim-mesh sim)
;;    (lambda (n)
;;      (when active)
;;      )
;;    )
;;   )

(defun solve-conjugant-gradients (A-operator b &key (tol 1d-9) (max-iters 10000))
  (declare (function a-operator))
  (let ((vector-size (magicl:nrows b))
        (b-norm (mag b)))
    (if (= b-norm 0d0)
        ;;Trivial case of 0 being the answer
        (cl-mpm/utils::arb-matrix vector-size 1)
        ;;Nontrivial case
        (let* ((x (cl-mpm/utils::arb-matrix vector-size 1))
               (r (fast-.- b (funcall a-operator x)))
               (p r)
               (ap (cl-mpm/utils::arb-matrix vector-size 1))
               (rs-old (cl-mpm/fastmaths::mag-squared r))
               (crit tol)
               (rs-new crit))
          (loop for i from 0 to max-iters
                while (>= rs-new crit)
                do
                   (progn
                     (setf ap (funcall a-operator p))
                     (let ((alpha
                             (/
                              rs-old
                              (dot p ap))))
                       (setf rs-new (cl-mpm/fastmaths::mag-squared r))
                       (fast-.+ x (fast-scale p alpha) x)
                       (fast-.- r (fast-scale p alpha) r)
                       (unless (< rs-new crit)
                         (setf p
                               (fast-.+
                                r
                                (fast-scale p (/ rs-new rs-old)))))
                       (setf rs-old rs-new))
                     ))
          x))))
