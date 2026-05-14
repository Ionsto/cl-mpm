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
                         &key (tol 1d-9) (max-iters 10000)
                           (jacobi-precondition nil))
  "Solve problems of the form Ax = b, where the computation of A(x) is supplied as a lambda function, and b is some vector"
  (let ((b-norm (cl-mpm/fastmaths::mag b))
        (vector-size (magicl:nrows b)))
    (format t "Start~%")
    (labels ((operator (x)
               (if jacobi-precondition
                   (fast-./ (funcall A-operator x) jacobi-precondition)
                   (funcall A-operator x))))
      (if (= b-norm 0d0)
          ;;Trivial case of 0 being the answer
          (cl-mpm/utils::arb-matrix vector-size 1)
          ;;Nontrivial case
          (let* ((x ;; (magicl:rand (list vector-size 1))
                   (magicl:zeros (list vector-size 1))
                   )
                 (xn (cl-mpm/utils::deep-copy x))
                 (r (cl-mpm/utils::deep-copy x))
                 (crit tol)
                 (err crit)
                 (eigen-max (* 1d0 (estimate-max-eigenvalue #'operator vector-size))))
            (setf r (cl-mpm/fastmaths::fast-.-
                     ;; (funcall a-operator x)
                     (operator x)
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
                                (operator x)
                                b))
                       (setf err (cl-mpm/fastmaths:mag r))
                       (when (= (mod i 1000) 0)
                         (format t "Iter ~D - ~E~%" i err))))
            (when (> err crit)
              (error "Richardson didn't converge"))
            x)))))




;; (defun project-x (sim x)
;;   (cl-mpm::iterate-over-nodes
;;    (cl-mpm:sim-mesh sim)
;;    (lambda (n)
;;      (when active)
;;      )
;;    )
;;   )

(defun make-preconditioner (sparse-mat)
  (let* ((nrows (cl-mpm/utils::sparse-matrix-nrows sparse-mat))
         (pre (cl-mpm/utils::arb-matrix nrows 1)))
    (dotimes (r nrows)
      (setf (varef pre r) (/ 1d0 (cl-mpm/utils::sparse-matrix-aref sparse-mat r r))))
    pre))

(defun make-lumped-preconditioner (mat)
  (let* ((rowindex (cl-mpm/utils::sparse-matrix-rowindex mat))
         (cols (cl-mpm/utils::sparse-matrix-cols mat))
         (values (cl-mpm/utils::sparse-matrix-values mat))
         (nrows (cl-mpm/utils::sparse-matrix-nrows mat))
         (pre (cl-mpm/utils::arb-matrix nrows 1)))
    (dotimes (r nrows)
      (let ((col-0 (aref rowindex r))
            (col-1 (aref rowindex (1+ r))))
        (loop for c from col-0 below col-1 do
          (incf (varef pre r)
                (the double-float
                     (abs (aref values c)))))
        (setf (varef pre r) (varef pre r))))
    pre))

(defun solve-conjugant-gradients (A-operator b &key
                                                 (tol 1d-9)
                                                 (max-iters 10000)
                                                 (mask nil))
  (declare (function a-operator))
  (let (;; (pre nil)
        )
    ;; (when jacobi-precondition
    ;;   (setf pre (make-preconditioner A)))
  (labels ((mask-op (x)
             (if nil;mask
                 (cl-mpm/fastmaths:fast-.* x mask)
                 x))
           (mask-inplace (x)
               (when mask
                   (cl-mpm/fastmaths:fast-.* x mask x)))
           (operation (x)
             (funcall a-operator x))
           )

    (mask-inplace b)
    (let ((vector-size (magicl:nrows b))
          (b-norm (cl-mpm/fastmaths::mag-squared b)))
      ;; (pprint b-norm)
      (if (= b-norm 0d0)
          ;;Trivial case of 0 being the answer
          (cl-mpm/utils::arb-matrix vector-size 1)
          ;;Nontrivial case
          (let* ((x (cl-mpm/utils::arb-matrix vector-size 1))
                 (r (fast-.- b (operation x)))
                 (p (cl-mpm/utils::deep-copy r))
                 (ap (cl-mpm/utils::arb-matrix vector-size 1))
                 (rs-old (cl-mpm/fastmaths::mag-squared r))
                 (crit tol)
                 (rs-new crit)
                 (residual crit))
            (loop for i from 0 to max-iters
                  while (>= residual crit)
                  do
                     (progn
                       ;; (mask-inplace p)
                       (setf ap (operation p))
                       ;; (mask-inplace ap)
                       (let ((alpha
                               (/
                                rs-old
                                (dot p ap))))
                         (fast-.+ x (fast-scale p alpha) x)
                         (fast-.- r (fast-scale ap alpha) r)
                         ;; (setf r (mask-op (fast-.- b (operation x))))
                         (mask-inplace x)
                         ;; (mask-inplace r)
                         (setf rs-new (cl-mpm/fastmaths::mag-squared r))
                         ;(setf residual (/ rs-new b-norm))
                         (setf residual rs-new)
                         ;; (setf residual rs-new)
                         (unless (< residual crit)
                           (setf p
                                 (fast-.+
                                  r
                                  (fast-scale p (/ rs-new rs-old))))
                           ;; (mask-inplace p)
                           )
                         (when (= (mod (1+ i) (round (* max-iters 0.1d0))) 0)
                           (format t "CG Iter ~D ~E ~E ~E~%" i rs-old rs-new residual))
                         (setf rs-old rs-new)))
                  finally (when (> residual crit) 
                            (error "Conjugate gradients didn't converge")))
            ;; (format t "Solved in ~D iters~%" iters)
            (mask-inplace x)
            x))))))

(defun test-CG-squared ()
  (let ((A (cl-mpm/utils::matrix-from-list (list 1d0 8d0 3d0
                                                 2d0 6d0 9d0
                                                 0d0 2d0 1d0)))
        ;; (A (cl-mpm/utils::matrix-from-list (list 1d0 0d0 0d0
        ;;                                          0d0 6d0 0d0
        ;;                                          0d0 0d0 1d0)))
        (B (cl-mpm/utils:vector-from-list (list 0d0 2d0 9d0))))
    (pprint (magicl:linear-solve A B))
    (pprint (solve-conjugant-gradients-squared (lambda (v) (magicl:@ A v)) B))
    (pprint (solve-conjugant-gradients-squared (lambda (v) (magicl:@ A v)) B
                                               :jacobi-precondition (cl-mpm/utils::vector-from-list (list 0.001d0 1d0 100d0)))
            )
    (pprint (solve-preconditioned-conjugant-gradients-squared (lambda (v) (magicl:@ A v)) B
                                               :jacobi-precondition (cl-mpm/utils::vector-from-list (list 2d0 2d0 2d0)))
            )
    )
  )


(defun solve-conjugant-gradients-squared (A-operator b &key
                                                 (tol 1d-9)
                                                 (max-iters 10000)
                                                 (jacobi-precondition nil)
                                                 (mask nil))
  (declare (function a-operator))
  (let ()
  (labels ((mask-op (x)
             (if mask
                 (cl-mpm/fastmaths:fast-.* x mask)
                 x))
           (mask-inplace (x)
               (if mask
                   (cl-mpm/fastmaths:fast-.* x mask x)
                   x))
           (preconditioner-op (x)
             (if jacobi-precondition
                 (fast-./ x jacobi-precondition)
                 x))
           (preconditioner-into (x r)
             (if jacobi-precondition
                 (fast-./ x jacobi-precondition r)
                 (cl-mpm/utils::copy-into x r)))
           (preconditioner-inplace (x)
             (if jacobi-precondition
                 (fast-./ x jacobi-precondition x)
                 x))
           (operation (x)
             (mask-inplace
              (funcall a-operator x))))
    (mask-inplace b)
    (let ((vector-size (magicl:nrows b))
          (b-norm (cl-mpm/fastmaths::mag-squared b)))
      (assert (= (magicl:nrows b) vector-size))
      (when jacobi-precondition
        (assert (= (magicl:nrows jacobi-precondition) vector-size)))
      (if (= b-norm 0d0)
          ;;Trivial case of 0 being the answer
          (cl-mpm/utils::arb-matrix vector-size 1)
          ;;Nontrivial case
          (let* ((x (cl-mpm/utils::arb-matrix vector-size 1))
                 (r (fast-.- b (operation x)))
                 (r-tilde (cl-mpm/utils::deep-copy r))
                 (p (cl-mpm/utils::deep-copy r))
                 (u (cl-mpm/utils::deep-copy r))
                 (u-hat (cl-mpm/utils::deep-copy r))
                 (q (cl-mpm/utils::empty-copy r))
                 ;; (q-hat (cl-mpm/utils::empty-copy r))
                 (ap (cl-mpm/utils::arb-matrix vector-size 1))
                 (crit tol)
                 (rs-new crit)
                 (rs-old (cl-mpm/fastmaths:dot r-tilde r))
                 (residual crit)
                 (beta 0d0))
            ;; (mask-inplace r)
            ;; (mask-inplace r-tilde)

            (when (= (cl-mpm/fastmaths:dot r-tilde r) 0d0)
              (error "Conjugate gradients squared inital guess is orthoganal, rho=0"))
            (loop for i from 0 to max-iters
                  while (>= residual crit)
                  do
                     (progn
                       (setf rs-new (cl-mpm/fastmaths:dot r-tilde r))
                       (when (= rs-new 0d0)
                         (error "Conjugate gradients squared didn't converge, rho=0"))
                       (unless (= i 0)
                         (setf beta (/ rs-new rs-old))
                         (fast-.+ r (fast-scale q beta) u)
                         (fast-.+ u (fast-scale
                                     (fast-.+
                                      q
                                      (fast-scale p beta))
                                     beta) p))

                       ;; (mask-inplace r)
                       ;; (mask-inplace u)

                       ;;We then want to apply preconditioner to ap
                       (setf ap (operation (preconditioner-op p)))
                       (let ((alpha
                               (/
                                rs-new
                                (dot r-tilde ap))))
                         (fast-.- u (fast-scale ap alpha) q)
                         ;;Apply preconditioning to q
                         (preconditioner-into (fast-.+ u q) u-hat)
                         (fast-.+ x (fast-scale u-hat alpha) x)
                         ;; (mask-inplace x)
                         (fast-.- r (fast-scale (operation u-hat) alpha) r)
                         ;; (mask-inplace r)
                         ;; (setf residual (/ (cl-mpm/fastmaths::mag-squared r) b-norm))
                         (setf residual (cl-mpm/fastmaths::mag-squared r))
                         (when (= (mod (1+ i) (round (* max-iters 0.1d0))) 0)
                           (format t "Iter ~D ~E ~E ~E~%" i rs-old rs-new residual))
                         (setf rs-old rs-new)))
                  finally (when (> residual crit)
                       (error "Conjugate gradients didn't converge ~E" residual)))
            (mask-inplace x)
            x))))))
(defun solve-preconditioned-conjugant-gradients-squared (A-operator b &key
                                                                        (tol 1d-9)
                                                                        (max-iters 10000)
                                                                        (jacobi-precondition nil)
                                                                        (mask nil))
  ;; This is the improved preconditioned CGs from Formulation of a Preconditioned Algorithm
  ;; for the Conjugate Gradient Squared Method in Accordance with Its Logical Structure
  ;; Shoji Itoh1, Masaaki Sugihara (2015)
  (declare (function a-operator))
  (let ()
    ;; (when jacobi-precondition
    ;;   (setf pre (make-preconditioner A)))
  (labels ((mask-op (x)
             (if nil;mask
                 (cl-mpm/fastmaths:fast-.* x mask)
                 x))
           (mask-inplace (x)
               (if mask
                   (cl-mpm/fastmaths:fast-.* x mask x)
                   x))
           (preconditioner-op (x)
             (if jacobi-precondition
                 (fast-./ x jacobi-precondition)
                 x))
           (preconditioner-into (x r)
             (if jacobi-precondition
                 (fast-./ x jacobi-precondition r)
                 (cl-mpm/utils::copy-into x r)))
           (preconditioner-inplace (x)
             (if jacobi-precondition
                 (fast-./ x jacobi-precondition x)
                 x))
           (operation (x)
             (mask-inplace
              (funcall a-operator x))))
    (mask-inplace b)
    (let ((vector-size (magicl:nrows b))
          (b-norm (cl-mpm/fastmaths::mag-squared b)))
      (assert (= (magicl:nrows b) vector-size))
      (when jacobi-precondition
        (assert (= (magicl:nrows jacobi-precondition) vector-size)))
      (if (= b-norm 0d0)
          ;;Trivial case of 0 being the answer
          (cl-mpm/utils::arb-matrix vector-size 1)
          ;;Nontrivial case
          (let* ((x (cl-mpm/utils::arb-matrix vector-size 1))
                 (r (fast-.- b (operation x)))
                 (r-tilde (cl-mpm/utils::deep-copy r))
                 (p (cl-mpm/utils::deep-copy r))
                 (u (cl-mpm/utils::deep-copy r))
                 (u-hat (cl-mpm/utils::deep-copy r))
                 (q (cl-mpm/utils::empty-copy r))
                 ;; (q-hat (cl-mpm/utils::empty-copy r))
                 (ap (cl-mpm/utils::arb-matrix vector-size 1))
                 (crit tol)
                 (rs-new crit)
                 (rs-old (cl-mpm/fastmaths:dot r-tilde r))
                 (residual crit)
                 (beta 0d0)
                 (iters 0)
                 )
            (preconditioner-inplace r-tilde)
            ;; (preconditioner-inplace r)
            (loop for i from 0 to max-iters
                  while (>= residual crit)
                  do
                     (progn
                       (setf rs-new (cl-mpm/fastmaths:dot r-tilde (preconditioner-op r)))
                       (when (= rs-new 0d0)
                         (error "Conjugate gradients squared didn't converge, rho=0"))
                       (setf beta (/ rs-new rs-old))
                       (fast-.+ (preconditioner-op r) (fast-scale q beta) u)
                       (fast-.+ u (fast-scale
                                   (fast-.+
                                    q
                                    (fast-scale p beta))
                                   beta) p)

                       ;; (mask-inplace r)
                       ;; (mask-inplace u)

                       ;;We then want to apply preconditioner to ap
                       (setf ap (preconditioner-op (operation p)))
                       (let ((alpha
                               (/
                                rs-new
                                (dot r-tilde ap))))
                         (fast-.- u (fast-scale ap alpha) q)
                         ;;Apply preconditioning to q
                         (fast-.+ u q u-hat)

                         (fast-.+ x (fast-scale u-hat alpha) x)
                         (mask-inplace x)
                         (fast-.- r (fast-scale (operation u-hat) alpha) r)
                         ;; (mask-inplace r)
                         ;(setf residual (/ (cl-mpm/fastmaths::mag-squared r) b-norm))
                         (setf residual (cl-mpm/fastmaths::mag-squared r))
                         (when (= (mod (1+ i) (round (* max-iters 0.1d0))) 0)
                           (format t "Iter ~D ~E ~E ~E~%" i rs-old rs-new residual))
                         (incf iters)
                         (setf rs-old rs-new)))
                  finally (when (> residual crit)
                       (error "Conjugate gradients didn't converge ~E" residual)))
            (format t "CG solved in ~D iters ~E ~E ~E ~%" iters residual rs-new rs-old)
            (mask-inplace x)
            x))))))

(defun solve-conjugant-gradients-jacobi (A-operator b &key
                                                 (tol 1d-9)
                                                 (max-iters 10000)
                                                 (jacobi-precondition nil)
                                                 (mask nil))
  (declare (function a-operator))
  (let (;; (pre nil)
        )
    ;; (when jacobi-precondition
    ;;   (setf pre (make-preconditioner A)))
  (labels ((mask-op (x)
             (if nil;mask
                 (cl-mpm/fastmaths:fast-.* x mask)
                 x))
           (mask-inplace (x)
               (when mask
                   (cl-mpm/fastmaths:fast-.* x mask x)))
           (pre-inplace (x)
             (when jacobi-precondition
               (cl-mpm/fastmaths:fast-.* x jacobi-precondition x))
             x)
           (pre-into (source target)
             (if jacobi-precondition
               (cl-mpm/fastmaths:fast-.* source jacobi-precondition target)
               (cl-mpm/utils::copy-into source target)))
           (operation (x)
             (funcall a-operator x))
           (preconditioned-operation (x)
             (pre-inplace
              (funcall a-operator x))))

    (mask-inplace b)
    ;; (pre-inplace b)
    (let ((vector-size (magicl:nrows b))
          (b-norm (mag b)))
      (if (= b-norm 0d0)
          ;;Trivial case of 0 being the answer
          (cl-mpm/utils::arb-matrix vector-size 1)
          ;;Nontrivial case
          (let* ((x (cl-mpm/utils::arb-matrix vector-size 1))
                 (r (fast-.- b (operation x)))
                 (z (pre-inplace (cl-mpm/utils::deep-copy r)))
                 (p (cl-mpm/utils::deep-copy z))
                 (ap (cl-mpm/utils::arb-matrix vector-size 1))
                 (rs-old (cl-mpm/fastmaths::dot r z))
                 (crit tol)
                 (rs-new crit))
            (loop for i from 0 to max-iters
                  while (>= rs-new crit)
                  do
                     (progn
                       ;; (mask-inplace p)
                       (setf ap (operation p))
                       ;; (mask-inplace ap)
                       (let ((alpha
                               (/
                                rs-old
                                (dot p ap))))
                         (fast-.+ x (fast-scale p alpha) x)
                         ;; (fast-.- r (fast-scale ap alpha) r)
                         (setf r (fast-.- b (operation x)))
                         (pre-into r z)
                         ;; (mask-inplace x)
                         ;; (mask-inplace r)
                         (setf rs-new (cl-mpm/fastmaths::dot r z))
                         (unless (< rs-new crit)
                           (setf p
                                 (fast-.+
                                  z
                                  (fast-scale p (/ rs-new rs-old))))
                           )
                         (when (= (mod (1+ i) (round (* max-iters 0.1d0))) 0)
                           (format t "Iter ~D ~E ~E~%" i rs-old rs-new))
                         (setf rs-old rs-new)))
                  finally (when (> rs-new crit)
                       (error "Conjugate gradients didn't converge")))
            (mask-inplace x)
            x))))))


(defun test ()
  (let ((A (cl-mpm/utils::voigt-to-matrix (cl-mpm/utils:voigt-from-list (list 1d0 2d0 3d0 4d0 5d0 6d0))))
        (b (cl-mpm/utils:vector-from-list (list 1d0 2d0 3d0))))
    (pprint (magicl:linear-solve A B))
    (pprint (solve-conjugant-gradients (lambda (x) (magicl:@ A x)) b :tol 1d-15 :max-iters 1000))))
