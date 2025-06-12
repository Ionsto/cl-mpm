(in-package :cl-mpm)
;;All the various ways of iterating over the mesh
;; (declaim (optimize (debug 0) (safety 0) (speed 3)))
(declaim #.cl-mpm/settings:*optimise-setting*)
;; (declaim (optimize (debug 3) (safety 3) (speed 0)))


;; (require 'sb-concurrency)
;; (defparameter *workers* nil)
;; (defparameter *workers-array* nil)
;; (defparameter *work-queue* (sb-concurrency:make-queue))
;; (defparameter *workers-run* (sb-thread:make-semaphore))
;; (defparameter *workers-finish* (sb-thread:make-semaphore))
;; (defparameter *workers-kill* nil)
;; (defparameter *workers-func* (lambda ()))
;; (defun make-workers ()
;;   (unless *workers*
;;     (let ((thread-count (lparallel:kernel-worker-count)))
;;       (setf *workers*
;;             (loop for i fixnum from 0 below thread-count
;;                   collect
;;                   (sb-thread:make-thread
;;                    (lambda (thread-number thread-count)
;;                      (loop
;;                        while (not *workers-kill*)
;;                        do
;;                        (sb-thread:wait-on-semaphore *workers-run*)
;;                        (unless *workers-kill*
;;                          (let ((func (the function (sb-concurrency:dequeue *work-queue*))))
;;                            (funcall func)))
;;                        (sb-thread:signal-semaphore *workers-finish*))
;;                      (values))
;;                    :arguments (list i thread-count)
;;                    ))))))
;; (defun kill-workers ()
;;   (when *workers*
;;     (setf *workers-kill* t)
;;     (sb-thread:signal-semaphore *workers-run* (length *workers*))
;;     (loop for worker in *workers* do (sb-thread:join-thread worker))
;;     (setf *workers* nil)
;;     (setf *workers-kill* nil)
;;     ))
;; (defun better-pdotimes (array func)
;;   (declare (vector array)
;;            (function func))
;;   (make-workers)
;;   (setf *workers-array* array)
;;   (setf *workers-func* func)
;;   (let ((length (length *workers*)))
;;     (loop for i from 0 below length
;;           do (sb-concurrency:enqueue
;;               (let ((i-capt i)
;;                     (length-capt length))
;;                 (lambda () (funcall func i-capt length-capt))) *work-queue*)))
;;   ;; (pprint *work-queue*)
;;   (sb-thread:signal-semaphore *workers-run* (length *workers*))
;;   (sb-thread:wait-on-semaphore *workers-finish* :n (length *workers*))
;;   )

;; (defmacro omp (array func)
;;   (declare (vector array)
;;            (function func))
;;   (let ((job-sym (gensym)))
;;     `(progn
;;        (let ((,job-sym (lambda (thread-number thread-count)
;;                          (let* ((array cl-mpm::*workers-array*)
;;                                 (total-size (length array))
;;                                 (block-size (round (ceiling total-size thread-count))))
;;                            (declare (vector array))
;;                            (loop for j fixnum from (* thread-number block-size) below (min (* (+ thread-number 1) block-size) total-size)
;;                                  do (funcall ,func j))
;;                            (values))
;;                          )))
;;          (declare (function ,job-sym))
;;          (cl-mpm::better-pdotimes ,array ,job-sym)
;;          ))))

;; (defun shit-pdotimes (array func)
;;   (declare (vector array)
;;            (function func))
;;   (let* ((thread-count 12)
;;          (total-size (length array))
;;          (block-size (round (ceiling total-size thread-count))))
;;     (declare (fixnum block-size thread-count))
;;     (let ((threads (loop for i fixnum from 0 below thread-count
;;                          collect
;;                          (progn
;;                            (sb-thread:make-thread
;;                             (lambda (thread-number)
;;                               (loop for j fixnum from (* thread-number block-size) below (min (* (+ thread-number 1) block-size) total-size)
;;                                     do (funcall func j))
;;                               (values))
;;                             :arguments i
;;                             )))))
;;       (loop for thread in threads
;;             do (sb-thread:join-thread thread)))))


(declaim (inline iterate-over-nodes)
         (ftype (function (cl-mpm/mesh::mesh function) (values)) iterate-over-nodes))
(defun iterate-over-nodes (mesh func)
  "Helper function for iterating over all nodes in a mesh
   Calls func with only the node"
  (declare (type function func))
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
    (declare (type (array cl-mpm/particle:particle) nodes))
    ;; (dotimes (i (array-total-size nodes))
    ;;   (let ((node (row-major-aref nodes i)))
    ;;     (funcall func node)))
    (lparallel:pdotimes (i (array-total-size nodes))
      (let ((node (row-major-aref nodes i)))
        (when node
          (funcall func node))))
    ;; (let ((wrapper (make-array (array-total-size nodes) :displaced-to nodes)))
    ;;   (omp
    ;;    wrapper
    ;;    (lambda (i)
    ;;      (funcall func (aref wrapper i)))))
    )
  (values))

(declaim (inline iterate-over-nodes-serial)
         (ftype (function (cl-mpm/mesh::mesh function) (values)) iterate-over-nodes-serial)
         )
(defun iterate-over-nodes-serial (mesh func)
  "Helper function for iterating over all nodes in a mesh - in a serial way
Calls func with only the node"
  (declare (type function func))
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
    (declare (type (array cl-mpm/particle:particle) nodes))
    (dotimes (i (array-total-size nodes))
      (let ((node (row-major-aref nodes i)))
        (when node
          (funcall func node)))))
  (values))

(declaim
 (ftype (function (cl-mpm/mesh::mesh function) (values)) iterate-over-cells))
(defun iterate-over-cells (mesh func)
  "Helper function for iterating over all nodes in a mesh
   Calls func with only the node"
  (declare (type function func))
  (let ((cells (cl-mpm/mesh::mesh-cells mesh)))
    (declare (type (array (or cl-mpm/mesh::cell null)) cells))
    (lparallel:pdotimes (i (array-total-size cells))
      (let ((cell (row-major-aref cells i)))
        (when cell
          (funcall func cell)))))
  (values))

(defun iterate-over-cells-serial (mesh func)
  "Helper function for iterating over all nodes in a mesh
   Calls func with only the node"
  (declare (type function func))
  (let ((cells (cl-mpm/mesh::mesh-cells mesh)))
    (declare (type (array (or cl-mpm/mesh::cell null)) cells))
    (dotimes (i (array-total-size cells))
      (let ((cell (row-major-aref cells i)))
        (when cell
          (funcall func cell)))))
  (values))

(defun iterate-over-bcs (sim func)
  "Helper function for iterating over all nodes in a mesh
   Calls func with only the node"
  (declare (type function func))
  (let ((bcs (sim-bcs sim)))
    (lparallel:pdotimes (i (array-total-size bcs))
      (let ((bc (aref bcs i)))
        (when bc
          (funcall func bc)))))
  (values))

(defun iterate-over-bcs-serial (sim func)
  "Helper function for iterating over all nodes in a mesh
   Calls func with only the node"
  (declare (type function func))
  (let ((bcs (sim-bcs sim)))
    (dotimes (i (array-total-size bcs))
      (let ((bc (aref bcs i)))
        (when bc 
          (funcall func bc)))))
  (values))

(defun iterate-over-bcs-force-serial (sim func)
  "Helper function for iterating over all nodes in a mesh
   Calls func with only the node"
  (declare (type function func))
  (let ((bcs-f (sim-bcs-force-list sim)))
    (loop for bcs in bcs-f
          do (dotimes (i (array-total-size bcs))
               (let ((bc (aref bcs i)))
                 (when bc
                   (funcall func bc))))))
  (values))



(declaim (ftype (function
                 ((vector cl-mpm/particle:particle *)
                  function)
                 (values))
                iterate-over-mps))
(defun iterate-over-mps (mps func)
  "Helper function for iterating over all nodes in a mesh
   Calls func with only the node"
  (declare (type function func)
           (type (array cl-mpm/particle:particle) mps))
  ;; (dotimes (i (length mps))
  ;;   (funcall func (aref mps i)))
  (lparallel:pdotimes (i (length mps))
                      (funcall func (aref mps i)))
  ;; (omp
  ;;  mps
  ;;   (lambda (i)
  ;;     (funcall func (aref mps i))))
  ;; (better-pdotimes
  ;;  mps
  ;;  (lambda (i)
  ;;    (funcall func (aref mps i))))

  (values))

(defun reduce-over-nodes (mesh map reduce)
  (with-accessors ((nodes cl-mpm/mesh:mesh-nodes))
      mesh
    (lparallel:pmap-reduce
     map
     reduce
     (make-array (array-total-size nodes)
                 :displaced-to nodes))))

(defun reduce-over-mps (mps map reduce)
  (lparallel:pmap-reduce
   map
   reduce
   mps))

(defun iterate-over-mps-serial (mps func)
  "Helper function for iterating over all nodes in a mesh
   Calls func with only the node"
  (declare (type function func)
           (type (array cl-mpm/particle:particle) mps))
  (loop for mp across mps
        do (funcall func mp))
  (values))




(defgeneric iterate-over-neighbours-shape (mesh shape-func mp func)
  (:documentation "For a given shape function iterate over relevant nodes and call func with mesh, mp, node, weight and gradients"))

(declaim
(inline iterate-over-neighbours)
         (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle function) (values)) iterate-over-neighbours))
(defun iterate-over-neighbours (mesh mp func)
  "For a given mesh and mp, iterate over all the neighbours with
weight greater than 0, calling func with the mesh, mp, node, svp, and grad"
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (if (> (length (cl-mpm/particle::mp-cached-nodes mp)) 0)
      (iterate-over-neighbours-cached mesh mp func)
      (create-node-cache mesh mp func))

  ;; (if (= (the fixnum (cl-mpm/mesh:mesh-nd mesh)) 2)
  ;;     (iterate-over-neighbours-shape-gimp-simd
  ;;        mesh mp func)
  ;;     (iterate-over-neighbours-shape-gimp-3d
  ;;        mesh mp func))
  ;; (iterate-over-neighbours-shape-gimp-2d mesh mp func)
  ;; (iterate-over-neighbours-shape-gimp-simd mesh mp func)
  ;; (iterate-over-neighbours-shape-gimp-3d mesh mp func)
  ;; (iterate-over-neighbours-shape-gimp mesh mp func)
  ;; (iterate-over-neighbours-shape-linear mesh mp func)
  (values)
  )

(declaim (inline create-node-cache))
(defun create-node-cache (mesh mp func)
  "A function iterating over neighbours executing function, while also caching the relevent node, gradients and weights"
  (declare (function func))
  (with-accessors ((nodes cl-mpm/particle::mp-cached-nodes))
      mp
    ;;Simple if statement - we take the hit
    (if (= (the fixnum (cl-mpm/mesh:mesh-nd mesh)) 2)
        (iterate-over-neighbours-shape-gimp-simd
         mesh mp
         (lambda (mesh mp node svp grads fsvp fgrads)
           (destructuring-bind (gx gy gz) grads
             (declare (ignore gz))
             (destructuring-bind (gfx gfy gfz) fgrads
               (declare (ignore gfz))
               (vector-push-extend
                (cl-mpm/particle::make-node-cache
                 node
                 svp
                 gx
                 gy
                 0d0
                 fsvp
                 gfx
                 gfy
                 0d0
                 )
                nodes)))
           (funcall func mesh mp node svp grads fsvp fgrads)))
        (iterate-over-neighbours-shape-gimp-3d
         mesh mp
         (lambda (mesh mp node svp grads fsvp fgrads)
           (destructuring-bind (gx gy gz) grads
             (destructuring-bind (gfx gfy gfz) fgrads
               (vector-push-extend
                (cl-mpm/particle::make-node-cache
                 node
                 svp
                 gx
                 gy
                 gz
                 fsvp
                 gfx
                 gfy
                 gfz
                 )
                nodes)))
           (funcall func mesh mp node svp grads fsvp fgrads))))))

(declaim (inline iterate-over-neighbours-cached))
(defun iterate-over-neighbours-cached (mesh mp func)
  "If a node iteration cache has been generated we loop over the data list"
  (declare (function func))
  (loop for nc across (the (vector cl-mpm/particle::node-cache *) (cl-mpm/particle::mp-cached-nodes mp))
        do
         (funcall func mesh mp
                  (cl-mpm/particle::node-cache-node nc)
                  (cl-mpm/particle::node-cache-weight nc)
                  (list
                   (cl-mpm/particle::node-cache-grad-x nc)
                   (cl-mpm/particle::node-cache-grad-y nc)
                   (cl-mpm/particle::node-cache-grad-z nc))
                  (cl-mpm/particle::node-cache-weight-fbar nc)
                  (list
                   (cl-mpm/particle::node-cache-grad-fbar-x nc)
                   (cl-mpm/particle::node-cache-grad-fbar-y nc)
                   (cl-mpm/particle::node-cache-grad-fbar-z nc))
                  )))

;;This is one method of dispatching over different types of shape functions
(defmethod iterate-over-neighbours-shape (mesh (shape-func cl-mpm/shape-function:shape-function-linear) mp func)
  (iterate-over-neighbours-shape-linear mesh mp func))

(defmethod iterate-over-neighbours-shape (mesh (shape-func cl-mpm/shape-function:shape-function-bspline) mp func)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/shape-function::shape-function-bspline shape-func))
  (iterate-over-neighbours-shape-bspline mesh mp func))

(declaim (inline iterate-over-neighbours-shape-linear))
(defun iterate-over-neighbours-shape-linear (mesh mp func)
  "The simplest shape function implementation - linear"
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (pos-vec (cl-mpm/particle:mp-position mp))
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
           (pos-index (cl-mpm/mesh:position-to-index-floor mesh pos-vec)))
      (loop for dx from 0 to 1
            do (loop for dy from 0 to 1
                     do (let* ((id (mapcar #'+ pos-index (list dx dy 0))))
                          (when (cl-mpm/mesh:in-bounds mesh id)
                            (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                   (node (cl-mpm/mesh:get-node mesh id))
                                   (weights
                                     (mapcar (lambda (x) (cl-mpm/shape-function::shape-linear x h)) dist))
                                   (weight (reduce #'* weights))
                                   (grads (mapcar (lambda (d w) (* (cl-mpm/shape-function::shape-linear-dsvp d h) w))
                                                  dist (nreverse weights)))
                                   )
                              (when (< 0d0 weight)
                                (funcall func mesh mp node weight (append grads (list 0d0))
                                         1d0
                                         (list
                                          (* 0.5d0 (cl-mpm/shape-function::shape-linear-dsvp (nth 0 dist) h))
                                          (* 0.5d0 (cl-mpm/shape-function::shape-linear-dsvp (nth 1 dist) h))
                                          0d0
                                          )))))))))))

(declaim (inline iterate-over-neighbours-shape-linear))
(defun iterate-over-neighbours-shape-linear-3d (mesh mp func)
  "The simplest shape function implementation - linear"
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (pos-vec (cl-mpm/particle:mp-position mp))
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)  (tref pos-vec 2 0)))
           (pos-index (cl-mpm/mesh:position-to-index-floor mesh pos-vec)))
      (loop for dx from 0 to 1
            do (loop for dy from 0 to 1
                     do
                        (loop for dz from 0 to 1
                              do
                                 (let* ((id (mapcar #'+ pos-index (list dx dy dz))))
                                   (when (cl-mpm/mesh:in-bounds mesh id)
                                     (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                            (node (cl-mpm/mesh:get-node mesh id))
                                            (weights (mapcar (lambda (x) (cl-mpm/shape-function::shape-linear x h)) dist))
                                            (weight (reduce #'* weights))
                                            (lin-grads (mapcar (lambda (d) (cl-mpm/shape-function::shape-linear-dsvp d h)) dist))
                                            (grads (cl-mpm/shape-function::grads-3d weights lin-grads))
                                            )
                                       (when (< 0d0 weight)
                                         (funcall func mesh mp node weight grads
                                                  0d0 (list 0d0 0d0 0d0))))))))))))


(declaim (inline iterate-over-neighbours-point-linear)
         (ftype (function (cl-mpm/mesh::mesh magicl:matrix/double-float function) (values))
                iterate-over-neighbours-point-linear))

(defun iterate-over-neighbours-point-linear (mesh position func)
  "Iterate over neighbours of an arbitrary point - using FEM linear basis"
  (if (= (the fixnum (cl-mpm/mesh:mesh-nd mesh)) 2)
      (iterate-over-neighbours-point-linear-simd mesh position func)
      (iterate-over-neighbours-point-linear-3d mesh position func)))

(declaim (inline iterate-over-neighbours-point-linear-lisp)
         (ftype (function (cl-mpm/mesh::mesh magicl:matrix/double-float function) (values))
                iterate-over-neighbours-point-linear-lisp))

(defun iterate-over-neighbours-point-linear-lisp (mesh position func)
  "Iterating over basis functions in 2D"
  (declare (cl-mpm/mesh::mesh mesh))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (pos-vec position)
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
           (pos-index (cl-mpm/mesh:position-to-index-floor mesh pos-vec)))
      (declare (dynamic-extent pos pos-index pos-vec))
      (loop for dx from 0 to 1
            do (loop for dy from 0 to 1
                     do (let* ((id (mapcar #'+ pos-index (list dx dy))))
                          (declare (dynamic-extent id))
                          (when (cl-mpm/mesh:in-bounds mesh id)
                            (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                   (node (cl-mpm/mesh:get-node mesh id))
                                   (weights (mapcar (lambda (x) (cl-mpm/shape-function::shape-linear x h)) dist))
                                   (weight (reduce #'* weights))
                                   (grads (mapcar (lambda (d w) (* (the double-float (cl-mpm/shape-function::shape-linear-dsvp d h))
                                                                   (the double-float w)))
                                                  dist (nreverse weights))))
                              (declare (double-float weight)
                                       (dynamic-extent dist weights))
                              (when (< 0d0 weight)
                                (funcall func mesh node weight grads))))
                          ))))))
(defun iterate-over-neighbours-point-linear-3d (mesh position func)
  "Iterating over basis functions in 3D"
  (declare (cl-mpm/mesh::mesh mesh))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (pos-vec position)
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0) (tref pos-vec 2 0)))
           (pos-index (cl-mpm/mesh:position-to-index-floor mesh pos-vec)))
      (declare (dynamic-extent pos pos-index pos-vec))
      (loop for dx from 0 to 1
            do (loop for dy from 0 to 1
                     do (loop for dz from 0 to 1
                              do (let* ((id (mapcar #'+ pos-index (list dx dy dz))))
                                   (declare (dynamic-extent id))
                                   (when (cl-mpm/mesh:in-bounds mesh id)
                                     (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                            (node (cl-mpm/mesh:get-node mesh id))
                                            (weights (mapcar (lambda (x) (cl-mpm/shape-function::shape-linear x h)) dist))
                                            (weight (reduce #'* weights))
                                            (lin-grads (mapcar (lambda (d)
                                                                 (cl-mpm/shape-function::shape-linear-dsvp d h))
                                                               dist))
                                            (grads (cl-mpm/shape-function::grads-3d weights lin-grads)))
                                       (declare
                                        (double-float weight)
                                        (dynamic-extent dist weights))
                                       (when (< 0d0 weight)
                                         (funcall func mesh node weight grads)))))))))))

(defun iterate-over-neighbours-point-linear-2d (mesh position func)
  "Iterating over basis functions in 3D"
  (declare (cl-mpm/mesh::mesh mesh))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (pos-vec position)
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
           (pos-index (cl-mpm/mesh:position-to-index-floor mesh pos-vec)))
      (declare (dynamic-extent pos pos-index pos-vec))
      (loop for dx from 0 to 1
            do (loop for dy from 0 to 1
                     do (let* ((id (mapcar #'+ pos-index (list dx dy 0))))
                          (declare (dynamic-extent id))
                          (when (cl-mpm/mesh:in-bounds mesh id)
                            (destructuring-bind (dist-x dist-y) (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id))
                              (let* ((node (cl-mpm/mesh:get-node mesh id))
                                     (weight-x (cl-mpm/shape-function::shape-linear dist-x h))
                                     (weight-y (cl-mpm/shape-function::shape-linear dist-y h))
                                     (weight (* weight-x weight-y))
                                     (grad-x (* weight-y (cl-mpm/shape-function::shape-linear-dsvp dist-x h)))
                                     (grad-y (* weight-x (cl-mpm/shape-function::shape-linear-dsvp dist-y h)))
                                    )
                                (when (< 0d0 weight)
                                  (funcall func mesh node weight (list grad-x grad-y 0d0)))
                                ))
                            ;; (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                            ;;        (node (cl-mpm/mesh:get-node mesh id))
                            ;;        (weights (mapcar (lambda (x) (cl-mpm/shape-function::shape-linear x h)) dist))
                            ;;        (weight (reduce #'* weights))
                            ;;        (lin-grads (mapcar (lambda (d)
                            ;;                             (cl-mpm/shape-function::shape-linear-dsvp d h))
                            ;;                           dist))
                            ;;        (grads (cl-mpm/shape-function::grads-3d weights lin-grads)))
                            ;;   (declare
                            ;;    (double-float weight)
                            ;;    (dynamic-extent dist weights))
                            ;;   (when (< 0d0 weight)
                            ;;     (funcall func mesh node weight grads)))
                            )))))))



(declaim (inline iterate-over-neighbours-point-linear-simd)
         (ftype (function (cl-mpm/mesh::mesh magicl:matrix/double-float function) (values))
                iterate-over-neighbours-point-linear-simd)
         )
(defun iterate-over-neighbours-point-linear-simd (mesh position func)
  "A fast implemenntation of 2D linear basis function iteration"
  (declare (cl-mpm/mesh::mesh mesh)
           (magicl:matrix/double-float position)
           (function func))
  (labels ((simd-abs (vec)
           (sb-simd-avx:f64.2-and vec (sb-simd-avx:f64.2-not -0.0d0)))
         ;; (in-bounds-simd (pos)
         ;;   t)
         (linear-weight-simd (dist h)
           (declare (double-float h))
           ;;Add an abs
           (sb-simd-avx:f64.2- 1d0 (sb-simd-avx:f64.2/ (simd-abs dist) h)))
         (linear-grads-simd (dist h)
           (declare (double-float h))
           (sb-simd-avx:f64.2-if (sb-simd-avx:f64.2> dist 0d0) (/ 1d0 h) (/ -1d0 h))
           dist)
         (in-bounds-simd (mesh dist)
           t)
         )
    (progn
      (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
             (pos-vec (sb-simd-avx:f64.2-aref (magicl::matrix/double-float-storage position) 0))
             (pos-index (sb-simd-avx:f64.2-floor
                         (sb-simd-avx:f64.2/ pos-vec h))))
        (declare (sb-simd-avx:f64.2 pos-vec)
                 (double-float h))
        (loop for dx fixnum from 0 to 1
              do (loop for dy fixnum from 0 to 1
                       do (let* ((id-vec
                                   (sb-simd-avx:f64.2+ pos-index (sb-simd-avx:make-f64.2 dx dy)))
                                 ;; (id (mapcar (lambda (x) (truncate x))
                                 ;;             (multiple-value-list (sb-simd-avx:f64.2-values id-vec))))
                                 (id
                                   (append
                                    (mapcar (lambda (x) (truncate (the double-float x)))
                                            (multiple-value-list (sb-simd-avx:f64.2-values id-vec)))
                                    '(0)))
                                 )
                            (declare (dynamic-extent id))
                            (when (cl-mpm/mesh:in-bounds mesh id)
                                ;(in-bounds-simd mesh id-vec)
                              (let* ((dist (sb-simd-avx:f64.2-
                                            pos-vec
                                            (sb-simd-avx:f64.2* id-vec h)))
                                     (node (cl-mpm/mesh:get-node mesh id))
                                     (weights (linear-weight-simd dist h))
                                     (weight (sb-simd-avx::f64.2-horizontal* weights))
                                     (grads-vec (sb-simd-avx:f64.2*
                                                 (linear-grads-simd dist h)
                                                 (sb-simd-avx:f64.2-shuffle weights weights 1)))
                                     (grads (append (multiple-value-list (sb-simd-avx:f64.2-values grads-vec))
                                                    (list 0d0)
                                                    ))
                                     )
                                (declare (double-float weight))
                                (when (< 0d0 weight)
                                  (funcall func mesh node weight grads))))
                            )))))))


(declaim ;(inline iterate-over-neighbours-shape-gimp-simd)
         (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle function) (values))
                iterate-over-neighbours-shape-gimp-simd))
(defun iterate-over-neighbours-shape-gimp-simd (mesh mp func)
  "Iterate over a gimp domains neighbours in 2D using SIMD constructs"
  (declare (type cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp)
           (function func)
           (optimize (speed 3) (safety 0) (debug 0)))
  (progn
    ;; (let ((pos-vec (sb-mop:standard-instance-access mp 0))
    ;;       (d0 (sb-mop:standard-instance-access mp 1))
    ;;       ))
    (with-accessors ((pos-vec cl-mpm/particle:mp-position)
                     (d0 cl-mpm/particle::mp-domain-size))
        mp
      (let ((h (the double-float (cl-mpm/mesh:mesh-resolution mesh))))
        (let* ((pa (sb-simd-avx:f64.2-aref (magicl::matrix/double-float-storage pos-vec) 0))
               (d0a (sb-simd-avx:f64.2-aref (magicl::matrix/double-float-storage d0) 0))
               (ia (sb-simd-avx:f64.2-round (sb-simd-avx:f64.2*
                                             pa (/ 1d0 h))))
               (ca (sb-simd-avx:f64.2- pa (sb-simd-avx:f64.2* ia h)))
               ;;Lower bound of domain
               (dfa (sb-simd-avx:f64.2-floor
                     (sb-simd-avx:f64.2/
                      (sb-simd-avx:f64.2-
                       ca
                       (sb-simd-avx:f64.2*
                        d0a 0.5d0))
                      h)))
               ;;Upper bound of domain
               (dca (sb-simd-avx:f64.2-ceiling
                     (sb-simd-avx:f64.2/
                      (sb-simd-avx:f64.2+
                       ca
                       (sb-simd-avx:f64.2*
                        d0a 0.5d0)
                       )
                      h))))
          (multiple-value-bind (ix iy) (sb-simd-avx:f64.2-values ia)
            (multiple-value-bind (dox doy) (sb-simd-avx:f64.2-values d0a)
              (multiple-value-bind (cx cy) (sb-simd-avx:f64.2-values ca)
                (multiple-value-bind (dxf dyf) (sb-simd-avx:f64.2-values dfa)
                  (multiple-value-bind (dxc dyc) (sb-simd-avx:f64.2-values dca)
                    (declare (type double-float h cx cy dox doy)
                                        ;(type integer dxf dxc dyf dyc ix iy)
                             )
                    (loop for dx fixnum from (the fixnum (truncate dxf)) to (the fixnum (truncate dxc))
                     ;; loop for dx fixnum from -2 to 2
                           do (loop for dy fixnum from (the fixnum (truncate dyf)) to (the fixnum (truncate dyc))
                               ;; loop for dy fixnum from -2 to 2
                                   do
                                      (let* ((id (list (+ (the fixnum (truncate ix)) dx)
                                                       (+ (the fixnum (truncate iy)) dy)
                                                       0
                                                       )))
                                        (declare (dynamic-extent id))
                                        (when (cl-mpm/mesh:in-bounds mesh id)
                                          (let* ((dist (sb-simd-avx:f64.2-
                                                        ca
                                                        (sb-simd-avx:f64.2*
                                                         (sb-simd-avx:make-f64.2 dx dy)
                                                         h))))
                                            (multiple-value-bind (distx disty) (sb-simd-avx:f64.2-values dist)
                                              (declare (type double-float distx disty))
                                              (let* ((weightsx (the double-float (cl-mpm/shape-function::shape-gimp distx (* 0.5d0 dox) h)))
                                                     (weightsy (the double-float (cl-mpm/shape-function::shape-gimp disty (* 0.5d0 doy) h)))
                                                     (weight (* weightsx weightsy))
                                                     (weights-fbar-x (the double-float (cl-mpm/shape-function::shape-gimp-fbar distx (* 0.5d0 dox) h)))
                                                     (weights-fbar-y (the double-float (cl-mpm/shape-function::shape-gimp-fbar disty (* 0.5d0 doy) h)))
                                                     (weight-fbar (* weights-fbar-x weights-fbar-y))
                                                     ;; #+cl-mpm-fbar
                                                     )
                                                (when (or (< 0d0 weight) (< 0d0 weight-fbar))
                                                  (let* ((node (cl-mpm/mesh:get-node mesh id))
                                                         (gradx (* (cl-mpm/shape-function::shape-gimp-dsvp distx (* 0.5d0 dox) h)
                                                                   (the double-float weightsy)))
                                                         (grady (* (cl-mpm/shape-function::shape-gimp-dsvp disty (* 0.5d0 doy) h)
                                                                   (the double-float weightsx)))

                                                         (gradz 0d0)
                                                         (grads-fbar
                                                           (list (* weights-fbar-y (cl-mpm/shape-function::shape-gimp-dsvp distx (* 0.5d0 dox) h))
                                                                 (* weights-fbar-x (cl-mpm/shape-function::shape-gimp-dsvp disty (* 0.5d0 doy) h))
                                                                 0d0)))
                                                    (declare (double-float gradx grady))
                                                    (funcall func mesh mp node
                                                             weight (list gradx grady gradz)
                                                             weight-fbar grads-fbar))))))))))))))))))))
;;This is more consise but half as fast
(defun iterate-over-neighbours-shape-gimp (mesh mp func)
  "Iterate over a gimp domains neighbours in 3D using simple lisp constructs"
  (declare (type cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp)
           (function func))
  (progn
    (with-accessors ((pos-vec cl-mpm/particle:mp-position)
                     (d0 cl-mpm/particle::mp-domain-size))
        mp
      (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
             (pos (list (varef pos-vec 0) (varef pos-vec 1) (varef pos-vec 2)))
             (pos-index (cl-mpm/mesh:position-to-index-round mesh pos-vec))
             )
        (declare (dynamic-extent pos pos-index))
        (declare (type double-float h)
                 (type list pos pos-index))
        (loop for dx from -2 to 2
              do (loop for dy from -2 to 2
                       do (loop for dz from -2 to 2
                                do (let* ((id (mapcar #'+ pos-index (list dx dy dz))))
                                     (declare (dynamic-extent id))
                                     (when (cl-mpm/mesh:in-bounds mesh id)
                                       (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                              (domain (loop for x across (magicl::matrix/double-float-storage d0) collect (* 0.5d0 (the double-float x))))
                                              (weights (mapcar (lambda (x l)
                                                                 (cl-mpm/shape-function::shape-gimp x l h))
                                                               dist domain))
                                              (weight (reduce #'* weights)))
                                         (declare (dynamic-extent domain))
                                         (declare (type double-float h)
                                                  (type list domain))
                                         (when (< 0d0 weight)
                                           (let* ((node (cl-mpm/mesh:get-node mesh id))
                                                  (lin-grads (mapcar (lambda (d l)
                                                                       (cl-mpm/shape-function::shape-gimp-dsvp d l h))
                                                                     dist domain))
                                                  (grads (cl-mpm/shape-function::grads-3d weights lin-grads))
                                                  )
                                             (funcall func mesh mp node weight grads 0d0 (list 0d0 0d0 0d0))))))))))))))

(defun iterate-over-neighbours-shape-gimp-2d (mesh mp func)
  "Iterate over a gimp domains neighbours in 2D unrolled, this version is more performant but worse than SIMD"
  (declare (type cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp)
           (function func)
           (optimize (speed 3) (safety 0) (debug 0))
           )
  (progn
    (with-accessors ((pos-vec cl-mpm/particle:mp-position)
                     (d0 cl-mpm/particle::mp-domain-size))
        mp
      (let ((h (the double-float (cl-mpm/mesh:mesh-resolution mesh))))
        (flet ((center-id (x)
                 (round x h)))
          (let* ((pa (cl-mpm/utils:fast-storage pos-vec))
                 (da (cl-mpm/utils:fast-storage d0))
                 (px (the double-float (aref pa 0)))
                 (py (the double-float (aref pa 1)))
                 (ix (the fixnum (truncate (center-id px))))
                 (iy (the fixnum (truncate (center-id py))))
                 (cx (- px (* h ix)))
                 (cy (- py (* h iy)))
                 (dox (* 0.5d0 (the double-float (aref da 0))))
                 (doy (* 0.5d0 (the double-float (aref da 1))))
                 (dxf (the fixnum (truncate (ffloor   (- cx dox) h))))
                 (dxc (the fixnum (truncate (fceiling (+ cx dox) h))))
                 (dyf (the fixnum (truncate (ffloor   (- cy doy) h))))
                 (dyc (the fixnum (truncate (fceiling (+ cy doy) h))))
                 )
            (declare ((simple-array double-float (3)) pa))
            ;; (declare (dynamic-extent pa))
            (declare (type double-float h cx cy dox doy px py )
                     (type integer dxf dxc dyf dyc ix iy )
                     )
            (loop for dx from -2 to 2;dxf to dxc
                  do (loop for dy from -2 to 2;dyf to dyc
                           do
                              (let* ((id (list (+ ix dx)
                                               (+ iy dy)
                                               0)))
                                (declare (dynamic-extent id))
                                (when (cl-mpm/mesh:in-bounds mesh id)
                                  (let* ((distx (- cx (* h dx)))
                                         (disty (- cy (* h dy)))
                                         (weightsx (cl-mpm/shape-function::shape-gimp distx dox h))
                                         (weightsy (cl-mpm/shape-function::shape-gimp disty doy h))
                                         (weights-fbar-x (the double-float (cl-mpm/shape-function::shape-gimp-fbar distx dox h)))
                                         (weights-fbar-y (the double-float (cl-mpm/shape-function::shape-gimp-fbar disty doy h)))
                                         (weight-fbar (* weights-fbar-x weights-fbar-y))
                                         (weight (* weightsx weightsy)))
                                    (declare ;(type double-float h)
                                     (double-float weight weightsx weightsy distx disty))
                                    (when (or (> weight 0d0) (> weight-fbar 0d0))
                                      (let* ((node (cl-mpm/mesh:get-node mesh id))
                                             (gradx
                                               (* (cl-mpm/shape-function::shape-gimp-dsvp distx dox h)
                                                       weightsy))
                                             (grady
                                               (* (cl-mpm/shape-function::shape-gimp-dsvp disty doy h)
                                                       weightsx))
                                             (grads-fbar
                                               (list (* weights-fbar-y (cl-mpm/shape-function::shape-gimp-dsvp distx dox h))
                                                     (* weights-fbar-x (cl-mpm/shape-function::shape-gimp-dsvp disty doy h))
                                                     0d0))
                                             )
                                        (declare (double-float gradx grady))
                                        (funcall func mesh mp node
                                                 weight (list gradx grady 0d0)
                                                 weight-fbar
                                                 grads-fbar))))))))))))))

(declaim ;(inline iterate-over-neighbours-shape-gimp)
         (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle function) (values))
                iterate-over-neighbours-shape-gimp-3d))
(defun iterate-over-neighbours-shape-gimp-3d (mesh mp func)
  "Iterate over gimp neighbours in 3D, unrolled for speed"
  (declare (type cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp)
           (function func)
           (optimize (speed 3) (safety 0) (debug 0))
           )
  (progn
    (with-accessors ((pos-vec cl-mpm/particle:mp-position)
                     (d0 cl-mpm/particle::mp-domain-size))
        mp
      (let ((h (the double-float (cl-mpm/mesh:mesh-resolution mesh))))
        (flet ((center-id (x)
                 (round x h))
               )
          (let* ((pa (magicl::matrix/double-float-storage pos-vec))
                 (da (magicl::matrix/double-float-storage d0))
                 (px (the double-float (aref pa 0)))
                 (py (the double-float (aref pa 1)))
                 (pz (the double-float (aref pa 2)))
                 (ix (the fixnum (truncate (center-id px))))
                 (iy (the fixnum (truncate (center-id py))))
                 (iz (the fixnum (truncate (center-id pz))))
                 (cx (- px (* h ix)))
                 (cy (- py (* h iy)))
                 (cz (- pz (* h iz)))
                 (dox (* 0.5d0 (the double-float (aref da 0))))
                 (doy (* 0.5d0 (the double-float (aref da 1))))
                 (doz (* 0.5d0 (the double-float (aref da 2))))
                 (dxf (the fixnum (truncate (ffloor   (- cx dox) h))))
                 (dxc (the fixnum (truncate (fceiling (+ cx dox) h))))
                 (dyf (the fixnum (truncate (ffloor   (- cy doy) h))))
                 (dyc (the fixnum (truncate (fceiling (+ cy doy) h))))
                 (dzf (the fixnum (truncate (ffloor   (- cz doz) h))))
                 (dzc (the fixnum (truncate (fceiling (+ cz doz) h))))
                 )
            (declare ((simple-array double-float *) pa))
            ;; (declare (dynamic-extent pa))
            (declare (type double-float h cx cy dox doy px py pz doz cz)
                     (type integer dxf dxc dyf dyc ix iy dzf dzc iz)
                     )
            (loop for dx from dxf to dxc
                  do (loop for dy from dyf to dyc
                           do
                              (loop for dz from dzf to dzc
                                    do
                              (let* ((id (list (the fixnum (+ ix dx))
                                               (the fixnum (+ iy dy))
                                               (the fixnum (+ iz dz))
                                               )))
                                (declare (dynamic-extent id))
                                (when (cl-mpm/mesh:in-bounds mesh id)
                                  (let* ((distx (- cx (* h dx)))
                                         (disty (- cy (* h dy)))
                                         (distz (- cz (* h dz)))
                                         (weightsx (cl-mpm/shape-function::shape-gimp-fast distx dox h))
                                         (weightsy (cl-mpm/shape-function::shape-gimp-fast disty doy h))
                                         (weightsz (cl-mpm/shape-function::shape-gimp-fast distz doz h))
                                         (weight (* weightsx weightsy weightsz))
                                         (weights-fbar-x (the double-float (cl-mpm/shape-function::shape-gimp-fbar distx dox h)))
                                         (weights-fbar-y (the double-float (cl-mpm/shape-function::shape-gimp-fbar disty doy h)))
                                         (weights-fbar-z (the double-float (cl-mpm/shape-function::shape-gimp-fbar distz doz h)))
                                         (weight-fbar (* weights-fbar-x weights-fbar-y weights-fbar-z))
                                         ;; (weight-fbar 0d0)
                                         )
                                    (declare ;(type double-float h)
                                     (double-float weight weightsx weightsy weightsz distx disty distz)
                                     )
                                    (when (< 0d0 weight)
                                      (let* ((node (cl-mpm/mesh:get-node mesh id))
                                             (gradx (* (cl-mpm/shape-function::shape-gimp-dsvp distx dox h)
                                                       weightsy weightsz
                                                       ))
                                             (grady (* (cl-mpm/shape-function::shape-gimp-dsvp disty doy h)
                                                       weightsx weightsz
                                                       ))
                                             (gradz (* (cl-mpm/shape-function::shape-gimp-dsvp distz doz h)
                                                       weightsx weightsy
                                                       ))
                                             (grads-fbar
                                               (list (* weights-fbar-y weights-fbar-z
                                                        (cl-mpm/shape-function::shape-gimp-dsvp distx dox h))
                                                     (* weights-fbar-x weights-fbar-z
                                                        (cl-mpm/shape-function::shape-gimp-dsvp disty doy h))
                                                     (* weights-fbar-x weights-fbar-y
                                                        (cl-mpm/shape-function::shape-gimp-dsvp distz doz h)))))
                                        (declare (double-float gradx grady gradz))
                                        (funcall func mesh mp node
                                                 weight (list gradx grady gradz)
                                                 weight-fbar grads-fbar)))))))))))))))

(defun iterate-over-neighbours-shape-gimp-simd-3d (mesh mp func)
  "Iterate over gimp neighbours in 3D, unrolled for speed"
  (declare (type cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp)
           (function func)
           (optimize (speed 3) (safety 0) (debug 0))
           )
  (progn
    (with-accessors ((pos-vec cl-mpm/particle:mp-position)
                     (d0 cl-mpm/particle::mp-domain-size))
        mp
      (let ((h (the double-float (cl-mpm/mesh:mesh-resolution mesh))))
        (flet ((center-id (x)
                 (round x h))
               )
          (let* ((pa (magicl::matrix/double-float-storage pos-vec))
                 (da (magicl::matrix/double-float-storage d0))

                 (ia (let ((arr (make-array 3 :initial-element 0d0 :element-type 'double-float)))
                       (loop for i from 0 to 2 do (setf (aref arr i) (fround (aref arr i) h))) arr))
                 (ca (let ((arr (make-array 3 :initial-element 0d0 :element-type 'double-float)))
                       ;; (loop for i from 0 to 2
                       ;;       for pv across pa
                       ;;       for iv across ia
                       ;;       do (setf (aref arr i) (- pv (* h iv))))
                       (setf
                        (sb-simd-avx:f64.2-aref
                         arr 0)
                        (sb-simd-avx:f64.2-
                         (sb-simd-avx:f64.2-aref
                          pa 0)
                         (sb-simd-avx:f64.2*
                          (sb-simd-avx:f64.2-aref
                           ia 0)
                          h))
                        (aref arr 2)
                        (- (aref pa 2) (* h (aref ia 2))))
                       arr))

                 (doa (let ((arr (make-array 3 :initial-element 0d0 :element-type 'double-float)))
                        (setf
                         (sb-simd-avx:f64.2-aref
                          arr 0)
                         (sb-simd-avx:f64.2*
                          (sb-simd-avx:f64.2-aref
                           da 0)
                          0.5d0)
                         (aref arr 2)
                         (* 0.5d0 (aref da 2)))
                       ;; (loop for i from 0 to 2
                       ;;       for dv across da
                       ;;       do (setf (aref arr i) (* 0.5d0 dv)))
                        arr))
                 (dfa (let ((arr (make-array 3 :initial-element 0d0 :element-type 'double-float)))
                       (loop for i from 0 to 2
                             for cv across ca
                             for dov across doa
                             do (setf (aref arr i) (ffloor (- cv dov) h))) arr))
                 (dca (let ((arr (make-array 3 :initial-element 0d0 :element-type 'double-float)))
                        (loop for i from 0 to 2
                              for cv across ca
                              for dov across doa
                              do (setf (aref arr i) (fceiling (+ cv dov) h))) arr))
                 (dfia (let ((arr (make-array 3 :initial-element 0 :element-type 'fixnum)))
                         (loop for v across dfa
                               for i from 0
                               do (setf (aref arr i) (truncate v))) arr))
                 (dcia (let ((arr (make-array 3 :initial-element 0 :element-type 'fixnum)))
                         (loop for v across dca
                               for i from 0
                               do (setf (aref arr i) (truncate v))) arr))

                 (iia (let ((arr (make-array 3 :initial-element 0 :element-type 'fixnum)))
                       (loop for i from 0 to 2 do (setf (aref arr i) (truncate (aref ia i)))) arr))
                 )
            (declare ((simple-array double-float (3)) pa da ia ca doa dfa dca))
            (declare (dynamic-extent pa da ia ca doa dfa dca dfia dcia iia))
            (loop for dx fixnum from (aref dfia 0) to (aref dcia 0)
                  do (loop for dy fixnum from (aref dfia 1) to (aref dcia 1)
                           do
                              (loop for dz fixnum from (aref dfia 2) to (aref dcia 2)
                                    do
                                       (let* ((id (list (the fixnum (+ (the fixnum (aref iia 0)) dx))
                                                        (the fixnum (+ (the fixnum (aref iia 1)) dy))
                                                        (the fixnum (+ (the fixnum (aref iia 2)) dz))
                                                        )))
                                         (declare (dynamic-extent id))
                                         (when (cl-mpm/mesh:in-bounds mesh id)
                                           (let* ((dista
                                                    (let ((arr (make-array 3 :initial-element 0d0 :element-type 'double-float)))
                                                      (declare ((simple-array double-float (3))))
                                                           (setf
                                                            (sb-simd-avx:f64.2-aref
                                                             arr 0)
                                                            (sb-simd-avx:f64.2-
                                                             (sb-simd-avx:f64.2-aref
                                                              ca 0)
                                                             (sb-simd-avx:f64.2*
                                                              (sb-simd-avx:make-f64.2
                                                               dx
                                                               dy)
                                                              h))
                                                            ;; (aref arr 0) (- (aref ca 0) (* h dx))
                                                            ;; (aref arr 1) (- (aref ca 1) (* h dy))
                                                            (aref arr 2) (- (aref ca 2) (* h dz)))
                                                            arr))
                                                  (weightsx (cl-mpm/shape-function::shape-gimp-fast
                                                             (aref dista 0) (aref doa 0) h))
                                                  (weightsy (cl-mpm/shape-function::shape-gimp-fast
                                                             (aref dista 1) (aref doa 1) h))
                                                  (weightsz (cl-mpm/shape-function::shape-gimp-fast
                                                             (aref dista 2) (aref doa 2) h))
                                                  (weight (* weightsx weightsy weightsz))
                                                  )
                                             (declare ;(type double-float h)
                                              (double-float weight weightsx weightsy weightsz)
                                              )
                                             (when (< 0d0 weight)
                                               (let* ((node (cl-mpm/mesh:get-node mesh id))
                                                      (gradx (* (cl-mpm/shape-function::shape-gimp-dsvp (aref dista 0) (aref doa 0) h)
                                                                weightsy weightsz
                                                                ))
                                                      (grady (* (cl-mpm/shape-function::shape-gimp-dsvp (aref dista 1) (aref doa 1) h)
                                                                weightsx weightsz
                                                                ))
                                                      (gradz (* (cl-mpm/shape-function::shape-gimp-dsvp (aref dista 2) (aref doa 2) h)
                                                                weightsx weightsy
                                                                ))
                                                      )
                                                 (declare (double-float gradx grady gradz))
                                                 (funcall func mesh mp node
                                                          weight (list gradx grady gradz)
                                                          0d0 (list 0d0 0d0 0d0))))))))))))))))

(defun make-knot-list (mesh pos)
  "Function for maybe being able to make a bspline knot list"
  (list
   (let ((dy 0))
     (loop for dx from -4 to 4
           collect (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos (list dx dy)))))
   (let ((dx 0))
     (loop for dy from -4 to 4
           collect (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos (list dx dy)))))))
(declaim (inline iterate-over-neighbours-shape-bspline))
(defun iterate-over-neighbours-shape-bspline (mesh mp func)
  "Iterate over 2D bsplines, however boundary adjustment doesn't work"
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (pos-vec (cl-mpm/particle:mp-position mp))
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0) 0))
           ;; (pos-index (magicl:scale pos-vec (/ 1 h)))
           (pos-index (cl-mpm/mesh:position-to-index-round mesh pos-vec))
           (border (not (and (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list 2 2 0)))
                             (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list -2 -2 0)))
                             ;; (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list 2 -2 0)))
                             ;; (cl-mpm/mesh:in-bounds mesh (mapcar #'+ pos-index (list -2 2 0)))
                             )))
           )
      (if border
          (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
                 (pos-vec (cl-mpm/particle:mp-position mp))
                 (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0)))
                 (pos-index (cl-mpm/mesh:position-to-index-floor mesh pos-vec)))
            (loop for dx from 0 to 1
                  do (loop for dy from 0 to 1
                           do (let* ((id (mapcar #'+ pos-index (list dx dy 0))))
                                (when (cl-mpm/mesh:in-bounds mesh id)
                                  (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                         (node (cl-mpm/mesh:get-node mesh id))
                                         (weights (mapcar (lambda (x) (cl-mpm/shape-function::shape-linear x h)) dist))
                                         (weight (reduce #'* weights))
                                         (grads (mapcar (lambda (d w) (* (cl-mpm/shape-function::shape-linear-dsvp d h) w))
                                                        dist (nreverse weights))))
                                    (when (< 0d0 weight)
                                      (funcall func mesh mp node weight (append grads (list 0d0)) 0d0 (list 0d0 0d0 0d0)))))))))

          (loop for dx from -1 to 1
                do (loop for dy from -1 to 1
                         do (let* ((id (mapcar #'+ pos-index (list dx dy 0))))
                              (when (cl-mpm/mesh:in-bounds mesh id)
                                (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                       (node (cl-mpm/mesh:get-node mesh id))
                                       (weights (mapcar (lambda (x) (cl-mpm/shape-function:shape-bspline x h)) dist))
                                       (weight (reduce #'* weights))
                                       (gradx (* (cl-mpm/shape-function::shape-bspline-dsvp (nth 0 dist) h)
                                                 (the double-float (nth 1 weights))))
                                       (grady (* (cl-mpm/shape-function::shape-bspline-dsvp (nth 1 dist) h)
                                                 (the double-float (nth 0 weights))))
                                       (gradz 0d0)
                                       )
                                  (funcall func mesh mp node weight (list gradx grady gradz) 0d0 (list 0d0 0d0 0d0)))))))))))

;; (declaim (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle function) (values)) iterate-over-corners-2d))
(defun iterate-over-corners-2d-step (mesh mp step func)
  (declare (cl-mpm/particle::particle mp)
           (function func))
  ;; (array-operations/utilities:nested-loop (x y) '(2 2))
  (flet ((integrate (x y)
           (let ((domain (cl-mpm/particle::mp-domain-size mp))
                 (position (cl-mpm/particle:mp-position mp))
                 (corner (cl-mpm/utils:vector-zeros)))
             (declare (double-float x y))
             (cl-mpm/fastmaths::fast-.+-vector
              position
              (cl-mpm/fastmaths:fast-scale!
               (cl-mpm/fastmaths:fast-.*
                (vector-from-list
                 (list
                  x
                  y
                  0d0))
                domain
                ) 0.5d0) corner)
             (funcall func corner)))
         )
    (loop for v from -1d0 to 1d0 by step
          do
             (progn
               (integrate v -1d0)
               (integrate v  1d0)
               (integrate -1d0 v)
               (integrate  1d0 v)))))

(declaim (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle function) (values)) iterate-over-corners-2d))
(defun iterate-over-corners-2d (mesh mp func)
  (declare (cl-mpm/particle::particle mp)
           (function func))
  ;; (array-operations/utilities:nested-loop (x y) '(2 2))
  (loop for x from -1d0 to 1d0 by 2d0
        do
           (loop for y from -1d0 to 1d0 by 2d0
                 do
                    (let ((domain (cl-mpm/particle::mp-domain-size mp))
                          (position (cl-mpm/particle:mp-position mp))
                          (corner (cl-mpm/utils:vector-zeros)))
                      (declare (double-float x y))
                      (cl-mpm/fastmaths::fast-.+-vector
                       position
                       (cl-mpm/fastmaths:fast-scale!
                        (cl-mpm/fastmaths:fast-.*
                         (vector-from-list
                          (list
                           x
                           y
                           0d0))
                         domain
                         ) 0.5d0) corner)
                      (funcall func corner)))))

(declaim (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle function) (values)) iterate-over-corners-normal-2d))
(defun iterate-over-corners-normal-2d (mesh mp func)
  (declare (cl-mpm/particle::particle mp)
           (function func))
  (loop for x from -1d0 to 1d0 by 2d0
        do (loop for y from -1d0 to 1d0 by 2d0
                 do (let ((domain (cl-mpm/particle::mp-domain-size mp))
                          (position (cl-mpm/particle:mp-position mp))
                          (corner (cl-mpm/utils:vector-zeros))
                          )
                      (declare (double-float x y))
                      (cl-mpm/fastmaths::fast-.+-vector
                       position
                       (cl-mpm/fastmaths:fast-scale!
                        (cl-mpm/fastmaths:fast-.*
                         (vector-from-list
                          (list
                           x
                           y
                           0d0))
                         domain
                         ) 0.5d0) corner)
                      (funcall func corner 
                               (cl-mpm/utils:vector-from-list (list x y 0d0)))))))

(defun iterate-over-corners-normal-3d (mesh mp func)
  (declare (cl-mpm/particle::particle mp)
           (function func))
  (loop for x from -1d0 to 1d0 by 2d0
        do (loop for y from -1d0 to 1d0 by 2d0
                 do (loop for z from -1d0 to 1d0 by 2d0
                          do (let ((domain (cl-mpm/particle::mp-domain-size mp))
                                   (position (cl-mpm/particle:mp-position mp))
                                   (corner (cl-mpm/utils:vector-zeros)))
                               (declare (double-float x y z))
                               (cl-mpm/fastmaths::fast-.+-vector
                                position
                                (magicl:scale!
                                 (magicl:.*
                                  (vector-from-list (list x y z))
                                  domain
                                  ) 0.5d0) corner)
                               (funcall func corner (cl-mpm/utils:vector-from-list (list x y z))))))))

(declaim (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle function) (values)) iterate-over-corners-3d))
(defun iterate-over-corners-3d (mesh mp func)
  (declare (cl-mpm/particle::particle mp)
           (function func))
  (loop for x from -1d0 to 1d0 by 2d0
        do
           (loop for y from -1d0 to 1d0 by 2d0
                 do
                    (loop for z from -1d0 to 1d0 by 2d0
                          do
                             (let ((domain (cl-mpm/particle::mp-domain-size mp))
                                   (position (cl-mpm/particle:mp-position mp))
                                   (corner (cl-mpm/utils:vector-zeros)))
                               (cl-mpm/fastmaths::fast-.+-vector
                                position
                                (magicl:scale!
                                 (magicl:.*
                                  (vector-from-list (list x y z))
                                  domain
                                  ) 0.5d0) corner)
                               (funcall func corner))))))

(declaim (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle function) (values)) iterate-over-corners))
(defun iterate-over-corners (mesh mp func)
  "Iterates a function (corner) over a single MPs corners"
  (if (= (the fixnum (cl-mpm/mesh:mesh-nd mesh)) 2)
      (iterate-over-corners-2d mesh mp func)
      (iterate-over-corners-3d mesh mp func)))

(defun iterate-over-corners-step (mesh mp step func)
  (if (= (the fixnum (cl-mpm/mesh:mesh-nd mesh)) 2)
      (iterate-over-corners-2d-step mesh mp step func)
      (error "Not implemented")))

(defun iterate-over-midpoints-normal-2d (mesh mp func)
  (declare (cl-mpm/particle::particle mp)
           (function func))
  (flet ((update (x y)
           (let ((domain (cl-mpm/particle::mp-domain-size mp))
                 (position (cl-mpm/particle:mp-position mp))
                 (corner (cl-mpm/utils:vector-zeros))
                 )
             (declare (double-float x y))
             (cl-mpm/fastmaths::fast-.+-vector
              position
              (cl-mpm/fastmaths:fast-scale!
               (cl-mpm/fastmaths:fast-.*
                (vector-from-list
                 (list
                  x
                  y
                  0d0))
                domain
                ) 0.5d0) corner)
             (funcall func corner 
                      (cl-mpm/utils:vector-from-list (list x y 0d0))))))
    (loop for x from -1d0 to 1d0 by 2d0
          do (update x 0d0))
    (loop for y from -1d0 to 1d0 by 2d0
          do (update 0d0 y))))

(defun iterate-over-midpoints-2d (mesh mp func)
  (declare (cl-mpm/particle::particle mp)
           (function func))
  (flet ((update (x y)
           (let ((domain (cl-mpm/particle::mp-domain-size mp))
                 (position (cl-mpm/particle:mp-position mp))
                 (corner (cl-mpm/utils:vector-zeros)))
             (declare (double-float x y))
             (cl-mpm/fastmaths::fast-.+-vector
              position
              (cl-mpm/fastmaths:fast-scale!
               (cl-mpm/fastmaths:fast-.*
                (vector-from-list
                 (list
                  x
                  y
                  0d0))
                domain
                ) 0.5d0) corner)
             (funcall func corner))))
    (loop for x from -1d0 to 1d0 by 2d0
          do (update x 0d0))
    (loop for y from -1d0 to 1d0 by 2d0
          do (update 0d0 y))))
(defun iterate-over-midpoints-3d (mesh mp func)
  (declare (cl-mpm/particle::particle mp)
           (function func))
  (flet ((update (x y)
           (let ((domain (cl-mpm/particle::mp-domain-size mp))
                 (position (cl-mpm/particle:mp-position mp))
                 (corner (cl-mpm/utils:vector-zeros)))
             (declare (double-float x y))
             (cl-mpm/fastmaths::fast-.+-vector
              position
              (cl-mpm/fastmaths:fast-scale!
               (cl-mpm/fastmaths:fast-.*
                (vector-from-list
                 (list
                  x
                  y
                  0d0))
                domain
                ) 0.5d0) corner)
             (funcall func corner))))
    (loop for x from -1d0 to 1d0 by 2d0
          do (update x 0d0))
    (loop for y from -1d0 to 1d0 by 2d0
          do (update 0d0 y))
    (loop for z from -1d0 to 1d0 by 2d0
          do (update z 0d0))))

(defun weight-local-pos (x y z nd dx dy dz)
  (let ((dx (- (* dx 2) 1))
        (dy (- (* dy 2) 1))
        (dz (- (* dz 2) 1)))
    (case nd
      (2 (* 0.25d0  (+ 1d0 (* dx x)) (+ 1d0 (* dy y))))
      (3 (* 0.125d0 (+ 1d0 (* dx x)) (+ 1d0 (* dy y)) (+ 1d0 (* dz z)))))))

(defun weights-local-pos (xi nd dxi)
  (let ((dxi (- (* dxi 2) 1))
       )
    (case nd
      (2 (* 0.5d0 (+ 1d0 (* dxi xi))))
      (3 (* 0.5d0 (+ 1d0 (* dxi xi)))))))

(defun grads-local-pos (x nd dx h)
  (let ((dx (- (* dx 2) 1)))
    (/
     (case nd
       (2 (* -1d0 dx))
       (3 (* -1d0 dx)))
     (/ h 1))))

(defun iterate-over-cell-shape-local (mesh cell local-position func)
  "Iterating over a given cell's basis functions"
  (declare (cl-mpm/mesh::mesh mesh))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (nd (cl-mpm/mesh::mesh-nd mesh))
           (pos-vec local-position)
           (cell-pos (cl-mpm/mesh::cell-centroid cell))
           (cell-index (cl-mpm/mesh::cell-index cell)))
      (declare (dynamic-extent cell-index pos-vec))
      (loop for dx from 0 to 1
            do (loop for dy from 0 to 1
                     do (loop for dz from 0 to 1
                              do (let* ((id (mapcar #'+ cell-index (list dx dy dz))))
                                   (declare (dynamic-extent id))
                                   (when (cl-mpm/mesh:in-bounds mesh id)
                                     (let* ((node (cl-mpm/mesh:get-node mesh id))
                                            (dist-x (* 2d0 (/ (- (varef pos-vec 0) (varef cell-pos 0)) h)))
                                            (dist-y (* 2d0 (/ (- (varef pos-vec 1) (varef cell-pos 1)) h)))
                                            (dist-z (* 2d0 (/ (- (varef pos-vec 2) (varef cell-pos 2)) h)))
                                            (weight-x (weights-local-pos dist-x nd dx))
                                            (weight-y (weights-local-pos dist-y nd dy))
                                            (weight-z (if (= nd 3) (weights-local-pos dist-z nd dz) 1d0))
                                            )
                                       (let* ((weight (* weight-x weight-y weight-z))
                                              (grad-x (grads-local-pos dist-x nd dx h))
                                              (grad-y (grads-local-pos dist-y nd dy h))
                                              (grad-z (if (= nd 3) (grads-local-pos dist-z nd dz h) 0d0))
                                              (grads (cl-mpm/shape-function::grads-3d
                                                      (list weight-x weight-y weight-z)
                                                      (list grad-x grad-y grad-z)))
                                              )
                                         (declare
                                          (double-float weight))
                                         (when t
                                           (funcall func node weight grads))))))))))))


(defun iterate-over-midpoints (mesh mp func)
  "Iterates a function (corner) over a single MPs corners"
  (if (= (the fixnum (cl-mpm/mesh:mesh-nd mesh)) 2)
      (iterate-over-midpoints-2d mesh mp func)
      (iterate-over-midpoints-3d mesh mp func)))

(defun iterate-over-cell-linear-3d (mesh cell position func)
  "Iterating over a given cell's basis functions"
  (declare (cl-mpm/mesh::mesh mesh))
  (progn
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (cell-vec (cl-mpm/mesh::cell-centroid cell))
           (pos-vec position)
           (pos (list (tref pos-vec 0 0) (tref pos-vec 1 0) (tref pos-vec 2 0)))
           (pos-index (cl-mpm/mesh:position-to-index-floor mesh cell-vec)))
      (declare (dynamic-extent pos pos-index pos-vec))
      (loop for dx from 0 to 1
            do (loop for dy from 0 to 1
                     do (loop for dz from 0 to 1
                              do (let* ((id (mapcar #'+ pos-index (list dx dy dz))))
                                   (declare (dynamic-extent id))
                                   (when (cl-mpm/mesh:in-bounds mesh id)
                                     (let* ((dist (mapcar #'- pos (cl-mpm/mesh:index-to-position mesh id)))
                                            (node (cl-mpm/mesh:get-node mesh id))
                                            (weights (mapcar (lambda (x) (cl-mpm/shape-function::shape-linear x h)) dist))
                                            (weight (reduce #'* weights))
                                            (lin-grads (mapcar (lambda (d)
                                                                 (cl-mpm/shape-function::shape-linear-dsvp d h))
                                                               dist))
                                            (grads (cl-mpm/shape-function::grads-3d weights lin-grads)))
                                       (declare
                                        (double-float weight)
                                        (dynamic-extent dist weights))
                                       (when t;(< 0d0 weight)
                                         (funcall func node weight grads)))))))))))


