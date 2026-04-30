(defpackage :cl-mpm/mpi
  (:use
   :cl
   :cl-mpm
   :cl-mpm/particle
   :cl-mpm/mesh
   :cl-mpm/utils
   :cl-mpm/fastmaths
   )
  (:import-from
   :magicl tref .+ .-
   )
  (:import-from
   :trivial-with-current-source-form with-current-source-form
   )
  (:export
   :mpm-sim-mpi
   :mpm-sim-mpi-nodes
   :mpm-sim-mpi-nodes-damage
   #:domain-decompose
   #:mpi-average
   #:mpi-sum
   ))

(declaim (optimize (debug 0) (safety 0) (speed 3)))
(in-package :cl-mpm/mpi)

(defclass mpm-sim-mpi (cl-mpm::mpm-sim-usf)
  ((neighbour-node-list)
   (neighbour-ranks
    :initform '())
   (domain-index
    :accessor mpm-sim-mpi-domain-index
    :initform '(0 0 0))
   (domain-bounds
    :accessor mpm-sim-mpi-domain-bounds
    :initform nil)
   (halo-nd
    :accessor mpm-sim-mpi-nd
    :initform 3)
   (halo-depth
    :accessor mpm-sim-mpi-halo-depth
    :initform 3d0
    )
   (min-size
    :accessor mpm-sim-mpi-min-size
    :initform 2d0
    )
   (domain-count
    :accessor mpm-sim-mpi-domain-count
    :initform '(1 1 1)
    :initarg :domain-count
    )
   (work-time
    :accessor mpm-sim-mpi-work-time
    :initform 0d0)
   (load-metric
    :accessor mpm-sim-mpi-load-metric
    :initform 0d0)
   (halo-damage-size
    :accessor mpm-sim-mpi-halo-damage-size
    :initform 1d0)
   (damage-mps-cache
    :accessor mpm-sim-mpi-damage-mps-cache
    :initform (make-array 0 :element-type t :adjustable t :fill-pointer 0)))
  (:documentation "Damage sim with only stress update on mpi"))

(defclass mpm-sim-mpi-nodes-damage (mpm-sim-mpi-nodes cl-mpm/damage::mpm-sim-damage)
  ())
(defclass mpm-sim-usl-mpi-nodes-damage (mpm-sim-mpi-nodes-damage)
  ())

(defclass mpm-sim-mpi-nodes (mpm-sim-mpi)
  ((halo-node-list
    :accessor mpm-sim-mpi-halo-node-list
    :initform (loop for i from 0 to 2
                    collect (loop for i from 0 to 1 collect (make-array 0 :element-type t))))))

(defmacro rank-0-time (rank &rest body)
  `(if (= ,rank 0)
       (time
        (progn
          ,@body))
       (progn
         ,@body)))

(defun mpi-rank-to-index (sim rank)
  (let ((size (mpm-sim-mpi-domain-count sim)))
    (let ((index (list
                  (truncate (floor rank (* (nth 2 size) (nth 1 size))))
                  (truncate (floor (mod rank (* (nth 2 size) (nth 1 size))) (nth 2 size)))
                  (truncate (mod rank (nth 2 size)))
                  )))
      (if (mpi-index-in-bounds sim index)
          index
          (list -1 -1 -1))
      )))

(defun mpi-index-to-rank (sim index)
  (if (mpi-index-in-bounds sim index)
      (let ((size (mpm-sim-mpi-domain-count sim)))
        (destructuring-bind (s-1 s-2 s-3) size
          (destructuring-bind (i-1 i-2 i-3) index
            (declare (fixnum s-1 s-2 s-3)
                     (fixnum i-1 i-2 i-3))
            (+ i-3
               (* s-3
                  (+ i-2
                     (* s-2 i-1)))))
          ;; (+ (nth 2 index)
          ;;    (* (nth 2 size)
          ;;       (+ (nth 1 index)
          ;;          (* (nth 1 size) (nth 0 index)))))
          ))
      -1))
(defun mpi-index-in-bounds (sim index)
  (let ((size (mpm-sim-mpi-domain-count sim)))
    (destructuring-bind (s-1 s-2 s-3) size
      (destructuring-bind (i-1 i-2 i-3) index
        (declare (fixnum s-1 s-2 s-3)
                 (fixnum i-1 i-2 i-3))
        (and
         (and (>= i-1 0) (< i-1 s-1))
         (and (>= i-2 0) (< i-2 s-2))
         (and (>= i-3 0) (< i-3 s-3)))))))


(defun mpi-average (value mp-count)
  "Average a single pre-reduced value which represents the reduction of N local samples"
  (let ((sum (mpi-sum value))
        (total-mp-count (mpi-sum mp-count)))
    ;;Don't care about zero sums
    (setf sum (/ sum (max 1d0 total-mp-count)))
    sum))

(defun mpi-sum (value)
  "Sum a scalar over all mpi nodes"
  (cl-mpi:mpi-waitall)
  (let ((output 0d0))
    ;; (static-vectors:with-static-vector (source 1 :element-type 'double-float
    ;;                                              :initial-element (coerce value 'double-float))
    ;;   (static-vectors:with-static-vector (dest 1 :element-type 'double-float :initial-element 0d0)
    ;;     (cl-mpi:mpi-allreduce source dest cl-mpi:+mpi-sum+)
    ;;     (setf output (aref dest 0))))
    (cffi:with-foreign-objects ((in :double)
                                (out :double))
      (setf (cffi:mem-ref in :double) (coerce value 'double-float))
      (mpi::%mpi-allreduce in out 1 mpi:+mpi-double+ mpi:+mpi-sum+ mpi:*standard-communicator*)
      (setf output (cffi:mem-ref out :double)))
    output))

(defun mpi-max (value)
  (cl-mpi:mpi-waitall)
  ;; (static-vectors:with-static-vector (source 1 :element-type 'double-float :initial-element (coerce value 'double-float))
  ;;   (static-vectors:with-static-vector (dest 1 :element-type 'double-float :initial-element 0d0)
  ;;     (cl-mpi:mpi-allreduce source dest cl-mpi:+mpi-max+)
  ;;     (aref dest 0)))

  (cffi:with-foreign-objects ((in :double)
                              (out :double))
    (setf (cffi:mem-ref in :double) (coerce value 'double-float))
    (mpi::%mpi-allreduce in out 1 mpi:+mpi-double+ mpi:+mpi-max+ mpi:*standard-communicator*)
    (cffi:mem-ref out :double))
  )

(defun mpi-min (value)
  (cl-mpi:mpi-waitall)
  (static-vectors:with-static-vector (source 1 :element-type 'double-float :initial-element (coerce value 'double-float))
    (static-vectors:with-static-vector (dest 1 :element-type 'double-float :initial-element 0d0)
      (cl-mpi:mpi-allreduce source dest cl-mpi:+mpi-min+)
      (aref dest 0))))


(defun mpi-vector-integer-sum (values)
  "Sum a scalar over all mpi nodes"
  (cl-mpi:mpi-waitall)
  (let ((size (length values)))
    (static-vectors:with-static-vector (source size :element-type '(signed-byte 32)
                                                    :initial-contents values)
      (static-vectors:with-static-vector (dest size :element-type '(signed-byte 32)
                                                    :initial-element 0)
        (cl-mpi:mpi-allreduce source dest cl-mpi:+mpi-sum+)
        (cl-mpi:mpi-waitall)
        (aops:copy-into values dest)))
    values))

(defun mpi-vector-integer-max (values)
  "Sum a scalar over all mpi nodes"
  (cl-mpi:mpi-waitall)
  (let ((size (length values)))
    (static-vectors:with-static-vector (source size :element-type '(signed-byte 32)
                                                    :initial-contents values)
      (static-vectors:with-static-vector (dest size :element-type '(signed-byte 32)
                                                    :initial-element 0)
        (cl-mpi:mpi-allreduce source dest cl-mpi:+mpi-max+)
        (cl-mpi:mpi-waitall)
        (aops:copy-into values dest)))
    values))

(defun mpi-vector-sum (values)
  "Sum a scalar over all mpi nodes"
  (cl-mpi:mpi-waitall)
  (let ((size (length values)))
    (static-vectors:with-static-vector (source size :element-type 'double-float
                                                    :initial-contents values)
      (static-vectors:with-static-vector (dest size :element-type 'double-float :initial-element 0d0)
        (cl-mpi:mpi-allreduce source dest cl-mpi:+mpi-sum+)
        (cl-mpi:mpi-waitall)
        (aops:copy-into values dest)))
    values))


(defun mpi-vector-max (values)
  "Sum a scalar over all mpi nodes"
  (cl-mpi:mpi-waitall)
  (let ((size (length values)))
    (static-vectors:with-static-vector (source size :element-type 'double-float
                                                    :initial-contents values)
      (static-vectors:with-static-vector (dest size :element-type 'double-float :initial-element 0d0)
        (cl-mpi:mpi-allreduce source dest cl-mpi:+mpi-max+)
        (cl-mpi:mpi-waitall)
        (aops:copy-into values dest)))
    values))
