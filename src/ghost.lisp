(defpackage :cl-mpm/ghost
  (:use :cl
   :cl-mpm
   :cl-mpm/particle
   :cl-mpm/mesh
   :cl-mpm/utils
   :cl-mpm/fastmath
   )
  (:import-from
   :magicl tref .+ .-)
  (:export
   ))
(declaim (optimize (debug 3) (safety 3) (speed 0)))
;; (declaim (optimize (debug 0) (safety 0) (speed 3)))
                                        ;    #:make-shape-function
(in-package :cl-mpm/ghost)



(defun find-local-coord (cell pos)
  )
(defun det-element-derivitives (cell)
  )

(defun iterate-over-face-gps (mesh cell-a cell-b func)
  (let ((nd (cl-mpm/mesh:mesh-nd mesh)))
    (ecase nd
        (2
         (iterate-over-face-gps-2d mesh cell-a cell-b func))
        (3
         (error "Not implemented 3d GP face iterateion")))))

(defun iterate-over-face-gps-2d (mesh cell-a cell-b func)
  (declare (function func))
  (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
         (pos-a (cl-mpm/mesh::cell-centroid cell-a))
         (pos-b (cl-mpm/mesh::cell-centroid cell-b))
         (normal (cl-mpm/fastmath:norm (cl-mpm/fastmath:fast-.- pos-a pos-b)))
         (midpoint (cl-mpm/fastmath::fast-scale (cl-mpm/fastmath:fast-.+ pos-a pos-b) 0.5d0))
         (face-basis (cl-mpm/utils:vector-from-list (list (- (tref normal 1 0)) (tref normal 0 0) 0d0))))
    (loop for gp in (list -1d0 1d0)
          do (progn
               (let ((gp-loc (cl-mpm/fastmath:fast-.+
                              midpoint
                              (magicl:scale face-basis (/ (* 0.5d0 h) (sqrt 3)))))
                     (gp-weight (* (/ (expt h 3) 3) (/ h 2d0))))
                 (funcall func gp-loc gp-weight normal))))))


(defun iterate-over-cell-nodes (mesh cell func)
  (declare (function func)
           (cl-mpm/mesh::cell cell)
           (cl-mpm/mesh::mesh mesh))
  (with-accessors ((nodes cl-mpm/mesh::cell-nodes))
      cell
    (loop for node in nodes
          do (funcall func node))))
(defun iterate-over-face-nodes (mesh midpoint normal func)
  (iterate-over-face-nodes-2d
   mesh midpoint normal func))

(defun iterate-over-face-nodes-2d (mesh midpoint normal func)
  (let ((root (cl-mpm/mesh::position-to-index mesh midpoint #'floor))
        (face-basis (list (if (< 1d-15 (abs (tref normal 0 0))) 0 1)
                          (if (< 1d-15 (abs (tref normal 1 0))) 0 1)
                          0))
        (normal-basis (list (if (< 1d-15 (abs (tref normal 0 0))) 1 0)
                            (if (< 1d-15 (abs (tref normal 1 0))) 1 0)
                            0
                            ))
        )
    (loop for dface in '(0 1)
          do
             (loop for dnorm in '(-1 0 1)
                   do
                      (let ((index (mapcar #'+
                                           root
                                           (mapcar (lambda (x) (* x dface)) face-basis)
                                           (mapcar (lambda (x) (* x dnorm)) normal-basis))))
                        (funcall func (cl-mpm/mesh:get-node mesh index)))))))

(defun iterate-over-gp-nodes (mesh gp-point normal func)
  (let ((h (cl-mpm/mesh:mesh-resolution mesh)))
    (iterate-over-face-nodes
     mesh
     gp-point
     normal
     (lambda (node)
       (let* ((node-pos (cl-mpm/mesh::node-position node))
              (dist (cl-mpm/fastmath::fast-.- gp-point node-pos))
              (weights (list
                        (cl-mpm/shape-function::shape-linear (tref dist 0 0) h)
                        (cl-mpm/shape-function::shape-linear (tref dist 1 0) h)))
              (weight (reduce #'* weights))
              (grads (list
                      (* (the double-float (cl-mpm/shape-function::shape-linear-dsvp (tref dist 0 0) h))
                         (the double-float (nth 1 weights)))
                      (* (the double-float (cl-mpm/shape-function::shape-linear-dsvp (tref dist 1 0) h))
                         (the double-float (nth 0 weights)))
                      0d0
                      )))
         (funcall func node weight grads))))))

(defun apply-ghost-cells (mesh cell-a cell-b)
  (iterate-over-face-gps
   mesh
   cell-a
   cell-b
   (lambda (gp-loc gp-weight normal)
     (iterate-over-gp-nodes
      mesh
      gp-loc
      normal
      (lambda (node weight grads)
        ;;Add force to node
        (with-accessors ((ghost cl-mpm/mesh::node-ghost))
            node
          (let ((dsvp (magicl:transpose! (cl-mpm/shape-function::assemble-dsvp-3d grads))))
            (iterate-over-gp-nodes
             mesh
             gp-loc
             normal
             (lambda (node-b weight-b grads-b)
               (with-accessors ((acc cl-mpm/mesh::node-acceleration))
                   node-b
                 (let* ((grad-adjust (list (* (nth 0 grads-b) (tref normal 0 0))
                                          (* (nth 1 grads-b) (tref normal 1 0))
                                          0d0))
                        (dsvp-b (cl-mpm/shape-function::assemble-dsvp-3d grad-adjust)))
                   (cl-mpm/fastmath:fast-.+
                    ghost
                    (magicl:scale!
                     (magicl:@
                      dsvp
                      dsvp-b
                      acc)
                     weight)
                    ghost))))))))))))

(defun apply-ghost (sim)
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (with-accessors ((cells cl-mpm/mesh::mesh-cells))
        mesh
      (lparallel:pdotimes (i (array-total-size cells))
        (let ((cell (row-major-aref cells i)))
          (let* ((index (cl-mpm/mesh::cell-index cell))
                 (cell (cl-mpm/mesh::get-cell mesh index)))
            (loop for direction from 0 to 2
                  do
                     (let ((ind-dir (list 0 0 0)))
                       (setf (nth direction ind-dir) 1)
                       (let ((index-b (mapcar #'+ index ind-dir)))
                         (when (cl-mpm/mesh::in-bounds-cell mesh index-b)
                           (apply-ghost-cells mesh cell
                                              (cl-mpm/mesh::get-cell mesh index-b))))))))))))
