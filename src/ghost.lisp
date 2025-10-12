(defpackage :cl-mpm/ghost
  (:use :cl
   :cl-mpm
   :cl-mpm/particle
   :cl-mpm/mesh
   :cl-mpm/utils
   :cl-mpm/fastmaths
   )
  (:import-from
   :magicl tref .+ .-)
  (:export
   #:apply-ghost))
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

(declaim (notinline get-face-displacement))
(defun get-face-displacement (mesh normal point)
  (let ((not-normal (cl-mpm/fastmaths:fast-.-
                     (cl-mpm/utils:vector-from-list (list 1d0 1d0 1d0))
                     normal))
        (dF (cl-mpm/utils:voigt-zeros))
        )
    (cl-mpm::iterate-over-neighbours-point-linear
     mesh
     point
     (lambda (mesh node weight grads)
       (cl-mpm/fastmaths::fast-.+
        (magicl:@
         (cl-mpm/shape-function::assemble-dsvp-3d grads)
         (cl-mpm/fastmaths:fast-.*
          not-normal
          (cl-mpm/mesh::node-displacment node)))
        dF
        dF)))
    (* (cl-mpm/mesh::mesh-resolution mesh) (+ 1d0 (cl-mpm/fastmaths::voight-det dF)))))

(defun get-face-df (mesh normal point)
  (let ((not-normal (cl-mpm/fastmaths:fast-.-
                     (cl-mpm/utils:vector-from-list (list 1d0 1d0 1d0))
                     normal))
        (dF (cl-mpm/utils:matrix-eye 1d0)))
    (iterate-over-gp-nodes
     mesh
     point
     normal
     (lambda (node weight grads dface dnorm)
       (cl-mpm/shape-function::@-combi-assemble-dstretch-3d grads (cl-mpm/mesh::node-displacment node) dF)))
    dF))

(defun iterate-over-face-gps-2d (mesh cell-a cell-b func)
  (declare (function func))
  (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
         (pos-a (cl-mpm/mesh::cell-centroid cell-a))
         (pos-b (cl-mpm/mesh::cell-centroid cell-b))
         (normal (cl-mpm/fastmaths:norm (cl-mpm/fastmaths:fast-.- pos-a pos-b)))
         (pos-trial-a (cl-mpm/mesh::cell-trial-centroid cell-a))
         (pos-trial-b (cl-mpm/mesh::cell-trial-centroid cell-b))
         (normal-trial (cl-mpm/fastmaths:norm (cl-mpm/fastmaths:fast-.- pos-trial-a pos-trial-b)))
         (midpoint (cl-mpm/fastmaths::fast-scale! (cl-mpm/fastmaths:fast-.+ pos-a pos-b) 0.5d0))
         (face-basis (cl-mpm/utils:vector-from-list
                      (list
                       (- 1d0 (abs (varef normal 0)))
                       (- 1d0 (abs (varef normal 1)))
                       0d0)))
         (root (cl-mpm/mesh::position-to-index mesh midpoint #'floor))
         (end (mapcar #'+ root
                      (list
                       (round (varef face-basis 0))
                       (round (varef face-basis 1))
                       0)))
         (length (sqrt
                  (max 0d0
                       (cl-mpm/fastmaths::diff-norm
                        (cl-mpm/mesh::node-position (cl-mpm/mesh:get-node mesh root))
                        (cl-mpm/mesh::node-position (cl-mpm/mesh:get-node
                                                     mesh
                                                     end)))))))
    ;; (format t "~A ~A ~A~%" normal root end)
    (loop for gp in (list -1d0 1d0)
          do (progn
               (let ((gp-loc (cl-mpm/fastmaths:fast-.+
                              midpoint
                              (magicl:scale face-basis (/ (* gp 0.5d0 h) (sqrt 3)))))
                     (gp-weight 1d0))
                 (funcall func gp-loc gp-weight normal normal-trial length))))))


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
(defun iterate-over-face-nodes-positive (mesh midpoint normal func)
  (iterate-over-face-nodes-2d
   mesh midpoint normal func
   '(-1 0)))

(defun iterate-over-face-nodes-2d (mesh midpoint normal func &optional (dnormal '(-1 0 1)))
  (let* ((root (cl-mpm/mesh::position-to-index mesh midpoint #'floor))
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
             (loop for dnorm in dnormal
                   do
                      (let* ((index (mapcar #'+
                                           root
                                           (mapcar (lambda (x) (* x dface)) face-basis)
                                           (mapcar (lambda (x) (* x dnorm)) normal-basis))))
                        (when (cl-mpm/mesh:in-bounds mesh index)
                          (let ((node (cl-mpm/mesh:get-node mesh index)))
                            (when (cl-mpm/mesh::node-active node)
                              (funcall func node dface dnorm)))))))))

(defun iterate-over-gp-nodes (mesh gp-point normal func)
  (let ((h (cl-mpm/mesh:mesh-resolution mesh)))
    (iterate-over-face-nodes
     mesh
     gp-point
     normal
     (lambda (node dface dnorm)
       (let* ((node-pos (cl-mpm/mesh::node-position node))
              (dist (cl-mpm/fastmaths::fast-.- gp-point node-pos))
              (weight-x (cl-mpm/shape-function::shape-linear (varef dist 0) h))
              (weight-y (cl-mpm/shape-function::shape-linear (varef dist 1) h))
              (weight (* weight-x weight-y))
              (grads (list
                      (*
                       (signum (varef dist 0))
                       (the double-float (cl-mpm/shape-function::shape-linear-dsvp (varef dist 0) h))
                       (the double-float weight-y))
                      (*
                       (signum (varef dist 1))
                       (the double-float (cl-mpm/shape-function::shape-linear-dsvp (varef dist 1) h))
                       (the double-float weight-x))
                      0d0)))
         (funcall func node weight grads dface dnorm))))))
(defun iterate-over-gp-nodes-positive (mesh gp-point normal func)
  (let ((h (cl-mpm/mesh:mesh-resolution mesh)))
    (iterate-over-face-nodes-positive
     mesh
     gp-point
     normal
     (lambda (node dface dnorm)
       (let* ((node-pos (cl-mpm/mesh::node-position node))
              (dist (cl-mpm/fastmaths::fast-.- gp-point node-pos))
              (weight-x (cl-mpm/shape-function::shape-linear (varef dist 0) h))
              (weight-y (cl-mpm/shape-function::shape-linear (varef dist 1) h))
              (weight (* weight-x weight-y))
              (grads (list
                      (*
                       (signum (varef dist 0))
                       (the double-float (cl-mpm/shape-function::shape-linear-dsvp (varef dist 0) h))
                       (the double-float weight-y))
                      (*
                       (signum (varef dist 1))
                       (the double-float (cl-mpm/shape-function::shape-linear-dsvp (varef dist 1) h))
                       (the double-float weight-x))
                      0d0)))
         (funcall func node weight grads dface dnorm))))))



(defun populate-cell-mp-count-volume (mesh mps)
  (let ((cells (cl-mpm/mesh::mesh-cells mesh)))
    (cl-mpm::iterate-over-cells
     mesh
     (lambda (cell)
       (setf (cl-mpm/mesh::cell-mp-count cell) 0
             (cl-mpm/mesh::cell-ghost-element cell) nil)))
    (cl-mpm::iterate-over-cells
     mesh
     (lambda (cell)
       (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                        (neighbours cl-mpm/mesh::cell-neighbours)
                        (index cl-mpm/mesh::cell-index)
                        (nodes cl-mpm/mesh::cell-nodes)
                        (centroid cl-mpm/mesh::cell-centroid))
           cell
         (when (every
                (lambda (n)
                  (and (cl-mpm/mesh:node-active n)
                       (> (cl-mpm/mesh:node-mass n) 0d0)))
                nodes)
           (setf mp-count 1)))))
    ))

(defun locate-ghost-elements (sim)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (cl-mpm::iterate-over-cells
     mesh
     (lambda (cell)
       (with-accessors ((active cl-mpm/mesh::cell-active)
                        (partial cl-mpm/mesh::cell-partial)
                        (ghost cl-mpm/mesh::cell-ghost-element))
           cell
         (when active
           (setf ghost nil)))))
    (cl-mpm::iterate-over-cells
     mesh
     (lambda (cell)
       (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                        (neighbours cl-mpm/mesh::cell-neighbours)
                        (index cl-mpm/mesh::cell-index)
                        (nodes cl-mpm/mesh::cell-nodes)
                        (pruned cl-mpm/mesh::cell-pruned)
                        (active cl-mpm/mesh::cell-active)
                        (partial cl-mpm/mesh::cell-partial)
                        (ghost cl-mpm/mesh::cell-ghost-element)
                        (pos cl-mpm/mesh::cell-centroid))
           cell
         (when (and active (not partial))
           (loop for neighbour in neighbours
                 do
                    (when (and
                           (cl-mpm/mesh::cell-active neighbour)
                           (cl-mpm/mesh::cell-partial neighbour))
                      (setf ghost t
                            (cl-mpm/mesh::cell-ghost-element neighbour) t)))))))))

(defun iterate-over-bicell-nodes (mesh cell-a cell-b gp-loc func)
  (declare (function func))
  (cl-mpm::iterate-over-cell-shape-local
   mesh
   cell-a
   gp-loc
   (lambda (node weight grads) (funcall func cell-a node weight grads 1d0)))
  (cl-mpm::iterate-over-cell-shape-local
   mesh
   cell-b
   gp-loc
   (lambda (node weight grads) (funcall func cell-b node weight grads -1d0))))

(defparameter *write-lock* (sb-thread:make-mutex))
(defun apply-ghost-cells-new (mesh cell-a cell-b ghost-factor)
  (when
      (and
       (and
        (= (cl-mpm/mesh::cell-mp-count cell-a) 1)
        (= (cl-mpm/mesh::cell-mp-count cell-b) 1))
       (or
        (cl-mpm/mesh::cell-ghost-element cell-a)
        (cl-mpm/mesh::cell-ghost-element cell-b)))
    (with-accessors ((h cl-mpm/mesh:mesh-resolution))
        mesh
      (let ((h (cl-mpm/mesh:mesh-resolution mesh)))
        (iterate-over-face-gps
         mesh
         cell-a
         cell-b
         (lambda (gp-loc gp-weight normal normal-trial face-length)
           (let* ((nx (varef normal-trial 0))
                  (ny (varef normal-trial 1))
                  (nz (varef normal-trial 2))
                  (dsvp-adjuster (magicl:transpose! (cl-mpm/utils::arb-matrix-from-list
                                                     ;;3d case
                                                     ;; (list
                                                     ;;  nx  0d0 0d0
                                                     ;;  0d0 ny  0d0
                                                     ;;  0d0 0d0 nz
                                                     ;;  ny  nx  0d0
                                                     ;;  0d0 nz  ny
                                                     ;;  nz  0d0 nx)
                                                     ;;2d case
                                                     (list
                                                      nx  0d0 0d0
                                                      0d0 ny  0d0
                                                      0d0 0d0 0d0
                                                      0d0 0d0 0d0
                                                      0d0 0d0 0d0
                                                      0d0 nx  ny
                                                      )
                                                     3
                                                     6))))
             (cl-mpm/ghost::iterate-over-bicell-nodes
              mesh
              cell-a
              cell-b
              gp-loc
              (lambda (cell-p node weight grads fact-p)
                ;;Iterate over every node in A and B
                (with-accessors ((ghost cl-mpm/mesh::node-ghost-force)
                                 (f-int cl-mpm/mesh::node-internal-force)
                                 (disp cl-mpm/mesh::node-displacment))
                    node
                  (let ((dsvp-p
                          (fast-scale!
                           (magicl:transpose! (cl-mpm/shape-function::assemble-dsvp-3d
                                               grads
                                               ;; (cl-mpm::gradient-push-forwards
                                               ;;  grads
                                               ;;  (cl-mpm/mesh::cell-deformation-gradient cell-p))
                                               ))
                           (* fact-p)
                           ;; -1d0
                           )))
                    (iterate-over-bicell-nodes
                     mesh
                     cell-a
                     cell-b
                     gp-loc
                     (lambda (cell-n node-b weight-b grads-b fact-n)
                       (with-accessors ((disp-b cl-mpm/mesh::node-displacment)
                                        (ghost-b cl-mpm/mesh::node-ghost-force))
                           node-b
                         (let* ((dsvp-n
                                  (fast-scale!
                                   (magicl:transpose!
                                    (cl-mpm/shape-function::assemble-dsvp-3d
                                     grads-b
                                     ;; (cl-mpm::gradient-push-forwards
                                     ;;  grads-b
                                     ;;  (cl-mpm/mesh::cell-deformation-gradient cell-n))
                                     ))
                                   (* fact-n)
                                   ;; 1d0
                                   ))
                                (dsvp (cl-mpm/fastmaths:fast-.- dsvp-p dsvp-n)))
                           (let ((ghost-mat
                                   (cl-mpm/fastmaths:fast-scale!
                                    (magicl:@
                                     dsvp-p
                                     dsvp-adjuster
                                     (magicl:transpose dsvp-adjuster)
                                     (magicl:transpose dsvp-n))
                                    (*
                                     -1d0
                                     (/ (* ghost-factor (expt h 3)) 6)
                                     face-length
                                     gp-weight))
                                   ))
                             ;; (format t "GS~%")
                             ;; (format t "~A~%" dsvp-adjuster)
                             ;; (format t "~A~%" ghost-mat)
                             (cl-mpm/fastmaths:fast-.+
                              (magicl:@
                               ghost-mat
                               disp-b)
                              ghost
                              ghost)
                             (cl-mpm/fastmaths:fast-.+
                              (magicl:@
                               ghost-mat
                               disp-b)
                             f-int
                             f-int)
                             ;; (cl-mpm/fastmaths:fast-.-
                             ;;  ghost-b
                             ;;  (magicl:@
                             ;;   ghost-mat
                             ;;   disp)
                             ;;  ghost-b)
                             )

                           ))))
                    ))))
             )))))))



;; (defun apply-ghost-cells (mesh cell-a cell-b ghost-factor)
;;   (when
;;       ;;XOR active element
;;       (and
;;        (and
;;         (= (cl-mpm/mesh::cell-mp-count cell-a) 1)
;;         (= (cl-mpm/mesh::cell-mp-count cell-b) 1))
;;        (or
;;         (cl-mpm/mesh::cell-ghost-element cell-a)
;;         (cl-mpm/mesh::cell-ghost-element cell-b)))
;;     ;; (setf (cl-mpm/mesh::cell-ghost-element cell-a) t
;;     ;;       (cl-mpm/mesh::cell-ghost-element cell-b) t
;;     ;;       )
;;     (with-accessors ((h cl-mpm/mesh:mesh-resolution))
;;         mesh
;;       (iterate-over-face-gps
;;        mesh
;;        cell-a
;;        cell-b
;;        (lambda (gp-loc gp-weight normal normal-trial)
;;          (let (;; (face-length (get-face-displacement mesh normal gp-loc))
;;                (df (get-face-df mesh normal gp-loc))
;;                ;; (face-length h)
;;                (face-length (get-face-displacement mesh normal gp-loc))
;;                )
;;            ;; (pprint face-length)
;;            (iterate-over-gp-nodes;-positive
;;             mesh
;;             gp-loc
;;             normal
;;             (lambda (node weight grads dface dnorm)
;;               ;;Add force to node
;;               (with-accessors ((ghost cl-mpm/mesh::node-ghost-force)
;;                                (disp cl-mpm/mesh::node-displacment))
;;                   node
;;                 (let* ((nx (varef normal-trial 0))
;;                        (ny (varef normal-trial 1))
;;                        (nz (varef normal-trial 2))
;;                        (dsvp-adjuster (magicl:transpose! (cl-mpm/utils::arb-matrix-from-list
;;                                                           (list
;;                                                            nx  0d0 0d0
;;                                                            0d0 ny  0d0
;;                                                            0d0 0d0 nz
;;                                                            ny  nx  0d0
;;                                                            0d0 nz  ny
;;                                                            nz  0d0 nx)
;;                                                           3
;;                                                           6)))
;;                        (dsvp-p (magicl:transpose! (cl-mpm/shape-function::assemble-dsvp-3d grads))))
;;                   (iterate-over-gp-nodes
;;                    mesh
;;                    gp-loc
;;                    normal
;;                    (lambda (node-b weight-b grads-b dface-b dnorm-b)
;;                      (when (not (eq node node-b))
;;                        (with-accessors ((disp-b cl-mpm/mesh::node-displacment)
;;                                         (ghost-b cl-mpm/mesh::node-ghost-force))
;;                            node-b
;;                          (let* ((dsvp-n (magicl:transpose! (cl-mpm/shape-function::assemble-dsvp-3d grads-b)))
;;                                 (dsvp (cl-mpm/fastmaths:fast-.- dsvp-p dsvp-n))
;;                                 )
;;                            (let ((ghost-mat
;;                                    (cl-mpm/fastmaths:fast-scale!
;;                                     (magicl:@
;;                                      dsvp
;;                                      dsvp-adjuster
;;                                      (magicl:transpose dsvp-adjuster)
;;                                      (magicl:transpose dsvp)
;;                                      )
;;                                     (*
;;                                      -1d0
;;                                      (/ (* ghost-factor (expt h 3)) 6)
;;                                      ;; (* dface dface-b)
;;                                      face-length
;;                                      gp-weight))
;;                                    ))
;;                              (cl-mpm/fastmaths:fast-.+
;;                               ghost
;;                               (magicl:@
;;                                ghost-mat
;;                                disp-b)
;;                               ghost)
;;                              ;; (cl-mpm/fastmaths:fast-.+
;;                              ;;  ghost-b
;;                              ;;  (magicl:@
;;                              ;;   ghost-mat
;;                              ;;   disp)
;;                              ;;  ghost-b)
;;                              ))))))))))))))
;;     ))

;; (let* ((nx 1d0)
;;        (ny 0d0)
;;        (nz 0d0)
;;        (dsvp (cl-mpm/shape-function::assemble-dsvp-3d (list 1d0 2d0 3d0)))
;;        (n-adj
;;          (cl-mpm/utils::arb-matrix-from-list
;;           (list
;;            nx  0d0 0d0
;;            0d0 ny  0d0
;;            0d0 0d0 nz
;;            ny  nx  0d0
;;            0d0 nz  ny
;;            nz  0d0 nx) 3 6))
;;        )
;;   (pprint n-adj)
;;   (pprint dsvp)
;;   (pprint (magicl:@ (magicl:transpose dsvp) (magicl:transpose n-adj))))

(declaim (notinline apply-ghost))
(defun apply-ghost (sim ghost-factor)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps))
      sim
    (with-accessors ((cells cl-mpm/mesh::mesh-cells)
                     (h cl-mpm/mesh::mesh-resolution))
        mesh
      (declare (double-float h ghost-factor))
      (locate-ghost-elements sim)
      (cl-mpm::iterate-over-cells
       mesh
        (lambda (cell)
          (let* ((index (cl-mpm/mesh::cell-index cell))
                 (cell (cl-mpm/mesh::get-cell mesh index)))
            (when (cl-mpm/mesh::cell-ghost-element cell)
                (loop for direction from 0 to 2
                      do
                         (let ((ind-dir (list 0 0 0)))
                           (setf (nth direction ind-dir) 1)
                           (let ((index-b (mapcar #'+ index ind-dir)))
                             (when (cl-mpm/mesh::in-bounds-cell mesh index-b)
                               (apply-ghost-cells-new mesh cell
                                                      (cl-mpm/mesh::get-cell mesh index-b)
                                                      ghost-factor
                                                      )))))))))
      (let ((gf (* 4d0 ghost-factor (expt h 2))))
        (cl-mpm::iterate-over-cells
         mesh
         (lambda (cell)
           (when (and (cl-mpm/mesh::cell-active cell)
                      (cl-mpm/mesh::cell-ghost-element cell))
             (loop for n in (cl-mpm/mesh::cell-nodes cell)
                   do (when (cl-mpm/mesh::node-active n)
                        (sb-thread:with-mutex ((cl-mpm/mesh:node-lock n))
                          ;; (setf (cl-mpm/mesh::node-mass n)
                          ;;       (max (cl-mpm/mesh::node-mass n)
                          ;;            ghost-factor))
                          (incf (cl-mpm/mesh::node-mass n)
                                gf)
                          ))))))))))

(defun test-markup (sim)
  (let* ((index (list 10 10 0))
         (mesh (cl-mpm:sim-mesh sim))
         (cell (cl-mpm/mesh::get-cell mesh index)))
    (loop for direction from 0 to 0
          do
             (let ((ind-dir (list 0 0 0)))
               (setf (nth direction ind-dir) 1)
               (let ((index-b (mapcar #'+ index ind-dir)))
                 (when (cl-mpm/mesh::in-bounds-cell mesh index-b)
                   (iterate-over-face-gps
                    mesh
                    cell
                    (cl-mpm/mesh::get-cell mesh index-b)
                    (lambda (gp-loc gp-weight normal normal-trial)
                      (pprint gp-loc)
                      (iterate-over-gp-nodes
                       mesh
                       gp-loc
                       normal
                       (lambda (node weight grad dface dnorm)
                         (pprint node)
                         (setf (cl-mpm/mesh::node-active node) t))
                       )
                      )
                    )
                   ))))))

(defun calculate-forces-ghost (node damping dt mass-scale)
  "Update forces and nodal velocities with viscous damping"
  (when (cl-mpm/mesh:node-active node)
    (with-accessors ((mass node-mass)
                     (vel node-velocity)
                     (force node-force)
                     (force-ext cl-mpm/mesh::node-external-force)
                     (force-int cl-mpm/mesh::node-internal-force)
                     (force-damp cl-mpm/mesh::node-damping-force)
                     (force-ghost cl-mpm/mesh::node-ghost-force)
                     (residual cl-mpm/mesh::node-residual)
                     (residual-prev cl-mpm/mesh::node-residual-prev)
                     (acc node-acceleration))
        node
      (declare (double-float mass dt damping mass-scale))
      (progn
        ;; (cl-mpm/fastmaths:fast-zero acc)
        (cl-mpm/fastmaths:fast-fmacc acc force-ghost (/ 1d0 (* mass mass-scale)))
        (cl-mpm/fastmaths:fast-fmacc vel force-ghost (/ dt (* mass mass-scale)))
        ;; (cl-mpm/fastmaths:fast-fmacc vel acc dt)
        ;; (cl-mpm/utils::vector-copy-into residual residual-prev)
        ;; (cl-mpm/fastmaths::fast-.+-vector force-int force-ext residual)
        ;; (cl-mpm/fastmaths::fast-.+-vector force-ghost residual residual)
        )))
  (values))
(defun update-node-forces-ghost (sim)
  (with-accessors ((damping sim-damping-factor)
                   (mass-scale sim-mass-scale)
                   (mesh sim-mesh)
                   (damping-algo sim-damping-algorithm)
                   (dt sim-dt))
      sim
    (iterate-over-nodes
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (calculate-forces-ghost node damping dt mass-scale))))))

(defun apply-half-step-ghost (sim)
  (with-accessors ((ghost-factor cl-mpm::sim-ghost-factor)
                   (mesh cl-mpm::sim-mesh)
                   (bcs cl-mpm::sim-bcs)
                   (dt cl-mpm::sim-dt))
      sim
    (when ghost-factor

      (cl-mpm::reset-node-displacement sim)
      (cl-mpm::update-nodes sim)
      (cl-mpm::update-cells sim)
      (cl-mpm::apply-bcs mesh bcs dt)
      (cl-mpm/ghost::apply-ghost sim ghost-factor)
      (update-node-forces-ghost sim)
      (cl-mpm::apply-bcs mesh bcs dt)
      (cl-mpm::update-nodes sim)
      (cl-mpm::update-cells sim)
      )))
