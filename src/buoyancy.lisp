(defpackage :cl-mpm/buoyancy
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
(declaim (optimize (debug 0) (safety 0) (speed 3)))
;    #:make-shape-function
(in-package :cl-mpm/buoyancy)

;(defgeneric virtual-stress ())
(let ((rho-true (/ 1000.0d0 1d0)))
  (defun pressure-at-depth (z datum-true)
    (let* ((rho rho-true)
           (g -9.8)
           (datum datum-true)
           (h (- datum z))
           (f (* rho g h))
           )
      (if (> h 0d0)
          f
          0d0)
      ))
  (declaim (notinline buoyancy-virtual-stress))
  (defun buoyancy-virtual-stress (z datum-true rho)
    (let* (;(rho rho-true)
           (g -9.8d0)
           (datum datum-true)
           (h (- datum z))
           (f (* rho g h))
           (p -1d5)
           )
      (if (> h 0d0)
          (magicl:from-list (list f f 0d0) '(3 1) :type 'double-float)
          ;; (magicl:zeros '(3 1))
          (magicl:zeros '(3 1))
          )
      ;; (magicl:from-list (list 0d0 p 0d0)
      ;;                   '(3 1)
      ;;                   :type 'double-float)
      ))
  (defun buoyancy-virtual-div (z datum-true rho)
    (let* (;(rho rho-true)
           (g -9.8d0)
           (datum datum-true)
           (h (- datum z))
           (f (* -1d0 rho g))
           )
      (if (> h 0d0)
          (magicl:from-list (list 0 f) '(2 1) :type 'double-float)
          (magicl:zeros '(2 1))
          )
      )))
(defun pressure-virtual-stress (pressure-x pressure-y)
  (magicl:from-list (list pressure-x pressure-y 0d0)
                    '(3 1)
                    :type 'double-float)
  )
(defun pressure-virtual-div ()
  (magicl:zeros '(2 1)))

;; (defun pressure-virtual-stress ()
;;   (let* ((p -1d3)
;;          )
;;     (magicl:from-list (list p p 0d0)
;;                       '(3 1)
;;                       :type 'double-float)
;;     ))
;; (defun pressure-virtual-div ()
;;   (magicl:zeros '(2 1)))

(declaim (ftype (function (cl-mpm/particle:particle function) (values)) calculate-val-mp))
(defun calculate-val-mp (mp func)
  (with-accessors ((pos cl-mpm/particle:mp-position))
      mp
    (funcall func pos)))

(declaim (ftype (function (cl-mpm/mesh::cell function) (values)) calculate-val-cell))
(defun calculate-val-cell (cell func)
  (with-accessors ((pos cl-mpm/mesh::cell-centroid))
      cell
    (funcall func pos)))

(defun calculate-virtual-stress-mp (mp datum)
  (with-accessors ((pos cl-mpm/particle:mp-position))
      mp
    (pressure-virtual-stress)
    ;; (buoyancy-virtual-stress (tref pos 1 0) datum))
  ))

(defun calculate-virtual-stress-cell (cell datum)
  (with-accessors ((pos cl-mpm/mesh::cell-centroid))
      cell
    (pressure-virtual-stress)
    ;; (buoyancy-virtual-stress (tref pos 1 0) datum)
    ))

(defun calculate-virtual-divergance-mp (mp datum)
  (with-accessors ((pos cl-mpm/particle:mp-position))
      mp
    (pressure-virtual-div)
    ;; (buoyancy-virtual-div (tref pos 1 0) datum)
    ))

(defun calculate-virtual-divergance-cell (cell datum)
  (with-accessors ((pos cl-mpm/mesh::cell-centroid))
      cell
    (pressure-virtual-div)
    ;; (buoyancy-virtual-div (tref pos 1 0) datum)
    ))

(defun mps-in-cell (mesh mps)
  (loop for mp across mps
        do
           (let* ((id (cl-mpm/mesh::position-to-index
                      mesh
                      (cl-mpm/particle:mp-position mp)
                      #'floor))
                  (cell (cl-mpm/mesh::get-cell mesh id)))
             (incf (cl-mpm/mesh::cell-mp-count cell) 1))))

(defun apply-force-mps (mesh mps func-stress func-div)
  "Update force on nodes, with virtual stress field from mps"
  (declare (function func-stress func-div))
  (lparallel:pdotimes (i (length mps))
    (let ((mp (aref mps i)))

      (with-accessors ((volume cl-mpm/particle:mp-volume))
          mp
        ;;Iterate over neighbour nodes
        (cl-mpm::iterate-over-neighbours
         mesh mp
         (lambda (mesh mp node svp grads)
           (with-accessors ((node-force cl-mpm/mesh:node-force)
                            (node-lock  cl-mpm/mesh:node-lock)
                            (node-boundary cl-mpm/mesh::node-boundary-node)
                            (node-active  cl-mpm/mesh:node-active))
               node
             (when (and node-active node-boundary)
               ;;Lock node for multithreading
               (sb-thread:with-mutex (node-lock)

                 ;;Add the gradient of stress
                 (setf node-force (magicl:.+
                                   node-force
                                   (magicl:scale!
                                    (magicl:@
                                     (magicl:transpose
                                      (cl-mpm/shape-function::assemble-dsvp-2d grads))
                                     (funcall func-stress mp)
                                     )
                                    (* volume))))
                 ;;Add the divergance
                 (setf node-force (magicl:.+
                                   node-force
                                   (magicl:scale!
                                    (funcall func-div mp)
                                    (* svp volume))))
                 )))))))))

(defun apply-force-cells (mesh func-stress func-div)
  "Update force on nodes, with virtual stress field from cells"
  (declare (function func-stress func-div))
  (let ((cells (cl-mpm/mesh::mesh-cells mesh))
        (n-vol (expt (cl-mpm/mesh:mesh-resolution mesh) 2)))

    (lparallel:pdotimes (i (array-total-size cells))
             (let ((cell (row-major-aref cells i)))
               (let ((nodal-volume 0d0))
                 ;;Possibly clip ill posed cells
                 (when t;(> (/ nodal-volume (cl-mpm/mesh::cell-volume cell)) 1d-5)
                   ;;Iterate over a cells nodes
                   (cl-mpm/mesh::cell-iterate-over-neighbours
                    mesh cell
                    (lambda (mesh cell pos volume node svp grads)
                      (with-accessors ((node-force cl-mpm/mesh:node-force)
                                       (node-lock  cl-mpm/mesh:node-lock)
                                       (node-active cl-mpm/mesh:node-active)
                                       (node-boundary cl-mpm/mesh::node-boundary-node)
                                       (node-volume cl-mpm/mesh::node-volume)
                                       )
                          node
                        (when (and node-active node-boundary)
                          ;;Lock node
                          (sb-thread:with-mutex (node-lock)
                            ;;Subtract gradient of stress from node force
                            (setf node-force (magicl:.+
                                              node-force
                                              (magicl:scale!
                                               (magicl:@
                                                (magicl:transpose
                                                 (cl-mpm/shape-function::assemble-dsvp-2d grads))
                                                (funcall func-stress pos))
                                               (* -1d0 volume))))
                            ;;Subtract stress divergance from node force
                            (setf node-force (magicl:.+
                                              node-force
                                              (magicl:scale!
                                               (funcall func-div pos)
                                               (* -1d0 svp volume))))
                            )))
                      ))))))))
(defun direct-mp-enforcment (mesh mps datum)
  (lparallel:pdotimes (i (length mps))
    (let ((mp (aref mps i)))
      (with-accessors ((volume cl-mpm/particle:mp-volume)
                       (g cl-mpm/particle:mp-gravity)
                       (pos cl-mpm/particle:mp-position))
          mp
        (cl-mpm::iterate-over-neighbours
         mesh mp
         (lambda (mesh mp node svp grads)
           (with-accessors ((node-force cl-mpm/mesh:node-force)
                            (node-lock  cl-mpm/mesh:node-lock)
                            (node-active  cl-mpm/mesh:node-active))
               node
             (when node-active
               (sb-thread:with-mutex (node-lock)
                 ;;External force
                 (let* ((rho 1000)
                        (h (- datum (magicl:tref pos 1 0)))
                        (f (* rho g volume))
                        )
                   (when (> h 0d0)
                     (incf (magicl:tref node-force 1 0) (* -1 volume g rho svp))))
                 )))))))))

(defun check-neighbour-cell (cell)
  "Check if neighbour cell has nodes in"
  (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                   (nodes cl-mpm/mesh::cell-nodes))
      cell
    (if (> mp-count 0)
        (progn
          (loop for n in nodes
                do
                   (setf (cl-mpm/mesh::node-boundary-node n) t))
          t)
        nil)))

;;For the MPM case
(defun populate-cell-mp-count (mesh mps)
  (let ((cells (cl-mpm/mesh::mesh-cells mesh)))
    (lparallel:pdotimes (i (array-total-size cells))
      (let ((cell (row-major-aref cells i)))
        (setf (cl-mpm/mesh::cell-mp-count cell) 0)))
    (loop for mp across mps
          do (with-accessors ((pos cl-mpm/particle:mp-position))
                 mp
               (let* ((id (cl-mpm/mesh:position-to-index mesh pos #'floor)))
                 (when (cl-mpm/mesh::in-bounds-cell mesh id)
                   (let ((cell (cl-mpm/mesh::get-cell mesh id)))
                     (incf (cl-mpm/mesh::cell-mp-count cell)))))))))

;;For the MPM case
(defun populate-cell-mp-count-gimp (mesh mps)
  (let ((cells (cl-mpm/mesh::mesh-cells mesh))
        (h (cl-mpm/mesh:mesh-resolution mesh)))
    (lparallel:pdotimes (i (array-total-size cells))
      (let ((cell (row-major-aref cells i)))
        (setf (cl-mpm/mesh::cell-mp-count cell) 0)))
    (loop for mp across mps
          do (with-accessors ((pos cl-mpm/particle:mp-position)
                              (size cl-mpm/particle::mp-domain-size))
                 mp
               (let* ((pmin (cl-mpm/mesh:position-to-index mesh (magicl:.- pos (magicl:scale size 0.5d0)) #'floor))
                      (pmax (cl-mpm/mesh:position-to-index mesh (magicl:.+ pos (magicl:scale size 0.5d0)) #'floor)))
                 (loop for x from (first pmin) to (first pmax)
                       do
                          (loop for y from (second pmin) to (second pmax)
                                do
                                   (let ((id (list x y)))
                                     (when (cl-mpm/mesh::in-bounds-cell mesh id)
                                       (let ((cell (cl-mpm/mesh::get-cell mesh id)))
                                         (incf (cl-mpm/mesh::cell-mp-count cell))))))))))))
(defun populate-cell-mp-count-volume (mesh mps)
  (let ((cells (cl-mpm/mesh::mesh-cells mesh)))
    (lparallel:pdotimes (i (array-total-size cells))
      (let ((cell (row-major-aref cells i)))
        (setf (cl-mpm/mesh::cell-mp-count cell) 0)))
    (lparallel:pdotimes (i (array-total-size cells))
      (let ((cell (row-major-aref cells i)))
        (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                         (neighbours cl-mpm/mesh::cell-neighbours)
                         (index cl-mpm/mesh::cell-index)
                         (nodes cl-mpm/mesh::cell-nodes)
                         )
            cell
          (when (every (lambda (n)
                         (> (cl-mpm/mesh:node-mass n) 0d0)
                         ) nodes)
            (setf mp-count 1))
          )))))

(defun locate-mps-cells (mesh mps)
  "Mark boundary nodes based on neighbour MP inclusion"
  (let ((cells (cl-mpm/mesh::mesh-cells mesh)))
    ;;The aproach described by the paper
    ;; (populate-cell-mp-count mesh mps)
    ;; (populate-cell-mp-count-gimp mesh mps)
    (populate-cell-mp-count-volume mesh mps)
    (lparallel:pdotimes (i (array-total-size cells))
      (let ((cell (row-major-aref cells i)))
        (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                         (neighbours cl-mpm/mesh::cell-neighbours)
                         (index cl-mpm/mesh::cell-index)
                         (nodes cl-mpm/mesh::cell-nodes)
                         )
            cell
          (when (= mp-count 0)
              (loop for neighbour in neighbours
                    do
                       (when (check-neighbour-cell neighbour)
                         (loop for n in nodes
                               do
                                  (when (cl-mpm/mesh:node-active n)
                                    (sb-thread:with-mutex ((cl-mpm/mesh:node-lock n))
                                      (setf (cl-mpm/mesh::node-boundary-node n) t))))))
              ))))))

(defun find-active-nodes (mesh mps)
  (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
         (vtotal (expt h 2)))
    (cl-mpm::iterate-over-nodes mesh
                                (lambda (node)
                                  (with-accessors ((node-volume cl-mpm/mesh::node-volume)
                                                   (vtotal cl-mpm/mesh::node-volume-true)
                                                   (active cl-mpm/mesh::node-active)
                                                   (boundary cl-mpm/mesh::node-boundary-node))
                                      node
                                      (when active
                                        (setf boundary t)
                                        ))))))

(defun apply-bouyancy (sim datum-true)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm::sim-mps))
      sim
    (with-accessors ((h cl-mpm/mesh:mesh-resolution))
        mesh
      (let ((datum (- datum-true (* 0 0.5d0 h))))
        (locate-mps-cells mesh mps)
        (apply-force-mps mesh mps datum)
        (apply-force-cells mesh datum)
        ;; (direct-mp-enforcment mesh mps datum-true)
        )))
  )

(defun apply-non-conforming-nuemann (sim func-stress func-div)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm::sim-mps))
      sim
    (with-accessors ((h cl-mpm/mesh:mesh-resolution))
        mesh
      (locate-mps-cells mesh mps)
      (apply-force-mps mesh mps
                       (lambda (mp) (calculate-val-mp mp func-stress))
                       (lambda (mp) (calculate-val-mp mp func-div)))
      (apply-force-cells mesh
                         func-stress
                         func-div
                         ;; (lambda (cell) (calculate-val-cell cell func-stress))
                         ;; (lambda (cell) (calculate-val-cell cell func-div))
                         )
      )))

(defun make-bc-pressure (sim pressure-x pressure-y)
  (make-instance 'cl-mpm/bc::bc-closure
                 :index '(0 0)
                 :func (lambda ()
                         (apply-non-conforming-nuemann
                          sim
                          (lambda (pos)
                            (pressure-virtual-stress pressure-x pressure-y))
                          (lambda (pos)
                            (pressure-virtual-div))
                          )
                         )))

(defun make-bc-buoyancy (sim datum rho)
  (make-instance 'cl-mpm/bc::bc-closure
                 :index '(0 0)
                 :func (lambda ()
                         (progn
                           (apply-non-conforming-nuemann
                            sim
                            (lambda (pos)
                              (buoyancy-virtual-stress (tref pos 1 0) datum rho))
                            (lambda (pos)
                              (buoyancy-virtual-div (tref pos 1 0) datum rho)))
                           (with-accessors ((mesh cl-mpm:sim-mesh)
                                            (mps cl-mpm:sim-mps))
                               sim
                             (loop for mp across mps
                                   do
                                      (cl-mpm::iterate-over-neighbours
                                       mesh mp
                                       (lambda (mesh mp node svp grads)
                                         (when t;(cl-mpm/mesh::node-boundary-node node)
                                           (with-accessors ((pos cl-mpm/particle:mp-position)
                                                            (pressure cl-mpm/particle::mp-pressure))
                                               mp
                                             (setf pressure (pressure-at-depth (tref pos 1 0) datum)))
                                           )))))))))
