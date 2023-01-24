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
(let ((datum-true (- 200))
      (rho-true (/ 1000.0d0 2d0))
      )
  (defun buoyancy-virtual-stress (z)
    (let* ((rho rho-true)
           (g 9.8)
           (datum datum-true)
           (h (- datum z))
           (f (* -1d0 rho g h))
           )
      (if (>= h 0d0)
          (magicl:from-list (list f f 0d0) '(3 1) :type 'double-float)
          ;; (magicl:zeros '(3 1))
          (magicl:zeros '(3 1))
          )))
  (defun buoyancy-virtual-div (z)
    (let* ((rho rho-true)
           (g 9.8)
           (datum datum-true)
           (h (- datum z))
           (f (* rho g))
           )
      (if (>= h 0d0)
          (magicl:from-list (list 0 f) '(2 1) :type 'double-float)
          ;; (magicl:zeros '(2 1))
          (magicl:zeros '(2 1))
          ))))

(defun calculate-virtual-stress-mp (mp)
  (with-accessors ((pos cl-mpm/particle:mp-position))
      mp
    (buoyancy-virtual-stress (tref pos 1 0)))
  )

(defun calculate-virtual-stress-cell (cell)
  (with-accessors ((pos cl-mpm/mesh::cell-centroid))
      cell
    (buoyancy-virtual-stress (tref pos 1 0))))

(defun calculate-virtual-divergance-mp (mp)
  (with-accessors ((pos cl-mpm/particle:mp-position))
      mp
    (buoyancy-virtual-div (tref pos 1 0)))
  )

(defun calculate-virtual-divergance-cell (cell)
  (with-accessors ((pos cl-mpm/mesh::cell-centroid))
      cell
    (buoyancy-virtual-div (tref pos 1 0))))

(defun mps-in-cell (mesh mps)
  (loop for mp across mps
        do
           (let* ((id (cl-mpm/mesh::position-to-index
                      mesh
                      (cl-mpm/particle:mp-position mp)
                      #'floor))
                  (cell (cl-mpm/mesh::get-cell mesh id)))
             (incf (cl-mpm/mesh::cell-mp-count cell) 1))))
(defun apply-force-mps (mesh mps)
  (lparallel:pdotimes (i (length mps))
    (let ((mp (aref mps i)))
      (with-accessors ((volume cl-mpm/particle:mp-volume))
          mp
        (cl-mpm::iterate-over-neighbours ;-shape-linear
         mesh mp
         (lambda (mesh mp node svp grads)
           (with-accessors ((node-force cl-mpm/mesh:node-force)
                            (node-lock  cl-mpm/mesh:node-lock)
                            (node-active  cl-mpm/mesh:node-active))
               node
             (when node-active
               (sb-thread:with-mutex (node-lock)
                 ;;External force
                 (cl-mpm/fastmath:fast-add
                  node-force
                  (magicl:scale!
                   (magicl:@
                    (magicl:transpose!
                     (cl-mpm/shape-function::assemble-dsvp-2d grads))
                    (calculate-virtual-stress-mp mp))
                   volume))
                                        ;Internal force
                 (cl-mpm/fastmath:fast-add
                  node-force
                  (magicl:scale!
                   (calculate-virtual-divergance-mp mp)
                   (* svp volume))))))))))))
(defun apply-force-cells (mesh)
  (let ((cells (cl-mpm/mesh::mesh-cells mesh))
        (n-vol (expt (cl-mpm/mesh:mesh-resolution mesh) 2)))
    (lparallel:pdotimes (i (array-total-size cells))
             (let ((cell (row-major-aref cells i)))
               (let ((nodal-volume 0d0))
                 (cl-mpm/mesh::cell-iterate-over-neighbours
                  mesh cell
                  (lambda (mesh cell node svp grads)
                    (with-accessors ((node-volume cl-mpm/mesh::node-volume))
                        node
                      (incf nodal-volume node-volume))))
                 ;;Clip out nodes with poor conditied nodes
                 (when (> (/ nodal-volume (cl-mpm/mesh::cell-volume cell)) 1d-3)
                   (cl-mpm/mesh::cell-iterate-over-neighbours
                    mesh cell
                    (lambda (mesh cell node svp grads)
                      (with-accessors ((node-force cl-mpm/mesh:node-force)
                                       (node-lock  cl-mpm/mesh:node-lock)
                                       (node-active cl-mpm/mesh:node-active)
                                       (node-volume cl-mpm/mesh::node-volume)
                                       )
                          node
                        (with-accessors ((volume cl-mpm/mesh::cell-volume))
                            cell
                          (when node-active
                            (when t;(< node-volume (* 0.9d0 n-vol))
                              (sb-thread:with-mutex (node-lock)
                                        ;Internal force
                                (cl-mpm/fastmath:fast-add
                                 node-force
                                 (magicl:scale!
                                  (magicl:@
                                   (magicl:transpose! (cl-mpm/shape-function::assemble-dsvp-2d grads))
                                   (calculate-virtual-stress-cell cell))
                                  (* -1d0 volume)))
                                        ;External force
                                (cl-mpm/fastmath:fast-add
                                 node-force
                                 (magicl:scale!
                                  (calculate-virtual-divergance-cell cell)
                                  (* -1d0 svp volume)))
                                )))))
                      ))))))))

;(defun locate-mps-cells (mesh mps)
;  (loop for mp across mps
;        do (with-accessors ((pos cl-mpm/particle:mp-position))
;               mp
;             (let* ((id (cl-mpm/mesh:position-to-index mesh pos #'floor))
;                    (cell (cl-mpm/mesh::get-cell mesh id)))
;               (incf (cl-mpm/mesh::cell-mp-count cell))))))

(defun find-active-nodes (mesh mps)

  )

(defun apply-bouyancy (mesh mps)
  (apply-force-mps mesh mps)
  (apply-force-cells mesh))

