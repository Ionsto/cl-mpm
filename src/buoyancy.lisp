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
           (g 9.8)
           (datum datum-true)
           (h (- datum z))
           (f (* -1d0 rho g h))
           )
      (if (> h 0d0)
          f
          0d0)
      ))
  (defun buoyancy-virtual-stress (z datum-true)
    (let* ((rho rho-true)
           (g 9.8)
           (datum datum-true)
           (h (- datum z))
           (f (* -1d0 rho g h))
           (p -1d3)
           )
      (if (> h 0d0)
          (magicl:from-list (list f f 0d0) '(3 1) :type 'double-float)
          ;; (magicl:zeros '(3 1))
          (magicl:zeros '(3 1))
          )
      ;; (magicl:from-list (list p p 0d0)
      ;;                   '(3 1)
      ;;                   :type 'double-float)
      ))
  (defun buoyancy-virtual-div (z datum-true)
    (let* ((rho rho-true)
           (g 9.8)
           (datum datum-true)
           (h (- datum z))
           (f (* rho g))
           )
      (if (> h 0d0)
          (magicl:from-list (list 0 f) '(2 1) :type 'double-float)
          ;; (magicl:zeros '(2 1))
          (magicl:zeros '(2 1))
          )
      ;; (magicl:zeros '(2 1))

      )))

(defun calculate-virtual-stress-mp (mp datum)
  (with-accessors ((pos cl-mpm/particle:mp-position))
      mp
    (buoyancy-virtual-stress (tref pos 1 0) datum))
  )

(defun calculate-virtual-stress-cell (cell datum)
  (with-accessors ((pos cl-mpm/mesh::cell-centroid))
      cell
    (buoyancy-virtual-stress (tref pos 1 0) datum)))

(defun calculate-virtual-divergance-mp (mp datum)
  (with-accessors ((pos cl-mpm/particle:mp-position))
      mp
    (buoyancy-virtual-div (tref pos 1 0) datum))
  )

(defun calculate-virtual-divergance-cell (cell datum)
  (with-accessors ((pos cl-mpm/mesh::cell-centroid))
      cell
    (buoyancy-virtual-div (tref pos 1 0) datum)))

(defun mps-in-cell (mesh mps)
  (loop for mp across mps
        do
           (let* ((id (cl-mpm/mesh::position-to-index
                      mesh
                      (cl-mpm/particle:mp-position mp)
                      #'floor))
                  (cell (cl-mpm/mesh::get-cell mesh id)))
             (incf (cl-mpm/mesh::cell-mp-count cell) 1))))

(defun apply-force-mps (mesh mps datum)
  (lparallel:pdotimes (i (length mps))
    (let ((mp (aref mps i)))
      (with-accessors ((volume cl-mpm/particle:mp-volume))
          mp
        (cl-mpm::iterate-over-neighbours ;-shape-linear
         mesh mp
         (lambda (mesh mp node svp grads)
           (with-accessors ((node-force cl-mpm/mesh:node-force)
                            (node-lock  cl-mpm/mesh:node-lock)
                            (node-boundary cl-mpm/mesh::node-boundary-node)
                            (node-active  cl-mpm/mesh:node-active))
               node
             (when (and node-active node-boundary)
               (with-accessors ((pos cl-mpm/particle:mp-position)
                                (pressure cl-mpm/particle::mp-pressure))
                   mp
                 (setf pressure (pressure-at-depth (tref pos 1 0) datum)))
               (sb-thread:with-mutex (node-lock)
                 (setf node-force (magicl:.+
                                   node-force
                                   (magicl:scale!
                                    (magicl:@
                                     (magicl:transpose
                                      (cl-mpm/shape-function::assemble-dsvp-2d grads))
                                     (calculate-virtual-stress-mp mp datum))
                                    (* volume))))
                 (setf node-force (magicl:.+
                                   node-force
                                   (magicl:scale!
                                    (calculate-virtual-divergance-mp mp datum)
                                    (* svp volume))))
                 ;;External force
                 ;; (cl-mpm/fastmath:fast-add
                 ;;  node-force
                 ;;  (magicl:scale
                 ;;   (magicl:@
                 ;;    (magicl:transpose!
                 ;;     (cl-mpm/shape-function::assemble-dsvp-2d grads))
                 ;;    (calculate-virtual-stress-mp mp datum))
                 ;;   volume))
                                        ;Internal force
                 ;; (cl-mpm/fastmath:fast-add
                 ;;  node-force
                 ;;  (magicl:scale!
                 ;;   (calculate-virtual-divergance-mp mp datum)
                 ;;   (* svp volume)))
                 )))))))))

(defun apply-force-cells (mesh datum)
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
                 (when t;(> (/ nodal-volume (cl-mpm/mesh::cell-volume cell)) 1d-5)
                   (cl-mpm/mesh::cell-iterate-over-neighbours
                    mesh cell
                    (lambda (mesh cell node svp grads)
                      (with-accessors ((node-force cl-mpm/mesh:node-force)
                                       (node-lock  cl-mpm/mesh:node-lock)
                                       (node-active cl-mpm/mesh:node-active)
                                       (node-boundary cl-mpm/mesh::node-boundary-node)
                                       (node-volume cl-mpm/mesh::node-volume)
                                       )
                          node
                        (with-accessors ((volume cl-mpm/mesh::cell-volume))
                            cell
                          (when (and node-active node-boundary)
                            (when t;(< node-volume (* 0.99d0 n-vol))
                              (sb-thread:with-mutex (node-lock)
                                        ;Internal force
                                (setf node-force (magicl:.+
                                                  node-force
                                                  (magicl:scale!
                                                   (magicl:@
                                                    (magicl:transpose
                                                     (cl-mpm/shape-function::assemble-dsvp-2d grads))
                                                    (calculate-virtual-stress-cell cell datum))
                                                   (* -1d0 volume))))
                                (setf node-force (magicl:.+
                                                  node-force
                                                  (magicl:scale!
                                                   (calculate-virtual-divergance-cell cell datum)
                                                   (* -1d0 svp volume))))
                                ;; (cl-mpm/fastmath:fast-add
                                ;;  node-force
                                ;;  (magicl:scale!
                                ;;   (magicl:@
                                ;;    (magicl:transpose!
                                ;;     (cl-mpm/shape-function::assemble-dsvp-2d grads))
                                ;;    (calculate-virtual-stress-cell cell datum))
                                ;;   (* -1d0 volume)))
                                        ;External force
                                ;; (cl-mpm/fastmath:fast-add
                                ;;  node-force
                                ;;  (magicl:scale!
                                ;;   (calculate-virtual-divergance-cell cell datum)
                                ;;   (* -1d0 svp volume)))
                                )))))
                      ))))))))
(defun direct-mp-enforcment (mesh mps datum)
  (lparallel:pdotimes (i (length mps))
    (let ((mp (aref mps i)))
      (with-accessors ((volume cl-mpm/particle:mp-volume)
                       (g cl-mpm/particle:mp-gravity)
                       (pos cl-mpm/particle:mp-position))
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
                 (let* ((rho 1000)
                        (h (- datum (magicl:tref pos 1 0)))
                        (f (* rho g volume))
                        )
                   (when (> h 0d0)
                     (incf (magicl:tref node-force 1 0) (* -1 volume g rho svp))))
                 )))))))))

(defun check-neighbour-cell (cell)
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
(defun locate-mps-cells (mesh mps)
  (let ((cells (cl-mpm/mesh::mesh-cells mesh)))
    (lparallel:pdotimes (i (array-total-size cells))
      (let ((cell (row-major-aref cells i)))
        (setf (cl-mpm/mesh::cell-mp-count cell) 0)))
    (loop for mp across mps
          do (with-accessors ((pos cl-mpm/particle:mp-position))
                 mp
               (let* ((id (cl-mpm/mesh:position-to-index mesh pos #'floor))
                      (cell (cl-mpm/mesh::get-cell mesh id)))
                 (incf (cl-mpm/mesh::cell-mp-count cell)))))
    (lparallel:pdotimes (i (array-total-size cells))
      (let ((cell (row-major-aref cells i)))
        (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                         (neighbours cl-mpm/mesh::cell-neighbours)
                         (index cl-mpm/mesh::cell-index)
                         (nodes cl-mpm/mesh::cell-nodes)
                         )
            cell
          (when (= mp-count 0)
              ;; (loop for n in nodes
              ;;       do
              ;;          (sb-thread:with-mutex ((cl-mpm/mesh:node-lock n))
              ;;            (setf (cl-mpm/mesh::node-boundary-node n) t)))
            ;; (let ((patch 1))
            ;;   (loop for dx from (- patch) to patch
            ;;         do
            ;;            (loop for dy from (- patch) to patch
            ;;                  do
            ;;                     (unless (and (= dx 0) (= dy 0))
            ;;                       (let ((di (mapcar #'+ index (list dx dy))))
            ;;                         (when (cl-mpm/mesh::in-bounds-cell mesh di)
            ;;                           (check-neighbour-cell (apply #'aref cells di))))))))
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
                                        ;; (when (and ;(> (/ node-volume vtotal) 1e-3)
                                        ;;            (< (/ node-volume vtotal) 0.99))
                                        ;;   (setf boundary t))
                                        ))))))

(defun apply-bouyancy (sim datum-true)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm::sim-mps))
      sim
    (with-accessors ((h cl-mpm/mesh:mesh-resolution))
        mesh
      (let ((datum (- datum-true (* 0 0.5d0 h))))
        (locate-mps-cells mesh mps)
        ;; (find-active-nodes mesh mps)
        (apply-force-mps mesh mps datum)
        (apply-force-cells mesh datum)
        ;; (direct-mp-enforcment mesh mps datum-true)
        ))))

