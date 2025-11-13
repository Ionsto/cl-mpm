(defpackage :cl-mpm/buoyancy
  (:use :cl
        :cl-mpm
        :cl-mpm/particle
        :cl-mpm/mesh
        :cl-mpm/utils
        :cl-mpm/fastmaths)
  (:import-from
    :magicl tref .+ .-)
  (:export
   ))
;; (declaim (optimize (debug 3) (safety 3) (speed 0)))
(declaim (optimize (debug 0) (safety 0) (speed 3)))

;    #:make-shape-function
(in-package :cl-mpm/buoyancy)

;(defgeneric virtual-stress ())
(defun pressure-at-depth (z datum-true rho g)
  (let* (;; (g -9.8d0)
         (datum datum-true)
         (h (- datum z))
         (f (* rho g h))
         )
    (if (> h 0d0)
        f
        0d0)
    ))
(declaim (notinline buoyancy-virtual-stress))
(defun buoyancy-virtual-stress (z datum-true rho g)
  (let* (;; (g -9.8d0)
         (datum datum-true)
         (h (- datum z))
         (f (* 1d0 rho g h))
         )
    (if (> h 0d0)
        (voigt-from-list (list f f f 0d0 0d0 0d0))
        (voigt-zeros))))

(defun buoyancy-virtual-div (z datum-true rho g)
  (let* (;; (g -9.8d0)
         (datum datum-true)
         (h (- datum z))
         (f (* -1d0 rho g))
         )
    ;; (vector-from-list (list 0d0 f))
    (if (> h 0d0)
        (vector-from-list (list 0d0 f 0d0))
        (vector-zeros)
        )
    ;; (vector-zeros)
    ))

(defun pressure-virtual-stress (pressure-x pressure-y)
  (voigt-from-list (list pressure-x pressure-y 0d0 0d0 0d0 0d0))
  ;; (voigt-from-list (list 0d0 0d0 0d0 0d0 0d0 0d0))
  )

(defun pressure-virtual-div ()
  (vector-zeros)
  ;; (vector-from-list (list 0d0 0d0 0d0))
  )

;; (defun pressure-virtual-stress ()
;;   (let* ((p -1d3)
;;          )
;;     (magicl:from-list (list p p 0d0)
;;                       '(3 1)
;;                       :type 'double-float)
;;     ))
;; (defun pressure-virtual-div ()
;;   (magicl:zeros '(2 1)))

(defun compute-mp-displacement (mesh mp)
  ;; (cl-mpm/fastmaths:fast-.+ corner (cl-mpm/particle::mp-displacement-increment mp))
  (with-accessors ((disp-inc cl-mpm/particle::mp-displacement-increment))
      mp
    (fast-zero disp-inc)
    (cl-mpm:iterate-over-neighbours
     mesh
     mp
     (lambda (mesh mp node svp grads fsvp fgrads)
       (declare
        (ignore mp mesh fsvp fgrads)
        (cl-mpm/mesh::node node)
        (cl-mpm/particle:particle mp)
        (double-float svp))
       (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
                        (node-acc cl-mpm/mesh:node-acceleration)
                        (node-disp cl-mpm/mesh::node-displacment)
                        (node-scalar cl-mpm/mesh::node-boundary-scalar)
                        (node-active cl-mpm/mesh:node-active)
                        ) node
         (declare (double-float node-scalar)
                  (boolean node-active))
         (when node-active
           (cl-mpm/fastmaths::fast-fmacc disp-inc node-disp svp)))))
    (fast-.+ (cl-mpm/particle::mp-position mp)
            disp-inc
            (cl-mpm/particle::mp-position-trial mp))))

(declaim (ftype (function (cl-mpm/particle:particle function) (values)) calculate-val-mp))
(defun calculate-val-mp (mp func)
  (let ((pos (fast-.+
              (cl-mpm/particle::mp-position mp)
              (cl-mpm/particle::mp-displacement-increment mp)
              )))
    (funcall func pos)))

(declaim (ftype (function (cl-mpm/mesh::cell function) (values)) calculate-val-cell))
(defun calculate-val-cell (cell func)
  (with-accessors ((pos cl-mpm/mesh::cell-centroid))
      cell
    (funcall func pos)))

;; (defun calculate-virtual-stress-mp (mp datum)
;;   (with-accessors ((pos cl-mpm/particle:mp-position))
;;       mp
;;     (pressure-virtual-stress)
;;     ;; (buoyancy-virtual-stress (tref pos 1 0) datum))
;;   ))

;; (defun calculate-virtual-stress-cell (cell datum)
;;   (with-accessors ((pos cl-mpm/mesh::cell-centroid))
;;       cell
;;     (pressure-virtual-stress)
;;     ;; (buoyancy-virtual-stress (tref pos 1 0) datum)
;;     ))

;; (defun calculate-virtual-divergance-mp (mp datum)
;;   (with-accessors ((pos cl-mpm/particle:mp-position))
;;       mp
;;     (pressure-virtual-div)
;;     ;; (buoyancy-virtual-div (tref pos 1 0) datum)
;;     ))

;; (defun calculate-virtual-divergance-cell (cell datum)
;;   (with-accessors ((pos cl-mpm/mesh::cell-centroid))
;;       cell
;;     (pressure-virtual-div)
;;     ;; (buoyancy-virtual-div (tref pos 1 0) datum)
;;     ))

(defun mps-in-cell (mesh mps)
  (loop for mp across mps
        do
           (let* ((id (cl-mpm/mesh::position-to-index
                      mesh
                      (cl-mpm/particle:mp-position mp)
                      #'floor))
                  (cell (cl-mpm/mesh::get-cell mesh id)))
             (incf (cl-mpm/mesh::cell-mp-count cell) 1))))

(defun melt-rate (pos)
  1d0
  ;; (if (< (magicl:tref pos 1 0) 300)
  ;;     ;(exp (* 0.1 (- (magicl:tref pos 1 0) 300)))
  ;;     ;; (+ 1 (exp (* 0.05 (- (magicl:tref pos 1 0) 300))))
  ;;     1d0
  ;;     ;; (/ 1 (+ 1 (expt
  ;;     ;;            (abs (min 0 (- (magicl:tref pos 1 0) 300))) 2)))
  ;;     0d0)
  )



(defun direct-mp-enforcment (mesh mps datum)
  (lparallel:pdotimes (i (length mps))
    (let ((mp (aref mps i)))
      (with-accessors ((volume cl-mpm/particle:mp-volume)
                       (g cl-mpm/particle:mp-gravity)
                       (pos cl-mpm/particle:mp-position))
          mp
        (cl-mpm::iterate-over-neighbours
         mesh mp
         (lambda (mesh mp node svp grads fsvp fgrads)
           (with-accessors ((node-force cl-mpm/mesh::node-external-force)
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
                   (active cl-mpm/mesh::cell-active)
                   (partial cl-mpm/mesh::cell-partial)
                   (nodes cl-mpm/mesh::cell-nodes))
      cell
    (if (and active (not partial))
     ;(> mp-count 0)
        (progn
          ;; (loop for n in nodes
          ;;       do (setf (cl-mpm/mesh::node-boundary-node n) t))
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
    (cl-mpm::iterate-over-cells
     mesh
     (lambda (cell)
       (setf (cl-mpm/mesh::cell-mp-count cell) 0)))
    (loop for mp across mps
          do (with-accessors ((pos cl-mpm/particle:mp-position)
                              (size cl-mpm/particle::mp-domain-size))
                 mp
               (let* ((pmin (cl-mpm/mesh:position-to-index mesh (magicl:.- pos (magicl:scale size 0.5d0)) #'floor))
                      (pmax (cl-mpm/mesh:position-to-index mesh (cl-mpm/fastmaths::fast-.+ pos (magicl:scale size 0.5d0)) #'floor)))
                 (loop for x from (first pmin) to (first pmax)
                       do
                          (loop for y from (second pmin) to (second pmax)
                                do
                                   (let ((id (list x y 0)))
                                     (when (cl-mpm/mesh::in-bounds-cell mesh id)
                                       (let ((cell (cl-mpm/mesh::get-cell mesh id)))
                                         (incf (cl-mpm/mesh::cell-mp-count cell))))))))))))
(defun populate-cell-mp-count-volume (mesh mps clip-func)
  (let ((cells (cl-mpm/mesh::mesh-cells mesh)))
    (cl-mpm::iterate-over-cells
     mesh
     (lambda (cell)
       (setf (cl-mpm/mesh::cell-mp-count cell) 0)))
    (cl-mpm::iterate-over-cells
     mesh
     (lambda (cell)
       (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                        (neighbours cl-mpm/mesh::cell-neighbours)
                        (index cl-mpm/mesh::cell-index)
                        (nodes cl-mpm/mesh::cell-nodes)
                        (centroid cl-mpm/mesh::cell-centroid))
           cell
         (when (funcall clip-func centroid)
           (when (every
                  (lambda (n)
                    (and (cl-mpm/mesh:node-active n)
                         (> (cl-mpm/mesh:node-mass n) 0d0)))
                  nodes)
             (setf mp-count 1))))))))

(defun cell-clipping (pos datum)
  (<= (magicl:tref pos 1 0) datum)
  ;; (<= (magicl:tref pos 1 0) 300)
  ;; t
  )
(defgeneric populate-cells-volume (sim clip-function))

(defmethod populate-cells-volume ((sim cl-mpm:mpm-sim) clip-function)
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (let ((cells (cl-mpm/mesh::mesh-cells mesh)))
      (cl-mpm::iterate-over-cells
       mesh
       (lambda (cell)
         (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                          (neighbours cl-mpm/mesh::cell-neighbours)
                          (index cl-mpm/mesh::cell-index)
                          (nodes cl-mpm/mesh::cell-nodes)
                          (pruned cl-mpm/mesh::cell-pruned)
                          (boundary cl-mpm/mesh::cell-boundary)
                          (pos cl-mpm/mesh::cell-centroid)
                          (vt cl-mpm/mesh::cell-volume)
                          (active cl-mpm/mesh::cell-active)
                          )
             cell
           (setf boundary nil)
           (when t
             (flet ((check-cell (c)
                      (with-accessors ((pos cl-mpm/mesh::cell-centroid)
                                       (neighbours cl-mpm/mesh::cell-neighbours)
                                       (vt cl-mpm/mesh::cell-volume)
                                       (nns cl-mpm/mesh::cell-nodes)
                                       )
                          c
                        (when (and (funcall clip-function pos))
                          (let ((vest 0d0))
                            (loop for n in nns
                                  do
                                     (when (cl-mpm/mesh:node-active n)
                                       (incf vest
                                             (* 0.25d0 (/
                                                        (cl-mpm/mesh::node-volume n)
                                                        (cl-mpm/mesh::node-volume-true n))))))
                            (when (< vest 0.5d0)
                              (setf boundary t)
                              (loop for n in nodes
                                    do
                                       (when (cl-mpm/mesh:node-active n)
                                         (sb-thread:with-mutex ((cl-mpm/mesh:node-lock n))
                                           (setf (cl-mpm/mesh::node-boundary-node n) t))))
                              ))))))
               (check-cell cell)
               (loop for neighbour in neighbours
                     do (check-cell neighbour)))))))
      )))

(defmethod populate-cells-volume ((sim cl-mpm/mpi::mpm-sim-mpi) clip-function)
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (let ((cells (cl-mpm/mesh::mesh-cells mesh)))
      (cl-mpm::iterate-over-cells
       mesh
       (lambda (cell)
         (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                          (neighbours cl-mpm/mesh::cell-neighbours)
                          (index cl-mpm/mesh::cell-index)
                          (nodes cl-mpm/mesh::cell-nodes)
                          (pruned cl-mpm/mesh::cell-pruned)
                          (boundary cl-mpm/mesh::cell-boundary)
                          (pos cl-mpm/mesh::cell-centroid)
                          (vt cl-mpm/mesh::cell-volume)
                          (active cl-mpm/mesh::cell-active)
                          )
             cell
           (setf boundary nil)
           (setf active (cl-mpm/mpi::in-computational-domain sim pos))
           (when active
             (flet ((check-cell (c)
                      (with-accessors ((pos cl-mpm/mesh::cell-centroid)
                                       (neighbours cl-mpm/mesh::cell-neighbours)
                                       (vt cl-mpm/mesh::cell-volume)
                                       (nns cl-mpm/mesh::cell-nodes)
                                       )
                          c
                        (when (and (funcall clip-function pos))
                          (let ((vest 0d0))
                            (loop for n in nns
                                  do
                                     (when (cl-mpm/mesh:node-active n)
                                       (incf vest
                                             (* 0.25d0 (/
                                                        (cl-mpm/mesh::node-volume n)
                                                        (cl-mpm/mesh::node-volume-true n))))))
                            (when (< vest 0.5d0)
                              (setf boundary t)
                              (loop for n in nodes
                                    do
                                       (when (cl-mpm/mesh:node-active n)
                                         (sb-thread:with-mutex ((cl-mpm/mesh:node-lock n))
                                           (setf (cl-mpm/mesh::node-boundary-node n) t))))
                              ))))))
               (check-cell cell)
               (loop for neighbour in neighbours
                     do (check-cell neighbour)))))))
      )))

;; (defmethod populate-cells-volume ((sim cl-mpm/mpi:mpm-sim-mpi) clip-function)
;;   (with-accessors ((mesh cl-mpm:sim-mesh))
;;       sim
;;     (let ((cells (cl-mpm/mesh::mesh-cells mesh)))
;;       (cl-mpm::iterate-over-cells
;;        mesh
;;        (lambda (cell)
;;          (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
;;                           (neighbours cl-mpm/mesh::cell-neighbours)
;;                           (index cl-mpm/mesh::cell-index)
;;                           (nodes cl-mpm/mesh::cell-nodes)
;;                           (pruned cl-mpm/mesh::cell-pruned)
;;                           (boundary cl-mpm/mesh::cell-boundary)
;;                           (pos cl-mpm/mesh::cell-centroid)
;;                           (vt cl-mpm/mesh::cell-volume)
;;                           )
;;              cell
;;            (setf boundary nil)
;;            (flet ((check-cell (c)
;;                     (with-accessors ((pos cl-mpm/mesh::cell-centroid)
;;                                      (neighbours cl-mpm/mesh::cell-neighbours)
;;                                      (vt cl-mpm/mesh::cell-volume)
;;                                      (nns cl-mpm/mesh::cell-nodes)
;;                                      )
;;                         c
;;                       (when (and (funcall clip-function pos))
;;                         (let ((vest 0d0))
;;                           (loop for n in nns
;;                                 do
;;                                    (when (cl-mpm/mesh:node-active n)
;;                                      (incf vest
;;                                            (* 0.25d0 (/
;;                                                       (cl-mpm/mesh::node-volume n)
;;                                                       (cl-mpm/mesh::node-volume-true n))))))
;;                           (when (< vest 0.5d0)
;;                             (setf boundary t)
;;                             (loop for n in nodes
;;                                   do
;;                                      (when (cl-mpm/mesh:node-active n)
;;                                        (sb-thread:with-mutex ((cl-mpm/mesh:node-lock n))
;;                                          (setf (cl-mpm/mesh::node-boundary-node n) t))))
;;                             ))

;;                         ))))
;;              (check-cell cell)
;;              (loop for neighbour in neighbours
;;                    do (check-cell neighbour))))))
;;       )))
(defun populate-nodes-volume (mesh clip-function)
  (cl-mpm::iterate-over-nodes
   mesh
   (lambda (node)
     (with-accessors ((node-lock cl-mpm/mesh::node-lock)
                      (volume cl-mpm/mesh::node-volume)
                      (vt cl-mpm/mesh::node-volume-true)
                      (boundary cl-mpm/mesh::node-boundary-node)
                      (active cl-mpm/mesh::node-active)
                      (pos cl-mpm/mesh::node-position)
                      )
         node
       (when (< volume (* 0.95d0 vt))
         (when (and (funcall clip-function pos))
           (when active
             (setf boundary t))))))))

(defun populate-nodes-volume-damage (mesh clip-function)
  (cl-mpm::iterate-over-nodes
   mesh
   (lambda (node)
     (with-accessors ((node-lock cl-mpm/mesh::node-lock)
                      (volume cl-mpm/mesh::node-volume)
                      (damage cl-mpm/mesh::node-damage)
                      (vt cl-mpm/mesh::node-volume-true)
                      (boundary cl-mpm/mesh::node-boundary-node)
                      (active cl-mpm/mesh::node-active)
                      (pos cl-mpm/mesh::node-position)
                      )
         node
       (when (< (* volume (- 1d0 damage)) (* 0.95d0 vt))
         (when (and (funcall clip-function pos))
           (when active
             (setf boundary t))))))))

(defun populate-nodes-domain (mesh clip-function)
  (cl-mpm::iterate-over-nodes
   mesh
   (lambda (node)
     (with-accessors ((node-lock cl-mpm/mesh::node-lock)
                      (volume cl-mpm/mesh::node-volume)
                      (vt cl-mpm/mesh::node-volume-true)
                      (boundary cl-mpm/mesh::node-boundary-node)
                      (active cl-mpm/mesh::node-active)
                      (pos cl-mpm/mesh::node-position)
                      )
         node
       (when (and active
                  (funcall clip-function pos))
         (setf boundary t))))))

(defun set-boundary (cell)
  (with-accessors ((nodes cl-mpm/mesh::cell-nodes)
                   (boundary cl-mpm/mesh::cell-boundary))
      cell
    (setf boundary t)
    (loop for n in nodes
          do
             (when (cl-mpm/mesh:node-active n)
               (sb-thread:with-mutex ((cl-mpm/mesh:node-lock n))
                 (setf (cl-mpm/mesh::node-boundary-node n) t))))))

(defgeneric locate-mps-cells (sim clip-function))
(defmethod locate-mps-cells (sim clip-function)
  "mark boundary nodes based on neighbour mp inclusion"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (let ((cells (cl-mpm/mesh::mesh-cells mesh)))
        ;;the aproach described by the paper
        ;; (populate-cell-mp-count mesh mps)
        ;; (populate-cell-mp-count-gimp mesh mps)
        ;; (populate-cell-mp-count-volume mesh mps clip-function)
        ;; (populate-cell-nodes mesh mps)
        ;; (prune-buoyancy-nodes mesh '(0 0) 300)
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
                            (boundary cl-mpm/mesh::cell-boundary)
                            (pos cl-mpm/mesh::cell-centroid))
               cell
             (setf boundary nil)
             ;; (when (and
             ;;        (funcall clip-function pos)
             ;;        active partial)
             ;;   (set-boundary cell))
             (when (and
                    active
                    partial
                    (funcall clip-function pos)
                    (not pruned))
               (set-boundary cell)
               ;; (setf boundary t)
               ;; (loop for n in (cl-mpm/mesh::cell-nodes cell)
               ;;       do (setf (cl-mpm/mesh::node-boundary-node n) t))
               (loop for neighbour in neighbours
                     do
                        (when (funcall clip-function (cl-mpm/mesh::cell-centroid neighbour))
                          (when (check-neighbour-cell neighbour)
                            ;; (setf (cl-mpm/mesh::cell-boundary neighbour) t)
                            (set-boundary neighbour)
                            ;; (setf boundary t)
                            ;; (loop for n in nodes
                            ;;       do
                            ;;          (when (cl-mpm/mesh:node-active n)
                            ;;            (sb-thread:with-mutex ((cl-mpm/mesh:node-lock n))
                            ;;              (setf (cl-mpm/mesh::node-boundary-node n) t))))
                            ))))))))))

(defmethod locate-mps-cells ((sim cl-mpm/mpi:mpm-sim-mpi) clip-function)
  "mark boundary nodes based on neighbour mp inclusion"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (let ((cells (cl-mpm/mesh::mesh-cells mesh)))
        ;;the aproach described by the paper
        ;; (populate-cell-mp-count mesh mps)
        ;; (populate-cell-mp-count-gimp mesh mps)
        ;; (populate-cell-mp-count-volume mesh mps clip-function)
        ;; (populate-cell-nodes mesh mps)
        ;; (prune-buoyancy-nodes mesh '(0 0) 300)
        (cl-mpm::iterate-over-cells
         mesh
         (lambda (cell)
           (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                            (neighbours cl-mpm/mesh::cell-neighbours)
                            (index cl-mpm/mesh::cell-index)
                            (centroid cl-mpm/mesh::cell-centroid)
                            (nodes cl-mpm/mesh::cell-nodes)
                            (pruned cl-mpm/mesh::cell-pruned)
                            (active cl-mpm/mesh::cell-active)
                            (boundary cl-mpm/mesh::cell-boundary)
                            (pos cl-mpm/mesh::cell-centroid))
               cell
             (setf boundary nil)
             (setf active (cl-mpm/mpi::in-computational-domain sim centroid))
             (when (and (= mp-count 0)
                        (funcall clip-function pos)
                        (not pruned))
               ;; (setf boundary t)
               ;; (loop for n in (cl-mpm/mesh::cell-nodes cell)
               ;;       do (setf (cl-mpm/mesh::node-boundary-node n) t))
               (loop for neighbour in neighbours
                     do
                        (when (funcall clip-function (cl-mpm/mesh::cell-centroid neighbour))
                          (when (check-neighbour-cell neighbour)
                            (setf boundary t)
                            (loop for n in nodes
                                  do
                                     (when (cl-mpm/mesh:node-active n)
                                       (sb-thread:with-mutex ((cl-mpm/mesh:node-lock n))
                                         (setf (cl-mpm/mesh::node-boundary-node n) t))))))))))))))




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
                                        (setf boundary t)))))))


(declaim (notinline apply-non-conforming-nuemann))
(defun apply-non-conforming-nuemann (sim func-stress func-div clip-function)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm::sim-mps))
      sim
    (with-accessors ((h cl-mpm/mesh:mesh-resolution))
        mesh
      (locate-mps-cells sim clip-function)
      ;; (populate-cells-volume sim clip-function)
      ;; (populate-nodes-volume mesh clip-function)
      ;; (populate-nodes-volume-damage mesh clip-function)
      ;; (populate-nodes-domain mesh clip-function)
      (apply-force-mps-3d mesh mps
                       (lambda (mp) (calculate-val-mp mp func-stress))
                       (lambda (mp) (calculate-val-mp mp func-div))
                       clip-function
                       )
      (apply-force-cells-3d mesh
                         func-stress
                         func-div
                         clip-function
                         )
      )))

(defclass bc-non-conforming-neumann (bc)
  ()
  (:documentation "A nonconforming neumann bc"))

(defclass bc-pressure (cl-mpm/bc::bc)
  ((sim
    :accessor bc-pressure-sim
    :initarg :sim)
   (pressures
    :accessor bc-pressure-pressures
    :initarg :pressures)
   (clip-func
    :accessor bc-pressure-clip-func
    :type function
    :initarg :clip-func
    :initform (lambda (&rest args) t)))
  (:documentation "A nonconforming pressure bc"))


(defmethod cl-mpm/bc::apply-bc ((bc bc-pressure) node mesh dt)
  "Arbitrary closure BC"
  (with-accessors ((pressures bc-pressure-pressures)
                   (clip-func bc-pressure-clip-func)
                   (sim bc-pressure-sim))
      bc
    (apply-non-conforming-nuemann
     sim
     (lambda (pos)
       (pressure-virtual-stress (first pressures) (second pressures)))
     (lambda (pos)
       (pressure-virtual-div))
     clip-func)))

(defclass bc-buoyancy (cl-mpm/bc::bc)
  ((sim
    :accessor bc-buoyancy-sim
    :initarg :sim
    :initform nil
    )
   (datum
    :accessor bc-buoyancy-datum
    :initarg :datum
    :initform 0d0
    )
   (rho
    :accessor bc-buoyancy-rho
    :initarg :rho
    :initform 0d0
    )
   (clip-func
    :accessor bc-buoyancy-clip-func
    :type function
    :initarg :clip-func
    :initform (lambda (&rest args) t)
    )
   (viscous-damping
    :accessor bc-viscous-damping
    :initarg :visc-damping
    :initform 0d0))
  (:documentation "A nonconforming buoyancy bc"))

(defclass bc-scalar (cl-mpm/bc::bc)
  ((sim
    :accessor bc-buoyancy-sim
    :initarg :sim
    :initform nil
    )
   (datum
    :accessor bc-buoyancy-datum
    :initarg :datum
    :initform 0d0
    )
   (scalar-func
    :accessor bc-scalar-func
    :initarg :scalar-func
    :initform (lambda (pos) 1d0)
    )
   (damage-volume
    :accessor bc-damage-volume
    :initarg :damage-volume
    :initform nil)
   (bc-enable
    :accessor bc-enable
    :initarg :enable
    :initform t)
   (clip-func
    :accessor bc-buoyancy-clip-func
    :type function
    :initarg :clip-func
    :initform (lambda (&rest args) t)))
  (:documentation "A nonconforming buoyancy bc"))

(defun make-bc-pressure (sim pressure-x pressure-y &key (clip-func (lambda (pos) t)))
  (make-instance 'bc-pressure
                 :index nil
                 :sim sim
                 :pressures (list (float pressure-x 0d0) (float pressure-y 0d0))
                 :clip-func clip-func))

(defun prune-buoyancy-nodes (mesh origin datum)
    (with-accessors ((size cl-mpm/mesh::mesh-count))
        mesh
      (let ((cell-size (mapcar #'round (mapcar #'- size '(1 1)))))
        (let ((fill-array (make-array cell-size :initial-element nil))
              (search-queue (list origin)))
          ;; (setf (apply #'aref fill-array origin) nil)
          (loop while search-queue
                do
                 (let* ((s (pop search-queue)))
                   (when (and (cl-mpm/mesh::in-bounds-cell mesh s)
                              (not (apply #'aref fill-array s))
                              )
                     (let ((cell (cl-mpm/mesh::get-cell mesh s)))
                       (when (and
                              (eq (cl-mpm/mesh::cell-mp-count cell) 0)
                              (< (magicl:tref (cl-mpm/mesh::cell-centroid cell) 1 0) datum))
                         (setf (apply #'aref fill-array s) t)
                         (push (mapcar #'+ s (list 1 0)) search-queue)
                         (push (mapcar #'+ s (list -1 0)) search-queue)
                         (push (mapcar #'+ s (list 0 1)) search-queue)
                         (push (mapcar #'+ s (list 0 -1)) search-queue)
                         )))))
          (loop for x from 0 below (nth 0 cell-size)
                do (loop for y from 0 below (nth 1 cell-size)
                         do
                            ;; (unless (aref fill-array x y))
                                        ;(cl-mpm/mesh::cell-nodes (cl-mpm/mesh::get-cell mesh (list x y)))
                            ;; (print (cl-mpm/mesh::get-cell mesh (list x y)))
                            (setf (cl-mpm/mesh::cell-pruned (cl-mpm/mesh::get-cell mesh (list x y)))
                                  (not (aref fill-array x y)))
                            ;; (loop for n in (cl-mpm/mesh::cell-nodes (cl-mpm/mesh::get-cell mesh (list x y)))
                            ;;       do (when (cl-mpm/mesh::node-boundary-node n)
                            ;;            ;; (format t "~F ~F ~%" x y)
                            ;;            (setf (cl-mpm/mesh::node-boundary-node n) nil)))
                         ))))))

(defun make-bc-buoyancy (sim datum rho)
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (with-accessors ((h cl-mpm/mesh:mesh-resolution))
        mesh
      (make-instance 'bc-buoyancy
                     :index nil
                     :sim sim
                     :clip-func (lambda (pos datum)
                                  (and
                                   t
                                   (cell-clipping pos datum)
                                   ;(>= (magicl:tref pos 1 0) (* 1 h))
                                   ))
                     :rho rho
                     :datum datum))))

(defun make-bc-buoyancy-clip (sim datum rho clip-func &key (visc-damping 0d0))
  (declare (function clip-func))
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (with-accessors ((h cl-mpm/mesh:mesh-resolution))
        mesh
      (make-instance 'bc-buoyancy
                     :index nil
                     :sim sim
                     :clip-func (lambda (pos datum)
                                  (and (funcall clip-func pos datum)
                                       (< (cl-mpm/utils::varef pos 1) datum)))
                     :rho rho
                     :visc-damping visc-damping
                     :datum datum))))

(defgeneric markup-cells-nodes (sim bc))
(defmethod markup-cells-nodes (sim (bc bc))
  )
(defmethod markup-cells-nodes (sim (bc cl-mpm/buoyancy::bc-buoyancy))
  (with-accessors ((datum bc-buoyancy-datum)
                   (clip-function bc-buoyancy-clip-func))
      bc
    (declare (function clip-function))
    ;; (populate-cells-volume sim (lambda (pos) (funcall clip-function pos datum)))
    (locate-mps-cells sim (lambda (pos) (funcall clip-function pos datum)))
    ))

(defgeneric apply-buoyancy (sim func-stress func-div clip-function datum))

(defmethod apply-buoyancy (sim func-stress func-div clip-function datum)
  (declare (function func-stress func-div clip-function))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm::sim-mps))
      sim
    (with-accessors ((h cl-mpm/mesh:mesh-resolution))
        mesh
      ;; (locate-mps-cells sim clip-function)
      ;; (populate-cells-volume sim (lambda (pos) (funcall clip-function pos datum)))
      ;; (locate-mps-cells sim (lambda (pos) (funcall clip-function pos datum)))
      ;; (markup-cells-nodes sim)
      ;; (locate-mps-cells sim clip-function)
      ;; (populate-nodes-volume mesh clip-function)
      ;; (populate-nodes-volume-damage mesh clip-function)
      ;; (populate-nodes-domain mesh clip-function)

      (apply-force-mps-3d mesh mps
                       (lambda (mp) (calculate-val-mp mp func-stress))
                       (lambda (mp) (calculate-val-mp mp func-div))
                       (lambda (pos) (funcall clip-function pos datum))
                       )
      (apply-force-cells-3d mesh
                         func-stress
                         func-div
                         (lambda (pos) (funcall clip-function pos datum)))
      )))


(in-package :cl-mpm/mpi)
(make-mpi-ser
 node-buoyancy
 ((index index cl-mpm/mesh::node-index)
  (int boundary-node (lambda (mp) (if (cl-mpm/mesh::node-boundary-node mp) 1 0)))
  (float scalar cl-mpm/mesh::node-boundary-scalar)
  ))
(in-package :cl-mpm/buoyancy)
(defgeneric exchange-bc-data (sim bc))
(defmethod exchange-bc-data ((sim cl-mpm:mpm-sim) bc))
(defmethod exchange-bc-data ((sim cl-mpm/mpi:mpm-sim-mpi) bc)
  (with-accessors ((mesh cl-mpm:sim-mesh))
   sim
   (cl-mpm/mpi::exchange-node-like
    sim
    #'cl-mpm/mpi::serialise-node-buoyancy
    #'cl-mpm/mpi::deserialise-node-buoyancy
    (lambda (node-list)
      (lparallel:pdotimes (i (length node-list))
        (let* ((mpi-node (aref node-list i))
               (index (cl-mpm/mpi::mpi-object-node-buoyancy-index mpi-node))
               (mpi-boundary (= 1 (cl-mpm/mpi::mpi-object-node-buoyancy-boundary-node mpi-node)))
               (node (cl-mpm/mesh:get-node mesh index)))
          (if node
              (with-accessors ((active cl-mpm/mesh:node-active)
                               (boundary cl-mpm/mesh::node-boundary-node)
                               (scalar cl-mpm/mesh::node-boundary-scalar)
                               )
                  node
                (declare (double-float scalar))
                (if (cl-mpm/mpi::in-computational-domain-buffer sim (cl-mpm/mesh::node-position node) 1)
                    (progn
                      (setf boundary (or
                                      mpi-boundary
                                      boundary))
                      (incf scalar (the double-float (cl-mpm/mpi::mpi-object-node-buoyancy-scalar mpi-node))))
                    (progn
                      (setf boundary mpi-boundary)
                      (setf scalar (the double-float (cl-mpm/mpi::mpi-object-node-buoyancy-scalar mpi-node))))
                    ))
              (error "Buoancy MPI exchange touched invalid node?"))))))))

(defmethod cl-mpm/bc::apply-bc ((bc bc-buoyancy) node mesh dt)
  "Arbitrary closure BC"
  (with-accessors ((datum bc-buoyancy-datum)
                   (rho bc-buoyancy-rho)
                   (clip-func bc-buoyancy-clip-func)
                   (sim bc-buoyancy-sim))
      bc
    (declare (function clip-func))

    (markup-cells-nodes sim bc)
    (let ((datum-rounding nil))
      (if datum-rounding
          (progn
            (let ((h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh sim))))
              (setf datum (* (round datum h) h)))
            (apply-buoyancy
             sim
             (lambda (pos)
               (buoyancy-virtual-stress (tref pos 1 0) datum rho (cl-mpm:sim-gravity sim)))
             (lambda (pos)
               (buoyancy-virtual-div (tref pos 1 0) datum rho (cl-mpm:sim-gravity sim)))
             (lambda (pos datum)
               (and
                (cell-clipping pos datum)
                (funcall clip-func pos datum)))
             datum)
            )
          (apply-buoyancy
           sim
           (lambda (pos)
             (buoyancy-virtual-stress (tref pos 1 0) datum rho (cl-mpm:sim-gravity sim)))
           (lambda (pos)
             (buoyancy-virtual-div (tref pos 1 0) datum rho (cl-mpm:sim-gravity sim)))
           (lambda (pos datum)
             (and
              (funcall clip-func pos datum)))
           datum)))
    (exchange-bc-data sim bc)
    ;;Reset pressure on MPs
    (with-accessors ((mesh cl-mpm:sim-mesh)
                     (mps cl-mpm:sim-mps))
        sim
      (cl-mpm:iterate-over-mps
       mps
       (lambda (mp)
         (with-accessors ((pos cl-mpm/particle::mp-position-trial)
                          (pressure cl-mpm/particle::mp-pressure)
                          (mp-datum cl-mpm/particle::mp-pressure-datum)
                          (mp-pfunc cl-mpm/particle::mp-pressure-func)
                          (mp-head cl-mpm/particle::mp-pressure-head)
                          (mp-boundary cl-mpm/particle::mp-boundary)
                          )
             mp
           (when ;(and (cell-clipping (cl-mpm/mesh::node-position node) datum)
                                        ;     (funcall clip-func (cl-mpm/mesh::node-position node) datum))

               (setf pressure 0d0)
             ;; (setf mp-pfunc
             ;;       (lambda (p)
             ;;         0d0))
             (setf mp-datum datum
                   mp-head rho
                   mp-boundary 0d0
                   )))))

      ;;Populate pressure on MPs
      (cl-mpm:iterate-over-mps
       mps
       (lambda (mp)
         (with-accessors ((pos cl-mpm/particle::mp-position-trial)
                          (pressure cl-mpm/particle::mp-pressure)
                          (mp-datum cl-mpm/particle::mp-pressure-datum)
                          (mp-head cl-mpm/particle::mp-pressure-head)
                          (mp-pfunc cl-mpm/particle::mp-pressure-func)
                          (mp-boundary cl-mpm/particle::mp-boundary)
                          (mp-volume cl-mpm/particle::mp-volume)
                          (mp-body-force cl-mpm/particle::mp-body-force))
             mp
           (cl-mpm::iterate-over-neighbours
            mesh mp
            (lambda (mesh mp node svp grads fsvp fgrad)
              (declare (double-float mp-boundary svp))
              (when t;(cl-mpm/mesh::node-boundary-node node)
                (when node
                  ;; (setf pressure (pressure-at-depth (tref pos 1 0) datum rho))
                  (when (and (cell-clipping (cl-mpm/mesh::node-position node) datum)
                             (funcall clip-func (cl-mpm/mesh::node-position node) datum))
                    (setf pressure (pressure-at-depth (tref pos 1 0) datum rho (cl-mpm:sim-gravity sim)))
                    ;; (setf mp-pfunc
                    ;;       (lambda (p)
                    ;;         (pressure-at-depth (tref p 1 0) datum rho)))
                    (setf mp-datum datum
                          mp-head rho)
                    (incf mp-boundary (* -1d0 svp (cl-mpm/mesh::node-boundary-scalar node)))
                    ;; (setf mp-boundary (cl-mpm/mesh:mesh-resolution mesh))
                    ;; (setf mp-boundary 1d3)
                    ))))))))

      (let ((damping (bc-viscous-damping bc)))
        (cl-mpm:iterate-over-nodes
         mesh
         (lambda (node)
           (when (cl-mpm/mesh:node-active node)
             (with-accessors ((force cl-mpm/mesh::node-damping-force)
                              (active cl-mpm/mesh:node-active)
                              (mass cl-mpm/mesh:node-mass)
                              (velocity cl-mpm/mesh:node-velocity)
                              (volume cl-mpm/mesh::node-volume)
                              (boundary cl-mpm/mesh::node-boundary-node)
                              (lock cl-mpm/mesh::node-lock)
                              (boundary-scalar cl-mpm/mesh::node-boundary-scalar))
                 node
               (sb-thread:with-mutex (lock)
                 (cl-mpm/fastmaths:fast-.-
                  force
                  (cl-mpm/fastmaths:fast-scale-vector
                   velocity
                   (*
                    1/2
                    damping
                    rho
                    (sqrt (max 0d0 (- boundary-scalar)))))
                  force))))))
        ;; (cl-mpm:iterate-over-mps
        ;;  mps
        ;;  (lambda (mp)
        ;;    (with-accessors ((pos cl-mpm/particle:mp-position)
        ;;                     (pressure cl-mpm/particle::mp-pressure)
        ;;                     (mp-boundary cl-mpm/particle::mp-boundary)
        ;;                     (mp-volume cl-mpm/particle::mp-volume)
        ;;                     (mp-velocity cl-mpm/particle::mp-velocity)
        ;;                     )
        ;;        mp
        ;;      (when (> mp-boundary 0d0)
        ;;        (cl-mpm::iterate-over-neighbours
        ;;         mesh mp
        ;;         (lambda (mesh mp node svp grads fsvp fgrad)
        ;;           (when (cl-mpm/mesh:node-active node)
        ;;             (with-accessors ((force cl-mpm/mesh::node-damping-force)
        ;;                              (active cl-mpm/mesh:node-active)
        ;;                              (mass cl-mpm/mesh:node-mass)
        ;;                              (velocity cl-mpm/mesh:node-velocity)
        ;;                              (volume cl-mpm/mesh::node-volume)
        ;;                              (boundary cl-mpm/mesh::node-boundary-node)
        ;;                              (lock cl-mpm/mesh::node-lock)
        ;;                              (boundary-scalar cl-mpm/mesh::node-boundary-scalar))
        ;;                 node
        ;;               (sb-thread:with-mutex (lock)
        ;;                 (cl-mpm/fastmaths:fast-.-
        ;;                  force
        ;;                  (cl-mpm/fastmaths:fast-scale-vector
        ;;                   ;; mp-velocity
        ;;                   velocity
        ;;                   (*
        ;;                    ;; (cl-mpm/fastmaths:mag vel)
        ;;                    1/2
        ;;                    damping
        ;;                    svp
        ;;                    ;; (sqrt mp-volume)
        ;;                    mp-volume
        ;;                    rho
        ;;                    (sqrt (max 0d0 (- boundary-scalar)))))
        ;;                  force)))
        ;;           )))))))
        ))))

(defun apply-viscous-damping ())

(defun set-pressure-all (sim bc)
  (with-accessors ((datum bc-buoyancy-datum)
                   (rho bc-buoyancy-rho)
                   (clip-func bc-buoyancy-clip-func)
                   ;(sim bc-buoyancy-sim)
                   )
      bc
    (with-accessors ((mesh cl-mpm:sim-mesh)
                     (mps cl-mpm:sim-mps))
        sim
      (loop for mp across mps
            do
               (with-accessors ((pos cl-mpm/particle:mp-position)
                                (pressure cl-mpm/particle::mp-pressure)
                                (mp-datum cl-mpm/particle::mp-pressure-datum)
                                (mp-head cl-mpm/particle::mp-pressure-head)
                                (mp-pfunc cl-mpm/particle::mp-pressure-func)
                                )
                   mp
                 (setf pressure (pressure-at-depth (tref pos 1 0) datum rho (cl-mpm:sim-gravity sim)))
                 ;; (setf pressure -1d0)
                 ;; (setf mp-pfunc
                 ;;       (lambda (p)
                 ;;         (pressure-at-depth (tref p 1 0) datum rho)))
                 (setf mp-datum datum
                       mp-head rho)
                 )))))

(defun get-cell-df (mesh point)
  (let ((dF (cl-mpm/utils:matrix-eye 1d0)))
    (cl-mpm::iterate-over-neighbours-point-linear
     mesh
     point
     (lambda (mesh node weight grads)
       (when (cl-mpm/mesh::node-active node)
         (cl-mpm/shape-function::@-combi-assemble-dstretch-3d grads (cl-mpm/mesh::node-displacment node) df))))
    dF))

(defun apply-force-cells-3d (mesh func-stress func-div clip-func)
  "Update force on nodes, with virtual stress field from cells"
  (declare (function func-stress func-div))
  (cl-mpm::iterate-over-cells
   mesh
   (lambda (cell)
     ;;Iterate over a cells nodes
     (with-accessors ((pos cl-mpm/mesh::cell-centroid)
                      (trial-pos cl-mpm/mesh::cell-trial-centroid)
                      (cell-active cl-mpm/mesh::cell-active)
                      (cell-buoyancy cl-mpm/mesh::cell-buoyancy)
                      (cell-pressure cl-mpm/mesh::cell-pressure)
                      (df cl-mpm/mesh::cell-deformation-gradient)
                      (disp cl-mpm/mesh::cell-displacement))
         cell
       (when t;cell-active
         (let* (
                ;; (pos (cl-mpm/fastmaths:fast-.+ pos disp))
                (cell-stress (funcall func-stress trial-pos))
                ;; (cell-stress (cl-mpm/utils::pull-back-voigt-stress cell-stress df))
                (cell-div (funcall func-div trial-pos))
                (f-stress (cl-mpm/utils:vector-zeros))
                (f-div (cl-mpm/utils:vector-zeros))
                ;; (df (cl-mpm/fastmaths:fast-.+ (cl-mpm/utils:matrix-eye 1d0)
                ;;                               ))
                )
           (setf cell-pressure (varef cell-stress 0))
           (cl-mpm/mesh::cell-iterate-over-neighbours
            mesh cell
            (lambda (mesh cell p volume node svp grads)
              (with-accessors ((node-force cl-mpm/mesh::node-force)
                               (node-force-ext cl-mpm/mesh::node-external-force)
                               (node-force-int cl-mpm/mesh::node-internal-force)
                               (node-pos cl-mpm/mesh::node-position)
                               (node-buoyancy-force cl-mpm/mesh::node-buoyancy-force)
                               (node-lock  cl-mpm/mesh:node-lock)
                               (node-active cl-mpm/mesh:node-active)
                               (node-boundary cl-mpm/mesh::node-boundary-node)
                               (node-volume cl-mpm/mesh::node-volume)
                               (node-boundary-scalar cl-mpm/mesh::node-boundary-scalar))
                  node
                (declare (double-float volume svp))
                (when (and node-active
                           node-boundary
                           ;; (funcall clip-func node-pos)
                           )
                  ;;Lock node
                  (cl-mpm/fastmaths:fast-zero f-stress)
                  (cl-mpm/forces::det-stress-force-unrolled cell-stress grads (- volume) f-stress)
                  (cl-mpm/fastmaths:fast-scale-vector
                   cell-div
                   (* volume svp)
                   f-div)
                  (let* ((f-total (cl-mpm/fastmaths::fast-.+ f-stress f-div)))
                    (sb-thread:with-mutex (node-lock)
                      (cl-mpm/fastmaths:fast-.- node-force-ext f-stress node-force-ext)
                      (cl-mpm/fastmaths:fast-.- node-force-ext f-div node-force-ext)
                      (cl-mpm/fastmaths:fast-.- node-buoyancy-force f-total node-buoyancy-force)
                      (incf node-boundary-scalar
                            (* -1d0 volume svp (the double-float (calculate-val-cell cell #'melt-rate)))))
                    )))))))))))

(defun apply-force-mps-3d (mesh mps func-stress func-div clip-func)
  "Update force on nodes, with virtual stress field from mps"
  (declare (function func-stress func-div clip-func))
  (cl-mpm:iterate-over-mps
   mps
    (lambda (mp)
      (compute-mp-displacement mesh mp)
      (with-accessors ((volume cl-mpm/particle:mp-volume)
                       ;(pos cl-mpm/particle::mp-position-trial)
                       (pos cl-mpm/particle::mp-position)
                       (disp cl-mpm/particle::mp-displacement-increment)
                       (df cl-mpm/particle::mp-deformation-gradient-increment)
                       (damage cl-mpm/particle::mp-damage))
          mp
        (when t;(funcall clip-func pos)
          (let* (
                 (mp-stress (funcall func-stress mp))
                 ;; (mp-stress (cl-mpm/utils::pull-back-voigt-stress mp-stress df))
                 (mp-div (funcall func-div mp))
                 (f-stress (cl-mpm/utils:vector-zeros))
                 (f-div (cl-mpm/utils:vector-zeros)))
            ;;Iterate over neighbour nodes
            (cl-mpm::iterate-over-neighbours
             mesh mp
             (lambda (mesh mp node svp grads fsvp fgrads)
               (when (cl-mpm/mesh:node-active node)
                 (with-accessors ((node-force cl-mpm/mesh::node-force)
                                  (node-force-ext cl-mpm/mesh::node-external-force)
                                  (node-force-int cl-mpm/mesh::node-internal-force)
                                  (node-pos cl-mpm/mesh::node-position)
                                  (node-buoyancy-force cl-mpm/mesh::node-buoyancy-force)
                                  (node-lock  cl-mpm/mesh:node-lock)
                                  (node-boundary cl-mpm/mesh::node-boundary-node)
                                  (node-boundary-scalar cl-mpm/mesh::node-boundary-scalar)
                                  (node-active  cl-mpm/mesh:node-active))
                     node
                   (declare (double-float volume svp))
                   (when (and node-boundary
                              ;; (funcall clip-func node-pos)
                              )
                     (cl-mpm/fastmaths:fast-zero f-stress)
                     (let (;(grads (cl-mpm::gradient-push-forwards grads df))
                           )
                       (cl-mpm/forces::det-stress-force-unrolled mp-stress grads (- volume) f-stress))
                     (cl-mpm/fastmaths:fast-scale-vector
                      mp-div
                      (* volume svp)
                      f-div)
                     (let* ((f-total (cl-mpm/fastmaths::fast-.+ f-stress f-div)))
                       (sb-thread:with-mutex (node-lock)
                         (cl-mpm/fastmaths:fast-.+ node-force-ext f-stress node-force-ext)
                         (cl-mpm/fastmaths:fast-.+ node-force-ext f-div node-force-ext)
                         (cl-mpm/fastmaths:fast-.+ node-buoyancy-force f-total node-buoyancy-force)
                         (incf node-boundary-scalar
                               (* volume svp (calculate-val-mp mp #'melt-rate))))))))))))))))

(defmethod cl-mpm/bc::apply-bc ((bc bc-scalar) node mesh dt)
  "Arbitrary closure BC"
  (with-accessors ((datum bc-buoyancy-datum)
                   (clip-func bc-buoyancy-clip-func)
                   (scalar-func bc-scalar-func)
                   (damage-volume bc-damage-volume)
                   (sim bc-buoyancy-sim))
      bc
    (declare (function clip-func))
    (apply-scalar
     sim
     scalar-func
     (lambda (pos)
       (and
        ;; (cell-clipping pos datum)
        ;; t
        (funcall clip-func pos)
        ))
     datum
     :damage-volume damage-volume
     )
    ))

(defun apply-scalar (sim func-scalar clip-function datum &key (damage-volume nil))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm::sim-mps))
      sim
    (with-accessors ((h cl-mpm/mesh:mesh-resolution))
        mesh
      (locate-mps-cells sim clip-function)
      ;; (populate-cells-volume mesh clip-function)
      ;; (locate-mps-cells sim clip-function)
      ;; (populate-nodes-volume mesh clip-function)
      ;; (populate-nodes-volume-damage mesh clip-function)
      ;; (populate-nodes-domain mesh clip-function)

      (apply-scalar-mps-3d mesh mps
                          (lambda (mp)
                            (calculate-val-mp mp func-scalar))
                          clip-function
                          :damage-volume damage-volume
                          )
      (apply-scalar-cells-3d mesh
                            func-scalar
                            clip-function)
      )))



(defun apply-scalar-cells-3d (mesh func-scalar clip-func)
  "Update force on nodes, with virtual stress field from cells"
  (declare (function func-scalar))
  (cl-mpm::iterate-over-cells
   mesh
   (lambda (cell)
     (when t;(loop for n in (cl-mpm/mesh::cell-nodes cell) thereis (cl-mpm/mesh::node-boundary-node n))
       (let ((nodal-volume 0d0))
         ;;Possibly clip ill posed cells
         (when t;(> (/ nodal-volume (cl-mpm/mesh::cell-volume cell)) 1d-5)
           ;;Iterate over a cells nodes
           (let ((dsvp (cl-mpm/utils::dsvp-3d-zeros)))
             (cl-mpm/mesh::cell-iterate-over-neighbours
              mesh cell
              (lambda (mesh cell pos volume node svp grads)
                (with-accessors ((node-force cl-mpm/mesh::node-force)
                                 (node-pos cl-mpm/mesh::node-position)
                                 (node-buoyancy-force cl-mpm/mesh::node-buoyancy-force)
                                 (node-lock  cl-mpm/mesh:node-lock)
                                 (node-active cl-mpm/mesh:node-active)
                                 (node-boundary cl-mpm/mesh::node-boundary-node)
                                 (node-volume cl-mpm/mesh::node-volume-true)
                                 (node-boundary-scalar cl-mpm/mesh::node-boundary-scalar))
                    node
                  (declare (double-float volume svp))
                  (when (and node-active
                             node-boundary
                             (funcall clip-func node-pos)
                             )
                    ;;Lock node
                    (sb-thread:with-mutex (node-lock)
                      (incf node-boundary-scalar
                            (* -1d0 volume svp (the double-float (calculate-val-cell cell func-scalar))))
                      ))))))))))))

(defun apply-scalar-mps-3d (mesh mps func-scalar clip-func &key (damage-volume nil))
  "Update force on nodes, with virtual stress field from mps"
  (declare (function func-scalar))
  (cl-mpm::iterate-over-mps
   mps
   (lambda (mp)
     (when t;(< (cl-mpm/particle::mp-damage mp) 1d0)
       (with-accessors ((volume cl-mpm/particle:mp-volume)
                        (pos cl-mpm/particle::mp-position)
                        (damage cl-mpm/particle::mp-damage)
                        )
           mp
         (let ((dsvp (cl-mpm/utils::dsvp-3d-zeros)))
           ;;Iterate over neighbour nodes
           (cl-mpm::iterate-over-neighbours
            mesh mp
            (lambda (mesh mp node svp grads fsvp fgrads)
              (with-accessors ((node-buoyancy-force cl-mpm/mesh::node-buoyancy-force)
                               (node-lock  cl-mpm/mesh:node-lock)
                               (node-pos  cl-mpm/mesh::node-position)
                               (node-boundary cl-mpm/mesh::node-boundary-node)
                               (node-boundary-scalar cl-mpm/mesh::node-boundary-scalar)
                               (node-active  cl-mpm/mesh:node-active))
                  node
                (declare (double-float volume svp))
                (when (and node-active
                           node-boundary
                           (funcall clip-func node-pos)
                           )
                  (sb-thread:with-mutex (node-lock)
                    (incf node-boundary-scalar (* (if damage-volume (- 1d0 damage) 1d0)
                                                  volume svp (funcall func-scalar mp))))))))))))))

(defclass bc-buoyancy-body (bc-buoyancy)
  ()
  (:documentation "A nonconforming buoyancy bc"))

(defun make-bc-buoyancy-body (sim datum rho clip-func &key (datum-rounding nil))
  (declare (function clip-func))
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (with-accessors ((h cl-mpm/mesh:mesh-resolution))
        mesh
      (when datum-rounding
        (let ((h (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh sim))))
          (setf datum (* (round datum h) h))))
      (make-instance 'bc-buoyancy-body
                     :index nil
                     :sim sim
                     :clip-func (lambda (pos)
                                  (and (funcall clip-func pos datum)
                                       (< (cl-mpm/utils::varef pos 1) datum)))
                     :rho rho
                     :datum datum))))

(defmethod cl-mpm/bc::apply-bc ((bc bc-buoyancy-body) node mesh dt)
  "Arbitrary closure BC"
  (with-accessors ((datum bc-buoyancy-datum)
                   (rho bc-buoyancy-rho)
                   (clip-func bc-buoyancy-clip-func)
                   (sim bc-buoyancy-sim))
      bc
    (declare (function clip-func))
    (with-accessors ((mesh cl-mpm:sim-mesh)
                     (mps cl-mpm:sim-mps))
        sim
      ;;Apply body force
      (markup-cells-nodes sim bc)
      (apply-buoyancy-body
       mesh mps
       (lambda (mp) (calculate-val-mp mp (lambda (pos) (buoyancy-virtual-div (tref pos 1 0) datum rho (cl-mpm:sim-gravity sim)))))
       ;; (lambda (mp)
       ;;   (cl-mpm/fastmaths:fast-scale-vector
       ;;    (cl-mpm/particle::mp-gravity-axis mp)
       ;;    (* -1d0 rho (cl-mpm/particle:mp-gravity mp))))
       clip-func)
      ;;Set all hydrostatic pressures
      (cl-mpm:iterate-over-mps
       mps
       (lambda (mp)
         (with-accessors ((pos cl-mpm/particle:mp-position)
                          (pressure cl-mpm/particle::mp-pressure)
                          (mp-datum cl-mpm/particle::mp-pressure-datum)
                          (mp-pfunc cl-mpm/particle::mp-pressure-func)
                          (mp-head cl-mpm/particle::mp-pressure-head)
                          )
             mp
           (when t;(and (cell-clipping (cl-mpm/mesh::node-position node) datum)
             (setf pressure 0d0)
             (setf mp-datum datum
                   mp-head rho)))))
      (cl-mpm:iterate-over-mps
       mps
       (lambda (mp)
         (with-accessors ((pos cl-mpm/particle:mp-position)
                          (pressure cl-mpm/particle::mp-pressure)
                          (mp-datum cl-mpm/particle::mp-pressure-datum)
                          (mp-head cl-mpm/particle::mp-pressure-head)
                          (mp-pfunc cl-mpm/particle::mp-pressure-func)
                          )
             mp
           (cl-mpm::iterate-over-neighbours
            mesh mp
            (lambda (mesh mp node svp grads fsvp fgrad)
              (when t;(cl-mpm/mesh::node-boundary-node node)
                (when node
                  (when (and (cell-clipping (cl-mpm/mesh::node-position node) datum)
                             (funcall clip-func (cl-mpm/mesh::node-position node) datum))
                    (setf pressure (pressure-at-depth (tref pos 1 0) datum rho (cl-mpm:sim-gravity sim)))
                    (setf mp-datum datum
                          mp-head rho))))))))))))

(defun apply-buoyancy-body (mesh mps func-div clip-func)
  "Update force on nodes, with virtual stress field from mps"
  (declare (function func-div))
  (cl-mpm:iterate-over-mps
   mps
    (lambda (mp)
      (when t
        (with-accessors ((volume cl-mpm/particle:mp-volume)
                         (pos cl-mpm/particle::mp-position)
                         (damage cl-mpm/particle::mp-damage)
                         )
            mp
          (let ((dsvp (cl-mpm/utils::dsvp-3d-zeros)))
            ;;Iterate over neighbour nodes
            (cl-mpm::iterate-over-neighbours
             mesh mp
             (lambda (mesh mp node svp grads fsvp fgrads)
               (with-accessors ((node-force cl-mpm/mesh::node-force)
                                (node-force-ext cl-mpm/mesh::node-external-force)
                                (node-pos cl-mpm/mesh::node-position)
                                (node-buoyancy-force cl-mpm/mesh::node-buoyancy-force)
                                (node-lock  cl-mpm/mesh:node-lock)
                                (node-boundary cl-mpm/mesh::node-boundary-node)
                                (node-boundary-scalar cl-mpm/mesh::node-boundary-scalar)
                                (node-active  cl-mpm/mesh:node-active))
                   node
                 (declare (double-float volume svp))
                 (when (and node-active
                            (funcall clip-func node-pos)
                            ;; (funcall clip-func pos)
                            )
                   ;;Lock node for multithreading
                   (sb-thread:with-mutex (node-lock)
                     (setf node-boundary t)
                     ;; Add gradient of stress
                     ;; Add divergance of stress
                     (let ((buoyancy-force (cl-mpm/fastmaths::fast-scale!
                                            (funcall func-div mp)
                                            (* svp volume))))
                       ;; (cl-mpm/fastmaths::fast-.+ node-force
                       ;;                           buoyancy-force
                       ;;                           node-force)
                       (cl-mpm/fastmaths::fast-.+ node-force-ext
                                                  buoyancy-force
                                                  node-force-ext)
                       (cl-mpm/fastmaths::fast-.+ node-buoyancy-force
                                                  buoyancy-force
                                                  node-buoyancy-force)
                       )
                     ;;Add boundary scalar (can be used for other ops)
                     (incf node-boundary-scalar
                           (* volume svp (calculate-val-mp mp #'melt-rate))))))))))))))

(defun apply-buoyancy-body-damage (mesh mps func-div clip-func datum)
  "Update force on nodes, with virtual stress field from mps"
  (declare (function func-div))
  (cl-mpm:iterate-over-mps
   mps
    (lambda (mp)
      (when t
        (with-accessors ((volume cl-mpm/particle:mp-volume)
                         (pos cl-mpm/particle::mp-position)
                         (damage cl-mpm/particle::mp-damage)
                         )
            mp
          (when (< (magicl:tref pos 1 0) datum)
            (let ((bf (cl-mpm/fastmaths::fast-scale!
                       (funcall func-div mp)
                       (* damage volume))))
              ;;Iterate over neighbour nodes
              (cl-mpm::iterate-over-neighbours
               mesh mp
               (lambda (mesh mp node svp grads fsvp fgrads)
                 (with-accessors ((node-force cl-mpm/mesh::node-force)
                                  (node-force-ext cl-mpm/mesh::node-external-force)
                                  (node-pos cl-mpm/mesh::node-position)
                                  (node-buoyancy-force cl-mpm/mesh::node-buoyancy-force)
                                  (node-lock  cl-mpm/mesh:node-lock)
                                  (node-boundary cl-mpm/mesh::node-boundary-node)
                                  (node-boundary-scalar cl-mpm/mesh::node-boundary-scalar)
                                  (node-active  cl-mpm/mesh:node-active))
                     node
                   (declare (double-float volume svp))
                   (when (and node-active
                              (funcall clip-func pos))
                     ;;Lock node for multithreading
                     (sb-thread:with-mutex (node-lock)
                       (let ((buoyancy-force (cl-mpm/fastmaths:fast-scale-vector bf svp)))
                         ;; (cl-mpm/fastmaths::fast-.+ node-force
                         ;;                           buoyancy-force
                         ;;                           node-force)
                         (cl-mpm/fastmaths::fast-.+ node-force-ext
                                                   buoyancy-force
                                                   node-force-ext)
                         (cl-mpm/fastmaths::fast-.+ node-buoyancy-force
                                                   buoyancy-force
                                                   node-buoyancy-force)
                         )
                       ;;Add boundary scalar (can be used for other ops)
                       (incf node-boundary-scalar
                             (* volume svp (calculate-val-mp mp #'melt-rate)))))))))))))))


(defmethod cl-mpm/bc::assemble-bc-stiffness (sim (bc bc-buoyancy))
  (with-accessors ((datum bc-buoyancy-datum)
                   (clip-func bc-buoyancy-clip-func)
                   (rho bc-buoyancy-rho))
      bc
    (declare (function clip-func))
    ;; (locate-mps-cells sim (lambda (pos) (funcall clip-func pos datum)))
    (markup-cells-nodes sim bc)
    (compute-stiffness-cells-3d
     (cl-mpm:sim-mesh sim)
     (lambda (pos) (abs (* (max 0d0 (- datum (varef pos 1))) rho (cl-mpm:sim-gravity sim))))
     clip-func))
  ;; (apply-force-mps-3d mesh mps
  ;;                     (lambda (mp) (calculate-val-mp mp func-stress))
  ;;                     (lambda (mp) (calculate-val-mp mp func-div))
  ;;                     (lambda (pos) (funcall clip-function pos datum))
  ;;                     )

  )

(defun compute-stiffness-cells-3d (mesh stiffness-func clip-func)
  "Update force on nodes, with virtual stress field from cells"
  (declare (function stiffness-func))
  (cl-mpm::iterate-over-cells
   mesh
   (lambda (cell)
     (with-accessors ((pos cl-mpm/mesh::cell-centroid)
                      (trial-pos cl-mpm/mesh::cell-trial-centroid)
                      (cell-active cl-mpm/mesh::cell-active)
                      (cell-buoyancy cl-mpm/mesh::cell-buoyancy)
                      (cell-pressure cl-mpm/mesh::cell-pressure)
                      (df cl-mpm/mesh::cell-deformation-gradient)
                      (disp cl-mpm/mesh::cell-displacement))
         cell
       (when cell-active
         (let* ()
           (cl-mpm/mesh::cell-iterate-over-neighbours
            mesh cell
            (lambda (mesh cell p volume node svp grads)
              (with-accessors ((node-pos cl-mpm/mesh::node-position)
                               (node-mass cl-mpm/mesh::node-mass)
                               (node-lock  cl-mpm/mesh:node-lock)
                               (node-active cl-mpm/mesh:node-active)
                               (node-boundary cl-mpm/mesh::node-boundary-node)
                               (node-volume cl-mpm/mesh::node-volume)
                               )
                  node
                (declare (double-float volume svp))
                (when (and node-active
                           node-boundary)
                  (sb-thread:with-mutex (node-lock)
                    (incf
                     node-mass
                     (max 0d0 (*
                               0d0
                               (funcall
                                stiffness-func
                                (cl-mpm/fastmaths::fast-.+ node-pos (cl-mpm/mesh::node-displacment node))) volume svp)))
                    )))))))))))
