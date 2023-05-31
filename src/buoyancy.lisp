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
(defun pressure-at-depth (z datum-true rho)
  (let* ((g -9.8)
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
  (let* ((g -9.8d0)
         (datum datum-true)
         (h (- datum z))
         (f (* rho g h))
         )
    (if (> h 0d0)
        (magicl:from-list (list f f 0d0) '(3 1) :type 'double-float)
        (magicl:zeros '(3 1))
        )
    ))
(defun buoyancy-virtual-div (z datum-true rho)
  (let* ((g -9.8d0)
         (datum datum-true)
         (h (- datum z))
         (f (* -1d0 rho g))
         )
    (if (> h 0d0)
        (magicl:from-list (list 0 f) '(2 1) :type 'double-float)
        (magicl:zeros '(2 1))
        )
    ))
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
        (let ((dsvp (cl-mpm/utils::stretch-dsvp-zeros)))
          ;;Iterate over neighbour nodes
          (cl-mpm::iterate-over-neighbours
           mesh mp
           (lambda (mesh mp node svp grads fsvp fgrads)
             (with-accessors ((node-force cl-mpm/mesh:node-force)
                              (node-lock  cl-mpm/mesh:node-lock)
                              (node-boundary cl-mpm/mesh::node-boundary-node)
                              (node-boundary-scalar cl-mpm/mesh::node-boundary-scalar)
                              (node-active  cl-mpm/mesh:node-active))
                 node
               (declare (double-float volume svp))
               (when (and node-active
                          node-boundary
                          )
                 ;;Lock node for multithreading
                 (sb-thread:with-mutex (node-lock)

                   (cl-mpm/shape-function::assemble-dsvp-2d-prealloc grads dsvp)
                   (cl-mpm/fastmath::mult-transpose-accumulate dsvp
                                                               (funcall func-stress mp)
                                                               (* volume)
                                                               node-force);
                   (cl-mpm/fastmath::fast-fmacc node-force
                                                (funcall func-div mp)
                                                (* svp volume))
                   ;;Add the gradient of stress
                   ;; (setf node-force (magicl:.+
                   ;;                   node-force
                   ;;                   (magicl:scale!
                   ;;                    (magicl:@
                   ;;                     (magicl:transpose
                   ;;                      (cl-mpm/shape-function::assemble-dsvp-2d grads))
                   ;;                     (funcall func-stress mp)
                   ;;                     )
                   ;;                    (* volume))))
                   ;; ;;Add the divergance
                   ;; (setf node-force (magicl:.+
                   ;;                   node-force
                   ;;                   (magicl:scale!
                   ;;                    (funcall func-div mp)
                   ;;                    (* svp volume))))
                   (incf node-boundary-scalar
                         (* volume svp (calculate-val-mp mp #'melt-rate)))
                   ))))))))))

(defun melt-rate (pos)
  (if (< (magicl:tref pos 1 0) 300)
      ;(exp (* 0.1 (- (magicl:tref pos 1 0) 300)))
      (+ 1 (exp (* 0.05 (- (magicl:tref pos 1 0) 300))))
      ;; (/ 1 (+ 1 (expt
      ;;            (abs (min 0 (- (magicl:tref pos 1 0) 300))) 2)))
      0d0))

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
                   (let ((dsvp (cl-mpm/utils::stretch-dsvp-zeros)))
                     (cl-mpm/mesh::cell-iterate-over-neighbours
                      mesh cell
                      (lambda (mesh cell pos volume node svp grads)
                        (with-accessors ((node-force cl-mpm/mesh:node-force)
                                         (node-lock  cl-mpm/mesh:node-lock)
                                         (node-active cl-mpm/mesh:node-active)
                                         (node-boundary cl-mpm/mesh::node-boundary-node)
                                         (node-volume cl-mpm/mesh::node-volume)
                                         (node-boundary-scalar cl-mpm/mesh::node-boundary-scalar)
                                         )
                            node
                          (declare (double-float volume svp))
                          (when (and node-active
                                     node-boundary
                                     )
                            ;;Lock node
                            (sb-thread:with-mutex (node-lock)
                              ;;Subtract gradient of stress from node force
                              (cl-mpm/shape-function::assemble-dsvp-2d-prealloc grads dsvp)
                              (cl-mpm/fastmath::mult-transpose-accumulate dsvp
                                                                          (funcall func-stress pos)
                                                                          (* -1d0 volume)
                                                                          node-force);
                              (cl-mpm/fastmath::fast-fmacc node-force
                                                           (funcall func-div pos)
                                                           (* -1d0 svp volume))
                              ;; (setf node-force (magicl:.+
                              ;;                   node-force
                              ;;                   (magicl:scale!
                              ;;                    (magicl:@
                              ;;                     (magicl:transpose
                              ;;                      (cl-mpm/shape-function::assemble-dsvp-2d grads))
                              ;;                     (funcall func-stress pos))
                              ;;                    (* -1d0 volume))))
                              ;; ;;Subtract stress divergance from node force
                              ;; (setf node-force (magicl:.+
                              ;;                   node-force
                              ;;                   (magicl:scale!
                              ;;                    (funcall func-div pos)
                              ;;                    (* -1d0 svp volume))))
                              (incf node-boundary-scalar
                                    (* -1d0 volume svp (the double-float (calculate-val-cell cell #'melt-rate))))
                              )))
                        )))))))))
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
                         (and (cl-mpm/mesh:node-active n) (> (cl-mpm/mesh:node-mass n) 0d0))
                         ) nodes)
            (setf mp-count 1))
          )))))

(defun cell-clipping (pos)
  (<= (magicl:tref pos 1 0) 300)
  ;; (<= (magicl:tref pos 1 0) 300)
  ;; t
  )
(defun populate-nodes-volume (mesh)
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
         (when t;(cell-clipping pos)
           (when active
             (setf boundary t))))))))
(defun populate-nodes-domain (mesh)
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
       (setf boundary t)))))

(defun locate-mps-cells (mesh mps)
  "Mark boundary nodes based on neighbour MP inclusion"
  (let ((cells (cl-mpm/mesh::mesh-cells mesh)))
    ;;The aproach described by the paper
    ;; (populate-cell-mp-count mesh mps)
    (populate-cell-mp-count-gimp mesh mps)
    ;; (populate-cell-mp-count-volume mesh mps)
    ;; (populate-cell-nodes mesh mps)
    ;; (prune-buoyancy-nodes mesh '(0 0) 300)
    (lparallel:pdotimes (i (array-total-size cells))
      (let ((cell (row-major-aref cells i)))
        (with-accessors ((mp-count cl-mpm/mesh::cell-mp-count)
                         (neighbours cl-mpm/mesh::cell-neighbours)
                         (index cl-mpm/mesh::cell-index)
                         (nodes cl-mpm/mesh::cell-nodes)
                         (pruned cl-mpm/mesh::cell-pruned)
                         (boundary cl-mpm/mesh::cell-boundary)
                         (pos cl-mpm/mesh::cell-centroid)
                         )
            cell
          (setf boundary nil)
          (when (and (= mp-count 0)
                     (cell-clipping pos) (not pruned))
              (loop for neighbour in neighbours
                    do
                       (when (cell-clipping (cl-mpm/mesh::cell-centroid neighbour))
                           (check-neighbour-cell neighbour)
                         (setf boundary t)
                         (loop for n in nodes
                               do
                                  (when (cl-mpm/mesh:node-active n)
                                    (sb-thread:with-mutex ((cl-mpm/mesh:node-lock n))
                                      (setf (cl-mpm/mesh::node-boundary-node n) t))))))
              ))))
    )
  )

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
(defun apply-non-conforming-nuemann (sim func-stress func-div)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm::sim-mps))
      sim
    (with-accessors ((h cl-mpm/mesh:mesh-resolution))
        mesh
      ;; (locate-mps-cells mesh mps)
      (populate-nodes-volume mesh)
      ;; (populate-nodes-domain mesh)
      (apply-force-mps mesh mps
                       (lambda (mp) (calculate-val-mp mp func-stress))
                       (lambda (mp) (calculate-val-mp mp func-div)))
      (apply-force-cells mesh
                         func-stress
                         func-div
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
    :initarg :pressures))
  (:documentation "A nonconforming pressure bc"))


(defmethod cl-mpm/bc::apply-bc ((bc bc-pressure) node mesh dt)
  "Arbitrary closure BC"
  (with-accessors ((pressures bc-pressure-pressures)
                   (sim bc-pressure-sim))
      bc
    (apply-non-conforming-nuemann
     sim
     (lambda (pos)
       (pressure-virtual-stress (first pressures) (second pressures)))
     (lambda (pos)
       (pressure-virtual-div))
     )))

(defclass bc-nc-buoyancy (bc)
  ()
  (:documentation "A nonconforming buoyancy bc"))

(defun make-bc-pressure (sim pressure-x pressure-y)
  (make-instance 'bc-pressure
                 :index '(0 0)
                 :sim sim
                 :pressures (list pressure-x pressure-y)
                 )
  ;; (make-instance 'cl-mpm/bc::bc-closure
  ;;                :index '(0 0)
  ;;                :func (lambda ()
  ;;                        (apply-non-conforming-nuemann
  ;;                         sim
  ;;                         (lambda (pos)
  ;;                           (pressure-virtual-stress pressure-x pressure-y))
  ;;                         (lambda (pos)
  ;;                           (pressure-virtual-div))
  ;;                         )
  ;;                        ))
  )
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
                                       (lambda (mesh mp node svp grads fsvp fgrad)
                                         (when t;(cl-mpm/mesh::node-boundary-node node)
                                           (with-accessors ((pos cl-mpm/particle:mp-position)
                                                            (pressure cl-mpm/particle::mp-pressure))
                                               mp
                                             (setf pressure (pressure-at-depth (tref pos 1 0) datum rho)))
                                           )))))))))
