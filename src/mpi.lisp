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
;; (declaim (optimize (debug 3) (safety 3) (speed 0)))
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
    :initform nil
    ;; '((0 0)
    ;;             (0 0)
    ;;             (0 0))
    )
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

(defmethod sim-add-mp ((sim mpm-sim-mpi) mp)
  (with-accessors ((uid-counter cl-mpm::sim-unique-index-counter)
                   (lock cl-mpm::sim-unique-index-lock))
      sim
    (let* ((size (cl-mpi:mpi-comm-size))
           (rank (cl-mpi:mpi-comm-rank))
           (shift (ceiling (log size 2))))
      (sb-thread:with-mutex (lock)
        (setf (cl-mpm/particle::mp-unique-index mp) (+ (ash uid-counter shift) rank))
        (incf uid-counter)))
    (call-next-method)))


(defgeneric update-min-domain-size (sim))
(defmethod update-min-domain-size ((sim mpm-sim-mpi))
  (setf (mpm-sim-mpi-min-size sim)
        (* (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) (mpm-sim-mpi-halo-depth sim))))
(defmethod update-min-domain-size ((sim mpm-sim-mpi-nodes-damage))
  (setf (mpm-sim-mpi-min-size sim)
        (max
         (* (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) (mpm-sim-mpi-halo-depth sim))
         (mpm-sim-mpi-halo-damage-size sim))))
(defmethod (setf cl-mpm::sim-mesh) :after (value (sim mpm-sim-mpi))
  (update-min-domain-size sim))
(defmethod (setf mpm-sim-mpi-halo-depth) :after (value (sim mpm-sim-mpi))
  (update-min-domain-size sim))
(defmethod (setf mpm-sim-mpi-halo-damage-size) :after (value (sim mpm-sim-mpi-nodes-damage))
  (update-min-domain-size sim))

(defclass mpm-sim-mpi-nodes (mpm-sim-mpi)
  ((halo-node-list
    :accessor mpm-sim-mpi-halo-node-list
    :initform (loop for i from 0 to 2
                    collect (loop for i from 0 to 1 collect (make-array 0 :element-type t))))
   )
  )

(defmacro push-bytes (array bytes index)
  `(progn
    ;(declare ((simple-array (unsigned-byte 8) *) array bytes))
     (loop for b across (the (simple-array (unsigned-byte 8) *) ,bytes)
          do
             (progn
               (setf (aref (the (simple-array (unsigned-byte 8) *) ,array) ,index) b)
               (incf ,index)))))

(defmacro pull-bytes (array index length)
  `(progn
                                        ;(declare ((simple-array (unsigned-byte 8) *) array bytes))
     (let ((o (make-array ,length :element-type '(unsigned-byte 8) :displaced-to ,array :displaced-index-offset ,index)))
       (declare (dynamic-extent o))
       (incf ,index ,length)
       o
       )))


(defmacro push-int (output value inc)
  (let ((bytes-per-int 2))
    `(push-bytes ,output (cl-intbytes:int->octets ,value ,bytes-per-int) ,inc)))

(defmacro push-index (output value inc)
  `(progn
     (push-int ,output (nth 0 ,value) ,inc)
     (push-int ,output (nth 1 ,value) ,inc)
     (push-int ,output (nth 2 ,value) ,inc)
     ))

(defmacro push-float (output value inc)
  `(push-bytes ,output (cl-intbytes:int64->octets (ieee-floats:encode-float64 ,value)) ,inc))
(defmacro push-vector (output vec inc)
  `(let ((v-s (magicl::matrix/double-float-storage ,vec)))
     (push-float ,output (aref v-s 0) ,inc)
     (push-float ,output (aref v-s 1) ,inc)
     (push-float ,output (aref v-s 2) ,inc)))


(eval-when
    (:compile-toplevel
     :load-toplevel
     :execute)
  (defun serialise-length (type)
    (ecase type
      (int 2)
      (index (* 2 3))
      (float 8)
      (vector (* 8 3))
      (t 0))))

(defmacro pull-int (array inc)
  (let ((bytes-per-int 2))
    `(cl-intbytes:octets->int (pull-bytes ,array ,inc ,bytes-per-int) ,bytes-per-int)))

(defmacro pull-float (array inc)
  (let ((bytes-per-float 8))
    `(ieee-floats:decode-float64 (cl-intbytes:octets->uint64 (pull-bytes ,array ,inc ,bytes-per-float)))))

(defmacro pull-vector (array inc)
  `(cl-mpm/utils:vector-from-list
    (list (pull-float ,array ,inc)
          (pull-float ,array ,inc)
          (pull-float ,array ,inc))))

(defmacro pull-index (array inc)
  `(list (pull-int ,array ,inc)
         (pull-int ,array ,inc)
         (pull-int ,array ,inc)))


(defmacro make-mpi-ser (name mapping-list)
  "Creates an MPI structure, with a list of variables and a mapping set to a CLOS object "
  (let ((ser-name (intern (format nil "SERIALISE-~:@(~A~)" name)))
        (deser-name (intern (format nil "DESERIALISE-~:@(~A~)" name)))
        (mpi-object-name (intern (format nil "MPI-OBJECT-~:@(~A~)" name)))
        (mpi-object-constructor (intern (format nil "MAKE-MPI-OBJECT-~:@(~A~)" name)))
        (packet-size (loop for map in mapping-list
                           sum (serialise-length (first map)))))
    ;; (format t "~A~%" packet-size)
    `(progn
       ;; (declaim (inline ,mpi-object-constructor))
       (defstruct ,mpi-object-name
         ,@(mapcar (lambda (slot-entry)
                     (destructuring-bind (type var-name accessor-name)
                         slot-entry
                       var-name)
                     ) mapping-list)
         )
       (defun ,ser-name (objects)
         (let* ((node-count (length objects))
                (packet-size ,packet-size)
                (output (static-vectors:make-static-vector (* node-count packet-size) :element-type '(unsigned-byte 8))))
           (declare ((simple-array (unsigned-byte 8) *) output))
           (lparallel:pdotimes (i node-count)
             (let* ((inc (* i packet-size))
                    (obj (aref objects i)))
               (with-accessors ,(mapcar (lambda (slot-entry)
                                          (with-current-source-form (slot-entry mapping-list)
                                            (unless (sb-int::proper-list-of-length-p slot-entry 3)
                                              (error "Malformed slot entry: ~s, should ~
                                  be (type variable-name accessor-name)"
                                                     slot-entry))
                                            (destructuring-bind (type var-name accessor-name)
                                                slot-entry
                                              `(,var-name ,accessor-name))
                                            ))
                                 mapping-list)
                   obj
                 ,@(mapcar (lambda (slot-entry)
                             (destructuring-bind (type var-name accessor-name)
                                 slot-entry
                               (let ((push-name (intern (format nil "PUSH-~:@(~A~)" type))))
                                 `(,push-name output ,var-name inc)
                                 ))
                             ) mapping-list)
                 )
               ))
           output
           )
         )
       (defun ,deser-name (array)
         (let* ((node-count (floor (length array) ,packet-size))
                (output (make-array node-count :element-type ',mpi-object-name))
                (inc 0))
           (declare (fixnum inc)
                    ((simple-array (unsigned-byte 8) *) array))
           (lparallel:pdotimes (i node-count)
             (let ((inc (* i ,packet-size)))
               (setf (aref output i)
                     (,mpi-object-constructor
                      ,@(apply #'append
                               (mapcar (lambda (slot-entry)
                                         (destructuring-bind (type var-name accessor-name)
                                             slot-entry
                                           (let ((pull-name (intern (format nil "PULL-~:@(~A~)" type)))
                                                 (var-keyword (intern (format nil "~:@(~A~)" var-name) "KEYWORD")))
                                             `(,var-keyword (,pull-name array inc))))) mapping-list))
                      ))
               ))
           output))
         )
       ))

(make-mpi-ser
 damage-mp
 ((vector position cl-mpm/particle::mp-position)
  (float volume cl-mpm/particle::mp-volume)
  (float damage cl-mpm/particle::mp-damage)
  (float y cl-mpm/particle::mp-damage-y-local)
  (float local-length cl-mpm/particle::mp-true-local-length)
  ))

(make-mpi-ser
 node
 (
  (index index cl-mpm/mesh::node-index)
  (float mass cl-mpm/mesh::node-mass)
  (float pmod cl-mpm/mesh::node-pwave)
  (float svp cl-mpm/mesh::node-svp-sum)
  (float vol cl-mpm/mesh::node-volume)
  (float j-inc cl-mpm/mesh::node-jacobian-inc)
  (vector velocity cl-mpm/mesh::node-velocity)
  (vector force cl-mpm/mesh::node-force)
  (vector force-int cl-mpm/mesh::node-internal-force)
  (vector force-ext cl-mpm/mesh::node-external-force)
  (vector force-damping cl-mpm/mesh::node-damping-force)
  (vector force-buoyancy cl-mpm/mesh::node-buoyancy-force)
  ))

(defun exchange-node-like (sim
                           serialise
                           deserialise
                           func)
  (declare (function func))
  (let* ((rank (cl-mpi:mpi-comm-rank))
         (size (cl-mpi:mpi-comm-size)))
    (with-accessors ((mesh cl-mpm:sim-mesh))
        sim
      (let* ((nd-nodes (cl-mpm/mesh:mesh-nodes mesh))
             (index (mpi-rank-to-index sim rank))
             (bounds-list (mpm-sim-mpi-domain-bounds sim))
             (h (cl-mpm/mesh:mesh-resolution mesh))
             (halo-depth (mpm-sim-mpi-halo-depth sim))
             (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
             )
        (loop for i from 0 below nd
              do
                 (let ((id-delta (list 0 0 0)))
                   (setf (nth i id-delta) 1)
                   (let ((left-neighbor (mpi-index-to-rank sim (mapcar (lambda (a b) (declare (fixnum a b)) (- a b)) index id-delta)))
                         (right-neighbor (mpi-index-to-rank sim (mapcar (lambda (a b) (declare (fixnum a b)) (+ a b)) index id-delta))))
                     (declare (fixnum left-neighbor right-neighbor))
                     (destructuring-bind (bl bu) (nth i bounds-list)
                       (when (not
                              (and
                               (= left-neighbor -1)
                               (= right-neighbor -1)))
                         (labels
                             ((active-filter (nodes)
                                (let ((res
                                        (lparallel:premove-if-not
                                         (lambda (mp)
                                           (and
                                            (cl-mpm/mesh:node-active mp)))
                                         nodes)))
                                  (make-array (length res) :initial-contents res))))
                           (let ((left-filter (nth 0 (nth i (mpm-sim-mpi-halo-node-list sim))))
                                 (right-filter (nth 1 (nth i (mpm-sim-mpi-halo-node-list sim))))
                                 )
                             (declare (fixnum left-neighbor right-neighbor))
                             ;; (format t "Rank ~D - left ~A~%" rank (length left-filter))
                             ;; (format t "Rank ~D - righ ~A~%" rank (length right-filter))
                             (let* ((cl-mpi-extensions::*standard-encode-function* serialise)
                                    (cl-mpi-extensions::*standard-decode-function* deserialise)
                                    (recv
                                      (cond
                                        ((and (not (= left-neighbor -1))
                                              (not (= right-neighbor -1)))
                                         (let ((l (active-filter left-filter))
                                               (r (active-filter right-filter)))
                                           (cl-mpi-extensions:mpi-waitall-anything
                                            (cl-mpi-extensions:mpi-irecv-anything right-neighbor :tag 1)
                                            (cl-mpi-extensions:mpi-irecv-anything left-neighbor :tag 2)
                                            (cl-mpi-extensions:mpi-isend-anything
                                             l
                                             left-neighbor :tag 1)
                                            (cl-mpi-extensions:mpi-isend-anything
                                             r
                                             right-neighbor :tag 2)
                                            )))
                                        ((and
                                          (= left-neighbor -1)
                                          (not (= right-neighbor -1)))
                                         (let ((r (active-filter right-filter)))
                                           (cl-mpi-extensions:mpi-waitall-anything
                                            (cl-mpi-extensions:mpi-irecv-anything right-neighbor :tag 1)
                                            (cl-mpi-extensions:mpi-isend-anything
                                             r
                                             right-neighbor :tag 2)
                                            )))
                                        ((and
                                          (not (= left-neighbor -1))
                                          (= right-neighbor -1))
                                         (let ((l (active-filter left-filter))
                                               )
                                           (cl-mpi-extensions:mpi-waitall-anything
                                            (cl-mpi-extensions:mpi-irecv-anything left-neighbor :tag 2)
                                            (cl-mpi-extensions:mpi-isend-anything
                                             l
                                             left-neighbor :tag 1)
                                            )))
                                        (t nil))))
                               (loop for packet in recv
                                     do
                                        (destructuring-bind (rank tag object) packet
                                          (when object
                                            (funcall func object))))))))))))))))

(defun exchange-nodes (sim func)
  (declare (function func))
  (exchange-node-like
   sim
   #'serialise-node
   #'deserialise-node
   func))

(defun exchange-domain-bounds (sim)
  (with-accessors ((mesh sim-mesh)
                   (domain-bounds mpm-sim-mpi-domain-bounds ))
      sim
      (with-accessors ((h mesh-resolution))
          mesh
          (setf domain-bounds
                (mapcar (lambda (v) (mapcar (lambda (x) (* (round x h) h)) v))
                        domain-bounds)))))

;; (defun exchange-nodes-nonblocking (sim func)
;;   (declare (function func))
;;   (let* ((rank (cl-mpi:mpi-comm-rank))
;;          (size (cl-mpi:mpi-comm-size)))
;;     (with-accessors ((mesh cl-mpm:sim-mesh))
;;         sim
;;       (let* ((nd-nodes (cl-mpm/mesh:mesh-nodes mesh))
;;             (index (mpi-rank-to-index sim rank))
;;             (bounds-list (mpm-sim-mpi-domain-bounds sim))
;;             (h (cl-mpm/mesh:mesh-resolution mesh))
;;             (halo-depth 2))


;;         (let* ((recv (cl-mpi-extensions:mpi-waitall-anything
;;                      (array-operations/utilities:nested-loop
;;                          (dx dy dz) '(2 2 2)
;;                        (let* ((tag (+ dx (* dy 2) (* 4 dz)))
;;                               (send-neighbor (mpi-index-to-rank sim (mapcar (lambda (a b) (declare (fixnum a b)) (+ a b)) index (list dx dy dz)))))
;;                          (declare (fixnum send-neighbor))
;;                          (destructuring-bind (bl bu) (nth i bounds-list)
;;                            (when (not
;;                                   (and
;;                                    (= left-neighbor -1)
;;                                    (= right-neighbor -1)))
;;                              (labels
;;                                  ((active-filter (nodes)
;;                                     (let ((res
;;                                             (lparallel:premove-if-not
;;                                              (lambda (mp)
;;                                                (and
;;                                                 (cl-mpm/mesh:node-active mp)))
;;                                              nodes)))
;;                                       (make-array (length res) :initial-contents res))))
;;                                (let ((left-filter (nth 0 (nth i (mpm-sim-mpi-halo-node-list sim))))
;;                                      (right-filter (nth 1 (nth i (mpm-sim-mpi-halo-node-list sim))))
;;                                      )
;;                                  (declare (fixnum left-neighbor right-neighbor))
;;                                  ;; (format t "Rank ~D - left ~A~%" rank (length left-filter))
;;                                  ;; (format t "Rank ~D - righ ~A~%" rank (length right-filter))
;;                                  (let* ((cl-mpi-extensions::*standard-encode-function* #'serialise-node)
;;                                         (cl-mpi-extensions::*standard-decode-function* #'deserialise-node))
;;                                    (cond
;;                                      ((and (not (= left-neighbor -1))
;;                                            (not (= right-neighbor -1))
;;                                            )
;;                                       (cl-mpi-extensions:mpi-irecv-anything right-neighbor :tag 1)
;;                                       (cl-mpi-extensions:mpi-irecv-anything left-neighbor :tag 2)
;;                                       (cl-mpi-extensions:mpi-isend-anything
;;                                        (active-filter left-filter)
;;                                        left-neighbor :tag 1)
;;                                       (cl-mpi-extensions:mpi-isend-anything
;;                                        (active-filter right-filter)
;;                                        right-neighbor :tag 2)
;;                                       )
;;                                      ((and
;;                                        (= left-neighbor -1)
;;                                        (not (= right-neighbor -1))
;;                                        )
;;                                       (cl-mpi-extensions:mpi-irecv-anything right-neighbor :tag 1)
;;                                       (cl-mpi-extensions:mpi-isend-anything
;;                                        (active-filter right-filter)
;;                                        right-neighbor :tag 2)
;;                                       )
;;                                      ((and
;;                                        (not (= left-neighbor -1))
;;                                        (= right-neighbor -1))
;;                                       (cl-mpi-extensions:mpi-irecv-anything left-neighbor :tag 2)
;;                                       (cl-mpi-extensions:mpi-isend-anything
;;                                        (active-filter left-filter)
;;                                        left-neighbor :tag 1)
;;                                       )
;;                                      (t nil))
;;                                    )))))

;;                          )
                       
;;                        (incf tag)
;;                        )
;;                     (loop for i from 0 to 2
;;                           do
;;                              (let ((id-delta (list 0 0 0)))
;;                                (setf (nth i id-delta) 1)
;;                                (let (
;;                                      (right-neighbor (mpi-index-to-rank sim (mapcar (lambda (a b) (declare (fixnum a b)) (+ a b)) index id-delta))))
;;                                  ;; (format t "Bounds list ~A~%" bounds-list)
;;                                  ))))))
;;           (loop for packet in recv
;;                 do
;;                    (destructuring-bind (rank tag object) packet
;;                      (when object
;;                        (funcall func object))))
;;           )))))





(declaim (notinline mpi-sync-momentum))
(defun mpi-sync-momentum (sim)
  ;; (format t "Sync momentum ~%")
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
      (exchange-nodes
       sim
       (lambda (node-list)
         (lparallel:pdotimes (i (length node-list))
           (let* ((mpi-node (aref node-list i))
                  (index (mpi-object-node-index mpi-node))
                  (node (cl-mpm/mesh:get-node mesh index)))
             (if node
                 (progn
                   (with-accessors ((active cl-mpm/mesh:node-active)
                                    (mass cl-mpm/mesh:node-mass)
                                    (velocity cl-mpm/mesh:node-velocity)
                                    (pmod cl-mpm/mesh::node-pwave)
                                    (vol cl-mpm/mesh::node-volume)
                                    (svp cl-mpm/mesh::node-svp-sum))
                       node
                     (declare (double-float mass svp vol pmod))
                     (setf active t)
                     (incf mass (the double-float (mpi-object-node-mass mpi-node)))
                     (incf svp (the double-float (mpi-object-node-svp mpi-node)))
                     (incf vol (the double-float (mpi-object-node-vol mpi-node)))
                     (incf pmod (the double-float (mpi-object-node-pmod mpi-node)))
                     (cl-mpm/fastmaths::fast-.+ velocity (mpi-object-node-velocity mpi-node) velocity)
                     ))
                 (error "MPI exchange touched invalid node?"))))))))

(defun mpi-sync-j-inc (sim)
  ;; (format t "Sync momentum ~%")
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
      (exchange-nodes
       sim
       (lambda (node-list)
         (lparallel:pdotimes (i (length node-list))
           (let* ((mpi-node (aref node-list i))
                  (index (mpi-object-node-index mpi-node))
                  (node (cl-mpm/mesh:get-node mesh index)))
             (if node
                 (progn
                   (with-accessors ((j-inc cl-mpm/mesh::node-jacobian-inc))
                       node
                     (declare (double-float j-inc))
                     (incf j-inc (the double-float (mpi-object-node-j-inc mpi-node)))))
                 (error "MPI exchange touched invalid node?"))))))))


(declaim (notinline mpi-sync-force))
(defun mpi-sync-force (sim)
  ;; (format t "Sync force ~%")
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (exchange-nodes
     sim
     (lambda (node-list)
       (lparallel:pdotimes (i (length node-list))
         (let* ((mpi-node (aref node-list i))
                (index (mpi-object-node-index mpi-node))
                (node (cl-mpm/mesh:get-node mesh index)))
           (if node
               (when (cl-mpm/mesh:node-active node)
                 (with-accessors ((mass cl-mpm/mesh:node-mass)
                                  (force cl-mpm/mesh:node-force)
                                  (force-int cl-mpm/mesh::node-internal-force)
                                  (force-ext cl-mpm/mesh::node-external-force)
                                  (force-damping cl-mpm/mesh::node-damping-force)
                                  (force-buoyancy cl-mpm/mesh::node-buoyancy-force))
                     node
                   (cl-mpm/fastmaths::fast-.+ force-int (mpi-object-node-force-int mpi-node) force-int)
                   (cl-mpm/fastmaths::fast-.+ force-ext (mpi-object-node-force-ext mpi-node) force-ext)
                   (cl-mpm/fastmaths::fast-.+ force-damping (mpi-object-node-force-damping mpi-node) force-damping)
                   (cl-mpm/fastmaths::fast-.+ force-buoyancy (mpi-object-node-force-buoyancy mpi-node) force-buoyancy)
                   ))
               (error "MPI force exchange touched invalid node ~A" index))))))))

(defparameter *damage-mp-send-cache* (make-array 0 :element-type 'cl-mpm::particle :adjustable t :fill-pointer 0))
(declaim (notinline mpi-sync-damage-mps))
(defun mpi-sync-damage-mps (sim &optional halo-depth)
  (let* ((rank (cl-mpi:mpi-comm-rank))
         (size (cl-mpi:mpi-comm-size)))
    (with-accessors ((mps cl-mpm:sim-mps)
                     (mesh cl-mpm:sim-mesh))
        sim
      (let ((all-mps mps)
            (index (mpi-rank-to-index sim rank))
            (bounds-list (mpm-sim-mpi-domain-bounds sim))
            (halo-depth (if halo-depth
                            halo-depth
                            1d0))
            (nd (cl-mpm/mesh:mesh-nd mesh))
            (damage-mps (mpm-sim-mpi-damage-mps-cache sim)))
        ;; (pprint halo-depth)
        (setf (fill-pointer damage-mps) 0)
        (setf (fill-pointer *damage-mp-send-cache*) 0)
        (loop for i from 0 below nd
              do
                 (let ((id-delta (list 0 0 0)))
                   (setf (nth i id-delta) 1)
                   (let ((left-neighbor (mpi-index-to-rank sim (mapcar #'- index id-delta)))
                         (right-neighbor (mpi-index-to-rank sim (mapcar #'+ index id-delta))))
                     (declare (double-float halo-depth))
                     (destructuring-bind (bl bu) (nth i bounds-list)
                       (declare (double-float bl bu))
                       (labels
                           ((halo-filter (test)
                              (let ((res
                                      (lparallel:premove-if-not
                                       (lambda (mp)
                                         (funcall test (cl-mpm/utils:varef (cl-mpm/particle:mp-position mp) i)))
                                       all-mps))
                                    (res-corners
                                      (lparallel:premove-if-not
                                       (lambda (mp)
                                         (funcall test (cl-mpm/utils:varef (cl-mpm/particle:mp-position mp) i)))
                                       damage-mps)))
                                (concatenate `(vector ,(array-element-type res)) res res-corners)))
                            (left-filter ()
                              (halo-filter (lambda (pos)
                                             (declare (double-float pos))
                                             (and
                                              (<= pos (+ bl halo-depth))))))
                            (right-filter ()
                              (halo-filter (lambda (pos)
                                             (declare (double-float pos))
                                             (and
                                              (> pos (- bu halo-depth)))))))
                         (let ((l (left-filter))
                               (r (right-filter)))
                           ;; (format t "Rank ~D - ~D - Sending ~D damage mps left~%"  rank i (length l))
                           ;; (format t "Rank ~D - ~D - Sending ~D damage mps right~%" rank i (length r))
                           ;; (unless (= left-neighbor -1)
                           ;;   (iterate-over-mps
                           ;;    l (lambda (mp) (setf (cl-mpm/particle::mp-index mp) 1))))
                           ;; (unless (= right-neighbor -1)
                           ;;   (iterate-over-mps
                           ;;    r (lambda (mp) (setf (cl-mpm/particle::mp-index mp) 2))))
                           )
                         (when t
                           (let* ((cl-mpi-extensions::*standard-encode-function* #'serialise-damage-mp)
                                  (cl-mpi-extensions::*standard-decode-function* #'deserialise-damage-mp)
                                  (recv
                                    (cond
                                      ((and (not (= left-neighbor -1))
                                            (not (= right-neighbor -1)))
                                       (let ((l (left-filter))
                                             (r (right-filter)))
                                         ;; (format t "Rank ~D - ~D - Sending ~D damage mps left~%"  rank i (length l))
                                         ;; (format t "Rank ~D - ~D - Sending ~D damage mps right~%" rank i (length r))
                                         (cl-mpi-extensions:mpi-waitall-anything
                                          (cl-mpi-extensions:mpi-irecv-anything right-neighbor :tag 1)
                                          (cl-mpi-extensions:mpi-irecv-anything left-neighbor :tag 2)
                                          (cl-mpi-extensions:mpi-isend-anything
                                           l
                                           left-neighbor :tag 1)
                                          (cl-mpi-extensions:mpi-isend-anything
                                           r
                                           right-neighbor :tag 2)
                                          )))
                                      ((and
                                        (= left-neighbor -1)
                                        (not (= right-neighbor -1)))
                                       (let ((r (right-filter)))
                                         ;; (format t "Rank ~D - ~D - Sending ~D damage mps left~%"  rank i 0)
                                         ;; (format t "Rank ~D - ~D - Sending ~D damage mps right~%" rank i (length r))
                                         (cl-mpi-extensions:mpi-waitall-anything
                                          (cl-mpi-extensions:mpi-irecv-anything right-neighbor :tag 1)
                                          (cl-mpi-extensions:mpi-isend-anything
                                           r
                                           right-neighbor :tag 2))))
                                      ((and
                                        (not (= left-neighbor -1))
                                        (= right-neighbor -1))
                                       (let ((l (left-filter)))
                                         ;; (format t "Rank ~D - ~D - Sending ~D damage mps left~%"  rank i (length l))
                                         ;; (format t "Rank ~D - ~D - Sending ~D damage mps right~%" rank i 0)
                                         (cl-mpi-extensions:mpi-waitall-anything
                                          (cl-mpi-extensions:mpi-irecv-anything left-neighbor :tag 2)
                                          (cl-mpi-extensions:mpi-isend-anything
                                           l
                                           left-neighbor :tag 1))))
                                      (t nil))))
                             (loop for packet in recv
                                   do
                                      (destructuring-bind (rank tag object) packet
                                        (when object
                                          (loop for mp across object
                                                do (progn
                                                     (vector-push-extend
                                                      (make-instance 'cl-mpm/particle::particle-damage-mpi
                                                                     :damage (mpi-object-damage-mp-damage mp)
                                                                     :volume (mpi-object-damage-mp-volume mp)
                                                                     :position (mpi-object-damage-mp-position mp)
                                                                     :damage-y (mpi-object-damage-mp-y mp)
                                                                     :local-length (mpi-object-damage-mp-local-length mp))
                                                      damage-mps)))))))))))))
        damage-mps))))



;; (defun test-mpi-mp-sync (sim)
;;   (let ((mpi-objects (cl-mpm/mpi::deserialise-damage-mp
;;                       (cl-mpm/mpi::serialise-damage-mp
;;                        (remove-if (lambda (mp) (= 0d0 (cl-mpm/particle::mp-damage mp))) (cl-mpm:sim-mps sim)))))
;;         (output (make-array 0 :adjustable t :fill-pointer 0)))
;;     (loop for mp across mpi-objects
;;           do (progn
;;                (vector-push-extend
;;                 (make-instance 'cl-mpm/particle::particle-damage
;;                                :nd 2
;;                                :volume (mpi-object-damage-mp-volume mp)
;;                                :position (mpi-object-damage-mp-position mp)
;;                                :damage-y (mpi-object-damage-mp-y mp)
;;                                :local-length-t (mpi-object-damage-mp-local-length mp))
;;                 output)
;;                ))
;;     output))

(defun update-stress-nodal-fbar (sim)
  "Update all stresses, with optional f-bar"
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps)
                   (dt cl-mpm:sim-dt)
                   (fbar cl-mpm::sim-enable-fbar))
      sim
    (declare ((array cl-mpm/particle:particle) mps) (cl-mpm/mesh::mesh mesh))
    (iterate-over-mps
     mps
     (lambda (mp)
       (cl-mpm::calculate-strain-rate mesh mp dt)
       (cl-mpm::map-jacobian mesh mp dt)))
    (mpi-sync-j-inc sim)
    (iterate-over-mps
     mps
     (lambda (mp)
       (cl-mpm::update-stress-mp mesh mp dt fbar)
       (cl-mpm::post-stress-step mesh mp dt)
       )))
  (values))

(defmethod cl-mpm::update-sim ((sim mpm-sim-mpi-nodes))
  (with-slots ((mesh cl-mpm::mesh)
               (mps  cl-mpm::mps)
               (bcs  cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (bcs-force-list cl-mpm::bcs-force-list)
               (dt cl-mpm::dt)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (nonlocal-damage cl-mpm::nonlocal-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (fbar cl-mpm::enable-fbar)
               (vel-algo cl-mpm::velocity-algorithm)
               )
                sim
    (declare (type double-float mass-filter))
                (progn
                    ;; (exchange-mps sim)
                  (when t;(> (length mps) 0)
                    (cl-mpm::reset-grid mesh)
                    (cl-mpm::p2g mesh mps)
                    (mpi-sync-momentum sim)
                    (when (> mass-filter 0d0)
                      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
                    (cl-mpm::update-node-kinematics mesh dt)
                    (cl-mpm::apply-bcs mesh bcs dt)
                    (cl-mpm::update-stress mesh mps dt fbar)
                    (cl-mpm::p2g-force mesh mps)
                    (loop for bcs-f in bcs-force-list
                          do (cl-mpm::apply-bcs mesh bcs-f dt))
                    (mpi-sync-force sim)
                    (cl-mpm::update-node-forces sim)
                    (cl-mpm::apply-bcs mesh bcs dt)
                    (cl-mpm::g2p mesh mps dt vel-algo)
                    (cl-mpm::update-particles sim)
                    (when remove-damage
                      (cl-mpm::remove-material-damaged sim))
                    (when split
                      (cl-mpm::split-mps sim))
                    (cl-mpm::check-mps sim)
                    (set-mp-mpi-index sim)
                    )
                  (exchange-mps sim 0d0)
                  (set-mp-mpi-index sim)
                  (clear-ghost-mps sim)
                    )))



(defmethod cl-mpm::update-sim ((sim mpm-sim-mpi-nodes-damage))
  (with-slots ((mesh cl-mpm::mesh)
               (mps  cl-mpm::mps)
               (bcs  cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (bcs-force-list cl-mpm::bcs-force-list)
               (dt cl-mpm::dt)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (nonlocal-damage cl-mpm::nonlocal-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (fbar cl-mpm::enable-fbar)
               (vel-algo cl-mpm::velocity-algorithm)
               )
                sim
    (declare (type double-float mass-filter))
                (progn
                  (cl-mpm::reset-grid mesh)
                  (cl-mpm::p2g mesh mps)
                  (mpi-sync-momentum sim)
                  (when (> mass-filter 0d0)
                    (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
                  (cl-mpm::update-node-kinematics mesh dt)
                  (cl-mpm::apply-bcs mesh bcs dt)
                  (cl-mpm::update-stress mesh mps dt fbar)
                  ;; (update-stress-nodal-fbar sim)
                  (cl-mpm/damage::calculate-damage sim)
                  (cl-mpm::p2g-force mesh mps)
                  (loop for bcs-f in bcs-force-list
                        do (cl-mpm::apply-bcs mesh bcs-f dt))
                  (mpi-sync-force sim)
                  (cl-mpm::update-node-forces sim)
                  (cl-mpm::apply-bcs mesh bcs dt)
                  (cl-mpm::update-dynamic-stats sim)
                  (cl-mpm::g2p mesh mps dt vel-algo)
                  (cl-mpm::update-particles sim)
                  (when remove-damage
                    (cl-mpm::remove-material-damaged sim))
                  (when split
                    (cl-mpm::split-mps sim))
                  (cl-mpm::check-mps sim)
                  ;; (cl-mpm::check-single-mps sim)
                  ;; (set-mp-mpi-index sim)
                  (exchange-mps sim 0d0)
                  (set-mp-mpi-index sim)
                  (clear-ghost-mps sim)
                    )))
(defmacro rank-0-time (rank &rest body)
  `(if (= ,rank 0)
       (time
        (progn
          ,@body))
       (progn
         ,@body)))

(defmethod cl-mpm::update-sim ((sim mpm-sim-usl-mpi-nodes-damage))
  (with-slots ((mesh cl-mpm::mesh)
               (mps  cl-mpm::mps)
               (bcs  cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (bcs-force-list cl-mpm::bcs-force-list)
               (dt cl-mpm::dt)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (nonlocal-damage cl-mpm::nonlocal-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (fbar cl-mpm::enable-fbar)
               (vel-algo cl-mpm::velocity-algorithm)
               )
                sim
    (declare (type double-float mass-filter))
                (progn
                  (cl-mpm::reset-grid mesh)
                  (cl-mpm::p2g mesh mps)
                  (mpi-sync-momentum sim)
                  (when (> mass-filter 0d0)
                    (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
                  (cl-mpm::update-node-kinematics mesh dt)
                  (cl-mpm::apply-bcs mesh bcs dt)
                  (cl-mpm::p2g-force mesh mps)
                  (loop for bcs-f in bcs-force-list
                        do
                           (cl-mpm::apply-bcs mesh bcs-f dt))
                  (mpi-sync-force sim)
                  (cl-mpm::update-node-forces sim)
                  (cl-mpm::apply-bcs mesh bcs dt)
                  (cl-mpm::g2p mesh mps dt vel-algo)


                  ;;2nd round of mapping for USL
                  (cl-mpm::reset-grid-velocity mesh)
                  (cl-mpm::p2g mesh mps)
                  (mpi-sync-momentum sim)
                  (when (> mass-filter 0d0)
                    (cl-mpm::filter-grid-velocity mesh (cl-mpm::sim-mass-filter sim)))
                  (cl-mpm::update-node-kinematics mesh dt)
                  (cl-mpm::apply-bcs mesh bcs dt)

                  (cl-mpm::update-stress mesh mps dt fbar)
                  (cl-mpm/damage::calculate-damage sim)
                  (cl-mpm::update-dynamic-stats sim)
                  (cl-mpm::update-particles sim)

                  (when remove-damage
                    (cl-mpm::remove-material-damaged sim))
                  (when split
                    (cl-mpm::split-mps sim))

                  ;;Must be done here
                  (set-mp-mpi-index sim)
                  (exchange-mps sim 0d0)
                  (set-mp-mpi-index sim)
                  (clear-ghost-mps sim)
                  (cl-mpm::check-mps sim)
                  (cl-mpm::check-single-mps sim)

                    )))

(defmethod cl-mpm::update-sim ((sim mpm-sim-mpi))
  (with-slots ((mesh cl-mpm::mesh)
               (mps  cl-mpm::mps)
               (bcs  cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (bcs-force-list cl-mpm::bcs-force-list)
               (dt cl-mpm::dt)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (nonlocal-damage cl-mpm::nonlocal-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (fbar cl-mpm::enable-fbar)
               (vel-algo cl-mpm::velocity-algorithm)
               )
                sim
    (declare (type double-float mass-filter))
                (progn
                    (exchange-mps sim)
                  (when (> (length mps) 0)
                    (cl-mpm::reset-grid mesh)
                    (cl-mpm::p2g mesh mps)
                    (when (> mass-filter 0d0)
                      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
                    (cl-mpm::update-node-kinematics mesh dt)
                    (cl-mpm::apply-bcs mesh bcs dt)
                    ;; (cl-mpm::update-stress mesh mps dt)
                    (cl-mpm::update-stress mesh mps dt fbar)
                                        ;(exchange-mps sim)
                    (cl-mpm/damage::calculate-damage sim)
                                        ;(exchange-mps sim)
                    (cl-mpm::p2g-force mesh mps)
                    (loop for bcs-f in bcs-force-list
                          do
                             (cl-mpm::apply-bcs mesh bcs-f dt))
                    (cl-mpm::update-node-forces sim)
                    (cl-mpm::apply-bcs mesh bcs dt)
                                        ;Also updates mps inline
                    (cl-mpm::g2p mesh mps dt vel-algo)
                    ;;MPI reduce new velocities
                                        ;(exchange-mps sim)
                    ;;Get new MPS

                    (cl-mpm::update-particles sim)
                    (when remove-damage
                      (cl-mpm::remove-material-damaged sim))
                    (when split
                      (cl-mpm::split-mps sim))
                    (cl-mpm::check-mps sim)
                    (set-mp-mpi-index sim))
                    ;; (clear-ghost-mps sim)
                    ;;Update mp list between processors
                    )))

(defun clear-ghost-mps (sim)
  (let ((rank (cl-mpi:mpi-comm-rank)))
     (cl-mpm::remove-mps-func
      sim
      (lambda (mp)
        (not (= rank (cl-mpm/particle::mp-mpi-index mp)))))))

;; (defun store-array (obj stream)
;;   (declare (optimize speed (safety 0) (debug 0))
;;            (type array obj))
;;   (cl-store:output-type-code cl-store::+array-code+ stream)
;;   (if (and (= (array-rank obj) 1)
;;            (array-has-fill-pointer-p obj))
;;       (cl-store:store-object (fill-pointer obj) stream)
;;       (cl-store:store-object nil stream))
;;   (cl-store:store-object (array-element-type obj) stream)
;;   (cl-store:store-object (adjustable-array-p obj) stream)
;;   (cl-store:store-object (array-dimensions obj) stream)
;;   (dolist (x (multiple-value-list (array-displacement obj)))
;;     (cl-store:store-object x stream))
;;   (cl-store:store-object (array-total-size obj) stream)
;;   (loop for x from 0 below (array-total-size obj) do
;;     (cl-store:store-object (row-major-aref obj x) stream)))

;; (cl-store:defstore-cl-store (obj array stream)
;;   (store-array))

(defun ser-part (mps &optional parts)
 ; (check-type size fixnum)
  (let ((size (length mps)))
    (when (plusp size)
      (let ()
        (flet ((compute-part (part-offset part-size)
                 (let ((index part-offset)
                       (end (+ part-offset part-size)))
                   (declare (type fixnum index end part-offset part-size))
                   (let ((cl-store:*current-backend* cl-store:*default-backend*)
                         (cl-store:*check-for-circs* nil))
                     (flexi-streams:with-output-to-sequence (stream :element-type '(unsigned-byte 8))
                       (loop while (< index end)
                             do (cl-store:store-object (aref mps index) stream)
                                ;(hu.dwim.serializer:serialize (aref mps index) :output stream)
                                (incf index))))
                   ;; (loop while (< index end)
                   ;;       do (funcall fn index)
                   ;;          (incf index))
                   )))
          (let ((parts (lparallel.cognate::get-parts-hint parts))
                (channel (lparallel.cognate::make-channel)))
            (lparallel.cognate::with-parts size parts
              (loop while (lparallel.cognate::next-part)
                    do (lparallel.cognate::submit-task
                        channel #'compute-part
                        (lparallel.cognate::part-offset) (lparallel.cognate::part-size)))
              (loop repeat (lparallel.cognate::num-parts)
                    collect (lparallel.cognate::receive-result channel))
              ;; (lparallel.cognate::repeat (lparallel.cognate::num-parts)
              ;;   (lparallel.cognate::receive-result channel))
              )))))))







;; (defstruct mpi-node
;;   (index nil :type list)
;;   (mass 0d0 :type double-float)
;;   (pmod 0d0 :type double-float)
;;   (svp  0d0 :type double-float)
;;   (vol  0d0 :type double-float)
;;   (velocity (vector-zeros) :type magicl:matrix/double-float)
;;   (force  (vector-zeros) :type magicl:matrix/double-float)
;;    )
;; (declaim (inline make-mpi-node))



(defun serialise-mps (mps)
  (cl-store-encoder mps))

(defun cl-store-decoder (x)
    (flexi-streams:with-input-from-sequence (stream x)
      (cl-store:restore stream)))

(defun cl-store-encoder (x)
  (let ((res (flexi-streams:with-output-to-sequence (stream)
               (cl-store:store x stream))))
    (static-vectors:make-static-vector (length res)
                                       :element-type '(unsigned-byte 8)
                                       :initial-contents res
                                       )))
(defun deserialise-mps (x)
  (when x
    (flexi-streams:with-input-from-sequence (stream x)
      (cl-store:restore stream))))


(defun test-ser (mps)
  (let ((res (flexi-streams:with-output-to-sequence (stream)
             (cl-store:store mps stream))))
  (static-vectors:make-static-vector (length res)
                                     :element-type '(unsigned-byte 8)
                                     :initial-contents res
                                     )))

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

(declaim (notinline exchange-mps))
(defun exchange-mps (sim &optional halo-depth)
  (let* ((rank (cl-mpi:mpi-comm-rank))
         (size (cl-mpi:mpi-comm-size))
         )
    (with-accessors ((mps cl-mpm:sim-mps)
                     (mesh cl-mpm:sim-mesh))
        sim
      (let ((all-mps mps)
            (index (mpi-rank-to-index sim rank))
            (bounds-list (mpm-sim-mpi-domain-bounds sim))
            (halo-depth (if halo-depth
                            halo-depth
                            (* 2 (cl-mpm/mesh:mesh-resolution mesh))))
            (nd (cl-mpm/mesh:mesh-nd mesh))
            )
        (loop for i from 0 below nd
              do
                 (let ((id-delta (list 0 0 0)))
                   (setf (nth i id-delta) 1)
                   (let ((left-neighbor (mpi-index-to-rank sim (mapcar #'- index id-delta)))
                         (right-neighbor (mpi-index-to-rank sim (mapcar #'+ index id-delta))))
                     (destructuring-bind (bl bu) (nth i bounds-list)
                       (when (not (= bl bu))
                         (labels
                             ((halo-filter (test)
                                (let ((res
                                        (lparallel:premove-if-not
                                         (lambda (mp) (funcall test (magicl:tref (cl-mpm/particle:mp-position mp) i 0)))
                                         all-mps)))
                                  res))
                              (left-filter ()
                                (halo-filter (lambda (pos)
                                               (and
                                                (< pos (+ bl halo-depth))))))
                              (right-filter ()
                                (halo-filter (lambda (pos)
                                               (and
                                                (>= pos (- bu halo-depth)))))))
                           (let* ((cl-mpi-extensions::*standard-encode-function* #'cl-store-encoder)
                                  (cl-mpi-extensions::*standard-decode-function* #'cl-store-decoder)
                                  (recv
                                    (cond
                                      ((and (not (= left-neighbor -1))
                                            (not (= right-neighbor -1)))
                                       (cl-mpi-extensions:mpi-waitall-anything
                                        (cl-mpi-extensions:mpi-irecv-anything right-neighbor :tag 1)
                                        (cl-mpi-extensions:mpi-irecv-anything left-neighbor :tag 2)
                                        (cl-mpi-extensions:mpi-isend-anything
                                         (left-filter)
                                         left-neighbor :tag 1)
                                        (cl-mpi-extensions:mpi-isend-anything
                                         (right-filter)
                                         right-neighbor :tag 2)
                                        ))
                                      ((and
                                        (= left-neighbor -1)
                                        (not (= right-neighbor -1))
                                        )
                                       (cl-mpi-extensions:mpi-waitall-anything
                                        (cl-mpi-extensions:mpi-irecv-anything right-neighbor :tag 1)
                                        (cl-mpi-extensions:mpi-isend-anything
                                         (right-filter)
                                         right-neighbor :tag 2)
                                        ))
                                      ((and
                                        (not (= left-neighbor -1))
                                        (= right-neighbor -1))
                                       (cl-mpi-extensions:mpi-waitall-anything
                                        (cl-mpi-extensions:mpi-irecv-anything left-neighbor :tag 2)
                                        (cl-mpi-extensions:mpi-isend-anything
                                         (left-filter)
                                         left-neighbor :tag 1)
                                        ))
                                      (t nil))))
                             (loop for packet in recv
                                   do
                                      (destructuring-bind (rank tag object) packet
                                        (when (> (length object) 0)
                                          (when object
                                            (loop for mp across object
                                                  do (progn
                                                       (setf (fill-pointer (cl-mpm/particle::mp-cached-nodes mp)) 0)
                                                       (when (slot-exists-p mp 'cl-mpm/particle::damage-position)
                                                         (setf (cl-mpm/particle::mp-damage-position mp) nil))
                                                       (vector-push-extend mp mps)
                                                       )))))))))))))))))
  ;; (cl-mpi:mpi-barrier)
;  )

(defvar *mutex-code* (cl-store:register-code 110 'sb-thread:mutex))
(cl-store:defstore-cl-store (obj sb-thread:mutex stream)
   (cl-store:output-type-code *mutex-code* stream))
(cl-store:defrestore-cl-store (sb-thread:mutex stream)
   (sb-thread:make-mutex))

;; (defvar *node-cache-code* (cl-store:register-code 112 'cl-mpm/particle::node-cache))
;; (cl-store:defstore-cl-store (obj cl-mpm/particle::node-cache stream)
;;   (cl-store:output-type-code *node-cache-code* stream)
;;   )
;; (cl-store:defrestore-cl-store (sb-thread:mutex stream)
;;   ;(sb-thread:make-mutex)
;;   )


(defvar *mesh-code* (cl-store:register-code 111 'cl-mpm/mesh::mesh))
(cl-store:defstore-cl-store (obj cl-mpm/mesh::mesh stream)
    (cl-store:output-type-code *mesh-code* stream)
  (cl-store:store-object (cl-mpm/mesh::mesh-nd obj) stream)
  (cl-store:store-object (cl-mpm/mesh::mesh-count obj) stream)
  (cl-store:store-object (cl-mpm/mesh::mesh-mesh-size obj) stream)
  (cl-store:store-object (cl-mpm/mesh::mesh-resolution obj) stream)
  (cl-store:store-object (cl-mpm/mesh::mesh-nodes obj) stream)
  (cl-store:store-object (cl-mpm/mesh::mesh-cells obj) stream)
  )

(cl-store:defrestore-cl-store (cl-mpm/mesh::mesh stream)
    (let ((obj (make-instance 'cl-mpm/mesh::mesh)))
      (setf (cl-mpm/mesh::mesh-nd obj) (cl-store:restore-object stream))
      (setf (cl-mpm/mesh::mesh-count obj) (cl-store:restore-object stream))
      (setf (cl-mpm/mesh::mesh-mesh-size obj)  (cl-store:restore-object stream))
      (setf (cl-mpm/mesh::mesh-resolution obj) (cl-store:restore-object stream))
      (setf (cl-mpm/mesh::mesh-nodes obj) (cl-store:restore-object stream))
      (setf (cl-mpm/mesh::mesh-cells obj) (cl-store:restore-object stream))
      obj))

(defvar *sim-code* (cl-store:register-code 112 'cl-mpm::mpm-sim))
(cl-store:defstore-cl-store (obj cl-mpm::mpm-sim stream)
    (cl-store:output-type-code *sim-code* stream)
  (cl-store:store-object (cl-mpm::sim-dt obj) stream)
  (cl-store:store-object (cl-mpm::sim-mesh obj) stream)
  )
(cl-store:defrestore-cl-store (cl-mpm::mpm-sim stream)
    (let ((obj (make-instance 'cl-mpm::mpm-sim)))
      (setf (cl-mpm::sim-dt obj) (cl-store:restore-object stream))
      (setf (cl-mpm::sim-mesh obj) (cl-store:restore-object stream))
      obj))

;; (defmethod cl-store:serializable-slots-using-class ((object t) (class cl-mpm/mesh::node))
;;   (delete 'cl-mpm/mesh::local-list (call-next-method) :key 'c2mop:slot-definition-name))

(defvar *node-code* (cl-store:register-code 113 'cl-mpm/mesh::node))
(cl-store:defstore-cl-store (obj cl-mpm/mesh::node stream)
  (cl-store:output-type-code *node-code* stream)
  (loop for slot in (cl-store:serializable-slots obj)
        do
           (cond
             ((eq 'cl-mpm/mesh::local-list (sb-mop:slot-definition-name slot))
              nil)
             (t
              (cl-store:store-object (sb-mop:slot-value-using-class 'cl-mpm/mesh::node obj slot) stream)))))

(cl-store:defrestore-cl-store (cl-mpm/mesh::node stream)
  (let ((obj (make-instance 'cl-mpm/mesh::node)))
    (loop for slot in (cl-store:serializable-slots obj)
          do
             (cond
               ((eq 'cl-mpm/mesh::local-list (sb-mop:slot-definition-name slot))
                nil)
               (t
                (setf (sb-mop:slot-value-using-class 'cl-mpm/mesh::node obj slot) (cl-store:restore-object stream)))))
  ;; (setf (cl-mpm::sim-dt obj) (cl-store:restore-object stream))
    obj))

;; (defun mpi-run (setup-func run-func)
;;   (lfarm:broadcast-task (lambda ()
;;                           (progn
;;                             ;; (asdf:compile-system :magicl :force t)
;;                             (ql:quickload :cl-mpm/examples/slump)
;;                             (ql:quickload :cl-mpm/mpi)
;;                             (in-package :cl-mpm/examples/slump)
;;                             (defparameter *local-sim* nil)
;;                             (defparameter *global-dt* 1d0)
;;                             (setf lparallel:*kernel* (lparallel:make-kernel 4))
;;                             t)))

;;   )

(defun collect-servers (n &optional (dir (uiop:getcwd)))
  (let ((servers (with-open-file (s "lfarm_connections") (read s))))
    (format t "~S ~%" servers)
    (defparameter *open-servers* servers)
    (setf lfarm:*kernel* (lfarm:make-kernel servers))
    )
  ;; (asdf:compile-system :magicl)
  (print "Broadcasting setup info")
  (lfarm:broadcast-task (lambda ()
                          (progn
                            ;; (asdf:compile-system :magicl :force t)
                            (ql:quickload :cl-mpm/examples/slump)
                            (ql:quickload :cl-mpm/mpi)
                            (in-package :cl-mpm/examples/slump)
                            (defparameter *local-sim* nil)
                            (defparameter *global-dt* 1d0)
                            (setf lparallel:*kernel* (lparallel:make-kernel 4))
                            t))))

(defun in-computational-domain (sim pos)
  (with-accessors ((bounds mpm-sim-mpi-domain-bounds))
      sim
      (let ((in-bounds t)
            (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim))))
        (loop for i from 0 below nd
              for bound in bounds
              while in-bounds
              do
                 (destructuring-bind (bl bu) bound
                   (declare (double-float bl bu))
                   (when (not (= bu bl))
                     (setf in-bounds
                           (and
                            in-bounds
                            (and
                             (> bu (cl-mpm/utils:varef pos i))
                             (<= bl (cl-mpm/utils:varef pos i)))
                            )))))
        in-bounds
        )))

(defun in-computational-domain-buffer (sim pos node-buffer)
  (let ((in-bounds t)
        (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim))))
    (loop for i from 0 below 3;(cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)) 
          do
             (destructuring-bind (bl bu) (nth i (mpm-sim-mpi-domain-bounds sim))
               (when (not (= bu bl))
                 (setf in-bounds
                       (and
                        in-bounds
                        (and
                         (> (+ bu (* node-buffer h)) (cl-mpm/utils:varef pos i))
                         (<= (- bl (* node-buffer h)) (cl-mpm/utils:varef pos i)))
                        )))))
    in-bounds
    ))

(defun calculate-domain-sizes (sim &optional size)
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (divs mpm-sim-mpi-domain-count))
      sim
    (when size
      (setf divs size))
    (let ((mc (cl-mpm/mesh:mesh-mesh-size mesh)))
      ;; (loop for i from 0 below (- rank-count 1)
      ;;       do (incf (nth (mod i 3) divs)))
      ;; divs
      (mapcar #'/ mc divs))))

(declaim (notinline set-mp-mpi-index))
(defun set-mp-mpi-index (sim)
  (let* ((rank (cl-mpi:mpi-comm-rank)))
    (cl-mpm:iterate-over-mps
     (cl-mpm:sim-mps sim)
     (lambda (mp)
       (setf (cl-mpm/particle::mp-mpi-index mp)
             (if (in-computational-domain sim (cl-mpm/particle:mp-position mp))
                 rank
                 -1))))))

;; (defun cl-mpi:mpi-comm-rank ()
;;   0)
;; (defun cl-mpi:mpi-comm-size ()
;;  4)


(defun setup-domain-bounds (sim &key (domain-scaler (lambda (x) x)))

  (with-accessors ((mesh-size cl-mpm/mesh::mesh-mesh-size)
                   (mesh-count cl-mpm/mesh::mesh-count)
                   (h cl-mpm/mesh::mesh-resolution))
      (cl-mpm:sim-mesh sim)
    (let* ((rank (cl-mpi:mpi-comm-rank))
           (count (cl-mpi:mpi-comm-size))
           (index (mpi-rank-to-index sim rank))
           (comp-size (mpm-sim-mpi-domain-count sim)))
      (setf (mpm-sim-mpi-domain-bounds sim)
            (mapcar (lambda (domain size) (mapcar (lambda (x) (* x size)) domain))
                    (funcall domain-scaler
                             (loop for i from 0 to 2
                                   collect
                                   (if (> (nth i comp-size) 1)
                                       (list (/ (* (/ (nth i mesh-size) (nth i comp-size)) (nth i index)) (nth i mesh-size))
                                             (/ (* (/ (nth i mesh-size) (nth i comp-size)) (+ (nth i index) 1)) (nth i mesh-size)))
                                       (list 0d0 1d0))))
                    mesh-size))
      (format t "Rank ~D: Domain bounds ~A~%" rank (mpm-sim-mpi-domain-bounds sim))
      (exchange-domain-bounds sim)
      (format t "Rank ~D: Exchanged domain bounds ~A~%" rank (mpm-sim-mpi-domain-bounds sim))))
  )


(defun origin-domain-scaler (sim real-x real-y real-z)
  "Build a domain scalar that re-scales all mpi domains to be tight packed from (0 - real-x) (0 - real-y) (0 - real-z)
leaves a hanging mpi domain at the back"
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (with-accessors ((size cl-mpm/mesh:mesh-mesh-size))
        mesh
      (let ((target-size (list real-x real-y real-z)))
        (lambda (domain)
          (loop for i from 0
                for dim in domain
                collect
                (mapcar (lambda (x)
                          (if (and (> x 0d0)
                                   (< x 1d0)
                                   (nth i target-size)
                                   )
                              (* x (/ (nth i target-size) (nth i size)))
                              x))
                        dim)))))))

;; (defun collect-in-bound-nodes (sim)
;;   (let (()))
;;   (cl-mpm:iterate-over-nodes-serial
;;    (cl-mpm:sim-mesh sim)

;;    )
;;   )

(defun domain-decompose (sim &key (domain-scaler (lambda (x) x)))
  (%domain-decompose sim domain-scaler)
  )
(defgeneric %domain-decompose (sim domain-scaler)
  )
(defmethod %domain-decompose (sim domain-scaler)
  "The aim of domain decomposition is to take a full simulation and cut it into subsections for MPI"

  (with-accessors ((mesh-size cl-mpm/mesh::mesh-mesh-size)
                   (mesh-count cl-mpm/mesh::mesh-count)
                   (h cl-mpm/mesh::mesh-resolution))
      (cl-mpm:sim-mesh sim)
    (let* ((rank (cl-mpi:mpi-comm-rank))
           (count (cl-mpi:mpi-comm-size))
           (index (mpi-rank-to-index sim rank))
           (x-length (first mesh-size))
           (slice-size (/ x-length count))
           (slice-count count)
           (bound-lower (* rank slice-size))
           (bound-upper (* (+ 1 rank) slice-size))
           (comp-size (mpm-sim-mpi-domain-count sim)))
      (unless (mpm-sim-mpi-domain-bounds sim)
        (setup-domain-bounds sim :domain-scaler domain-scaler))
      (setf (mpm-sim-mpi-domain-index sim) index)
      (set-mp-mpi-index sim)
      (clear-ghost-mps sim)
      ;; (format t "Sim MPs: ~a~%" (length (cl-mpm:sim-mps sim)))
      (let ((min-mps (mpi-min (length (cl-mpm:sim-mps sim))))
            (max-mps (mpi-max (length (cl-mpm:sim-mps sim)))))
        (when (and (= rank 0)
                   (> min-mps 0))
          (format t "Occupancy ratio : ~F%~%" (* 100d0 (/ max-mps min-mps))))))
    ))
(defparameter *prune-nodes* t)
(defmethod %domain-decompose :after ((sim cl-mpm/mpi::mpm-sim-mpi-nodes) domain-scaler)
  (let* ((rank (cl-mpi:mpi-comm-rank))
         (size (cl-mpi:mpi-comm-size)))
    (with-accessors ((mesh cl-mpm:sim-mesh))
        sim
      (let* ((nd-nodes (cl-mpm/mesh:mesh-nodes mesh))
             (all-nodes (make-array (array-total-size nd-nodes) :displaced-to nd-nodes :displaced-index-offset 0))
             (index (mpi-rank-to-index sim rank))
             (bounds-list (mpm-sim-mpi-domain-bounds sim))
             (h (cl-mpm/mesh:mesh-resolution mesh))
             (halo-depth (cl-mpm/mpi::mpm-sim-mpi-halo-depth sim)))
        (loop for i from 0 to 2
              do
                 (let ((id-delta (list 0 0 0)))
                   (setf (nth i id-delta) 1)
                   (let ((left-neighbor (mpi-index-to-rank sim (mapcar #'- index id-delta)))
                         (right-neighbor (mpi-index-to-rank sim (mapcar #'+ index id-delta))))
                     (destructuring-bind (bl bu) (nth i bounds-list)
                       (labels
                           ((halo-filter (test)
                              (let ((res
                                      (lparallel:premove-if-not
                                       (lambda (node)
                                         (if node
                                             (funcall test (nth i (cl-mpm/mesh:node-index node)))
                                             nil))
                                       all-nodes)))
                                (make-array (length res) :initial-contents res)))
                            (left-filter ()
                              (halo-filter (lambda (pos)
                                             (and
                                              (>= pos (- (/ bl h) halo-depth))
                                              (<= pos (+ (/ bl h) halo-depth))
                                              )
                                             )))
                            (right-filter ()
                              (halo-filter (lambda (pos)
                                             (and
                                              (>= pos (- (/ bu h) halo-depth))
                                              (<= pos (+ (/ bu h) halo-depth)))
                                             ))))
                         (when t;(not (= left-neighbor -1))
                           (setf (nth 0 (nth i (mpm-sim-mpi-halo-node-list sim)))
                                 (left-filter)))

                         (when t;(not (= right-neighbor -1))
                           (setf (nth 1 (nth i (mpm-sim-mpi-halo-node-list sim)))
                                 (right-filter))))))))

        (when *prune-nodes*
          (with-accessors ((nodes cl-mpm/mesh:mesh-nodes)
                           (cells cl-mpm/mesh::mesh-cells))
              mesh
            (with-accessors ((bcs cl-mpm:sim-bcs))
                sim
              (let* ((domain-sizes (mapcar (lambda (x) (abs (reduce #'- x)))
                                           (cl-mpm/mpi::mpm-sim-mpi-domain-bounds sim)))
                     (h (cl-mpm/mesh:mesh-resolution mesh))
                     (min-domain-length (reduce #'max (remove 0d0 domain-sizes)))
                     (buffer-size (+ 1 (max halo-depth (ceiling min-domain-length h)))))
                (format t "Pruning nodes and BCs up to ~D nodes away~%" buffer-size)
                (format t "Domain size ~A~%" domain-sizes)
                ;;Remove nil bcs

                (setf bcs (delete nil bcs))
                (let ((prune-count 0))
                  (dotimes (i (length bcs))
                    (let ((bc (aref bcs i)))
                      (when bc
                        (let ((index (cl-mpm/bc:bc-index bc)))
                          ;;nil indexes indiciate global bcs
                          (when index
                            (progn
                              ;;SBCL is very sure that the result cannot be nil!
                              ;; (when (not (equal (cl-mpm/mesh:get-node mesh index) nil))
                              ;;   (let ((node (cl-mpm/mesh:get-node mesh index)))
                                  ;; (format t "Pruning BC at index: ~A node: ~A~%" index node)
                                  (when (not (in-computational-domain-buffer
                                              sim
                                              ;; (cl-mpm/mesh::node-position node)
                                              (cl-mpm/utils:vector-from-list  (cl-mpm/mesh:index-to-position mesh index))
                                              buffer-size))
                                    (incf prune-count)
                                    (setf (aref bcs i) nil))
                                        ;))
                              ))))))

                  (format t "Rank ~D - Pruned ~D bcs~%" rank prune-count))
                (setf (cl-mpm:sim-bcs sim) (delete nil bcs))

                ;;Trim out all nodes that we can get rid of
                (let ((prune-count 0))
                  (dotimes (i (array-total-size nodes))
                    (let ((node (row-major-aref nodes i)))
                      (when node
                        (when (not (in-computational-domain-buffer
                                    sim
                                    (cl-mpm/mesh::node-position node)
                                    buffer-size))
                          (incf prune-count)
                          (setf (row-major-aref nodes i) nil)))))
                  (format t "Rank ~D - Pruned ~D nodes~%" rank prune-count))

                ;;Prune orphan bcs
                (let ((prune-count 0))
                  (dotimes (i (length bcs))
                    (let ((bc (aref bcs i)))
                      (when bc
                        (let ((index (cl-mpm/bc:bc-index bc)))
                          ;;nil indexes indiciate global bcs
                          (when index
                            ;;SBCL is very sure that the result cannot be nil!
                            (when (equal (cl-mpm/mesh:get-node mesh index) nil)
                              (incf prune-count)
                              (setf (aref bcs i) nil)))))))
                  (setf (cl-mpm:sim-bcs sim) (delete nil bcs))
                  (format t "Rank ~D - Pruned ~D orphan bcs~%" rank prune-count))

                (loop for bc across (cl-mpm:sim-bcs sim)
                      do
                        (when bc (when (equal (cl-mpm/mesh:get-node mesh (cl-mpm/bc:bc-index bc)) nil)
                                   (error "How on earth has bc ~A got a nil node rank ~D" bc rank))))
                ;;Cells
                (let ((prune-count 0))
                  (dotimes (i (array-total-size cells))
                    (let ((cell (row-major-aref cells i)))
                      (when cell
                        (when (not (in-computational-domain-buffer sim (cl-mpm/mesh::cell-centroid cell)
                                                                   0
                                                                   ;; buffer-size
                                                                   ))
                          (incf prune-count)
                          (setf (row-major-aref cells i) nil)))))
                  (format t "Rank ~D - Pruned ~D cells~%" rank prune-count))))))

        ))))
;; (defmethod %domain-decompose :after ((sim cl-mpm/mpi::mpm-sim-mpi-nodes-damage) domain-scaler)
;;   (let* ((rank (cl-mpi:mpi-comm-rank))
;;          (size (cl-mpi:mpi-comm-size)))
;;     (with-accessors ((mesh cl-mpm:sim-mesh)
;;                      )
;;         sim
;;       (let* ((domain-sizes (mapcar (lambda (x) (abs (reduce #'- x)))
;;                                     (cl-mpm/mpi::mpm-sim-mpi-domain-bounds sim)))
;;               (h (cl-mpm/mesh:mesh-resolution mesh))
;;               (min-domain-length (reduce #'max (remove 0d0 domain-sizes)))
;;               (buffer-size (ceiling min-domain-length h)))

;;          ))))

;; (defun kill-servers ()
;;     (dolist (server *open-servers*)
;;       (lfarm-admin:end-server (first server) (second server)))
;;   (setf *open-servers* nil))

;; (defparameter *global-dt* 1d0)
;; (lfarm:deftask uls (mp)
;;   (update-stress-mp mp *global-dt*)
;;   (setf (fill-pointer (cl-mpm/particle::mp-cached-nodes mp)) 0)
;;   mp)

;; (setf cl-mpi-extensions::*standard-encode-function*
;;       (lambda (x)
;;         (let ((res (flexi-streams:with-output-to-sequence (stream)
;;                      (cl-store:store x stream))))
;;           (static-vectors:make-static-vector (length res)
;;                       :element-type '(unsigned-byte 8)
;;                       :initial-contents res
;;                       ))))

;; (setf cl-mpi-extensions::*standard-decode-function*
;;       (lambda (x)
;;         (flexi-streams:with-input-from-sequence (stream x)
;;           (cl-store:restore stream))))


(defmethod cl-mpm::calculate-min-dt-mps ((sim cl-mpm/mpi::mpm-sim-mpi))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mass-scale cl-mpm::sim-mass-scale))
      sim
    (let ((inner-factor 
            ;;MPI not a fan of most-positive-double-float
            1d50
            ;most-positive-double-float
                        ))
      (iterate-over-nodes-serial
       mesh
       (lambda (node)
         (with-accessors ((node-active  cl-mpm/mesh:node-active)
                          (node-pos cl-mpm/mesh::node-position)
                          (pmod cl-mpm/mesh::node-pwave)
                          (mass cl-mpm/mesh::node-mass)
                          (svp-sum cl-mpm/mesh::node-svp-sum)
                          (vol cl-mpm/mesh::node-volume)
                          ) node
           (when (and node-active
                      ;(in-computational-domain sim node-pos)
                      (> vol 0d0)
                      (> pmod 0d0)
                      (> svp-sum 0d0))
             (let ((nf (/ mass (* vol (/ pmod svp-sum)))))
                 (when (< nf inner-factor)
                   (setf inner-factor nf)))))))
      (let ((rank (cl-mpi:mpi-comm-rank))
            (size (cl-mpi:mpi-comm-size)))
        (static-vectors:with-static-vector (source 1 :element-type 'double-float :initial-element inner-factor)
          (static-vectors:with-static-vector (dest 1 :element-type 'double-float :initial-element 0d0)
            (cl-mpi:mpi-allreduce source dest cl-mpi:+mpi-min+)
            (setf inner-factor (aref dest 0))
            (if (< inner-factor most-positive-double-float)
                (progn
                  ;; (format t "Rank ~D: dt - ~F~%" rank (* (sqrt mass-scale) (sqrt inner-factor) (cl-mpm/mesh:mesh-resolution mesh)))
                  ;; (format t "global : dt - ~F~%" (* (sqrt mass-scale) (sqrt inner-factor) (cl-mpm/mesh:mesh-resolution mesh)))
                  (* (sqrt mass-scale) (sqrt inner-factor) (cl-mpm/mesh:mesh-resolution mesh)))
                (cl-mpm::sim-dt sim))))))))



;; (defmacro define-constitutive (name lambda-list  &rest body)
;;   (let ((actual-name (intern (format nil "CONSTITUTIVE-~:@(~A~)" name)))
;;         (const-name (intern (format nil "MAKE-~:@(~A~)-DATA" name)))
;;         (iter 0)
;;         (new-ll (append (list 'constitutive-array) lambda-list))
;;         )
;;     `(progn
;;        (defun ,actual-name ,new-ll
;;          ,@(subst `(aref constitutive-array ,iter) 'give_me_zeros  body))
;;        (defun ,const-name ()
;;          (make-array ,iter)
;;          )
;;        )))

;; (defun give_me_zeros ())

;; (define-constitutive my-consitutive (strain)
;;   (let ((a (give_me_zeros))
;;         (b (give_me_zeros))
;;         )
;;     ;(magicl:@ strain a b)
;;     b
;;     ))



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


(defun save-damage-vtk (filename mps)
  (with-open-file (fs filename :direction :output :if-exists :supersede)
    (format fs "# vtk DataFile Version 2.0~%")
    (format fs "Lisp generated vtk file, SJVS~%")
    (format fs "ASCII~%")
    (format fs "DATASET UNSTRUCTURED_GRID~%")
    (format fs "POINTS ~d double~%" (length mps))
    (loop for mp across mps
          do (format fs "~E ~E ~E ~%"
                     (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) 'single-float)
                     (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 1 0) 'single-float)
                     (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 2 0) 'single-float)
                     ))
    (format fs "~%")

    ;; (cl-mpm/output::with-parameter-list fs mps
    ;;   ("mass" 'cl-mpm/particle:mp-mass)
    ;;   ("density" (lambda (mp) (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp))))
    ;;   )
    (let ((id 1))
      (declare (special id))
      (format fs "POINT_DATA ~d~%" (length mps))
      (cl-mpm/output::save-parameter "damage-y"
                                     (if (slot-exists-p mp 'cl-mpm/particle::damage-y-local)
                                         (cl-mpm/particle::mp-damage-y-local mp)
                                         0d0))
      )))
(defmethod cl-mpm/damage::update-localisation-lengths ((sim cl-mpm/mpi::mpm-sim-mpi-nodes-damage))
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (let ((damage-mps (cl-mpm/mpi::mpi-sync-damage-mps
                       sim
                       (cl-mpm/mpi::mpm-sim-mpi-halo-damage-size sim))))
      (lparallel:pdotimes (i (length damage-mps))
        (cl-mpm/damage::local-list-add-particle mesh (aref damage-mps i)))
      (call-next-method)
      (lparallel:pdotimes (i (length damage-mps))
        (cl-mpm/damage::local-list-remove-particle mesh (aref damage-mps i))))
    (values))
  )
(defmethod cl-mpm/damage::delocalise-damage ((sim cl-mpm/mpi::mpm-sim-mpi-nodes-damage))
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (let ((damage-mps (cl-mpm/mpi::mpi-sync-damage-mps
                       sim
                       (cl-mpm/mpi::mpm-sim-mpi-halo-damage-size sim))))
      (lparallel:pdotimes (i (length damage-mps))
        (cl-mpm/damage::local-list-add-particle mesh (aref damage-mps i)))
      (call-next-method)
      (lparallel:pdotimes (i (length damage-mps))
        (cl-mpm/damage::local-list-remove-particle mesh (aref damage-mps i))))
    (values)))

;;MPI stubs for setup
(defmethod cl-mpm/setup::%estimate-elastic-dt ((sim cl-mpm/mpi:mpm-sim-mpi))
  (cl-mpm/mpi::mpi-min
   (call-next-method)))
(defmethod cl-mpm/setup::estimate-critical-damping ((sim cl-mpm/mpi:mpm-sim-mpi))
  (cl-mpm/mpi::mpi-min
   (call-next-method)))




(defgeneric load-balance-metric (sim))

(defmethod load-balance-metric (sim)
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (set-mp-mpi-index sim)
    (lparallel:pcount-if (lambda (mp) (= (cl-mpm/particle::mp-mpi-index mp) rank))
                         (cl-mpm:sim-mps sim))))

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

(defun load-balance-setup (sim)
  (setf (mpm-sim-mpi-load-metric sim) (float (load-balance-metric sim) 0d0)))

(defun load-balance-dimension (sim dim &key (step-size 1d-2)
                                         (exchange nil))
  (let* ((rank (cl-mpi::mpi-comm-rank))
         (index-bottom-rank (mpi-rank-to-index sim rank))
         (dim-length (nth dim (mpm-sim-mpi-domain-count sim)))
         (dim-index (nth dim (mpi-rank-to-index sim rank)))
         (dim-size (abs (- (apply #'- (nth dim (mpm-sim-mpi-domain-bounds sim))))))
         (metric-array (make-array dim-length :element-type '(signed-byte 32) :initial-element 0))
         (size-array (make-array dim-length :element-type 'double-float :initial-element 0d0))
         (increment-length (- dim-length 1))
         (increment-array (make-array increment-length :element-type 'double-float :initial-element 0d0)))
    ;; (format t "Dims ~D length ~D~%" dim dim-length)
    (setf (nth dim index-bottom-rank) 0)
    (let ((bottom-rank (every #'zerop index-bottom-rank))
          (min-size (* 2 (mpm-sim-mpi-min-size sim))))
      ;; (when (= rank 0)
      ;;   (format t "Min size ~F~%" min-size)
      ;;   )
      (when (> dim-length 1)
        (setf (aref metric-array dim-index)
              (coerce (round (mpm-sim-mpi-load-metric sim)) 'integer)
              )
        ;; (format t "Metric ~E~%" (float (mpm-sim-mpi-load-metric sim) 0d0))
        (mpi-vector-integer-sum metric-array)
        ;; (format t "Dim index ~D~%" dim-index)
        ;; (when (= rank 0)
        ;;   (format t "Metric ~A~%" metric-array))
        (setf (aref size-array dim-index)
              (float dim-size 0d0))
        (mpi-vector-max size-array)
        (when bottom-rank
          (when (< dim-index increment-length)
            (setf (aref increment-array dim-index)
                  (float
                   (* step-size
                      (min (aref size-array dim-index)
                           (aref size-array (+ dim-index 1)))
                      (float
                       (if (or (> (aref metric-array dim-index) 0d0)
                               (> (aref metric-array (1+ dim-index)) 0d0))
                           (/ (- (aref metric-array dim-index)
                                 (aref metric-array (+ dim-index 1)))
                              (max (aref metric-array dim-index)
                                   (aref metric-array (+ dim-index 1))))
                           (signum
                            (-
                             (- (aref metric-array dim-index)
                                (aref metric-array (+ dim-index 1))))))
                       0d0))
                   0d0))
            (when (and
                   (<= (aref size-array dim-index) min-size)
                   (> (aref increment-array dim-index) 0d0))
              (setf (aref increment-array dim-index) 0d0)
              (format t "Rank ~D: Dim ~D hitting min size limit ~%" rank dim))
            (when (and
                   (<= (aref size-array (1+ dim-index)) min-size)
                   (< (aref increment-array dim-index) 0d0))
              (setf (aref increment-array dim-index) 0d0)
              (format t "Rank ~D: Dim ~D hitting min size limit ~%" (1+ rank) dim))
            ))
        (mpi-vector-sum increment-array)
        (let ((left-index (- dim-index 1)))
          (when (>= left-index 0)
            (incf (first (nth dim (mpm-sim-mpi-domain-bounds sim)))
                  (- (aref increment-array left-index)))
            ;; (format t "Rank ~D: Dim ~D moved lower by ~F ~%" (1+ rank) dim (- (aref increment-array left-index)))
            ))
        (let ((right-index dim-index))
          (when (< right-index increment-length)
            (incf (second (nth dim (mpm-sim-mpi-domain-bounds sim)))
                  (- (aref increment-array right-index)))
            ;; (format t "Rank ~D: Dim ~D moved upper by ~F ~%" (1+ rank) dim (- (aref increment-array right-index)))
            ))))
    (every (lambda (x) (= x 0d0)) increment-array)))

(defun load-balance-value (sim)
  (let ((rank (cl-mpi::mpi-comm-rank))
        (size (cl-mpi::mpi-comm-size)))
    (let* ((metric (float (mpm-sim-mpi-load-metric sim) 0d0))
           (sum-mps (mpi-sum metric))
           (max-mps (mpi-max metric))
           (balance nil))
      (when (> sum-mps 0)
        (setf balance
              (/
               max-mps
               (/ sum-mps
                  (cl-mpi:mpi-comm-size)))))
      ;; (when t;(> min-mps 0)
      ;;   ;; (setf balance (mpi-max (/ max-mps min-mps)))
      ;;   (when (= rank 0)
      ;;     (format t "Occupancy ratio : ~F%~%" (* 100d0 balance))))
      balance)
    ;; (let ((min-mps (mpi-min (float (load-balance-metric sim) 0d0)))
    ;;       (max-mps (mpi-max (float (load-balance-metric sim) 0d0)))
    ;;       (sum-mps (mpi-sum (float (load-balance-metric sim) 0d0)))
    ;;       (balance nil))
    ;;   ;; (when (= rank 0)
    ;;   ;;   (format t "Min : ~F%~%" min-mps)
    ;;   ;;   (format t "Max : ~F%~%" max-mps)
    ;;   ;;   )
    ;;   (when (> min-mps 0)
    ;;     (setf balance (mpi-max (/ max-mps min-mps)))
    ;;     (when (= rank 0)
    ;;       (format t "Occupancy ratio : ~F%~%" (* 100d0 balance))))
    ;;   ;; (format t "Rank ~D: Domain bounds ~A~%" rank (mpm-sim-mpi-domain-bounds sim))
    ;;   balance)
    ))


(defun load-balance (sim &key (substeps 10)
                           (step-size 1d-1)
                           (exchange-mps t)
                           (dims (list :x :y :z))
                           )
  (let ((rank (cl-mpi::mpi-comm-rank))
        (stagnent nil))
    (load-balance-setup sim)
    (loop for j from 0 to substeps
          while (not stagnent)
          do
             (progn
               (setf stagnent t)
               (loop for dim in dims ;(i (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
                     do
                        (let ((i (position dim (list :x :y :z))))
                          (when t;(< i (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
                            (setf
                             stagnent
                             (and
                              stagnent
                              (load-balance-dimension sim i :step-size step-size)))
                            (when exchange-mps
                              (set-mp-mpi-index sim)
                              (exchange-mps sim 0d0)
                              (set-mp-mpi-index sim)
                              (clear-ghost-mps sim))
                            (load-balance-setup sim))))))
    (let* ((metric (float (mpm-sim-mpi-load-metric sim) 0d0))
           (min-mps (mpi-min metric))
           (max-mps (mpi-max metric))
           (balance nil))
      (setf balance
            (/
             (mpi-max metric)
             (/ (mpi-sum metric)
                (cl-mpi:mpi-comm-size))))
      (when t;(> min-mps 0)
        ;; (setf balance (mpi-max (/ max-mps min-mps)))
        (when (= rank 0)
          (format t "Occupancy ratio : ~F%~%" (* 100d0 balance))))
      (format t "Rank ~D: Domain bounds ~A~%" rank (mpm-sim-mpi-domain-bounds sim))
      (values balance stagnent))))

(defun load-balance-algo (sim &key (substeps 10)
                                (max-iters 50)
                                (min-bounds 1.1d0)
                                (max-bounds 1.5d0)
                                (step-size 1d-1)
                                (dims (list :x :y :z))
                                )
  (when (typep sim 'cl-mpm/mpi:mpm-sim-mpi)
    (load-balance-setup sim)
    (let ((balance (load-balance-value sim))
          (rank (cl-mpi::mpi-comm-rank)))
      (when (= rank 0)
        (format t "Check balance~%"))
      (when (and balance (> balance max-bounds))
        (let ((balance nil)
              (stag nil))
          (loop repeat max-iters
                while (and (if balance (> balance min-bounds) t)
                           (not stag))
                do
                   (multiple-value-bind (balance stagnent) (cl-mpm/mpi::load-balance sim
                                                                                     :exchange-mps t
                                                                                     :step-size step-size
                                                                                     :substeps substeps
                                                                                     :dims dims)
                     (setf balance balance
                           stag stagnent))))
        (domain-decompose sim)))))
