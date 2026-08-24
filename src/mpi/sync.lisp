(in-package :cl-mpm/mpi)


;; (defparameter *exception-comm* (cl-mpi::mpi-comm-dup))

(defun check-exception ()
  (static-vectors:with-static-vector (source 1 :element-type 'double-float :initial-element -1d0)
    (static-vectors:with-static-vector (dest 1 :element-type 'double-float :initial-element 0d0)
      (cl-mpi:mpi-allreduce source dest cl-mpi:+mpi-max+)
      (unless (= (aref dest 0) -1d0)
        (format t "MPI error propogated ~D~%" (cl-mpi:mpi-comm-rank))
        (let ((cl-mpi-extensions::*standard-encode-function* #'cl-store-encoder)
              (cl-mpi-extensions::*standard-decode-function* #'cl-store-decoder))
          (let ((er
                  (cl-mpi-extensions:mpi-broadcast-anything
                   (round (aref dest 0)))))
            (error er)))))))

(defun throw-mpi-error (err)
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (format t "Throwing error over MPI from ~D~%" (cl-mpi:mpi-comm-rank))
    (format t "~A~%" err)
    (static-vectors:with-static-vector (source 1 :element-type 'double-float :initial-element (float (cl-mpi:mpi-comm-rank) 0d0))
      (static-vectors:with-static-vector (dest 1 :element-type 'double-float :initial-element 0d0)
        (cl-mpi:mpi-allreduce source dest cl-mpi:+mpi-max+)
        (format t "All threads reached checkpoint~%")
        (let ((cl-mpi-extensions::*standard-encode-function* #'cl-store-encoder)
              (cl-mpi-extensions::*standard-decode-function* #'cl-store-decoder))
          (if (= rank (round (aref dest 0)))
              (progn
                (cl-mpi-extensions:mpi-broadcast-anything
                 (cl-mpi:mpi-comm-rank)
                 :object err)
                (format t "Broadcast error, now rethrow~%")
                (error err))
              (let ((er (cl-mpi-extensions:mpi-broadcast-anything
                         (cl-mpi:mpi-comm-rank)
                         )))
                (format t "Multiple errors thrown, ~D defers to ~D~%" rank (round (aref dest 0)))
                (error er))))))))

;; (defmacro with-mpi-errors (&body body)
;;   `;; (handler-bind
;;    ;;     ((error (lambda (e)
;;    ;;               (format t "Throwing err~%")
;;    ;;               (throw-mpi-error e))))
;;    ;;   (progn
;;    ;;     ,@body
;;    ;;     (check-exception)))
;;   ;; (let ((block-name (gensym)))
;;   ;;   `(block ,block-name
;;   ;;      (handler-case
;;   ;;          (error (lambda (e)
;;   ;;                   (format t "Throwing err~%")
;;   ;;                   (throw-mpi-error e)))
;;   ;;          (progn
;;   ;;            ,@body
;;   ;;            (check-exception)))))
;;   )

(defmacro with-mpi-errors (&body body)
  `(progn
     (check-exception)
     (handler-case
         (progn
           ,@body)
       (error (e)
         (format t "Throwing err~%")
         (throw-mpi-error e)))
     (check-exception)))


(defun exchange-node-like (sim
                           serialise
                           deserialise
                           func)
  (declare (function func))
  (let* ((rank (cl-mpi:mpi-comm-rank)))
    (with-accessors ((mesh cl-mpm:sim-mesh))
        sim
      (let* ((index (mpi-rank-to-index sim rank))
             (bounds-list (mpm-sim-mpi-domain-bounds sim))
             (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim))))
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
                                             right-neighbor :tag 2))))
                                        ((and
                                          (= left-neighbor -1)
                                          (not (= right-neighbor -1)))
                                         (let ((r (active-filter right-filter)))
                                           (cl-mpi-extensions:mpi-waitall-anything
                                            (cl-mpi-extensions:mpi-irecv-anything right-neighbor :tag 1)
                                            (cl-mpi-extensions:mpi-isend-anything
                                             r
                                             right-neighbor :tag 2))))
                                        ((and
                                          (not (= left-neighbor -1))
                                          (= right-neighbor -1))
                                         (let ((l (active-filter left-filter))
                                               )
                                           (cl-mpi-extensions:mpi-waitall-anything
                                            (cl-mpi-extensions:mpi-irecv-anything left-neighbor :tag 2)
                                            (cl-mpi-extensions:mpi-isend-anything
                                             l
                                             left-neighbor :tag 1))))
                                        (t nil))))
                               (loop for packet in recv
                                     do (destructuring-bind (rank tag object) packet
                                          (when object
                                            (funcall func object))))))))))))))))

(defun exchange-nodes (sim func)
  (declare (function func))
  (exchange-node-like
   sim
   #'serialise-node
   #'deserialise-node
   func))






(declaim (notinline mpi-sync-momentum))
(defun mpi-sync-momentum (sim)
  ;; (format t "Sync momentum ~%")
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
      (exchange-nodes
       sim
       (lambda (node-list)
         (cl-mpm/utils::bpdotimes (i (length node-list))
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

(declaim (notinline mpi-sync-displacement))
(defun mpi-sync-displacement (sim)
  ;; (format t "Sync momentum ~%")
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (cl-mpm::iterate-over-nodes
     (cl-mpm:sim-mesh sim)
     (lambda (n)
       (when (not (in-computational-domain sim (cl-mpm/mesh::node-position n)))
         (cl-mpm::fast-zero (cl-mpm/mesh::node-acceleration n))
         (cl-mpm::fast-zero (cl-mpm/mesh::node-displacment n)))))
    (exchange-nodes
     sim
     (lambda (node-list)
       (cl-mpm/utils::bpdotimes (i (length node-list))
         (let* ((mpi-node (aref node-list i))
                (index (mpi-object-node-index mpi-node))
                (node (cl-mpm/mesh:get-node mesh index)))
           (if node
               (progn
                 (with-accessors ((disp cl-mpm/mesh::node-displacment)
                                  (acc cl-mpm/mesh::node-acceleration)
                                  )
                     node
                   (cl-mpm/fastmaths::fast-.+ disp (mpi-object-node-displacement mpi-node) disp)
                   (cl-mpm/fastmaths::fast-.+ acc (mpi-object-node-acc mpi-node) acc)))
               (error "MPI exchange touched invalid node?"))))))))

(declaim (notinline mpi-sync-mass))
(defun mpi-sync-mass (sim)
  ;; (format t "Sync momentum ~%")
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (exchange-node-like
     sim
     #'serialise-node-mass
     #'deserialise-node-mass
     (lambda (node-list)
       (cl-mpm/utils::bpdotimes (i (length node-list))
         (let* ((mpi-node (aref node-list i))
                (index (mpi-object-node-mass-index mpi-node))
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
                   (incf mass (the double-float (mpi-object-node-mass-mass mpi-node)))
                   ;; (incf svp (the double-float  (mpi-object-node-mass-svp mpi-node)))
                   ;; (incf vol (the double-float  (mpi-object-node-mass-vol mpi-node)))
                   ;; (incf pmod (the double-float (mpi-object-node-mass-pmod mpi-node)))
                   ))
               (error "MPI exchange touched invalid node?"))))))))

(defun mpi-sync-j-inc (sim)
  ;; (format t "Sync momentum ~%")
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
      (exchange-nodes
       sim
       (lambda (node-list)
         (cl-mpm/utils::bpdotimes (i (length node-list))
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
       (cl-mpm/utils::bpdotimes (i (length node-list))
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
                   (cl-mpm/fastmaths::fast-.+ force-buoyancy (mpi-object-node-force-buoyancy mpi-node) force-buoyancy)))
               (error "MPI force exchange touched invalid node ~A" index))))))))



(defun ser-part (mps &optional parts)
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


(defun clear-ghost-mps (sim)
  (let ((rank (cl-mpi:mpi-comm-rank)))
     (cl-mpm::remove-mps-func
      sim
      (lambda (mp)
        (not (= rank (cl-mpm/particle::mp-mpi-index mp)))))))


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


(declaim (notinline exchange-mps))
(defun exchange-mps (sim &optional halo-depth)
  (let* ((rank (cl-mpi:mpi-comm-rank))
         (size (cl-mpi:mpi-comm-size)))
    (with-accessors ((mps cl-mpm:sim-mps)
                     (mesh cl-mpm:sim-mesh))
        sim
      (cl-mpm:iterate-over-mps
       (cl-mpm:sim-mps sim)
       (lambda (mp)
         (setf (fill-pointer (cl-mpm/particle::mp-cached-nodes mp)) 0)))
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
                                        (when object
                                          (when (> (length object) 0)
                                            (cl-mpm/utils::bpdotimes
                                             (i (length object))
                                             (let ((mp (aref object i)))
                                               (setf (fill-pointer (cl-mpm/particle::mp-cached-nodes mp)) 0)
                                               (when (slot-exists-p mp 'cl-mpm/particle::damage-position)
                                                 (setf (cl-mpm/particle::mp-damage-position mp) nil))
                                               (cl-mpm::sim-add-mp sim mp)
                                               ;; (vector-push-extend mp mps)
                                               ))))))))))))))))
  (cl-mpi:mpi-barrier))
;  )

(defvar *mutex-code* (cl-store:register-code 110 'sb-thread:mutex))
(cl-store:defstore-cl-store (obj sb-thread:mutex stream)
   (cl-store:output-type-code *mutex-code* stream))
(cl-store:defrestore-cl-store (sb-thread:mutex stream)
   (sb-thread:make-mutex))


(defvar *particle-code* (cl-store:register-code 100 'cl-mpm/particle::particle))
(cl-store:defstore-cl-store (obj cl-mpm/particle::particle stream)
  (cl-store:output-type-code *particle-code* stream)
  (loop for slot in (cl-store:serializable-slots obj)
        do
           (cond
             ((eq 'cl-mpm/particle::cached-nodes (sb-mop:slot-definition-name slot))
              nil)
             (t
              (cl-store:store-object (sb-mop:slot-value-using-class 'cl-mpm/mesh::node obj slot) stream)))))

(cl-store:defrestore-cl-store (cl-mpm/particle::particle stream)
  (let ((obj (make-instance 'cl-mpm/particle::particle)))
    (loop for slot in (cl-store:serializable-slots obj)
          do
             (cond
               ((eq 'cl-mpm/particle::cached-nodes (sb-mop:slot-definition-name slot))
                (setf (cl-mpm/particle::mp-cached-nodes obj) (make-array 8 :fill-pointer 0 :element-type 'cl-mpm/particle::node-cache :initial-element (cl-mpm/particle::make-empty-node-cache)))
                )
               (t
                (setf (sb-mop:slot-value-using-class 'cl-mpm/mesh::node obj slot) (cl-store:restore-object stream)))))
    obj))

;; (defvar *node-cache-code* (cl-store:register-code 113 'cl-mpm/particle::node-cache))
;; (cl-store:defstore-cl-store (obj cl-mpm/particle::node-cache stream)
;;   (cl-store:output-type-code *node-cache-code* stream)
;;   )
;; (cl-store:defrestore-cl-store (cl-mpm/particle::node-cache stream)
;;   (let ((obj (cl-mpm/particle::make-node-cache )))
;;     obj)
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
             ((eq 'cl-mpm/mesh::agg-interior-cell (sb-mop:slot-definition-name slot))
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
