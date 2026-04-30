(in-package :cl-mpm/mpi)

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
 node-mass
 ((index index cl-mpm/mesh::node-index)
  (float mass cl-mpm/mesh::node-mass)
  (float pmod cl-mpm/mesh::node-pwave)
  (float svp cl-mpm/mesh::node-svp-sum)
  (float vol cl-mpm/mesh::node-volume)))

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
  (vector displacement cl-mpm/mesh::node-displacment)
  ))
