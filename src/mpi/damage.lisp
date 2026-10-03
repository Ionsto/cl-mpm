(in-package :cl-mpm/mpi)
(declaim #.cl-mpm/settings:*optimise-setting*)

;; (defparameter *damage-mp-send-cache* (make-array 0 :element-type 'cl-mpm::particle :adjustable t :fill-pointer 0))
(declaim (notinline mpi-sync-damage-mps))
(defun mpi-sync-damage-mps (sim &key
                                  (halo-depth nil)
                                  (update-mps nil))
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

        (setf (fill-pointer damage-mps) 0)

        (loop for i from 0 below nd
              do
                 (let ((id-delta (list 0 0 0)))
                   (setf (nth i id-delta) 1)
                   (let ((left-neighbor (mpi-index-to-rank sim (mapcar #'- index id-delta)))
                         (right-neighbor (mpi-index-to-rank sim (mapcar #'+ index id-delta))))

                     (declare (double-float halo-depth))
                     ;; (format t "Transfer ~D ~D~%" i rank)
                     (destructuring-bind (bl bu) (nth i bounds-list)
                       (declare (double-float bl bu))
                       (labels
                           ((halo-filter (test)
                              (let ((res
                                      (lparallel:premove-if-not
                                       ;remove-if-not
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
                         (let* ((cl-mpi-extensions::*standard-encode-function* #'serialise-damage-mp)
                                (cl-mpi-extensions::*standard-decode-function* #'deserialise-damage-mp)
                                (recv
                                  (cond
                                    ((and (not (= left-neighbor -1))
                                          (not (= right-neighbor -1)))
                                     ;; (format t "Build filters~%")
                                     (let ((l (left-filter))
                                           (r (right-filter)))
                                       ;; (format t "Built~%")
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
                                     ;; (format t "Build filters~%")
                                     (let ((r (right-filter)))
                                       ;; (format t "Built - ~D~%" rank)
                                       (cl-mpi-extensions:mpi-waitall-anything
                                        (cl-mpi-extensions:mpi-irecv-anything right-neighbor :tag 1)
                                        (cl-mpi-extensions:mpi-isend-anything
                                         r
                                         right-neighbor :tag 2))))
                                    ((and
                                      (not (= left-neighbor -1))
                                      (= right-neighbor -1))
                                     ;; (format t "Build filters~%")
                                     (let ((l (left-filter)))
                                       ;; (format t "Built - ~D~%" rank)
                                       (cl-mpi-extensions:mpi-waitall-anything
                                        (cl-mpi-extensions:mpi-irecv-anything left-neighbor :tag 2)
                                        (cl-mpi-extensions:mpi-isend-anything
                                         l
                                         left-neighbor :tag 1))))
                                    (t nil))))
                           ;; (format t "Received~%")
                           (let ((current-list (sim-mpi-damage-mps-list sim)))
                             (declare ((vector t *) current-list))
                             (when (and (not update-mps)
                                        (not (= (length current-list) 0)))
                               (format t "We've got MPS but not updating?"))
                             (loop for packet in recv
                                   do
                                      (destructuring-bind (rank tag object) packet
                                        (when object
                                          (loop for mp across object
                                                do (progn
                                                     (let ((dummy-mp nil))
                                                       (if update-mps
                                                           (let ((index (position (mpi-object-damage-mp-unique-id mp) current-list :key (lambda (mp-o) (cl-mpm/particle::mp-unique-index mp-o)))))
                                                             (unless index
                                                               (format t "~A~%" current-list)
                                                               (error "MP with unique id ~A not found in current list" (mpi-object-damage-mp-unique-id mp)))
                                                             (setf dummy-mp (aref current-list index)))
                                                           (setf dummy-mp (allocate-instance (find-class 'cl-mpm/particle::particle-damage))))
                                                       (setf (slot-value dummy-mp 'cl-mpm/particle::damage) (mpi-object-damage-mp-damage mp)
                                                             (slot-value dummy-mp 'cl-mpm/particle::volume) (mpi-object-damage-mp-volume mp)
                                                             (slot-value dummy-mp 'cl-mpm/particle::volume-n) (mpi-object-damage-mp-volume-n mp)
                                                             (slot-value dummy-mp 'cl-mpm/particle::unique-index) (mpi-object-damage-mp-unique-id mp)
                                                             (slot-value dummy-mp 'cl-mpm/particle::position) (mpi-object-damage-mp-position mp)
                                                             (slot-value dummy-mp 'cl-mpm/particle::position-trial) (mpi-object-damage-mp-position-trial mp)
                                                             (slot-value dummy-mp 'cl-mpm/particle::damage-y-local) (mpi-object-damage-mp-y mp)
                                                             (slot-value dummy-mp 'cl-mpm/particle::true-local-length) (mpi-object-damage-mp-local-length mp)
                                                             (slot-value dummy-mp 'cl-mpm/particle::average-damage) (mpi-object-damage-mp-average-damage mp))
                                                       ;; (unless update-mps
                                                       ;;   (slot-value dummy-mp 'cl-mpm/particle::damage-position) nil)
                                                       ;; (when (and (slot-boundp dummy-mp 'cl-mpm/particle::damage-position)
                                                       ;;            (not update-mps))
                                                       ;;   (cl-mpm/damage::local-list-remove-particle mesh dummy-mp))
                                                       (unless update-mps
                                                         (setf (slot-value dummy-mp 'cl-mpm/particle::damage-position) nil))
                                                       ;; (unless update-mps)
                                                       (vector-push-extend
                                                        dummy-mp
                                                        damage-mps))))))))))))))
        damage-mps))))




(defun partial-rebuild-mp-local-list (sim)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)
                   (dhalo cl-mpm/mpi::mpm-sim-mpi-halo-damage-size))
      sim
    (cl-mpm::iterate-over-mps
     mps
     (lambda (mp)
       (when (in-computational-domain-buffer sim (cl-mpm/particle::mp-position mp)
                                             (/ dhalo (cl-mpm/mesh::mesh-resolution mesh)))
         (cl-mpm/damage::build-mp-local-list mesh mp))))))

(defmethod cl-mpm/damage::update-localisation-lengths ((sim cl-mpm/mpi::mpm-sim-mpi-damage))
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (let ((damage-mps (cl-mpm/mpi::mpi-sync-damage-mps
                       sim
                       :halo-depth (cl-mpm/mpi::mpm-sim-mpi-halo-damage-size sim)
                       :update-mps t)))
      ;; (partial-rebuild-mp-local-list sim)
      ;; (format t "Call next~%")
      (call-next-method)
      ;; (partial-rebuild-mp-local-list sim)
      )
    (values))
  )
(defmethod cl-mpm/damage::delocalise-damage ((sim cl-mpm/mpi::mpm-sim-mpi-damage))
  ;; (format t "Sync delocalise length~%")
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (let ((damage-mps (cl-mpm/mpi::mpi-sync-damage-mps
                       sim
                       :halo-depth (cl-mpm/mpi::mpm-sim-mpi-halo-damage-size sim)
                       :update-mps t)))
      (call-next-method)
      )
    ;; (cl-mpm/dynamic-relaxation::save-vtks-dr-step sim "./" 0 0 0)
    ;; (cl-mpi::mpi-waitall)
    ;; (break)
    (values)))

;; (defun full-reset (sim)
;;   (cl-mpm::iterate-over-mps
;;    (cl-mpm::sim-mps sim)
;;    (lambda (mp)
;;      (setf (cl-mpm/particle::mp-damage-position mp) nil)))
;;   (cl-mpm::iterate-over-nodes
;;    (cl-mpm::sim-mesh sim)
;;    (lambda (node)
;;      (setf (cl-mpm/mesh::node-local-list node)
;;            (make-array 0 :element-type 'cl-mpm::particle :adjustable t :fill-pointer 0)))))

(defparameter *lliters* 0)
(defmethod cl-mpm/damage::update-delocalisation-list ((sim cl-mpm/mpi::mpm-sim-mpi-damage))
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mps cl-mpm:sim-mps))
      sim
      (with-accessors ((nodes cl-mpm/mesh:mesh-nodes)
                       (h cl-mpm/mesh:mesh-resolution))
          mesh
        (cl-mpm/utils::bpdotimes
         (i (length (sim-mpi-damage-mps-list sim)))
         (progn
           (unless (cl-mpm/particle::mp-damage-position (aref (sim-mpi-damage-mps-list sim) i))
             (format t "Damage position unknown?~%"))
           (cl-mpm/damage::local-list-remove-particle mesh (aref (sim-mpi-damage-mps-list sim) i))))
        (setf (fill-pointer (sim-mpi-damage-mps-list sim)) 0)
        ;; (full-reset sim)

        (cl-mpm:iterate-over-mps
         mps
         (lambda (mp)
           (when (typep mp 'cl-mpm/particle:particle-damage)
             (if (eq (cl-mpm/particle::mp-damage-position mp) nil)
                 (cl-mpm/damage::local-list-add-particle mesh mp)
                 (let* ((delta (cl-mpm/fastmaths::diff-mag (cl-mpm/particle:mp-position mp)
                                                           (cl-mpm/particle::mp-damage-position mp))))
                   (declare (double-float delta h))
                   (when (> delta (/ h 16d0))
                     (when (not (equal
                                 (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle:mp-position mp))
                                 (cl-mpm/mesh:position-to-index mesh (cl-mpm/particle::mp-damage-position mp))))
                       (cl-mpm/damage::local-list-remove-particle mesh mp)
                       (cl-mpm/damage::local-list-add-particle mesh mp))))))))


        (let ((damage-mps (cl-mpm/mpi::mpi-sync-damage-mps
                           sim
                           :halo-depth (cl-mpm/mpi::mpm-sim-mpi-halo-damage-size sim)
                           :update-mps nil)))
          (setf (fill-pointer (sim-mpi-damage-mps-list sim)) 0)
          (loop for mp across damage-mps
                do (vector-push-extend mp (sim-mpi-damage-mps-list sim)))
          ;; (setf (sim-mpi-damage-mps-list sim) (copy-seq damage-mps))
          (cl-mpm/utils::bpdotimes
           (i (length (sim-mpi-damage-mps-list sim)))
           (cl-mpm/damage::local-list-add-particle mesh (aref (sim-mpi-damage-mps-list sim) i))))
        (cl-mpm/damage::setup-mp-local-list sim)
        ;; (cl-mpm/dynamic-relaxation::save-vtks-dr-step sim "./" 0 0 *lliters*)
        ;; (incf *lliters*)
        ;; (cl-mpi::mpi-waitall)
        ;; (break)
        )))

(in-package :cl-mpm/output)
(defun save-mpi-damage-vtk (filename sim)
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (let ((mps (cl-mpm/mpi::sim-mpi-damage-mps-list sim)))
      (with-open-file (fs filename :direction :output :if-exists :supersede)
        (format fs "# vtk DataFile Version 2.0~%")
        (format fs "Lisp generated vtk file, SJVS~%")
        ;; (format fs "ASCII~%")
        (format fs "BINARY~%")
        (format fs "DATASET UNSTRUCTURED_GRID~%")
        (format fs "POINTS ~d double~%" (length mps)))
      (with-open-file (fs filename :direction :output :if-exists :append)
        (with-open-file (fs-bin filename :direction :output :if-exists :append :element-type '(unsigned-byte 8))
          (force-output fs)
          (loop for mp across mps
                do
                   (let ((pos (cl-mpm/particle::mp-position-trial mp)))
                     (write-binary-float (cl-mpm/utils:varef pos 0) fs-bin)
                     (write-binary-float (cl-mpm/utils:varef pos 1) fs-bin)
                     (write-binary-float (cl-mpm/utils:varef pos 2) fs-bin)
                     ;; (format fs "~E ~E ~E ~%"
                     ;;         (coerce (cl-mpm/utils:varef pos 0) 'single-float)
                     ;;         (coerce (cl-mpm/utils:varef pos 1) 'single-float)
                     ;;         (coerce (cl-mpm/utils:varef pos 2) 'single-float))
                     ))
          (force-output fs-bin)
          (format fs "~%")
          (let ((id 1)
                (nd (cl-mpm/mesh:mesh-nd mesh)))
            (declare (special id))
            (format fs "POINT_DATA ~d~%" (length mps))
            (let ((output-list (list (list :SCALAR "damage" #'cl-mpm/particle::mp-damage)
                                     (list :SCALAR "volume" #'cl-mpm/particle::mp-volume)
                                     (list :SCALAR "volume-n" #'cl-mpm/particle::mp-volume-n)
                                     (list :SCALAR "y" #'cl-mpm/particle::mp-damage-y-local)
                                     (list :SCALAR "local-length" #'cl-mpm/particle::mp-true-local-length)
                                     (list :SCALAR "average-damage" #'cl-mpm/particle::mp-av-damage))))
              (dolist (f output-list)
                (destructuring-bind (type name accessor) f
                  (case type
                    (:BOOL
                     (cl-mpm/output::save-parameter name (if (funcall accessor mp) 1d0 0d0)))
                    (:SCALAR
                     (cl-mpm/output::save-parameter name (funcall accessor mp)))
                    (:VECTOR
                     (cl-mpm/output::save-parameter (format nil "~A_mag" name) (cl-mpm/fastmaths::mag (funcall accessor mp)))
                     (cl-mpm/output::save-parameter (format nil "~A_x" name) (varef (funcall accessor mp) 0))
                     (cl-mpm/output::save-parameter (format nil "~A_y" name) (varef (funcall accessor mp) 1))
                     (when (= nd 3)
                       (cl-mpm/output::save-parameter (format nil "~A_z" name) (varef (funcall accessor mp) 2))))
                    (:VOIGT
                     (cl-mpm/output::save-parameter (format nil "~A_xx" name) (varef (funcall accessor mp) 0))
                     (cl-mpm/output::save-parameter (format nil "~A_yy" name) (varef (funcall accessor mp) 1))
                     (when (= nd 3)
                       (cl-mpm/output::save-parameter (format nil "~A_zz" name) (varef (funcall accessor mp) 2))
                       (cl-mpm/output::save-parameter (format nil "~A_yz" name) (varef (funcall accessor mp) 3))
                       (cl-mpm/output::save-parameter (format nil "~A_xz" name) (varef (funcall accessor mp) 4)))
                     (cl-mpm/output::save-parameter (format nil "~A_xy" name) (varef (funcall accessor mp) 5)))
                    (:MATRIX
                     (cl-mpm/output::save-parameter (format nil "~A_xx" name) (mtref (funcall accessor mp) 0 0))
                     (cl-mpm/output::save-parameter (format nil "~A_yy" name) (mtref (funcall accessor mp) 1 1))
                     (cl-mpm/output::save-parameter (format nil "~A_xy" name) (mtref (funcall accessor mp) 0 1))
                     (cl-mpm/output::save-parameter (format nil "~A_yx" name) (mtref (funcall accessor mp) 1 0)))))))))))))
