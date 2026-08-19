(in-package :cl-mpm/mpi)

;; (defparameter *damage-mp-send-cache* (make-array 0 :element-type 'cl-mpm::particle :adjustable t :fill-pointer 0))
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

        (setf (fill-pointer damage-mps) 0)
        ;; (format t "Start transfer~%")
        ;; (setf damage-mps (make-array 0 :element-type t :adjustable t :fill-pointer 0))
        ;; (if damage-mps
        ;;   (setf (fill-pointer damage-mps) 0)
        ;;   (setf damage-mps (make-array 0 :element-type t :adjustable t :fill-pointer 0)))
        ;; (pprint damage-mps)
        ;; (setf (fill-pointer *damage-mp-send-cache*) 0)
        ;; (format t "Begin transfer~%")
        (cl-mpi::mpi-barrier)
        ;; (format t "Synced~%")
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
                                       ;remove-if-not
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
                           (loop for packet in recv
                                 do
                                    (destructuring-bind (rank tag object) packet
                                      (when object
                                        (loop for mp across object
                                              do (progn
                                                   (vector-push-extend
                                                    (make-instance 'cl-mpm/particle::particle-damage
                                                                   :damage (mpi-object-damage-mp-damage mp)
                                                                   :volume (mpi-object-damage-mp-volume mp)
                                                                   :position (mpi-object-damage-mp-position mp)
                                                                   :damage-y (mpi-object-damage-mp-y mp)
                                                                   :local-length (mpi-object-damage-mp-local-length mp))
                                                    damage-mps))))))))))))
        damage-mps))))

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
(defmethod cl-mpm/damage::update-localisation-lengths ((sim cl-mpm/mpi::mpm-sim-mpi-damage))
  ;; (format t "Sync update-localisation length~%")
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (let ((damage-mps (cl-mpm/mpi::mpi-sync-damage-mps
                       sim
                       (cl-mpm/mpi::mpm-sim-mpi-halo-damage-size sim))))
      ;; (format t "Update local list~%")
      (cl-mpm/utils::bpdotimes (i (length damage-mps))
        (cl-mpm/damage::local-list-add-particle mesh (aref damage-mps i)))
      ;; (format t "Call next~%")
      (call-next-method)
      (cl-mpm/utils::bpdotimes (i (length damage-mps))
        (cl-mpm/damage::local-list-remove-particle mesh (aref damage-mps i)))
      )
    (values))
  )
(defmethod cl-mpm/damage::delocalise-damage ((sim cl-mpm/mpi::mpm-sim-mpi-damage))
  ;; (format t "Sync delocalise length~%")
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (let ((damage-mps (cl-mpm/mpi::mpi-sync-damage-mps
                       sim
                       (cl-mpm/mpi::mpm-sim-mpi-halo-damage-size sim))))
      (cl-mpm/utils::bpdotimes (i (length damage-mps))
        (cl-mpm/damage::local-list-add-particle mesh (aref damage-mps i)))
      ;; (format t "Call next~%")
      (call-next-method)
      (cl-mpm/utils::bpdotimes (i (length damage-mps))
        (cl-mpm/damage::local-list-remove-particle mesh (aref damage-mps i)))
      )
    (values)))
