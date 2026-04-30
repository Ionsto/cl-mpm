(in-package :cl-mpm/mpi)
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
