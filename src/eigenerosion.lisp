(defpackage :cl-mpm/eigenerosion
  (:use :cl
        :cl-mpm
        :cl-mpm/utils)
  (:export
    #:update-fracture
    )
  )
(in-package :cl-mpm/eigenerosion)
(declaim (optimize (debug 0) (safety 0) (speed 3)))
(declaim (ftype (function (cl-mpm/particle:particle) (values)) calculate-strain-energy-mp))
(defun calculate-strain-energy-mp (mp)
  (declare (cl-mpm/particle:particle-damage mp))
  (with-accessors ((stress cl-mpm/particle:mp-stress)
                   (strain cl-mpm/particle:mp-strain)
                   (damage cl-mpm/particle:mp-damage)
                   (volume cl-mpm/particle:mp-volume)
                   (strain-energy-density cl-mpm/particle::mp-strain-energy-density)
                   ) mp
    (progn
      (multiple-value-bind (l v) (magicl:eig (voight-to-matrix strain))
        ;; (loop for i from 0 to 1
        ;;       do (let ((sii (nth i l)))
        ;;            (when (< sii 0) (setf (nth i l) 0))))
        (map-into l (lambda (sii) (max sii 0)) l)

        (let* ((driving-strain (matrix-to-voight (magicl:@ v (magicl:from-diag l :type 'double-float) (magicl:transpose v))))
               (ss (magicl:.* stress driving-strain))
               ;; (energy (* 0.5d0 (+ (magicl::sum ss) (magicl:tref ss 2 0))))
               (energy (* 0.5d0 (+ (magicl:tref ss 0 0) (magicl:tref ss 1 0) (magicl:tref ss 2 0) (magicl:tref ss 2 0)))))
          ;; (print driving-strain)
          (setf strain-energy-density energy)
          ))))
  (values))
(defun calculate-strain-energy (mps)
  ;; (lparallel:pmapcar #'calculate-strain-energy-mp mps)
  (lparallel:pdotimes (i (length mps)) 
    (calculate-strain-energy-mp (aref mps i)))
  )
(defun weight-func (distance length-parameter)
  ;(- 1 (/ (apply #'+ (mapcar #'* distance distance)) (*  length-parameter length-parameter)))
  (apply #'* (mapcar (lambda (x) (cl-mpm::shape-bspline x length-parameter)) distance))
  )
(defun reset-strain (mesh)
  (declare (cl-mpm/mesh::mesh mesh))
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
  (dotimes (i (array-total-size nodes))
    (let ((node (row-major-aref nodes i)))
      (when (> (cl-mpm/mesh::node-strain-energy-density node) 0)
        (setf (cl-mpm/mesh::node-strain-energy-density node) 0))))))
(defun scatter-local-energy (sim mps length-parameter order)
  (let ((mesh (cl-mpm:sim-mesh sim)))
    (lparallel:pdotimes (i (length mps))
        (cl-mpm::iterate-over-neighbours-shape-bspline mesh (aref mps i)
                                                       (lambda (mesh mp node weight dsvp)
                                                   (with-accessors ((strain cl-mpm/mesh::node-strain-energy-density)
                                                                   (lock cl-mpm/mesh::node-lock) 
                                                                    )
                                                       node
                                                     (let* ((mp-strain (cl-mpm/particle::mp-strain-energy-density mp))
                                                          (mass (cl-mpm/particle:mp-mass mp))
                                                          (volume (cl-mpm/particle:mp-volume mp))
                                                          (strain-inc (* mp-strain weight volume)))
                                                     (sb-thread:with-mutex (lock)
                                                       (progn
                                                         (incf strain strain-inc)))))))
      ;; (cl-mpm::iterate-over-neighbours-general mesh (aref mps i) order
      ;;                                          (lambda (mesh mp id distance)
      ;;                                            (with-accessors ((strain cl-mpm/mesh::node-strain-energy-density)
      ;;                                                             (lock cl-mpm/mesh::node-lock))
      ;;                                                (cl-mpm/mesh::get-node mesh id)
      ;;                                              (let* ((mp-strain (cl-mpm/particle::mp-strain-energy-density mp))
      ;;                                                     (weight (weight-func distance length-parameter))
      ;;                                                     (mass (cl-mpm/particle:mp-mass mp))
      ;;                                                     (volume (cl-mpm/particle:mp-mass mp))
      ;;                                                     (strain-inc (* mp-strain weight volume)))
      ;;                                                (sb-thread:with-mutex (lock)
      ;;                                                  (progn
      ;;                                                    (incf strain strain-inc)))))))
      )))
(defun gather-local-energy (sim mps length-parameter order)
  (let ((mesh (cl-mpm:sim-mesh sim)))
    (lparallel:pdotimes (i (length mps))
      (let ((total-energy 0)
            (mp (aref mps i)))
        (cl-mpm::iterate-over-neighbours-shape-bspline mesh mp
                                                       (lambda (mesh mp node weight dsvp)
                                                   (with-accessors ((strain cl-mpm/mesh::node-strain-energy-density)
                                                                    (mass cl-mpm/mesh::node-mass)
                                                                    )
                                                       node
                                                     (when (> mass (cl-mpm:sim-mass-filter sim))
                                                       (incf total-energy (/ (* strain weight) 1))))))
        ;; (cl-mpm::iterate-over-neighbours-general mesh mp order
        ;;                                          (lambda (mesh mp id distance)
        ;;                                            (with-accessors ((strain cl-mpm/mesh::node-strain-energy-density)
        ;;                                                             (mass cl-mpm/mesh::node-mass)
        ;;                                                             )
        ;;                                                (cl-mpm/mesh::get-node mesh id)
        ;;                                              (when (> mass (cl-mpm:sim-mass-filter sim))
        ;;                                                (let* ((weight (weight-func distance length-parameter)))
        ;;                                                  (incf total-energy (/ (* strain weight) 1))
        ;;                                                  )))))
        (setf (cl-mpm/particle::mp-strain-energy-density mp) total-energy)))))
(defun delocalise-strain-energy (sim mps)
  "Take some strain energy value and generate a delocalised variation"
  (lparallel:pmapcar (lambda (mp)
                       (setf (cl-mpm/particle::mp-strain-energy-density-local mp)
                             (cl-mpm/particle::mp-strain-energy-density mp))) mps)
  (let ((order 1)
        (length-parameter (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh sim)))
        )
    (scatter-local-energy sim mps length-parameter order)
    (gather-local-energy sim mps length-parameter order)
    )
  ;; (let ((potential-mps (lparallel:premove-if-not
  ;;                       (lambda (mp) (> (cl-mpm/particle::mp-strain-energy-density mp) 0))
  ;;                       mps)))
  ;;   (multiple-value-bind (energy mp) (loop for mp in potential-mps
  ;;                                          maximizing (values (cl-mpm/particle::mp-strain-energy-density mp) mp))
  ;;     (progn
  ;;         (setf (cl-mpm/particle::mp-damage mp) 1)
  ;;       )))
  )

(defun errode-material (mps)
  "Damages mps that exceeed fracture threshold"
    (lparallel:pdotimes (i (length mps))
      (with-accessors ((gc cl-mpm/particle::mp-fracture-toughness)
                       (damage cl-mpm/particle:mp-damage)
                       (volume cl-mpm/particle:mp-volume)
                       (strain-energy cl-mpm/particle::mp-strain-energy-density))
          (aref mps i)
        (when (> (* 1 strain-energy) gc)
          (setf damage 1)
          ))))

(defun remove-material (mps)
  "Remove material points that have strain energy density exceeding fracture toughness"
  (delete-if (lambda (mp)
               (with-accessors ((gc cl-mpm/particle::mp-fracture-toughness)
                                (damage cl-mpm/particle:mp-damage)
                                (volume cl-mpm/particle:mp-volume)
                                (strain-energy cl-mpm/particle::mp-strain-energy-density))
                   mp
                 (> strain-energy gc)))
             mps))
(defun update-fracture (sim)
  (with-accessors ((mps sim-mps)
                   (mesh sim-mesh))
      sim
    ;; Caclulate local strain value
    (calculate-strain-energy mps)
    ;; (delocalise-strain-energy sim mps)
    (errode-material mps)
    ;; (remove-material mps)
    ;; (reset-strain mesh)
    ;; Apply some nonlocal smoothing?
    ))









