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
(defun calculate-strain-energy-mp (mp)
  (with-accessors ((stress cl-mpm/particle:mp-stress)
                   (strain cl-mpm/particle:mp-strain)
                   (damage cl-mpm/particle:mp-damage)
                   (volume cl-mpm/particle:mp-volume)
                   (strain-energy-density cl-mpm/particle::mp-strain-energy-density)
                   ) mp
    (progn
      (multiple-value-bind (l v) (magicl:eig (voight-to-matrix strain))
        (loop for i from 0 to 1
              do (let ((sii (nth i l)))
                   (when (< sii 0) (setf (nth i l) 0))))
        ;; (print l)
        ;; (print (magicl:@ (magicl:zeros '(2 2)) (magicl:from-list l '(2 1))))
        (let* ((driving-strain (matrix-to-voight (magicl:@ v (magicl:from-list (list (first l) 0 0 (second l)) '(2 2) :type 'double-float) (magicl:inv v))))
               ;(ss (magicl:.* stress driving-strain))
               (energy (* 0.5d0 volume (magicl::sum (magicl:.* stress driving-strain)))))
          ;; (print driving-strain)
          (setf strain-energy-density energy)
          )))))
(defun calculate-strain-energy (mps)
  (lparallel:pdotimes (i (length mps)) 
    (calculate-strain-energy-mp (aref mps i))))
(defun weight-func (distance length-parameter)
  (- 1 (/ (apply #'+ (mapcar #'* distance distance)) (*  length-parameter length-parameter))))
(defun reset-strain (mesh)
  (declare (cl-mpm/mesh::mesh mesh))
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
  (dotimes (i (array-total-size nodes))
    (let ((node (row-major-aref nodes i)))
      (when (> (cl-mpm/mesh::node-strain-energy-density node) 0)
        (setf (cl-mpm/mesh::node-strain-energy-density node) 0))))))
(defun delocalise-strain-energy (sim mesh mps)
  (let* ((length-parameter (* (cl-mpm/mesh::mesh-resolution mesh) 1))
         (order (floor (/ length-parameter (cl-mpm/mesh::mesh-resolution mesh)))))
    (lparallel:pdotimes (i (length mps))
      (cl-mpm::iterate-over-neighbours-general mesh (aref mps i) order
                                               (lambda (mesh mp id distance)
                                                 (with-accessors ((strain cl-mpm/mesh::node-strain-energy-density)
                                                                  (lock cl-mpm/mesh::node-lock))
                                                     (cl-mpm/mesh::get-node mesh id)
                                                   (let* ((mp-strain (cl-mpm/particle::mp-strain-energy-density mp))
                                                          (weight (weight-func distance length-parameter))
                                                          (strain-inc (* mp-strain weight (cl-mpm/particle:mp-mass mp))))
                                                     (sb-thread:with-mutex (lock)
                                                       (progn
                                                         (incf strain strain-inc))))))))
    (lparallel:pdotimes (i (length mps))
      (let ((total-energy 0)
            (mp (aref mps i)))
        (cl-mpm::iterate-over-neighbours-general mesh mp order
                                                 (lambda (mesh mp id distance)
                                                   (with-accessors ((strain cl-mpm/mesh::node-strain-energy-density)
                                                                    (mass cl-mpm/mesh::node-mass)
                                                                    )
                                                       (cl-mpm/mesh::get-node mesh id)
                                                     (when (> mass (cl-mpm:sim-mass-filter sim))
                                                       (let* ((weight (weight-func distance length-parameter)))
                                                         (incf total-energy (/ (* strain weight) mass))
                                                         )))))
        (setf (cl-mpm/particle::mp-strain-energy-density mp) total-energy)))))
(defun errode-material (mps)
    (lparallel:pdotimes (i (length mps))
      (with-accessors ((gc cl-mpm/particle::mp-fracture-toughness)
                       (damage cl-mpm/particle:mp-damage)
                       (strain-energy cl-mpm/particle::mp-strain-energy-density))
          (aref mps i)
        (when (> strain-energy gc)
          (setf damage 1)
          ))))
(defun update-fracture (sim)
  (with-accessors ((mps sim-mps)
                   (mesh sim-mesh))
      sim
    ;; Caclulate local strain value
    (calculate-strain-energy mps)
    ;; (delocalise-strain-energy sim mesh mps)
    (errode-material mps)
    ;; (reset-strain mesh)
    ;; Apply some nonlocal smoothing?

    ))









