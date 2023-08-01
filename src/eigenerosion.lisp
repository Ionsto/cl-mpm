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
(declaim (ftype (function (cl-mpm/particle:particle) (values double-float)) calculate-strain-energy-mp))
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
        (map-into l (lambda (sii) (max sii 0)) l)
        (let* ((driving-strain (matrix-to-voight (magicl:@ v (magicl:from-diag l :type 'double-float) (magicl:transpose v))))
               (ss (magicl.simd::.*-simd stress driving-strain))
               (energy (* 0.5d0 (+ (magicl:tref ss 0 0) (magicl:tref ss 1 0) (magicl:tref ss 2 0) (magicl:tref ss 2 0)))))
          (values energy))))))
(defun calculate-strain-energy (mps)
  ;; (lparallel:pmapcar #'calculate-strain-energy-mp mps)
  (lparallel:pdotimes (i (length mps)) 
    (let ((mp (aref mps i)))
      (setf (cl-mpm/particle::mp-strain-energy-density mp) (calculate-strain-energy-mp mp)))))
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

(defparameter *priority-queue* (lparallel.queue:make-queue))
(defun delocalise-strain-energy (sim mps)
  "Run the delocalised strain energy eigenerosion"
  ;Find the local strain of each MP
  (lparallel:pmapcar (lambda (mp)
                       (with-accessors ((strain-energy cl-mpm/particle::mp-strain-energy-density-local)) mp
                         (setf strain-energy
                               (max strain-energy (calculate-strain-energy-mp mp))))) mps)
  ;Using a general length-param
  (let ((length-parameter (* 1.5d0 (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh sim)))))
    (let* (
          (pq *priority-queue*)
          (nmps (length mps))
          ;; (neighbours (find-neighbours mps length-parameter))
          (neighbours (ensure-neighbour-array nmps))
          (neighbours-count (ensure-neighbour-count nmps))
          (fracture-count 3))
      ;We only allow a finite amount of fracture per timestep
      (loop for f from 0 to fracture-count
            do
               ;Update each mps local strain
               (lparallel:pdotimes (i nmps)
                 (let* ((mp (aref mps i))
                       (mpneighbour-count (aref neighbours-count i))
                       (local-energy (with-accessors ((mass cl-mpm/particle:mp-mass)
                                                       (seng cl-mpm/particle::mp-strain-energy-density-local)) mp
                                        (* mass seng)))
                       (mass-total (cl-mpm/particle:mp-mass mp)))
                   (loop for mp-i from 0 below mpneighbour-count
                         ;Get the index of the mp-ith neighbour of mp i
                         do (let* ((mp-id (aref neighbours i mp-i))
                                   (mp-b (aref mps mp-id))
                                   (mp-mass (cl-mpm/particle::mp-mass mp)))
                              (incf local-energy (* mp-mass (cl-mpm/particle::mp-strain-energy-density-local mp-b)))
                              (incf mass-total mp-mass)))
                   (when (> mass-total 0)
                     (setf local-energy (/ local-energy mass-total))
                     (with-accessors ((actual-strain-energy cl-mpm/particle::mp-strain-energy-density))
                         mp
                       (setf actual-strain-energy (max actual-strain-energy local-energy)))
                     (when (< (cl-mpm/particle::mp-damage mp) 1)
                       (lparallel.queue:push-queue (list local-energy i) pq)))
                   ))
               (let* ((mp-fracture-list (loop while (not (lparallel.queue:queue-empty-p pq))
                                              collect (lparallel.queue:pop-queue pq)))
                      (max-e 0)
                      (max-mp nil)
                      (to-fracture (loop for fl in mp-fracture-list
                                         do (when (> (first fl) max-e)
                                              (setf max-e (first fl))
                                              (setf max-mp (second fl))))))
                 (when max-mp
                   (when (> max-e (cl-mpm/particle::mp-fracture-toughness (aref mps max-mp))) 
                     (with-accessors ((damage cl-mpm/particle:mp-damage)
                                      (local-strain cl-mpm/particle::mp-strain-energy-density-local))
                         (aref mps max-mp)
                       (setf damage 1)
                       ;; (setf local-strain 0)
                       )))
                 ))
      ;; (errode-material)
      (remove-material-damaged sim)
      )))

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

(defun remove-material-damaged (sim)
  "Remove material points that have strain energy density exceeding fracture toughness"
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (setf mps (remove-if (lambda (mp)
               (with-accessors ((gc cl-mpm/particle::mp-fracture-toughness)
                                (damage cl-mpm/particle:mp-damage)
                                (volume cl-mpm/particle:mp-volume))
                   mp
                 (>= damage 1))) mps)))
  ;; (delete-if (lambda (mp)
  ;;              (with-accessors ((gc cl-mpm/particle::mp-fracture-toughness)
  ;;                               (damage cl-mpm/particle:mp-damage)
  ;;                               (volume cl-mpm/particle:mp-volume))
  ;;                  mp
  ;;                (>= damage 1)))
  ;;            mps)
  )
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
(declaim
 (inline diff-squared)
 (ftype (function (cl-mpm/particle:particle cl-mpm/particle:particle) double-float) diff-squared))
(defun diff-squared (mp-a mp-b)
  (let ((dist (magicl:.- (cl-mpm/particle:mp-position mp-a)
                         (cl-mpm/particle:mp-position mp-b)))
        )
    (values (the double-float (magicl::sum (magicl.simd::.*-simd dist dist))))))

(defparameter *neighbour-queue* nil)
(defparameter *neighbour-array* nil)
(defparameter *neighbour-count* nil)
(defun ensure-neighbour-queue (nmps)
  (if *neighbour-queue*
      (if (>= (length *neighbour-queue*) nmps)
          *neighbour-queue*
          (adjust-array *neighbour-queue* nmps :initial-contents
                        (loop repeat (- (length *neighbour-queue*) nmps) collect (lparallel.queue:make-queue))))
      (setf *neighbour-queue* (make-array nmps :adjustable t
                                               :initial-contents
                                               (loop repeat nmps collect (lparallel.queue:make-queue))))))

(defun make-neighbour-array (nmps)
  (make-array (list nmps (max 10 (floor (/ nmps 2)))) :element-type 'fixnum :initial-element 0))

(defun ensure-neighbour-array (nmps)
  (if *neighbour-array*
      (if (>= (array-dimension *neighbour-array* 0) nmps)
          *neighbour-array*
          (setf *neighbour-array* (make-neighbour-array nmps)))
      (setf *neighbour-array* (make-neighbour-array nmps))))

(defun make-neighbour-count (nmps)
  (make-array nmps :element-type '(unsigned-byte 64) :initial-element 0))
(defun ensure-neighbour-count (nmps)
  (if *neighbour-count*
      (if (>= (length *neighbour-count*) nmps)
          *neighbour-count*
          (setf *neighbour-count* (make-neighbour-count nmps)))
      (setf *neighbour-count* (make-neighbour-count nmps))))

(declaim (ftype (function ((array cl-mpm/particle:particle) double-float) (values)) find-neighbours))
(defun find-neighbours (mps length-scale)
  "Find the neighbours and stuff it into the global array"
  (declare ((array cl-mpm/particle:particle) mps)
           (double-float length-scale))
  (let* ((nmps (length mps))
         ;; (neighbours (ensure-neighbour-queue nmps))
         (lsquared (* length-scale length-scale)))
    (let ((neighbour-array (ensure-neighbour-array nmps))
          (neighbour-count (ensure-neighbour-count nmps)))
      (declare ((simple-array (unsigned-byte 64)) neighbour-count)
               ((simple-array fixnum) neighbour-array)
               (double-float lsquared))
      ;Reset neighbour counter
      (lparallel:pdotimes (i nmps)
                          (setf (aref neighbour-count i) 0))
      ;Iterate over every mp except last
      (lparallel:pdotimes (i (- nmps 1))
        (progn
          (let ((mp-i (aref mps i)))
            ;Compare every mp ahead of us in the list
            (loop for j from (+ i 1) below nmps
                  do (let ((mp-j (aref mps j)))
                       (when (< (diff-squared mp-i mp-j) lsquared)
                                        ;Atomic update the list with less distance comparisons 
                                                  (let ((pos (sb-ext:atomic-incf (aref neighbour-count i))))
                                                    (setf (aref neighbour-array i pos) j))
                                                  (let ((pos (sb-ext:atomic-incf (aref neighbour-count j))))
                                                    (setf (aref neighbour-array j pos) i))
                                                  ))))))))
  (values))

(defun print-neighbours (mps)
  (loop for i from 0 below (length mps)
        do
           (progn
             (format t "MP: ~D : " i)
             (loop for n from 0 below (aref *neighbour-count* i)
                                      do (format t "~D " (aref *neighbour-array* i n)))
             (format t "~%"))))
(defparameter *neighbour-search-counter* 0)
(defun update-fracture (sim)
  (with-accessors ((mps sim-mps)
                   (mesh sim-mesh))
      sim
    ;; Caclulate local strain value
    ;; (calculate-strain-energy mps)
    ;; (when (<= (incf *neighbour-search-counter* -1) 0)
    ;;   (progn
    ;;     (find-neighbours mps (* 1.5d0 (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh sim))))
    ;;     (setf *neighbour-search-counter* 100)))
    (find-neighbours mps (* 1.5d0 (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh sim))))
    (delocalise-strain-energy sim mps)
    ;; (errode-material mps)
    ;; (remove-material mps)
    ;; (reset-strain mesh)
    ;; Apply some nonlocal smoothing?
    ))








