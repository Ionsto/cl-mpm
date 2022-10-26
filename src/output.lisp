(defpackage :cl-mpm/output
  (:use :cl)
  (:export
   #:save-vtk-mesh
   #:save-vtk))

(in-package :cl-mpm/output)

(defun save-sim (sim)
  "Save a simulation such that it can be restarted from that point")
(defun load-sim (sim))

(defun sample-line (sim start end steps get-value)
	"Sample over a line from start to end with fixed steps count, get-value function gets given a sampling mp and a mesh and should return a scalar"
    (let (  (sample-mp (cl-mpm/particle:make-particle 2 :pos start))
            (direction (magicl:from-list (mapcar (lambda (s e) (/ (- e s) steps)) start end) '(2 1))))
      (loop :repeat (+ steps 1) collect 
              (prog1
                (funcall get-value (cl-mpm:sim-mesh sim) sample-mp)
                (setf (cl-mpm/particle:mp-position sample-mp) (magicl:.+ (cl-mpm/particle:mp-position sample-mp) direction))))))

(defun sample-point (sim point get-value)
    (let ((sample-mp (cl-mpm/particle:make-particle 2 :pos point)))
        (funcall get-value (cl-mpm:sim-mesh sim) sample-mp)))

(defmacro scalar-average (accessor)
    `(lambda (mesh mp)
        (let ((av 0))
          (progn
            (cl-mpm::iterate-over-neighbours mesh mp 
                 (lambda (mesh mp node svp dsvp) 
                   (setf av (+ av (* svp (,accessor node))))))
            av))))

(defmacro matrix-average (accessor shape)
    `(lambda (mesh mp)
        (let ((av (magicl:zeros ,shape)))
          (progn
            (cl-mpm::iterate-over-neighbours mesh mp 
             (lambda (mesh mp node svp dsvp) 
               (setf av (magicl:.+ av (magicl:scale (,accessor node) svp)))))
            av))))

(defun sample-line-mass (sim start end steps)
    (sample-line sim start end steps 
        (scalar-average cl-mpm/mesh:node-mass)))

(defun sample-line-velocity (sim start end steps)
    (sample-line sim start end steps 
        (matrix-average cl-mpm/mesh:node-velocity '(2 1))))

(defun sample-point-mass (sim point)
    (sample-point sim point 
        (scalar-average cl-mpm/mesh:node-mass)))

(defun sample-point-velocity (sim point)
    (sample-point sim point
        (matrix-average cl-mpm/mesh:node-velocity '(2 1))))
(defun format-scalar (stream name id mps accessor)
  (format stream "SCALARS ~a FLOAT ~d~%" name 1)
  (format stream "LOOKUP_TABLE default~%")
    (loop for mp across mps
        do (format stream "~f ~%" (funcall accessor mp)))
  (format stream "~%")
  )
(defmacro save-parameter (name accessor)
  `(progn
     (format-scalar fs ,name id mps (lambda (mp) ,accessor))
     (incf id)))
(defun save-vtk-mesh (filename sim)
  (with-accessors ((mesh cl-mpm:sim-mesh)) sim
    (with-open-file (fs filename :direction :output :if-exists :supersede)
      (format fs "# vtk DataFile Version 2.0~%")
      (format fs "Lisp generated vtk file, SJVS~%")
      (format fs "ASCII~%")
      (format fs "DATASET UNSTRUCTURED_GRID~%")
      (with-accessors ((nodes cl-mpm/mesh::mesh-nodes)
                       (size cl-mpm/mesh::mesh-count)
                       (h cl-mpm/mesh::mesh-resolution)) mesh
        (format fs "POINTS ~d double~%" (floor (apply #'* size)))
        (loop for x from 0 to (- (first size) 1)
              do
                 (loop for y from 0 to (- (second size) 1)
                       do
                          (format fs "~f ~f ~f ~%" (* h x) (* h y)
                                  ;(magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                                  ;(magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                                  0)))
        (let ((nels (floor (* (- (first size) 1) (- (second size) 1)))))
          (format fs "CELLS ~D ~D~%"
                  nels
                  (floor (* 5 nels)))
          (flet ((id (x y) (floor (+ y (* x (second size))))))
            (loop for x from 0 to (- (first size) 2)
                  do
                     (loop for y from 0 to (- (second size) 2)
                           do
                              (format fs "~D ~D ~D ~D ~D ~%"
                                      4
                                      (id x y)
                                      (id x (+ y 1))
                                      (id (+ x 1) y)
                                      (id (+ x 1) (+ y 1))
                                      ))))

          (format fs "CELL_TYPES ~d~%" nels)
          (loop repeat nels
                do (format fs "~D~%" 8))
          ;; (loop for x from 0 to (- (first size) 1)
          ;;       do
          ;;          (loop for y from 0 to (- (second size) 1)
          ;;                do
          ;;                   (format fs "~D~%" 8)))
          )))))
(defun save-vtk (filename sim)
  (with-accessors ((mps cl-mpm:sim-mps)) sim
    (with-open-file (fs filename :direction :output :if-exists :supersede)
      (format fs "# vtk DataFile Version 2.0~%")
      (format fs "Lisp generated vtk file, WMC~%")
      (format fs "ASCII~%")
      (format fs "DATASET UNSTRUCTURED_GRID~%")
      (format fs "POINTS ~d double~%" (length mps))
      (loop for mp across mps
            do (format fs "~f ~f ~f ~%"
                       (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                       (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                       0))
      (format fs "~%")
      (let ((id 1))
        (declare (special id))
        (format fs "POINT_DATA ~d~%" (length mps))
        (save-parameter "damage" (cl-mpm/particle:mp-damage mp))
        (save-parameter "mass" (cl-mpm/particle:mp-mass mp))
        (save-parameter "stress_xx" (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0))
        (save-parameter "stress_yy" (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0))
        (save-parameter "stress_xy" (magicl:tref (cl-mpm/particle:mp-stress mp) 2 0))
        (save-parameter "strain_energy" (cl-mpm/particle::mp-strain-energy-density mp))
        ;; (format-scalar fs "damage" id mps (lambda (mp) (cl-mpm/particle:mp-damage mp)))
        ;; (format-scalar fs "stress_x" id mps (lambda (mp) (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0)))
        ;; (format-scalar fs "stress_y" id mps (lambda (mp) (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0)))
        ;; (format-scalar fs "stress_xy" id mps (lambda (mp) (magicl:tref (cl-mpm/particle:mp-stress mp) 2 0)))
        )
      )))
