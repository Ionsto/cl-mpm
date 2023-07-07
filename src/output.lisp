(defpackage :cl-mpm/output
  (:use :cl)
  (:export
   #:save-vtk-mesh
   #:save-vtk
   #:save-csv
   ))

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
          do (format stream "~E ~%"
                     (coerce (funcall accessor mp) 'single-float)))
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
                       (h cl-mpm/mesh::mesh-resolution)
                       (border cl-mpm/mesh::mesh-boundary-order)) mesh
        (format fs "POINTS ~d double~%" (floor (apply #'* size)))
        (loop for x from 0 to (- (first size) 1)
              do
                 (loop for y from 0 to (- (second size) 1)
                       do
                          (format fs "~E ~E ~E ~%"
                                  (coerce (* h (- x border)) 'single-float)
                                  (coerce (* h (- y border)) 'single-float)
                                  0e0)))
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
            do (format fs "~E ~E ~E ~%"
                       (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 1 0) 'single-float)
                       0e0))
      (format fs "~%")
      (let ((id 1))
        (declare (special id))
        (format fs "POINT_DATA ~d~%" (length mps))

        (save-parameter "mass" (cl-mpm/particle:mp-mass mp))
        (save-parameter "density" (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)))
        (save-parameter "index" (cl-mpm/particle::mp-index mp))
        (save-parameter "vel_x" (magicl:tref (cl-mpm/particle:mp-velocity mp) 0 0))
        (save-parameter "vel_y" (magicl:tref (cl-mpm/particle:mp-velocity mp) 1 0))
        (save-parameter "acc_x" (magicl:tref (cl-mpm/particle::mp-acceleration mp) 0 0))
        (save-parameter "acc_y" (magicl:tref (cl-mpm/particle::mp-acceleration mp) 1 0))
        (save-parameter "disp_x" (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0))
        (save-parameter "disp_y" (magicl:tref (cl-mpm/particle::mp-displacement mp) 1 0))
        (save-parameter "sig_xx" (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0))
        (save-parameter "sig_yy" (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0))
        (save-parameter "sig_xy" (magicl:tref (cl-mpm/particle:mp-stress mp) 2 0))

        (save-parameter "e_xx" (magicl:tref (cl-mpm/particle::mp-strain mp) 0 0))
        (save-parameter "e_yy" (magicl:tref (cl-mpm/particle::mp-strain mp) 1 0))
        (save-parameter "e_xy" (magicl:tref (cl-mpm/particle::mp-strain mp) 2 0))
        (save-parameter "temp" (magicl:tref (cl-mpm/particle::mp-velocity-rate mp) 2 0))
        ;; (save-parameter "viscosity" (cl-mpm/particle::mp-true-visc mp))

        (save-parameter "strain_rate"
                        (multiple-value-bind (l v)
                            (magicl:eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle::mp-velocity-rate mp)))
                          (reduce #'+ (mapcar #'* l l))))
        (save-parameter "pressure" (cl-mpm/particle::mp-pressure mp))
        ;; (save-parameter "pressure" (/ (+ (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0)
        ;;                                 (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0)) 2d0))
        (labels ((dot (a b) (magicl::sum (magicl:.* a b)))
                 (norm (a) (magicl:scale a (/ 1d0 (sqrt (dot a a)))))
                 (radial-stress (mp)
                   (with-accessors ((stress cl-mpm/particle:mp-stress)
                                    (pos cl-mpm/particle:mp-position))
                       mp
                     (let ((normal (norm (magicl:.- pos (magicl:from-list '(250d0 250d0) '(2 1))))))
                       (dot normal
                            (magicl:@ (cl-mpm/utils:voight-to-matrix stress)
                                      normal))))))
          (save-parameter "s_rr" (radial-stress mp)))

        (save-parameter "s_1"
                        (multiple-value-bind (l v) (magicl:eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))
                          (loop for sii in l maximize sii)))
        (save-parameter "s_vm"
                        (multiple-value-bind (l v) (magicl:eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))

                          (let* ((l (sort l #'>))
                                 (s_1 (max 0 (- (first l) (cl-mpm/particle::mp-pressure mp))))
                                 (s_2 (max 0 (- (second l) (cl-mpm/particle::mp-pressure mp)))))

                            (* (sqrt (/ 3 4)) (- s_1 s_2))
                            )
                          ))

        (save-parameter "EPS"
                        (multiple-value-bind (l v) (magicl:eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))
                          (- (loop for sii in l maximize sii) (cl-mpm/particle::mp-pressure mp))))
        (save-parameter "size_x" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 0 0))
        (save-parameter "size_y" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0))
        (save-parameter "damage"
                        (if (slot-exists-p mp 'cl-mpm/particle::damage)
                            (cl-mpm/particle:mp-damage mp)
                            0d0))
        (save-parameter "damage_inc"
                        (if (slot-exists-p mp 'cl-mpm/particle::damage)
                            (cl-mpm/particle::mp-damage-increment mp)
                            0d0))
        (save-parameter "damage_ybar"
                        (if (slot-exists-p mp 'cl-mpm/particle::damage)
                            (cl-mpm/particle::mp-damage-ybar mp)
                            0d0))
        (save-parameter "local_length"
                        (if (slot-exists-p mp 'cl-mpm/particle::damage)
                            (cl-mpm/particle::mp-true-local-length mp)
                            0d0))
        )
      )))

(defun save-csv (filename sim)
  (with-accessors ((mps cl-mpm:sim-mps)) sim
    (with-open-file (fs filename :direction :output :if-exists :supersede)
      (format fs "coord_x,coord_y,stress_xx,stress_yy,tau_xy,velocity_x,velocity_y,stress_1,eps,damage,lx,ly~%")
      (loop for mp across mps
            do (format fs "~E, ~E, ~F, ~F, ~F, ~F, ~F, ~F, ~F, ~F, ~F, ~F ~%"
                       (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle:mp-position mp) 1 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle:mp-stress mp) 2 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle:mp-velocity mp) 0 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle:mp-velocity mp) 1 0) 'single-float)
                       (coerce (multiple-value-bind (l v)
                                   (magicl:eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))
                                 (loop for sii in l maximize sii)) 'single-float)
                       (coerce
                        (multiple-value-bind (l v) (magicl:eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))
                          (- (apply #'max l) (cl-mpm/particle::mp-pressure mp))) 'single-float)
                       (coerce (if (slot-exists-p mp 'cl-mpm/particle::damage)
                            (cl-mpm/particle:mp-damage mp)
                            0d0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle::mp-domain-size mp) 0 0) 'single-float)
                       (coerce (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0) 'single-float)
                       ))
      )))
