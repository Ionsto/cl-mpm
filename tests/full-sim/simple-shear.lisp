(in-package :cl-mpm-tests)

(defun apply-simple-shear (sim shear-rate)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (node)
       (with-accessors ((pos cl-mpm/mesh::node-position)
                        (vel cl-mpm/mesh::node-velocity)
                        (active cl-mpm/mesh::node-active)
                        (acc cl-mpm/mesh::node-acceleration)
                        )
           node
         (when active
           (setf (cl-mpm/utils:varef vel 0) (* shear-rate (magicl:tref pos 1 0))
                 (cl-mpm/utils:varef vel 1) 0d0)))))
    ))

(defun test-simple-shear ()
  (let ((h 1d0))
    (let ((sim (cl-mpm/setup::make-simple-sim
                h
                (list 1 1)
                :sim-type 'cl-mpm:mpm-sim-usf)))
      (cl-mpm::add-mps
       sim
       (cl-mpm/setup::make-block-mps
        (list 0 0)
        (list 1 1)
        (list 2 2)
        1d0
        'cl-mpm/particle::particle-elastic
        :E 1d0
        :nu 0.2d0
        ))
      (cl-mpm/setup::setup-bcs
       sim
       :left   '(nil nil nil)
       :right  '(nil nil nil)
       :top    '(nil nil nil)
       :bottom '(0 0 nil))
      (setf (cl-mpm:sim-dt sim) 1d-3)
      (format t "dt ~E ~%" (cl-mpm:sim-dt sim))
      (cl-mpm:add-bcs
       sim
       (cl-mpm/bc::make-bc-closure
        '(0 0 0)
        (lambda ()
          (apply-simple-shear
           sim
           1d0))))
      (dotimes (i 100)
        (cl-mpm:update-sim sim)
        (pprint (cl-mpm/particle:mp-deformation-gradient (aref (cl-mpm:sim-mps sim) 0)))
        (cl-mpm/output::save-vtk-nodes (merge-pathnames "./output/" (format nil "sim_nodes_~5,'0d.vtk" i)) sim)
        (format t "~%")
        )
      (pprint (cl-mpm/particle:mp-stress (aref (cl-mpm:sim-mps sim) 0)))
      (format t "~%")
      (pprint (cl-mpm/particle::mp-elastic-matrix (aref (cl-mpm:sim-mps sim) 0)))
      (format t "~%")
      )))
