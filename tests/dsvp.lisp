(in-package :cl-mpm-tests)


(defun test-dsvp ()
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
       :bottom '(nil nil nil))
      (cl-mpm:update-sim sim)
      ;; (cl-mpm::iterate-over-mps-serial
      ;;  (cl-mpm:sim-mps sim))
      (let ((mp (aref (cl-mpm:sim-mps sim) 0)))
        (format t "MP at ~A ~%" (cl-mpm/utils:fast-storage (cl-mpm/particle::mp-position mp)))
        (cl-mpm::iterate-over-neighbours
         (cl-mpm:sim-mesh sim)
         mp
         (lambda (mesh mp node svp grads fsvp fgrads)
           (format t "Node ~A - svp ~E - grads ~E - dsvp~%"
                   (cl-mpm/mesh:node-index node) svp grads)
           (pprint (cl-mpm/shape-function::assemble-dsvp-3d grads))
           (format t "~%")
           )))
      (setf (cl-mpm/utils:varef (cl-mpm/particle::mp-stress-kirchoff (aref (cl-mpm:sim-mps sim) 0)) 1) 1d0)
      (setf (cl-mpm/utils:varef (cl-mpm/particle::mp-stress (aref (cl-mpm:sim-mps sim) 0)) 1) 1d0)
      ;; (cl-mpm:update-sim sim)
      (cl-mpm::p2g-force (cl-mpm:sim-mesh sim) (cl-mpm:sim-mps sim))
      (cl-mpm::update-node-forces sim)
      (cl-mpm/output::save-vtk-nodes "./output/nodes_0.vtk" sim)
      (cl-mpm/output::save-vtk-nodes "./output/nodes_1.vtk" sim)
      (cl-mpm/output::save-vtk "./output/mps_0.vtk" sim)
      (cl-mpm/output::save-vtk "./output/mps_1.vtk" sim)
      )))
