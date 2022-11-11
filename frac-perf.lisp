(ql:quickload "cl-mpm/examples/fracture")
(cl-mpm/examples/fracture::setup)
(time (progn
        (dotimes (i 50)
          (progn
            (format t "~D~%" i)
            (cl-mpm::update-sim cl-mpm/examples/fracture::*sim*)
            ))))
