(defpackage :cl-mpm/dynamic-relaxation-mpi
  (:use :cl
   :cl-mpm/dynamic-relaxation)
  (:export))
(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)

(defclass mpm-sim-dr-mpi (cl-mpm/dynamic-relaxation::mpm-sim-dr-ul
                          cl-mpm/mpi::mpm-sim-mpi-nodes)
  ()
  (:documentation "DR psudo-linear step with update stress last update"))

(defclass mpm-sim-quasi-static-mpi (mpm-sim-dr-mpi)
  ()
  (:documentation "DR psudo-linear step with update stress last update"))

(defclass mpm-sim-damage-quasi-static-mpi (mpm-sim-quasi-static-mpi cl-mpm/mpi::mpm-sim-mpi-damage)
  ()
  (:default-initargs
   :vel-algo :QUASI-STATIC)
  (:documentation "DR psudo-linear step with update stress last update"))

(defmethod cl-mpm/dynamic-relaxation::save-vtks ((sim cl-mpm/mpi::mpm-sim-mpi) output-dir step &optional prefix)
  (let ((pre (if prefix (format nil "_~A" prefix) "")))
    (let ((rank (cl-mpi:mpi-comm-rank)))
      (when (= rank 0)
        (format t "Save vtks ~D~%" step))
      (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim~A_~5,'0d_~5,'0d.vtk" pre rank step)) sim)
      (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim~A_nodes_~5,'0d_~5,'0d.vtk" pre rank step)) sim)
      (cl-mpm/output::save-vtk-cells (merge-pathnames output-dir (format nil "sim~A_cells_~5,'0d_~5,'0d.vtk" pre rank step)) sim)
      (cl-mpm/penalty:save-vtk-penalties (uiop:merge-pathnames* output-dir (format nil "sim~A_p_~5,'0d_~5,'0d.vtk" pre rank step)) sim ))))

(defmethod cl-mpm/dynamic-relaxation::save-vtks-dr-step ((sim cl-mpm/mpi::mpm-sim-mpi) output-dir trial-solve step iter)
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (when (= rank 0)
      (format t "Save vtks ~D~%" step))
    (let ((post (format nil "~5,'0d_~5,'0d_~5,'0d_~5,'0d.vtk" rank trial-solve step iter)))
      (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_step_~A.vtk" post)) sim)
      (cl-mpm/output:save-vtk-nodes (merge-pathnames output-dir (format nil "sim_step_nodes_~A.vtk" post)) sim)
      ;; (cl-mpm/output:save-vtk-cells (merge-pathnames output-dir (format nil "sim_step_cells_~A.vtk" post)) sim)
      (cl-mpm/penalty:save-vtk-penalties (merge-pathnames output-dir (format nil "sim_step_p_~A.vtk" post)) sim))))
