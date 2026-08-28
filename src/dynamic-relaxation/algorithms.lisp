(in-package :cl-mpm/dynamic-relaxation)
(declaim #.cl-mpm/settings:*optimise-setting*)


(defun save-test-vtks (sim &key (output-dir "./output/"))
  (cl-mpm/output:save-vtk (merge-pathnames "test.vtk" output-dir) sim)
  (cl-mpm/output:save-vtk-nodes (merge-pathnames "test_nodes_0.vtk" output-dir) sim)
  (cl-mpm/output:save-vtk-nodes (merge-pathnames "test_nodes_1.vtk" output-dir) sim)
  (cl-mpm/output:save-vtk-cells (merge-pathnames "test_cells.vtk" output-dir) sim)
  )


















