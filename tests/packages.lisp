(fiasco:define-test-package #:cl-mpm-tests
    (:use :cl-mpm/utils)
    (:export
     #:run-cl-mpm-tests))

(in-package :cl-mpm-tests)
(defun run-cl-mpm-tests ()
  (in-package :cl-mpm-tests)
  (run-package-tests))
