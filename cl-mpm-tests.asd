(in-package :asdf-user)
(defsystem "cl-mpm-tests"
  :depends-on ("cl-mpm/all"
               "lisp-stat"
               "fiasco")
  :description "cl-mpm tests"
  :serial t
  :components ((:module "tests"
                :serial t
                :components
                ((:file "packages")
                 (:file "constitutive")
                 (:file "cpp")
                 ))))
