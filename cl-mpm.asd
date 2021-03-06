(in-package :asdf-user)

(defsystem "cl-mpm"
  :class :package-inferred-system
  :depends-on ("magicl"
               "alexandria"
               "array-operations"
               "lparallel"
               "symbolic-derivation"
               "cl-mpm/constitutive"
               "cl-mpm/particle"
               ;"cl-mpm/shape-function"
               "cl-mpm/bc"
               "cl-mpm/mesh"
               )
  :description ""
  :in-order-to ((test-op (load-op "test/all")))
  :perform (test-op (o c) (symbol-call :test/all :test-suite))
  :serial t
  :components (
               (:file "src/shape-function")
               (:file "src/forces")
               (:file "src/core")
               ))
(defsystem "symbolic-derivation"
  :description "Symbolic derivation library from #TODO"
  :author "Aleksander Ksiazek"
  :serial t
  :components ((:file "src/symbolic-derivation")))

(defsystem "cl-mpm/mesh"
  :depends-on ("magicl")
  :description "MPM boundary conditions"
  :serial t
  :components ((:file "src/mesh")))

(defsystem "cl-mpm/bc"
  :depends-on ("magicl"
               "cl-mpm/mesh"
               )
  :description "MPM boundary conditions"
  :serial t
  :components ((:file "src/bc")))

(defsystem "cl-mpm/setup"
  :depends-on ("cl-mpm"
               "cl-mpm/mesh"
               )
  :description "MPM test simulations"
  :serial t
  :components ((:file "src/setup")))


(defsystem "cl-mpm/constitutive"
  :depends-on ("magicl")
  :description "Various constitutive models"
  :serial t
  :components ((:file "src/constitutive")))
;
(defsystem "cl-mpm/particle"
  :depends-on ("magicl")
  :description "MPM particle definitions"
  :serial t
  :components ((:file "src/particle")))

(defsystem "cl-mpm/output"
  :depends-on ("magicl")
  :description "MPM output helper functions"
  :serial t
  :components ((:file "src/output")))
;
;(defsystem "cl-mpm/shape-function"
;  :depends-on ("magicl")
;  :description "MPM shape function definitions"
;  :serial t
;  :components ((:file "src/shape-function")))
;
;(defsystem "cl-mpm/mesh"
;  :depends-on ("magicl")
;  :description "MPM mesh definitions"
;  :serial t
;  :components ((:file "src/mesh")))
;
(defsystem "cl-mpm/test"
  :depends-on ("cl-mpm"))
