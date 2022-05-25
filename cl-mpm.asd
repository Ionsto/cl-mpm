(in-package :asdf-user)

(defsystem "cl-mpm"
  :class :package-inferred-system
  :depends-on ("magicl"
               "alexandria"
               "array-operations"
               "cl-autodiff"
               ;"cl-mpm/constitutive"
               ;"cl-mpm/particle"
               ;"cl-mpm/shape-function"
               ;"cl-mpm/mesh"
               )
  :description ""
  :in-order-to ((test-op (load-op "test/all")))
  :perform (test-op (o c) (symbol-call :test/all :test-suite))
  :serial t
  :components (
               (:file "src/symbolic-derivation")
               (:file "src/shape-function")
               (:file "src/constitutive")
               (:file "src/forces")
               (:file "src/particle")
               (:file "src/mesh")
               (:file "src/core")
               ))

;(defsystem "cl-mpm/constitutive"
;  :depends-on ("magicl")
;  :description "Various constitutive models"
;  :serial t
;  :components ((:file "src/constitutive")))
;
;(defsystem "cl-mpm/particle"
;  :depends-on ("magicl")
;  :description "MPM particle definitions"
;  :serial t
;  :components ((:file "src/particle")))
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
