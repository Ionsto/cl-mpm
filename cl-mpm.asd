(in-package :asdf-user)

(defsystem "symbolic-derivation"
  :description "Symbolic derivation library from #TODO"
  :author "Aleksander Ksiazek"
  :serial t
  :components ((:file "src/symbolic-derivation")))

(defsystem "cl-mpm/utils"
  :depends-on ("magicl")
  :description "MPM utility functiosn definitions"
  :serial t
  :components ((:file "src/utils")))

(defsystem "cl-mpm/fastmath"
  :depends-on ("magicl"
               :sb-simd
               )
  :description "MPM fast maths operations definitions"
  :serial t
  :components ((:file "src/fastmath")))

(defsystem "cl-mpm/mesh"
  :depends-on ("magicl"
               "cl-mpm/shape-function")
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
  :depends-on ("magicl"
               "cl-mpm/utils")
  :description "Various constitutive models"
  :serial t
  :components ((:file "src/constitutive")))
;
(defsystem "cl-mpm/particle"
  :depends-on ("magicl")
  :description "MPM particle definitions"
  :serial t
  :components ((:file "src/particle")))

(defsystem "cl-mpm/damage"
  :depends-on ("magicl"
               "cl-mpm/utils"
               "cl-mpm/particle")
  :description "MPM smeared damage mechanics"
  :serial t
  :components ((:file "src/damage")))
(defsystem "cl-mpm/eigenerosion"
  :depends-on ("magicl"
               "cl-mpm"
               "cl-mpm/utils"
               "cl-mpm/particle")
  :description "MPM eigenerosion damage mechanics"
  :serial t
  :components ((:file "src/eigenerosion")))

(defsystem "cl-mpm/output"
  :depends-on ("magicl")
  :description "MPM output helper functions"
  :serial t
  :components ((:file "src/output")))
;
(defsystem "cl-mpm/shape-function"
 :depends-on ("magicl")
 :description "MPM shape function definitions"
 :serial t
 :components ((:file "src/shape-function")))
;
;(defsystem "cl-mpm/mesh"
;  :depends-on ("magicl")
;  :description "MPM mesh definitions"
;  :serial t
;  :components ((:file "src/mesh")))
;
(defsystem "cl-mpm"
  ;; :class :package-inferred-system
  :depends-on ("magicl"
               "cl-mpm/fastmath"
               "cl-mpm/utils"
               "alexandria"
               "array-operations"
               "lparallel"
               "symbolic-derivation"
               "cl-mpm/constitutive"
               "cl-mpm/particle"
               "cl-mpm/shape-function"
               "cl-mpm/bc"
               "cl-mpm/mesh"
               "cl-mpm/damage"
               )
  :description ""
  :in-order-to ((test-op (load-op "test/all")))
  :perform (test-op (o c) (symbol-call :test/all :test-suite))
  :serial t
  :components (
               ;; (:file "src/shape-function")
               (:file "src/forces")
               (:file "src/core")
               ))
(defsystem "cl-mpm/buoyancy"
  :depends-on ("cl-mpm"
               "cl-mpm/bc")
  :description ""
  :components (
               (:file "src/buoyancy")
               ))

(defsystem "cl-mpm/test"
  :depends-on ("cl-mpm"))
(defsystem "cl-mpm/example"
  :class :package-inferred-system
  :depends-on ("cl-mpm")
  :serial t
  :components ((:file "examples/bounce")))

(defsystem "cl-mpm/examples/column"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/eigenerosion"
               "vgplot"
               "swank.live"
               "magicl")
  :serial t
  :components ((:file "examples/column")))

(defsystem "cl-mpm/examples/fracture"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/eigenerosion"
               "vgplot"
               "swank.live"
               "magicl")
  :serial t
  :components ((:file "examples/fracture")))

(defsystem "cl-mpm/examples/plate-with-hole"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/eigenerosion"
               "vgplot"
               "swank.live"
               "magicl")
  :serial t
  :components ((:file "examples/plate-with-hole")))

(defsystem "cl-mpm/examples/slump"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/buoyancy"
               "cl-mpm/eigenerosion"
               "vgplot"
               "swank.live"
               "magicl")
  :serial t
  :components ((:file "examples/slump")))
(defsystem "cl-mpm/examples/notch"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/buoyancy"
               "vgplot"
               "swank.live"
               "magicl")
  :serial t
  :components ((:file "examples/notch")))
(defsystem "cl-mpm/examples/float"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/buoyancy"
               "vgplot"
               "swank.live"
               "magicl")
  :serial t
  :components ((:file "examples/float")))
(defsystem "cl-mpm/examples/flow"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/eigenerosion"
               "vgplot"
               "swank.live"
               "magicl")
  :serial t
  :components ((:file "examples/flow")))

(defsystem "cl-mpm/examples/beam"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/eigenerosion"
               "vgplot"
               "swank.live"
               "magicl")
  :serial t
  :components ((:file "examples/beam")))

(defsystem "cl-mpm/examples/shear"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "vgplot"
               "swank.live"
               "magicl")
  :serial t
  :components ((:file "examples/shear")))

(defsystem "cl-mpm/examples/pullout"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/buoyancy"
               "cl-mpm/eigenerosion"
               "vgplot"
               "swank.live"
               "magicl")
  :serial t
  :components ((:file "examples/pullout")))
