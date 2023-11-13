(in-package :asdf-user)

(defsystem "symbolic-derivation"
  :description "Symbolic derivation library from #TODO"
  :author "Aleksander Ksiazek"
  :serial t
  :components ((:file "src/symbolic-derivation")))

(defsystem "cl-mpm/utils"
  :depends-on ("magicl"
               "array-operations"
               )
  :description "MPM utility functions definitions"
  :serial t
  :components ((:file "src/utils")))

(defsystem "cl-mpm/fastmath"
  :depends-on ("magicl"
               "cl-mpm/utils"
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

(defsystem "cl-mpm/forces"
  :depends-on ("magicl"
               "cl-mpm/utils"
               "cl-mpm/particle"
               "cl-mpm/mesh"
               )
  :description "MPM forces calculations"
  :serial t
  :components ((:file "src/forces")))

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
  :description "MPM setup system"
  :serial t
  :components ((:file "src/setup")))


(defsystem "cl-mpm/constitutive"
  :depends-on ("magicl"
               "cl-mpm/utils"
               "cl-mpm/fastmath")
  :description "Various constitutive models"
  :serial t
  :components ((:file "src/constitutive")))
;
(defsystem "cl-mpm/particle"
  :depends-on ("magicl"
               "cl-mpm/constitutive"
               )
  :description "MPM particle definitions"
  :serial t
  :components ((:file "src/particle")))

(defsystem "cl-mpm/damage"
  :depends-on ("magicl"
               "cl-mpm/utils"
               "cl-mpm/constitutive"
               "cl-mpm/output"
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
  :depends-on ("magicl"
               "jonathan"
               "str")
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
               ;"cl-mpm/damage"
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
(defsystem "cl-mpm/mpi"
  :depends-on ("cl-mpm"
               "lfarm-client"
               "lfarm-server"
               "lfarm-admin"
               "cl-mpm/fastmath"
               "cl-mpm/utils"
               "alexandria"
               "array-operations"
               "lparallel"
               "cl-intbytes"
               "ieee-floats"
               "symbolic-derivation"
               "cl-mpm/constitutive"
               "cl-mpm/particle"
               "cl-mpm/shape-function"
               "cl-mpm/bc"
               "cl-mpm/mesh"
               "cl-mpm/damage"
               "cl-mpi"
               "cl-mpi-extensions"
               "trivial-with-current-source-form"
               )
  :description ""
  :components (
               (:file "src/mpi")
               ))
(defsystem "cl-mpm/buoyancy"
  :depends-on ("cl-mpm"
               "cl-mpm/bc")
  :description ""
  :components (
               (:file "src/buoyancy")
               ))

(defsystem "cl-mpm/penalty"
  :depends-on ("cl-mpm"
               "cl-mpm/bc")
  :description ""
  :components (
               (:file "src/penalty")
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
               "cl-mpm/forces"
               "cl-mpm/output"
               "cl-mpm/eigenerosion"
               "vgplot"
               "swank.live"
               "magicl")
  :serial t
  :components ((:file "examples/column")))


(defsystem "cl-mpm/plotter"
  :depends-on ("cl-mpm"
               "cl-mpm/particle"
               "vgplot")
  :serial t
  :components ((:file "src/plotter")))

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
               "cl-mpm/penalty"
               "cl-mpm/eigenerosion"
               "cl-mpm/damage"
               "cl-mpm/mpi"
               "lfarm-client"
               "vgplot"
               "swank.live"
               "lisp-stat"
               "magicl")
  :serial t
  :components ((:file "examples/slump")))
(defsystem "cl-mpm/examples/notch"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/buoyancy"
               "cl-mpm/penalty"
               "cl-mpm/damage"
               "lfarm-client"
               "vgplot"
               "lisp-stat"
               "swank.live"
               "magicl")
  :serial t
  :components ((:file "examples/notch")))

(defsystem "cl-mpm/examples/single-crack"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/buoyancy"
               "cl-mpm/penalty"
               "cl-mpm/damage"
               "lfarm-client"
               "vgplot"
               "lisp-stat"
               "swank.live"
               "magicl")
  :serial t
  :components ((:file "examples/single-crack")))

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
               "cl-mpm/damage"
               "cl-mpm/penalty"
               "cl-mpm/plotter"
               "cl-mpm/eigenerosion"
               "vgplot"
               "lisp-stat"
               "swank.live"
               "magicl")
  :serial t
  :components ((:file "examples/pullout")))

(defsystem "cl-mpm/examples/split-test"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/buoyancy"
               "cl-mpm/damage"
               "cl-mpm/penalty"
               "cl-mpm/plotter"
               "cl-mpm/eigenerosion"
               "vgplot"
               "lisp-stat"
               "swank.live"
               "magicl")
  :serial t
  :components ((:file "examples/split-test")))


(defsystem "cl-mpm/examples/creep"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/buoyancy"
               "cl-mpm/plotter"
               "cl-mpm/damage"
               "cl-mpm/eigenerosion"
               "vgplot"
               "lisp-stat"
               "swank.live"
               "magicl")
  :serial t
  :components ((:file "examples/creep")))

(defsystem "cl-mpm/examples/sliding"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/buoyancy"
               "cl-mpm/penalty"
               "cl-mpm/damage"
               "lfarm-client"
               "vgplot"
               "swank.live"
               "lisp-stat"
               "magicl")
  :serial t
  :components ((:file "examples/sliding")))

(defsystem "cl-mpm/examples/collapse"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/buoyancy"
               "cl-mpm/penalty"
               "cl-mpm/damage"
               "cl-mpm/plotter"
               "vgplot"
               "swank.live")
  :serial t
  :components ((:file "examples/collapse")))



(defsystem "cl-mpm/examples/chalk"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/buoyancy"
               "cl-mpm/penalty"
               "cl-mpm/damage"
               "cl-mpm/plotter"
               "vgplot"
               "swank.live")
  :serial t
  :components ((:file "examples/chalk")))


(defsystem "cl-mpm/examples/tpb"
  :depends-on ("cl-mpm"
               "magicl"
               "lisp-stat"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/buoyancy"
               "cl-mpm/penalty"
               "cl-mpm/damage"
               "cl-mpm/plotter"
               "vgplot"
               "swank.live")
  :serial t
  :components ((:file "examples/tpb")))

(defsystem "cl-mpm/examples/lbar"
  :depends-on ("cl-mpm"
               "magicl"
               "lisp-stat"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/buoyancy"
               "cl-mpm/penalty"
               "cl-mpm/damage"
               "cl-mpm/plotter"
               "vgplot"
               "swank.live")
  :serial t
  :components ((:file "examples/lbar")))
