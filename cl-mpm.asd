(in-package :asdf-user)

(defsystem "symbolic-derivation"
  :description "Symbolic derivation library from #TODO"
  :author "Aleksander Ksiazek"
  :serial t
  :components ((:file "src/core/symbolic-derivation")))


(defsystem "cl-mpm/settings"
  :components ((:file "src/settings"))
  :description "Dummy package that can be used to prune magicl dependancies, i.e. lisp only implementation")

(defsystem "cl-mpm/magicl"
  :depends-on ("magicl"
               "magicl/ext-expokit")
  :description "Dummy package that can be used to prune magicl dependancies, i.e. lisp only implementation")

(defsystem "cl-mpm/utils"
  :depends-on ("cl-mpm/magicl"
               "array-operations"
               )
  :description "MPM utility functions definitions"
  :serial t
  :components ((:file "src/utils")))

(defsystem "cl-mpm/fastmaths"
  :depends-on ("cl-mpm/magicl"
               "cl-mpm/utils"
               :sb-simd)
  :description "MPM fast maths operations definitions"
  :serial t
  :components ((:file "src/fastmaths")))

(defsystem "cl-mpm/mesh"
  :depends-on ("cl-mpm/magicl"
               "cl-mpm/utils"
               "cl-mpm/fastmaths"
               "cl-mpm/shape-function")
  :description "MPM boundary conditions"
  :serial t
  :components ((:file "src/core/mesh")))

(defsystem "cl-mpm/forces"
  :depends-on ("cl-mpm/magicl"
               "cl-mpm/utils"
               "cl-mpm/particle"
               "cl-mpm/mesh")
  :description "MPM forces calculations"
  :serial t
  :components ((:file "src/core/forces")))

(defsystem "cl-mpm/bc"
  :depends-on ("cl-mpm/magicl"
               "cl-mpm/mesh"
               )
  :description "MPM boundary conditions"
  :serial t
  :components ((:file "src/bc/bc")))

(defsystem "cl-mpm/setup"
  :depends-on ("cl-mpm"
               "cl-mpm/mesh"
               )
  :description "MPM setup system"
  :serial t
  :components ((:file "src/setup")))


(defsystem "cl-mpm/constitutive"
  :depends-on ("cl-mpm/magicl"
               "cl-mpm/utils"
               "cl-mpm/fastmaths")
  :description "Various constitutive models"
  :serial t
  :components ((:file "src/constitutive")))
;
(defsystem "cl-mpm/particle"
  :depends-on ("cl-mpm/magicl"
               "cl-mpm/constitutive"
               ;; "cl-mpm/damage"
               "cl-mpm/mesh")
  :description "MPM particle definitions"
  :serial t
  :components ((:module "src"
                :serial t
                :components
                ((:file "particle")
                 (:module "models"
                  :serial t
                  :components
                  (;(:file "damage")
                   (:file "plastic")
                   (:file "concrete")
                   ;; (:file "limestone")
                   ;; (:file "ice")
                   ;(:file "chalk")
                   ))))))

(defsystem "cl-mpm/damage"
  :depends-on ("cl-mpm/magicl"
               "cl-mpm/utils"
               "cl-mpm/constitutive"
               "cl-mpm/output"
               "cl-mpm/particle")
  :description "MPM damage mechanics"
  :serial t
  :components ((:module "src"
                :serial t
                :components
                ((:module "damage"
                  :serial t
                  :components
                  ((:file "package")
                   (:file "softening")
                   (:file "criteria")
                   (:file "damage")
                   (:file "delay-damage")
                   ))
                 (:file "models/damage"))
               )))
(defsystem "cl-mpm/eigenerosion"
  :depends-on ("cl-mpm/magicl"
               "cl-mpm"
               "cl-mpm/utils"
               "cl-mpm/particle")
  :description "MPM eigenerosion damage mechanics"
  :serial t
  :components ((:file "src/eigenerosion")))

(defsystem "cl-mpm/output"
  :depends-on ("cl-mpm/magicl"
               "jonathan"
               "cl-mpm"
               "str")
  :description "MPM output helper functions"
  :serial t
  :components ((:file "src/output")))
;
(defsystem "cl-mpm/shape-function"
 :depends-on ("cl-mpm/magicl"
              "symbolic-derivation")
 :description "MPM shape function definitions"
 :serial t
 :components ((:file "src/core/shape-function")))

(defsystem "cl-mpm/dynamic-relaxation"
  :depends-on ("cl-mpm"
               "cl-mpm/penalty"
               "cl-mpm/mpi"
               "swank.live"
               "cl-mpm/utils")
  :description "MPM dynamic relaxation helpers - allows for modelling quasi-static problems"
  :serial t
  :components ((:file "src/dynamic-relaxation")))
(defsystem "cl-mpm/ext"
  :depends-on ("cl-mpm/magicl"
               "cffi"
               "cl-mpm/utils"
               "cl-mpm/constitutive"
               "array-operations")
  :description "External libraries"
  :serial t
  :components ((:file "src/cpp"))
  )
(defsystem "cl-mpm"
  ;; :class :package-inferred-system
  :depends-on ("cl-mpm/magicl"
               "cl-mpm/fastmaths"
               "cl-mpm/ext"
               "cl-mpm/utils"
               "cl-mpm/forces"
               "alexandria"
               "serapeum"
               "local-time"
               "array-operations"
               "lparallel"
               "symbolic-derivation"
               "cl-mpm/constitutive"
               "cl-mpm/particle"
               "cl-mpm/shape-function"
               "cl-mpm/bc"
               "cl-mpm/mesh"
               "cl-mpm/ext")
  :description "An explicit Material Point Method implementation"
  :in-order-to ((test-op (load-op "test/all")))
  :perform (test-op (o c) (symbol-call :test/all :test-suite))
  :serial t
  :components ((:file "src/cl-mpm")
               (:file "src/core/iterate")
               (:file "src/core/domain-update")
               (:file "src/core/stress-update")
               (:file "src/core")
               (:file "src/solvers/usf")
               (:file "src/solvers/usl")
               ))
(defsystem "cl-mpm/mpi"
  :depends-on ("cl-mpm"
               "lfarm-client"
               "lfarm-server"
               "lfarm-admin"
               "cl-mpm/fastmaths"
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
               "cl-mpm/setup"
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
               "cl-mpm/bc"
               "cl-mpm/mpi"
               )
  :description ""
  :components (
               (:file "src/bc/buoyancy")
               ))

(defsystem "cl-mpm/erosion"
  :depends-on ("cl-mpm"
               "cl-mpm/bc"
               "cl-mpm/damage")
  :description "Erosion specific BCs"
  :components ((:file "src/bc/erosion")))

(defsystem "cl-mpm/penalty"
  :depends-on ("cl-mpm"
               "cl-mpm/bc")
  :description ""
  :components (
               (:file "src/bc/penalty")
               ))

(defsystem "cl-mpm/ghost"
  :depends-on ("cl-mpm"
               "cl-mpm/bc")
  :description ""
  :components (
               (:file "src/ghost")
               ))

(defsystem "cl-mpm/test"
  :depends-on ("cl-mpm"))

(defsystem "cl-mpm/example"
  :depends-on ("cl-mpm/all")
  :serial t
  :components ((:file "examples/core")))



(defsystem "cl-mpm/example/bounce"
  :class :package-inferred-system
  :depends-on ("cl-mpm")
  :serial t
  :components ((:file "examples/bounce")))

(defsystem "cl-mpm/examples/column"
  :depends-on ("cl-mpm"
               "cl-mpm/all"
               "vgplot"
               "swank.live"
               "cl-mpm/magicl")
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
               "cl-mpm/magicl")
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
               "cl-mpm/magicl")
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
               "cl-mpm/dynamic-relaxation"
               "cl-mpm/plotter"
               "cl-mpm/models/ice"
               "lfarm-client"
               "vgplot"
               "swank.live"
               "lisp-stat"
               "cl-mpm/magicl")
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
               "cl-mpm/magicl")
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
               "cl-mpm/magicl")
  :serial t
  :components ((:file "examples/single-crack")))

(defsystem "cl-mpm/examples/float"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/buoyancy"
               "cl-mpm/plotter"
               "cl-mpm/models/ice"
               "vgplot"
               "swank.live"
               "cl-mpm/magicl")
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
               "cl-mpm/magicl")
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
               "cl-mpm/magicl")
  :serial t
  :components ((:file "examples/beam")))

(defsystem "cl-mpm/examples/shear"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "vgplot"
               "swank.live"
               "cl-mpm/magicl")
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
               "cl-mpm/magicl")
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
               "cl-mpm/magicl")
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
               "cl-mpm/magicl")
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
               "cl-mpm/magicl")
  :serial t
  :components ((:file "examples/sliding")))

(defsystem "cl-mpm/examples/collapse"
  :depends-on ("cl-mpm"
               "cl-mpm/all"
               "vgplot"
               "cl-mpm/models/visco"
               "swank.live")
  :serial t
  :components ((:file "examples/collapse")))

(defsystem "cl-mpm/examples/rotate"
  :depends-on ("cl-mpm/all"
               "vgplot"
               "swank.live")
  :serial t
  :components ((:file "examples/rotate")))



(defsystem "cl-mpm/examples/penalty-friction"
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
  :components ((:file "examples/penalty-friction")))





(defsystem "cl-mpm/examples/chalk"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/buoyancy"
               "cl-mpm/penalty"
               "cl-mpm/damage"
               "cl-mpm/plotter"
               "cl-mpm/mpi"
               "parse-float"
               "vgplot"
               "swank.live")
  :serial t
  :components ((:file "examples/chalk")))

(defsystem "cl-mpm/examples/joss"
  :depends-on ("cl-mpm"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/buoyancy"
               "cl-mpm/erosion"
               "cl-mpm/penalty"
               "cl-mpm/damage"
               "cl-mpm/plotter"
               "cl-mpm/mpi"
               "cl-mpm/models/chalk"
               "cl-mpm/dynamic-relaxation"
               "cl-mpm/examples/chalk")
  :serial t
  :components ((:file "examples/joss")))


(defsystem "cl-mpm/examples/tpb"
  :depends-on ("cl-mpm"
               "cl-mpm/magicl"
               "lisp-stat"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/buoyancy"
               "cl-mpm/penalty"
               "cl-mpm/damage"
               "cl-mpm/plotter"
               "cl-mpm/dynamic-relaxation"
               "cl-mpm/models/limestone"
               "vgplot"
               "swank.live")
  :serial t
  :components ((:file "examples/tpb")))

(defsystem "cl-mpm/examples/lbar"
  :depends-on ("cl-mpm"
               "cl-mpm/magicl"
               "lisp-stat"
               "cl-mpm/setup"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/buoyancy"
               "cl-mpm/penalty"
               "cl-mpm/damage"
               "cl-mpm/plotter"
               "cl-mpm/dynamic-relaxation"
               "cl-mpm/models/limestone"
               "parse-float"
               ;; "cl-mpm/models/limestone"
               "vgplot"
               "swank.live")
  :serial t
  :components ((:file "examples/lbar")))
(defsystem "cl-mpm/all"
  :depends-on ("cl-mpm"
               "cl-mpm/forces"
               "cl-mpm/fastmaths"
               "cl-mpm/setup"
               "cl-mpm/bc"
               "cl-mpm/constitutive"
               "cl-mpm/particle"
               "cl-mpm/output"
               "cl-mpm/buoyancy"
               "cl-mpm/penalty"
               "cl-mpm/damage"
               "cl-mpm/plotter"
               "cl-mpm/mpi"
               "cl-mpm/dynamic-relaxation"
               "cl-mpm/ghost"))


(defsystem "cl-mpm/examples/uniaxial"
  :depends-on ("cl-mpm/examples/tpb"
               "cl-mpm/models/chalk"
               )
  :serial t
  :components ((:file "examples/uniaxial")))
(defsystem "cl-mpm/examples/brazilian"
  :depends-on ("cl-mpm/examples/tpb")
  :serial t
  :components ((:file "examples/brazilian")))

(defsystem "cl-mpm/examples/shear-box"
  :depends-on ("cl-mpm/all"
               "parse-float"
               "lisp-stat"
               "cl-mpm/models/chalk"
               "vgplot"
               "swank.live")
  :serial t
  :components ((:file "examples/shear-box")))

(defsystem "cl-mpm/examples/erode"
  :depends-on ("cl-mpm/example")
  :components ((:file "examples/erode")))


(defsystem "cl-mpm/examples/ice-buoyancy"
  :depends-on ("cl-mpm/example"
               "cl-mpm/erosion"
               "cl-mpm/models/chalk"
               "cl-mpm/models/visco"
               "cl-mpm/models/ice"
               )
  :components ((:file "examples/ice-buoyancy")))

;; (defsystem "cl-mpm/models/all"
;;   :depends-on ("cl-mpm/particle"
;;                "cl-mpm/models/damage"
;;                "cl-mpm/models/ice"
;;                "cl-mpm/models/chalk")
;;   :description "MPM plastic particle definitions"
;;   :serial t)

;; (defsystem "cl-mpm/models/plastic"
;;   :depends-on ("cl-mpm/particle")
;;   :description "MPM plastic particle definitions"
;;   :serial t
;;   :components ((:file "src/models/plastic")))

;; (defsystem "cl-mpm/models/damage"
;;   :depends-on ("cl-mpm/particle"
;;                "cl-mpm/damage")
;;   :description "MPM plastic particle definitions"
;;   :serial t
;;   :components ((:file "src/models/damage")))

;; (defsystem "cl-mpm/models/ice"
;;   :depends-on ("cl-mpm/particle"
;;                "cl-mpm/models/plastic"
;;                "cl-mpm/models/damage"
;;                )
;;   :description "MPM ice particle definitions"
;;   :serial t
;;   :components ((:file "src/models/ice")))

(defsystem "cl-mpm/models/chalk"
  :depends-on ("cl-mpm/particle"
               ;; "cl-mpm/models/plastic"
               ;; "cl-mpm/models/damage"
               "cl-mpm/damage"
               )
  :description "MPM chalk particle definitions"
  :serial t
  :components ((:file "src/models/chalk")))


(defsystem "cl-mpm/models/ice"
  :depends-on ("cl-mpm/particle"
               "cl-mpm/models/chalk"
               "cl-mpm/damage"
               )
  :description "MPM ice particle definitions"
  :serial t
  :components ((:file "src/models/ice")))

(defsystem "cl-mpm/models/limestone"
  :depends-on ("cl-mpm/particle"
               "cl-mpm/damage")
  :description "MPM limestone definitions"
  :serial t
  :components ((:file "src/models/limestone")))


(defsystem "cl-mpm/models/visco"
  :depends-on ("cl-mpm/particle"
               "cl-mpm/ext"
               )
  :description "MPM limestone definitions"
  :serial t
  :components ((:file "src/models/visco")))

(defsystem "cl-mpm/implicit"
  :depends-on ("cl-mpm"
               "cl-mpm/utils"
               "cl-mpm/constitutive"
               "cl-mpm/output"
               "cl-mpm/particle")
  :description "Implicit quasi-static implementation"
  :serial t
  :components ((:file "src/implicit")))

(defsystem "cl-mpm/examples/ice-visco"
  :depends-on ("cl-mpm/example"
               "cl-mpm/erosion"
               "cl-mpm/models/chalk"
               "cl-mpm/models/visco"
               )
  :components ((:file "examples/ice-visco")))

