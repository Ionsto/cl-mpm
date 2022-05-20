(in-package :asdf-user)

(defsystem "cl-mpm"
  :class :package-inferred-system
  :depends-on ( "src/symbolic-derivation"
                "src/main"
                "l-math")
  :description ""
  :in-order-to ((test-op (load-op "test/all")))
  :perform (test-op (o c) (symbol-call :test/all :test-suite)))

(defsystem "cl-mpm/test"
  :depends-on ("test/all"))

(register-system-packages "cl-mpm/src/main" '(:main))
(register-system-packages "cl-mpm/test/all" '(:test/all))
