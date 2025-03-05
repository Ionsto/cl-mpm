(defpackage :cl-mpm/damage
  (:use :cl
   :cl-mpm/utils)
  (:export
   #:mpm-sim-damage
   #:calculate-damage
   #:update-damage
   #:post-stress-step
   #:update-damage
   #:damage-model-calculate-y))
