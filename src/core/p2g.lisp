(in-package :cl-mpm)
(declaim #.cl-mpm/settings:*optimise-setting*)

(declaim
 (inline p2g-mp)
 (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle) (values))
                p2g-mp))
(defun p2g-mp (mesh mp)
  "P2G transfer for one MP"
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (with-accessors ((mp-vel  cl-mpm/particle:mp-velocity)
                   (mp-mass cl-mpm/particle:mp-mass)
                   (mp-volume cl-mpm/particle:mp-volume)
                   (mp-pmod cl-mpm/particle::mp-p-modulus)
                   (mp-disp cl-mpm/particle::mp-displacement-increment)
                   (mp-damage cl-mpm/particle::mp-damage))
      mp
    (let ((mp-mass mp-mass)
          (mp-vel mp-vel)
          (mp-volume mp-volume)
          (mp-pmod mp-pmod)
          (mp-damage mp-damage))
      (declare (type double-float mp-mass mp-volume))
      (iterate-over-neighbours
       mesh mp
       (lambda (mesh mp node svp grads fsvp fgrads)
         (declare
          (cl-mpm/particle:particle mp)
          (cl-mpm/mesh::node node)
          (double-float svp)
          (ignore mesh))
         (with-accessors ((node-vel   cl-mpm/mesh:node-velocity)
                          (node-active  cl-mpm/mesh::node-active)
                          (node-mass  cl-mpm/mesh:node-mass)
                          (node-volume  cl-mpm/mesh::node-volume)
                          (node-volume-true  cl-mpm/mesh::node-volume-true)
                          (node-svp-sum  cl-mpm/mesh::node-svp-sum)
                          (node-force cl-mpm/mesh:node-force)
                          (node-p-wave cl-mpm/mesh::node-pwave)
                          (node-damage cl-mpm/mesh::node-damage)
                          (node-disp cl-mpm/mesh::node-displacment)
                          (node-lock  cl-mpm/mesh:node-lock))
             node
           (declare (type double-float node-mass node-volume mp-volume mp-pmod mp-damage node-svp-sum svp node-p-wave node-damage)
                    (type sb-thread:mutex node-lock))
           (sb-thread:with-mutex (node-lock)
             (setf node-active t)
             (incf node-mass (* mp-mass svp))
             (incf node-volume (* mp-volume svp))
             (incf node-p-wave (* mp-pmod svp))
             (incf node-damage (* mp-damage svp))
             (incf node-svp-sum svp)
             (cl-mpm/fastmaths::fast-fmacc node-vel mp-vel (* mp-mass svp))))))))
  (values))


(declaim (inline p2g-force-mp)
         (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle) (values)) p2g-force-mp)
         )
(defun p2g-force-mp (mesh mp)
  "Map particle forces to the grid for one mp"
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (iterate-over-neighbours
     mesh mp
     (lambda (mesh mp node svp grads fsvp fgrads)
       (declare
        (cl-mpm/particle:particle mp)
        (cl-mpm/mesh::node node)
        (double-float svp)
        (ignore mesh fsvp fgrads))
       (with-accessors ((node-active  cl-mpm/mesh:node-active)
                        (node-int-force cl-mpm/mesh::node-internal-force)
                        (node-ext-force cl-mpm/mesh::node-external-force)
                        (node-lock  cl-mpm/mesh:node-lock))
           node
         (declare (boolean node-active)
                  (sb-thread:mutex node-lock)
                  (magicl:matrix/double-float node-int-force node-ext-force))
         (when node-active
           (sb-thread:with-mutex (node-lock)
             (det-ext-force mp node svp node-ext-force)
             (det-int-force-unrolled mp grads node-int-force))))))
  (values))

(declaim (notinline p2g-force-mp-2d)
         (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle) (values)) p2g-force-mp-2d)
         )
;; (defun p2g-force-mp-2d (mesh mp)
;;   "Map particle forces to the grid for one mp"
;;   (declare (cl-mpm/mesh::mesh mesh)
;;            (cl-mpm/particle:particle mp))
;;   (let (;; (df-inv (cl-mpm/fastmaths::fast-inv-3x3 (cl-mpm/particle::mp-deformation-gradient-increment mp)))
;;         (df-inv (magicl:inv (cl-mpm/particle::mp-deformation-gradient-increment mp)))
;;         (dsvp (cl-mpm/utils::dsvp-3d-zeros))
;;         (dsvp-push (cl-mpm/utils::dsvp-3d-zeros))
;;         )
;;     (iterate-over-neighbours
;;      mesh mp
;;      (lambda (mesh mp node svp grads fsvp fgrads)
;;        (declare
;;         (cl-mpm/particle:particle mp)
;;         (cl-mpm/mesh::node node)
;;         (double-float svp))
;;        (with-accessors ((node-active  cl-mpm/mesh:node-active)
;;                         (node-int-force cl-mpm/mesh::node-internal-force)
;;                         (node-ext-force cl-mpm/mesh::node-external-force)
;;                         (node-lock  cl-mpm/mesh:node-lock)) node
;;          (declare (boolean node-active)
;;                   (sb-thread:mutex node-lock)
;;                   (magicl:matrix/double-float node-int-force node-ext-force))
;;          (when node-active
;;            (cl-mpm/shape-function::assemble-dsvp-3d-prealloc grads dsvp)
;;            (magicl:mult dsvp df-inv :target dsvp-push)
;;            (sb-thread:with-mutex (node-lock)
;;              (det-ext-force-2d mp node svp node-ext-force)
;;              (det-int-force mp dsvp-push node-int-force)))))))
;;   (values))

(defun p2g-force-mp-2d (mesh mp)
  "Map particle forces to the grid for one mp"
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (let ((df-inv (magicl:inv (cl-mpm/particle::mp-deformation-gradient-increment mp))))
    (iterate-over-neighbours
     mesh mp
     (lambda (mesh mp node svp grads fsvp fgrads)
       (declare
        (cl-mpm/particle:particle mp)
        (cl-mpm/mesh::node node)
        (double-float svp))
       (with-accessors ((node-active  cl-mpm/mesh:node-active)
                        (node-int-force cl-mpm/mesh::node-internal-force)
                        (node-ext-force cl-mpm/mesh::node-external-force)
                        (node-lock  cl-mpm/mesh:node-lock)) node
         (declare (boolean node-active)
                  (sb-thread:mutex node-lock)
                  (magicl:matrix/double-float node-int-force node-ext-force))
         (let ((grads (cl-mpm::gradient-push-forwards-cached grads df-inv)))
           (when node-active
             (sb-thread:with-mutex (node-lock)
               (det-ext-force-2d mp node svp node-ext-force)
               (det-int-force-unrolled-2d mp grads node-int-force))))))))
  (values))

(defun p2g-force-mp-2d-ss (mesh mp)
  "Map particle forces to the grid for one mp"
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (iterate-over-neighbours
   mesh mp
   (lambda (mesh mp node svp grads fsvp fgrads)
     (declare
      (cl-mpm/particle:particle mp)
      (cl-mpm/mesh::node node)
      (double-float svp))
     (with-accessors ((node-active  cl-mpm/mesh:node-active)
                      (node-int-force cl-mpm/mesh::node-internal-force)
                      (node-ext-force cl-mpm/mesh::node-external-force)
                      (node-lock  cl-mpm/mesh:node-lock)) node
       (declare (boolean node-active)
                (sb-thread:mutex node-lock)
                (magicl:matrix/double-float node-int-force node-ext-force))
       (when node-active
         (sb-thread:with-mutex (node-lock)
           (det-ext-force-2d mp node svp node-ext-force)
           (det-int-force-unrolled-2d mp grads node-int-force))))))
  (values))

(defgeneric special-p2g (mp node svp dsvp)
  (:documentation "P2G behaviour for specific features")
  (:method (mp node svp dsvp)))

(defmethod special-p2g ((mp cl-mpm/particle::particle-thermal) node svp dsvp)
  (with-accessors ((node-mass  cl-mpm/mesh:node-mass)
                   (node-temp  cl-mpm/mesh:node-temperature)
                   (node-dtemp  cl-mpm/mesh:node-dtemp)
                   (node-volume  cl-mpm/mesh::node-volume)
                   (node-lock  cl-mpm/mesh:node-lock)) node
    (with-accessors ((mp-mass cl-mpm/particle:mp-mass)
                     (mp-temp cl-mpm/particle::mp-temperature)
                     (mp-heat-capaciy cl-mpm/particle::mp-heat-capacity)
                     (mp-thermal-conductivity cl-mpm/particle::mp-thermal-conductivity)
                     (mp-volume cl-mpm/particle:mp-volume)
                     ) mp
            (progn
              (let* ((weighted-temp (* mp-temp mp-mass svp))
                     (weighted-dtemp (* (/ mp-volume
                                           (* mp-mass
                                              mp-heat-capaciy))
                                        mp-thermal-conductivity
                                        mp-temp
                                        (magicl::sum
                                         (cl-mpm/fastmaths::fast-.* dsvp dsvp)))))
                (sb-thread:with-mutex (node-lock)
                    (setf node-temp
                          (+ node-temp weighted-temp))
                    (setf node-dtemp
                          (+ node-dtemp weighted-dtemp))))))))


(declaim (notinline p2g))
(defun p2g (mesh mps)
  "Map particle momentum to the grid"
  (declare (type (vector cl-mpm/particle:particle) mps) (cl-mpm/mesh::mesh mesh))
  (iterate-over-mps
   mps
   (lambda (mp)
     (p2g-mp mesh mp))))


(declaim (inline p2g-force))
(defun p2g-force (mesh mps)
  "Map particle forces to the grid"
  (declare (type (array cl-mpm/particle:particle) mps) (cl-mpm/mesh::mesh mesh))
  (if (= (the fixnum (cl-mpm/mesh:mesh-nd mesh)) 2)
      (iterate-over-mps
       mps
       (lambda (mp)
         (p2g-force-mp-2d mesh mp)))
      (iterate-over-mps
       mps
       (lambda (mp)
         (p2g-force-mp mesh mp)))))
