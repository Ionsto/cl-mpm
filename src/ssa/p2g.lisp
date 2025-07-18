
(in-package :cl-mpm/ssa)

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
      (cl-mpm::iterate-over-neighbours
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

(defun p2g (mesh mps)
  "Map particle momentum to the grid"
  (declare (type (vector cl-mpm/particle:particle) mps) (cl-mpm/mesh::mesh mesh))
  (cl-mpm::iterate-over-mps
   mps
   (lambda (mp)
     (p2g-mp mesh mp))))


(defun det-ext-force (mp node svp grads node-ext-force)
  (with-accessors ((mass cl-mpm/particle::mp-mass))
      mp
      (let* (;; (dsvp (cl-mpm/shape-function::assemble-dsvp-3d grads))
             ;; (h-v (cl-mpm/particle::mp-height mp))
             (gravity 9.8d0)
             (gravity 9.8d0)
             ;; (h-v (* gravity mass))
             (h-v (/ mass volume))
             (height-vec (cl-mpm/utils:vector-from-list (list h-v h-v h-v)))
             )
        ;; (@-dsvp-vec-simd dsvp stress volume res)
        (cl-mpm/fastmaths::fast-.+
         ;; (magicl:@ dsvp height-vec)
         (cl-mpm/fastmaths:fast-.* (cl-mpm/utils:vector-from-list grads) height-vec)
         node-ext-force node-ext-force))))

(defun p2g-force-mp (mesh mp)
  "Map particle forces to the grid for one mp"
  (declare (cl-mpm/mesh::mesh mesh)
    (cl-mpm/particle:particle mp))
  (cl-mpm::iterate-over-neighbours
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
           (det-ext-force mp node svp grads node-ext-force)
           (cl-mpm::det-int-force-unrolled-2d mp grads node-int-force)
           ;; (det-ext-force-2d mp node svp node-ext-force)
           ;; (det-int-force-unrolled-2d mp grads node-int-force)
           )))))
  (values))
(defun p2g-force (mesh mps)
  "Map particle forces to the grid"
  (declare (type (array cl-mpm/particle:particle) mps) (cl-mpm/mesh::mesh mesh))
  (cl-mpm::iterate-over-mps
   mps
   (lambda (mp)
     (p2g-force-mp mesh mp))))
