(in-package :cl-mpm)
(declaim #.cl-mpm/settings:*optimise-setting*)

(defgeneric special-g2p (mesh mp node svp grads)
  (:documentation "G2P behaviour for specific features")
  (:method (mesh mp node svp grads)))

(defmethod special-g2p (mesh (mp cl-mpm/particle::particle-thermal) node svp grads)
  (declare (ignore mesh))
  "Map grid to particle for one mp-node pair"
  (with-accessors ((node-temp cl-mpm/mesh:node-temperature)) node
    (with-accessors ((temp cl-mpm/particle::mp-temperature)) mp
      (progn
        (setf temp (+ temp (* node-temp svp)))))))

(defmethod special-g2p (mesh (mp cl-mpm/particle::particle-damage) node svp grads)
  (declare (ignore mesh))
  "Map grid to particle for one mp-node pair"
  (with-accessors ((node-damage cl-mpm/mesh::node-damage)
                   (node-temp cl-mpm/mesh::node-temp)) node
    (with-accessors ((temp cl-mpm/particle::mp-damage)) mp
      (progn
        (setf temp (+ temp (* node-temp svp)))
        ))))

(defgeneric reset-mps-g2p (mp)
  (:method (mp)))

(defmethod reset-mps-g2p ((mp cl-mpm/particle::particle-thermal))
  (with-accessors ((temp cl-mpm/particle::mp-temperature)) mp
      (setf temp 0d0)))

(declaim
 (notinline g2p-mp)
 (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle double-float) (values))
                g2p-mp))


(macrolet ((def-g2p-mp (name &body update)
             `(defun ,name (mesh mp dt damping)
                (declare (cl-mpm/mesh::mesh mesh)
                         (cl-mpm/particle:particle mp)
                         (double-float dt))
                "Map one MP from the grid"
                (with-accessors ((vel mp-velocity)
                                 (pos mp-position)
                                 (pos-trial cl-mpm/particle::mp-position-trial)
                                 (disp cl-mpm/particle::mp-displacement)
                                 (disp-inc cl-mpm/particle::mp-displacement-increment))
                    mp
                  (let* ((mapped-vel (cl-mpm/utils:vector-zeros))
                         (acc (cl-mpm/utils:vector-zeros)))
                    (cl-mpm/fastmaths:fast-zero disp-inc)
                    ;; Map variables
                    (iterate-over-neighbours
                     mesh mp
                     (lambda (mesh mp node svp grads fsvp fgrads)
                       (declare
                        (ignore mp mesh fsvp fgrads)
                        (cl-mpm/mesh::node node)
                        (cl-mpm/particle:particle mp)
                        (double-float svp))
                       (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
                                        (node-acc cl-mpm/mesh:node-acceleration)
                                        (node-disp cl-mpm/mesh::node-displacment)
                                        (node-scalar cl-mpm/mesh::node-boundary-scalar)
                                        (node-active cl-mpm/mesh:node-active)
                                        ) node
                         (declare (double-float node-scalar)
                                  (boolean node-active))
                         (when node-active
                           (cl-mpm/fastmaths::fast-fmacc mapped-vel node-vel svp)
                           (cl-mpm/fastmaths::fast-fmacc disp-inc node-disp svp)
                           (cl-mpm/fastmaths::fast-fmacc acc node-acc svp)
                           ;; (incf temp (* svp node-scalar))
                           ;;With special operations we want to include this operation
                           #+cl-mpm-special (special-g2p mesh mp node svp grads)
                           ))))
                    ;;Update particle
                    (progn
                      ;;Invalidate shapefunction/gradient cache
                      ;; (update-particle mesh mp dt)
                      ,@update
                      ))
                  ))

             ))

  (def-g2p-mp g2p-mp-flip
      (progn
        (cl-mpm/fastmaths:fast-.+ pos disp-inc pos-trial)
        (cl-mpm/fastmaths:fast-fmacc vel acc dt)))
  (def-g2p-mp g2p-mp-pic
      (progn
        (cl-mpm/utils::vector-copy-into mapped-vel vel)
        (cl-mpm/fastmaths:fast-.+ pos disp-inc pos-trial)))
  (def-g2p-mp g2p-mp-quasi-static
      (progn
        ;; (declare (ignore mapped-vel acc))
        (cl-mpm/fastmaths:fast-.+ pos disp-inc pos-trial)))
  (def-g2p-mp g2p-mp-blend
      (let ((pic-value 1d-3))
        (cl-mpm/fastmaths:fast-.+ pos disp-inc pos-trial)
        (cl-mpm/fastmaths:fast-.+
         (cl-mpm/fastmaths:fast-scale-vector
          ;; FLIP value
          (cl-mpm/fastmaths:fast-.+ vel (cl-mpm/fastmaths:fast-scale-vector acc dt))
          (- 1d0 pic-value))
         ;; PIC update
         (cl-mpm/fastmaths:fast-scale-vector mapped-vel pic-value)
         vel)))
  (def-g2p-mp g2p-mp-blend-2nd-order
      (let* ((pic-value (/ 1d-3 dt))
             (vel-inc
               (cl-mpm/fastmaths::fast-scale!
                (cl-mpm/fastmaths::fast-.--vector
                 acc
                 (cl-mpm/fastmaths:fast-scale!
                  (cl-mpm/fastmaths::fast-.--vector vel mapped-vel) pic-value))
                dt)))
        (cl-mpm/fastmaths::fast-.+-vector vel vel-inc vel)
        (let ((dx (cl-mpm/fastmaths::fast-.-
                   (cl-mpm/fastmaths:fast-scale-vector mapped-vel dt)
                   (cl-mpm/fastmaths:fast-scale-vector vel-inc (* 0.5d0 dt)))))
          (cl-mpm/fastmaths::fast-scale! disp-inc dt)))))

(declaim (notinline g2p))

(defun g2p (mesh mps dt damping &optional (update-type :FLIP))
  (declare (double-float dt damping))
  (ecase update-type
    (:QUASI-STATIC
     (iterate-over-mps
      mps
      (lambda (mp)
        (g2p-mp-quasi-static mesh mp dt damping))))
    (:FLIP
     (iterate-over-mps
      mps
      (lambda (mp)
        (g2p-mp-flip mesh mp dt damping))))
    (:PIC
     (iterate-over-mps
      mps
      (lambda (mp)
        (g2p-mp-pic mesh mp dt damping))))
    (:BLEND
     (iterate-over-mps
      mps
      (lambda (mp)
        (g2p-mp-blend mesh mp dt damping))))
    (:BLEND-2ND-ORDER
     (iterate-over-mps
      mps
      (lambda (mp)
        (g2p-mp-blend-2nd-order mesh mp dt damping))))))
