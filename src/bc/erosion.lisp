(defpackage :cl-mpm/erosion
  (:use :cl
   :cl-mpm
        :cl-mpm/utils
   :cl-mpm/fastmaths)
  (:import-from
   :magicl tref .+ .-)
  (:export
   ))

;; (declaim (optimize (debug 3) (safety 3) (speed 0)))
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(in-package :cl-mpm/particle)
(defclass particle-erosion (particle-elastic-damage)
  ((eroded-volume
    :initform 0d0
    :accessor mp-eroded-volume)
   (erosion-modulus
    :accessor mp-erosion-modulus
    :initform 1d0
    :initarg :erosion-modulus)))
(in-package :cl-mpm/erosion)

(defclass bc-erode (cl-mpm/buoyancy::bc-scalar)
  ((damage-rate
    :initform 1d0
    :initarg :damage-rate
    :accessor bc-water-damage-damage-rate)))

(defun make-bc-erode (sim &key (rate 1d0)
                            (clip-func (lambda (pos) t))
                            (scalar-func (lambda (pos) 1d0))
                            )
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (make-instance 'bc-erode
                   :index nil
                   :damage-rate rate
                   :damage-volume nil
                   :clip-func clip-func
                   :scalar-func scalar-func
                   :sim sim)))
(defmethod cl-mpm/bc::apply-bc ((bc bc-erode) node mesh dt)
  (call-next-method)
  (with-accessors ((sim cl-mpm/buoyancy::bc-buoyancy-sim)
                   (rate bc-water-damage-damage-rate))
      bc
    (when (cl-mpm/buoyancy::bc-enable bc)
      (cl-mpm:iterate-over-mps
       (cl-mpm:sim-mps sim)
       (lambda (mp)
         (with-accessors ((volume cl-mpm/particle::mp-volume)
                          (mass cl-mpm/particle::mp-mass)
                          (erode cl-mpm/particle::mp-eroded-volume)
                          (erosion-modulus cl-mpm/particle::mp-erosion-modulus)
                          )
             mp
           (declare (double-float volume mass erode erosion-modulus))
           (let ((weathering 0d0))
             (cl-mpm:iterate-over-neighbours
              mesh
              mp
              (lambda (mesh mp node svp grads fsvp fgrads)
                (with-accessors ((node-boundary-scalar cl-mpm/mesh::node-boundary-scalar)
                                 (node-volume cl-mpm/mesh::node-volume)
                                 (boundary-node cl-mpm/mesh::node-boundary-node)
                                 (node-active cl-mpm/mesh::node-active)
                                 )
                    node
                  (when (and node-active
                             boundary-node)
                    (incf weathering (* svp (/ node-boundary-scalar node-volume)))))))

             (setf weathering (min weathering 0d0))
             (setf weathering (* rate (min weathering 0d0) volume))
             ;; (setf weathering (- (sqrt (abs weathering))))
             (setf weathering (* weathering (+ 1d0 (* 8 (cl-mpm/particle:mp-damage mp)))))
             (setf (cl-mpm/particle::mp-boundary mp) weathering)
             (incf erode (* (/ (- weathering) erosion-modulus) dt))
               ;; (let ((density (/ mass volume)))
               ;;   (setf
               ;;    mass
               ;;    (max
               ;;     0d0
               ;;     (-
               ;;      mass
               ;;      (abs (*
               ;;            (bc-water-damage-damage-rate bc)
               ;;            weathering dt)))))
               ;;   ;; (setf mass (* density volume))
               ;;   )
               ))))

      (cl-mpm::remove-mps-func sim (lambda (mp) (>= (cl-mpm/particle::mp-eroded-volume mp)
                                                    (cl-mpm/particle::mp-mass mp)))))))
