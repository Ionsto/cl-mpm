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
   (eroded-volume-n
    :initform 0d0
    :accessor mp-eroded-volume-n)
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
                            (enable t)
                            )
  (with-accessors ((mesh cl-mpm:sim-mesh))
      sim
    (make-instance 'bc-erode
                   :index nil
                   :damage-rate rate
                   :damage-volume nil
                   :clip-func clip-func
                   :scalar-func scalar-func
                   :enable enable
                   :sim sim)))

(defgeneric mp-erosion-enhancment (mp))
(defmethod mp-erosion-enhancment ((mp cl-mpm/particle::particle))
  1d0)

(defun apply-erosion (bc mesh dt)
  (with-accessors ((sim cl-mpm/buoyancy::bc-buoyancy-sim)
                   (clip-func cl-mpm/buoyancy::bc-buoyancy-clip-func)
                   (scalar-func cl-mpm/buoyancy::bc-scalar-func)
                   (datum cl-mpm/buoyancy::bc-buoyancy-datum)
                   (rate bc-water-damage-damage-rate))
      bc
    (when (cl-mpm/buoyancy::bc-enable bc)
      (cl-mpm/buoyancy::apply-scalar
       sim
       scalar-func
       (lambda (pos)
         (and
          (funcall clip-func pos)))
       datum
       :damage-volume nil;damage-volume
       )
      (cl-mpm:iterate-over-mps
       (cl-mpm:sim-mps sim)
       (lambda (mp)
         (with-accessors ((volume cl-mpm/particle::mp-volume)
                          (mass cl-mpm/particle::mp-mass)
                          (erode cl-mpm/particle::mp-eroded-volume)
                          (erode-n cl-mpm/particle::mp-eroded-volume-n)
                          (erosion-modulus cl-mpm/particle::mp-erosion-modulus))
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
                                 (node-active cl-mpm/mesh::node-active))
                    node
                  (when (and node-active
                             boundary-node)
                    (incf weathering (* svp (/ node-boundary-scalar node-volume)))))))

             (setf weathering (min weathering 0d0))
             (setf weathering (* rate (min weathering 0d0) volume))
             (setf weathering (* weathering (mp-erosion-enhancment mp)))
             ;; (setf weathering (- (sqrt (abs weathering))))
             ;(setf weathering (* weathering (+ 1d0 (* 8 (cl-mpm/particle:mp-damage mp)))))
             ;(setf weathering (* weathering (+ 1d0 (* 10 (cl-mpm/particle::mp-strain-plastic-vm mp)))))
             ;; (setf weathering (* weathering (+ 1d0 (* 10 (cl-mpm/particle::mp-strain-plastic-vm mp)))))

             (setf weathering (* (/ (- weathering) erosion-modulus) dt))
             (setf (cl-mpm/particle::mp-boundary mp) weathering)
             (setf erode (+ erode-n weathering))
             ;; (pprint erode)
             ;; (let ((density (/ mass volume)))
             ;;   (setf
             ;;    mass
             ;;    (max
             ;;     0d0
             ;;     (-
             ;;      mass
             ;;      weathering
             ;;      )))
             ;;   (setf volume (/ mass density))
               ))))

      ;; (cl-mpm::remove-mps-func sim (lambda (mp) (>= (cl-mpm/particle::mp-eroded-volume mp)
      ;;                                               (cl-mpm/particle::mp-mass mp))))
      ))
  )

(defmethod cl-mpm/bc::apply-bc ((bc bc-erode) node mesh dt)
  (apply-erosion bc mesh dt))

(defmethod cl-mpm::new-loadstep :after ((sim cl-mpm::mpm-sim))
  (cl-mpm::remove-mps-func
   sim
   (lambda (mp)
     (when (typep mp 'cl-mpm/particle::particle-erosion)
       (>= (cl-mpm/particle::mp-eroded-volume mp)
           (cl-mpm/particle::mp-mass mp))))))

(defmethod cl-mpm/particle::reset-loadstep-mp ((mp cl-mpm/particle::particle-erosion))
  (with-accessors ((k    cl-mpm/particle::mp-eroded-volume)
                   (k-n    cl-mpm/particle::mp-eroded-volume-n))
      mp
    (setf k k-n)
    (call-next-method)))

(defmethod cl-mpm/particle::new-loadstep-mp ((mp cl-mpm/particle::particle-erosion))
  (with-accessors ((k    cl-mpm/particle::mp-eroded-volume)
                   (k-n    cl-mpm/particle::mp-eroded-volume-n))
      mp
    (setf k-n k)
    (call-next-method)))
