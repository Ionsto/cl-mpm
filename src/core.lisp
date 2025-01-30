;    #:make-shape-function
(in-package :cl-mpm)

;; (declaim (optimize (debug 3) (safety 3) (speed 0)))
(declaim (optimize (debug 0) (safety 0) (speed 3)))
(eval-when (:compile-toplevel :load-toplevel :execute)
  #+cl-mpm-special (print "Compiled with special hooks")
  #-cl-mpm-special (print "Compiled without special hooks")
  (format t "Compiled with optimize ~A~%" (sb-ext:restrict-compiler-policy)))

(defun FLIP-status ()
  #+cl-mpm-pic (print "Compiled with PIC")
  #-cl-mpm-pic (print "Compiled with FLIP")
  )

(defun check-mps (sim)
  "Function to check that stresses and positions are sane, deleting mps moving very fast"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
  (with-accessors ((h cl-mpm/mesh::mesh-resolution))
      mesh
    ;; loop for mp across mps do
    ;; (cl-mpm::iterate-over-mps
    ;;  mps
    ;;  (lambda (mp)
    ;;    (progn
    ;;      (with-accessors ((pos cl-mpm/particle::mp-position)
    ;;                       (vel cl-mpm/particle::mp-velocity)
    ;;                       (domain cl-mpm/particle::mp-domain-size)
    ;;                       ) mp
    ;;        (loop for i from 0 to 2
    ;;              do (progn
    ;;                   (when (sb-ext:float-nan-p (varef pos i)) 
    ;;                     (pprint mp)
    ;;                     (error "NaN location found for ~A" mp))
    ;;                   (when (equal (abs (varef vel i)) #.sb-ext:double-float-positive-infinity)
    ;;                     (pprint mp)
    ;;                     (error "Infinite velocity found"))
    ;;                   (when (> (abs (varef vel i)) 1e10)
    ;;                     (pprint mp)
    ;;                     (error "High velocity found"))))))))
    (remove-mps-func
       sim
       (lambda (mp)
         (with-accessors ((damage cl-mpm/particle:mp-damage)
                          (def cl-mpm/particle::mp-deformation-gradient))
             mp
           (or
            (gimp-removal-criteria mp h)
            ;; (cl-mpm/mesh::in-bounds-array mesh (magicl::matrix/double-float-storage
            ;;                                     (cl-mpm/particle::mp-position mp)))
            ))))
    )))
;; (sb-ext:restrict-compiler-policy 'debug  0 0)
;; (sb-ext:restrict-compiler-policy 'speed  3 3)
(defun check-single-mps (sim)
  "Function to check and remove single material points"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
  (with-accessors ((h cl-mpm/mesh::mesh-resolution))
      mesh
    (iterate-over-mps
     mps
     (lambda (mp)
       (setf (cl-mpm/particle::mp-single-particle mp)
             (single-particle-criteria mesh mp))))
      (remove-mps-func
       sim
       #'cl-mpm/particle::mp-single-particle))))


(defgeneric update-sim (sim)
  (:documentation "Update an mpm simulation by one timestep"))

;; (defmethod update-sim ((sim mpm-sim))
;;   (declare (cl-mpm::mpm-sim sim))
;;   (with-slots ((mesh mesh)
;;                (mps mps)
;;                (bcs bcs)
;;                (bcs-force bcs-force)
;;                (dt dt)
;;                (mass-filter mass-filter)
;;                (split allow-mp-split)
;;                (enable-damage enable-damage)
;;                (nonlocal-damage nonlocal-damage)
;;                (remove-damage allow-mp-damage-removal)
;;                (fbar enable-fbar)
;;                (bcs-force-list bcs-force-list)
;;                (vel-algo velocity-algorithm)
;;                (time time))
;;                 sim
;;     (declare (type double-float mass-filter))
;;                 (progn

;;                   (reset-grid mesh)
;;                   (p2g mesh mps)
;;                   (when (> mass-filter 0d0)
;;                     (filter-grid mesh (sim-mass-filter sim)))
;;                   (update-node-kinematics mesh dt )
;;                   (apply-bcs mesh bcs dt)
;;                   (update-stress mesh mps dt fbar)
;;                   ;; ;; Map forces onto nodes
;;                   (p2g-force mesh mps)
;;                   (loop for bcs-f in bcs-force-list
;;                         do (apply-bcs mesh bcs-f dt))
;;                   (update-node-forces sim)
;;                   ;; Reapply velocity BCs
;;                   (apply-bcs mesh bcs dt)
;;                   ;; Also updates mps inline
;;                   (g2p mesh mps dt vel-algo)
;;                   (when split
;;                     (split-mps sim))
;;                   (check-mps sim)
;;                   (incf time dt))))
(defmethod update-sim ((sim mpm-sim-usf))
  "Update stress first algorithm"
  (declare (cl-mpm::mpm-sim-usf sim))
  (with-slots ((mesh mesh)
               (mps mps)
               (bcs bcs)
               (bcs-force bcs-force)
               (bcs-force-list bcs-force-list)
               (dt dt)
               (mass-filter mass-filter)
               (split allow-mp-split)
               (enable-damage enable-damage)
               (nonlocal-damage nonlocal-damage)
               (remove-damage allow-mp-damage-removal)
               (fbar enable-fbar)
               (time time)
               (vel-algo velocity-algorithm)
               )
                sim
    (declare (double-float mass-filter dt time))
    (progn
      (reset-grid mesh)
      (when (> (length mps) 0)
        (p2g mesh mps)
        (when (> mass-filter 0d0)
          (filter-grid mesh (sim-mass-filter sim)))
        (update-node-kinematics mesh dt )
        (apply-bcs mesh bcs dt)
        (update-stress mesh mps dt fbar)
        ;; Map forces onto nodes
        (p2g-force mesh mps)
        (loop for bcs-f in bcs-force-list
              do (apply-bcs mesh bcs-f dt))
        (update-node-forces sim)
        ;; Reapply velocity BCs
        (apply-bcs mesh bcs dt)
        ;; Also updates mps inline
        (g2p mesh mps dt vel-algo)
        (when split
          (split-mps sim)))
      ;; (check-mps sim)
      (incf time dt))))

(defmethod update-sim ((sim mpm-sim-usl))
  "Update stress last algorithm"
  (declare (cl-mpm::mpm-sim sim))
  (with-slots ((mesh mesh)
               (mps mps)
               (bcs bcs)
               (bcs-force bcs-force)
               (dt dt)
               (mass-filter mass-filter)
               (split allow-mp-split)
               (enable-damage enable-damage)
               (nonlocal-damage nonlocal-damage)
               (remove-damage allow-mp-damage-removal)
               (fbar enable-fbar)
               (bcs-force-list bcs-force-list)
               (vel-algo velocity-algorithm)
               )
                sim
    (declare (type double-float mass-filter))
                (progn
                    (reset-grid mesh)
                    ;; Map momentum to grid
                    (p2g mesh mps)
                    ;;Reset nodes below our mass-filter
                    (when (> mass-filter 0d0)
                      (filter-grid mesh (sim-mass-filter sim)))
                    ;;Turn momentum into velocity
                    (update-node-kinematics mesh dt)
                    (p2g-force mesh mps)
                    (loop for bcs-f in bcs-force-list
                          do (apply-bcs mesh bcs-f dt))
                    ;; (apply-bcs mesh bcs-force dt)
                    ;;Update our nodes after force mapping
                    (update-node-forces sim)
                    ;;Apply velocity bcs
                    (apply-bcs mesh bcs dt)
                    ;;Grid to particle mapping
                    (g2p mesh mps dt vel-algo)

                    ;;2nd round of momentum mapping
                    (reset-grid-velocity mesh)
                    (p2g mesh mps)
                    (when (> mass-filter 0d0)
                     (filter-grid-velocity mesh (sim-mass-filter sim)))
                    (update-node-kinematics mesh dt)
                    (apply-bcs mesh bcs dt)
                    ;;Update stress last
                    (update-stress mesh mps dt fbar)

                    (when remove-damage
                      (remove-material-damaged sim))
                    (when split
                      (split-mps sim))
                    (check-mps sim)
                    )))




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

;; (defmacro def-p2g (name body-func)
;;   `(defun ,name (mesh mp)
;;      (declare (cl-mpm/mesh::mesh mesh)
;;               (cl-mpm/particle:particle mp))
;;      (iterate-over-neighbours
;;       mesh mp
;;       (lambda (mesh mp node svp grads fsvp fgrads)
;;         (declare
;;          (cl-mpm/particle:particle mp)
;;          (cl-mpm/mesh::node node)
;;          (double-float svp))
;;         (funcall ,body-func (mesh mp node svp))))))

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
          )
         (with-accessors ((node-vel   cl-mpm/mesh:node-velocity)
                          (node-active  cl-mpm/mesh::node-active)
                          (node-mass  cl-mpm/mesh:node-mass)
                          (node-volume  cl-mpm/mesh::node-volume)
                          (node-volume-true  cl-mpm/mesh::node-volume-true)
                          (node-svp-sum  cl-mpm/mesh::node-svp-sum)
                          (node-force cl-mpm/mesh:node-force)
                          (node-p-wave cl-mpm/mesh::node-pwave)
                          (node-damage cl-mpm/mesh::node-damage)
                          (node-lock  cl-mpm/mesh:node-lock)) node
           (declare (type double-float node-mass node-volume mp-volume mp-pmod mp-damage node-svp-sum svp node-p-wave node-damage)
                    (type sb-thread:mutex node-lock))
           (sb-thread:with-mutex (node-lock)
             (setf node-active t)
             (incf node-mass
                   (* mp-mass svp))
             (incf node-volume
                   (* mp-volume svp))
             (incf node-p-wave
                   (* mp-pmod svp))
             (incf node-damage
                   (* mp-damage svp))
             (incf node-svp-sum svp)
             (cl-mpm/fastmaths::fast-fmacc node-vel mp-vel (* mp-mass svp))
             )
           ;;Ideally we include these generic functions for special mapping operations, however they are slow
           ;; (special-p2g mp node svp dsvp)
           )))
      ))
  (values))

(declaim (notinline p2g))
(defun p2g (mesh mps)
  "Map particle momentum to the grid"
  (declare (type (vector cl-mpm/particle:particle) mps) (cl-mpm/mesh::mesh mesh))
  (iterate-over-mps
   mps
   (lambda (mp)
     (p2g-mp mesh mp))))

(declaim (notinline p2g-force-mp)
         (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle) (values)) p2g-force-mp)
         )
(defun p2g-force-mp (mesh mp)
  "Map particle forces to the grid for one mp"
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (with-accessors ((mp-vel  cl-mpm/particle:mp-velocity)
                   (mp-mass cl-mpm/particle:mp-mass)
                   ) mp
    (declare (type double-float mp-mass))
    (let ((dsvp (cl-mpm/utils::dsvp-3d-zeros)))
      ;; (declare (dynamic-extent dsvp))
      (iterate-over-neighbours
       mesh mp
       (lambda (mesh mp node svp grads fsvp fgrads)
         (declare
          (cl-mpm/particle:particle mp)
          (cl-mpm/mesh::node node)
          (double-float svp))
         (with-accessors ((node-vel   cl-mpm/mesh:node-velocity)
                          (node-active  cl-mpm/mesh:node-active)
                          (node-mass  cl-mpm/mesh:node-mass)
                          (node-force cl-mpm/mesh:node-force)
                          (node-int-force cl-mpm/mesh::node-internal-force)
                          (node-ext-force cl-mpm/mesh::node-external-force)
                          (node-lock  cl-mpm/mesh:node-lock))
             node
           (declare (double-float node-mass)
                    (boolean node-active)
                    (sb-thread:mutex node-lock)
                    (magicl:matrix/double-float node-vel node-force node-int-force node-ext-force))
           (when node-active
             (cl-mpm/shape-function::assemble-dsvp-3d-prealloc grads dsvp)
             (sb-thread:with-mutex (node-lock)
               (det-ext-force mp node svp node-ext-force)
               (det-int-force-unrolled mp grads node-int-force)
               ;; (det-int-force mp dsvp node-int-force)
               )))))))
  (values))
(declaim (notinline p2g-force-mp-2d)
         (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle) (values)) p2g-force-mp-2d)
         )
(defun p2g-force-mp-2d (mesh mp)
  "Map particle forces to the grid for one mp"
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp))
  (with-accessors ((mp-vel  cl-mpm/particle:mp-velocity)
                   (mp-mass cl-mpm/particle:mp-mass)
                   ) mp
    (declare (type double-float mp-mass))
    (iterate-over-neighbours
     mesh mp
     (lambda (mesh mp node svp grads fsvp fgrads)
       (declare
        (cl-mpm/particle:particle mp)
        (cl-mpm/mesh::node node)
        (double-float svp))
       (with-accessors ((node-vel   cl-mpm/mesh:node-velocity)
                        (node-active  cl-mpm/mesh:node-active)
                        (node-mass  cl-mpm/mesh:node-mass)
                        (node-force cl-mpm/mesh:node-force)
                        (node-int-force cl-mpm/mesh::node-internal-force)
                        (node-ext-force cl-mpm/mesh::node-external-force)
                        (node-lock  cl-mpm/mesh:node-lock)) node
         (declare (double-float node-mass)
                  (boolean node-active)
                  (sb-thread:mutex node-lock)
                  (magicl:matrix/double-float node-vel node-force node-int-force node-ext-force))
         (when node-active
           (sb-thread:with-mutex (node-lock)
             (det-ext-force-2d mp node svp node-ext-force)
             (det-int-force-unrolled-2d mp grads node-int-force)
             ;; (det-ext-force mp node svp node-ext-force)
             ;; (det-int-force mp grads node-int-force)
             ))))))
  (values))

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
             `(defun ,name (mesh mp dt)
                (declare (cl-mpm/mesh::mesh mesh)
                         (cl-mpm/particle:particle mp)
                         (double-float dt))
                "Map one MP from the grid"
                (with-accessors ((vel mp-velocity)
                                 (pos mp-position)
                                 (disp cl-mpm/particle::mp-displacement)
                                 (acc cl-mpm/particle::mp-acceleration)
                                 (nc cl-mpm/particle::mp-cached-nodes))
                    mp
                  (let* ((mapped-vel (cl-mpm/utils:vector-zeros)))
                    ;; (declare (dynamic-extent mapped-vel))
                    (progn
                      ;;With special operations we need to reset some params for g2p
                      ;; (reset-mps-g2p mp)
                      ;(setf temp 0d0)
                      )
                    (cl-mpm/fastmaths::fast-zero acc)
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
                                        (node-scalar cl-mpm/mesh::node-boundary-scalar)
                                        (node-active cl-mpm/mesh:node-active)
                                        ) node
                         (declare (double-float node-scalar)
                                  (boolean node-active))
                         (when node-active
                           (cl-mpm/fastmaths::fast-fmacc mapped-vel node-vel svp)
                           (cl-mpm/fastmaths::fast-fmacc acc node-acc svp)
                           ;; (incf temp (* svp node-scalar))
                           ;;With special operations we want to include this operation
                           #+cl-mpm-special (special-g2p mesh mp node svp grads)
                           )
                         )
                       ;; (g2p-mp-node mp node svp grads)
                       ))
                    ;;Update particle
                    (progn
                      ;;Invalidate shapefunction/gradient cache
                      (update-particle mesh mp dt)
                      (setf (fill-pointer nc) 0)
                      (setf (cl-mpm/particle::mp-penalty-contact-step mp) (cl-mpm/particle::mp-penalty-contact mp))
                      (unless (cl-mpm/particle::mp-penalty-contact mp)
                        (cl-mpm/fastmaths:fast-zero (cl-mpm/particle:mp-penalty-frictional-force mp))
                        (setf (cl-mpm/particle::mp-penalty-normal-force mp) 0d0))
                      (setf (cl-mpm/particle::mp-penalty-contact mp) nil)
                      ,@update
                      ))
                  ))

             ))

  (def-g2p-mp g2p-mp-flip
      (progn
        (cl-mpm/fastmaths::fast-scale! mapped-vel dt)
        (cl-mpm/fastmaths::fast-.+-vector pos mapped-vel pos)
        (cl-mpm/fastmaths::fast-.+-vector disp mapped-vel disp)
        (cl-mpm/fastmaths:fast-fmacc vel acc dt)))
  (def-g2p-mp g2p-mp-pic
      (progn
        (cl-mpm/utils::vector-copy-into mapped-vel vel)
        (cl-mpm/fastmaths::fast-scale! mapped-vel dt)
        (cl-mpm/fastmaths::fast-.+-vector pos mapped-vel pos)
        (cl-mpm/fastmaths::fast-.+-vector disp mapped-vel disp)))
  (def-g2p-mp g2p-mp-blend
      (let ((pic-value 1d-3)
            (pic-vel (cl-mpm/utils:vector-copy mapped-vel)))
        (cl-mpm/fastmaths::fast-scale! mapped-vel dt)
        (cl-mpm/fastmaths::fast-.+-vector pos mapped-vel pos)
        (cl-mpm/fastmaths::fast-.+-vector disp mapped-vel disp)
        (cl-mpm/fastmaths:fast-.+
         (cl-mpm/fastmaths:fast-scale-vector
          ;; FLIP value
          (cl-mpm/fastmaths:fast-.+ vel (cl-mpm/fastmaths:fast-scale-vector acc dt))
          (- 1d0 pic-value))
         ;; PIC update
         (cl-mpm/fastmaths:fast-scale-vector pic-vel pic-value)
         vel))
      ;; (let* ((pic-value 1d-2)
      ;;        ;(pic-vel (cl-mpm/utils:vector-copy mapped-vel))
      ;;        (vel-inc
      ;;          (cl-mpm/fastmaths::fast-scale!
      ;;           (cl-mpm/fastmaths::fast-.--vector
      ;;            acc
      ;;            (cl-mpm/fastmaths:fast-scale!
      ;;             (cl-mpm/fastmaths::fast-.--vector vel mapped-vel) pic-value))
      ;;           dt)))
      ;;   (cl-mpm/fastmaths::fast-.+-vector vel vel-inc vel)
      ;;   (let ((dx (cl-mpm/fastmaths::fast-.-
      ;;              (cl-mpm/fastmaths:fast-scale-vector vel dt)
      ;;              (cl-mpm/fastmaths:fast-scale-vector vel-inc (* 0.5d0 dt)))))
      ;;     (cl-mpm/fastmaths::fast-.+-vector pos dx pos)
      ;;     (cl-mpm/fastmaths::fast-.+-vector disp dx disp))
      ;;   )
    )
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
          (cl-mpm/fastmaths::fast-.+-vector pos dx pos)
          (cl-mpm/fastmaths::fast-.+-vector disp dx disp))
        )
  )
  )

(defun g2p-mp (mesh mp dt)
  (declare (cl-mpm/mesh::mesh mesh)
           (cl-mpm/particle:particle mp)
           (double-float dt))
  "Map one MP from the grid"
  (with-accessors ((mass mp-mass)
                   (vel mp-velocity)
                   (pos mp-position)
                   (disp cl-mpm/particle::mp-displacement)
                   (acc cl-mpm/particle::mp-acceleration)
                   (temp cl-mpm/particle::mp-boundary)
                   (strain-rate mp-strain-rate)
                   (vorticity cl-mpm/particle:mp-vorticity)
                   (nc cl-mpm/particle::mp-cached-nodes)
                   (fixed-velocity cl-mpm/particle::mp-fixed-velocity)
                   )
    mp
    (let* ((mapped-vel (cl-mpm/utils:vector-zeros)))
      (progn
        ;;With special operations we need to reset some params for g2p
        ;; (reset-mps-g2p mp)
        (setf temp 0d0))
      (cl-mpm/fastmaths::fast-zero acc)
      ;; Map variables
      (iterate-over-neighbours
       mesh mp
       (lambda (mesh mp node svp grads fsvp fgrads)
         (declare
          (ignore mp mesh fsvp fgrads)
          (cl-mpm/mesh::node node)
          (cl-mpm/particle:particle mp)
          ( double-float svp))
         (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
                          (node-acc cl-mpm/mesh:node-acceleration)
                          (node-scalar cl-mpm/mesh::node-boundary-scalar)
                          (node-active cl-mpm/mesh:node-active)
                          ) node
           (declare (double-float node-scalar temp)
                    (boolean node-active))
           (when node-active
             (cl-mpm/fastmaths::fast-fmacc mapped-vel node-vel svp)
             (cl-mpm/fastmaths::fast-fmacc acc node-acc svp)
             (incf temp (* svp node-scalar))
             ;;With special operations we want to include this operation
             #+cl-mpm-special (special-g2p mesh mp node svp grads)
             )
           )
         ;; (g2p-mp-node mp node svp grads)
         ))
      ;;Update particle
      (progn
        ;;Invalidate shapefunction/gradient cache
        (setf (fill-pointer nc) 0)

        (let ((pic-value 0d-3)
              (pic-vel (cl-mpm/utils:vector-copy mapped-vel)))

          ;;PIC
          ;; #+ :cl-mpm-pic (cl-mpm/utils::vector-copy-into mapped-vel vel)
          (cl-mpm/fastmaths::fast-scale! mapped-vel dt)
          (cl-mpm/fastmaths::simd-add pos mapped-vel)
          (cl-mpm/fastmaths::simd-add disp mapped-vel)
          (setf (cl-mpm/particle::mp-penalty-contact-step mp) (cl-mpm/particle::mp-penalty-contact mp))
          (unless (cl-mpm/particle::mp-penalty-contact mp)
            (cl-mpm/fastmaths:fast-zero (cl-mpm/particle:mp-penalty-frictional-force mp))
            (setf (cl-mpm/particle::mp-penalty-normal-force mp) 0d0))
          (setf (cl-mpm/particle::mp-penalty-contact mp) nil)
          ;;FLIP
          ;; #- :cl-mpm-pic (cl-mpm/fastmaths:fast-fmacc vel acc dt)
          ;;Blended pic-flip
          (cl-mpm/fastmaths:fast-.+
           (cl-mpm/fastmaths:fast-scale-vector
            ;;New FLIP value
            (cl-mpm/fastmaths:fast-.+ vel (cl-mpm/fastmaths:fast-scale-vector acc dt))
            (- 1d0 pic-value))
           (cl-mpm/fastmaths:fast-scale-vector
            ;;New FLIP value
           (cl-mpm/fastmaths:fast-scale-vector pic-vel) pic-value)
           vel))))))
(defgeneric pre-particle-update-hook (particle dt)
  )
(defmethod pre-particle-update-hook (particle dt))

(declaim (notinline g2p))
;; (defun g2p-flip (mesh mps dt)
;;   (declare (cl-mpm/mesh::mesh mesh) (array mps))
;;   "Map grid values to all particles"
;;   (iterate-over-mps
;;    mps
;;    (lambda (mp)
;;      (g2p-mp mesh mp dt))))

(defun g2p (mesh mps dt &optional (update-type :FLIP))
  (ecase update-type
    (:FLIP
     (iterate-over-mps
      mps
      (lambda (mp)
        (g2p-mp-flip mesh mp dt))))
    (:PIC
     (iterate-over-mps
      mps
      (lambda (mp)
        (g2p-mp-pic mesh mp dt))))
    (:BLEND
     (iterate-over-mps
      mps
      (lambda (mp)
        (g2p-mp-blend mesh mp dt))))
    (:BLEND-2ND-ORDER
     (iterate-over-mps
      mps
      (lambda (mp)
        (g2p-mp-blend-2nd-order mesh mp dt))))
    ))

(defgeneric special-update-node (mesh dt node damping)
  (:documentation "Update node method")
  (:method (mesh dt node damping)))

(defmethod special-update-node (mesh dt (node cl-mpm/mesh::node-thermal) damping)
  (with-accessors ((mass  node-mass)
                   (temp node-temperature)
                   (dtemp node-dtemp))
      node
    (progn
      (setf temp (+ (/ temp mass) (* dtemp dt)))
      )))

(defun calculate-kinematics (node)
  "Calculate velocity from momentum on a single node"
  (when (cl-mpm/mesh:node-active node)
      (with-accessors ( (mass  node-mass)
                        (vel   node-velocity))
        node
        (declare (double-float mass))
        (progn
          (cl-mpm/fastmaths::fast-scale! vel (/ 1.0d0 mass))))))

(declaim (notinline calculate-forces)
         (ftype (function (cl-mpm/mesh::node double-float double-float double-float) (vaules)) calculate-forces))
(defun calculate-forces (node damping dt mass-scale)
  "Update forces and nodal velocities with viscous damping"
  (when (cl-mpm/mesh:node-active node)
    (with-accessors ((mass  node-mass)
                     (vel   node-velocity)
                     (force node-force)
                     (force-ext cl-mpm/mesh::node-external-force)
                     (force-int cl-mpm/mesh::node-internal-force)
                     (acc   node-acceleration))
        node
      (declare (double-float mass dt damping mass-scale))
      (progn
        (cl-mpm/fastmaths:fast-zero acc)
        ;;Set acc to f/m
        ;; (cl-mpm/fastmaths:fast-fmacc force-ext vel (* -2d0 damping mass))
        (cl-mpm/fastmaths::fast-.+-vector force-int force-ext force)
        ;; (cl-mpm/fastmaths:fast-fmacc acc vel (/ (* damping -1d0) 1d0))
        (cl-mpm/fastmaths:fast-fmacc acc force (/ 1d0 (* mass mass-scale)))
        (cl-mpm/fastmaths:fast-fmacc acc vel (/ (* damping -1d0) mass-scale))
        ;; (cl-mpm/fastmaths:fast-fmacc acc vel (/ (* damping -1d0) (sqrt mass-scale)))
        ;; (cl-mpm/fastmaths:fast-fmacc acc vel (* damping -1d0))
        (cl-mpm/fastmaths:fast-fmacc vel acc dt)
        )))
  (values))
(defun calculate-forces-psudo-viscous (node damping dt mass-scale)
  "Update forces and nodal velocities with viscous damping - except without scaling by mass
This allows for a non-physical but viscous damping scheme that is robust to GIMP domains "
  (when (cl-mpm/mesh:node-active node)
    (with-accessors ( (mass  node-mass)
                     (vel   node-velocity)
                     (force node-force)
                     (force-ext cl-mpm/mesh::node-external-force)
                     (force-int cl-mpm/mesh::node-internal-force)
                     (acc   node-acceleration))
        node
      (declare (double-float mass dt damping mass-scale))
      (progn
        (magicl:scale! acc 0d0)
        ;;Set acc to f/m
        (cl-mpm/fastmaths::fast-.+-vector force-int force-ext force)
        (cl-mpm/fastmaths:fast-fmacc acc force (/ 1d0 (* mass mass-scale)))
        (cl-mpm/fastmaths:fast-fmacc acc vel (* damping -1d0))
        (cl-mpm/fastmaths:fast-fmacc vel acc dt)
        )))
  (values))

(defun calculate-forces-cundall (node damping dt mass-scale)
  "Apply cundall damping to the system"
  (when (cl-mpm/mesh:node-active node)
    (with-accessors ((mass  node-mass)
                     (vel   node-velocity)
                     (force node-force)
                     (force-ext cl-mpm/mesh::node-external-force)
                     (force-int cl-mpm/mesh::node-internal-force)
                     (acc   node-acceleration)
                     )
        node
      (declare (double-float mass dt damping mass-scale))
      (progn
        (magicl:scale! acc 0d0)
        ;;Set acc to f/m
        (cl-mpm/fastmaths::fast-.+-vector force-int force-ext force)
        (let* ((vel-sign (cl-mpm/utils:vector-zeros))
               (f-s (magicl::matrix/double-float-storage force))
               (v-s (magicl::matrix/double-float-storage vel))
               (fnorm (cl-mpm/fastmaths::mag force))
               (vs-s (magicl::matrix/double-float-storage vel-sign)))
          (loop for i from 0 to 2
                do
                   ;;Possible form of damping
                   (when t;(> (* (aref f-s i) (aref v-s i)) 0)
                     (incf (aref f-s i)
                           (*
                            (signum (aref v-s i))
                            (abs (aref f-s i))
                            damping
                            -1d0))))
          (cl-mpm/fastmaths:fast-fmacc acc force (/ 1d0 (* mass mass-scale))))
        (cl-mpm/fastmaths:fast-fmacc vel acc dt))))
  (values))
(defun calculate-forces-cundall-conservative (node damping dt mass-scale)
  "Apply cundall damping to the system only when doing negative work"
  (when (cl-mpm/mesh:node-active node)
    (with-accessors ((mass  node-mass)
                     (vel   node-velocity)
                     (force node-force)
                     (acc   node-acceleration)
                     )
        node
      (declare (double-float mass dt damping))
      (progn
        (magicl:scale! acc 0d0)
        ;;Set acc to f/m

        (let* ((vel-sign (cl-mpm/utils:vector-zeros))
               (f-s (magicl::matrix/double-float-storage force))
               (v-s (magicl::matrix/double-float-storage vel))
               (fnorm (cl-mpm/fastmaths::mag force))
               (vs-s (magicl::matrix/double-float-storage vel-sign)))
          (loop for i from 0 to 2
                do
                   ;;Possible form of damping
                   (when (> (* (aref f-s i) (aref v-s i)) 0)
                     (incf (aref f-s i)
                           (*
                            (signum (aref v-s i))
                            (abs (aref f-s i))
                            damping
                            -1d0)))
                )
          (cl-mpm/fastmaths:fast-fmacc acc force (/ 1d0 (* mass mass-scale))))
        (cl-mpm/fastmaths:fast-fmacc vel acc dt))))
  (values))

(defun update-node-kinematics (mesh dt)
  (iterate-over-nodes mesh
                      (lambda (node)
                        (calculate-kinematics node))))
(defgeneric update-node-forces (sim)
  (:documentation "Update the acceleration from forces and apply any damping"))
(defmethod update-node-forces ((sim mpm-sim))
  (with-accessors ((damping sim-damping-factor)
                   (mass-scale sim-mass-scale)
                   (mesh sim-mesh)
                   (dt sim-dt))
      sim
    (iterate-over-nodes
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (calculate-forces node damping dt mass-scale))))))

(defmethod update-node-forces ((sim mpm-sim-quasi-static))
  (with-accessors ((damping sim-damping-factor)
                   (mass-scale sim-mass-scale)
                   (mesh sim-mesh)
                   (dt sim-dt))
      sim
    ;; (break)
    (iterate-over-nodes
     mesh
     (lambda (node)
       (when (cl-mpm/mesh:node-active node)
         (calculate-forces-cundall node damping dt mass-scale))))))



(defun apply-bcs (mesh bcs dt)
  "Apply all normal bcs onto the mesh"
  (declare (cl-mpm/mesh::mesh mesh))
  (with-accessors ((nodes  mesh-nodes)
                   (nD     mesh-nD)
                   (mc     mesh-count)) mesh
                                        ;each bc is a list (pos value)
    (lparallel:pdotimes (i (length bcs))
      (let ((bc (aref bcs i)))
        (when bc
          (with-accessors ((node cl-mpm/bc::bc-node)
                           (index cl-mpm/bc::bc-index))
              bc
            (if (or node
                    (not index))
                (cl-mpm/bc:apply-bc bc node mesh dt)
                (progn
                  (setf node (cl-mpm/mesh:get-node mesh index))
                  (if node
                      (cl-mpm/bc:apply-bc bc node mesh dt)
                      (error "BC attempted to get a nil node ~A ~A" bc index))))))))))


;Could include this in p2g but idk
(declaim (notinline calculate-strain-rate)
         (ftype (function (cl-mpm/mesh::mesh  cl-mpm/particle:particle double-float)) calculate-strain-rate))
(defun calculate-strain-rate (mesh mp dt)
  "Calculate the strain rate, stretch rate and vorticity"
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt)
           (optimize (speed 3) (safety 1))
           )
  (with-accessors ((strain-rate cl-mpm/particle:mp-strain-rate)
                   (vorticity cl-mpm/particle:mp-vorticity)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (stretch-tensor-fbar cl-mpm/particle::mp-stretch-tensor-fbar)
                   (velocity-rate cl-mpm/particle::mp-velocity-rate)
                   ) mp
    (declare (magicl:matrix/double-float strain-rate vorticity stretch-tensor stretch-tensor-fbar velocity-rate))
        (progn
          ;; (cl-mpm/fastmaths::fast-zero strain-rate)
          ;; (cl-mpm/fastmaths::fast-zero vorticity)
          (cl-mpm/fastmaths::fast-zero stretch-tensor)
          (cl-mpm/fastmaths::fast-zero stretch-tensor-fbar)
          (let (
                ;; (stretch-dsvp (stretch-dsvp-3d-zeros))
                ;; (temp-mult (cl-mpm/utils::stretch-dsvp-voigt-zeros))
                ;; (temp-add (cl-mpm/utils::matrix-zeros))
                )
            ;; (declare (magicl:matrix/double-float stretch-dsvp temp-mult temp-add)
            ;;          ;; (dynamic-extent stretch-dsvp temp-mult temp-add)
            ;;          )
            (iterate-over-neighbours
             mesh mp
             (lambda (mesh mp node svp grads fsvp fgrads)
               (declare (ignore mp svp fsvp))
               (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
                                (node-active cl-mpm/mesh:node-active))
                   node
                 (declare (magicl:matrix/double-float node-vel)
                          (boolean node-active))
                 (when node-active

                   (cl-mpm/shape-function::@-combi-assemble-dstretch-3d grads node-vel stretch-tensor)
                   (cl-mpm/shape-function::@-combi-assemble-dstretch-3d fgrads node-vel stretch-tensor-fbar)
                   )))))

            ;; (cl-mpm/utils::stretch-to-sym stretch-tensor strain-rate)
            ;; (cl-mpm/utils::stretch-to-skew stretch-tensor vorticity)
            ;; (aops:copy-into (cl-mpm/utils::fast-storage velocity-rate) (cl-mpm/utils::fast-storage strain-rate))
            ;; (setf velocity-rate (magicl:scale strain-rate 1d0))
            (cl-mpm/fastmaths::fast-scale! stretch-tensor dt)
            (cl-mpm/fastmaths::fast-scale! stretch-tensor-fbar dt)
            ;; (cl-mpm/fastmaths::fast-scale! strain-rate dt)
            ;; (cl-mpm/fastmaths::fast-scale! vorticity dt)
            )))


(defun rotation-matrix (degrees)
  (let ((angle (/ (* pi degrees) 180)))
    (magicl:from-list (list (cos angle) (- (sin angle))
                            (sin angle) (cos angle))
                      '(2 2)
                      :type 'double-float
                      )))
(defun update-strain-linear (mesh mp dt fbar)
  "Linear small strain update algorithm"
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (volume-0 cl-mpm/particle::mp-volume-0)
                   (strain cl-mpm/particle:mp-strain)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   ) mp
    (let ((df (calculate-df mesh mp fbar))
          (dstrain (magicl:scale strain-rate dt)))
                   (progn
                     (setf def (cl-mpm/fastmaths::fast-@-matrix-matrix df def))
                     (setf strain (cl-mpm/fastmaths::fast-.+-voigt strain dstrain))
                     (setf volume (* volume-0 (cl-mpm/fastmaths:det-3x3 def)))
                     ;(update-domain-corner mesh mp dt)
                     (update-domain-midpoint mesh mp dt)))))

(declaim (inline update-strain-kirchoff))
(declaim (ftype (function (cl-mpm/mesh::mesh
                           cl-mpm/particle:particle
                           double-float
                           boolean) (values))
               update-strain-kirchoff))

(defun update-strain-kirchoff-noupdate (mesh mp dt fbar)
  "Finite strain kirchhoff strain update algorithm"
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (volume-0 cl-mpm/particle::mp-volume-0)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (df-inc cl-mpm/particle::mp-deformation-gradient-increment)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   (strain-rate-tensor cl-mpm/particle::mp-strain-rate-tensor)
                   (velocity-rate cl-mpm/particle::mp-velocity-rate)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0)
                   (eng-strain-rate cl-mpm/particle::mp-eng-strain-rate)
                   ) mp
    (declare (type double-float volume)
             (type magicl:matrix/double-float
                   domain))
    (progn
      (multiple-value-bind (df dj) (calculate-df mesh mp fbar)
        (progn
          (setf df-inc df)
          (setf def (cl-mpm/fastmaths::fast-@-matrix-matrix df def))
          ;; (cl-mpm/utils:voigt-copy-into strain strain-rate)
          (cl-mpm/ext:kirchoff-update strain df)
          ;; (cl-mpm/fastmaths:fast-.- strain strain-rate strain-rate)
          ;;Post multiply to turn to eng strain
                                        ;(setf volume (* volume (the double-float (cl-mpm/fastmaths:det-3x3 df))))
          ;; (setf volume (* volume-0 (the double-float (cl-mpm/fastmaths:det-3x3 def))))
          ;; (setf volume (* volume (the double-float dj)))
          ;; (when (<= volume 0d0)
          ;;   (error "Negative volume"))
          ))))
  (values))

(defun update-strain-kirchoff (mesh mp dt fbar)
  "Finite strain kirchhoff strain update algorithm"
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (volume-0 cl-mpm/particle::mp-volume-0)
                   (strain cl-mpm/particle:mp-strain)
                   (strain-n cl-mpm/particle:mp-strain-n)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (df-inc    cl-mpm/particle::mp-deformation-gradient-increment)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   ;; (strain-rate cl-mpm/particle:mp-strain-rate)
                   ;; (strain-rate-tensor cl-mpm/particle::mp-strain-rate-tensor)
                   (velocity-rate cl-mpm/particle::mp-velocity-rate)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0)
                   (eng-strain-rate cl-mpm/particle::mp-eng-strain-rate)
                   ) mp
    (declare (type double-float volume)
             (type magicl:matrix/double-float
                   domain))
    (progn
      (multiple-value-bind (df dj) (calculate-df mesh mp fbar)
        (progn
          (setf df-inc df)
          (setf def (cl-mpm/fastmaths::fast-@-matrix-matrix df def))
          (cl-mpm/utils:voigt-copy-into strain strain-n)
          (cl-mpm/ext:kirchoff-update strain df)
          ;; (cl-mpm/fastmaths:fast-.- strain strain-rate strain-rate)
          ;;Post multiply to turn to eng strain
                                        ;(setf volume (* volume (the double-float (cl-mpm/fastmaths:det-3x3 df))))
          ;(setf volume (* volume (the double-float (cl-mpm/fastmaths:det-3x3 df))))
          (setf volume (* volume (the double-float dj)))
          (when (<= volume 0d0)
            (error "Negative volume"))
          ;; ;;Stretch rate update
          ;; (update-domain-corner mesh mp dt)
          ;; (scale-domain-size mesh mp)
          ))))
  (values))
(declaim
 (inline update-strain-kirchoff-damaged)
 (ftype (function (cl-mpm/mesh::mesh
                           cl-mpm/particle:particle
                           double-float
                           boolean) (values))
                update-strain-kirchoff-damaged))
(defun update-strain-kirchoff-damaged (mesh mp dt fbar)
  "Finite strain kirchhoff strain update algorithm"
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   (domain cl-mpm/particle::mp-domain-size)
                   (eng-strain-rate cl-mpm/particle::mp-eng-strain-rate)
                   ) mp
    (declare (type double-float volume)
             (type magicl:matrix/double-float
                   domain))
    (progn
      (let ((df (calculate-df mesh mp fbar)))
        (progn
          (let ((def-n (cl-mpm/utils::matrix-zeros)))
            (magicl:mult df def :target def-n)
            (cl-mpm/utils::matrix-copy-into def-n def))
          ;; (aops:copy-into (magicl::matrix/double-float-storage eng-strain-rate)
          ;;                 (magicl::matrix/double-float-storage strain))
          ;; (when (> (abs (magicl:tref strain 4 0)) 0d0)
          ;;   (error "Pre Nonzero out of plane strain with ~A" (loop for v across (magicl::storage strain) collect v)))
          (cl-mpm/utils:voigt-copy-into strain strain-rate)
          (cl-mpm/ext:kirchoff-update strain df)
          (cl-mpm/fastmaths:fast-.- strain strain-rate strain-rate)
          ;; (when (> (abs (magicl:tref strain 4 0)) 0d0)
          ;;   (error "Post Nonzero out of plane strain with ~A" (loop for v across (magicl::storage strain) collect v)))
          ;; (magicl:.- eng-strain-rate strain eng-strain-rate)
          ;; (magicl:scale! eng-strain-rate (/ 1d0 dt))
          ;;Post multiply to turn to eng strain
          (setf volume (* volume (the double-float (cl-mpm/fastmaths:det-3x3 df))))
          (when (<= volume 0d0)
            (error "Negative volume"))
          (update-domain-stretch-rate-damage stretch-tensor (cl-mpm/particle::mp-damage mp) domain
                                             (cl-mpm/particle::mp-damage-domain-update-rate mp))
          (scale-domain-size mesh mp)
          )))
    )
  (values))
(defun update-strain-kirchoff-ugimp (mesh mp dt fbar)
  "Finite strain kirchhoff strain update algorithm"
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   (domain cl-mpm/particle::mp-domain-size)
                   (eng-strain-rate cl-mpm/particle::mp-eng-strain-rate)
                   ) mp
    (declare (type double-float volume)
             (type magicl:matrix/double-float
                   domain))
    (progn
      (let ((df (calculate-df mesh mp fbar)))
        (let ((def-n (cl-mpm/utils::matrix-zeros)))
          (magicl:mult df def :target def-n)
          (cl-mpm/utils::matrix-copy-into def-n def))
        (cl-mpm/utils:voigt-copy-into strain strain-rate)
        (cl-mpm/ext:kirchoff-update strain df)
        (cl-mpm/fastmaths:fast-.- strain strain-rate strain-rate)
        (setf volume (* volume (the double-float (cl-mpm/fastmaths:det-3x3 df))))
        (when (<= volume 0d0)
          (error "Negative volume")))))
  (values))

(defun update-strain-kirchoff-det (mesh mp dt fbar)
  "Finite strain kirchhoff strain update algorithm"
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   (domain cl-mpm/particle::mp-domain-size)
                   (eng-strain-rate cl-mpm/particle::mp-eng-strain-rate)
                   ) mp
    (declare (type double-float volume)
             (type magicl:matrix/double-float
                   domain))
    (progn
      (let ((df (calculate-df mesh mp fbar)))
        (let ((def-n (cl-mpm/utils::matrix-zeros)))
          (magicl:mult df def :target def-n)
          (cl-mpm/utils::matrix-copy-into def-n def))
        (cl-mpm/utils:voigt-copy-into strain strain-rate)
        (cl-mpm/ext:kirchoff-update strain df)
        (cl-mpm/fastmaths:fast-.- strain strain-rate strain-rate)
        (setf volume (* volume (the double-float (cl-mpm/fastmaths:det-3x3 df))))
        (update-domain-def mesh mp)
        (when (<= volume 0d0)
          (error "Negative volume")))))
  (values))



(defun update-stress-kirchoff (mesh mp dt fbar)
  "Update stress for a single mp"
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt)
           (optimize (speed 3) (safety 0) (debug 0)))
  (with-accessors ((stress cl-mpm/particle:mp-stress)
                   (stress-kirchoff cl-mpm/particle::mp-stress-kirchoff)
                   (volume cl-mpm/particle:mp-volume)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   ) mp
    (declare (magicl:matrix/double-float stress stress-kirchoff strain def strain-rate)
             (double-float volume))
    (progn
      (progn
        (calculate-strain-rate mesh mp dt)
        ;; Turn cauchy stress to kirchoff
        (cl-mpm/utils::voigt-copy-into stress-kirchoff stress)
        ;; Update our strains
        (update-strain-kirchoff mesh mp dt fbar)
        ;; Update our kirchoff stress with constitutive model
        (cl-mpm/utils::voigt-copy-into (cl-mpm/particle:constitutive-model mp strain dt) stress-kirchoff)
        ;; Check volume constraint!
        ;; (when (<= volume 0d0)
        ;;   (error "Negative volume"))
        ;; Turn kirchoff stress to cauchy
        (cl-mpm/utils::voigt-copy-into stress-kirchoff stress)
        (cl-mpm/fastmaths::fast-scale! stress (/ 1.0d0 (the double-float (cl-mpm/fastmaths:det-3x3 def))))
        ))))

(defun update-stress-kirchoff-mapped-jacobian (mesh mp dt fbar)
  "Update stress for a single mp"
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt)
           (optimize (speed 3) (safety 0) (debug 0)))
  (with-accessors ((stress cl-mpm/particle:mp-stress)
                   (stress-kirchoff cl-mpm/particle::mp-stress-kirchoff)
                   (volume cl-mpm/particle:mp-volume)
                   (volume-0 cl-mpm/particle::mp-volume-0)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   ) mp
    (declare (magicl:matrix/double-float stress stress-kirchoff strain def strain-rate)
             (double-float volume))
    (progn
      (progn
        ;; (calculate-strain-rate mesh mp dt)
        ;; Turn cauchy stress to kirchoff
        (cl-mpm/utils::voigt-copy-into stress-kirchoff stress)
        ;; Update our strains
        (update-strain-kirchoff-noupdate mesh mp dt fbar)
        ;; (scale-domain-size mesh mp)
        ;; Update our kirchoff stress with constitutive model
        (cl-mpm/utils::voigt-copy-into (cl-mpm/particle:constitutive-model mp strain dt) stress-kirchoff)
        ;; (cl-mpm/constitutive::linear-elastic-mat strain (cl-mpm/particle::mp-elastic-matrix mp) stress-kirchoff)
        ;; Check volume constraint!
        (when (<= volume 0d0)
          (error "Negative volume"))
        ;; Turn kirchoff stress to cauchy
        (cl-mpm/utils::voigt-copy-into stress-kirchoff stress)
        (cl-mpm/fastmaths::fast-scale! stress (/ 1.0d0 (the double-float (cl-mpm/fastmaths:det-3x3 def))))
        ))))

(defun update-stress-kirchoff-noscale (mesh mp dt fbar)
  "Update stress for a single mp"
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt)
           (optimize (speed 3) (safety 0) (debug 0)))
  (with-accessors ((stress cl-mpm/particle:mp-stress)
                   (stress-kirchoff cl-mpm/particle::mp-stress-kirchoff)
                   (volume cl-mpm/particle:mp-volume)
                   (volume-0 cl-mpm/particle::mp-volume-0)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   ) mp
    (declare (magicl:matrix/double-float stress stress-kirchoff strain def strain-rate)
             (double-float volume))
    (progn
      (progn
        (calculate-strain-rate mesh mp dt)
        ;; Turn cauchy stress to kirchoff
        (cl-mpm/utils::voigt-copy-into stress-kirchoff stress)
        ;; Update our strains
        (update-strain-kirchoff-noupdate mesh mp dt fbar)
        ;; (scale-domain-size mesh mp)
        ;; Update our kirchoff stress with constitutive model
        (cl-mpm/utils::voigt-copy-into (cl-mpm/particle:constitutive-model mp strain dt) stress-kirchoff)
        ;; (cl-mpm/constitutive::linear-elastic-mat strain (cl-mpm/particle::mp-elastic-matrix mp) stress-kirchoff)
        ;; Check volume constraint!
        (when (<= volume 0d0)
          (error "Negative volume"))
        ;; Turn kirchoff stress to cauchy
        (cl-mpm/utils::voigt-copy-into stress-kirchoff stress)
        (cl-mpm/fastmaths::fast-scale! stress (/ 1.0d0 (the double-float (cl-mpm/fastmaths:det-3x3 def))))
        ))))

(declaim (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle double-float boolean) (values)) update-stress-kirchoff-damaged))
(defun update-stress-kirchoff-damaged (mesh mp dt fbar)
  "Update stress for a single mp"
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt))
  (with-accessors ((stress cl-mpm/particle:mp-stress)
                   (stress-kirchoff cl-mpm/particle::mp-stress-kirchoff)
                   (volume cl-mpm/particle:mp-volume)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   ) mp
    (declare (magicl:matrix/double-float stress stress-kirchoff strain def strain-rate))
    (progn
      ;; (unless fbar)
      (calculate-strain-rate mesh mp dt)

      ;; Turn cauchy stress to kirchoff
      (cl-mpm/utils::voigt-copy-into stress-kirchoff stress)
      ;; Update our strains
      (update-strain-kirchoff-damaged mesh mp dt fbar)
      ;; Update our kirchoff stress with constitutive model
      (cl-mpm/utils::voigt-copy-into (cl-mpm/particle:constitutive-model mp strain dt) stress-kirchoff)
      ;; Check volume constraint!
      (when (<= volume 0d0)
        (error "Negative volume"))
      ;; Turn kirchoff stress to cauchy
      (cl-mpm/fastmaths::fast-scale! stress (/ 1.0d0 (the double-float (cl-mpm/fastmaths:det-3x3 def))))
      )))
(declaim (ftype (function (cl-mpm/mesh::mesh cl-mpm/particle:particle double-float boolean) (values)) update-stress-kirchoff-ugimp))
(defun update-stress-kirchoff-ugimp (mesh mp dt fbar)
  "Update stress for a single mp"
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt))
  (with-accessors ((stress cl-mpm/particle:mp-stress)
                   (stress-kirchoff cl-mpm/particle::mp-stress-kirchoff)
                   (volume cl-mpm/particle:mp-volume)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   ) mp
    (declare (magicl:matrix/double-float stress stress-kirchoff strain def strain-rate))
    (progn
      ;; (unless fbar)
      (calculate-strain-rate mesh mp dt)

      ;; Turn cauchy stress to kirchoff
      (cl-mpm/utils::voigt-copy-into stress-kirchoff stress)
      ;; Update our strains
      (update-strain-kirchoff-ugimp mesh mp dt fbar)
      ;; Update our kirchoff stress with constitutive model
      (cl-mpm/utils::voigt-copy-into (cl-mpm/particle:constitutive-model mp strain dt) stress-kirchoff)
      ;; Check volume constraint!
      (when (<= volume 0d0)
        (error "Negative volume"))
      ;; Turn kirchoff stress to cauchy
      (cl-mpm/fastmaths::fast-scale! stress (/ 1.0d0 (the double-float (cl-mpm/fastmaths:det-3x3 def))))
      )))

(defun update-stress-kirchoff-det (mesh mp dt fbar)
  "Update stress for a single mp"
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt))
  (with-accessors ((stress cl-mpm/particle:mp-stress)
                   (stress-kirchoff cl-mpm/particle::mp-stress-kirchoff)
                   (volume cl-mpm/particle:mp-volume)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   ) mp
    (declare (magicl:matrix/double-float stress stress-kirchoff strain def strain-rate))
    (progn
      (calculate-strain-rate mesh mp dt)
      ;; Turn cauchy stress to kirchoff
      (cl-mpm/utils::voigt-copy-into stress-kirchoff stress)
      ;; Update our strains
      ;(update-strain-kirchoff-det mesh mp dt fbar)
      (update-strain-kirchoff-det mesh mp dt fbar)
      ;; Update our kirchoff stress with constitutive model
      (cl-mpm/utils::voigt-copy-into (cl-mpm/particle:constitutive-model mp strain dt) stress-kirchoff)
      ;; Check volume constraint!
      (when (<= volume 0d0)
        (error "Negative volume"))
      ;; Turn kirchoff stress to cauchy
      (cl-mpm/fastmaths::fast-scale! stress (/ 1.0d0 (the double-float (cl-mpm/fastmaths:det-3x3 def))))
      )))

(defun update-stress-linear (mesh mp dt fbar)
  "Update stress for a single mp"
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt))
  (with-accessors ((stress cl-mpm/particle:mp-stress)
                   (stress-kirchoff cl-mpm/particle::mp-stress-kirchoff)
                   (volume cl-mpm/particle:mp-volume)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   ) mp
    (declare (magicl:matrix/double-float stress stress-kirchoff strain def strain-rate))
    (progn
      ;;   ;;For no FBAR we need to update our strains
      (progn
        (calculate-strain-rate mesh mp dt)

        ;; Update our strains
        (update-strain-linear mesh mp dt fbar)

        ;; Update our kirchoff stress with constitutive model
        (setf stress (cl-mpm/particle:constitutive-model mp strain dt))
        ;; Check volume constraint!
        (when (<= volume 0d0)
          (error "Negative volume"))
        ))))

(defun update-particle-kirchoff (mesh mp dt)
  ;; (with-accessors ((volume cl-mpm/particle:mp-volume)
  ;;                  (volume-0 cl-mpm/particle::mp-volume-0)
  ;;                  (def    cl-mpm/particle:mp-deformation-gradient)
  ;;                  (df    cl-mpm/particle::mp-deformation-gradient-increment)
  ;;                  )
  ;;     mp
  ;;   ;; (setf volume (* volume-0 (the double-float (cl-mpm/fastmaths:det-3x3 def))))
  ;;   (setf volume (* volume (the double-float (cl-mpm/fastmaths:det-3x3 df))))
  ;;   (when (<= volume 0d0)
  ;;     (error "Negative volume")))
  )

(defgeneric update-particle (mesh mp dt)
  (:documentation "End of step update"))
(defmethod update-particle (mesh (mp cl-mpm/particle:particle) dt)
  (update-particle-kirchoff mesh mp dt)
  (update-domain-corner mesh mp dt)
  ;; (update-domain-stretch mesh mp dt)
  (cl-mpm::scale-domain-size mesh mp)
  )

(defgeneric update-stress-mp (mesh mp dt fbar)
  (:documentation "A mp dependent stress update scheme"))

(defmethod update-stress-mp (mesh (mp cl-mpm/particle::particle) dt fbar)
  (update-stress-kirchoff mesh mp dt fbar)
  ;; (update-stress-linear mesh mp dt fbar)
  )

(defun map-jacobian (mesh mp dt)
  (with-accessors ((stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (volume cl-mpm/particle:mp-volume)
                   (def cl-mpm/particle:mp-deformation-gradient))
      mp
    (let* ((df (cl-mpm/utils::matrix-from-list '(1d0 0d0 0d0
                                                 0d0 1d0 0d0
                                                 0d0 0d0 1d0))))
      (cl-mpm/fastmaths::fast-.+ df stretch-tensor df)
      (let ((j-inc (cl-mpm/fastmaths:det-3x3 df))
            (j-n
              1d0
              ;; (cl-mpm/fastmaths:det-3x3 def)
              ))
        (iterate-over-neighbours
         mesh mp
         (lambda (mesh mp node svp grads fsvp fgrads)
           (with-accessors ((node-active cl-mpm/mesh:node-active)
                            (node-j-inc cl-mpm/mesh::node-jacobian-inc)
                            (node-volume cl-mpm/mesh::node-volume)
                            (node-svp cl-mpm/mesh::node-svp-sum)
                            (node-lock cl-mpm/mesh::node-lock))
               node
             (declare (double-float node-j-inc svp volume j-inc j-n node-volume))
             (when node-active
               (sb-thread:with-mutex (node-lock)
                 ;(incf node-j-inc (* svp j-inc j-n (/ volume node-volume)))
                 (incf node-j-inc (* j-inc j-n (/ svp node-svp)))
                 )))))))))

(declaim (inline calculate-df)
         (ftype (function (cl-mpm/mesh::mesh
                           cl-mpm/particle:particle
                           boolean) magicl:matrix/double-float)
                calculate-df))
(defun calculate-df (mesh mp fbar)
  (with-accessors ((dstrain cl-mpm/particle::mp-strain-rate)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (stretch-tensor-fbar cl-mpm/particle::mp-stretch-tensor-fbar)
                   (jfbar cl-mpm/particle::mp-j-fbar)
                   (def cl-mpm/particle:mp-deformation-gradient)
                   (pos cl-mpm/particle:mp-position)
                   )
      mp
    (let* ((df (cl-mpm/utils::matrix-from-list '(1d0 0d0 0d0
                                                 0d0 1d0 0d0
                                                 0d0 0d0 1d0)))
           (dJ 1d0)
           )
      (cl-mpm/fastmaths::fast-.+-matrix df stretch-tensor df)
      (setf dJ (cl-mpm/fastmaths:det-3x3 df))
      ;;Explicit fbar
      (when fbar
        (if nil;;t exp: nil Coobs
            (progn
              (let ((j-inc (cl-mpm/fastmaths:det-3x3 df))
                    (j-n
                      1d0
                      ;; (cl-mpm/fastmaths:det-3x3 def)
                         )
                    (gather-j 0d0)
                    (nd (cl-mpm/mesh::mesh-nd mesh))
                    (svp-sum 0d0)
                    )
                (iterate-over-neighbours
                 mesh mp
                 (lambda (mesh mp node svp grads fsvp fgrads)
                   (with-accessors ((node-active cl-mpm/mesh:node-active)
                                    (node-j-inc cl-mpm/mesh::node-jacobian-inc)
                                    (node-volume cl-mpm/mesh::node-volume))
                       node
                     (when node-active
                       (incf svp-sum svp)
                       ;; (incf gather-j (/ (* svp node-j-inc) node-volume))
                       (incf gather-j (* svp node-j-inc))
                       ))))
                ;; (setf gather-j (/ gather-j svp-sum))

                (setf (cl-mpm/particle::mp-debug-j mp) (/ gather-j (* j-inc j-n))
                      (cl-mpm/particle::mp-debug-j-gather mp) gather-j)

                (cl-mpm/fastmaths:fast-scale!
                 df
                 (expt
                  (the double-float (/ gather-j (* j-inc j-n)))
                  (/ 1 nd)))
                (when (= nd 2)
                  (setf (magicl:tref df 2 2) 1d0))
                )
              )
            ;;Coombs fbar
            (progn
              (let* ((df-fbar (cl-mpm/utils::matrix-from-list '(1d0 0d0 0d0
                                                                0d0 1d0 0d0
                                                                0d0 0d0 1d0)))
                     (nd (cl-mpm/mesh::mesh-nd mesh)))
                (cl-mpm/fastmaths::fast-.+-matrix df-fbar stretch-tensor-fbar df-fbar)
                (setf (cl-mpm/particle::mp-debug-j mp) (cl-mpm/fastmaths:det-3x3 df)
                      (cl-mpm/particle::mp-debug-j-gather mp) (cl-mpm/fastmaths:det-3x3 df-fbar))
                (cl-mpm/fastmaths::fast-scale!
                 df
                 (expt
                  (the double-float (/ (cl-mpm/fastmaths:det-3x3 df-fbar)
                                       (cl-mpm/fastmaths:det-3x3 df)))
                  (the double-float (/ 1d0 nd))))
                (when (= nd 2)
                  (setf (magicl:tref df 2 2) 1d0))
                ))))
      (values df dJ))))
(defgeneric post-stress-step (mesh mp dt))
(defmethod post-stress-step (mesh mp dt))
(defmethod post-stress-step (mesh (mp cl-mpm/particle::particle) dt)
  (declare (ignore mesh mp dt)))

(declaim (ftype (function (cl-mpm/mesh::mesh
                           (array cl-mpm/particle:particle)
                           double-float
                           &optional boolean) (values)) update-stress))
(defun update-stress (mesh mps dt &optional (fbar nil))
  "Update all stresses, with optional f-bar"
  (declare ((array cl-mpm/particle:particle) mps) (cl-mpm/mesh::mesh mesh))
  ;; (iterate-over-mps
  ;;  mps
  ;;  (lambda (mp)
  ;;    (calculate-strain-rate mesh mp dt)
  ;;    (map-jacobian mesh mp dt)))
  (iterate-over-mps
   mps
   (lambda (mp)
     (update-stress-mp mesh mp dt fbar)
     (post-stress-step mesh mp dt)
     ))
  (values))

(declaim (notinline reset-grid))
(defun reset-grid (mesh)
  "Reset all nodes on the grid"
  (declare (cl-mpm/mesh::mesh mesh))
  (iterate-over-nodes
   mesh
   (lambda (node)
     (when (cl-mpm/mesh:node-active node)
       (cl-mpm/mesh:reset-node node)))))
(declaim (notinline reset-grid-velocity))
(defun reset-grid-velocity (mesh)
  "Reset all velocity map on grid for MUSL"
  (declare (cl-mpm/mesh::mesh mesh))
  (iterate-over-nodes
   mesh
   (lambda (node)
     (when (cl-mpm/mesh:node-active node)
       (cl-mpm/mesh::reset-node-velocity node)))))


(defun filter-grid (mesh mass-thresh)
  "Filter out nodes with very small massess"
  (declare (cl-mpm/mesh::mesh mesh) (double-float mass-thresh))
  (iterate-over-nodes
   mesh
   (lambda (node)
     (when (and (cl-mpm/mesh:node-active node)
                (<= (cl-mpm/mesh:node-mass node) mass-thresh))
       (cl-mpm/mesh:reset-node node)))))

(defun filter-grid-velocity (mesh mass-thresh)
  "Filter out nodes with very small massess"
  (declare (cl-mpm/mesh::mesh mesh) (double-float mass-thresh))
  (iterate-over-nodes
   mesh
   (lambda (node)
     (when (and (cl-mpm/mesh:node-active node)
                (<= (cl-mpm/mesh:node-mass node) mass-thresh))
       (cl-mpm/mesh::reset-node-velocity node)))))

(defun filter-grid-volume (mesh volume-ratio)
  (declare (cl-mpm/mesh::mesh mesh) (double-float volume-ratio))
  (let ((nodes (cl-mpm/mesh:mesh-nodes mesh)))
    (dotimes (i (array-total-size nodes))
      (let ((node (row-major-aref nodes i)))
        (when (and (cl-mpm/mesh:node-active node)
                   (< (the double-float (cl-mpm/mesh:node-mass node))
                      (* volume-ratio (the double-float (cl-mpm/mesh::node-volume-true node)))))
          (setf (cl-mpm/mesh::node-active node) nil)
          (cl-mpm/mesh:reset-node node)
          )))))

;; (defun delete-mps-func (mps func)
;;   (declare (function func)
;;            (vector mps))
;;   (let* ((mp-count (length mps))
;;          (bit-vector (make-array mp-count :element-type 'bit :initial-element 0)))
;;     (declare (dynamic-extent bit-vector))
;;     (lparallel:pdotimes (i mp-count)
;;       (setf (aref bit-vector i) (if (funcall func (aref mps i)) 1 0)))
;;     (loop for i fixnum from 0 below mp-count
;;           do (when (= 1 (aref bit-vector i))
;;               (setf (aref mps i) (aref mps (- mp-count 1)))
;;               (setf (aref bit-vector i) (aref bit-vector (- mp-count 1)))
;;               (decf mp-count)
;;               (decf i)
;;               (decf (fill-pointer mps))
;;               )))
;;   (when (and (not (adjustable-array-p mps))
;;              (= (length mps) 0))
;;     (setf mps (make-array 0 :adjustable t :fill-pointer 0))))

(defgeneric sim-add-mp (sim mp)
  (:documentation "A function to append an mp to a simulation"))

(defmethod sim-add-mp (sim mp)
  (with-accessors ((uid-counter sim-unique-index-counter)
                   (lock sim-unique-index-lock))
    sim
      (sb-thread:with-mutex (lock)
        (setf (cl-mpm/particle::mp-unique-index mp) uid-counter)
        (incf uid-counter))
    (vector-push-extend mp (cl-mpm:sim-mps sim))))

(defgeneric remove-mps-func (sim func)
  (:documentation "A function for removing mps from a sim"))

(defmethod remove-mps-func (sim func)
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (declare ((vector cl-mpm/particle::particle) mps))
    (when (> (length mps) 0)
      ;; (delete-mps-func mps func)
      (setf mps
            (delete-if func mps)))
    ;;Sometimes when compacting the array; sbcl will just discard make and unadjustable array in place which is a bit wild
    (when (and (not (adjustable-array-p mps))
               (= (length mps) 0))
      (setf mps (make-array 0 :adjustable t :fill-pointer 0)))
    )
  (values))

(defun remove-material-damaged (sim)
  "Remove material points that have strain energy density exceeding fracture toughness"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)
                   (instant-damage-removal cl-mpm::sim-mp-damage-removal-instant)
                   )
      sim
    (let ((h (cl-mpm/mesh:mesh-resolution mesh)))
      (remove-mps-func
       sim
       (lambda (mp)
         (with-accessors ((damage cl-mpm/particle:mp-damage)
                          (def cl-mpm/particle::mp-deformation-gradient))
             mp
           (and (>= damage 0.9d0)
                (or instant-damage-removal
                    (damage-removal-criteria mp h) )
                )))))
    ))

(defun damage-removal-criteria (mp h)
  "Some criteria for when we should remove fully damaged MPs"
  (with-accessors ((def cl-mpm/particle:mp-deformation-gradient)
                   (lens cl-mpm/particle::mp-domain-size)
                   (strain cl-mpm/particle::mp-strain)
                   (lens-0 cl-mpm/particle::mp-domain-size-0)
                   ;; (ft cl-mpm/particle::mp-initiation-stress)
                   ;; (ductility cl-mpm/particle::mp-ductility)
                   ;; (E cl-mpm/particle::mp-e)
                   )
      mp
    ;; (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix strain)))
    (let* ((l-factor 2.00d0)
           ;; (e0 (/ ft E))
           ;; (ef (/ (* ft (+ ductility 1d0)) (* 2d0 E)))
           ;; (ef (* ef 50d0))
           ;; (e_max (reduce #' max l))
           )
      (cond
        ;; ((< ef e_max) :x)
        ;; ((< ef (tref strain 0 0)) :x)
        ;; ((< ef (tref strain 1 0)) :y)
        ;; ((< ef (tref strain 2 0)) :z)
        ((< (* l-factor (varef lens-0 0)) (varef lens 0)) :x)
        ((< (* l-factor (varef lens-0 1)) (varef lens 1)) :y)
        ;; ((and (< (* l-factor (tref lens-0 2 0)) (tref lens 2 0))
        ;;       (> (tref lens-0 2 0) 0d0)
        ;;       ) :z)
        (t nil)
        ))))

(defparameter *max-split-depth* 3)
(defun split-criteria (mp h)
  "Some numerical splitting estimates"
  (with-accessors ((def cl-mpm/particle:mp-deformation-gradient)
                   (lens cl-mpm/particle::mp-domain-size)
                   (lens-0 cl-mpm/particle::mp-domain-size-0)
                   (split-depth cl-mpm/particle::mp-split-depth)
                   )
      mp
    (when (< split-depth *max-split-depth*)
      (let ((l-factor 1.00d0)
            (h-factor (* 0.6d0 h))
            (s-factor 1.5d0))
        (cond
          ((< h-factor (varef lens 0)) :x)
          ((< h-factor (varef lens 1)) :y)
          ((< h-factor (varef lens 2)) :z)
          ;; ((< (* l-factor (tref lens-0 0 0)) (tref lens 0 0)) :x)
          ;; ((< (* l-factor (tref lens-0 1 0)) (tref lens 1 0)) :y)
          ;; ((< (* l-factor (tref lens-0 2 0)) (tref lens 2 0)) :z)
          ;; ((< l-factor (/ (tref lens 0 0) (tref lens-0 0 0))) t)
          ;; ((< l-factor (/ (tref lens 1 0) (tref lens-0 1 0))) t)
                                        ;((< 2.0d0 (tref def 1 1)) t)
                                        ;((> 1.5d0 (tref def 0 0)) t)
          (t nil)
          )))))

(defun single-particle-criteria (mesh mp)
  "Criteria for checking if material point is unconnected to another MP"
  (let ((svp-sum 0d0)
        (alone t))
    (iterate-over-neighbours
      mesh mp
      (lambda (mesh mp node svp grads fsvp fgrads)
        (declare
          (cl-mpm/mesh::node node)
          (cl-mpm/particle:particle mp)
          (type double-float svp))
        (with-accessors ((node-svp cl-mpm/mesh::node-svp-sum)
                         (node-active cl-mpm/mesh:node-active)
                         ) node
          (when node-active
            (incf svp-sum node-svp)))))
    (and (< svp-sum 2d0) (not (= svp-sum 0d0)))
    ;; (setf alone t)
    ;; alone
    ))
(defun gimp-removal-criteria (mp h)
  "Criteria for removal of gimp mps based on domain length"
  (with-accessors ((lens cl-mpm/particle::mp-domain-size))
      mp
    (declare (double-float h))
    (let ((h-factor (* 0.7d0 h))
          (aspect 0.01d0))
      (cond
        ((< h-factor (the double-float (varef lens 0))) :x)
        ((< h-factor (the double-float (varef lens 1))) :y)
        ((< h-factor (the double-float (varef lens 2))) :z)
        ((> aspect (/ (the double-float (varef lens 0))
                      (the double-float (varef lens 1))
                      )) :x)
        ((> aspect (/ (the double-float (varef lens 1))
                      (the double-float (varef lens 0))
                      )) :x)
        (t nil)
        ))))
(defun copy-particle (original &rest initargs &key &allow-other-keys)
  "Function for copying a particle, allowing for re-initialising members"
  (let* ((class (class-of original))
         (copy (allocate-instance class)))
    (dolist (slot (mapcar #'sb-mop:slot-definition-name (sb-mop:class-slots class)))
      (when (slot-boundp original slot)
        (cond
          ((typep (slot-value original slot) 'magicl::abstract-tensor)
           (setf (slot-value copy slot)
                 (magicl:scale (slot-value original slot) 1d0)))
          (t (setf (slot-value copy slot)
                  (slot-value original slot))))))
    (apply #'reinitialize-instance copy initargs)))

(defun split-vector (mp split-vec)
  "Helper macro for single splitting along cartesian directions "
  (with-accessors ((lens cl-mpm/particle::mp-domain-size)
                   (lens-0 cl-mpm/particle::mp-domain-size-0)
                   (mass cl-mpm/particle::mp-mass)
                   (pos cl-mpm/particle:mp-position)
                   (volume cl-mpm/particle:mp-volume)
                   (true-domain cl-mpm/particle::mp-true-domain)
                   )
      mp
    (let ((new-size (vector-from-list (list 1d0 1d0 1d0)))
          (new-size-0 (vector-from-list (list 1d0 1d0 1d0)))
          (pos-offset (vector-zeros))
          (new-split-depth (+ (cl-mpm/particle::mp-split-depth mp) 1))
          (new-domain nil)
          (vec-scaler (cl-mpm/fastmaths:fast-scale-vector (cl-mpm/fastmaths:norm split-vec) 0.5d0))
          )
      ;; (break)
      ;; (setf new-size (cl-mpm/fastmaths:fast-.- new-size vec-scaler))
      (setf new-size-0 (cl-mpm/fastmaths:fast-.- new-size-0 vec-scaler))
      (let ((domain-scaler (magicl:eye 3)))
        (setf (tref domain-scaler 0 0) (- 1d0 (abs (varef vec-scaler 0))))
        (setf (tref domain-scaler 1 1) (- 1d0 (abs (varef vec-scaler 1))))
        (setf (tref domain-scaler 2 2) (- 1d0 (abs (varef vec-scaler 2))))
        (setf new-domain (magicl:@ domain-scaler true-domain)))

      (setf pos-offset (magicl:@
                        true-domain
                        (cl-mpm/fastmaths::fast-scale-vector vec-scaler 0.5d0)))

      (setf
       (varef new-size 0)
       (cl-mpm/fastmaths:mag
        (cl-mpm/fastmaths::fast-@-matrix-vector
         new-domain
         (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0))
         ))
       (varef new-size 1)
       (cl-mpm/fastmaths:mag
        (cl-mpm/fastmaths::fast-@-matrix-vector
         new-domain
         (cl-mpm/utils:vector-from-list (list 0d0 1d0 0d0)))
        ))
      ;; (cl-mpm/fastmaths::fast-.* lens new-size new-size)
      (cl-mpm/fastmaths::fast-.* lens-0 new-size-0 new-size-0)
      ;; (cl-mpm/fastmaths::fast-.* lens pos-offset pos-offset)
      ;; (setf pos-offset (magicl:@ true-domain vec-scaler))
      ;; (break)
      (list
       (copy-particle mp
                      :mass (/ mass 2)
                      :volume (/ volume 2)
                      :size (cl-mpm/utils::vector-copy new-size)
                      :size-0 (cl-mpm/utils::vector-copy new-size-0)
                      :position (cl-mpm/fastmaths::fast-.+-vector pos pos-offset)
                      :nc (make-array 8 :fill-pointer 0 :element-type 'node-cache)
                      :split-depth new-split-depth
                      :true-domain (cl-mpm/utils:matrix-copy new-domain)
                      )
       (copy-particle mp
                      :mass (/ mass 2)
                      :volume (/ volume 2)
                      :size (cl-mpm/utils::vector-copy new-size)
                      :size-0 (cl-mpm/utils::vector-copy new-size-0)
                      :position (magicl:.- pos pos-offset)
                      :nc (make-array 8 :fill-pointer 0 :element-type 'node-cache)
                      :split-depth new-split-depth
                      :true-domain (cl-mpm/utils:matrix-copy new-domain)
                      )))))

(defmacro split-linear (dir direction dimension)
  "Helper macro for single splitting along cartesian directions "
  `((eq ,dir ,direction)
    (let ((new-size (vector-from-list (list 1d0 1d0 1d0)))
          (new-size-0 (vector-from-list (list 1d0 1d0 1d0)))
          (pos-offset (vector-zeros))
          (new-split-depth (+ (cl-mpm/particle::mp-split-depth mp) 1))
          (true-domain (cl-mpm/particle::mp-true-domain mp))
          (new-domain nil))
      (setf (tref new-size ,dimension 0) 0.5d0)
      (setf (tref new-size-0 ,dimension 0) 0.5d0)
      (setf (tref pos-offset ,dimension 0) 0.25d0)
      (let ((domain-scaler (magicl:eye 3)))
        (setf (tref domain-scaler ,dimension ,dimension) 0.5d0)
        (setf new-domain (magicl:@ true-domain domain-scaler)))
      (cl-mpm/fastmaths::fast-.* lens new-size new-size)
      (cl-mpm/fastmaths::fast-.* lens-0 new-size-0 new-size-0)
      (cl-mpm/fastmaths::fast-.* lens pos-offset pos-offset)
      (list
       (copy-particle mp
                      :mass (/ mass 2)
                      :volume (/ volume 2)
                      :size (cl-mpm/utils::vector-copy new-size)
                      :size-0 (cl-mpm/utils::vector-copy new-size-0)
                      :position (cl-mpm/fastmaths::fast-.+-vector pos pos-offset)
                      :nc (make-array 8 :fill-pointer 0 :element-type 'node-cache)
                      :split-depth new-split-depth
                      :true-domain (cl-mpm/utils:matrix-copy new-domain)
                       )
       (copy-particle mp
                      :mass (/ mass 2)
                      :volume (/ volume 2)
                      :size (cl-mpm/utils::vector-copy new-size)
                      :size-0 (cl-mpm/utils::vector-copy new-size-0)
                      :position (cl-mpm/fastmaths::fast-.--vector pos pos-offset)
                      :nc (make-array 8 :fill-pointer 0 :element-type 'node-cache)
                      :split-depth new-split-depth
                      :true-domain (cl-mpm/utils:matrix-copy new-domain)
                      )))))
(defmacro split-cases (direction)
  "Another helper macro for splitting mps"
  `(cond
     ,(macroexpand-1 '(split-linear direction :x 0))
     ,(macroexpand-1 '(split-linear direction :y 1))
     ,(macroexpand-1 '(split-linear direction :z 2))
     (t nil)
     ))
(defun split-mp (mp h direction)
  "Function to split an mp across a direction
   Directions should be :x,:y,:z "
  (with-accessors ((def cl-mpm/particle:mp-deformation-gradient)
                   (lens cl-mpm/particle::mp-domain-size)
                   (lens-0 cl-mpm/particle::mp-domain-size-0)
                   (pos cl-mpm/particle:mp-position)
                   (mass cl-mpm/particle:mp-mass)
                   (volume cl-mpm/particle:mp-volume)
                   ;; (volume cl-mpm/particle::mp-volume-0)
                   )
      mp
    (split-cases direction)))

(defun split-mps-eigenvalue (sim)
  (let* ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (crit (* h 0.25d0)))
    (cl-mpm::split-mps-vector
     sim
     (lambda (mp)
       (let ((split-dir nil))
         (multiple-value-bind (l v)
             (cl-mpm/utils:eig (cl-mpm/particle::mp-true-domain mp))
           (loop ;for i from 0 to 2
                 for lv in l
                 for i from 0
                 while (not split-dir)
                 do
                    (progn
                      ;; (pprint lv)
                      (when (> (abs lv) crit)
                        ;; (break)
                        (setf split-dir (magicl:column v i))
                        ))))
         split-dir)))))

(defun split-mps (sim)
  "Split mps that match the split-criteria"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (mps-to-split (remove-if-not (lambda (mp) (split-criteria mp h)) mps))
           (split-direction (map 'list (lambda (mp) (split-criteria mp h)) mps-to-split)))
      (remove-mps-func sim (lambda (mp) (split-criteria mp h)))
      (loop for mp across mps-to-split
            for direction in split-direction
            do (loop for new-mp in (split-mp mp h direction)
                     do (sim-add-mp sim new-mp)))
      )))

(defun split-mps-criteria (sim criteria)
  "Split mps that fail an arbritary criteria"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (mps-to-split (remove-if-not (lambda (mp) (funcall criteria mp h)) mps))
           (split-direction (map 'list (lambda (mp) (funcall criteria mp h)) mps-to-split)))
      (remove-mps-func sim (lambda (mp) (funcall criteria mp h)))
      (loop for mp across mps-to-split
            for direction in split-direction
            do (loop for new-mp in (split-mp mp h direction)
                     do (sim-add-mp sim new-mp))))))

(defun split-mps-vector (sim criteria)
  "Split mps that fail an arbritary criteria"
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (split-directions (lparallel:pmap 'vector criteria mps))
           (mps-to-split (delete-if-not #'identity (lparallel:pmap 'vector (lambda (mp crit) (when crit mp)) mps split-directions)))
           ;; (mps-to-split (remove-if-not (lambda (mp) (funcall criteria mp h)) mps))
           ;; (split-direction (map 'list (lambda (mp) (funcall criteria mp h)) mps-to-split))
           )
      (remove-mps-func sim (lambda (mp) (position mp mps-to-split)))
      (loop for mp across mps-to-split
            for direction across (remove-if-not #'identity split-directions)
            do (loop for new-mp in (split-vector mp direction)
                     do (sim-add-mp sim new-mp))))))

(defgeneric calculate-min-dt (sim)
  (:documentation "A function for calculating an approximate stable timestep"))
(defmethod calculate-min-dt ((sim mpm-sim))
  "Estimate minimum p-wave modulus"
  (with-accessors ((mesh cl-mpm:sim-mesh)
                   (mass-scale cl-mpm::sim-mass-scale))
      sim
    (let ((inner-factor most-positive-double-float))
      (declare (double-float inner-factor mass-scale))
      (iterate-over-nodes-serial
       mesh
       (lambda (node)
         (with-accessors ((node-active  cl-mpm/mesh:node-active)
                          (pmod cl-mpm/mesh::node-pwave)
                          (mass cl-mpm/mesh::node-mass)
                          (svp-sum cl-mpm/mesh::node-svp-sum)
                          (vol cl-mpm/mesh::node-volume)
                          (vel cl-mpm/mesh::node-velocity)
                          ) node
           (when (and node-active
                      (> vol 0d0)
                      (> pmod 0d0)
                      (> svp-sum 0d0))
             (let ((nf (+ (/ mass (* vol (+ (/ pmod svp-sum) ;; (* svp-sum (cl-mpm/fastmaths::mag-squared vel))
                                            )))))) 
                 (when (< nf inner-factor)
                   (setf inner-factor nf)))))))
      (if (< inner-factor most-positive-double-float)
          (* (sqrt mass-scale) (sqrt inner-factor) (cl-mpm/mesh:mesh-resolution mesh))
          (cl-mpm:sim-dt sim)))))
(defun calculate-adaptive-time (sim target-time &key (dt-scale 1d0))
  "Given that the system has previously been mapped to, caluclate an estimated dt and substep for a given target time
The value dt-scale allows for the estimated dt to be scaled up or down
This modifies the dt of the simulation in the process
"
  (let* ((dt-e (* dt-scale (calculate-min-dt sim)))
         (substeps-e (floor target-time dt-e)))
    ;; (format t "CFL dt estimate: ~f~%" dt-e)
    ;; (format t "CFL step count estimate: ~D~%" substeps-e)
    (setf (cl-mpm:sim-dt sim) dt-e)
    ;; (setf substeps substeps-e)
    (values dt-e substeps-e)
    )
  )


(defgeneric add-mps (sim mps-array))
(defmethod add-mps (sim mps-array)
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
      (if (> (length mps) 0)
          (progn
            (loop for mp across mps-array
                  do (sim-add-mp sim mp)))
          (progn
            (setf (cl-mpm:sim-mps sim) (make-array (length mps-array) :adjustable t :fill-pointer 0))
            (loop for mp across mps-array
                  do (sim-add-mp sim mp))))))

(defun add-bcs (sim bcs-array)
  "Add nodal essential bcs"
  (with-accessors ((bcs cl-mpm:sim-bcs))
      sim
    (if (> (length bcs) 0)
        (progn
          (loop for bc across bcs-array
                do (vector-push-extend bc bcs)))
        (setf bcs bcs-array))))
(defun add-bcs-force-list (sim new-bcs)
  "Add bcs that apply forces, ordered in a FILO"
  (with-accessors ((bcs-force-list cl-mpm:sim-bcs-force-list))
      sim
    ;; (when (= (length bcs-force-list 0)))
    (typecase new-bcs
      (list
       ;;Common case, we have a list, or arrays of BCs, running bcs in serial down the list, but in parrallel over the array
       (setf bcs-force-list (nconc new-bcs bcs-force-list)))
      (array
       ;;We have an array of bcs, that we want to add to the list -> add it to the back in its own list
       (setf bcs-force-list (nconc (list new-bcs) bcs-force-list)))
      (cl-mpm/bc::bc
       ;;Uncommon case we have a singluar bc we wish to add, so double wrap it
       (setf bcs-force-list (nconc (list (cl-mpm/bc:make-bcs-from-list (list new-bcs))) bcs-force-list)))
      (t (error "BCs must be a list, array or singluar bc")))))

(defun get-mp (sim index)
  (aref (sim-mps sim) index))
