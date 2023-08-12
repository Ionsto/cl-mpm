(defpackage :cl-mpm/mpi
  (:use
   :cl
   :cl-mpm
   :cl-mpm/particle
   :cl-mpm/mesh
   :cl-mpm/utils
   :cl-mpm/fastmath
   )
  (:import-from
   :magicl tref .+ .-
   )
  )
(declaim (optimize (debug 0) (safety 0) (speed 3)))
                                        ;    #:make-shape-function
(in-package :cl-mpm/mpi)


(declaim (inline calculate-strain-rate)
         (ftype (function (cl-mpm/particle:particle double-float)) calculate-strain-rate))
(defun calculate-strain-rate (mp dt)
  (declare (cl-mpm/particle:particle mp) (double-float dt))
  (with-accessors ((strain-rate cl-mpm/particle:mp-strain-rate)
                   (vorticity cl-mpm/particle:mp-vorticity)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (jfbar cl-mpm/particle::mp-j-fbar)
                   (velocity-rate cl-mpm/particle::mp-velocity-rate)
                   ) mp
        (progn
          (magicl:scale! strain-rate 0d0)
          (magicl:scale! vorticity 0d0)
          (magicl:scale! stretch-tensor 0d0)
          (let (;(stretch-dsvp (magicl:zeros '(4 2)))
                #+cl-mpm-fbar (stretch-tensor-fbar (magicl:zeros '(2 2) :type 'double-float))
                ;(v-s (magicl:zeros '(2 2)))
                (v-s (matrix-zeros))
                (stretch-dsvp (stretch-dsvp-zeros))
                )
            (declare (magicl:matrix/double-float stretch-dsvp v-s)
                     (dynamic-extent stretch-dsvp v-s))
            (cl-mpm::iterate-over-neighbours-cached
             nil mp
             (lambda (mesh mp node svp grads fsvp fgrads)
               (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
                                (node-active cl-mpm/mesh:node-active))
                   node
                 (declare (double-float))
                 (when node-active
                   (magicl.simd::.+-simd
                    stretch-tensor
                    (cl-mpm/utils::voight-to-stretch-prealloc
                     (magicl:@ (cl-mpm/shape-function::assemble-dstretch-2d-prealloc
                                grads stretch-dsvp) node-vel) v-s)
                    stretch-tensor)
                   #+cl-mpm-fbar (magicl.simd::.+-simd
                                  stretch-tensor-fbar
                                  (cl-mpm/utils::voight-to-stretch-prealloc
                                   (magicl:@ (cl-mpm/shape-function::assemble-dstretch-2d-prealloc
                                              fgrads stretch-dsvp) node-vel) v-s)
                                  stretch-tensor-fbar)
                   ;; (mult (cl-mpm/shape-function::assemble-dsvp-2d grads) node-vel strain-rate)
                   ;; (mult (cl-mpm/shape-function::assemble-vorticity-2d grads) node-vel vorticity)
                   ))))
            #+cl-mpm-fbar (setf jfbar (magicl:det (magicl.simd::.+-simd (magicl:eye 2) stretch-tensor-fbar)))
            )
          #+cl-mpm-fbar (when (<= jfbar 0d0)
            (error "FBAR volume non-positive"))
            (cl-mpm/fastmath::stretch-to-sym stretch-tensor strain-rate)
            (cl-mpm/fastmath::stretch-to-skew stretch-tensor vorticity)
            (aops:copy-into (magicl::matrix/double-float-storage velocity-rate) (magicl::matrix/double-float-storage strain-rate))
            ;; (setf velocity-rate (magicl:scale strain-rate 1d0))
            (magicl:scale! stretch-tensor dt)
            (magicl:scale! strain-rate dt)
            (magicl:scale! vorticity dt)
            )))

(declaim (inline update-strain-kirchoff))
(declaim (ftype (function (cl-mpm/particle:particle
                           double-float) (values))
               update-strain-kirchoff))
(defun update-strain-kirchoff (mp dt)
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (volume-0 cl-mpm/particle::mp-volume-0)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
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
                   domain
                   ))
    (progn
      (let ((df (cl-mpm::calculate-df mp)))
        (progn
          ;; (magicl:mult df def :target def)
          (setf def (magicl:@ df def))
          (let ((initial-strain (magicl:scale strain 1d0))
                )
            (multiple-value-bind (l v) (cl-mpm/utils::eig (voigt-to-matrix strain))
              (let (;(trail-lgs (cl-mpm/utils::matrix-zeros))
                    (trial-lgs (magicl:@ df
                                         v
                                         (cl-mpm/utils::matrix-from-list
                                          (list
                                           (the double-float (exp (* 2d0 (the double-float (nth 0 l)))))
                                           0d0 0d0
                                           (the double-float (exp (* 2d0 (the double-float (nth 1 l)))))))
                                         (magicl:transpose v)
                                         (magicl:transpose df)))
                    )
                (multiple-value-bind (lf vf) (cl-mpm/utils::eig
                                              (magicl:scale! (magicl:.+ trial-lgs (magicl:transpose trial-lgs)) 0.5d0))
                  (setf strain (magicl:scale!
                                (matrix-to-voigt
                                 (magicl:@
                                  vf
                                  (cl-mpm/utils::matrix-from-list
                                   (list
                                    (the double-float (log (the double-float (nth 0 lf))))
                                    0d0 0d0
                                    (the double-float (log (the double-float (nth 1 lf))))))
                                  (magicl:transpose vf)))
                                0.5d0))
                  ;; (print strain)
                  ;; (setf (magicl:tref strain 2 0) (* 2d0 (the double-float (magicl:tref strain 2 0))))
                  )
                ))
            (magicl:.- initial-strain strain initial-strain)
            (setf eng-strain-rate initial-strain)
            (magicl:scale! eng-strain-rate (/ 1d0 dt))

            ;; (setf eng-strain-rate (magicl:scale! (magicl:.- strain initial-strain) (/ 1d0 dt)))
            )

          ;;Post multiply to turn to eng strain
          ;; (setf volume (* volume (magicl:det df)))
          (setf volume (* volume-0 (magicl:det def)))
          (when (<= volume 0d0)
            (error "Negative volume"))
          ;;Stretch rate update
          (let ((F (cl-mpm/utils::matrix-zeros)))
            (magicl:mult def def :target F :transb :t)
            (multiple-value-bind (l v) (cl-mpm/utils::eig F)
              (let ((stretch
                      (magicl:@
                       v
                       (cl-mpm/utils::matrix-from-list
                        (list (the double-float (sqrt (the double-float (nth 0 l))))
                              0d0 0d0
                              (the double-float (sqrt (the double-float (nth 1 l))))))
                       (magicl:transpose v)))
                    )
                (declare (type magicl:matrix/double-float stretch))
                (setf (tref domain 0 0) (* (the double-float (tref domain-0 0 0))
                                           (the double-float (tref stretch 0 0))))
                (setf (tref domain 1 0) (* (the double-float (tref domain-0 1 0))
                                           (the double-float (tref stretch 1 1))))
                )))
          ;; (update-domain-corner mp dt)
          ;; (setf (tref domain 0 0) (* (the double-float (tref domain 0 0))
          ;;                             (the double-float (tref df 0 0))))
          ;; (setf (tref domain 1 0) (* (the double-float (tref domain 1 0))
          ;;                             (the double-float (tref df 1 1))))
          )
          )))
  (values))


(declaim (inline update-stress-mp)
         (ftype (function (cl-mpm/particle:particle double-float) (values)) update-stress-mp))
(defun update-stress-mp (mp dt)
  ;; (declare ((cl-mpm/particle::particle mp)
  ;;           (double-float dt)))
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
            (calculate-strain-rate mp dt)

            ;;; Turn cauchy stress to kirchoff
            (setf stress stress-kirchoff)

            ;;; Update our strains
            (update-strain-kirchoff mp dt)

            ;;;Update our kirchoff stress with constitutive model
            (setf stress-kirchoff (cl-mpm/particle:constitutive-model mp strain dt))

            ;; Check volume constraint!
            (when (<= volume 0d0)
              (error "Negative volume"))
            ;; Turn kirchoff stress to cauchy
            (setf stress (magicl:scale stress-kirchoff (/ 1.0d0 (the double-float (magicl:det def)))))
            ))))

(defclass mpm-sim-mpi-stress (cl-mpm/damage::mpm-sim-damage)
  ()
  (:documentation "Damage sim with only stress update on mpi"))
(defmethod cl-mpm::update-sim ((sim mpm-sim-mpi-stress))
  (with-slots ((mesh cl-mpm::mesh)
               (mps  cl-mpm::mps)
               (bcs  cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (bcs-force-list cl-mpm::bcs-force-list)
               (dt cl-mpm::dt)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (nonlocal-damage cl-mpm::nonlocal-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               )
                sim
    (declare (type double-float mass-filter))
                (progn
                    (cl-mpm::reset-grid mesh)
                    (cl-mpm::p2g mesh mps)
                    (when (> mass-filter 0d0)
                      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
                    (cl-mpm::update-node-kinematics mesh dt )
                    (cl-mpm::apply-bcs mesh bcs dt)
                    ;; ;(cl-mpm::update-stress mesh mps dt)
                    ;; (loop for mp across mps
                    ;;       do (cl-mpm/mpi::update-stress-mp mp dt))

                    (let ((dt-e dt))
                      (lfarm:broadcast-task (lambda ()
                                              (progn
                                                (setf *global-dt* dt-e)
                                                t))))
                    (lfarm:pmap-into (cl-mpm::sim-mps sim)
                                     'uls
                                     (cl-mpm::sim-mps sim)
                                     )
                    (when enable-damage
                     (cl-mpm/damage::calculate-damage mesh
                                                      mps
                                                      dt
                                                      50d0
                                                      nonlocal-damage
                                                      ))
                    ;Map forces onto nodes
                    (cl-mpm::p2g-force mesh mps)
                    ;(cl-mpm::apply-bcs mesh bcs-force dt)
                    (loop for bcs-f in bcs-force-list
                          do
                             (cl-mpm::apply-bcs mesh bcs-f dt))
                    (cl-mpm::update-node-forces mesh (cl-mpm::sim-damping-factor sim) dt (cl-mpm::sim-mass-scale sim))
                    ;Reapply velocity BCs
                    (cl-mpm::apply-bcs mesh bcs dt)
                    ;Also updates mps inline
                    (cl-mpm::g2p mesh mps dt)

                    (when remove-damage
                      (cl-mpm::remove-material-damaged sim))
                    (when split
                      (cl-mpm::split-mps sim))
                    (cl-mpm::check-mps mps)
                    )))


(defvar *mutex-code* (cl-store:register-code 110 'sb-thread:mutex))
(cl-store:defstore-cl-store (obj sb-thread:mutex stream)
   (cl-store:output-type-code *mutex-code* stream))
(cl-store:defrestore-cl-store (sb-thread:mutex stream)
   (sb-thread:make-mutex))

;; (defvar *node-cache-code* (cl-store:register-code 112 'cl-mpm/particle::node-cache))
;; (cl-store:defstore-cl-store (obj cl-mpm/particle::node-cache stream)
;;   (cl-store:output-type-code *node-cache-code* stream)
;;   )
;; (cl-store:defrestore-cl-store (sb-thread:mutex stream)
;;   ;(sb-thread:make-mutex)
;;   )


(defvar *mesh-code* (cl-store:register-code 111 'cl-mpm/mesh::mesh))
(cl-store:defstore-cl-store (obj cl-mpm/mesh::mesh stream)
    (cl-store:output-type-code *mesh-code* stream)
  (cl-store:store-object (cl-mpm/mesh::mesh-nd obj) stream)
  (cl-store:store-object (cl-mpm/mesh::mesh-count obj) stream)
  (cl-store:store-object (cl-mpm/mesh::mesh-mesh-size obj) stream)
  (cl-store:store-object (cl-mpm/mesh::mesh-resolution obj) stream)
  (cl-store:store-object (cl-mpm/mesh::mesh-nodes obj) stream)
  (cl-store:store-object (cl-mpm/mesh::mesh-cells obj) stream)
  )

(cl-store:defrestore-cl-store (cl-mpm/mesh::mesh stream)
    (let ((obj (make-instance 'cl-mpm/mesh::mesh)))
      (setf (cl-mpm/mesh::mesh-nd obj) (cl-store:restore-object stream))
      (setf (cl-mpm/mesh::mesh-count obj) (cl-store:restore-object stream))
      (setf (cl-mpm/mesh::mesh-mesh-size obj)  (cl-store:restore-object stream))
      (setf (cl-mpm/mesh::mesh-resolution obj) (cl-store:restore-object stream))
      (setf (cl-mpm/mesh::mesh-nodes obj) (cl-store:restore-object stream))
      (setf (cl-mpm/mesh::mesh-cells obj) (cl-store:restore-object stream))
      obj))

(defvar *sim-code* (cl-store:register-code 112 'cl-mpm::mpm-sim))
(cl-store:defstore-cl-store (obj cl-mpm::mpm-sim stream)
    (cl-store:output-type-code *sim-code* stream)
  (cl-store:store-object (cl-mpm::sim-dt obj) stream)
  (cl-store:store-object (cl-mpm::sim-mesh obj) stream)
  )
(cl-store:defrestore-cl-store (cl-mpm::mpm-sim stream)
    (let ((obj (make-instance 'cl-mpm::mpm-sim)))
      (setf (cl-mpm::sim-dt obj) (cl-store:restore-object stream))
      (setf (cl-mpm::sim-mesh obj) (cl-store:restore-object stream))
      obj))

;; (defmethod cl-store:serializable-slots-using-class ((object t) (class cl-mpm/mesh::node))
;;   (delete 'cl-mpm/mesh::local-list (call-next-method) :key 'c2mop:slot-definition-name))

(defvar *node-code* (cl-store:register-code 113 'cl-mpm/mesh::node))
(cl-store:defstore-cl-store (obj cl-mpm/mesh::node stream)
  (cl-store:output-type-code *node-code* stream)
  (loop for slot in (cl-store:serializable-slots obj)
        do
           (cond
             ((eq 'cl-mpm/mesh::local-list (sb-mop:slot-definition-name slot))
              nil)
             (t
              (cl-store:store-object (sb-mop:slot-value-using-class 'cl-mpm/mesh::node obj slot) stream)))))
(cl-store:defrestore-cl-store (cl-mpm/mesh::node stream)
  (let ((obj (make-instance 'cl-mpm/mesh::node)))
    (loop for slot in (cl-store:serializable-slots obj)
          do
             (cond
               ((eq 'cl-mpm/mesh::local-list (sb-mop:slot-definition-name slot))
                nil)
               (t
                (setf (sb-mop:slot-value-using-class 'cl-mpm/mesh::node obj slot) (cl-store:restore-object stream)))))
  ;; (setf (cl-mpm::sim-dt obj) (cl-store:restore-object stream))
    obj))


(defun collect-servers (n)
  (let ((servers (with-open-file (s "lfarm_connections") (read s))))
    (format t "~S ~%" servers)
    (defparameter *open-servers* servers)
    (setf lfarm:*kernel* (lfarm:make-kernel servers))
    )
  (asdf:compile-system :magicl)
  (print "Broadcasting setup info")
  (lfarm:broadcast-task (lambda ()
                          (progn
                            ;; (asdf:compile-system :magicl :force t)
                            (ql:quickload :cl-mpm/examples/slump)
                            (ql:quickload :cl-mpm/mpi)
                            (in-package :cl-mpm/examples/slump)
                            (defparameter *local-sim* nil)
                            (defparameter *global-dt* 1d0)
                            (setf lparallel:*kernel* (lparallel:make-kernel 4))
                            t))))

(defun kill-servers ()
    (dolist (server *open-servers*)
      (lfarm-admin:end-server (first server) (second server))))

(defparameter *global-dt* 1d0)
(lfarm:deftask uls (mp)
  (update-stress-mp mp *global-dt*)
  (setf (fill-pointer (cl-mpm/particle::mp-cached-nodes mp)) 0)
  mp)
