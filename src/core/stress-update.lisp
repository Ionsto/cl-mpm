(in-package :cl-mpm)
;;Various types of stress update
(declaim #.cl-mpm/settings:*optimise-setting*)

(declaim (notinline calculate-strain-rate)
         (ftype (function (cl-mpm/mesh::mesh  cl-mpm/particle:particle double-float)) calculate-strain-rate))
(defun calculate-strain-rate (mesh mp dt)
  "Calculate the strain rate, stretch rate and vorticity"
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt)
           (optimize (speed 3) (safety 0)))
  (with-accessors ((stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (stretch-tensor-fbar cl-mpm/particle::mp-stretch-tensor-fbar)
                   ) mp
    (declare (magicl:matrix/double-float stretch-tensor stretch-tensor-fbar)
             (double-float dt))
        (progn
          (cl-mpm/fastmaths::fast-zero stretch-tensor)
          (cl-mpm/fastmaths::fast-zero stretch-tensor-fbar)
          (iterate-over-neighbours
           mesh mp
           (lambda (mesh mp node svp grads fsvp fgrads)
             (declare (ignore mp svp fsvp))
             (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
                              ;; (node-disp cl-mpm/mesh::node-displacment)
                              (node-active cl-mpm/mesh:node-active))
                 node
               (declare (magicl:matrix/double-float node-vel)
                        (boolean node-active))
               (when node-active
                 (cl-mpm/shape-function::@-combi-assemble-dstretch-3d grads node-vel stretch-tensor)
                 (cl-mpm/shape-function::@-combi-assemble-dstretch-3d fgrads node-vel stretch-tensor-fbar)))))
          (cl-mpm/fastmaths:fast-scale! stretch-tensor dt)
          (cl-mpm/fastmaths:fast-scale! stretch-tensor-fbar dt)
          )))

(declaim (notinline calculate-strain-rate)
         (ftype (function (cl-mpm/mesh::mesh  cl-mpm/particle:particle double-float)) calculate-strain-rate))
(defun calculate-strain-rate-disp (mesh mp dt)
  "Calculate the strain rate, stretch rate and vorticity"
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt)
           (optimize (speed 3) (safety 0)))
  (with-accessors ((stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (stretch-tensor-fbar cl-mpm/particle::mp-stretch-tensor-fbar)
                   ) mp
    (declare (magicl:matrix/double-float stretch-tensor stretch-tensor-fbar)
             (double-float dt))
        (progn
          (cl-mpm/fastmaths::fast-zero stretch-tensor)
          (cl-mpm/fastmaths::fast-zero stretch-tensor-fbar)
          (let ()
            (iterate-over-neighbours
             mesh mp
             (lambda (mesh mp node svp grads fsvp fgrads)
               (declare (ignore mp svp fsvp))
               (with-accessors ((node-disp cl-mpm/mesh::node-displacment)
                                (node-active cl-mpm/mesh:node-active))
                   node
                 (declare (magicl:matrix/double-float node-disp)
                          (boolean node-active))
                 (when node-active
                   (cl-mpm/shape-function::@-combi-assemble-dstretch-3d grads node-disp stretch-tensor)
                   (cl-mpm/shape-function::@-combi-assemble-dstretch-3d fgrads node-disp stretch-tensor-fbar)))))))))

;Could include this in p2g but idk


(defun calculate-strain-rate-nofbar (mesh mp dt)
  "Calculate the strain rate, stretch rate and vorticity"
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt)
           (optimize (speed 3) (safety 0)))
  (with-accessors ((stretch-tensor cl-mpm/particle::mp-stretch-tensor)) mp
    (declare (magicl:matrix/double-float stretch-tensor)
             (double-float dt))
        (progn
          (cl-mpm/fastmaths::fast-zero stretch-tensor)
          (iterate-over-neighbours
           mesh mp
           (lambda (mesh mp node svp grads fsvp fgrads)
             (declare (ignore mesh mp svp fsvp fgrads))
             (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
                              (node-active cl-mpm/mesh:node-active))
                 node
               (declare (magicl:matrix/double-float node-vel)
                        (boolean node-active))
               (when node-active
                 (cl-mpm/shape-function::@-combi-assemble-dstretch-3d grads node-vel stretch-tensor)))))
          (cl-mpm/fastmaths::fast-scale! stretch-tensor dt))))

(declaim (notinline calculate-strain-rate)
         (ftype (function (cl-mpm/mesh::mesh  cl-mpm/particle:particle double-float)) calculate-strain-rate))
(defun calculate-strain-rate-DR (mesh mp dt)
  "Calculate the strain rate, stretch rate and vorticity"
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt)
           (optimize (speed 3) (safety 0)))
  (with-accessors ((stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (stretch-tensor-fbar cl-mpm/particle::mp-stretch-tensor-fbar)
                   (df cl-mpm/particle::mp-deformation-gradient-increment)
                   ) mp
    (declare (magicl:matrix/double-float stretch-tensor stretch-tensor-fbar)
             (double-float dt))
        (progn
          (cl-mpm/fastmaths::fast-zero stretch-tensor)
          (cl-mpm/fastmaths::fast-zero stretch-tensor-fbar)
          (let ()
            (iterate-over-neighbours
             mesh mp
             (lambda (mesh mp node svp grads fsvp fgrads)
               (declare (ignore mp svp fsvp))
               (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
                                (node-disp cl-mpm/mesh::node-displacment)
                                (node-active cl-mpm/mesh:node-active))
                   node
                 (declare (magicl:matrix/double-float node-vel)
                          (boolean node-active))
                 (when node-active
                   (let* ((grads (cl-mpm::gradient-push-forwards grads df)))
                     (cl-mpm/shape-function::@-combi-assemble-dstretch-3d grads node-disp stretch-tensor)
                     (cl-mpm/shape-function::@-combi-assemble-dstretch-3d fgrads node-disp stretch-tensor-fbar))))))))))

(defun make-df (inc)
  (let ((df (cl-mpm/utils::matrix-from-list '(1d0 0d0 0d0
                                              0d0 1d0 0d0
                                              0d0 0d0 1d0))))
    (cl-mpm/fastmaths::fast-.+-matrix df inc df)
    df))

(declaim (notinline calculate-df)
         (ftype (function (cl-mpm/mesh::mesh
                           cl-mpm/particle:particle
                           boolean
                           &optional (or null magicl:matrix/double-float)
                           ) magicl:matrix/double-float)
                calculate-df))
(defun calculate-df (mesh mp fbar &optional (result nil))
  (with-accessors ((dstrain cl-mpm/particle::mp-strain-rate)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (stretch-tensor-fbar cl-mpm/particle::mp-stretch-tensor-fbar)
                   (jfbar cl-mpm/particle::mp-j-fbar)
                   (def cl-mpm/particle:mp-deformation-gradient)
                   (pos cl-mpm/particle:mp-position))
      mp
    (let* ((df (if result
                   (cl-mpm/fastmaths::matrix-reset-identity result)
                   (cl-mpm/utils:matrix-eye 1d0)))
           (dJ 1d0))
      (cl-mpm/fastmaths::fast-.+-matrix df stretch-tensor df)
      (setf dJ (cl-mpm/fastmaths:det-3x3 df))
      (when (< dJ 0d0)
        (error 'cl-mpm/errors:error-dF-negative))
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
              (let* ((df-fbar (cl-mpm/utils::matrix-eye 1d0))
                     (nd (cl-mpm/mesh::mesh-nd mesh)))
                (cl-mpm/fastmaths::fast-.+-matrix df-fbar stretch-tensor-fbar df-fbar)
                (when (< (cl-mpm/fastmaths:det-3x3 df-fbar) 0d0)
                  (error 'cl-mpm/errors:error-dF-negative)
                  )
                (cl-mpm/fastmaths::fast-scale!
                 df
                 (expt
                  (the double-float (/ (cl-mpm/fastmaths:det-3x3 df-fbar)
                                       (cl-mpm/fastmaths:det-3x3 df)))
                  (the double-float (/ 1d0 nd))))
                (when (= nd 2)
                  (setf (magicl:tref df 2 2) 1d0))))))
      (values df dJ))))

;;; Strain updates

(defun update-strain-linear (mesh mp dt fbar)
  "Linear small strain update algorithm"
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (volume-n cl-mpm/particle::mp-volume-n)
                   (volume-0 cl-mpm/particle::mp-volume-0)
                   (strain cl-mpm/particle:mp-strain)
                   (strain-n cl-mpm/particle:mp-strain-n)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (def-0    cl-mpm/particle::mp-deformation-gradient-0)
                   (df-inc    cl-mpm/particle::mp-deformation-gradient-increment)
                   (df-inc-inv    cl-mpm/particle::mp-deformation-gradient-increment-inverse)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   (delta-vol cl-mpm/particle::mp-delta-volume)
                   (strain-rate cl-mpm/particle:mp-strain-rate))
      mp
    (declare (double-float volume volume-0))
    (multiple-value-bind (df dj) (calculate-df mesh mp fbar)
      (progn
        (setf def (cl-mpm/fastmaths::fast-@-matrix-matrix df-inc def-0 def))
        ;; (cl-mpm/utils:voigt-copy-into strain-n strain)
        (cl-mpm/fastmaths:fast-.+
         strain-n
         (cl-mpm/utils::stretch-to-sym stretch-tensor)
         strain)
        (setf delta-vol dj)
        (setf volume (* volume-n (the double-float delta-vol)))
        (setf df-inc-inv (cl-mpm/fastmaths::fast-inv-3x3 df-inc))
        (when (<= volume 0d0)
          (error 'cl-mpm/errors:error-volume-negative))))))

(declaim (inline update-strain-kirchoff))
(declaim (ftype (function (cl-mpm/mesh::mesh
                           cl-mpm/particle:particle
                           double-float
                           boolean) (values))
               update-strain-kirchoff))


(defun update-strain-kirchoff (mesh mp dt fbar)
  "Finite strain kirchhoff strain update algorithm"
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (volume-n cl-mpm/particle::mp-volume-n)
                   (volume-0 cl-mpm/particle::mp-volume-0)
                   (strain cl-mpm/particle:mp-strain)
                   (strain-n cl-mpm/particle:mp-strain-n)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (df-inc    cl-mpm/particle::mp-deformation-gradient-increment)
                   (eng-strain-rate cl-mpm/particle::mp-eng-strain-rate)
                   ) mp
    (declare (type double-float volume))
    (progn
      (multiple-value-bind (df dj) (calculate-df mesh mp fbar)
        (progn
          (cl-mpm/utils:matrix-copy-into df df-inc)
          (let ((temp (cl-mpm/utils:matrix-zeros)))
            (cl-mpm/fastmaths::fast-@-matrix-matrix df def temp)
            (cl-mpm/utils:matrix-copy-into temp def))
          (cl-mpm/utils:voigt-copy-into strain strain-n)
          (cl-mpm/ext:kirchoff-update strain df)
          (setf volume (* volume-n (the double-float dj)))
          ;; (setf volume (* volume-0 (magicl:det def)))
          (when (<= volume 0d0)
            (error 'cl-mpm/errors:error-volume-negative))))))
  (values))

(defun update-strain-kirchoff-dynamic-relaxation (mesh mp dt fbar)
  "Finite strain kirchhoff strain update algorithm"
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (volume-n cl-mpm/particle::mp-volume-n)
                   (volume-0 cl-mpm/particle::mp-volume-0)
                   (strain cl-mpm/particle:mp-strain)
                   (strain-n cl-mpm/particle:mp-strain-n)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (def-0    cl-mpm/particle::mp-deformation-gradient-0)
                   (df-inc    cl-mpm/particle::mp-deformation-gradient-increment)
                   (delta-vol cl-mpm/particle::mp-delta-volume)
                   (df-inc-inv    cl-mpm/particle::mp-deformation-gradient-increment-inverse)
                   ) mp
    (declare (type double-float volume))
    (progn
      (multiple-value-bind (df dj) (calculate-df mesh mp fbar df-inc)
        (progn
          ;; (cl-mpm/utils:matrix-copy-into df df-inc)
          (setf def (cl-mpm/fastmaths::fast-@-matrix-matrix df-inc def-0 def))
          ;; (setf def (cl-mpm/fastmaths::fast-@-matrix-matrix df-inc def-0))
          (cl-mpm/utils:voigt-copy-into strain-n strain)
          (cl-mpm/ext:kirchoff-update strain df-inc)
          (setf delta-vol dj)
          (setf volume (* volume-n
                          delta-vol
                          ;; (the double-float (cl-mpm/fastmaths:det-3x3 df))
                          ))
          (setf df-inc-inv (cl-mpm/fastmaths::fast-inv-3x3 df-inc))
          (when (<= volume 0d0)
            (error 'cl-mpm/errors:error-volume-negative))))))
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
        (update-strain-kirchoff mesh mp dt fbar)
        (cl-mpm/utils::voigt-copy-into (cl-mpm/particle:constitutive-model mp strain dt) stress-kirchoff)
        (cl-mpm/utils::voigt-copy-into stress-kirchoff stress)
        ;; (setf stress-kirchoff (cl-mpm/utils::voigt-copy stress))
        (cl-mpm/fastmaths::fast-scale! stress (/ 1.0d0 (the double-float (cl-mpm/fastmaths:det-3x3 def))))
        ))))

(defun update-stress-kirchoff-dynamic-relaxation (mesh mp dt fbar)
  "Update stress for a single mp"
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt))
  (with-accessors ((stress cl-mpm/particle:mp-stress)
                   (stress-kirchoff cl-mpm/particle::mp-stress-kirchoff)
                   (volume cl-mpm/particle:mp-volume)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (def-0    cl-mpm/particle::mp-deformation-gradient-0)
                   (df-inc    cl-mpm/particle::mp-deformation-gradient-increment)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   ) mp
    (declare (magicl:matrix/double-float stress stress-kirchoff strain def strain-rate)
             (double-float volume))
    (progn
      (progn
        (calculate-strain-rate-disp mesh mp dt)
        ;; ;; Update our strains
        (update-strain-kirchoff-dynamic-relaxation mesh mp dt fbar)
        (cl-mpm/utils::voigt-copy-into (cl-mpm/particle:constitutive-model mp strain dt) stress-kirchoff)
        ;; ;; Turn kirchoff stress to cauchy
        (cl-mpm/utils::voigt-copy-into stress-kirchoff stress)
        (cl-mpm/fastmaths::fast-scale! stress (/ 1.0d0 (the double-float (cl-mpm/fastmaths:det-3x3 def))))
        ))))

(defun update-strain-kirchoff-dynamic-relaxation-incremental (mesh mp dt fbar)
  "Finite strain kirchhoff strain update algorithm"
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (volume-0 cl-mpm/particle::mp-volume-0)
                   (strain cl-mpm/particle:mp-strain)
                   (strain-n cl-mpm/particle:mp-strain-n)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (def-0    cl-mpm/particle::mp-deformation-gradient-0)
                   (df-inc    cl-mpm/particle::mp-deformation-gradient-increment)
                   (df-inc-inv    cl-mpm/particle::mp-deformation-gradient-increment-inverse)
                   ) mp
    (declare (type double-float volume))
    (progn
      (multiple-value-bind (df dj) (calculate-df mesh mp fbar)
        (progn
          (setf df-inc (cl-mpm/fastmaths::fast-@-matrix-matrix df df-inc))
          (setf def (cl-mpm/fastmaths::fast-@-matrix-matrix df-inc def-0 def))
          (cl-mpm/utils:voigt-copy-into strain-n strain)
          (cl-mpm/ext:kirchoff-update strain df-inc)
          (setf volume (* volume-0 (the double-float (cl-mpm/fastmaths:det-3x3 def))))
          (setf df-inc-inv (cl-mpm/fastmaths::fast-inv-3x3 df-inc))
          (when (<= volume 0d0)
            (error 'cl-mpm/errors::error-volume-negative))))))
  (values))



(defun update-stress-kirchoff-dynamic-relaxation-incremental (mesh mp dt fbar)
  "Update stress for a single mp"
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt))
  (with-accessors ((stress cl-mpm/particle:mp-stress)
                   (stress-kirchoff cl-mpm/particle::mp-stress-kirchoff)
                   (volume cl-mpm/particle:mp-volume)
                   (strain cl-mpm/particle:mp-strain)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (def-0    cl-mpm/particle::mp-deformation-gradient-0)
                   (df-inc    cl-mpm/particle::mp-deformation-gradient-increment)
                   (strain-rate cl-mpm/particle:mp-strain-rate)
                   (stretch-tensor cl-mpm/particle::mp-stretch-tensor)
                   ) mp
    (declare (magicl:matrix/double-float stress stress-kirchoff strain def strain-rate)
             (double-float volume))
    (progn
      (progn
        ;; (calculate-strain-rate-dr mesh mp dt)
        ;; (calculate-strain-rate mesh mp dt)
        (calculate-strain-rate mesh mp dt)
        ;; Turn cauchy stress to kirchoff
        (cl-mpm/utils::voigt-copy-into stress-kirchoff stress)
        ;; Update our strains
        (update-strain-kirchoff-dynamic-relaxation-incremental mesh mp dt fbar)
        (cl-mpm/utils::voigt-copy-into (cl-mpm/particle:constitutive-model mp strain dt) stress-kirchoff)
        ;; Turn kirchoff stress to cauchy
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
          (error 'cl-mpm/errors:error-volume-negative))
        ;; Turn kirchoff stress to cauchy
        (cl-mpm/utils::voigt-copy-into stress-kirchoff stress)
        (cl-mpm/fastmaths::fast-scale! stress (/ 1.0d0 (the double-float (cl-mpm/fastmaths:det-3x3 def))))
        ))))


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
        (calculate-strain-rate-disp mesh mp dt)
        ;; Update our strains
        (update-strain-linear mesh mp dt fbar)
        ;; Update our kirchoff stress with constitutive model
        (setf stress (cl-mpm/particle:constitutive-model mp strain dt))
        ;; Check volume constraint!
        (when (<= volume 0d0)
          (error 'cl-mpm/errors:error-volume-negative))))))


(defun update-stress-kirchoff-p (mesh mp dt fbar)
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
        (calculate-strain-rate-nofbar mesh mp dt)
        ;; (calculate-strain-rate mesh mp dt)
        ;; Turn cauchy stress to kirchoff
        (cl-mpm/utils::voigt-copy-into stress-kirchoff stress)
        ;; Update our strains
        (update-strain-kirchoff mesh mp dt fbar)
        ;; Update our kirchoff stress with constitutive model
        (cl-mpm/utils::voigt-copy-into (cl-mpm/particle:constitutive-model mp strain dt) stress-kirchoff)
        ;; Turn kirchoff stress to cauchy
        (cl-mpm/utils::voigt-copy-into stress-kirchoff stress)
        (cl-mpm/fastmaths::fast-scale! stress (/ 1.0d0 (the double-float (cl-mpm/fastmaths:det-3x3 def))))
        ))))

(defun calculate-strain-rate-p (mesh mp dt)
  "Calculate the strain rate, stretch rate and vorticity"
  (declare (cl-mpm/mesh::mesh mesh) (cl-mpm/particle:particle mp) (double-float dt)
           (optimize (speed 3) (safety 0)))
  (declare (double-float dt))
  (with-accessors ((stretch-tensor cl-mpm/particle::mp-stretch-tensor-fbar))
      mp
    (cl-mpm/fastmaths::fast-zero stretch-tensor)
    (iterate-over-neighbours
     mesh mp
     (lambda (mesh mp node svp grads fsvp fgrads)
       (declare (ignore mesh mp svp fsvp fgrads))
       (with-accessors ((node-vel cl-mpm/mesh:node-velocity)
                        (node-active cl-mpm/mesh:node-active))
           node
         (declare (magicl:matrix/double-float node-vel)
                  (boolean node-active))
         (when node-active
           (cl-mpm/shape-function::@-combi-assemble-dstretch-3d grads node-vel stretch-tensor)))))
    (cl-mpm/fastmaths::fast-scale! stretch-tensor dt)))

(defun update-stress-p (mesh mps dt)
  "Update all stresses, with optional f-bar"
  (declare ((array cl-mpm/particle:particle) mps) (cl-mpm/mesh::mesh mesh))
  (iterate-over-mps
   mps
   (lambda (mp)
     (calculate-strain-rate-p mesh mp dt)))
  (values))


(defun calculate-df-p (mesh mp fbar)
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
           (dJ 1d0))
      (cl-mpm/fastmaths::fast-.+-matrix df stretch-tensor df)
      (setf dJ (cl-mpm/fastmaths:det-3x3 df))
      ;;Explicit fbar
      (when fbar
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
              1d0
              ;; (the double-float (/ 1d0 nd))
              ))
            (when (= nd 2)
              (setf (magicl:tref df 2 2) 1d0)))))
      (values df dJ))))
