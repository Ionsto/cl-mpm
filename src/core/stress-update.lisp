(in-package :cl-mpm)
;;Various types of stress update
(declaim (optimize (debug 0) (safety 0) (speed 3)))


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
           (dJ 1d0))
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
                ;; (setf (cl-mpm/particle::mp-debug-j mp) (cl-mpm/fastmaths:det-3x3 df)
                ;;       (cl-mpm/particle::mp-debug-j-gather mp) (cl-mpm/fastmaths:det-3x3 df-fbar))
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
        ))))

(declaim (inline update-strain-kirchoff))
(declaim (ftype (function (cl-mpm/mesh::mesh
                           cl-mpm/particle:particle
                           double-float
                           boolean) (values))
               update-strain-kirchoff))


(defun update-strain-kirchoff (mesh mp dt fbar)
  "Finite strain kirchhoff strain update algorithm"
  (with-accessors ((volume cl-mpm/particle:mp-volume)
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
          (setf volume (* volume (the double-float dj)))
          (when (<= volume 0d0)
            (error "Negative volume"))))))
  (values))

(defun update-strain-kirchoff-dynamic-relaxation (mesh mp dt fbar)
  "Finite strain kirchhoff strain update algorithm"
  (with-accessors ((volume cl-mpm/particle:mp-volume)
                   (volume-0 cl-mpm/particle::mp-volume-0)
                   (strain cl-mpm/particle:mp-strain)
                   (strain-n cl-mpm/particle:mp-strain-n)
                   (def    cl-mpm/particle:mp-deformation-gradient)
                   (def-0    cl-mpm/particle::mp-deformation-gradient-0)
                   (df-inc    cl-mpm/particle::mp-deformation-gradient-increment)
                   ) mp
    (declare (type double-float volume))
    (progn
      (multiple-value-bind (df dj) (calculate-df mesh mp fbar)
        (progn
          ;(setf df-inc df)
          (cl-mpm/utils:matrix-copy-into df df-inc)
          ;; (setf df-inc (cl-mpm/fastmaths::fast-@-matrix-matrix df df-inc))
          (setf def (cl-mpm/fastmaths::fast-@-matrix-matrix df-inc def-0 def))
          (cl-mpm/utils:voigt-copy-into strain-n strain)
          (cl-mpm/ext:kirchoff-update strain df-inc)
          ;(setf volume (* volume (the double-float dj)))
          ;; (setf volume (* volume (the double-float (cl-mpm/fastmaths:det-3x3 df))))
          ;(setf volume (* volume-0 (the double-float (cl-mpm/fastmaths:det-3x3 def))))
          (setf volume (* volume-0 (the double-float (cl-mpm/fastmaths:det-3x3 def-0))))
          (when (<= volume 0d0)
            (error "Negative volume"))))))
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
        ;; Turn kirchoff stress to cauchy
        (cl-mpm/utils::voigt-copy-into stress-kirchoff stress)
        (cl-mpm/fastmaths::fast-scale! stress (/ 1.0d0 (the double-float (cl-mpm/fastmaths:det-3x3 def))))
        ))))

(defun update-stress-kirchoff-dynamic-relaxation (mesh mp dt fbar)
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
        (update-strain-kirchoff-dynamic-relaxation mesh mp dt fbar)
        (cl-mpm/utils::voigt-copy-into (cl-mpm/particle:constitutive-model mp strain dt) stress-kirchoff)
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
