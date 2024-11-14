(in-package :cl-mpm)
;;All the various ways of iterating over the mesh
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(defun update-domain-deformation-rate (domain df)
  "Update the domain length based on the increment of defomation rate"
  (setf (varef domain 0) (* (the double-float (varef domain 0))
                             (the double-float (mtref df 0 0))))
  (setf (varef domain 1) (* (the double-float (varef domain 1))
                             (the double-float (mtref df 1 1))))
  )
(defun update-domain-stretch-rate (df domain)
  "Update the domain length based on the increment of the stretch rate"
  (let ((F (cl-mpm/utils::matrix-zeros)))
    (magicl:mult df df :target F :transb :t)
    (multiple-value-bind (l v) (cl-mpm/utils::eig F)
      (let* ((stretch
              (magicl:@
               v
               (cl-mpm/utils::matrix-from-list
                (list (the double-float (sqrt (the double-float (nth 0 l)))) 0d0 0d0
                      0d0 (the double-float (sqrt (the double-float (nth 1 l)))) 0d0
                      0d0 0d0 (the double-float (sqrt (the double-float (nth 2 l))))))
               (magicl:transpose v)))
            )
        (declare (type magicl:matrix/double-float stretch))
        (setf (varef domain 0) (* (the double-float (varef domain 0))
                                   (the double-float (mtref stretch 0 0))))
        (setf (varef domain 1) (* (the double-float (varef domain 1))
                                   (the double-float (mtref stretch 1 1))))
        (setf (varef domain 2) (* (the double-float (varef domain 2))
                                   (the double-float (mtref stretch 2 2))))
        ))))


(declaim (ftype (function (magicl:matrix/double-float
                           double-float
                           magicl:matrix/double-float
                           double-float) (values))
                update-domain-stretch-rate-damage))
(defun update-domain-stretch-rate-damage (stretch-rate damage domain
                                          damage-domain-rate)
  "Update the domain length based on the increment of the stretch rate"
  (declare (double-float damage damage-domain-rate))
  (let ((df (cl-mpm/utils::matrix-from-list '(1d0 0d0 0d0
                                              0d0 1d0 0d0
                                              0d0 0d0 1d0)))
        (degredation (- 1d0 (* damage-domain-rate damage)))
        (domain-array (cl-mpm/utils:fast-storage domain))
        )

    (cl-mpm/fastmaths:fast-.+ df (magicl:scale stretch-rate degredation) df)
    ;; (cl-mpm/fastmaths::fast-.+ df (magicl:scale stretch-rate degredation) df)
    (let ((F (cl-mpm/utils::matrix-zeros)))
      (magicl:mult df df :target F :transb :t)
      (multiple-value-bind (l v) (cl-mpm/utils::eig F)
        (destructuring-bind (l1 l2 l3) l
          (declare (double-float l1 l2 l3))
          (let* ((stretch
                   (magicl:@
                    v
                    (cl-mpm/utils::matrix-from-list
                     (list (the double-float (sqrt l1)) 0d0 0d0
                           0d0 (the double-float (sqrt l2)) 0d0
                           0d0 0d0 (the double-float (sqrt l3))))
                    (magicl:transpose v)))
                 )
            (declare (type magicl:matrix/double-float stretch))
            (setf (aref domain-array 0) (* (the double-float (varef domain 0))
                                       (the double-float (mtref stretch 0 0))))
            (setf (aref domain-array 1) (* (the double-float (varef domain 1))
                                       (the double-float (mtref stretch 1 1))))
            (setf (aref domain-array 2) (* (the double-float (varef domain 2))
                                       (the double-float (mtref stretch 2 2))))
            ))))))

(defun update-domain-stretch (def domain domain-0)
  "Update the domain length based on the total stretch rate"
  (let ((F (cl-mpm/utils::matrix-zeros)))
    (magicl:mult def def :target F :transb :t)
    (multiple-value-bind (l v) (cl-mpm/utils::eig F)
      (let* ((stretch
              (magicl:@
               v
               (cl-mpm/utils::matrix-from-list
                (list (the double-float (sqrt (the double-float (nth 0 l)))) 0d0 0d0
                      0d0 (the double-float (sqrt (the double-float (nth 1 l)))) 0d0
                      0d0 0d0 (the double-float (sqrt (the double-float (nth 2 l))))))
               (magicl:transpose v)))
            )
        (declare (type magicl:matrix/double-float stretch))
        (setf (varef domain 0) (* (the double-float (varef domain-0 0))
                                   (the double-float (mtref stretch 0 0))))
        (setf (varef domain 1) (* (the double-float (varef domain-0 1))
                                   (the double-float (mtref stretch 1 1))))
        (setf (varef domain 2) (* (the double-float (varef domain-0 2))
                                   (the double-float (mtref stretch 2 2))))
        ))))
(defun update-domain-midpoint (mesh mp dt)
  "Use a corner tracking scheme to update domain lengths"
  (with-accessors ((position cl-mpm/particle::mp-position)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0)
                   )
      mp
    (let ((nd (the fixnum (cl-mpm/mesh:mesh-nd mesh))))
      (with-accessors ((mesh-size cl-mpm/mesh::mesh-mesh-size))
          mesh
        (let ((diff (make-array 3 :initial-element 0d0))
              (domain-storage (magicl::matrix/double-float-storage domain)))
          (loop for d from 0 below nd
                do
                   (let ((disp (cl-mpm/utils:vector-zeros)))
                     (loop for direction in (list 1d0 -1d0)
                           do
                              (let ((corner (cl-mpm/utils::vector-copy position)))
                                (incf (varef corner d)
                                      (* direction 0.5d0 (aref domain-storage d)))
                                (iterate-over-neighbours-point-linear
                                 mesh corner
                                 (lambda (mesh node svp grads)
                                   (declare (double-float dt svp))
                                   (with-accessors ((vel cl-mpm/mesh:node-velocity))
                                       node
                                     (cl-mpm/fastmaths:fast-fmacc corner vel (* dt svp))
                                     )))
                                (cl-mpm/fastmaths:fast-fmacc disp corner direction)
                                ))
                     (setf (aref diff d) (varef disp d))))
          (setf
           (aref domain-storage 0) (aref diff 0)
           (aref domain-storage 1) (aref diff 1)
           (aref domain-storage 2) (aref diff 2))
          ;; (incf (aref domain-storage 0) (aref diff 0))
          ;; (incf (aref domain-storage 1) (aref diff 1))
          ;; (incf (aref domain-storage 2) (aref diff 2))
          (if (= 2 nd)
              (let* ((jf  (magicl:det def))
                     (jl  (* (varef domain 0) (varef domain 1)))
                     (jl0 (* (varef domain-0 0) (varef domain-0 1)))
                     (scaling (expt (/ (* jf jl0) jl) (/ 1d0 2d0))))
                (setf (varef domain 0) (* (varef domain 0) scaling)
                      (varef domain 1) (* (varef domain 1) scaling)))
              (let* ((jf  (magicl:det def))
                     (jl  (* (varef domain 0) (varef domain 1) (varef domain 2)))
                     (jl0 (* (varef domain-0 0) (varef domain-0 1) (varef domain-0 2)))
                     (scaling (expt (/ (* jf jl0) jl) (/ 1d0 3d0))))
                (setf (varef domain 0) (* (varef domain 0) scaling)
                      (varef domain 1) (* (varef domain 1) scaling)
                      (varef domain 2) (* (varef domain 2) scaling)
                      ))))))))

(defun update-domain-corner-2d (mesh mp dt)
  "Use a corner tracking scheme to update domain lengths"
  (with-accessors ((position cl-mpm/particle::mp-position)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0)
                   )
      mp
    (with-accessors ((mesh-size cl-mpm/mesh::mesh-mesh-size))
        mesh
        (let ((diff (make-array 2 :initial-element 0d0))
              (domain-storage (magicl::matrix/double-float-storage domain)))
          (array-operations/utilities:nested-loop (x y) '(2 2)
            (let ((corner (cl-mpm/utils:vector-zeros))
                  (disp (cl-mpm/utils:vector-zeros)))
              (cl-mpm/fastmaths::fast-.+-vector
               position
               (cl-mpm/fastmaths:fast-scale!
                (magicl:.*
                 (vector-from-list (mapcar (lambda (x) (- (* 2d0 (coerce x 'double-float)) 1d0)) (list x y 0)))
                 domain
                 ) 0.5d0) corner)
              (loop for i from 0 to 1
                    do (setf (the double-float (varef corner i))
                             (max 0d0 (min
                                       (the double-float (coerce (nth i mesh-size) 'double-float))
                                       (the double-float (varef corner i))))))
              (iterate-over-neighbours-point-linear
               mesh corner
               (lambda (mesh node svp grads)
                 (declare (double-float dt svp))
                 (with-accessors ((vel cl-mpm/mesh:node-velocity))
                     node
                   (cl-mpm/fastmaths:fast-fmacc disp vel (* dt svp)))))

              (incf (the double-float (aref diff 0)) (* 1.0d0 (the double-float (varef disp 0)) (- (* 2d0 (coerce x 'double-float)) 1d0)))
              (incf (the double-float (aref diff 1)) (* 1.0d0 (the double-float (varef disp 1)) (- (* 2d0 (coerce y 'double-float)) 1d0)))
              ))
          (incf (the double-float (aref domain-storage 0)) (* 0.5d0 (the double-float (aref diff 0))))
          (incf (the double-float (aref domain-storage 1)) (* 0.5d0 (the double-float (aref diff 1))))
          (let* ((jf  (the double-float (magicl:det def)))
                 (jl  (* (the double-float (varef domain 0))
                         (the double-float (varef domain 1))))
                 (jl0 (* (the double-float (varef domain-0 0))
                         (the double-float (varef domain-0 1))))
                 (scaling (the double-float
                               (expt (the double-float (/ (the double-float (* (the double-float jf) (the double-float jl0))) (the double-float jl)))
                                     (the double-float (/ 1d0 2d0))))))
            (setf (varef domain 0) (* (the double-float (varef domain 0)) scaling)
                  (varef domain 1) (* (the double-float (varef domain 1)) scaling)
                  ))
          ))))

(defun update-domain-corner-3d (mesh mp dt)
  "Use a corner tracking scheme to update domain lengths"
  (with-accessors ((position cl-mpm/particle::mp-position)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0)
                   )
      mp
    (with-accessors ((mesh-size cl-mpm/mesh::mesh-mesh-size))
        mesh
        (let ((diff (make-array 3 :initial-element 0d0))
              (domain-storage (magicl::matrix/double-float-storage domain)))
          (array-operations/utilities:nested-loop (x y z) '(2 2 2)
            (let ((corner (cl-mpm/utils:vector-zeros))
                  (disp (cl-mpm/utils:vector-zeros)))
              (cl-mpm/fastmaths::fast-.+-vector position
                                    (magicl:scale!
                                     (magicl:.*
                                      (vector-from-list (mapcar (lambda (x) (- (* 2d0 (coerce x 'double-float)) 1d0)) (list x y z)))
                                      domain
                                      ) 0.5d0) corner)
              (loop for i from 0 to 2
                    do (setf (the double-float (varef corner i))
                             (max 0d0 (min
                                       (the double-float (coerce (nth i mesh-size) 'double-float))
                                       (the double-float (varef corner i))))))
              (iterate-over-neighbours-point-linear
               mesh corner
               (lambda (mesh node svp grads)
                 (declare (double-float dt svp))
                 (with-accessors ((vel cl-mpm/mesh:node-velocity))
                     node
                   (cl-mpm/fastmaths:fast-fmacc disp vel (* dt svp)))))

              (incf (the double-float (aref diff 0)) (* (the double-float (varef disp 0)) (- (* 2d0 (coerce x 'double-float)) 1d0)))
              (incf (the double-float (aref diff 1)) (* (the double-float (varef disp 1)) (- (* 2d0 (coerce y 'double-float)) 1d0)))
              (incf (the double-float (aref diff 2)) (* (the double-float (varef disp 2)) (- (* 2d0 (coerce z 'double-float)) 1d0)))
              ))
          (let ((nd (the fixnum (cl-mpm/mesh:mesh-nd mesh))))
            (incf (the double-float (aref domain-storage 0)) (* 0.5d0 (the double-float (aref diff 0))))
            (incf (the double-float (aref domain-storage 1)) (* 0.5d0 (the double-float (aref diff 1))))
            (incf (the double-float (aref domain-storage 2)) (* 0.5d0 (the double-float (aref diff 2))))
            (let* ((jf  (magicl:det def))
                   (jl  (* (varef domain 0) (varef domain 1) (varef domain 2)))
                   (jl0 (* (varef domain-0 0) (varef domain-0 1) (varef domain-0 2)))
                   (scaling (expt (/ (* jf jl0) jl) (/ 1d0 3d0))))
              (setf (varef domain 0) (* (varef domain 0) scaling)
                    (varef domain 1) (* (varef domain 1) scaling)
                    (varef domain 2) (* (varef domain 2) scaling)
                    ))
            )
          ))))
(defun update-domain-corner (mesh mp dt)
  (if (= (the fixnum (cl-mpm/mesh:mesh-nd mesh)) 2)
      (update-domain-corner-2d mesh mp dt)
      (update-domain-corner-3d mesh mp dt)))

