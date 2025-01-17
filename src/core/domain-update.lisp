(in-package :cl-mpm)
;;All the various ways of iterating over the mesh
(declaim (optimize (debug 0) (safety 0) (speed 3)))


(defun scale-domain-size-2d (mp)
  (with-accessors ((def cl-mpm/particle::mp-deformation-gradient)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0))
      mp
    (declare (magicl::matrix/double-float domain domain-0 def))
    (let* ((jf  (cl-mpm/fastmaths:det-3x3 def))
           (jl  (* (varef domain 0) (varef domain 1)))
           (jl0 (* (varef domain-0 0) (varef domain-0 1)))
           (scaling (the double-float (expt (the double-float (/ (* jf jl0) jl)) (/ 1d0 2d0)))))
      (declare (double-float jf jl0 jl))
      (setf (varef domain 0) (* (varef domain 0) scaling)
            (varef domain 1) (* (varef domain 1) scaling)))))

(defun scale-domain-size-3d (mp)
  (with-accessors ((def cl-mpm/particle::mp-deformation-gradient)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0))
      mp
    (declare (magicl::matrix/double-float domain domain-0 def))
    (let* ((jf  (cl-mpm/fastmaths:det-3x3 def))
           (jl  (* (varef domain 0) (varef domain 1) (varef domain 2)))
           (jl0 (* (varef domain-0 0) (varef domain-0 1) (varef domain-0 2)))
           (scaling (expt (/ (* jf jl0) jl) (/ 1d0 3d0))))
      (declare (double-float scaling jf jl0 jl))
      (setf (varef domain 0) (* (varef domain 0) scaling)
            (varef domain 1) (* (varef domain 1) scaling)
            (varef domain 2) (* (varef domain 2) scaling)
            )))
  )
(defun scale-domain-size (mesh mp)
  (if (= 2 (cl-mpm/mesh:mesh-nd mesh))
      (scale-domain-size-2d mp)
      (scale-domain-size-3d mp)))

(defun update-domain-deformation-rate (domain df)
  "Update the domain length based on the increment of defomation rate"
  (setf (varef domain 0) (* (the double-float (varef domain 0))
                             (the double-float (mtref df 0 0))))
  (setf (varef domain 1) (* (the double-float (varef domain 1))
                             (the double-float (mtref df 1 1))))
  )
(defun update-domain-stretch-rate (mesh mp dt)
  "Update the domain length based on the increment of the stretch rate"
  (with-accessors ((df cl-mpm/particle::mp-deformation-gradient-increment)
                   (domain cl-mpm/particle::mp-domain-size))
      mp
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
          )))))


(declaim (ftype (function (magicl:matrix/double-float
                           double-float
                           magicl:matrix/double-float
                           double-float) (values))
                update-domain-stretch-rate-damage))
(defun update-domain-stretch-rate-damage (stretch-rate damage domain damage-domain-rate)
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
            )))))
  )
(defun update-domain-F (def domain domain-0)
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
                                   (the double-float (mtref F 0 0))))
        (setf (varef domain 1) (* (the double-float (varef domain-0 1))
                                   (the double-float (mtref F 1 1))))
        (setf (varef domain 2) (* (the double-float (varef domain-0 2))
                                   (the double-float (mtref F 2 2))))
        ))))

(defun update-domain-stretch (mesh mp dt)
  "Update the domain length based on the increment of the stretch rate"
  (with-accessors ((def cl-mpm/particle::mp-deformation-gradient)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0)
                   )
      mp
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
          )))))

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
              (let* ((jf  (cl-mpm/fastmaths:det-3x3 def))
                     (jl  (* (varef domain 0) (varef domain 1)))
                     (jl0 (* (varef domain-0 0) (varef domain-0 1)))
                     (scaling (expt (/ (* jf jl0) jl) (/ 1d0 2d0))))
                (setf (varef domain 0) (* (varef domain 0) scaling)
                      (varef domain 1) (* (varef domain 1) scaling)))
              (let* ((jf  (cl-mpm/fastmaths:det-3x3 def))
                     (jl  (* (varef domain 0) (varef domain 1) (varef domain 2)))
                     (jl0 (* (varef domain-0 0) (varef domain-0 1) (varef domain-0 2)))
                     (scaling (expt (/ (* jf jl0) jl) (/ 1d0 3d0))))
                (setf (varef domain 0) (* (varef domain 0) scaling)
                      (varef domain 1) (* (varef domain 1) scaling)
                      (varef domain 2) (* (varef domain 2) scaling)
                      ))))))))

(defun clamp-point-domain (mesh point)
  (let ((nd (cl-mpm/mesh:mesh-nd mesh))
        (ms-list (cl-mpm/mesh:mesh-mesh-size mesh)))
    (declare (fixnum nd)
             (list ms-list))
    (loop for i from 0 below nd
          for ms in ms-list
          do (setf (the double-float (varef point i))
                   (max 0d0 (min
                             (the double-float ms)
                             (the double-float (varef point i))))))))

(defun generate-rotation-matrix (omega)
  ;;Finite rotation case
  (let ((empty t))
    (loop for v across (cl-mpm/utils:fast-storage omega)
          while empty
          do (setf empty (= v 0d0)))
    (if empty
        (cl-mpm/utils:matrix-eye 1d0)
        (magicl:expm omega)))
  ;;Infinitesmal rotation case
  ;; (cl-mpm/fastmaths:fast-.+ (cl-mpm/utils:matrix-eye 1d0) omega)
  )

(declaim (notinline co-domain-corner-2d))
(defun co-domain-corner-2d (mesh mp dt)
  "Use a corner tracking scheme to update domain lengths"
  (let ((inc (cl-mpm/utils:vector-zeros)))
    (with-accessors ((position cl-mpm/particle::mp-position)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     (domain cl-mpm/particle::mp-domain-size)
                     (domain-0 cl-mpm/particle::mp-domain-size-0)
                     )
        mp
      (with-accessors ((mesh-size cl-mpm/mesh::mesh-mesh-size))
          mesh
        (let ((diff (make-array 2 :initial-element 0d0 :element-type 'double-float))
              (domain-storage (magicl::matrix/double-float-storage domain)))
          (iterate-over-corners-normal-2d
           mesh mp
           (lambda (corner normal)
             (let ((disp (cl-mpm/utils:vector-zeros)))
               (iterate-over-neighbours-point-linear-simd
                mesh corner
                (lambda (mesh node svp grads)
                  (declare (double-float dt svp))
                  (with-accessors ((vel cl-mpm/mesh:node-velocity))
                      node
                    (cl-mpm/fastmaths:fast-fmacc disp vel (* dt svp)))))
               (incf (the double-float (aref diff 0)) (* 1d0
                                                         (the double-float (varef disp 0))
                                                         (varef normal 0)))
               (incf (the double-float (aref diff 1)) (* 1d0
                                                         (the double-float (varef disp 1))
                                                         (varef normal 1))))))

          (incf (the double-float (cl-mpm/utils:varef inc 0)) (* 0.5d0 (the double-float (aref diff 0))))
          (incf (the double-float (cl-mpm/utils:varef inc 1)) (* 0.5d0 (the double-float (aref diff 1))))
          )))

    (with-accessors ((domain cl-mpm/particle::mp-domain-size)
                     (true-domain cl-mpm/particle::mp-true-domain)
                     (D cl-mpm/particle::mp-stretch-tensor))
        mp
      (let* ((omega (magicl:scale
                     (magicl:.-
                      D
                      (magicl:transpose D)) 0.5d0))
             (dom
               true-domain)
             (dom-inc (cl-mpm/utils:matrix-from-list
                   (list
                    (cl-mpm/utils:varef inc 0) 0d0 0d0
                    0d0 (cl-mpm/utils:varef inc 1) 0d0
                    0d0 0d0 0d0)))
             ;; (co-inc (cl-mpm/fastmaths:fast-.+
             ;;          dom-inc
             ;;          (cl-mpm/fastmaths::fast-.-
             ;;           (magicl:@ omega dom)
             ;;           (magicl:@ dom omega)
             ;;           )))
             ;(R (magicl:.+ (magicl:eye 3) omega))
             (corner-stretch (cl-mpm/utils:matrix-from-list (list
                                                             (+ 1d0 (/ (varef inc 0) (varef domain 0))) 0d0 0d0
                                                             0d0 (+ 1d0 (/ (varef inc 1) (varef domain 1))) 0d0
                                                             0d0 0d0 1d0
                                                             )))
             (R (generate-rotation-matrix omega))
             )
        ;; (multiple-value-bind (u s v) (magicl:svd (cl-mpm/particle::mp-deformation-gradient-increment mp))
        ;;   (let* ((Rdef (magicl:@ u (magicl:transpose v))))
        ;;     (setf R Rdef)
        ;;     ))
        ;; (cl-mpm/fastmaths:fast-.+
        ;;  true-domain
        ;;  dom-inc
        ;;  true-domain)
        (setf true-domain (magicl:@ R true-domain corner-stretch (magicl:transpose R)))

        ;; (cl-mpm/fastmaths::fast-.+
        ;;  true-domain
        ;;  co-inc
        ;;  true-domain)
        (setf
         (varef domain 0)
         ;; (magicl:tref true-domain 0 0)
         (;cl-mpm/fastmaths:mag
          varef
          (cl-mpm/fastmaths::fast-@-matrix-vector
           true-domain
           (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0)))
          0
          )
         (varef domain 1)
         ;; (magicl:tref true-domain 1 1)
         (;cl-mpm/fastmaths:mag
          varef
          (cl-mpm/fastmaths::fast-@-matrix-vector
           true-domain
           (cl-mpm/utils:vector-from-list (list 0d0 1d0 0d0)))
          0)
         ;; (varef domain 2)
         ;; (cl-mpm/fastmaths:mag
         ;;  (cl-mpm/fastmaths::fast-@-matrix-vector
         ;;   true-domain
         ;;   (cl-mpm/utils:vector-from-list (list 0d0 0d0 1d0)))
         ;;  )
         )
        ;; (incf (cl-mpm/utils:varef domain 0)
        ;;       (magicl:tref co-inc 0 0))
        ;; (incf (cl-mpm/utils:varef domain 1)
        ;;       (magicl:tref co-inc 1 1))
        ))
    ))
(defun update-domain-polar-2d (mesh mp dt)
  "Use a corner tracking scheme to update domain lengths"
  (with-accessors ((dF cl-mpm/particle::mp-deformation-gradient-increment)
                   (domain cl-mpm/particle::mp-domain-size)
                   (true-domain cl-mpm/particle::mp-true-domain)
                   (D cl-mpm/particle::mp-stretch-tensor))
      mp
    (multiple-value-bind (u s vt) (magicl:svd dF)
      (let* ((R (magicl:@ u vt))
             (U (magicl:@ (magicl:transpose vt) s vt))
             )
        ;; (setf true-domain (magicl:@ R (magicl:@ true-domain U) (magicl:transpose R)))
        (setf true-domain (magicl:@ R (magicl:@ true-domain U) (magicl:transpose R)))
        ;; (setf true-domain (magicl:@ true-domain U))
        ;; (setf true-domain (magicl:@ R true-domain (magicl:transpose R)))
        ;; (setf true-domain (magicl:@ U true-domain))
        ))
    ;; (setf true-domain (magicl:@ R true-domain U (magicl:transpose R)))
    ;; (setf true-domain (magicl:@ R (magicl:@ true-domain U) (magicl:transpose R)))
    ;; (setf true-domain (magicl:@ true-domain dF))

    (setf
     (varef domain 0)
     (cl-mpm/fastmaths:mag
      (cl-mpm/fastmaths::fast-@-matrix-vector
       true-domain
       (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0))))
     (varef domain 1)
     (cl-mpm/fastmaths:mag
      (cl-mpm/fastmaths::fast-@-matrix-vector
       true-domain
       (cl-mpm/utils:vector-from-list (list 0d0 1d0 0d0)))))))

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
        (let ((diff (make-array 2 :initial-element 0d0 :element-type 'double-float))
              (domain-storage (magicl::matrix/double-float-storage domain)))
          (iterate-over-corners-normal-2d
           mesh mp
           (lambda (corner normal)
             (let ((disp (cl-mpm/utils:vector-zeros)))
               (iterate-over-neighbours-point-linear-simd
                mesh corner
                (lambda (mesh node svp grads)
                  (declare (double-float dt svp))
                  (with-accessors ((vel cl-mpm/mesh:node-velocity))
                      node
                    (cl-mpm/fastmaths:fast-fmacc disp vel (* dt svp)))))
               (incf (the double-float (aref diff 0)) (* 1d0
                                                         (the double-float (varef disp 0))
                                                         (varef normal 0)))
               (incf (the double-float (aref diff 1)) (* 1d0
                                                         (the double-float (varef disp 1))
                                                         (varef normal 1))))))

          (incf (the double-float (aref domain-storage 0)) (* 0.5d0 (the double-float (aref diff 0))))
          (incf (the double-float (aref domain-storage 1)) (* 0.5d0 (the double-float (aref diff 1))))
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
          (iterate-over-corners-normal-3d
           mesh mp
           (lambda (corner normal)
             (let ((disp (cl-mpm/utils:vector-zeros)))
               (iterate-over-neighbours-point-linear-3d
                mesh corner
                (lambda (mesh node svp grads)
                  (declare (double-float dt svp))
                  (with-accessors ((vel cl-mpm/mesh:node-velocity))
                      node
                    (cl-mpm/fastmaths:fast-fmacc disp vel (* dt svp)))))
               (incf (the double-float (aref diff 0)) (* (the double-float (varef disp 0))
                                                         (varef normal 0)))
               (incf (the double-float (aref diff 1)) (* (the double-float (varef disp 1))
                                                         (varef normal 1)))
               (incf (the double-float (aref diff 2)) (* (the double-float (varef disp 2))
                                                         (varef normal 2))))))))))
(defun update-domain-corner (mesh mp dt)
  (if (= (the fixnum (cl-mpm/mesh:mesh-nd mesh)) 2)
      (update-domain-corner-2d mesh mp dt)
      (update-domain-corner-3d mesh mp dt)))


(defun update-domain-det (mesh mp)
  "Update the domain length based on the increment of the stretch rate"
  (with-accessors ((domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0)
                   (def cl-mpm/particle::mp-deformation-gradient))
      mp
    (let ((det (cl-mpm/fastmaths:det-3x3 def)))
      (declare (double-float det))
      (if (= (the fixnum (cl-mpm/mesh:mesh-nd mesh)) 2)
          (let ((scale (expt det 1/2)))
            (declare (double-float scale))
            (setf (varef domain 0) (* (the double-float (varef domain-0 0)) scale)
                  (varef domain 1) (* (the double-float (varef domain-0 1)) scale)))
          (let ((scale (expt det 1/3)))
            (declare (double-float scale))
            (setf (varef domain 0) (* (the double-float (varef domain-0 0)) scale)
                  (varef domain 1) (* (the double-float (varef domain-0 1)) scale)
                  (varef domain 2) (* (the double-float (varef domain-0 2)) scale)))))))


(defun update-domain-max-corner-2d (mesh mp dt)
  "Use a corner tracking scheme to update domain lengths"
  (with-accessors ((position cl-mpm/particle::mp-position)
                   (def cl-mpm/particle::mp-deformation-gradient)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0)
                   )
      mp
    (with-accessors ((mesh-size cl-mpm/mesh::mesh-mesh-size))
        mesh
        (let ((max-diff (make-array 2 :initial-element 0d0 :element-type 'double-float))
              (min-diff (make-array 2 :initial-element 0d0 :element-type 'double-float))
              (domain-storage (magicl::matrix/double-float-storage domain)))
          (iterate-over-midpoints-normal-2d
           mesh mp
           (lambda (corner normal)
             (let ((disp (cl-mpm/utils:vector-zeros)))
               (iterate-over-neighbours-point-linear-simd
                mesh corner
                (lambda (mesh node svp grads)
                  (declare (double-float dt svp))
                  (with-accessors ((vel cl-mpm/mesh:node-velocity))
                      node
                    (cl-mpm/fastmaths:fast-fmacc disp vel (* dt svp)))))
               (loop for i from 0 to 1
                     do
                        (progn
                          (setf (the double-float (aref max-diff i))
                                (max (aref max-diff i)
                                     (+ (varef corner i) (varef disp i))))
                          (setf (the double-float (aref min-diff i))
                                (min (aref min-diff i)
                                     (+ (varef corner i) (varef disp i)))))))))

          (setf (the double-float (aref domain-storage 0))
                (* 0.5d0 (- (aref max-diff 0) (aref min-diff 0))))
          (setf (the double-float (aref domain-storage 1))
                (* 0.5d0 (- (aref max-diff 1) (aref min-diff 1))))
          ))))

(defun update-domain-deformation (mesh mp dt)
  "Update the domain length based on the increment of the stretch rate"
  (with-accessors ((def cl-mpm/particle::mp-deformation-gradient)
                   (domain cl-mpm/particle::mp-domain-size)
                   (domain-0 cl-mpm/particle::mp-domain-size-0)
                   )
      mp
    
    (let* ((d-0 (cl-mpm/utils::matrix-diag (list (expt (varef domain-0 0) 2)
                                                 (expt (varef domain-0 1) 2)
                                                 (expt (varef domain-0 2) 2))))
           (l (magicl:@ def d-0 (magicl:transpose def))))
      (multiple-value-bind (l v) (cl-mpm/utils::eig l)
        (let* ((stretch
                 (magicl:@
                  v
                  (cl-mpm/utils::matrix-from-list
                   (list (the double-float (sqrt (the double-float (nth 0 l)))) 0d0 0d0
                         0d0 (the double-float (sqrt (the double-float (nth 1 l)))) 0d0
                         0d0 0d0 (the double-float (sqrt (the double-float (nth 2 l))))))
                  (magicl:transpose v)))
               )
          (setf (varef domain 0) (abs (magicl:tref stretch 0 0)))
          (setf (varef domain 1) (abs (magicl:tref stretch 1 1)))))
      ;; (setf
      ;;  (varef domain 0)
      ;;  ;; (magicl:tref true-domain 0 0)
      ;;  (cl-mpm/fastmaths:mag
      ;;   (cl-mpm/fastmaths::fast-@-matrix-vector
      ;;    l
      ;;    (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0))))
      ;;  (varef domain 1)
      ;;  ;; (magicl:tref true-domain 1 1)
      ;;  (cl-mpm/fastmaths:mag
      ;;   (cl-mpm/fastmaths::fast-@-matrix-vector
      ;;    l
      ;;    (cl-mpm/utils:vector-from-list (list 0d0 1d0 0d0))))
      ;;  ;; (varef domain 2)
      ;;  ;; (cl-mpm/fastmaths:mag
      ;;  ;;  (cl-mpm/fastmaths::fast-@-matrix-vector
      ;;  ;;   true-domain
      ;;  ;;   (cl-mpm/utils:vector-from-list (list 0d0 0d0 1d0)))
      ;;  ;;  )
      ;;  )
      ;; (setf (varef domain 2) (* (the double-float (varef domain-0 2))
      ;;                           ))
      )))

;; (defun update-domain-deformation (mesh mp dt)
;;   "Update the domain length based on the increment of the stretch rate"
;;   (with-accessors ((def cl-mpm/particle::mp-deformation-gradient)
;;                    (domain cl-mpm/particle::mp-domain-size)
;;                    (domain-0 cl-mpm/particle::mp-domain-size-0)
;;                    )
;;       mp
    
;;     (let* ((d-0 (cl-mpm/utils::matrix-diag (list (expt (varef domain-0 0) 1)
;;                                                  (expt (varef domain-0 1) 1)
;;                                                  (expt (varef domain-0 2) 1))))
;;            ;(l (magicl:@ def d-0 (magicl:transpose def)))
;;            (l (magicl:@ def d-0 (magicl:transpose def)))
;;            )
;;       ;; (setf (varef domain 0) (abs (magicl:tref stretch 0 0)))
;;       ;; (setf (varef domain 1) (abs (magicl:tref stretch 1 1)))
;;       (setf
;;        (varef domain 0)
;;        ;; (magicl:tref true-domain 0 0)
;;        (cl-mpm/fastmaths:mag
;;         (cl-mpm/fastmaths::fast-@-matrix-vector
;;          l
;;          (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0))))
;;        (varef domain 1)
;;        ;; (magicl:tref true-domain 1 1)
;;        (cl-mpm/fastmaths:mag
;;         (cl-mpm/fastmaths::fast-@-matrix-vector
;;          l
;;          (cl-mpm/utils:vector-from-list (list 0d0 1d0 0d0))))
;;        ;; (varef domain 2)
;;        ;; (cl-mpm/fastmaths:mag
;;        ;;  (cl-mpm/fastmaths::fast-@-matrix-vector
;;        ;;   true-domain
;;        ;;   (cl-mpm/utils:vector-from-list (list 0d0 0d0 1d0)))
;;        ;;  )
;;        )
;;       ;; (setf (varef domain 2) (* (the double-float (varef domain-0 2))
;;       ;;                           ))
;;       )))
