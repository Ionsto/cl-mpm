(in-package :cl-mpm)
;;MP splitting code
(declaim #.cl-mpm/settings:*optimise-setting*)


(defun copy-particle (original &rest initargs &key &allow-other-keys)
  "Function for copying a particle, allowing for re-initialising members"
  (let* ((class (class-of original))
         (copy (allocate-instance class)))
    (dolist (slot (mapcar #'sb-mop:slot-definition-name (sb-mop:class-slots class)))
      (when (slot-boundp original slot)
        (cond
          ((typep (slot-value original slot) 'magicl::abstract-tensor)
           (setf (slot-value copy slot)
                 ;; (magicl:deep-copy-tensor (slot-value original slot))
                 ;; (magicl:scale (slot-value original slot) 1d0)
                 (cl-mpm/utils::deep-copy (slot-value original slot))
                 ))
          (t (setf (slot-value copy slot)
                  (slot-value original slot))))))
    (apply #'reinitialize-instance copy initargs)))

(defun make-scaling-matrix (vec scale)
  (let* ((vec (cl-mpm/fastmaths:norm vec))
         (a (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0)))
         (dotprod (cl-mpm/fastmaths:dot vec a)))
    (if (= dotprod 1d0)
        (progn
          (cl-mpm/utils::matrix-from-diag (list scale 1d0 1d0)))
        (progn
          (let* ((b (cl-mpm/fastmaths:norm (cl-mpm/fastmaths::cross-product vec a)))
                 (angle (acos dotprod))
                 (bx (varef b 0))
                 (by (varef b 1))
                 (bz (varef b 2))
                 (K (cl-mpm/utils:matrix-from-list (list 0d0 (- bz) by
                                                         bz 0d0 (- bx)
                                                         (- by) bx 0d0)))
                 (R (cl-mpm/fastmaths:fast-.+ (cl-mpm/utils:matrix-eye 1d0)
                                              (cl-mpm/fastmaths:fast-.+
                                               (cl-mpm/fastmaths:fast-scale K (sin angle))
                                               (magicl:@ K (cl-mpm/fastmaths:fast-scale K (- 1d0 (cos angle))))))))
            (magicl:@ R (cl-mpm/utils::matrix-from-diag (list scale 1d0 1d0)) (magicl:transpose R)))))))


(defun split-vector (mp split-vec)
  "Helper macro for single splitting along cartesian directions "
  (with-accessors ((lens cl-mpm/particle::mp-domain-size)
                   (lens-0 cl-mpm/particle::mp-domain-size-0)
                   (mass cl-mpm/particle::mp-mass)
                   (pos cl-mpm/particle:mp-position)
                   (volume cl-mpm/particle:mp-volume)
                   (volume-0 cl-mpm/particle::mp-volume-0)
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
      (let ((domain-scaler (make-scaling-matrix split-vec 0.5d0)))
        (setf new-domain (magicl:@ domain-scaler true-domain))
        (setf pos-offset (magicl:@
                          true-domain
                          (cl-mpm/fastmaths::fast-scale-vector vec-scaler 0.5d0)))
        (setf new-size-0 (magicl:@ domain-scaler lens-0))
        (setf new-size (magicl:@ domain-scaler lens))
        (list
         (copy-particle mp
                        :mass (/ mass 2)
                        :volume (/ volume 2)
                        :volume-0 (/ volume-0 2)
                        :size (cl-mpm/utils::vector-copy new-size)
                        :size-0 (cl-mpm/utils::vector-copy new-size-0)
                        :position (cl-mpm/fastmaths::fast-.+-vector pos pos-offset)
                        :nc (make-array 8 :fill-pointer 0 :element-type 'cl-mpm/particle::node-cache
                                          :initial-element (cl-mpm/particle::make-empty-node-cache))
                        :split-depth new-split-depth
                        :true-domain (cl-mpm/utils:matrix-copy new-domain)
                        )
         (copy-particle mp
                        :mass (/ mass 2)
                        :volume (/ volume 2)
                        :volume-0 (/ volume-0 2)
                        :size (cl-mpm/utils::vector-copy new-size)
                        :size-0 (cl-mpm/utils::vector-copy new-size-0)
                        :position (magicl:.- pos pos-offset)
                        :nc (make-array 8 :fill-pointer 0 :element-type 'cl-mpm/particle::node-cache
                                          :initial-element (cl-mpm/particle::make-empty-node-cache))
                        :split-depth new-split-depth
                        :true-domain (cl-mpm/utils:matrix-copy new-domain)
                        ))))))

;; (defmacro split-linear (dir direction dimension)
;;   "Helper macro for single splitting along cartesian directions "
;;   `((eq ,dir ,direction)
;;     (let ((new-size (vector-from-list (list 1d0 1d0 1d0)))
;;           (new-size-0 (vector-from-list (list 1d0 1d0 1d0)))
;;           (pos-offset (vector-zeros))
;;           (new-split-depth (+ (cl-mpm/particle::mp-split-depth mp) 1))
;;           (true-domain (cl-mpm/particle::mp-true-domain mp))
;;           (new-domain nil))
;;       (setf (tref new-size ,dimension 0) 0.5d0)
;;       (setf (tref new-size-0 ,dimension 0) 0.5d0)
;;       (setf (tref pos-offset ,dimension 0) 0.25d0)
;;       (let ((domain-scaler (magicl:eye 3)))
;;         (setf (tref domain-scaler ,dimension ,dimension) 0.5d0)
;;         (setf new-domain (magicl:@ true-domain domain-scaler)))

;;       (cl-mpm/fastmaths::fast-.* lens new-size new-size)
;;       (cl-mpm/fastmaths::fast-.* lens-0 new-size-0 new-size-0)
;;       (cl-mpm/fastmaths::fast-.* lens pos-offset pos-offset)
;;       (list
;;        (copy-particle mp
;;                       :mass (/ mass 2)
;;                       :volume (/ volume 2)
;;                       :volume-0 (/ volume-0 2)
;;                       :size (cl-mpm/utils::vector-copy new-size)
;;                       :size-0 (cl-mpm/utils::vector-copy new-size-0)
;;                       :position (cl-mpm/fastmaths::fast-.+-vector pos pos-offset)
;;                       :nc (make-array 8 :fill-pointer 0 :element-type 'node-cache)
;;                       :split-depth new-split-depth
;;                       :true-domain (cl-mpm/utils:matrix-copy new-domain)
;;                       )
;;        (copy-particle mp
;;                       :mass (/ mass 2)
;;                       :volume (/ volume 2)
;;                       :volume-0 (/ volume-0 2)
;;                       :size (cl-mpm/utils::vector-copy new-size)
;;                       :size-0 (cl-mpm/utils::vector-copy new-size-0)
;;                       :position (cl-mpm/fastmaths::fast-.--vector pos pos-offset)
;;                       :nc (make-array 8 :fill-pointer 0 :element-type 'node-cache)
;;                       :split-depth new-split-depth
;;                       :true-domain (cl-mpm/utils:matrix-copy new-domain)
;;                       )))))

;; (defmacro split-cases (direction)
;;   "Another helper macro for splitting mps"
;;   `(cond
;;      ,(macroexpand-1 '(split-linear direction :x 0))
;;      ,(macroexpand-1 '(split-linear direction :y 1))
;;      ,(macroexpand-1 '(split-linear direction :z 2))
;;      (t nil)
;;      ))
(defun split-criteria-variable (mp h factor)
  "Some numerical splitting estimates"
  (with-accessors ((def cl-mpm/particle:mp-deformation-gradient)
                   (lens cl-mpm/particle::mp-domain-size)
                   (lens-0 cl-mpm/particle::mp-domain-size-0))
      mp
    (when t
      (let ((h-factor (* factor h)))
        (cond
          ((< h-factor (varef lens 0)) :x)
          ((< h-factor (varef lens 1)) :y)
          ((< h-factor (varef lens 2)) :z)
          (t nil))))))
(defun split-mp (mp h direction)
  "Function to split an mp across a direction
   Directions should be :x,:y,:z "
  (with-accessors ((def cl-mpm/particle:mp-deformation-gradient)
                   (lens cl-mpm/particle::mp-domain-size)
                   (lens-0 cl-mpm/particle::mp-domain-size-0)
                   (pos cl-mpm/particle:mp-position)
                   (mass cl-mpm/particle:mp-mass)
                   (volume cl-mpm/particle:mp-volume)
                   (volume-0 cl-mpm/particle::mp-volume-0)
                   ;; (volume cl-mpm/particle::mp-volume-0)
                   )
      mp
    (let ((dir (cl-mpm/utils:vector-zeros)))
      ;; (setf (varef dir direction) 1d0)
      (ecase direction
        (:x (setf (varef dir 0) 1d0))
        (:y (setf (varef dir 1) 1d0))
        (:z (setf (varef dir 2) 1d0)))
      (split-vector mp dir))))

(defun abs-vec (v)
  (let ((v-d (cl-mpm/utils::vector-zeros)))
    (loop for i from 0 below (length (fast-storage v-d))
          do (setf (varef v-d i) (abs (varef v i))))
    v-d))

(defun split-mps-eigenvalue (sim)
  (declare (optimize (speed 0) (debug 3) (safety 3)))
  (let* ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim)))
         (split-factor (cl-mpm::sim-split-factor sim))
         (crit (* h split-factor)))
    (pprint crit)
    (split-mps-vector
     sim
     (lambda (mp)
       (let ((split-dir nil))
         (let ((td (cl-mpm/particle::mp-true-domain mp)))
           (multiple-value-bind (l v) (cl-mpm/utils:eig
                                       ;; td
                                       (cl-mpm/utils::slice-matrix-nd td nd)
                                                        )
             (let* ((v (cl-mpm/utils::pad-matrix-nd v nd))
                    (abs-l (mapcar #'abs l))
                    (max-l (reduce #'max abs-l))
                    (min-l (reduce #'min (remove 0d0 abs-l))))
               (when (> max-l crit)
                 (let* ((pos (position max-l abs-l))
                        (vec (magicl::vector->column-matrix (magicl:column v pos))))
                   (setf split-dir (cl-mpm/fastmaths:norm
                                    (cl-mpm/utils:vector-from-list (list (varef vec 0) (varef vec 1) (varef vec 2))))))))))
         split-dir)))))


(defun split-mps-cartesian (sim)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)
                   (max-split-depth cl-mpm::sim-max-split-depth)
                   (split-factor cl-mpm::sim-split-factor))
      sim
    (declare (fixnum max-split-depth))
    (let* ((h (cl-mpm/mesh:mesh-resolution mesh))
           (mps-to-split (remove-if-not (lambda (mp)
                                          (and
                                           (< (cl-mpm/particle::mp-split-depth mp) max-split-depth)
                                           (split-criteria-variable mp h split-factor))) mps))
           (split-direction (map 'list (lambda (mp) (split-criteria-variable mp h split-factor)) mps-to-split)))
      (remove-mps-func sim (lambda (mp) (split-criteria-variable mp h split-factor)))
      (loop for mp across mps-to-split
            for direction in split-direction
            do (loop for new-mp in (split-mp mp h direction)
                     do (sim-add-mp sim new-mp)))))
  )

(defun split-mps (sim)
  "Split mps that match the split-criteria"
  (split-mps-eigenvalue sim)
  ;; (split-mps-cartesian sim)
  )

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
                   (mesh cl-mpm:sim-mesh)
                   (max-split-depth cl-mpm::sim-max-split-depth))
      sim
    (let* ((split-directions (lparallel:pmap '(vector t *) criteria mps))
           (mps-to-delete (delete-if-not #'identity (lparallel:pmap '(vector t *)
                                                                   (lambda (mp crit)
                                                                     (if crit mp nil)) mps split-directions)))
           (mps-to-split (remove-if-not (lambda (mp) (< (cl-mpm/particle::mp-split-depth mp) max-split-depth)) mps-to-delete)))
      (remove-mps-func sim (lambda (mp) (position mp mps-to-delete)))
      (loop for mp across mps-to-split
            do (loop for new-mp in (split-vector mp (funcall criteria mp))
                     do (sim-add-mp sim new-mp))))))
