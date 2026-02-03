(in-package :cl-mpm/damage)
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(defmethod cl-mpm/output::save-vtk (filename (sim cl-mpm/damage::mpm-sim-damage))
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh)) sim
    (with-open-file (fs filename :direction :output :if-exists :supersede)
      (format fs "# vtk DataFile Version 2.0~%")
      (format fs "Lisp generated vtk file, SJVS~%")
      (format fs "ASCII~%")
      (format fs "DATASET UNSTRUCTURED_GRID~%")
      (format fs "POINTS ~d double~%" (length mps))
      (loop for mp across mps
            do
               (let ((pos (cl-mpm/particle::mp-position-trial mp)))
                 (format fs "~E ~E ~E ~%"
                         (coerce (cl-mpm/utils:varef pos 0) 'single-float)
                         (coerce (cl-mpm/utils:varef pos 1) 'single-float)
                         (coerce (cl-mpm/utils:varef pos 2) 'single-float))))
      (format fs "~%")
      ;; (cl-mpm/output::with-parameter-list fs mps
      ;;   ("mass" 'cl-mpm/particle:mp-mass)
      ;;   ("density" (lambda (mp) (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp))))
      ;;   )
      (let ((id 1)
            (nd (cl-mpm/mesh:mesh-nd mesh)))
        (declare (special id))
        (format fs "POINT_DATA ~d~%" (length mps))

        (cl-mpm/output::save-parameter "viscosity" (if (slot-exists-p mp 'cl-mpm/particle::viscosity) (cl-mpm/particle::mp-viscosity mp) 0d0))
        (cl-mpm/output::save-parameter "mass" (cl-mpm/particle:mp-mass mp))
        (cl-mpm/output::save-parameter "density" (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)))
        (cl-mpm/output::save-parameter "unique-id" (cl-mpm/particle::mp-unique-index mp))
        (cl-mpm/output::save-parameter "index" (cl-mpm/particle::mp-index mp))
        (cl-mpm/output::save-parameter "mpi-index" (cl-mpm/particle::mp-mpi-index mp))
        (cl-mpm/output::save-parameter "j" (magicl:det (cl-mpm/particle::mp-deformation-gradient mp)))
        (cl-mpm/output::save-parameter "volume" (cl-mpm/particle::mp-volume mp))
        (cl-mpm/output::save-parameter "vel_x" (magicl:tref (cl-mpm/particle:mp-velocity mp) 0 0))
        (cl-mpm/output::save-parameter "vel_y" (magicl:tref (cl-mpm/particle:mp-velocity mp) 1 0))
        ;; (cl-mpm/output::save-parameter "acc_x" (magicl:tref (cl-mpm/particle::mp-acceleration mp) 0 0))
        ;; (cl-mpm/output::save-parameter "acc_y" (magicl:tref (cl-mpm/particle::mp-acceleration mp) 1 0))

        (cl-mpm/output::save-parameter "disp_x" (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0))
        (cl-mpm/output::save-parameter "disp_y" (magicl:tref (cl-mpm/particle::mp-displacement mp) 1 0))
        (when (= 3 nd)
          (cl-mpm/output::save-parameter "disp_z" (magicl:tref (cl-mpm/particle::mp-displacement mp) 2 0)))

        (cl-mpm/output::save-parameter "sig_xx" (magicl:tref (cl-mpm/particle:mp-stress mp) 0 0))
        (cl-mpm/output::save-parameter "sig_yy" (magicl:tref (cl-mpm/particle:mp-stress mp) 1 0))
        (cl-mpm/output::save-parameter "sig_xy" (magicl:tref (cl-mpm/particle:mp-stress mp) 5 0))
        (when (= 3 nd)
          (cl-mpm/output::save-parameter "sig_zz" (magicl:tref (cl-mpm/particle:mp-stress mp) 2 0))
          (cl-mpm/output::save-parameter "sig_yz" (magicl:tref (cl-mpm/particle:mp-stress mp) 3 0))
          (cl-mpm/output::save-parameter "sig_zx" (magicl:tref (cl-mpm/particle:mp-stress mp) 4 0)))

        (cl-mpm/output::save-parameter "eps_xx" (magicl:tref (cl-mpm/particle:mp-strain mp) 0 0))
        (cl-mpm/output::save-parameter "eps_yy" (magicl:tref (cl-mpm/particle:mp-strain mp) 1 0))
        (cl-mpm/output::save-parameter "eps_xy" (magicl:tref (cl-mpm/particle:mp-strain mp) 5 0))
        (cl-mpm/output::save-parameter "eps_1"
                        (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-strain mp)))
                          (loop for sii in l maximize sii)))
        (when (= 3 nd)
          (cl-mpm/output::save-parameter "eps_zz" (magicl:tref (cl-mpm/particle:mp-strain mp) 2 0))
          (cl-mpm/output::save-parameter "eps_yz" (magicl:tref (cl-mpm/particle:mp-strain mp) 3 0))
          (cl-mpm/output::save-parameter "eps_zx" (magicl:tref (cl-mpm/particle:mp-strain mp) 4 0)))

        ;; (cl-mpm/output::save-parameter "temp" (magicl:tref (cl-mpm/particle::mp-velocity-rate mp) 2 0))

        (cl-mpm/output::save-parameter "damage-inc-average"
                                       (if (slot-exists-p mp 'cl-mpm/particle::time-averaged-damage-inc)
                                           (let ((v (/ (cl-mpm/particle::mp-time-averaged-damage-inc mp)
                                                       (max 1d0
                                                            (cl-mpm/particle::mp-time-averaged-counter mp)))))
                                             (setf (cl-mpm/particle::mp-time-averaged-damage-inc mp) 0d0)
                                             v)
                                           0d0
                                           ))
        (cl-mpm/output::save-parameter "damage-ybar-average"
                                       (if (slot-exists-p mp 'cl-mpm/particle::time-averaged-damage-inc)
                                           (let ((v (/ (cl-mpm/particle::mp-time-averaged-ybar mp)
                                                       (max 1d0
                                                            (cl-mpm/particle::mp-time-averaged-counter mp)))))
                                             (setf (cl-mpm/particle::mp-time-averaged-counter mp) 0d0
                                                   (cl-mpm/particle::mp-time-averaged-ybar mp) 0d0)
                                             v)
                                           0))
        ;; (cl-mpm/output::save-parameter "erosion"
        ;;                                (if (slot-exists-p mp 'cl-mpm/particle::eroded-volume)
        ;;                                    (cl-mpm/particle::mp-eroded-volume mp)
        ;;                                    0))
        (cl-mpm/output::save-parameter "pressure" (cl-mpm/particle::mp-pressure mp))
        (cl-mpm/output::save-parameter "boundary" (cl-mpm/particle::mp-boundary mp))

        (cl-mpm/output::save-parameter "s_1"
                                       (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))

                                         (nth 0 (sort l #'>))))
        (cl-mpm/output::save-parameter "s_3"
                                       (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))

                                         (nth 2 (sort l #'>))))

        (cl-mpm/output::save-parameter "su_1"
                                       (if (typep mp 'cl-mpm/particle::particle-damage)
                                           (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle::mp-undamaged-stress mp)))

                                             (nth 0 (sort l #'>)))
                                           0d0))
        (cl-mpm/output::save-parameter "su_3"
                                       (if (typep mp 'cl-mpm/particle::particle-damage)
                                           (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle::mp-undamaged-stress mp)))

                                             (nth 2 (sort l #'>)))
                                           0d0))
        ;; (cl-mpm/output::save-parameter "e_1"
        ;;                                (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-strain mp)))
        ;;                                  (loop for sii in l maximize sii)))


        (cl-mpm/output::save-parameter "EPS"
                                       (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voight-to-matrix (cl-mpm/particle:mp-stress mp)))
                                         (- (loop for sii in l maximize sii) (cl-mpm/particle::mp-pressure mp))))
        (cl-mpm/output::save-parameter "size_x" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 0 0))
        (cl-mpm/output::save-parameter "size_y" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0))
        (when (= 3 nd)
          (cl-mpm/output::save-parameter "size_z" (magicl:tref (cl-mpm/particle::mp-domain-size mp) 2 0)))

        (cl-mpm/output::save-parameter "damage"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage)
                                           (cl-mpm/particle:mp-damage mp)
                                           0d0))

        (cl-mpm/output::save-parameter "damage-shear"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-shear)
                                           (cl-mpm/particle::mp-damage-shear mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-compression"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-compression)
                                           (cl-mpm/particle::mp-damage-compression mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-tension"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-tension)
                                           (cl-mpm/particle::mp-damage-tension mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-inc"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-increment)
                                           (cl-mpm/particle::mp-damage-increment mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-ybar"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-ybar)
                                           (cl-mpm/particle::mp-damage-ybar mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-ybar-prev"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-ybar-prev)
                                           (cl-mpm/particle::mp-damage-ybar-prev mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-ybar-scaled"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-ybar)
                                           (*
                                            (- 1d0 (cl-mpm/particle::mp-damage mp))
                                            (cl-mpm/particle::mp-damage-ybar mp)
                                            )
                                           0d0))
        (cl-mpm/output::save-parameter "damage-y"
                                       (if (slot-exists-p mp 'cl-mpm/particle::damage-y-local)
                                           (cl-mpm/particle::mp-damage-y-local mp)
                                           0d0))
        (cl-mpm/output::save-parameter "damage-k"
                                       (if (slot-exists-p mp 'cl-mpm/particle::history-stress)
                                           (cl-mpm/particle::mp-history-stress mp)
                                           0d0))
        (cl-mpm/output::save-parameter "local-length"
                                       (if (slot-exists-p mp 'cl-mpm/particle::true-local-length)
                                           (cl-mpm/particle::mp-true-local-length mp)
                                           0d0))
        (cl-mpm/output::save-parameter "p-undamaged"
                                       (if (slot-exists-p mp 'cl-mpm/particle::undamaged-stress)
                                           (cl-mpm/utils::trace-voigt (cl-mpm/particle::mp-undamaged-stress mp))
                                           (cl-mpm/utils::trace-voigt (cl-mpm/particle::mp-stress mp))))
        (cl-mpm/output::save-parameter "p-damaged"
                                       (if (slot-exists-p mp 'cl-mpm/particle::undamaged-stress)
                                           (cl-mpm/utils::trace-voigt (cl-mpm/particle::mp-stress mp))
                                           0d0))

        (cl-mpm/output::save-parameter "q-undamaged"
                                       (if (slot-exists-p mp 'cl-mpm/particle::undamaged-stress)
                                           (sqrt (cl-mpm/constitutive::voigt-j2 (cl-mpm/utils:deviatoric-voigt (cl-mpm/particle::mp-undamaged-stress mp))))
                                           (sqrt (cl-mpm/constitutive::voigt-j2 (cl-mpm/utils:deviatoric-voigt (cl-mpm/particle::mp-stress mp))))))
        (cl-mpm/output::save-parameter "q-damaged"
                                       (if (slot-exists-p mp 'cl-mpm/particle::undamaged-stress)
                                           (sqrt (cl-mpm/constitutive::voigt-j2 (cl-mpm/utils:deviatoric-voigt (cl-mpm/particle::mp-stress mp))))
                                           0d0))

        (cl-mpm/output::save-parameter "split-depth"
                                       (cl-mpm/particle::mp-split-depth mp))

        (cl-mpm/output::save-parameter "plastic-iterations"
                                       (if (slot-exists-p mp 'cl-mpm/particle::plastic-iterations)
                                           (cl-mpm/particle::mp-plastic-iterations mp) 0d0))
        (cl-mpm/output::save-parameter
         "plastic_strain"
         (if (slot-exists-p mp 'cl-mpm/particle::yield-func)
             (cl-mpm/particle::mp-strain-plastic-vm mp)
             0d0)
         )
        (cl-mpm/output::save-parameter
         "yield-func"
         (if (slot-exists-p mp 'cl-mpm/particle::yield-func)
             (cl-mpm/particle::mp-yield-func mp)
             0d0))

        (cl-mpm/output::save-parameter
         "current-effective-angle"
         (if (typep mp 'cl-mpm/particle::particle-ice-delayed)
             (* (/ 180 pi) (atan (* (/ (- 1d0 (cl-mpm/particle::mp-damage-shear mp))
                                       (- 1d0 (cl-mpm/particle::mp-damage-compression mp)))
                                    (tan (cl-mpm/particle::mp-phi mp)))))
             0d0)
         )
        (cl-mpm/output::save-parameter
         "current-effective-coheasion"
         (if (typep mp 'cl-mpm/particle::particle-ice-delayed)
             (*
              (- 1d0 (cl-mpm/particle::mp-damage-compression mp))
              (cl-mpm/particle::mp-c mp))
             0d0)
         )
        ;; (cl-mpm/output::save-parameter
        ;;  "plastic-c"
        ;;  (if (slot-exists-p mp 'cl-mpm/particle::c)
        ;;      (cl-mpm/particle::mp-c mp)
        ;;      0d0))
        ;; (cl-mpm/output::save-parameter
        ;;  "plastic-phi"
        ;;  (if (slot-exists-p mp 'cl-mpm/particle::phi)
        ;;      (* (cl-mpm/particle::mp-phi mp) (/ 180 pi))
        ;;      0d0))

        (cl-mpm/output::save-parameter
         "energy"
         ;; (cl-mpm/particle::mp-penalty-energy mp)
         (* (cl-mpm/particle::mp-mass mp)
            (cl-mpm/fastmaths::mag-squared (cl-mpm/particle::mp-velocity mp)))
         )

        ;; (cl-mpm/output::save-parameter
        ;;  "fbar-j"
        ;;  (cl-mpm/particle::mp-debug-j mp)
        ;;  )
        ;; (cl-mpm/output::save-parameter
        ;;  "fbar-j-gather"
        ;;  (cl-mpm/particle::mp-debug-j-gather mp)
        ;;  )
        ;; (cl-mpm/output::save-parameter
        ;;  "fbar-j-diff"
        ;;  (if (> (cl-mpm/particle::mp-debug-j mp) 0d0) 
        ;;      (/ (- (cl-mpm/particle::mp-debug-j-gather mp) (cl-mpm/particle::mp-debug-j mp)) 
        ;;         (cl-mpm/particle::mp-debug-j mp)
        ;;         )
        ;;      0d0)
        ;;  )
        (cl-mpm/output::save-parameter "erosion"
                                       (/
                                        (cl-mpm/output::optional-slot-access 'cl-mpm/particle::eroded-volume mp)
                                        (cl-mpm/particle::mp-mass mp)))
        (cl-mpm/output::save-parameter "fric-contact" (if (cl-mpm/particle::mp-penalty-contact-step mp) 1 0))
        (cl-mpm/output::save-parameter "fric-contact-stick" (if (cl-mpm/particle::mp-penalty-friction-stick mp) 1 0))
        (cl-mpm/output::save-parameter "fric-normal" (cl-mpm/particle::mp-penalty-normal-force mp))
        (cl-mpm/output::save-parameter "fric-x" (magicl:tref (cl-mpm/particle::mp-penalty-frictional-force mp) 0 0))
        (cl-mpm/output::save-parameter "fric-y" (magicl:tref (cl-mpm/particle::mp-penalty-frictional-force mp) 1 0))
        (cl-mpm/output::save-parameter "fric-k" (cl-mpm/particle::mp-penalty-stiffness mp))
        )
      )))
