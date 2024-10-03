(in-package :cl-mpm-tests)

(defun test-plastic-dp (&rest strain-list)
  (let* ((E 1d0)
         (nu 0.1d0)
         (phi 0.1d0)
         (psi 0.05d0)
         (c 1d0)
         (strain (cl-mpm/constitutive::swizzle-coombs->voigt (voigt-from-list strain-list)))
         (de (cl-mpm/constitutive:linear-elastic-matrix E nu))
         (stress (magicl:@ de strain)))
    (multiple-value-bind (sig eps f) (cl-mpm/constitutive::plastic-dp stress de strain E nu phi psi c)
      (cl-mpm/constitutive::swizzle-voigt->coombs eps))))
(defun test-plastic-mc (&rest strain-list)
  (let* ((E 1d0)
         (nu 0.1d0)
         (phi 0.1d0)
         (psi 0.05d0)
         (c 0.5d0)
         (strain (cl-mpm/constitutive::swizzle-coombs->voigt (voigt-from-list strain-list)))
         (de (cl-mpm/constitutive:linear-elastic-matrix E nu))
         (stress (magicl:@ de strain)))
    (multiple-value-bind (sig eps f) (cl-mpm/constitutive::mc-plastic stress de strain E nu phi psi c)
      (cl-mpm/constitutive::swizzle-voigt->coombs eps))))

(defun test-plastic-vm (&rest strain-list)
  (let* ((E 1d0)
         (nu 0.1d0)
         (rho 10d0)
         (strain (cl-mpm/constitutive::swizzle-coombs->voigt (voigt-from-list strain-list)))
         (de (cl-mpm/constitutive:linear-elastic-matrix E nu))
         (stress (magicl:@ de strain)))
    (multiple-value-bind (sig eps f) (cl-mpm/constitutive::vm-plastic stress de strain rho)
      (cl-mpm/constitutive::swizzle-voigt->coombs eps))))

(defun test (strain-list result-list)
  (let* ((out (apply #'test-plastic-dp strain-list))
         (err
           (cl-mpm/fastmaths:fast-.-
            (voigt-from-list result-list)
            out))
         (norm
           (cl-mpm/fastmaths:mag err)))
    (pprint norm)
    (if (< norm 1d-4)
        (format t "~%Pass:~%~A~%~A~%~A~%" result-list (loop for i from 0 below 6 collect (varef out i))
                (loop for i from 0 below 6 collect (varef err i)))
        (format t "~%Fail:~%~A~%~A~%~A~%" result-list (loop for i from 0 below 6 collect (varef out i))
                (loop for i from 0 below 6 collect (varef err i))))))

(defun voigt-test (a b)
  (let* ((err (cl-mpm/fastmaths:fast-.-
               a b))
         (norm (cl-mpm/fastmaths:mag err)))
    (< norm 1d-4)))

(defun voigt-list-test (a b)
  (let* ((err (cl-mpm/fastmaths:fast-.-
               a (cl-mpm/utils:voigt-from-list b)))
         (norm (cl-mpm/fastmaths:mag err)))
    (< norm 1d-4)))

(defun dp-test (strain-list result-list)
  (let* ((out (apply #'test-plastic-dp strain-list))
         (err (cl-mpm/fastmaths:fast-.-
               (cl-mpm/constitutive::swizzle-voigt->coombs (voigt-from-list result-list)) out))
         (norm (cl-mpm/fastmaths:mag err)))
    (< norm 1d-4)))
(defun mc-test (strain-list result-list)
  (let* ((out (apply #'test-plastic-mc strain-list))
         (err (cl-mpm/fastmaths:fast-.-
               (cl-mpm/constitutive::swizzle-voigt->coombs (voigt-from-list result-list))
               out))
         (norm (cl-mpm/fastmaths:mag err)))
    (< norm 1d-4)))

(defun vm-test (strain-list result-list)
  (let* ((out (apply #'test-plastic-vm strain-list))
         (err (cl-mpm/fastmaths:fast-.-
               (cl-mpm/constitutive::swizzle-voigt->coombs (voigt-from-list result-list))
               out))
         (norm (cl-mpm/fastmaths:mag err)))
    (< norm 1d-4)))

(deftest test-drucker-prager ()
  ;;Not a good test
  (is (voigt-list-test (test-plastic-dp 1d0 0d0 0d0 0d0 0d0 0d0) (list 1d0 0d0 0d0 0d0 0d0 0d0)))
  (is (voigt-list-test (test-plastic-dp -1d0 0d0 0d0 0d0 0d0 0d0) (list -1d0 0d0 0d0 0d0 0d0 0d0)))
  (is (voigt-list-test (test-plastic-dp 1d0 2d0 3d0 1d0 2d0 3d0) (list 1.475498236391475d0 1.954875103787714d0 2.434251971183953d0 0.479376867396238d0 0.958753734792476d0 1.438130602188714d0)))
  (is (voigt-list-test (test-plastic-dp 0d0 0d0 0d0 3d0 0d0 0d0)
               (list -0.006199645592901d0 -0.006199645592901d0 -0.006199645592899d0 2.696533775890750d0 0d0 0d0)
               ;; (list 0d0 0d0 0d0 2.6944d0 0d0 0d0)
               ))
  )

(deftest test-mohr-coloumb ()
  ;;Not a good test
  (is (voigt-list-test (test-plastic-mc 1d0 0d0 0d0 0d0 0d0 0d0)
                       (list
                        0.979934879811914d0
                        0.009077457299942d0
                        0.009077457299942d0
                        0d0
                        0d0
                        0d0
                        )
                      ))
  (is (voigt-list-test (test-plastic-mc -1d0 0d0 0d0 0d0 0d0 0d0) (list -1d0 0d0 0d0 0d0 0d0 0d0)))
  (is (voigt-test (test-plastic-mc 1d0 2d0 3d0 1d0 2d0 3d0)
       (cl-mpm/utils:voigt-from-list
        (list
         1.851075535530768d0
         1.855075170381020d0
         2.099777263559730d0
         0.198854005878051d0
         0.374783968105772d0
         0.367321581130859d0
         )))
      )
  )

(deftest test-von-mises ()
  ;;Not a good test
  (is (vm-test (list 1d0 0d0 0d0 0d0 0d0 0d0) (list 1d0 0d0 0d0 0d0 0d0 0d0)))
  (is (vm-test (list -1d0 0d0 0d0 0d0 0d0 0d0) (list -1d0 0d0 0d0 0d0 0d0 0d0)))
  (is (vm-test (list 15d0 0d0 0d0 0d0 0d0 0d0) (list
                                                13.981462390204987d0
                                                0.509268804897507d0
                                                0.509268804897507d0
                                                0d0
                                                0d0
                                                0d0)))
  (is (voigt-test
       (test-plastic-vm 10d0 20d0 30d0 10d0 20d0 30d0)
       (cl-mpm/utils:voigt-from-list
        (list
         16.333333333333332d0
         20.000000000000000d0
         23.666666666666671d0
         3.666666666666669d0
         7.333333333333337d0
         11.000000000000004d0)))))
