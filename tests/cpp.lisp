(in-package :cl-mpm-tests)

(defun test-ext-plastic-dp (&rest strain-list)
  (let* ((E 1d0)
         (nu 0.1d0)
         (phi 0.1d0)
         (psi 0.05d0)
         (c 1d0)
         (strain (cl-mpm/constitutive::swizzle-coombs->voigt (voigt-from-list strain-list)))
         (de (cl-mpm/constitutive:linear-elastic-matrix E nu)))
    (multiple-value-bind (sig eps f) (cl-mpm/ext::constitutive-drucker-prager strain de E nu phi psi c)
      (cl-mpm/constitutive::swizzle-voigt->coombs eps))))

(defun test-ext-plastic-mc (&rest strain-list)
  (let* ((E 1d0)
         (nu 0.1d0)
         (phi 0.1d0)
         (psi 0.05d0)
         (c 1d0)
         (strain (cl-mpm/constitutive::swizzle-coombs->voigt (voigt-from-list strain-list)))
         (de (cl-mpm/constitutive:linear-elastic-matrix E nu)))
    (multiple-value-bind (sig eps f) (cl-mpm/ext::constitutive-drucker-prager strain de E nu phi psi c)
      (cl-mpm/constitutive::swizzle-voigt->coombs eps))))

(deftest test-ext-drucker-prager ()
  ;;Not a good test
  (is (voigt-list-test (test-ext-plastic-dp 0d0 0d0 0d0 0d0 0d0 0d0) (list 0d0 0d0 0d0 0d0 0d0 0d0)))
  (is (voigt-list-test (test-ext-plastic-dp 1d0 0d0 0d0 0d0 0d0 0d0) (list 1d0 0d0 0d0 0d0 0d0 0d0)))
  (is (voigt-list-test (test-ext-plastic-dp -1d0 0d0 0d0 0d0 0d0 0d0) (list -1d0 0d0 0d0 0d0 0d0 0d0)))
  (is (voigt-list-test (test-ext-plastic-dp 1d0 2d0 3d0 1d0 2d0 3d0) (list 1.475498236391475d0 1.954875103787714d0 2.434251971183953d0 0.479376867396238d0 0.958753734792476d0 1.438130602188714d0)))
  (is (voigt-list-test (test-ext-plastic-dp 0d0 0d0 0d0 3d0 0d0 0d0)
                       (list -0.006199645592901d0 -0.006199645592901d0 -0.006199645592899d0 2.696533775890750d0 0d0 0d0))))



(deftest test-ext-mohr-coloumb ()
  ;;Not a good test
  (is (voigt-list-test (test-ext-plastic-mc 0d0 0d0 0d0 0d0 0d0 0d0) (list 0d0 0d0 0d0 0d0 0d0 0d0)))
  (is (voigt-list-test (test-ext-plastic-mc 1d0 0d0 0d0 0d0 0d0 0d0) (list 1d0 0d0 0d0 0d0 0d0 0d0)))
  (is (voigt-list-test (test-ext-plastic-mc -1d0 0d0 0d0 0d0 0d0 0d0) (list -1d0 0d0 0d0 0d0 0d0 0d0)))
  (is (voigt-list-test (test-ext-plastic-mc -5d0 0d0 0d0 0d0 0d0 0d0) (list -3.5522d0 -0.8001d0 -0.8001d0 0d0 0d0 0d0)))
  (is (voigt-list-test (test-ext-plastic-mc 1d0 2d0 3d0 1d0 2d0 3d0)
                       (list
                        1.6919d0
                        1.7525d0
                        2.4287d0
                        0.5589d0
                        1.0592d0
                        1.0905d0)))

  (is (voigt-list-test (test-ext-plastic-mc 0d0 0d0 0d0 3d0 0d0 0d0)
                       (list
                        -0.0201d0
                        -0.0201d0
                        0.0000d0
                        2.1940d0
                        0d0
                        0d0
                        )
                       )))
