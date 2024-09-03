(in-package :cl-mpm-tests)

(defun test-plastic-dp (&rest strain-list)
  (let* ((E 1d0)
         (nu 0.1d0)
         (phi 0.1d0)
         (psi 0d0)
         (c 1d0)
         (strain (voigt-from-list strain-list))
         (de (cl-mpm/constitutive:linear-elastic-matrix E nu))
         (stress (magicl:@ de strain)))
    (multiple-value-bind (sig eps f) (cl-mpm/constitutive::plastic-dp stress de strain E nu phi psi c)
      eps)))
(defun test-plastic-mc (&rest strain-list)
  (let* ((E 1d0)
         (nu 0.1d0)
         (phi 0.1d0)
         (psi 0d0)
         (c 1d0)
         (strain (voigt-from-list strain-list))
         (de (cl-mpm/constitutive:linear-elastic-matrix E nu))
         (stress (magicl:@ de strain)))
    (multiple-value-bind (sig eps f) (cl-mpm/constitutive::mc-plastic stress de strain E nu phi psi c)
      eps)))

(defun test (strain-list result-list)
  (let* ((out (apply #'test-plastic-dp strain-list))
         (err
           (cl-mpm/fastmath:fast-.-
            (voigt-from-list result-list)
            out))
         (norm
           (cl-mpm/fastmath:mag err)))
    (pprint norm)
    (if (< norm 1d-4)
        (format t "~%Pass:~%~A~%~A~%~A~%" result-list (loop for i from 0 below 6 collect (varef out i))
                (loop for i from 0 below 6 collect (varef err i)))
        (format t "~%Fail:~%~A~%~A~%~A~%" result-list (loop for i from 0 below 6 collect (varef out i))
                (loop for i from 0 below 6 collect (varef err i))))))

(defun dp-test (strain-list result-list)
  (let* ((out (apply #'test-plastic-dp strain-list))
         (err (cl-mpm/fastmath:fast-.-
               (voigt-from-list result-list)
               out))
         (norm (cl-mpm/fastmath:mag err)))
    (< norm 1d-4)))
(defun mc-test (strain-list result-list)
  (let* ((out (apply #'test-plastic-mc strain-list))
         (err (cl-mpm/fastmath:fast-.-
               (voigt-from-list result-list)
               out))
         (norm (cl-mpm/fastmath:mag err)))
    (< norm 1d-4)))

;; (pprint (test-plastic-dp 1d0 0d0 0d0 0d0 0d0 0d0))
;; (pprint (test-plastic-dp -1d0 0d0 0d0 0d0 0d0 0d0))
;; (pprint (test-plastic-dp 1d0 2d0 3d0 0d0 0d0 0d0))
;; (pprint (test-plastic-dp 1d0 2d0 3d0 1d0 2d0 3d0))

;; (test (list 1d0 0d0 0d0 0d0 0d0 0d0) (list 1d0 0d0 0d0 0d0 0d0 0d0))
;; (test (list -1d0 0d0 0d0 0d0 0d0 0d0) (list -1d0 0d0 0d0 0d0 0d0 0d0))
;; (test (list 1d0 2d0 3d0 1d0 2d0 3d0) (list
;;                                       1.5242d0
;;                                       2.0000d0
;;                                       2.4758d0
;;                                       0.4758d0
;;                                       0.9516d0
;;                                       1.4273d0))

;; (test (list 0d0 0d0 0d0 3d0 0d0 0d0) (list
;;                                       0d0
;;                                       0d0
;;                                       0d0
;;                                       2.6944d0
;;                                       0d0
;;                                       0d0
;;                                       ))


(deftest test-drucker-prager ()
  ;;Not a good test
  (is (dp-test (list 1d0 0d0 0d0 0d0 0d0 0d0) (list 1d0 0d0 0d0 0d0 0d0 0d0)))
  (is (dp-test (list -1d0 0d0 0d0 0d0 0d0 0d0) (list -1d0 0d0 0d0 0d0 0d0 0d0)))
  (is (dp-test (list 1d0 2d0 3d0 1d0 2d0 3d0) (list 1.5242d0 2.0000d0 2.4758d0 0.4758d0 0.9516d0 1.4273d0)))
  (is (dp-test (list 0d0 0d0 0d0 3d0 0d0 0d0) (list 0d0 0d0 0d0 2.6944d0 0d0 0d0)))
  )

(deftest test-mohr-coloumb ()
  ;;Not a good test
  (is (mc-test (list 1d0 0d0 0d0 0d0 0d0 0d0) (list 1d0 0d0 0d0 0d0 0d0 0d0)))
  (is (mc-test (list -1d0 0d0 0d0 0d0 0d0 0d0) (list -1d0 0d0 0d0 0d0 0d0 0d0)))
  (is (mc-test (list 1d0 2d0 3d0 1d0 2d0 3d0) (list 1.5242d0 2.0000d0 2.4758d0 0.4758d0 0.9516d0 1.4273d0)))
  (is (mc-test (list 0d0 0d0 0d0 3d0 0d0 0d0) (list 0d0 0d0 0d0 2.6944d0 0d0 0d0)))
  )
