
(defpackage :cl-mpm/models/visco
  (:use :cl :cl-mpm/utils :cl-mpm/fastmaths)
  (:export
   )
  )
(in-package :cl-mpm/particle)
(defclass particle-finite-viscoelastic (particle-elastic)
  ((viscosity
    :accessor mp-viscosity
    :initarg :viscosity
    :initform 1d0
    )))
(defmethod constitutive-model ((mp particle-finite-viscoelastic) strain dt)
  (with-accessors ((de mp-elastic-matrix)
                   (stress mp-stress)
                   (strain mp-strain)
                   (strain-n mp-strain-n)
                   (viscosity mp-viscosity)
                   )
      mp
    ;;Train elastic strain - plus trail kirchoff stress
    ;; (setf stress (cl-mpm/constitutive::linear-elastic-mat strain de))
    (let* ((etr strain)
           (en strain)
           (f-tol 1d-5)
           (tau nil)
           (a
             (magicl:inv
              (magicl:.+
               (magicl:inv de)
               (magicl:scale
                (magicl:.- (magicl:eye 6)
                           (magicl:scale
                            (magicl:from-list
                             (list 1d0 1d0 1d0 0d0 0d0 0d0
                                   1d0 1d0 1d0 0d0 0d0 0d0
                                   1d0 1d0 1d0 0d0 0d0 0d0
                                   0d0 0d0 0d0 0d0 0d0 0d0
                                   0d0 0d0 0d0 0d0 0d0 0d0
                                   0d0 0d0 0d0 0d0 0d0 0d0
                                   )
                             '(6 6)) (/ 1d0 3d0)))
                (/ dt (* 2d0 viscosity))))))
           (f f-tol))
      (loop for i from 0 to 50
            while (>= f f-tol)
            do
               (progn
                 ;; (pprint i)
                 (setf tau (cl-mpm/constitutive::linear-elastic-mat en de tau))
                 (let* ((r (magicl:.-
                            (magicl:.+ en
                                       (magicl:scale
                                        (cl-mpm/utils:deviatoric-voigt tau)
                                        (/ dt (* 2d0 viscosity))))
                            etr)))
                   (setf f (magicl:norm r))
                   ;; (format t "f ~F ~%"f)
                   ;; (pprint r)
                   (when (>= f f-tol)
                     ;; (pprint  r)
                     ;; (pprint (magicl:linear-solve a r))
                     ;; (format t "~%")
                     (setf
                      en
                      (magicl:.+
                       en
                       (magicl:.*
                        (magicl:linear-solve a (magicl:scale r -1))
                        (cl-mpm/utils:voigt-from-list (list 1d0 1d0 1d0 0.5d0 0.5d0 0.5d0))))))
                   ))
               )
      (cl-mpm/utils:voigt-copy-into tau stress)
      (cl-mpm/utils:voigt-copy-into en strain)
      )))


;; (let ((test (cl-mpm/utils:voigt-from-list (list 0d0 2d0 3d0 2d0 0d0 3d0))))
;;   (pprint (cl-mpm/utils:deviatoric-voigt test))
;;   (pprint (magicl:@ (magicl:.- (magicl:eye 6)
;;                                (magicl:scale
;;                                 (magicl:from-list
;;                                  (list 1d0 1d0 1d0 0d0 0d0 0d0
;;                                        1d0 1d0 1d0 0d0 0d0 0d0
;;                                        1d0 1d0 1d0 0d0 0d0 0d0
;;                                        0d0 0d0 0d0 0d0 0d0 0d0
;;                                        0d0 0d0 0d0 0d0 0d0 0d0
;;                                        0d0 0d0 0d0 0d0 0d0 0d0
;;                                        )
;;                                  '(6 6)) (/ 1d0 3d0))
;;                                ) test))
;;   )

(let* ((de (cl-mpm/constitutive:linear-elastic-matrix 1d0 0d0))
       (strain (cl-mpm/utils:voigt-from-list (list 1d0 0d0 0d0 0d0 2d0 0d0)))
       (viscosity 1d0)
       (dt 0.1d0)
       (etr strain)
       (en strain)
       (f-tol 1d-5)
       (tau nil)
       (a
         (magicl:inv
          (magicl:.+
           (magicl:inv de)
           (magicl:scale
            (magicl:.- (magicl:eye 6)
                       (magicl:scale
                        (magicl:from-list
                         (list 1d0 1d0 1d0 0d0 0d0 0d0
                               1d0 1d0 1d0 0d0 0d0 0d0
                               1d0 1d0 1d0 0d0 0d0 0d0
                               0d0 0d0 0d0 0d0 0d0 0d0
                               0d0 0d0 0d0 0d0 0d0 0d0
                               0d0 0d0 0d0 0d0 0d0 0d0
                               )
                         '(6 6)) (/ 1d0 3d0)))
            (/ dt (* 2d0 viscosity))))))
       (f f-tol))
  (pprint "start")
  (pprint a)
  (loop for i from 0 to 5
        while (>= f f-tol)
        do
           (progn
             ;; (pprint i)
             (setf tau (cl-mpm/constitutive::linear-elastic-mat en de tau))
             (let* ((r (magicl:.-
                        (magicl:.+ en
                                   (magicl:scale
                                    (cl-mpm/utils:deviatoric-voigt tau)
                                    (/ dt (* 2d0 viscosity))))
                        etr)))
               (setf f (magicl:norm r))
               (format t "f ~F ~%"f)
               ;; (pprint r)
               (when (>= f f-tol)
                 (pprint  r)
                 (pprint (magicl:linear-solve a r))
                 (format t "~%")
                 (setf
                  en
                  (magicl:.+
                   en
                   (magicl:.*
                    (magicl:linear-solve a (magicl:scale r -1))
                    (cl-mpm/utils:voigt-from-list (list 1d0 1d0 1d0 0.5d0 0.5d0 0.5d0))))))
               )))
  (when (< f f-tol)
    (format t "Converged!~%"))
  (pprint en)
  ;; (cl-mpm/utils:voigt-copy-into tau stress)
  )
