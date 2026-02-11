(defpackage :cl-mpm/inflow
  (:use :cl
        :cl-mpm
        :cl-mpm/particle
        :cl-mpm/mesh
        :cl-mpm/utils
        :cl-mpm/fastmaths
        )
  (:import-from
    :magicl tref .+ .-)
  (:import-from
   :cl-mpm/fastmaths  fast-.+ fast-.- fast-.*)
  (:export
   #:make-bc-inflow))
(declaim #.cl-mpm/settings:*optimise-setting*)
(in-package :cl-mpm/inflow)
(defclass bc-true-inflow ()
  ((sim
    :accessor bc-sim
    :initarg :sim
    )
   (particle-constructor
    :accessor bc-inflow-particle-constructor
    :initarg :particle-constructor
    ))
  (:documentation "A true inflow algorithm"))



(defclass bc-mesh-inflow (bc-true-inflow)
   ((cell-list
     :accessor bc-mesh-inflow-cell-list)
    (current-disp
     :accessor bc-current-disp
     :initform 0d0)
    (velocity
     :accessor bc-inflow-velocity
     :initform (cl-mpm/utils:vector-zeros)))
  (:documentation "A mesh based inflow condition"))

(defun make-bc-mesh-inflow (sim cell-list inflow-velocity constructor)
  (make-instance 'bc-mesh-inflow
                 :index nil
                 :sim sim
                 :cell-list cell-list
                 :velocity inflow-velocity
                 :particle-constructor constructor)
  (dolist (cell cell-list)
    (cl-mpm/bc::make-bc-fixed )
    )
  )

(defmethod cl-mpm/bc::apply-bc ((bc bc-mesh-inflow) node mesh dt)
"Nodal inflow"
  (with-accessors ((sim bc-sim)
                   ;(volume-threshold volume-threshold)
                   (pcons bc-inflow-particle-constructor)
                   (cell-list bc-mesh-inflow-cell-list)
               )
      bc
    (loop for c in cell-list
          )
    (when t
      (print "Created particle")
      ;; Add particle?
      (funcall pcons (list (magicl:tref (cl-mpm/mesh::node-position node) 0 0)
                           (magicl:tref (cl-mpm/mesh::node-position node) 1 0)
                           ))))
  )

(defclass bc-mps-inflow (bc-true-inflow)
  ((mps-list
    :accessor bc-inflow-mps-list
    :mps)
   (current-disp
    :accessor bc-current-disp
    :initform 0d0)
   (velocity
    :accessor bc-inflow-velocity
    :initform (cl-mpm/utils:vector-zeros)))
  (:documentation "An MP based inflow condition"))

(defun make-bc-mps-inflow (sim mps-list inflow-velocity constructor)
  (make-instance 'bc-mesh-inflow
                 :index nil
                 :sim sim
                 :mps (make-array (length mps-list) :initial-contents mps-list)
                 :velocity inflow-velocity
                 :particle-constructor constructor)
  ;; (dolist (cell cell-list)
  ;;   (cl-mpm/bc::make-bc-fixed
  ;;    )
  ;;   )
  )

(defmethod cl-mpm/bc::apply-bc ((bc bc-mps-inflow) node mesh dt)
  "Nodal inflow"
  (with-accessors ((sim bc-sim)
                   (pcons bc-inflow-particle-constructor)
                   (mps-list bc-inflow-mps-list)
                   )
      bc
    (cl-mpm::iterate-over-mps
     mps-list
     (lambda (mp)

       ))
    (when t
      ;; (print "Created particle")
      ;; Add particle?
      ;; (funcall pcons (list (magicl:tref (cl-mpm/mesh::node-position node) 0 0)
      ;;                      (magicl:tref (cl-mpm/mesh::node-position node) 1 0)
      ;;                      ))
      ))
  )
