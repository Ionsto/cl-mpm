(in-package :cl-mpm/penalty)

(format t "Load octree-penalty package~%")

(defmethod cl-mpm/bc::apply-sim-bc ((sim cl-mpm/dynamic-relaxation::mpm-sim-octree) (bc bc-penalty-structure) dt)
  (with-accessors ((epsilon bc-penalty-epsilon)
                   (friction bc-penalty-friction)
                   (damping bc-penalty-damping)
                   (sub-bcs bc-penalty-structure-sub-bcs)
                   (debug-mutex bc-penalty-load-lock)
                   (debug-force bc-penalty-load)
                   ;; (sim bc-penalty-sim)
                   )
      bc
    (reset-load bc)
    (with-accessors ((mps cl-mpm:sim-mps)
                     (mesh-list cl-mpm::sim-mesh-list))
        sim
      (cl-mpm:iterate-over-mps
       mps
       (lambda (mp)
         (let ((mesh (aref mesh-list (cl-mpm/particle::mp-mesh-index mp))))
           (when t;(early-sweep-intersection bc mp)
             (cl-mpm::iterate-over-midpoints
              ;; cl-mpm::iterate-over-corners
              mesh
              mp
              (lambda (corner-trial)
                (let* ((disp (compute-corner-displacement mesh corner-trial))
                       (corner (cl-mpm/fastmaths:fast-.+ corner-trial disp)))
                  (cl-mpm/mesh::clamp-point-to-bounds mesh corner)
                  (multiple-value-bind (in-contact pen closest-point) (resolve-closest-contact bc corner)
                    ;; (cl-mpm/fastmaths:fast-.- (contact-point closest-point) disp-inc (contact-point closest-point))
                    (when in-contact
                      (setf (contact-trial-point closest-point) corner-trial
                             (contact-disp closest-point) disp)
                      (let ((load (apply-penalty-point
                                   mesh
                                   (contact-sub-bc closest-point)
                                   mp
                                   (contact-point closest-point)
                                   (contact-trial-point closest-point)
                                   (contact-disp closest-point)
                                   dt)))
                        (sb-thread:with-mutex (debug-mutex)
                          (incf debug-force load)))))))))))))))

(defmethod cl-mpm/bc::apply-sim-bc ((sim cl-mpm/dynamic-relaxation::mpm-sim-octree) (bc cl-mpm/penalty::bc-penalty) dt)
  (with-accessors ((epsilon bc-penalty-epsilon)
                   (friction bc-penalty-friction)
                   (normal bc-penalty-normal)
                   (datum bc-penalty-datum)
                   (damping bc-penalty-damping)
                   (sub-bcs bc-penalty-structure-sub-bcs)
                   (debug-mutex bc-penalty-load-lock)
                   (debug-force bc-penalty-load)
                   (mp-stiffness bc-mp-stiffness)
                   (sim bc-penalty-sim))
      bc
    (reset-load bc)
    (setf mp-stiffness nil)
    (with-accessors ((mps cl-mpm:sim-mps)
                     (mesh-list cl-mpm::sim-mesh-list))
        sim
      (cl-mpm:iterate-over-mps
       mps
       (lambda (mp)
         (when t;(early-sweep-intersection bc mp)
           (let ((in-contact nil)
                 (closest-point (make-contact
                                 nil
                                 0d0
                                 0d0
                                 nil
                                 nil)))
             (let ((mesh (aref mesh-list (cl-mpm/particle::mp-mesh-index mp))))
               ;; (pprint (cl-mpm/particle::mp-mesh-index mp))
               (;cl-mpm::iterate-over-corners
                cl-mpm::iterate-over-midpoints
                mesh
                mp
                (lambda (corner-trial)
                  (let* ((disp (compute-corner-displacement mesh corner-trial))
                         (corner (cl-mpm/fastmaths:fast-.+ corner-trial disp)))
                    (cl-mpm/mesh::clamp-point-to-bounds mesh corner)
                    (let* ((penetration-dist (penetration-distance-point corner datum normal)))
                      (declare (double-float penetration-dist))
                      (when (and
                             (>= penetration-dist 0d0)
                             (penalty-contact-valid bc corner))
                        (apply-penalty-point mesh
                                             bc
                                             mp
                                             corner
                                             corner-trial
                                             disp
                                             dt))))))))))))))


(defmethod cl-mpm/bc::assemble-bc-stiffness ((sim cl-mpm/dynamic-relaxation::mpm-sim-octree) (bc cl-mpm/penalty::bc-penalty))
  (let* (;; (mesh (cl-mpm:sim-mesh sim))
         (mps (cl-mpm:sim-mps sim))
         (stiffness-scale (bc-penalty-stiffness-scale bc))
         )
    (let ((contacts (get-contact-points bc)))
      ;; (format t "Contacts ~D~%" (length contacts))
      (when (> (length contacts) 0)
        ;; (format t "Contacts ~D~%" (length (bc-penalty-contact-points bc)))
        (cl-mpm/utils::bpdotimes (i (length contacts))
          (let* ((contact (aref contacts i))
                 (mesh (dr-contact-point-mesh contact))
                 (nd (cl-mpm/mesh::mesh-nd mesh))
                )
            (iterate-over-neighbours-point-linear
             mesh
             (dr-contact-point-position contact)
             (lambda (mesh node svp grads)
               (with-accessors ((node-active cl-mpm/mesh:node-active)
                                (node-volume cl-mpm/mesh::node-volume)
                                (node-mass cl-mpm/mesh::node-mass)
                                (node-lock cl-mpm/mesh::node-lock))
                   node
                 (declare (double-float node-mass node-volume svp))
                 (when node-active
                   (sb-thread:with-mutex (node-lock)
                     (setf
                      node-mass
                      (+
                       node-mass
                       (*
                        stiffness-scale
                        svp
                        (dr-contact-point-stiffness contact)))))))))))))))
