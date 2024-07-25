(defpackage :cl-mpm/setup
  (:use :cl
   :cl-mpm/utils
   )
  (:export
   #:make-column
   #:make-block-mps
   #:remove-sdf
   #:ellipse-sdf
   #:circle-sdf
   #:rectangle-sdf
   #:estimate-elastic-dt
   #:estimate-critical-damping
  ))

(in-package :cl-mpm/setup)




(defun make-block (res element-count &key (shape-maker #'cl-mpm/shape-function::make-shape-function-linear)
                                          (sim-type 'cl-mpm::mpm-sim-usf)
                                       )
  "Make a 2D column of height size, and width 1 - filled with elements"
  (let ((nd (length element-count)))
    (let* ((nD nd)
           (size (mapcar (lambda (x) (* x res)) element-count))
           (sim (cl-mpm:make-mpm-sim size res 1d-3 (funcall shape-maker nD res) :sim-type sim-type)))
      (progn
        (setf (cl-mpm:sim-mps sim) #())
        (setf (cl-mpm:sim-bcs sim) (cl-mpm/bc:make-outside-bc (cl-mpm:sim-mesh sim)))
        sim))))

(defun make-column (height element-count &optional (shape-maker #'cl-mpm::make-shape-function-linear))
  "Make a 2D column of heigh size, and width 1 - filled with elements"
  (let* ((nD 2)
         (mp-spacing (/ height element-count))
         (sim (cl-mpm:make-mpm-sim (list mp-spacing height) mp-spacing 1e-3
                                   nil
                                   ;(funcall shape-maker nD mp-spacing)
                                   )))
    (progn
          (setf (cl-mpm:sim-mps sim) #())
          (setf (cl-mpm:sim-bcs sim) (cl-mpm/bc:make-outside-bc (cl-mpm/mesh:mesh-count (cl-mpm:sim-mesh sim)))) 
           sim)))

(defun make-block-mps-list (offset size mps density constructor &rest args &key (angle 0) (clip-func (lambda (x y z) t))&allow-other-keys)
  (declare (function clip-func))
  "Construct a block of mxn (mps) material points real size (size) with a density (density)"
  (if (= (length size) 2)
      (setf mps (append mps '(1))
            size (append size '(0))
            offset (append offset '(0))))
  (format t "Angle: ~E~%" angle)
  (let*  ((nD 3)
          (args (alexandria:remove-from-plist args :angle :clip-func))
          (spacing (mapcar (lambda (s m) (/ (coerce s 'double-float)
                                            (coerce m 'double-float)
                                            )) size mps))
          (offset (mapcar #'+ offset
                          (mapcar (lambda (x) (* x 0.5d0)) spacing)
                          ))
          (volume (reduce #'* (remove-if (lambda (x) (= 0d0 x)) spacing)))
          (data (loop for x from 0 to (- (nth 0 mps) 1)
                      append
                      (loop
                        for y from 0 to (- (nth 1 mps) 1)
                        append
                        (loop
                          for z from 0 to (- (nth 2 mps) 1)
                          when (funcall clip-func
                                        (+ (* (nth 0 spacing) x) (nth 0 offset))
                                        (+ (* (nth 1 spacing) y) (nth 1 offset))
                                        (+ (* (nth 2 spacing) z) (nth 2 offset)))
                            collect
                            (let* ((rot (cl-mpm/utils::rotation-matrix angle))
                                   (origin-vec (cl-mpm/utils:vector-from-list offset))
                                   (position-vec (magicl:from-list (list (* (nth 0 spacing) x)
                                                                         (* (nth 1 spacing) y)
                                                                         (* (nth 2 spacing) z))
                                                               '(3 1) :type 'double-float))
                               (size-vec (magicl:from-list spacing '(3 1) :type 'double-float))
                               (position-vec (cl-mpm/fastmath::fast-.+ origin-vec
                                                        (magicl:@ rot position-vec)))
                               )
                          (flet ((lisp-list (m) (loop for i from 0 to 2 collect (magicl:tref m i 0))))
                            (apply #'cl-mpm::make-particle
                                   (append (list 3 constructor)
                                           args
                                           (list
                                            :position (lisp-list position-vec)
                                            :mass (* density volume)
                                            :volume volume
                                            :size size-vec
                                            :size-0 size-vec
                                            )))))
                          )))))
    data))
(defun make-mps-from-list (mp-list)
  (make-array (length mp-list)
              :element-type t
              :adjustable t
              :fill-pointer (length mp-list)
              :initial-contents mp-list))

(defun make-block-mps (offset size mps constructor &rest args)
  (if (= (length size) 2)
      (append size '(1)))

  (let*  ((data (apply #'make-block-mps-list offset size mps constructor args)))
    (make-mps-from-list data)))

(defun make-block-mps-clipped (offset size mps clip-func constructor &rest args)
  (if (= (length size) 2)
      (append size '(1)))

  (let*  ((data (apply #'make-block-mps-list offset size mps constructor args)))
    (make-mps-from-list data)))


(defun make-column-mps (size mp-spacing constructor &rest args)
  (let* ((mp-spacing-x (first mp-spacing))
        (mp-spacing-y (second mp-spacing))
        (data (loop for i from 0 to (- size 1) collect 
                    (apply constructor (append '(2) args 
                                               (list
                                                 :position (list (/ mp-spacing-x 2d0)
                                                            (+ (/ mp-spacing-y 2d0) (* mp-spacing-y i)))
                                                 :volume (* mp-spacing-x mp-spacing-y)
                                                 :size (magicl:from-list (list mp-spacing-x mp-spacing-y) '(2 1))))))))
    (make-array (length data) :initial-contents data)))

(defun make-column-mps-elastic (element-count spacing E nu)
  (make-column-mps element-count spacing 'cl-mpm::make-particle-elastic E nu))

(defun node-to-mp-pos (position res mps)
  (mapcar #'+ (list (* res (+ (/ 1 (* 2 mps))))
                    (* res (+ (/ 1 (* 2 mps)))))
                      position))


(defun remove-sdf (sim sdf)
  (setf (cl-mpm:sim-mps sim)
        (delete-if (lambda (mp)
                     (with-accessors ((pos cl-mpm/particle:mp-position)) mp
                       (>= 0 (funcall sdf pos))
                       ))
                   (cl-mpm:sim-mps sim))))

(defun damage-sdf (sim sdf &optional (d 1d0))
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (loop for mp across mps
          do (with-accessors ((pos cl-mpm/particle:mp-position)
                              (damage cl-mpm/particle:mp-damage)) mp
               (when (>= 0 (funcall sdf pos))
                 (setf damage (coerce d 'double-float)))
               ))))
(defun apply-sdf (sim sdf func)
  (declare (function sdf func))
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (loop for mp across mps
          do (with-accessors ((pos cl-mpm/particle:mp-position)
                              (damage cl-mpm/particle:mp-damage)) mp
               (let ((sdf-v (funcall sdf pos)))
                 (when (>= 0 sdf-v)
                   (funcall func mp sdf-v)
                   ))))))

(defun rectangle-sdf (position size)
  (lambda (pos)
    (let* ((pos (magicl:from-list (list
                                   (magicl:tref pos 0 0)
                                   (magicl:tref pos 1 0)
                                   ) '(2 1) :type 'double-float))
           (position (magicl:from-list position '(2 1) :type 'double-float))
           (dist-vec (magicl:.- (magicl:map! #'abs (magicl:.- pos position))
                                (magicl:from-list size '(2 1) :type 'double-float))))

      (+ (sqrt (magicl::sum
                (magicl:map! (lambda (x) (* x x))
                             (magicl:map! (lambda (x) (max 0d0 x)) dist-vec))))
         (min (max (magicl:tref dist-vec 0 0)
                   (magicl:tref dist-vec 1 0)
                   ) 0d0)))))

(defun ellipse-sdf (position x-l y-l)
  (let ((aspect (/ x-l y-l)))
    (lambda (pos)
      (let* ((position (cl-mpm/utils:vector-from-list (append
                                                       position '(0d0))))
             (dist-vec (cl-mpm/fastmath::fast-.*
                        (magicl:.- position pos)
                        (magicl:from-list (list 1d0 aspect 1d0) '(3 1)
                                                                             :type 'double-float)))
             (distance (sqrt (magicl:tref (magicl:@ (magicl:transpose dist-vec)
                                                    dist-vec) 0 0))))
        (- distance x-l)))))
(defun circle-sdf (position radius)
  (lambda (pos)
    (let* ((position (cl-mpm/utils:vector-from-list (append
                                                     position '(0d0))))
           (dist-vec (magicl:.- position pos))
           (distance (sqrt (magicl:tref (magicl:@ (magicl:transpose dist-vec)
                                                  dist-vec) 0 0))))
      (- distance radius)))
  (ellipse-sdf position radius radius))

(defun line-sdf (position a b width)
  (let* ((start (vector-from-list a))
         (end (vector-from-list b))
         (pa (magicl:.- position start))
         (ba (magicl:.- end start))
         (h (min 1d0 (max 0d0 (/ (cl-mpm/fastmath::dot pa ba)
                                 (cl-mpm/fastmath::dot ba ba)
                                 ))))
         (v (magicl:.- pa (magicl:scale ba h))))
    (- (sqrt (cl-mpm/fastmath::dot v v)) width)))
(defun plane-sdf (position normal distance)
  (- distance (cl-mpm/fastmath::dot position (cl-mpm/fastmath::norm normal))))
(defun plane-point-sdf (position normal point)
  (let ((distance (cl-mpm/fastmath::dot point (cl-mpm/fastmath::norm normal))))
    (- distance (cl-mpm/fastmath::dot position (cl-mpm/fastmath::norm normal)))))



(defun make-block-mps-sloped-list (offset size mps density constructor &rest args &key (slope 0) &allow-other-keys)
  "Construct a block of mxn (mps) material points real size (size) with a density (density)"
  (if (= (length size) 2)
      (setf mps (append mps '(1))
            size (append size '(0))
            offset (append offset '(0))
            ))
  (let*  ((nD 3)
          (args (alexandria:remove-from-plist args :slope))
          (spacing-0 (mapcar #'/ size mps))
          (offset (mapcar (lambda (x) (* x 0.5d0))
                          (mapcar #'+ offset spacing-0)))
          (data (loop for x from 0 to (- (first mps) 1)
                      append
                      (loop
                        for y from 0 to (- (second mps) 1)
                        append 
                        (loop
                          for z from 0 to (- (nth 2 mps) 1)
                          collect
                          (let* ((i (+ y (* x (first mps))))
                                 (spacing (list (first spacing-0)
                                                (/ (+ (second size) (* slope (+ x 1) (first spacing-0)))
                                                   (second mps))
                                                (third spacing-0)
                                                ))

                                 (volume (reduce #'* (remove-if (lambda (x) (= 0d0 x)) spacing)))
                                 (position-vec (magicl:from-list (list (+ (first offset) (* (first spacing) x))
                                                                       (+ (second offset) (* (second spacing) y))
                                                                       (+ (nth 2 offset) (* (nth 2 spacing) z))
                                                                       )
                                                                 '(3 1) :type 'double-float))
                                 (size-vec (magicl:from-list spacing '(3 1) :type 'double-float))
                                 )
                            (flet ((lisp-list (m) (loop for i from 0 to 2 collect (magicl:tref m i 0))))
                              (apply #'cl-mpm::make-particle
                                     (append (list 2 constructor)
                                             args
                                             (list
                                              :position (lisp-list position-vec)
                                              :mass (* density volume)
                                              :volume volume
                                              :size size-vec
                                              :size-0 size-vec
                                              ))))
                            ))))))
    data))

(defun estimate-elastic-dt-mp (sim p-modulus density)
  "Calculate the estimated critical timestep for an elastic mp in a mesh"
  (declare (double-float p-modulus density))
  (*
   (the double-float (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
   (the double-float (sqrt (the double-float (cl-mpm::sim-mass-scale sim))))
   (the double-float (sqrt (the double-float (/ density p-modulus))))))

(defgeneric %estimate-elastic-dt (sim))
(defmethod %estimate-elastic-dt ((sim cl-mpm:mpm-sim))
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (if (> (length mps) 0)
        (loop for mp across mps
              minimize
              (estimate-elastic-dt-mp
               sim
               (cl-mpm/particle::mp-p-modulus mp)
               (/ (the double-float (cl-mpm/particle:mp-mass mp))
                  (the double-float (cl-mpm/particle:mp-volume mp)))))
        sb-ext:double-float-positive-infinity)))

(defun estimate-elastic-dt (sim &key (dt-scale 1d0))
  (* dt-scale
     (%estimate-elastic-dt sim)))

(defun estimate-critical-damping-mp (sim E density)
  (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
        (nd (cl-mpm/mesh:mesh-nd (cl-mpm:sim-mesh sim))))
    (*
     (/ pi 2)
     (sqrt (/ E (* (expt h nd)
                   density))))))

(defgeneric estimate-critical-damping (sim))
(defmethod estimate-critical-damping ((sim cl-mpm:mpm-sim))
  (with-accessors ((mps cl-mpm:sim-mps))
      sim
    (if (> (length mps) 0)
        (loop for mp across mps
              minimize
              (estimate-critical-damping-mp
               sim
               (cl-mpm/particle::mp-e mp)
               (/ (cl-mpm/particle:mp-mass mp)
                  (cl-mpm/particle:mp-volume mp))))
        sb-ext:double-float-positive-infinity)))
