(defpackage :cl-mpm/output
  (:use :cl))

(in-package :cl-mpm/output)

(defun save-sim (sim)
  "Save a simulation such that it can be restarted from that point")
(defun load-sim (sim))

(defun sample-line (sim start end steps get-value)
	"Sample over a line from start to end with fixed steps count, get-value function gets given a sampling mp and a mesh and should return a scalar"
    (let (  (sample-mp (cl-mpm/particle:make-particle 2 :pos start))
            (direction (magicl:from-list (mapcar (lambda (s e) (/ (- e s) steps)) start end) '(2 1))))
      (loop :repeat (+ steps 1) collect 
              (prog1
                (funcall get-value (cl-mpm:sim-mesh sim) sample-mp)
                (setf (cl-mpm/particle:mp-position sample-mp) (magicl:.+ (cl-mpm/particle:mp-position sample-mp) direction))))))

(defun sample-point (sim point get-value)
    (let ((sample-mp (cl-mpm/particle:make-particle 2 :pos point)))
        (funcall get-value (cl-mpm:sim-mesh sim) sample-mp)))

(defmacro scalar-average (accessor)
    `(lambda (mesh mp)
        (let ((av 0))
          (progn
            (cl-mpm::iterate-over-neighbours mesh mp 
                 (lambda (mesh mp node svp dsvp) 
                   (setf av (+ av (* svp (,accessor node))))))
            av))))

(defmacro matrix-average (accessor shape)
    `(lambda (mesh mp)
        (let ((av (magicl:zeros ,shape)))
          (progn
            (cl-mpm::iterate-over-neighbours mesh mp 
             (lambda (mesh mp node svp dsvp) 
               (setf av (magicl:.+ av (magicl:scale (,accessor node) svp)))))
            av))))

(defun sample-line-mass (sim start end steps)
    (sample-line sim start end steps 
        (scalar-average cl-mpm/mesh:node-mass)))

(defun sample-line-velocity (sim start end steps)
    (sample-line sim start end steps 
        (matrix-average cl-mpm/mesh:node-velocity '(2 1))))

(defun sample-point-mass (sim point)
    (sample-point sim point 
        (scalar-average cl-mpm/mesh:node-mass)))

(defun sample-point-velocity (sim point)
    (sample-point sim point
        (matrix-average cl-mpm/mesh:node-velocity '(2 1))))
