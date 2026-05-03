;;;; Standard instance access macro
;;;

#-(or SBCL)
(eval-when (:compile-toplevel :load-toplevel :execute)
  (error "Not supported in ~A" (lisp-implementation-type)))

#+org.tfeb.tools.require-module
(org.tfeb.tools.require-module:needs
 ((:org.tfeb.hax.utilities)
  :compile t))

;; Copyright 1989-2026 Tim Bradshaw

;; Permission is hereby granted, free of charge, to any person obtaining
;; a copy of this software and associated documentation files (the
;; "Software"), to deal in the Software without restriction, including
;; without limitation the rights to use, copy, modify, merge, publish,
;; distribute, sublicense, and/or sell copies of the Software, and to
;; permit persons to whom the Software is furnished to do so, subject to
;; the following conditions:

;; The above copyright notice and this permission notice shall be
;; included in all copies or substantial portions of the Software.

;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
;; EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
;; MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
;; NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
;; LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
;; OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
;; WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

(defpackage :org.tfeb.sbcl.sia-slots
  (:use :cl)
  (:use :org.tfeb.hax.utilities)
  #+SBCL
  (:use :sb-mop)
  (:export
   #:with-sia-slots))

(in-package :org.tfeb.sbcl.sia-slots)

(provide :org.tfeb.sbcl.sia-slots)

(defstruct csi
  slot
  macro
  lvar
  svar)

(defmacro with-sia-slots ((&rest slot-specs) instance &body decls/forms)
  "Like WITH-SLOTS but used STANDARD-INSTANCE-ACCESS

This means it has all the constraints implied by that.  It probably
is worth it if you are accessing a slot many times."
  (let ((csis (mapcar (lambda (ss)
                        (etypecase ss
                          (symbol
                           (make-csi :slot ss :macro ss
                                     :lvar (symbolify nil ss "-LOCATION")
                                     :svar (symbolify nil ss "-SEEN")))
                          (cons
                           (destructuring-bind (m s) ss
                             (make-csi :slot s :macro m
                                       :lvar (symbolify nil s "-LOCATION")
                                       :svar (symbolify nil s "-SEEN"))))))
                      slot-specs)))
    (with-names (<instance>)
      `(let ((,<instance> ,instance)
             ,@(mapcan (lambda (csi)
                         `((,(csi-lvar csi) 0)
                           (,(csi-svar csi) nil)))
                       csis))
         (dolist (sd (class-slots (class-of ,<instance>)))
           (case (slot-definition-name sd)
             ,@(mapcar (lambda (csi)
                         `(,(csi-slot csi)
                           (setf ,(csi-lvar csi) (slot-definition-location sd)
                                 ,(csi-svar csi) t)))
                csis)))
         (unless (and ,@(mapcar #'csi-svar csis))
           (error "Missed slots"))
         (symbol-macrolet ,(mapcar (lambda (csi)
                                     `(,(csi-macro csi)
                                       (standard-instance-access
                                        ,<instance> ,(csi-lvar csi))))
                            csis)
           ,@decls/forms)))))

#||
(org.tfeb.tools.require-module:needs
 (:org.tfeb.tools.timing :compile t :use t))

(defclass a ()
  ((i :initform 0)))

(defmethod inc-i/sv ((a a) n)
  (declare (type fixnum n)
           (optimize speed))
  (with-slots (i) a
    (dotimes (_ n)
      (incf (the fixnum i)))
    a))

(defmethod inc-i/sia ((a a) n)
  (declare (type fixnum n)
           (optimize speed))
  (with-sia-slots (i) a
    (dotimes (_ n)
      (incf (the fixnum i)))
    a))

(defun time-sv (&key (n 10000000))
  (let ((a (make-instance 'a)))
    (timing (:tries 10 :divider (* n 2))
      (inc-i/sv a n))))

(defun time-sia (&key (n 10000000))
  (let ((a (make-instance 'a)))
    (timing (:tries 10 :divider (* n 2))
      (inc-i/sia a n))))

(defun estimate-sia-speedup (&key (n 100000000))
  (/ (time-sv :n n)
     (time-sia :n n)))

||#
