
;; This buffer is for text that is not saved, and for Lisp evaluation.
;; To create a file, visit it with C-x C-f and enter text in its buffer.

(defun ql-dist::dist-name-pathname (name)
  "Return the pathname that would be used for an installed dist with
 the given NAME."
  (ql-dist::qmerge (make-pathname :directory (list* :relative "dists"
                                                    (uiop:split-string name :separator "/")))))
(defun digikar99-dist-enumeration-function ()
  "The default function used for producing a list of dist objects."
  (loop for file in (directory (ql-dist::qmerge "dists/digikar99/*/distinfo.txt"))
        collect (ql-dist::make-dist-from-file file)))
(push 'digikar99-dist-enumeration-function ql::*dist-enumeration-functions*)
(ql-dist:install-dist "http://dist.ultralisp.org/digikar99/specialized-array-dispatch.txt"
                      :prompt nil)
;;; If the install-dist step gives a "can't create directory" error, manually
;;; create the directory $QUICKLISP_HOME/dists/digikar99
(ql:update-dist "digikar99/specialized-array-dispatch")
(ql:quickload "dense-numericals") ; or "numericals"
(asdf:test-system "dense-numericals") ; or "numericals"
