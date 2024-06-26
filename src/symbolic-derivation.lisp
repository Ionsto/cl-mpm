;
;Copyright Aleksander Ksiazek 
;Made for the course Obliczenia Symboliczne II at AGH in 2014
;Developed and tested using Emacs with SLIME and SBCL 1.1


;; to compile into a standalone executable with sbcl

;; uncomment the following: (defun main () (sb-impl::toplevel-repl nil))

;; and use (sb-ext:save-lisp-and-die "symbolic" :executable t :toplevel 'main)

(defpackage :symbolic-derivation
  (:use :cl)
  (:export #:derive))
(in-package :symbolic-derivation)

(defun flatten (list)
  (cond
    ((null list) nil)
    ((atom list) (list list))
    (T (mapcan #'flatten list))
    ))

(defun not-in (var expr)
  (not (member var (flatten expr))))

(defun rewrite (expr)
; reduce redundant operations
  (if (atom expr) expr
      (let ((operator (car expr)) (arg1 (cadr expr)) (arg2 (caddr expr)))
	(cond
	  ((eq operator '+)
	   (cond
	     ((null arg2) (rewrite arg1)) ; +(- arg1)
	     ((eq arg2 0) arg1) ; +0
	     ((eq arg1 0) arg2) ; 0+
	     ((and (numberp arg1) (numberp arg2)) (+ arg1 arg2))
	     (T expr)
	     ))
	  ((eq operator '-)
	   (cond
	     ((and (eq arg1 0) (null arg2)) 0) ; (-0)
	     ((and (null arg2) (listp arg1) (eq (car arg1) '-)) (rewrite (cadr arg1))) ; -(- arg1)
	     ((eq arg1 0) (if (eq arg2 0) 0 (list '- arg2))) ; 0 - arg2
	     ((eq arg2 0) arg1) ; arg1 - 0
	     ((and (numberp arg1) (numberp arg2)) (- arg1 arg2))
	     (T expr)
	     ))
	  ((eq operator '*)
	   (cond
	     ((or (eq arg2 0) (eq arg1 0)) 0) ; 0*
	     ((eq arg2 1) arg1) ; 1*
	     ((eq arg1 1) arg2) ; *1
	     ((and (numberp arg1) (numberp arg2)) (* arg1 arg2))
	     (T expr)
	     ))
	  ((eq operator '/) 
	   (cond
	     ((eq arg2 0) (error "divide by zero"))
	     ((eq arg1 0) 0)
	     ((and (numberp arg1) (numberp arg2)) (/ arg1 arg2))
	     (T expr)
	     ))
	  ((eq operator 'expt) 
	   (cond
	     ((eq arg2 0) 1) ; f(x)expt0
	     ((eq arg2 1) arg1) ; f(x)expt1
	     ((and (numberp arg1) (numberp arg2)) (if (eq arg1 0) 0 (exp (* arg2 (log arg1)))))
	     (T expr)
	     ))
	  (T expr)
	  ))))     


(defmacro derivative (var expr)
  `(derive (quote ,var) (quote ,expr)))

(defun is-named-function (name)
  (if (member name '(expt abs sin cos tg ctg sqrt exp ln log asin acos atg actg)) T nil))

(defun named-function-derivative (name arg)
  (cond
    ;((eq name 'abs) `(signum ,(rewrite arg)))
    ((eq name 'abs) `(if (> ,(rewrite arg) 0) 1 -1))
    ((eq name 'sin) `(cos ,(rewrite arg)))
    ((eq name 'cos) `(- (sin ,(rewrite arg))))
    ((eq name 'tg) `(+ 1 (expt (tg ,(rewrite arg)) 2)))
    ((eq name 'ctg) `(- (+ 1 (expt (ctg ,(rewrite arg)) 2))))
    ((eq name 'sqrt) `(/ 1 (* 2 (sqrt ,(rewrite arg)))))
    ((eq name 'exp) (rewrite `(exp ,(rewrite arg))))
    ((eq name 'ln) `(/ 1 ,(rewrite arg)))
    ((eq name 'log) `(/ 1 (* ,(rewrite arg) (ln 10))))
    ((eq name 'asin) `(/ 1 (sqrt (- 1 (expt ,(rewrite arg) 2)))))
    ((eq name 'acos) `(- (/ 1 (sqrt (- 1 (expt ,(rewrite arg) 2))))))
    ((eq name 'atg) `(/ 1 (+ 1 (expt ,(rewrite arg) 2))))
    ((eq name 'actg) `(- (/ 1 (+ 1 (expt ,(rewrite arg) 2)))))
    ))

(defun pow (var arg1 arg2)
  (cond
	; aexptf(x)
    ((not-in var arg1)
     (rewrite (list '* 
		    (rewrite (list '* 
				   (rewrite (list 'expt arg1 arg2))
				   `(ln ,arg1)))
		    (derive var arg2))))
	; f(x) expt a
    ((not-in var arg2)
     (rewrite (list '* 
		    (rewrite (list '* arg2
				   (rewrite (list 'expt arg1 
						  (rewrite (list '- arg2 '1))))))
		    (derive var arg1))))
	; f(x) expt g(x)
    (T (derive var (list 'exp (rewrite (list '* arg2 (list 'ln arg1))))))
       ))

; processes already quoted expression by derivation rules and named-functions
(defun derive (var expr)
  (cond
    ((null expr) nil)
    ((and (atom expr) (eq var expr)) 1)
    ((not (member var (flatten expr))) 0)
    (T (let ((operator (car expr)) (arg1 (rewrite (cadr expr))) (arg2 (rewrite (caddr expr))))
	 (cond
	   ((member operator '(+ -)) 
	    (if (null arg2) (rewrite (list operator (derive var arg1)))
		(rewrite (list operator (derive var arg1) (derive var arg2)))))
	   ((eq operator '*)
	    (rewrite (list '+
			   (rewrite (list '* (derive var arg1) arg2))
			   (rewrite (list '* arg1 (derive var arg2)))
			   )))
	   ((eq operator '/)
	    (rewrite (list '/ 
			   (rewrite (list '- 
					  (rewrite (list '* (derive var arg1) arg2))
					  (rewrite (list '* arg1 (derive var arg2)))
					  ))
			   (rewrite (list 'expt arg2 2)))))
	   
	   ((eq operator 'expt) 
	    (pow var arg1 arg2))
	   ((eq operator 'if) 
        (let ((arg3 (rewrite (cadddr expr))))
          (rewrite (list 'if arg1 (derive var arg2) (derive var arg3)))))
	   (T 
	    (if (is-named-function operator) 
		(rewrite (list '* (named-function-derivative operator arg1) (derive var arg1)))
		(error "unknown function")))
	   )))))

(defmacro create-svp (arg form)
  `(lambda (,arg) ,form))

(defmacro create-dsvp (arg form)
  (let ((dx (derive arg form)))
    `(lambda (,arg) ,dx)))

