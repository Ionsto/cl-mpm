(defpackage :cl-mpm/implicit
  (:use :cl
   :cl-mpm
   :cl-mpm/particle
   :cl-mpm/mesh
   :cl-mpm/utils
   :cl-mpm/fastmaths
   )
  (:import-from
   :cl-mpm/utils varef)
  ;; (:export
  ;;  #:mpm-sim-agg-usf)
  )
(in-package :cl-mpm/implicit)

(defclass mpm-sim-implicit (cl-mpm/dynamic-relaxation::mpm-sim-dr-ul)
  ((sim-nodes-fd
    :initform nil
    :accessor sim-nodes-fd)
   (global-k
    :initform nil
    :accessor sim-global-k))
  (:documentation "NR implicit algorithm"))

(defmacro project-global-vec (sim vector accessor)
  `(progn
     (let* ((active-nodes (sim-nodes-fd ,sim))
            (proj-val ,vector)
            (nd (cl-mpm/mesh::mesh-nd (cl-mpm:sim-mesh ,sim))))
       (cl-mpm::iterate-over-nodes-array
        active-nodes
        (lambda (n)
          (loop for d from 0 below nd
                do (setf (cl-mpm/utils:varef (funcall ,accessor n) d)
                         (cl-mpm/utils:varef proj-val (+ (* nd (cl-mpm/mesh::node-stiffness-fd n)) d)))))))))

(defmacro increment-global-vec (sim vector accessor)
  `(progn
     (let* ((active-nodes (sim-nodes-fd ,sim))
            (proj-val ,vector)
            (nd (cl-mpm/mesh::mesh-nd (cl-mpm:sim-mesh ,sim))))
       (cl-mpm::iterate-over-nodes-array
        active-nodes
        (lambda (n)
          (loop for d from 0 below nd
                do (incf (cl-mpm/utils:varef (funcall ,accessor n) d)
                         (cl-mpm/utils:varef proj-val (+ (* nd (cl-mpm/mesh::node-stiffness-fd n)) d)))))))))

(defun assemble-global-vec (sim accessor)
  (declare (function accessor))
  (let* ((active-nodes (sim-nodes-fd sim))
         (nd (cl-mpm/mesh::mesh-nd (cl-mpm:sim-mesh sim)))
         (ndof (* nd (length active-nodes)))
         (v (cl-mpm/utils::arb-matrix ndof 1)))
    (cl-mpm::iterate-over-nodes-array
     active-nodes
     (lambda (n)
       (loop for d from 0 below nd
             do (setf (cl-mpm/utils:varef v (+ (* nd (cl-mpm/mesh::node-stiffness-fd n)) d))
                      (cl-mpm/utils:varef (funcall accessor n) d)))))
    (values v)))


;;Matlab version to translate
;; t = [1 2 3 4 4 5 5 6 6];                                                    % 6 to 9 component steering vector
;; J = det(F);                                                                 % volume ratio
;; [bV,bP] = eig(B); bP = [bP(1); bP(5); bP(9)];                               % eigen values/vector of the trial elastic strain tensor
;; L = parDerGen(B,bV,bP,log(bP),1./bP);                                       % derivative of the logarithmic strain
;; S = [xx   0    0    xy   0    0    0    0    xz;                          % matrix form of sigma_{il}delta_{jk}
;;      0    yy   0    0    xy   yz   0    0    0   ;
;;      0    0    zz   0    0    0    yz   xz   0   ;
;;      0    xy   0    0    xx   xz   0    0    0   ;
;;      xy   0    0    yy   0    0    0    0    yz;
;;      0    0    yz   0    0    0    yy   xy   0   ;
;;      0    yz   0    0    xz   zz   0    0    0   ;
;;      xz   0    0    yz   0    0    0    0    zz;
;;      0    0    xz   0    0    0    xy   xx   0   ];

;; T = [2*bxx  0      0      2*byx 0      0      2*bzx  0      0
;;      0      2*byy  0      0     2*bxy  2*bzy  0      0      0
;;      0      0      2*bzz  0     0      0      2*byz  2*bxz  0
;;      bxy    byx    0      byy   bxx    bzx    0      0      bzy
;;      bxy    byx    0      byy   bxx    bzx    0      0      bzy
;;      0      byz    bzy    0     bxz    bzz    byy    bxy    0
;;      0      byz    bzy    0     bxz    bzz    byy    bxy    0
;;      bxz    0      bzx    byz   0      0      byx    bxx    bzz
;;      bxz    0      bzx    byz   0      0      byx    bxx    bzz];
;; A = D(t,t)*L(t,t)*T/(2*J)-S;                                                % consistent tangent stiffness matrix


(defun 6x6-stress-to-9x9 (s)
  ;;;Original will mapping?
  ;; xx yy zz xy yx zy yz zx xz????
  (let ((mapping-list (list 0 1 2 3 3 4 4 5 5))
        ;;Potentially corrected listing?
        ;; (mapping-list (list 0 1 2 5 5 3 3 4 4))
        (result (cl-mpm/utils::arb-matrix 9 9)))
    (loop for i from 0
          for im in mapping-list
          do
      (loop for j from 0
            for jm in mapping-list
            do (setf (magicl:tref result i j) (magicl:tref s im jm))))
    result))

(defun form-G-will-form (gradients)

  )

(defun test-ul ()
  (let ((F (cl-mpm/utils::matrix-diag (list 1d0 1d0 1d0)))
        (D (cl-mpm/constitutive::linear-elastic-matrix 2d0 0d0))
        (s (cl-mpm/utils::voigt-from-list (list 0d0 0d0 0d0 0d0 0d0 0d0)))
        (b (cl-mpm/utils::matrix-from-list (list 2d0 1d0 0d0
                                                 1d0 1d0 0d0
                                                 0d0 0d0 1d0
                                                 )))
        )
    (pprint (form-ul-stiffness F D s b))
    )
  )
(defun swizzle-stiffness (d)
  ;;[TODO fix this guy]
  d)

(defun form-ul-stiffness (F D s b)
  (let ((J (magicl:det F))
        (L (cl-mpm/implicit::tensor-2nd-partial-deriv b #'log (lambda (x) (/ 1d0 x)))))
    (let* ((xx (varef s 0))
           (yy (varef s 1))
           (zz (varef s 2))
           (zx (varef s 3))
           (zy (varef s 4))
           (xy (varef s 5))
           (xz zx)
           (yz zy)
           (yx xy)
           (s-tensor
             (magicl:transpose
              (cl-mpm/utils::arb-matrix-from-list
               (list xx  0d0 0d0 xy  0d0 0d0 0d0 0d0 zx
                     0d0 yy  0d0 0d0 xy  zy  0d0 0d0 0d0
                     0d0 0d0 zz  0d0 0d0 0d0 zy  zx  0d0
                     0d0 xy  0d0 0d0 xx  xz  0d0 0d0 0d0
                     xy  0d0 0d0 yy  0d0 0d0 0d0 0d0 yz
                     0d0 0d0 yz  0d0 0d0 0d0 yy  xy  0d0
                     0d0 yz  0d0 0d0 xz  zz  0d0 0d0 0d0
                     xz  0d0 0d0 yz  0d0 0d0 0d0 0d0 zz
                     0d0 0d0 xz  0d0 0d0 0d0 xy  xx  0d0)
               9 9))))
      ;; (format t "L ~A~%" L)
      (let* ((bxx (cl-mpm/utils:mtref b 0 0))
             (bxy (cl-mpm/utils:mtref b 1 0))
             (bxz (cl-mpm/utils:mtref b 2 0))
             (byx (cl-mpm/utils:mtref b 0 1))
             (byy (cl-mpm/utils:mtref b 1 1))
             (byz (cl-mpm/utils:mtref b 2 1))
             (bzx (cl-mpm/utils:mtref b 0 2))
             (bzy (cl-mpm/utils:mtref b 1 2))
             (bzz (cl-mpm/utils:mtref b 2 2))
             (t-tensor
               (magicl:transpose
                (cl-mpm/utils::arb-matrix-from-list
                 (list
                  (* 2 bxx)  0d0        0d0        (* 2 byx) 0d0        0d0        (* 2 bzx)  0d0        0d0
                  0d0        (* 2 byy)  0d0        0d0       (* 2 bxy)  (* 2 bzy)  0d0        0d0        0d0
                  0d0        0d0        (* 2 bzz)  0d0       0d0        0d0        (* 2 byz)  (* 2 bxz)  0d0
                  bxy        byx        0d0        byy       bxx        bzx        0d0        0d0        bzy
                  bxy        byx        0d0        byy       bxx        bzx        0d0        0d0        bzy
                  0d0        byz        bzy        0d0       bxz        bzz        byy        bxy        0d0
                  0d0        byz        bzy        0d0       bxz        bzz        byy        bxy        0d0
                  bxz        0d0        bzx        byz       0d0        0d0        byx        bxx        bzz
                  bxz        0d0        bzx        byz       0d0        0d0        byx        bxx        bzz                       ) 9 9))))
        ;; (pprint (6x6-stress-to-9x9 D))
        ;; (pprint (6x6-stress-to-9x9 L))
        ;; (pprint t-tensor)
        ;; (pprint s-tensor)
        (cl-mpm/fastmaths:fast-.+
         (magicl:@
          (6x6-stress-to-9x9 (swizzle-stiffness D))
          (6x6-stress-to-9x9 L)
          (cl-mpm/fastmaths:fast-scale
           t-tensor
           (/ 1d0 (* J 2d0))))
         (cl-mpm/fastmaths:fast-scale! s-tensor -1d0)
         )))))

(defun upsize-agg ()
  )

(defun test ()
  (setup-2d :mps 3)
  ;; (change-class *sim* 'cl-mpm/dynamic-relaxation::mpm-sim-quasi-static)
  (let ((step (list ))
        (res (list )))
    (vgplot:close-all-plots)
    (setf (cl-mpm::sim-stats-oobf *sim*) 1d0)
    (loop for i from 0 below 50
          while (>= (cl-mpm::sim-stats-oobf *sim*) 1d-9)
          do
             (progn
               (cl-mpm:update-sim *sim*)
               (push i step)
               (push (cl-mpm::sim-stats-oobf *sim*) res)
               (vgplot:loglog step res)
               ;; (cl-mpm/plotter:simple-plot *sim* :plot :deformed :trial t)
               (pprint (cl-mpm::sim-stats-oobf *sim*))
               (cl-mpm/dynamic-relaxation::save-vtks *sim* "./output/" i)
               ))))

(defun test-partial ()
  (pprint (tensor-2nd-partial-deriv
           (cl-mpm/utils::matrix-from-list (list 1d0 0d0 0.1d0
                                                 0d0 1d0 0.1d0
                                                 0.1d0 0.1d0 1d0))
           #'log
           (lambda (x) (/ 1d0 x)))))
(declaim (notinline tensor-2nd-partial-deriv))
(defun tensor-2nd-partial-deriv (X func func-deriv)
  ;;; Matlab version we are translating
  ;; tol=1e-9;
  ;; Is=diag([1 1 1 0.5 0.5 0.5]); 
  ;; if (abs(eP(1))<tol && abs(eP(2))<tol && abs(eP(3))<tol)                     % all zero eigenvalues case
  ;;     L = Is;
  ;; elseif abs(eP(1)-eP(2))<tol && abs(eP(1)-eP(3))<tol                         % equal eigenvalues case
  ;;     L = ydash(1)*Is;
  ;; elseif abs(eP(1)-eP(2))<tol || abs(eP(2)-eP(3))<tol || abs(eP(1)-eP(3))<tol % repeated eigenvalues case
  ;;     if     abs(eP(1)-eP(2))<tol
  ;;         xa  = eP(3);    xc  = eP(1);
  ;;         ya  = yP(3);    yc  = yP(1);
  ;;         yda = ydash(3); ydc = ydash(1);
  ;;     elseif abs(eP(2)-eP(3))<tol
  ;;         xa  = eP(1);    xc  = eP(2);
  ;;         ya  = yP(1);    yc  = yP(2);
  ;;         yda = ydash(1); ydc = ydash(2);
  ;;     else
  ;;         xa  = eP(2);    xc  = eP(1);
  ;;         ya  = yP(2);    yc  = yP(1);
  ;;         yda = ydash(2); ydc = ydash(1);
  ;;     end
  ;;     x  = X([1 5 9 4 6 3].');
  ;;     s1 = (ya-yc)/(xa-xc)^2-ydc/(xa-xc);
  ;;     s2 = 2*xc*(ya-yc)/(xa-xc)^2-(xa+xc)/(xa-xc)*ydc;
  ;;     s3 = 2*(ya-yc)/(xa-xc)^3-(yda+ydc)/(xa-xc)^2;
  ;;     s4 = xc*s3;
  ;;     s5 = xc^2*s3;
  ;;     dX2dX=[2*xx 0      0      xy         0             xz   ;
  ;;     0      2*yy 0      xy         yz          0             ;
  ;;     0      0      2*zz 0            yz          xz          ;
  ;;     xy   xy   0     (xx+yy)/2 xz/2        yz/2        ;
  ;;     0      yz   yz   xz/2       (yy+zz)/2 xy/2        ;
  ;;     xz   0      X(3)   yz/2       xy/2        (xx+zz)/2];
  ;;     bm1  = [1 1 1 0 0 0].';
  ;;     bm11 = [1 1 1 0 0 0 ;
  ;;             1 1 1 0 0 0 ;
  ;;             1 1 1 0 0 0 ;
  ;;             0 0 0 0 0 0 ;
  ;;             0 0 0 0 0 0 ;
  ;;             0 0 0 0 0 0 ];
  ;;     L = s1*dX2dX-s2*Is-s3*(x*x.')+s4*(x*bm1.'+bm1*x.')-s5*bm11;
  ;; else                                                                        % general case (no repeated eigenvalues)
  ;;     D=[(eP(1)-eP(2))*(eP(1)-eP(3));
  ;;        (eP(2)-eP(1))*(eP(2)-eP(3));
  ;;        (eP(3)-eP(1))*(eP(3)-eP(2))];
  ;;     alfa=0; bta=0; gama=zeros(3,1); eDir=zeros(6,3);
  ;;     for i=1:3
  ;;         alfa = alfa+yP(i)*eP(i)/D(i);
  ;;         bta  = bta+yP(i)/D(i)*det(X);
  ;;         for j=1:3
  ;;             gama(i) = gama(i)+yP(j)*eP(j)/D(j)*(det(X)/eP(j)-eP(i)^2)*1/eP(i)^2;
  ;;         end
  ;;         esq = eV(:,i)*eV(:,i).';
  ;;         eDir(:,i) = [esq(1,1) esq(2,2) esq(3,3) esq(1,2) esq(2,3) esq(3,1)].';
  ;;     end
  ;;     y = inv(X);
  ;;     Ib=[y(1)^2    y(2)^2    y(7)^2     y(1)*y(2)               y(2)*y(7)               y(1)*y(7)              ;
  ;;         y(2)^2    y(5)^2    y(6)^2     y(5)*y(2)               y(5)*y(6)               y(2)*y(6)              ;
  ;;         y(7)^2    y(6)^2    y(9)^2     y(6)*y(7)               y(9)*y(6)               y(9)*y(7)              ;
  ;;         y(1)*y(2) y(5)*y(2) y(6)*y(7) (y(1)*y(5)+y(2)^2)/2    (y(2)*y(6)+y(5)*y(7))/2 (y(1)*y(6)+y(2)*y(7))/2 ;
  ;;         y(2)*y(7) y(5)*y(6) y(9)*y(6) (y(2)*y(6)+y(5)*y(7))/2 (y(9)*y(5)+y(6)^2)/2    (y(9)*y(2)+y(6)*y(7))/2 ;
  ;;         y(1)*y(7) y(2)*y(6) y(9)*y(7) (y(1)*y(6)+y(2)*y(7))/2 (y(9)*y(2)+y(6)*y(7))/2 (y(9)*y(1)+y(7)^2)/2   ];
  ;;     L = alfa*Is-bta*Ib;
  ;;     for i=1:3
  ;;         L = L+(ydash(i)+gama(i))*eDir(:,i)*eDir(:,i).';
  ;;     end
  ;; end
  ;; end
  ;;;
  (multiple-value-bind (eP eV) (cl-mpm/utils:eig X)
    (let* ((yP (mapcar func eP))
           (ydashP (mapcar func-deriv eP))
           (tol 1d-9))
      (labels ((approx-equal (x y)
                 (< (abs (- x y)) tol))
               (approx-zero (x)
                 (< (abs x) tol)))
        (let ((Is (cl-mpm/utils::tensor-voigt-4th-from-diag (list 1d0 1d0 1d0 0.5d0 0.5d0 0.5d0))))
          (destructuring-bind (e0 e1 e2) eP
            (destructuring-bind (y0 y1 y2) yP
              (destructuring-bind (yd0 yd1 yd2) ydashP
                (cond
                  ((and (approx-zero y0) (approx-zero y1) (approx-zero y2))
                   Is)
                  ((and (approx-equal y0 y1) (approx-equal y0 y2))
                   (cl-mpm/fastmaths:fast-scale! Is yd0))
                  ;;Repeated eigenvalue case
                  ((or (approx-equal y0 y1) (approx-equal y0 y2) (approx-equal y1 y2))
                   (let ((xa 0d0)
                         (xc 0d0)
                         (ya 0d0)
                         (yc 0d0)
                         (yda 0d0)
                         (ydc 0d0))
                     (cond
                       ((approx-equal y0 y1)
                        (setf xa e2
                              xc e0
                              ya y2
                              yc y0
                              yda yd2
                              ydc yd0))
                       ((approx-equal y1 y2)
                        (setf xa e0
                              xc e1
                              ya y0
                              yc y1
                              yda yd0
                              ydc yd1)
                        )
                       ((approx-equal y0 y2)
                        (setf xa e1
                              xc e0
                              ya y1
                              yc y0
                              yda yd1
                              ydc yd0)))
                     ;;TODO check order of opps here
                     (let* ((x (cl-mpm/utils:matrix-to-voight X))
                            (s1 (- (/ (- ya yc) (expt (- xa xc) 2)) (/ ydc (- xa xc))))
                            (s2 (-
                                 (* 2d0 xc (/ (- ya yc) (expt (- xa xc) 2)))
                                 (* ydc (/ (+ xa xc) (- xa xc)))))
                            (s3 (- (/ (* 2 (- ya yc)) (expt (- xa xc) 3)) (/ (+ yda ydc) (expt (- xa xc) 2))))
                            (s4 (* xc s3))
                            (s5 (* (expt xc 2) s3))
                            )
                       (destructuring-bind (xx yy zz zy zx xy) (list (varef x 0) (varef x 1) (varef x 2) (varef x 3) (varef x 4) (varef x 5))

                         ;;     dX2dX=[2*xx 0      0      xy         0             xz   ;
                         ;;     0      2*yy 0      xy         yz          0             ;
                         ;;     0      0      2*zz 0            yz          xz          ;
                         ;;     xy   xy   0     (xx+yy)/2 xz/2        yz/2        ;
                         ;;     0      yz   yz   xz/2       (yy+zz)/2 xy/2        ;
                         ;;     xz   0      X(3)   yz/2       xy/2        (xx+zz)/2];

                         ;; 2*X(1) 0      0      X(2)         0             X(3)          ;
                         ;; 0      2*X(5) 0      X(2)         X(6)          0             ;
                         ;; 0      0      2*X(9) 0            X(6)          X(3)          ;
                         ;; X(2)   X(2)   0     (X(1)+X(5))/2 X(3)/2        X(6)/2        ;
                         ;; 0      X(6)   X(6)   X(3)/2       (X(5)+X(9))/2 X(2)/2        ;
                         ;; X(3)   0      X(3)   X(6)/2       X(2)/2        (X(1)+X(9))/2];                        

                         ;;New mapping is xx yy zz zx zy xy
                         (let ((dX2dX (cl-mpm/utils::tensor-voigt-4th-from-list
                                       (list
                                        (* 2d0 xx) 0d0        0d0        xy   0d0 zx
                                        0d0        (* 2d0 yy) 0d0        xy   zy  0d0
                                        0d0        0d0        (* 2d0 zz) 0d0  zx  zx
                                        xy         xy         0d0        (/ (+ xx yy) 2)  (/ zy 2) (/ zx 2) 
                                        0d0        zy         zy         (/ zy 2)  (/ (+ yy zz) 2) (/ xy 2)
                                        zx         0d0        zx         (/ zx 2)  (/ xy 2) (/ (+ xx zz) 2)
                                        )))
                               (bm1 (cl-mpm/utils:voigt-from-list (list 1d0 1d0 1d0 0d0 0d0 0d0)))
                               (bm11 (cl-mpm/utils:tensor-voigt-4th-from-list
                                     (list
                                      1d0 1d0 1d0 0d0 0d0 0d0
                                      1d0 1d0 1d0 0d0 0d0 0d0
                                      1d0 1d0 1d0 0d0 0d0 0d0
                                      0d0 0d0 0d0 0d0 0d0 0d0
                                      0d0 0d0 0d0 0d0 0d0 0d0
                                      0d0 0d0 0d0 0d0 0d0 0d0
                                      ))))
                           ;; (pprint s1)
                           ;; (pprint s2)
                           ;; (pprint s3)
                           ;; (pprint s4)
                           ;; (pprint s5)
                           ;; (format t "dxdx ~A~%" dx2dx )
                           ;; (format t "Is ~A~%" Is)
                           ;; (format t "xx ~A~%" (magicl:@ x (magicl:transpose x)))
                           ;; (pprint
                           ;;  (cl-mpm/fastmaths:fast-.+
                           ;;   (cl-mpm/fastmaths:fast-.+
                           ;;    (cl-mpm/fastmaths:fast-scale!
                           ;;     (magicl:@ x (magicl:transpose x))
                           ;;     (- s3))
                           ;;    (cl-mpm/fastmaths:fast-scale!
                           ;;     (cl-mpm/fastmaths:fast-.+
                           ;;      (magicl:@ x (magicl:transpose bm1))
                           ;;      (magicl:@ bm1 (magicl:transpose x)))
                           ;;     s4))
                           ;;   (cl-mpm/fastmaths:fast-scale!
                           ;;    bm11
                           ;;    (- s5)))
                           ;;  )
                           ;; (pprint (cl-mpm/fastmaths:fast-.+
                           ;;          (cl-mpm/fastmaths:fast-scale
                           ;;           dX2dX
                           ;;           s1)
                           ;;          (cl-mpm/fastmaths:fast-scale
                           ;;           Is
                           ;;           (- s2))))
                           (cl-mpm/fastmaths:fast-.+
                            (cl-mpm/fastmaths:fast-.+
                             (cl-mpm/fastmaths:fast-scale
                              dX2dX
                              s1)
                             (cl-mpm/fastmaths:fast-scale
                              Is
                              (- s2)))
                            (cl-mpm/fastmaths:fast-.+
                             (cl-mpm/fastmaths:fast-.+
                              (cl-mpm/fastmaths:fast-scale!
                               (magicl:@ x (magicl:transpose x))
                               (- s3))
                              (cl-mpm/fastmaths:fast-scale!
                               (cl-mpm/fastmaths:fast-.+
                                (magicl:@ x (magicl:transpose bm1))
                                (magicl:@ bm1 (magicl:transpose x)))
                               s4))
                             (cl-mpm/fastmaths:fast-scale!
                              bm11
                              (- s5))))
                           )))))
                  ;;General case
                  (t
                   (let* ((D (cl-mpm/utils:vector-from-list (list
                                                             (* (- e0 e1) (- e0 e2))
                                                             (* (- e1 e0) (- e1 e2))
                                                             (* (- e2 e0) (- e2 e1)))))
                          (edir (cl-mpm/utils::arb-matrix 6 3))
                          )
                     (let ((alfa 0d0)
                           (beta 0d0)
                           (detx (cl-mpm/fastmaths:det-3x3 X))
                           (gama (cl-mpm/utils:vector-zeros)))
                       (dotimes (i 3)
                         (incf alfa (/ (* (nth i yP) (nth i eP)) (varef D i)))
                         (incf beta (* (/ (nth i yP) (varef D i)) detx))
                         (dotimes (j 3)
                           (incf (varef gama i)
                                 (* (nth j yP)
                                    (/ (nth j eP) (varef D j))
                                    (-
                                     (/ detx (nth j eP))
                                     (expt (nth i eP) 2))
                                    (/ 1d0 (expt (nth i eP) 2)))
                                 )
                           )
                         (let ((esq (cl-mpm/constitutive::swizzle-voigt->coombs (cl-mpm/utils:matrix-to-voight
                                                                                (magicl:@
                                                                                 (magicl:vector->column-matrix (magicl:column eV i))
                                                                                 (magicl:transpose (magicl:vector->column-matrix (magicl:column eV i))))))))
                           (dotimes (iter 6)
                             (setf (magicl:tref edir iter i) (varef esq iter)))))
                       (let* ((Ib (cl-mpm/utils::tensor-voigt-4th-zeros))
                              (y (cl-mpm/utils::matrix-copy (magicl:inv X))))
                         (flet ((s (ix iy ax ay bx by)
                                  (let ((ax (- ax 1))
                                        (bx (- bx 1))
                                        (ay (- ay 1))
                                        (by (- by 1)))
                                    (incf
                                     (magicl:tref Ib ix iy)
                                     (* 0.5d0
                                        (+ (* (mtref y ax bx) (mtref y ay by))
                                           (* (mtref y ax by) (mtref y ay bx))))))))
                           (s 0 0 1 1 1 1)
                           (s 0 1 1 1 2 2)
                           (s 0 2 1 1 3 3)
                           (s 0 3 1 1 1 2)
                           (s 0 4 1 1 2 3)
                           (s 0 5 1 1 3 1)

                           (s 1 0 2 2 1 1)
                           (s 1 1 2 2 2 2)
                           (s 1 2 2 2 3 3)
                           (s 1 3 2 2 1 2)
                           (s 1 4 2 2 2 3)
                           (s 1 5 2 2 3 1)

                           (s 2 0 3 3 1 1)
                           (s 2 1 3 3 2 2)
                           (s 2 2 3 3 3 3)
                           (s 2 3 3 3 1 2)
                           (s 2 4 3 3 2 3)
                           (s 2 5 3 3 3 1)

                           (s 3 0 1 2 1 1)
                           (s 3 1 1 2 2 2)
                           (s 3 2 1 2 3 3)
                           (s 3 3 1 2 1 2)
                           (s 3 4 1 2 2 3)
                           (s 3 5 1 2 3 1)

                           (s 4 0 2 3 1 1)
                           (s 4 1 2 3 2 2)
                           (s 4 2 2 3 3 3)
                           (s 4 3 2 3 1 2)
                           (s 4 4 2 3 2 3)
                           (s 4 5 2 3 3 1)

                           (s 5 0 3 1 1 1)
                           (s 5 1 3 1 2 2)
                           (s 5 2 3 1 3 3)
                           (s 5 3 3 1 1 2)
                           (s 5 4 3 1 2 3)
                           (s 5 5 3 1 3 1)
                           )
                         (let ((L (cl-mpm/fastmaths:fast-.-
                                   (cl-mpm/fastmaths:fast-scale Is alfa)
                                   (cl-mpm/fastmaths:fast-scale Ib beta))))
                           ;; (pprint Is)
                           ;; (pprint alfa)
                           ;; (pprint Ib)
                           ;; (pprint beta)
                           ;; (pprint ydashP)
                           ;; (pprint gama)
                           ;; (pprint edir)
                           ;; (break)
                           (dotimes (i 3)
                             (cl-mpm/fastmaths:fast-.+
                              (cl-mpm/fastmaths:fast-scale 
                               (magicl:@
                                (magicl:vector->column-matrix (magicl:column edir i))
                                (magicl:transpose (magicl:vector->column-matrix (magicl:column edir i)))
                                )
                               (* (+ (nth i ydashP) (varef gama i))))
                              L L))
                           L))))))))))))))



(defun assemble-b (strain)
  (multiple-value-bind (l v) (cl-mpm/utils::eig (cl-mpm/utils:voigt-to-matrix strain))
    (magicl:@ v
              (cl-mpm/utils::matrix-from-list
               (list
                (the double-float (exp (* 2d0 (the double-float (nth 0 l))))) 0d0 0d0
                0d0 (the double-float (exp (* 2d0 (the double-float (nth 1 l))))) 0d0
                0d0 0d0 (the double-float (exp (* 2d0 (the double-float (nth 2 l)))))
                ))
              (magicl:transpose v))))

(defgeneric assemble-mp-stiffness (mesh mp))
(defmethod assemble-mp-stiffness (mesh (mp cl-mpm/particle::particle-elastic))
  (let* ((F (cl-mpm/particle::mp-deformation-gradient mp))
         (D (cl-mpm/particle::mp-elastic-matrix mp))
         (s (cl-mpm/particle::mp-stress mp))
         (b (assemble-b (cl-mpm/particle::mp-strain mp))))
    ;; (pprint F)
    ;; (pprint D)
    ;; (pprint s)
    ;; (pprint b)
    ;; (pprint (cl-mpm/implicit::tensor-2nd-partial-deriv b #'log (lambda (x) (/ 1d0 x))))
    (form-ul-stiffness F D s b)))

(defun assemble-forces (sim)
  (assemble-global-vec sim #'cl-mpm/mesh::node-force))

(defun filter-nodes (sim filter)
  (let ((nodes (cl-mpm/mesh::mesh-nodes (cl-mpm:sim-mesh sim))))
    (remove-if-not filter (make-array (array-total-size nodes) :displaced-to nodes))))

(defun build-stiffness-matrix (sim)
  (with-accessors ((nodes-fd sim-nodes-fd)
                   (mesh cl-mpm::sim-mesh))
      sim
    (setf
     nodes-fd
     (filter-nodes sim #'cl-mpm/mesh::node-active))
    (let ()
      (lparallel:pdotimes (i (length nodes-fd))
        do (progn
             (setf (cl-mpm/mesh::node-stiffness-fd (aref nodes-fd i)) i))))
    (let* ((nd (cl-mpm/mesh:mesh-nd mesh))
           (len (* nd (length nodes-fd))))
      (setf (sim-global-k sim)
            (cl-mpm/utils::arb-matrix
             len len)))))


(defun assemble-g-3d (dsvp)
  (assemble-g-3d-prealloc dsvp (cl-mpm/utils::stretch-dsvp-3d-zeros)))

(declaim (inline assemble-g-3d-prealloc)
         (ftype (function (list magicl:matrix/double-float) magicl:matrix/double-float)
                assemble-dstretch-3d-prealloc))
(defun assemble-g-3d-prealloc (dsvp result)
  (declare (list dsvp)
           (magicl:matrix/double-float result)
           (optimize (speed 3) (safety 0) (debug 0)))
  "Assemble d/di to the strain-displacement matrix"
  (destructuring-bind (dx dy dz) dsvp
    (let* ((s (magicl::matrix/double-float-storage result)))
      (declare (double-float dx dy dz))
      (setf
       (aref s (+ 0 (* 9 0))) dx ;dx/dx
       (aref s (+ 1 (* 9 1))) dy ;dy/dy
       ;; (aref s (+ 2 (* 9 2))) dz ;dz/dz
       (aref s (+ 3 (* 9 0))) dy ;Dy/dx
       (aref s (+ 4 (* 9 1))) dx ;Dx/dy
       ;; (aref s (+ 5 (* 9 1))) dz ;dz/dx
       ;; (aref s (+ 6 (* 9 2))) dy ;dx/dz
       ;; (aref s (+ 7 (* 9 1))) dx ;Dz/dy
       ;; (aref s (+ 8 (* 9 2))) dz ;Dy/dz
       )
      result)))

(defun assemble-stiffness (sim)
  (with-accessors ((mesh cl-mpm::sim-mesh)
                   (mps cl-mpm::sim-mps)
                   (global-k sim-global-k))
      sim
    (with-accessors ((nd cl-mpm/mesh::mesh-nd))
        mesh
      (cl-mpm/fastmaths:fast-zero global-k)
      (cl-mpm::iterate-over-mps-serial
       mps
       (lambda (mp)
         (let ((stiffness (assemble-mp-stiffness mesh mp))
               (df-inv (cl-mpm/particle::mp-deformation-gradient-increment-inverse mp))
               (mp-volume (cl-mpm/particle::mp-volume mp))
               )
           (cl-mpm::iterate-over-neighbours
            mesh
            mp
            (lambda (mesh mpa node svp grads fsvp fgrads)
              (with-accessors ((node-active  cl-mpm/mesh:node-active)
                               (node-lock  cl-mpm/mesh:node-lock))
                  node
                (declare (boolean node-active)
                         (sb-thread:mutex node-lock))
                (when node-active
                  (let ((g-a (assemble-g-3d
                              (cl-mpm::gradient-push-forwards-cached
                               grads
                               df-inv))))
                    (cl-mpm::iterate-over-neighbours
                     mesh
                     mp
                     (lambda (mesh mpb node-b svp-b grads-b fsvp fgrads)
                       (when (cl-mpm::node-active node-b)
                         (let ((g-b (assemble-g-3d
                                     (cl-mpm::gradient-push-forwards-cached
                                      grads-b
                                      df-inv))))
                           (let ((stiff (magicl:@ (magicl:transpose g-a) stiffness g-b)))
                             ;; (pprint stiff)
                             ;; (sb-thread:with-mutex (node-lock)
                             ;;   (sb-thread:with-mutex ((cl-mpm/mesh::node-lock node-b))))
                             (dotimes (i nd)
                               (dotimes (j nd)
                                 (incf (magicl:tref global-k
                                                    (+ (* nd (cl-mpm/mesh::node-stiffness-fd node)) i)
                                                    (+ (* nd (cl-mpm/mesh::node-stiffness-fd node-b)) j))
                                       (* mp-volume
                                          (magicl:tref stiff i j)))))))))))
                  ))
              ))))))))

(declaim (notinline linear-solve-with-bcs))
(defun linear-solve-with-bcs (ma v bcs &optional (target-vi nil))
  (let ((target-vi (if target-vi
                       target-vi
                       (cl-mpm/utils::arb-matrix (magicl:nrows v) 1)
                       ))
        (bc-map (make-array (magicl:nrows bcs) :fill-pointer 0)))
    (cl-mpm/fastmaths:fast-zero target-vi)
    (loop for i from 0
          for bc across (magicl::storage bcs)
          do (when (> bc 0d0)
               (vector-push-extend i bc-map)))
    ;;TODO have a fallback if no bcs are applied (i.e. no resizing required)
    (let ((reduced-size (length bc-map)))
      ;;When we are solving a fully fixed system - i.e. out of plane dimensions
      (when (> reduced-size 0)
        (let* ((A-r (cl-mpm/utils::arb-matrix reduced-size reduced-size))
               (v-r (cl-mpm/utils::arb-matrix reduced-size 1)))
          (lparallel:pdotimes (i reduced-size)
            (dotimes (j reduced-size)
              (setf (mtref a-r i j) (mtref ma (aref bc-map i) (aref bc-map j))))
            (setf (varef v-r i) (varef v (aref bc-map i))))
          (let ((vs (magicl::linear-solve A-r v-r)))
            (lparallel:pdotimes (i reduced-size)
              (setf (varef target-vi (aref bc-map i)) (varef vs i))))
          )))
    target-vi))

(defun assemble-global-bcs (sim)
  (let* ((active-nodes (sim-nodes-fd sim))
         (nd (cl-mpm/mesh::mesh-nd (cl-mpm:sim-mesh sim)))
         (ndof (* nd (length active-nodes)))
         (v (cl-mpm/utils::arb-matrix ndof 1)))
    (cl-mpm::iterate-over-nodes-array
     active-nodes
     (lambda (n)
       (loop for d from 0 below nd
             do
                (let ((index (+ d (* nd (cl-mpm/mesh::node-stiffness-fd n)))))
                  (let ((bcs (cl-mpm/mesh::node-bcs n)))
                    (if bcs
                        (setf (cl-mpm/utils:varef v index) (cl-mpm/utils:varef bcs d))
                        (setf (cl-mpm/utils:varef v index) 1d0))))
                )))
    (values v)))

(declaim (notinline setup))
(defun setup-2d (&key(refine 1) (mps 2))
  (let* ((h (* 1d0 refine))
         (block-size (list 10d0 10d0))
         (domain-size (list 20d0 20d0)))
    (defparameter *sim* (cl-mpm/setup:make-simple-sim
                         h
                         (mapcar (lambda (x) (* x refine)) domain-size)
                         :sim-type 'mpm-sim-implicit
                         :args-list
                         (list
                          :enable-aggregate nil
                          :gravity -10d0
                          )))
    (cl-mpm/setup::setup-bcs
     *sim*
     :right (list nil nil nil)
     )
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup::make-block-mps (list 0d0 0d0)
                                   block-size
                                   (mapcar (lambda (e) (* (/ e h) mps)) block-size)
                                   1d3
                                   'cl-mpm/particle::particle-elastic
                                   :E 1d6
                                   :nu 0.3d0
                                   )))
  ;; (setf
  ;;  (cl-mpm:sim-gravity *sim*)

  ;;  )
  )
(defun setup ()
  (let* ((elements (expt 2 6))
         (L 50d0)
         (h (/ L elements))
         )
    (defparameter *sim* (cl-mpm/setup:make-simple-sim
                         h
                         (list elements 1)
                         :sim-type 'mpm-sim-implicit
                         :args-list
                         (list
                          :enable-aggregate nil)))
    (cl-mpm/setup::setup-bcs
     *sim*
     :right (list nil nil nil))
    (cl-mpm:add-mps
     *sim*
     (cl-mpm/setup::make-block-mps (list 0d0 h)
                                   (list L h)
                                   (list (* 2d0 (/ L h)) 1)
                                   80d0
                                   'cl-mpm/particle::particle-elastic
                                   :E 1d4
                                   :nu 0.2d0
                                   :gravity-axis (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0)))))
  (setf
   (cl-mpm:sim-gravity *sim*)
   (/ -10d0 10)
   )
  )

;; (defun setup-implicit (sim)
;;   (with-slots ((mesh cl-mpm::mesh)
;;                (mps cl-mpm::mps)
;;                (bcs cl-mpm::bcs)
;;                (mass-filter cl-mpm::mass-filter)
;;                (initial-setup initial-setup))
;;       sim
;;     (cl-mpm::reset-grid mesh)
;;     (cl-mpm::reset-node-displacement sim)
;;     (cl-mpm::p2g mesh mps vel-algo)
;;     (when (> mass-filter 0d0)
;;       (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
;;     (cl-mpm::filter-cells sim)
;;     (setf
;;      (cl-mpm/aggregate::sim-agg-nodes-fd sim)
;;      (cl-mpm/aggregate::filter-nodes sim #'cl-mpm/mesh::node-active))
;;     (let ((fd 0))
;;       (loop for n across (cl-mpm/aggregate::sim-agg-nodes-fd sim)
;;             do (progn
;;                  (setf (cl-mpm/mesh::node-agg-fd n) fd)
;;                  (incf fd))))))
(defun reduce-with-bcs (v bcs)
  (let ((bc-map (make-array (magicl:nrows bcs) :fill-pointer 0)))
    ;; (cl-mpm/fastmaths:fast-zero target-vi)
    (loop for i from 0
          for bc across (magicl::storage bcs)
          do (when (> bc 0d0)
               (vector-push-extend i bc-map)))
    (let ((reduced-size (length bc-map)))
      (if (> reduced-size 0)
        (let* ((v-r (cl-mpm/utils::arb-matrix reduced-size 1)))
          (lparallel:pdotimes (i reduced-size)
            (setf (cl-mpm/utils:varef v-r i) (cl-mpm/utils:varef v (aref bc-map i))))
          ;; (let ((vs (magicl::linear-solve A-r v-r)))
          ;;   (lparallel:pdotimes (i reduced-size)
          ;;     (setf (varef target-vi (aref bc-map i)) (varef vs i))))
          v-r
          )
        (error "no degrees of freedom")))))

(defun expand-with-bcs (v bcs)
  (let* ((expanded-size (magicl:nrows bcs))
         (target (cl-mpm/utils::arb-matrix expanded-size 1))
         (bc-map (make-array (magicl:nrows bcs) :fill-pointer 0)))
    (loop for i from 0
          for bc across (magicl::storage bcs)
          do (when (> bc 0d0)
               (vector-push-extend i bc-map)))
    (lparallel:pdotimes (i (magicl:nrows v))
      (setf (cl-mpm/utils:varef target (aref bc-map i)) (cl-mpm/utils:varef v i)))
    target))

(defun linear-solve-with-bcs (ma v bcs &optional (target-vi nil))
  (let ((target-vi (if target-vi
                       target-vi
                       (cl-mpm/utils::arb-matrix (magicl:nrows v) 1)
                       ))
        (bc-map (make-array (magicl:nrows bcs) :fill-pointer 0)))
    (cl-mpm/fastmaths:fast-zero target-vi)
    (loop for i from 0
          for bc across (magicl::storage bcs)
          do (when (> bc 0d0)
               (vector-push-extend i bc-map)))
    (let ((reduced-size (length bc-map)))
      (when (> reduced-size 0)
        (let* ((A-r (cl-mpm/utils::arb-matrix reduced-size reduced-size))
               (v-r (cl-mpm/utils::arb-matrix reduced-size 1)))
          (lparallel:pdotimes (i reduced-size)
            (dotimes (j reduced-size)
              (setf (mtref a-r i j) (mtref ma (aref bc-map i) (aref bc-map j))))
            (setf (varef v-r i) (varef v (aref bc-map i))))
          (let ((vs (magicl::linear-solve A-r v-r)))
            (lparallel:pdotimes (i reduced-size)
              (setf (varef target-vi (aref bc-map i)) (varef vs i))))
          )))
    target-vi
    ))

(defun test-imp ()
  (setup :refine 0.5)
  ;; (let ((d  (assemble-global-vec *sim* #'cl-mpm/mesh::node-displacment)))
  ;;   (setf (cl-mpm/utils:varef d 0) 1d0)
  ;;   (project-global-vec *sim* d cl-mpm/mesh::node-displacment)
  ;;   )
  (let ((lstps 50))
    (loop for l from 1 to lstps
          do
             (let ((sim *sim*)
                   (iter 0))

               (setf (cl-mpm::sim-ghost-factor *sim*) (* 1d6 1d-3))
               (setup-implicit *sim*)
               (setf (cl-mpm/aggregate::sim-enable-aggregate sim) nil)
               (with-slots ((mesh cl-mpm::mesh)
                            (mps cl-mpm::mps)
                            (bcs cl-mpm::bcs)
                            (dt-loadstep cl-mpm/dynamic-relaxation::dt-loadstep)
                            (dt cl-mpm::dt)
                            (ghost-factor cl-mpm::ghost-factor)
                            (fbar cl-mpm::enable-fbar)
                            (initial-setup initial-setup))
                   *sim*
                 (cl-mpm::p2g-force-fs *sim*)
                 (cl-mpm::apply-bcs mesh bcs dt)
                 (let* ((bcs-vec (assemble-global-bcs sim))
                        (f-ext-full (assemble-global-vec sim #'cl-mpm/mesh::node-external-force))
                        (f-ext
                          (reduce-with-bcs
                           f-ext-full
                           bcs-vec)))
                   (cl-mpm/fastmaths:fast-scale! f-ext (* -1d0 (/ (float l) lstps)))
                   (cl-mpm/linear-solver::solve-richardson
                    (lambda (d)
                      (project-global-vec sim
                                          (expand-with-bcs
                                           d
                                           bcs-vec)
                                          #'cl-mpm/mesh::node-displacment)
                      (cl-mpm::apply-bcs mesh bcs dt)
                      (cl-mpm::reset-nodes-force sim)
                      (cl-mpm::update-stress mesh mps dt-loadstep fbar)
                      (cl-mpm::p2g-force-fs sim)
                      (when ghost-factor
                        (cl-mpm/ghost::apply-ghost sim ghost-factor)
                        (cl-mpm::apply-bcs mesh bcs dt))
                      (cl-mpm::apply-bcs mesh bcs dt)
                      (incf iter)
                      (cl-mpm:iterate-over-nodes
                       mesh
                       (lambda (n)
                         (when (cl-mpm/mesh::node-active n)
                           (cl-mpm/fastmaths:fast-.+ (cl-mpm/mesh::node-internal-force n)
                                                     (cl-mpm/mesh::node-ghost-force n)
                                                     (cl-mpm/mesh::node-force n)))))
                      (reduce-with-bcs
                       (assemble-global-vec sim #'cl-mpm/mesh::node-force)
                       bcs-vec))
                    f-ext
                    :tol 1d-3
                    )))
               (cl-mpm::finalise-loadstep sim)
               (cl-mpm/output:save-vtk (merge-pathnames "./output/" (format nil "sim_~5,'0d.vtk" l)) sim)
               (cl-mpm/output:save-vtk-nodes (merge-pathnames "./output/" (format nil "sim_n_~5,'0d.vtk" l)) sim)
               (plot *sim*)
               (pprint iter)))))

(defmethod cl-mpm/dynamic-relaxation::pre-step ((sim mpm-sim-implicit))
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (dt cl-mpm::dt)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (initial-setup cl-mpm/dynamic-relaxation::initial-setup)
               (enable-aggregate cl-mpm/aggregate::enable-aggregate)
               (damping cl-mpm::damping-factor)
               (ghost-factor cl-mpm::ghost-factor)
               (vel-algo cl-mpm::velocity-algorithm))
      sim

    (setf (cl-mpm/dynamic-relaxation::sim-solve-count sim) 0)
    (cl-mpm::reset-grid mesh :reset-displacement t)
    (cl-mpm::reset-node-displacement sim)
    (cl-mpm::p2g mesh mps vel-algo)
    (when (> mass-filter 0d0)
      (cl-mpm::filter-grid mesh (cl-mpm::sim-mass-filter sim)))
    (cl-mpm::filter-cells sim)
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (n)
       (setf
        (cl-mpm/mesh::node-true-mass n) (cl-mpm/mesh:node-mass n))
       (cl-mpm/fastmaths:fast-zero (cl-mpm/mesh::node-true-velocity n))))

    (cl-mpm::zero-grid-velocity (cl-mpm:sim-mesh sim))
    (cl-mpm::apply-bcs mesh bcs dt)
    (build-stiffness-matrix sim)
    (setf initial-setup t)))


(defmethod cl-mpm::update-sim ((sim mpm-sim-implicit))
  (declare (cl-mpm::mpm-sim sim))
  (with-slots ((mesh cl-mpm::mesh)
               (mps cl-mpm::mps)
               (bcs cl-mpm::bcs)
               (bcs-force cl-mpm::bcs-force)
               (dt cl-mpm::dt)
               (dt-loadstep cl-mpm/dynamic-relaxation::dt-loadstep)
               (mass-filter cl-mpm::mass-filter)
               (split cl-mpm::allow-mp-split)
               (enable-damage cl-mpm::enable-damage)
               (nonlocal-damage cl-mpm::nonlocal-damage)
               (remove-damage cl-mpm::allow-mp-damage-removal)
               (fbar cl-mpm::enable-fbar)
               (bcs-force-list cl-mpm::bcs-force-list)
               (ghost-factor cl-mpm::ghost-factor)
               (initial-setup cl-mpm/dynamic-relaxation::initial-setup)
               (enable-aggregate cl-mpm/aggregate::enable-aggregate)
               (damping cl-mpm::damping-factor)
               (damping-scale cl-mpm/dynamic-relaxation::damping-scale)
               (solve-count cl-mpm/dynamic-relaxation::solve-count)
               (mass-update-iter cl-mpm/dynamic-relaxation::mass-update-count)
               (vel-algo cl-mpm::velocity-algorithm))
      sim
    (unless initial-setup
      (cl-mpm/dynamic-relaxation::pre-step sim))
    (setf dt 1d0)
    (cl-mpm::reset-nodes-force sim)
    (cl-mpm::update-stress mesh mps dt-loadstep fbar)
    (cl-mpm::p2g-force-fs sim)
    ;; (cl-mpm::update-node-forces sim)
    (iterate-over-nodes
     mesh
     (lambda (node)
       (when (and (cl-mpm/mesh:node-active node))
         (cl-mpm::calculate-forces node 0d0 0d0 1d0))))
    (incf solve-count)
    (cl-mpm::apply-bcs mesh bcs dt)
    (assemble-stiffness sim)
    ;;Assemble the global K matrix
    (let* ((K (sim-global-k sim))
           (f (assemble-forces sim))
           ;; (f-int (assemble-global-vec sim #'cl-mpm/mesh::node-internal-force))
           ;; (f-ext (assemble-global-vec sim #'cl-mpm/mesh::node-external-force))
           ;; (f (cl-mpm/fastmaths::fast-.+ f-ext f-int))
           (bcs (assemble-global-bcs sim)))
      (let ((ddisp (cl-mpm/fastmaths:fast-scale! (linear-solve-with-bcs K f bcs) 1d0)))
        ;; (pprint f)
        ;; (pprint ddisp)
        (increment-global-vec sim ddisp #'cl-mpm/mesh::node-displacment)))
    ;;Do 1 step of linear solver
    (cl-mpm::zero-grid-velocity mesh)
    (cl-mpm::apply-bcs mesh bcs dt)

    ;; (cl-mpm::update-nodes sim)
    ;; (cl-mpm::update-filtered-cells sim)
    (cl-mpm::apply-bcs mesh bcs dt)
    (cl-mpm::update-dynamic-stats sim)
    (cl-mpm::g2p mesh mps dt damping :TRIAL)
    (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC))
  )






(defun test-nr-imp ()
  (setup :refine 0.5)
  ;; (setup-implicit *sim*)
  ;; (let ((d  (assemble-global-vec *sim* #'cl-mpm/mesh::node-displacment)))
  ;;   (setf (cl-mpm/utils:varef d 0) 1d0)
  ;;   (project-global-vec *sim* d cl-mpm/mesh::node-displacment)
  ;;   )
  (let ((sim *sim*)
        (iter 0))
    (setf (cl-mpm/aggregate::sim-enable-aggregate sim) nil)

    (let* ((bcs-vec (assemble-global-bcs sim))
           (d-0
             (reduce-with-bcs
              (assemble-global-vec sim #'cl-mpm/mesh::node-displacment)
              bcs-vec))
           (nr-crit 1d-9)
           (nr-error nr-crit))
      ;;NR iter
      (loop for niter from 0 to 10
            while (>= nr-error nr-crit)
            do
               (progn
                 (with-slots ((mesh cl-mpm::mesh)
                              (mps cl-mpm::mps)
                              (bcs cl-mpm::bcs)
                              (dt-loadstep cl-mpm/dynamic-relaxation::dt-loadstep)
                              (dt cl-mpm::dt)
                              (fbar cl-mpm::enable-fbar)
                              (initial-setup initial-setup))
                     *sim*
                   (cl-mpm::p2g-force-fs *sim*)
                   (cl-mpm::apply-bcs mesh bcs dt)
                   (let* ((f-ext
                            (reduce-with-bcs
                             (assemble-global-vec sim #'cl-mpm/mesh::node-external-force)
                             bcs-vec)))
                     (cl-mpm::reset-node-displacement sim)
                     (cl-mpm/fastmaths:fast-scale! f-ext 1d-3)
                     (cl-mpm/linear-solver::solve-conjugant-gradients
                      (lambda (d)
                        (project-global-vec sim
                                            (expand-with-bcs
                                             d
                                             bcs-vec)
                                            #'cl-mpm/mesh::node-displacment)
                        ;; (cl-mpm::reset-nodes-force sim)
                        (cl-mpm::apply-bcs mesh bcs dt)
                        ;; (cl-mpm::update-nodes sim)
                        ;; (cl-mpm::update-cells sim)
                        (cl-mpm::reset-nodes-force sim)
                        (cl-mpm::update-stress mesh mps dt-loadstep fbar)
                        (cl-mpm::p2g-force-fs sim)
                        (cl-mpm::apply-bcs mesh bcs dt)

                        (incf iter)
                        (reduce-with-bcs
                         (assemble-global-vec sim #'cl-mpm/mesh::node-internal-force)
                         bcs-vec))
                      f-ext)))
                 (let ((r
                         (reduce-with-bcs
                          (cl-mpm/fastmaths::fast-.+
                           (assemble-global-vec sim #'cl-mpm/mesh::node-external-force)
                           (assemble-global-vec sim #'cl-mpm/mesh::node-internal-force))
                          bcs-vec)))
                   (setf nr-crit (cl-mpm/fastmaths::mag-squared r)))
                 (format t "NR iter ~D - ~E ~%" niter nr-crit)
                 )))))
