function [epsE] = MCconst(epsTr,mCst)

%Mohr-Coulomb linear elastic perfectly plastic constitutive model
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   29/01/2019
% Description:
% Mohr-Coulomb perfect plasticity constitutive model with non-associated
% plastic flow.  The stress integration routine is based on the
% following paper:
% Clausen, J., Damkilde, L. and Andersen, L. (2007). An efficient
% return algorithm for non-associated plasticity with linear yield criteria
% in principal stress space, Comput. Struct., 85(23-24), 1795-1807.
%
%--------------------------------------------------------------------------
% [Dalg,sig,epsE] = MCCONST(epsTr,mCst)
%--------------------------------------------------------------------------
% Input(s):
% epsTr  - trial elastic strain (6,1)
% mCst   - material constants 
%--------------------------------------------------------------------------
% Ouput(s);
% sig    - Cauchy stress (6,1)
% epsE   - elastic strain (6,1)
% Dalg   - algorithmic consistent tangent (6,6)
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------

%E=mCst(1); v=mCst(2); phi=mCst(4); psi=mCst(5); c=mCst(6);                  % material constants
E = 1;
v = 0.1;
phi=0.1;                                                                   % friction angle (opening angle of yield surface)
psi=0.05;                                                                  % dilation angle (set equal to phi for associated flow)
c=0.5;                                                                        % cohesion

tol=1e-12;
Ce=(-ones(3)*v+(1+v)*eye(3))/E;                                             % principal elastic compliance matrix
De=[inv(Ce) zeros(3); zeros(3) E/(2*(1+v))*eye(3)];                         % six-component elastic stiffness matric
De3=De(1:3,1:3);                                                            % principal elastic stiffness matric
[t,epsVal]=eig([epsTr(1)   epsTr(4)/2 epsTr(6)/2; 
                epsTr(4)/2 epsTr(2)   epsTr(5)/2; 
                epsTr(6)/2 epsTr(5)/2 epsTr(3)]);                           % eigen values and vectors of trial strain
epsTr=[epsVal(1); epsVal(5); epsVal(9)];                                    % principal values of strain
[epsTr,sO]=sort(epsTr,'descend');                                           % order principal values
t=t(:,sO);                                                                  % order eigen vectors 
Q=[t(1)*t(1)   t(2)*t(2)   t(3)*t(3)   t(1)*t(2)           t(2)*t(3)           t(3)*t(1)           ;     % six component principal mapping matrix
   t(4)*t(4)   t(5)*t(5)   t(6)*t(6)   t(4)*t(5)           t(5)*t(6)           t(6)*t(4)           ;
   t(7)*t(7)   t(8)*t(8)   t(9)*t(9)   t(7)*t(8)           t(8)*t(9)           t(9)*t(7)           ;
   2*t(1)*t(4) 2*t(2)*t(5) 2*t(3)*t(6) t(1)*t(5)+t(4)*t(2) t(2)*t(6)+t(5)*t(3) t(3)*t(4)+t(6)*t(1) ;
   2*t(4)*t(7) 2*t(5)*t(8) 2*t(6)*t(9) t(4)*t(8)+t(7)*t(5) t(5)*t(9)+t(8)*t(6) t(6)*t(7)+t(9)*t(4) ;
   2*t(7)*t(1) 2*t(8)*t(2) 2*t(9)*t(3) t(7)*t(2)+t(1)*t(8) t(8)*t(3)+t(2)*t(9) t(9)*t(1)+t(3)*t(7)];
sig=Ce\epsTr;                                                               % trial Cauchy stress
epsE=epsTr;                                                                 % initial elastic strain
Dalg=De;                                                                    % algorithmic tangent (elastic behaviour)
k=(1+sin(phi))/(1-sin(phi));  sigC=2*c*sqrt(k);                             % M-C yield parameters 
f=k*sig(1)-sig(3)-sigC;                                                     % M-C yield function
if f>tol                                                                    % elasto-plastic behaviour
  m=(1+sin(psi))/(1-sin(psi));                                              % M-C plastic potential parameter  
  sigA=sigC/(k-1)*ones(3,1);                                                % apex stress
  r1 =[1 1 k].'; r2 =[1 k k].';                                             % M-C yield surface tension & compression meridians
  rg1=[1 1 m].'; rg2=[1 m m].';                                             % M-C plastic potential surface tension & compression meridians
  df=[k 0 -1].';                                                            % yield surface normal
  dg=[m 0 -1].';                                                            % flow direction
  rp=De3*dg/(dg.'*De3*df);                                                  % elasto-plastic return direction (planar return)
  t1=(rg1.'*Ce*(sig-sigA))/(rg1.'*Ce*r1);                                   % return region equations
  t2=(rg2.'*Ce*(sig-sigA))/(rg2.'*Ce*r2);
  f12=(cross(rp,r1)).'*(sig-sigA); 
  f13=(cross(rp,r2)).'*(sig-sigA);
  if t1>tol && t2>tol                                                       % apex stress return
    sig=sigA;                                                               % return stress
    Dep=zeros(6);                                                           % elasto-plastic tangent
  elseif f12<tol && f13<tol                                                 % line 1 stress return
    sig=sigA+t1*r1;
    Dep=[r1*rg1'/(r1'*Ce*rg1) zeros(3); zeros(3) E/(2*(1+v))*eye(3)];
  elseif f12>tol && f13>tol                                                 % line 2 stress return 
    sig=sigA+t2*r2;
    Dep=[r2*rg2'/(r2'*Ce*rg2) zeros(3); zeros(3) E/(2*(1+v))*eye(3)];
  else                                                                      % main plane stress return
    sig=sig-f*rp;
    Dep=De3-(De3*(dg*df.')*De3)/(df.'*De3*dg);
    Dep=[Dep zeros(3); zeros(3) E/(2*(1+v))*eye(3)];
  end
  epsE=Ce*sig;                                                              % elastic strains
  Dalg=Dep;                                                                 % algorithmic tangent
  T=zeros(3);                                                               % shear modification matrix
  sigTr=Ce\epsTr;                                                           % trial stress
  if abs((sigTr(1)-sigTr(2)))>1e-3 
      T(1,1)=(sig(1)-sig(2))/(sigTr(1)-sigTr(2));                           % modify xy term
  end
  if abs((sigTr(2)-sigTr(3)))>1e-3
      T(2,2)=(sig(2)-sig(3))/(sigTr(2)-sigTr(3));                           % modify yz term
  end
  if abs((sigTr(1)-sigTr(3)))>1e-3
      T(3,3)=(sig(1)-sig(3))/(sigTr(1)-sigTr(3));                           % modify zx term
  end
  Dalg(4:6,4:6)=T*Dalg(4:6,4:6);                                            % modify shear terms
  Dalg=Q.'*Dalg*Q;                                                          % map algorithmic tangent to generalised space
end
sig=Q.'*[sig; zeros(3,1)];                                                  % map principal stresses to generalised stresses
epsE=Q\[epsE; zeros(3,1)];                                                  % map principal elastic strains to generalised elastic strains
