%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: William Coombs
% Date:   ???
% Description:
% Non-associated flow perfect plasticity Drucker-Prager model.  The stress
% update uses a NR proceedure combined with a identified region to catch
% points that should return to the apex of the yield surface.  The stress
% state is initially shifted according to the apex hydrostatic stress to
% make the numerics simpler so therefore the yield funciton of line 60
% assumes a zero apex hydrostatic stress. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% epsTr:    trial elastic strains (6 by 1)
% E:        Young's modulus (scalar)
% v:        Poisson's ratio (scalar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
% sig:      updated Cauchy stress (6 by 1)
% epsE:     updated elastic strains (6 by 1)
% Dalg:     consistent stiffness matrix (6 by 6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [sig,epsE,Dalg] = DPconst(epsTr)

E = 1;
v = 0.1;

tol=1e-12;                                                                  % convergence tolerance for NR iterations
tolf=1e-6;                                                                  % yield function tolerance
maxit=5;                                                                    % maximum NR iterations for convergence
phi=0.1;                                                                   % friction angle (opening angle of yield surface)
psi=0;                                                                  % dilation angle (set equal to phi for associated flow)
c=1;                                                                        % cohesion
alfa=-tan(phi);                                                             % yield surface slope (xi-rho space)
bta=-tan(psi);                                                              % plastic potential slope (xi-rho space)
xsic=sqrt(3)*cot(phi)*c;                                                    % apex hydrostatic stress

Ce=(-ones(3)*v+(1+v)*eye(3))/E;                                             % principal compliance matrix
De=E/(1+v)/(1-2*v)*[ones(3)*v+eye(3)*(1-2*v) zeros(3);                      % elastic stiffness matrix 
                    zeros(3)                 (1-2*v)/2*eye(3)];
De3=De(1:3,1:3)                                                            % principal elastic stiffness matrix
[t,epsVal]=eig([epsTr(1) epsTr(4)/2 epsTr(6)/2;                             % principal strains & directions
                epsTr(4)/2 epsTr(2) epsTr(5)/2; 
                epsTr(6)/2 epsTr(5)/2 epsTr(3)]);
epsTr=[epsVal(1); epsVal(5); epsVal(9)];                                    % principal strains
[epsTr,sO]=sort(epsTr,'descend');                                           % order principal strains
t=t(:,sO);                                                                  % order principal directions
t
epsTr

Q=[t(1)*t(1)   t(2)*t(2)   t(3)*t(3)   t(1)*t(2)           t(2)*t(3)           t(3)*t(1)           ;  % principal to 6 component mapping matrix 
   t(4)*t(4)   t(5)*t(5)   t(6)*t(6)   t(4)*t(5)           t(5)*t(6)           t(6)*t(4)           ;
   t(7)*t(7)   t(8)*t(8)   t(9)*t(9)   t(7)*t(8)           t(8)*t(9)           t(9)*t(7)           ;
   2*t(1)*t(4) 2*t(2)*t(5) 2*t(3)*t(6) t(1)*t(5)+t(4)*t(2) t(2)*t(6)+t(5)*t(3) t(3)*t(4)+t(6)*t(1) ;
   2*t(4)*t(7) 2*t(5)*t(8) 2*t(6)*t(9) t(4)*t(8)+t(7)*t(5) t(5)*t(9)+t(8)*t(6) t(6)*t(7)+t(9)*t(4) ;
   2*t(7)*t(1) 2*t(8)*t(2) 2*t(9)*t(3) t(7)*t(2)+t(1)*t(8) t(8)*t(3)+t(2)*t(9) t(9)*t(1)+t(3)*t(7)];
sig=De3*epsTr-xsic/sqrt(3)                                                 % shift stress state according to apex stress
epsE=Ce*sig;                                                                % estimate of elastic strains 
epsEtr=epsE;                                                                % trial elastic strain state
xi=sum(sig)/sqrt(3)                                                        % hydrostatic stress
s=sig-xi/sqrt(3)  
rho=sqrt(sum(s.^2))                                     % deviatoric stress
f=rho-alfa*xi;                                                              % yield function
if f>tol                                                                    % if yielded
  disp("Yeilding")
  fap=rho*sqrt(1+v)+xi*sqrt(1-2*v)/(bta*sqrt(1+v)/sqrt(1-2*v));             % apex retunr surface 
  if fap<tol                                                                % apex return
    sig=xsic/sqrt(3)*ones(3,1);                                             % updated stress
    Dalg=zeros(6);                                                          % Algorithmic tangent
  else                                                                      % standard surface return
    b=zeros(4,1); b(4)=f; itnum=0; dgam=0;                                  % set initial NR parameters
    df=s/rho-alfa/sqrt(3);                                                  % derivative of the YS wrt. stress
    dg=s/rho-bta/sqrt(3);                                                   % flow direction
    ddg=1/3*(3*eye(3)-ones(3))/rho-s*s.'/rho^3;                             % derivative of the flow direction wrt. stress 
    while (itnum<maxit) && ((norm(b(1:3))>tol) || (abs(b(4))>tolf))         % NR loop 
      A=[eye(3)+dgam*ddg*De3 dg; df'*De3 0]                                % Jacobian matrix for NR iterations
      dx=-inv(A)*b 
      epsE=epsE+dx(1:3)
      dgam=dgam+dx(4);                     % update unknowns
      
      sig=De3*epsE;                                                         % update stress 
      xi=sum(sig)/sqrt(3);                                                  % hydrostatic stress 
      s=sig-xi/sqrt(3); rho=sqrt(sum(s.^2));                                % deviatoric stress
      df=s/rho-alfa/sqrt(3);                                                % derivative of the YS wrt. stress
      dg=s/rho-bta/sqrt(3);                                                 % flow direction
      ddg=1/3*(3*eye(3)-ones(3))/rho-s*s.'/rho^3;                           % derivative of the flow direction wrt. stress 
      b=[epsE-epsEtr+dgam*dg; rho-alfa*xi];                                 % residuals
      itnum=itnum+1;                                                        % increment iteration number 
    end
    B=inv([Ce+dgam*ddg dg; df.' 0]);                                        % linearisation for Dalg
    Dalg=[B(1:3,1:3) zeros(3); zeros(3) E/(2*(1+v))*eye(3)];                % consistent tangent
    sig=sig+xsic/sqrt(3);                                                   % shift stress back according to apex stress
    sigTr=De3*epsTr;                                                        % trial stress state
    T=zeros(3);                                                             % shear modification of Dalg
    if abs((sigTr(1)-sigTr(2)))>1e-3; T(1,1)=(sig(1)-sig(2))/(sigTr(1)-sigTr(2)); end
    if abs((sigTr(2)-sigTr(3)))>1e-3; T(2,2)=(sig(2)-sig(3))/(sigTr(2)-sigTr(3)); end
    if abs((sigTr(1)-sigTr(3)))>1e-3; T(3,3)=(sig(1)-sig(3))/(sigTr(1)-sigTr(3)); end
    Dalg(4:6,4:6)=T*Dalg(4:6,4:6);                                          % modify shear terms
    Dalg=Q.'*Dalg*Q;                                                        % map Dalg to 6 component stress space
  end
  sig
  epsE=Ce*sig;                                                              % updated elastic strains
else                                                                        % elastic behaviour
  disp("elastic")
  sig=De3*epsTr; epsE=epsTr; Dalg=De;                                       % stress, elastic strains and Dalg for elastic behaviour
end
sig=Q.'*[sig; zeros(3,1)];
epsE=Q\[epsE; zeros(3,1)];                       % map stress and elastic strains to 6 component stress space
