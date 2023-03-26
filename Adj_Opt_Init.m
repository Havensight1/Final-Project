function [u_k,x_k]=Adj_Opt_Init(u_start,s)
%Named Example_21_1 in NR, Ref pg: 723, Algorithm 22.1.
%if ~exist('s.Q','var'), s.Q=diag([0 0 0 0 0 0]); else, s.Q=diag(s.Q); end, if ~exist('s.R', 'var'), s.R=0; else, s.R=s.R; end, if ~exist('s.QT','var'), s.QT=diag([5 35 15 .1 60 10]); else, s.QT=diag(s.QT); end%Cost funciton penalty weights
t=s.t; T=s.T; alpha=s.alpha;  u_k=zeros(s.N+1,1);
if size(u_start)==size(zeros(s.N+1,1)),u_k=u_start; end %warm start option
[m,~]=size(u_start);
if u_start~=0,u_k(1:m,1)=u_start; end   %partial start option
s.x0=transpose([0 pi pi 0 0 0]); x_k(1:6,1)=s.x0; res=0;
for k=0:150,k;
    u=u_k(1); x=s.x0;
    J=.25*s.h*(x'*s.Q*x+u'*s.R*u); c=.5; %Step 1: march/save state from t=0 ->T. compute cost
    for n=1:s.N,u=u_k(n); f1=RHS(x,u,s); f2=RHS(x+s.h*f1/2,u,s); f3=RHS(x+s.h*f2/2,u,s); f4=RHS(x+s.h*f3,u,s);
        x=x+s.h*(f1/6+(f2+f3)/3+f4/6); x_k(1:6,n+1)=x; u=u_k(n+1);
        x_k(7:9,n)=f1(4:6); if n==s.N, c=.25; end, J=J+c*s.h*(x'*s.Q*x+u'*s.R*u);
    end, f1=RHS(x,u,s); x_k(7:9,s.N+1)=f1(4:6); E=Compute_E(x,s);    J=J+.5*(x'*E'*s.QT*E*(x));
    r=s.QT*E*x; g(s.N+1,1)=s.B'*r+s.R*u_k(s.N+1);    %Steps 2&3: March adjoint from t=T->t=0, compute gradient
    for n=s.N:-1:1, xh=(x_k(:,n+1)+x_k(:,n))/2;
        f1=RHSa(r,x_k(:,n+1),s); f2=RHSa(r-s.h*f1/2,xh,s); f3=RHSa(r-s.h*f2/2,xh,s);
        f4=RHSa(r-s.h*f3,x_k(:,n+1),s); r=r-s.h*(f1/6+(f2+f3)/3+f4/6); g(n,1)=s.B'*r+s.R*u_k(n);
    end, res1=res; res=g'*g;            %Steps 4&5: update and repeat
    if (mod(k,4)==0|alpha<1e-4), p_k=-g; else, p_k=-g+p_k*res/res1; end %conjugate gradient (why do we jump in if k=4?)
%     figure(1); clf; subplot(2,1,1); plot(t,x_k(1,:),'r-',t,x_k(2,:),'b-', t,x_k(3,:),'g-');
%                     xlabel('Time from 0 to T (sec)');
%                     ylabel('Cart Position and Pendulum Angles');
%                     legend('Cart Pos','Pendulum 1 Angle', 'Pendulum 2 Angle');
%                     subplot(2,1,2); plot(t,u_k,'r--');
%                     xlabel('Time from 0 to T (sec)');
%                     ylabel('Actuator Output');
    [AA,AB,AC,JA,JB,JC]=Bracket(@Compute_J_Adj_Opt,0,alpha,J,u_k,p_k,s); %find triplet
    [alpha, J]=Brent(@Compute_J_Adj_Opt,AA,AB,AC,JA,JB,JC,1E-5,u_k,p_k,s); %refine triplet
    u_k=u_k+alpha*p_k; pause(0.01); if abs(alpha)<1E-12, break, end     %update u_k
end
s.mc=10; for n=1:s.N+1       %Compute u_k corresponding to different s.mc to give same x_k
E=Compute_E(x_k(1:6,n),s); N=Compute_N(x_k(1:6,n),0,s); u_k(n,1)=s.B'*(E*x_k(4:9,n)-N);
         end
x_k(:,s.N)
end %funciton Adj_Opt_Init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Sub Fucitons%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E=Compute_E(x,s); I=eye(3); Z=zeros(3);
E=[I Z;Z [s.mc+s.m1+s.m2 -s.m1*s.e111*cos(x(2)) -s.m2*s.e112*cos(x(3));
          -s.m1*s.e111*cos(x(2)) s.I1+s.m1*s.e111^2 0;
          -s.m2*s.e112*cos(x(3)) 0 s.I2+s.m2*s.e112^2]];
end %function Compute_E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N=Compute_N(x,u,s)
N=[x(4);x(5);x(6); -s.m1*s.e111*sin(x(2))*x(5)^2-s.m2*s.e112*sin(x(3))*x(6)^2+u;
    s.m1*9.8*s.e111*sin(x(2)); s.m2*9.8*s.e112*sin(x(3))];
end %function Compute_N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A=Compute_A(x,s); g=9.8;
a42=s.m1*s.e111*(x(8)*sin(x(2))+x(5)^2*cos(x(2))); a45=2*s.m1*s.e111*x(5)*sin(x(2));
a43=s.m2*s.e112*(x(9)*sin(x(3))+x(6)^2*cos(x(3))); a46=2*s.m2*s.e112*x(6)*sin(x(3));
a52=s.m1*s.e111*(g*cos(x(2))-x(7)*sin(x(2))); a63=s.m2*s.e112*(g*cos(x(3))-x(7)*sin(x(3)));
A=[zeros(3) eye(3); 0 -a42 -a43 0 -a45 -a46; 0 a52 0 0 0 0; 0 0 a63 0 0 0];
end %function Compute_A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J=Compute_J_Adj_Opt(u_trial,s)
x=s.x0; u=u_trial(1); J=.25*s.h*(x'*s.Q*x+u'*s.R*u); c=.5;
for n=1:s.N, u=u_trial(n); if n==s.N, c=.25; end
    f1=RHS(x,u,s); f2=RHS(x+s.h*f1/2,u,s); f3=RHS(x+s.h*f2/2,u,s); f4=RHS(x+s.h*f3,u,s);
    x=x+s.h*(f1/6+(f2+f3)/3+f4/6); J=J+c*s.h*(x'*s.Q*x+u'*s.R*u);
end, E=Compute_E(x,s); J=J+.5*(x'*E'*s.QT*E*x);
end % function Compute_J
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=RHS(x,u,s)
E=Compute_E(x,s); N=Compute_N(x,u,s); R=E\N;
end %function RHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=RHSa(r,x,s)
E=Compute_E(x,s); A=Compute_A(x,s); R=-E'\(A'*r+s.Q*x(1:6));
end %function RHSa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AA,AB,AC,JA,JB,JC]=Bracket(Compute_J,AA,AB,JA,X,P,V)
% INPUT: {AA,AB} are guesses of A near a minimum of J(A)=Compute_J(X+A*P), with JA=J(AA).
% OUTPUT: {AA,AB,AC} bracket the minimum of J(A), with values {JA,JB,JC}.
JB=Compute_J(X+AB*P,V);  if JB>JA, [AA,AB]=Swap(AA,AB); [JA,JB]=Swap(JA,JB); end
AC=AB+2*(AB-AA); JC=Compute_J(X+AC*P,V);
while (JB>JC)                                 % At this point, {AA,AB,AC} has JA>JB>JC.
  AN=AC+2.0*(AC-AB);  JN=Compute_J(X+AN*P,V);  % Compute new point AN outside of triplet
  AA=AB; AB=AC; AC=AN;  JA=JB; JB=JC; JC=JN;  % {AB,AC,AN} -> {AA,AB,AC}
end
end % function Bracket
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AB,JB]=Brent(Compute_J,AA,AB,AC,JA,JB,JC,TOL,X,P,V)
% INPUT: {AA,AB,AC} bracket a minimum of J(A)=Compute_J(X+A*P), with values {JA,JB,JC}.
% OUTPUT: AB locally minimizes J(A), with accuarcy TOL*abs(AB) and value JB.
AINC=0; AL=min(AC,AA); AR=max(AC,AA);
if (abs(AB-AA) > abs(AC-AB)); [AA,AC]=Swap(AA,AC); [JA,JC]=Swap(JA,JC); end
for ITER=1:50
  if ITER<3; AINT=2*(AR-AL); end; TOL1=TOL*abs(AB)+1E-25; TOL2=2*TOL1; FLAG=0; % Initialize
  AM=(AL+AR)/2; if (AR-AL)/2+abs(AB-AM)<TOL2; ITER; return; end         % Check convergence
  if (abs(AINT)>TOL1 | ITER<3)
    % Perform a parabolic fit based on points {AA,AB,AC} [see (15.2)]
    T=(AB-AA)*(JB-JC); D=(AB-AC)*(JB-JA); N=(AB-AC)*D-(AB-AA)*T; D=2*(T-D); 
    if D<0.; N=-N; D=-D; end; T=AINT; AINT=AINC;
    if (abs(N)<abs(D*T/2) & N>D*(AL-AB) & N<D*(AR-AB)) % AINC=N/D within reasonable range?
      AINC=N/D; AN=AB+AINC; FLAG=1;                    % Success! AINC is new increment.
      if (AN-AL<TOL2 | AR-AN<TOL2); AINC=abs(TOL1)*sign(AM-AB); end  % Fix if AN near ends
    end
  end
  % If parabolic fit unsuccessful, do golden section step based on bracket {AL,AB,AR}
  if FLAG==0; if AB>AM; AINT=AL-AB; else; AINT=AR-AB; end; AINC=0.381966*AINT; end
  if abs(AINC)>TOL1; AN=AB+AINC; else; AN=AB+abs(TOL1)*sign(AINC); end
  JN=Compute_J(X+AN*P,V);
  if JN<=JB                                   % Keep 6 (not necessarily distinct) points
    if AN>AB; AL=AB; else; AR=AB; end         % defining the interval from one iteration
    AC=AA; JC=JA; AA=AB; JA=JB; AB=AN; JB=JN; % to the next:
  else                                        % {AL,AB,AR} bracket the minimum
    if AN<AB; AL=AN; else; AR=AN; end         % AB=Lowest point, most recent if tied w/ AA
    if (JN<=JA | AA==AB)                      % AA=Second-to-lowest point.
      AC=AA; JC=JA; AA=AN; JA=JN;             % AC=Third-to-lowest point
    elseif (JN<=JC | AC==AB | AC==AA)         % AN=Newest point
      AC=AN; JC=JN;                           % Parabolic fit based on {AA,AB,AC}
    end                                       % Golden section search based on {AL,AB,AR}
  end
end
end % function Brent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_1,new_2]=Swap(old_1,old_2)
%switch the order of two variables
new_1=old_2;
new_2=old_1;
end %funciton Swap