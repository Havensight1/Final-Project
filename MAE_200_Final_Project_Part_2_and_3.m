%%Step 2 MAE 200 Final Project
%start by reading in the csv files 'x_optfin' and 'u_optfin' to a variable
%x_opt and u_opt respectively.
n_steps=size(x_opt,2);
for i=1:n_steps %create A(t), and E(t) from part 1 output. expand B to duplicates of size B(t) to simplify later calcs
    x=x_opt(:,i);
    A_opt(:,:,i)=RHS_NR(x,s);
    E_opt(:,:,i)=Compute_E(x,s);
    B_opt(:,:,i)=s.B;

end
X_T=100*diag([1 1 1 1 1 1]);            % set weights for K
R=50; Q_ric=diag([1 1 1 1 1 1]);
X_ric=DRE(A_opt,B_opt,E_opt,X_T,s.T,R,Q_ric,1);

P_0=eye(6);         %set weights for L
Q_2_ric=eye(3); Q_1_ric=diag([0 0 0 0 0 0]);
P_ric=DRE(A_opt,s.C,E_ones,P_0,s.T,Q_2_ric,Q_1_ric,0);
K_opt=zeros(size(s.B,2),size(s.B,1),size(X_ric,3));
L_opt=zeros(size(s.C,2),size(s.C,1),size(P_ric,3));
for i=1:size(X_ric,3)

    K_opt(:,:,i)=-inv(R)*transpose(B_opt(:,:,i))*X_ric(:,:,i)*E_opt(:, :, i);
    L_opt(:,:,i)=-P_ric(:,:,i)*transpose(s.C)*Q_2_ric;
end
eig(E_opt(:,:,end)\A_opt(:,:,end)+E_opt(:,:,end)\B_opt(:,:,end)*K_opt(:,:,end))
eig(E_opt(:,:,end)\A_opt(:,:,end)+L_opt(:,:,end)*s.C)
swingup_sim(A_opt,B_opt,C_opt,E_opt,K_opt,L_opt,u_opt,x_0,x_tilda_0,T,h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=RHS_NR(x,s);  g=9.8;
a42=s.m1*s.ell1*(x(8)*sin(x(2))+x(5)^2*cos(x(2))); a45=2*s.m1*s.ell1*x(5)*sin(x(2));
a43=s.m2*s.ell2*(x(9)*sin(x(3))+x(6)^2*cos(x(3))); a46=2*s.m2*s.ell2*x(6)*sin(x(3));
a52=s.m1*s.ell1*(g*cos(x(2))-x(7)*sin(x(2))); a63=s.m2*s.ell2*(g*cos(x(3))-x(7)*sin(x(3)));
A=[zeros(3) eye(3); 0 -a42 -a43 0 -a45 -a46; 0 a52 0 0 0 0; 0 0 a63 0 0 0];
R=A;
end % function RHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E=Compute_E(x,s); I=eye(3); Z=zeros(3);
E=[I Z; Z [s.mc+s.m1+s.m2         -s.m1*s.ell1*cos(x(2)) -s.m2*s.ell2*cos(x(3));
           -s.m1*s.ell1*cos(x(2))  s.I1+s.m1*s.ell1^2             0            ;
           -s.m2*s.ell2*cos(x(3))          0              s.I2+s.m2*s.ell2^2   ]];
end % function Compute_E