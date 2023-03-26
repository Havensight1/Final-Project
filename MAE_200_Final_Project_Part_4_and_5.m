%% Final Project for Controls (MAE 200)
    %% Isaac Need      03/16/2023
%define system properties and constants
g=9.81;      %acceleration of gravity @ sea level in m/s^2
m_c=10;     %mass of cart in kg
m_1=1;      %mass of pendulum 1 in kg
m_2=1/2;    %mass of pendulum 2 in kg
l_1=1;      %length of pendulum 1 in m
l_2=1/2;    %length of pendulum 2 in m
I_1=m_1*(l_1)^2*1/3;   %mass moment of inertia of pendulum 1
I_2=m_2*(l_2)^2*1/3;   %mass moment of inertia of pendulum 2

%intermediate variables to simplify input of the linearized dynamics
a_f=m_c+m_1+m_2;
b_f=-m_1*l_1;
c_f=-m_2*l_2;
d_f=I_1+m_1*l_1^2;
e_f=I_2+m_2*l_2^2;
f_f=-m_1*g*l_1;
h_f=-m_2*g*l_2;

%% infinite time horizon in the form E*dx/dt=A*x+B*u and y=C*x
E_f=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 a_f b_f c_f;0 0 0 b_f d_f 0; 0 0 0 c_f 0 e_f]; %matrix E based on the dynamics
N_f=[0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1;0 0 0 0 0 0;0 -f_f 0 0 0 0; 0 0 -h_f 0 0 0]; %matrix A based on the dynamics
B_f=[0;0;0;1;0;0];    %matrix B for our actuator only impacting the velocity of the base
C_f=[1 0 0 0 0 0;0 1 0 0 0 0; 0 0 1 0 0 0];%matrix C for our our sensors only observing position and angles of the cart and pendulums repsectively

A_f=E_f\N_f;            %multiply through by the inverse of E to simplify
B_f=E_f\B_f;            %multiply through by the inverse of E to simplify
Q_f=eye(6);R_f=1;  %set Q and R to eye for starters
Q_1=eye(6); Q_2=eye(size(C_f,1));
[X_f,K_f,]= icare(A_f,B_f,Q_f,R_f);  %solve the continuous time algebraic Ricarri equation for K setting Q and R to eye and 1 respectively and letting S and G default to zero and E default to 1
[P_f,L_f,]= icare(transpose(A_f),transpose(C_f),Q_1,Q_2);  %solve the continuous time algebraic Ricarri equation for L setting Q and R to eye and 1 respectively and letting S and G default to zero and E default to 1
% %[cont_f,obs_f]=Cont_Obs_Matrix(A_f,B_f,C_f);    %check that the solution we just got stabilizes the system
% rank(cont_f);
% eig(A_f-B_f*K_f);
% rank(obs_f);
K_f=-K_f;
L_f;
L_f=P_f*C_f'*Q_2;
Plot_inf_resp(A_f,B_f,K_f,x_0,T,h);
