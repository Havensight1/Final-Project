x = x_opt(1:6, 1);
x_prime=zeros(size(x_opt,1)-3,size(x_opt,2));
h=0.01;
num_step=size(x_opt,2);
time=zeros(1,num_step);
for i=1:num_step-1
    time(i+1)=time(i)+h;
end
x_sim_opt=zeros(size(x_opt,1)-3,size(x_opt,2));
x_sim_opt(:,1)=x_opt(1:6,1);
u_sim=zeros(size(x_opt,2),1);
u_sim(1)=u_opt(1);
for i=1:size(x_opt,2)-1
    %%using derrived A(t), and E(t) to simulate
    x_bar=x_opt(1:6,i); u_bar=u_opt(i); K=K_opt(:,:,i);
    
    A=A_opt(:,:,i); B=B_opt(:,:,i); E=E_opt(:,:,i);
    f1=RHS_sys(A,B,E,K,x_temp,u_bar,x_bar);
    f2=RHS_sys(A,B,E,K,x_temp+h*f1/2,u_bar,x_bar);
    f3=RHS_sys(A,B,E,K,x_temp+h*f2/2,u_bar,x_bar);
    f4=RHS_sys(A,B,E,K,x_temp+h*f3,u_bar,x_bar);
%% Muhan's way: calulate RHS based on x_bar+x_prime
%     f1=RHS_time_march(x_temp+x_bar,u_sim(i),s);
%     f2=RHS_time_march((x_temp+h*f1/2)+x_bar,u_sim(i),s);
%     f3=RHS_time_march((x_temp+h*f2/2)+x_bar,u_sim(i),s);
%     f4=RHS_time_march((x_temp+h*f3)+x_bar,u_sim(i),s);

    x_prime(:,i+1)=x_prime(:,i)+h*(f1/6+(f2+f3)/3+f4/6);
    u_sim(i)=K*x_prime(1:6,i)+u_opt(i);
    x_sim_opt(:,i+1)=x_opt(1:6,i)+h*x_prime(1:6,i);
end


%Plot various states
figure(1)
plot(time, x_opt(1, :), '-..b'); grid on; hold on;
plot(time, x_sim_opt(1, :), '-.b');

figure(1)
plot(time, x_opt(2, :), '-..r'); grid on; hold on;
plot(time, x_sim_opt(2, :), '-.r');

figure(1)
plot(time, x_opt(3, :), '-..g'); grid on; hold on;
plot(time, x_sim_opt(3, :), '-.g');

% figure(1)
% plot(time, x_opt(5, :), '-b'); grid on; hold on;
% plot(time, x_sim_opt(5, :), '-.b');
% 
% figure(6)
% plot(time, x_opt(6, :), '-b'); grid on; hold on;
% plot(time, x_sim_opt(6, :), '-.b');
% 
% figure(7)
% plot(time, u_opt, '-b'); grid on; hold on;
% plot(time, u_sim, '-.b');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=RHS_time_march(x,u,s)
E=Compute_E(x,s); N=Compute_N(x,u,s); R=E\N;
end % function RHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N=Compute_N(x,u,s)
N=[x(4); x(5); x(6); -s.m1*s.ell1*sin(x(2))*x(5)^2-s.m2*s.ell2*sin(x(3))*x(6)^2+u;
 s.m1*9.8*s.ell1*sin(x(2)); s.m2*9.8*s.ell2*sin(x(3)) ];
end % function Compute_N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E=Compute_E(x,s); I=eye(3); Z=zeros(3);
E=[I Z; Z [s.mc+s.m1+s.m2         -s.m1*s.ell1*cos(x(2)) -s.m2*s.ell2*cos(x(3));
           -s.m1*s.ell1*cos(x(2))  s.I1+s.m1*s.ell1^2             0            ;
           -s.m2*s.ell2*cos(x(3))          0              s.I2+s.m2*s.ell2^2   ]];
end % function Compute_E

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dX=RHS_sys(A,B,E,K,x_prime,u_bar,x_bar)
    dX=E\(A+B*K)*x_prime-E\B*K*x_bar+E\B*u_bar;
end