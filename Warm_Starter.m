%% Loop Adj_Opt_Init to minimize theta 1 and theta 2
clear n_loops 
n_loops=5;
%initialize variables (~.15rad and ~0.1rad/sec)
%u_k=zeros(301,n_loops);
%x_k=zeros(9,301,n_loops);
small_angle_limit=.3; %limit for small angel approx in rad (0.1 rad convets to ~5deg)
vel_limit=.2; %limits the maximum angular velocity at the final state to 0.5rad/sec

%maintain options to pick between both stock provided code and my code modified to take the struct s as
%a funciton input. ***dynamics are wrong in the stock code, to fix, pull
%from Github and then edit with correct magnitudes.
[u_k(:,1),x_k(:,:,1)]=Adj_Opt_Init(u_k(:,end),s);
%[u_k(:,1),x_k(:,:,1)]=NR_Dual_Pendulum_Swingup(5,u_k(:,end));
noise_factor=0.10;
%figure(1); plot(1,x_k(2,301,1)); hold on; plot(1,x_k(3,301,1)); hold on; plot(1,x_k(5,301,1)); hold on; plot(1,x_k(2,301,1)); hold on;
for n=1:n_loops     %control loops based on state thresholds. This feature never triggered when state thrsholds were even remotely desireable
    if abs(x_k(2,end,n))<small_angle_limit && abs(x_k(3,end,n))<small_angle_limit && abs(x_k(5,end,n))<theta_1_dot_crit && abs(x_k(6,end,n))<theta_2_dot_crit, break % abs(x_k(5,301,n))<vel_limit, abs(x_k(6,301,n))<vel_limit
    else
        [u_k(:,n+1),x_k(:,:,n+1)]=Adj_Opt_Init(u_k(:,n),s);
    %[u_k(:,n+1),x_k(:,:,n+1)]=NR_Dual_Pendulum_Swingup(T,u_k(:,n));
    
    n       %print iteration number and states each loop because Adj_Opt_Init has the printing and the plots turned off.
    x_k(1:6,end,n)
    end

end


    