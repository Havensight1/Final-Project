function X=DRE(A,BC,E,XP_0,T,R,Q,ricatti)
%%solves the differential Ricatti equaiton for the provided matricies and
%%weights over the given interval T using a timestep calculated from the size of A and returns
% either X or P based on if the string argument is 'K' or 'L' respectively.
%ricatti sets the ricatti equation being solved, 1=controller for X,
%0=kalman filter for P
X=zeros(size(A,1),size(A,2),size(A,3));
X(:,:,end)=XP_0;
h=T/size(A,3)/2;
if ricatti~=1 
X(:,:,1)=XP_0; X(:,:,end)=X(:,:,2);
C=zeros(size(R,2),size(A,2),size(A,3));
if size(BC,3)>1
    for i=1:size(A,3)
        C(:,:,i)=BC(:,:,i);
    end
else
        for i=1:size(A,3)
        C(:,:,i)=BC(:,:);
        end
end
end
% keyboard
if ricatti~=1 
for i=1:size(A,3)-1
        f1 = RHS_P(X(:, :, i), E(:, :, i), A(:, :, i), C(:, :, i), R, Q);
        f2 = RHS_P(X(:, :, i)+f1/2*h, E(:, :, i), A(:, :, i), C(:, :, i), R, Q);
        f3 = RHS_P(X(:, :, i)+f2/2*h, E(:, :, i), A(:, :, i), C(:, :, i), R, Q);
        f4 = RHS_P(X(:, :, i)+f3*h, E(:, :, i), A(:, :, i), C(:, :, i), R, Q);
        X(:, :, i+1) = X(:, :, i) + h*(f1/6+(f2+f3)/3+f4/6);
end
else
    for i=size(A,3):-1:2
        f1 = RHS(X(:, :, i), E(:, :, i), A(:, :, i), BC(:, :, i), R, Q);
        f2 = RHS(X(:, :, i)-f1/2*h, E(:, :, i), A(:, :, i), BC(:, :, i), R, Q);
        f3 = RHS(X(:, :, i)-f2/2*h, E(:, :, i), A(:, :, i), BC(:, :, i), R, Q);
        f4 = RHS(X(:, :, i)-f3*h, E(:, :, i), A(:, :, i), BC(:, :, i), R, Q);
        X(:, :, i-1) = X(:, :, i) - h*(f1/6+(f2+f3)/3+f4/6);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dX=RHS(X, E,A,B,R,Q) %equation 22.13a
% keyboard
dX = - (E')^-1 * A' * X - X * A * (E)^-1 + X * B * R^-1 * B' * X - (E')^-1 * Q * E^-1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dX=RHS_P(X,E,A,C,R,Q) %equation 22.30
% keyboard
dX =(E')^-1 * A' * X + X * A * (E)^-1 - X * C' * R^-1 * C * X + (E')^-1 * Q * E^-1;
end