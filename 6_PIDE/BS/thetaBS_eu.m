function Price = ThetaEu_BS(S0,T,r,K,  sigma,theta, flag, M, N,D,U )
%THETA Price an European Option with finite difference using Theta Method
%   Detailed explanation goes here

%% Input




%% Grids
dt=T/M;
if D==-inf    
    xmin=max(D, (r-sigma^2/2)*T-6*sigma*sqrt(T));
else
    xmin = log(D/S0); 
end

if U==inf
    xmax=min(U, (r-sigma^2/2)*T+6*sigma*sqrt(T));
else
    xmax = log(U/S0);   
end

dx=(xmax-xmin)/N;
x=linspace(xmin, xmax, N+1)';

%% Construct the matrix matA and matB
Ad=(1-theta)*(-(r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));
A=-1/dt+(1-theta)*(-sigma^2/(dx^2)-r);
Au=(1-theta)*((r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));
% 1nd way
Mhat=spdiags( [Ad*ones(N-1,1) A*ones(N-1,1) Au*ones(N-1,1)],[-1 0 1],N-1,N-1);
% 2nd way
MatA=sparse(N+1,N+1); MatA(2:end-1,2:end-1)=Mhat; clear Mhat
MatA(1,1)=1; MatA(2,1)=Ad; MatA(end,end)=1; MatA(end-1,end)=Au;


Bd=-theta*(-(r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));
B=-1/dt-theta*(-sigma^2/(dx^2)-r);
Bu=-theta*((r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));
% 1nd way
Mhat=spdiags( [Bd*ones(N-1,1) B*ones(N-1,1) Bu*ones(N-1,1)],[-1 0 1],N-1,N-1);
% 2nd way
MatB=sparse(N+1,N+1); MatB(2:end-1,2:end-1)=Mhat; clear Mhat
MatB(2,1)=Bd; MatB(end-1,end)=Bu;

%% Backward-in-time procedure
V=max( flag*( S0*exp(x)-K ),0).*( S0*exp(x)>D).*( S0*exp(x)<U); % add minus for put option
BC=zeros(N+1,1);
for j=M-1:-1:0
    if flag==1 %Boundary Conditions
        if U==inf %only if not up and out
            BC(end)=(S0*exp(xmax)-K*exp(-r*(T-j*dt)) ); % For Call option
        end
    elseif flag==-1
        if D==-inf % Only if not down and out
            BC(1)=-(S0*exp(xmin)-K*exp(-r*(T-j*dt)) ); % For Put option
        end
    end
    rhs=MatB*V+BC;
    V=MatA\rhs;
end

%% Plot
% figure
% plot(S0*exp(x),V); title('Price'); xlabel('S - spot price');
Price=interp1( S0*exp(x),V, S0,'spline');
end

