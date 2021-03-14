function Price =  explicitEuler(S0,Strike,Rate,Time,Volatility, flag, M, N )
% Price a European Call (or Put) Option in the B&S model using Explicit
% Euler
%% Input
if nargin <8
    N = 5e2;
end
if nargin <7
    M=5e5;
end
if nargin <6
    flag = 1;
end

%% Grids
dt=Time/M; 
xmin=(Rate-Volatility^2/2)*Time-6*Volatility*sqrt(Time);
xmax=(Rate-Volatility^2/2)*Time+6*Volatility*sqrt(Time);
dx=(xmax-xmin)/N;
x=linspace(xmin, xmax, N+1)';

%% Backward-in-time procedure
V=max(flag*( S0*exp(x)-Strike) ,0); Vnew=zeros(size(V));
for j=M-1:-1:0
    % BC
    Vnew(1)=0;
    % PDE
    Vnew(2:end-1)=V(2:end-1)+dt*( (Rate-Volatility^2/2)/(2*dx)*(V(3:end)-V(1:end-2))... 
    +Volatility^2/2/(dx^2)*(V(3:end)-2*V(2:end-1)+V(1:end-2))...
    -Rate*V(2:end-1) );
    % BC
    Vnew(end)=S0*exp(xmax)-Strike*exp(-Rate*(M-j)*dt);
    % Replace the old vector with the new one
    V=Vnew;
        
end

%% Plot & Return

% figure
% plot(S0*exp(x),V); title('Price'); xlabel('S - spot price');

Price=interp1( S0*exp(x),V, S0,'spline');

end