function Price =  implicitEuler(S0,Strike,Rate,Time,Volatility, flag, M, N )
% Price a European Call (or Put) Option in the B&S model using Explicit
% Euler
%% Input
if nargin <8
    N=1e3;
end
if nargin <7
    M = 1e3;
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

%% Construct the matrix
A=-(Rate-Volatility^2/2)/(2*dx)+Volatility^2/(2*dx^2);
B=-1/dt-Volatility^2/(dx^2)-Rate;
C=(Rate-Volatility^2/2)/(2*dx)+Volatility^2/(2*dx^2);
% 1nd way
Mhat=spdiags( [A*ones(N-1,1) B*ones(N-1,1) C*ones(N-1,1)],[-1 0 1],N-1,N-1);
% 2nd way
Mat=sparse(N+1,N+1); Mat(2:end-1,2:end-1)=Mhat; clear Mhat
Mat(1,1)=1; Mat(2,1)=A; Mat(end,end)=1; Mat(end-1,end)=C;

%% Backward-in-time procedure
V=max( flag*(S0*exp(x)-Strike),0); 
BC=zeros(N+1,1);
for j=M-1:-1:0
    if flag==1 %Boundary Conditions
        BC(end)=(S0*exp(xmax)-Strike*exp(-Rate*(Time-j*dt)) ); % For Call option
    elseif flag==-1
        BC(1)=-(S0*exp(xmax)-K*exp(-r*(T-j*dt)) ); % For Put option
    end
    rhs=-1/dt*V+BC;
    V=Mat\rhs;
end


%% Plot
% figure
% plot(S0*exp(x),V); title('Price'); xlabel('S - spot price');
Price=interp1( S0*exp(x),V, S0,'spline');

end