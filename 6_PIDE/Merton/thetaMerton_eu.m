function Price =  thetaMerton_eu(S0,T,r,K,  pMer,theta, flag, M, N,D,U )
% Price a European Call (or Put) Option using FD SCHEME - THETA
% Model: Mer
% Framework: General Levy - LogPrice

%% For Merton Model
sigma=pMer(1); lambdaK=pMer(2); mu=pMer(3); delta=pMer(4);

nu=@(y) lambdaK/(delta*sqrt(2*pi))*exp(-(y-mu).^2/(2*delta^2));
%% Grids
dt=T/M; 
% xmin=(r-sigma^2/2)*T-6*sigma*sqrt(T);
% xmax=(r-sigma^2/2)*T+6*sigma*sqrt(T);
Smin = max(0.2*S0,D);
Smax = min(3*S0,U);
xmin=log(Smin/S0); xmax=log(Smax/S0);
dx=(xmax-xmin)/N;
pMer=linspace(xmin, xmax, N+1)';

%% Computing alpha, lambda, and truncating the integral
% Truncate the integral
tol=1e-14;
ymin=-0.1;
while nu(ymin)>tol
    ymin=ymin-0.5;
end
ymax=0.1;
while nu(ymax)>tol
    ymax=ymax+0.5;
end
ynodes=linspace(ymin,ymax,2*N);
% figure
% plot(ynodes,nu(ynodes));
alpha=trapz(ynodes, (exp(ynodes)-1).*nu(ynodes));
lambdaNum=trapz(ynodes, nu(ynodes));
% COMMENT: If we use the actual lambdaNum and not lambda for the FD 
% equation we can get away with using a much smaller N, discretization on
% price

%% Construct the matrix


Ad=(1-theta)*(-(r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2));
A=-1/dt+(1-theta)*(-sigma^2/(dx^2)-(r+lambdaNum));
Au=(1-theta)*((r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2));
% 1nd way
Mhat=spdiags( [Ad*ones(N-1,1) A*ones(N-1,1) Au*ones(N-1,1)],[-1 0 1],N-1,N-1);
% 2nd way
MatA=sparse(N+1,N+1); MatA(2:end-1,2:end-1)=Mhat; clear Mhat
MatA(1,1)=1; MatA(2,1)=Ad; MatA(end,end)=1; MatA(end-1,end)=Au;


Bd=-theta*(-(r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2));
B=-1/dt-theta*(-sigma^2/(dx^2)-(r+lambdaNum));
Bu=-theta*((r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2));
% 1nd way
Mhat=spdiags( [Bd*ones(N-1,1) B*ones(N-1,1) Bu*ones(N-1,1)],[-1 0 1],N-1,N-1);
% 2nd way
MatB=sparse(N+1,N+1); MatB(2:end-1,2:end-1)=Mhat; clear Mhat
MatB(2,1)=Bd; MatB(end-1,end)=Bu;

%% Backward-in-time procedure
V=max( flag*(S0*exp(pMer)-K),0).*( S0*exp(pMer)>D).*( S0*exp(pMer)<U); 
BC=zeros(N+1,1);
for j=M-1:-1:0
    I=integral(nu,pMer,V,ynodes,S0,K*exp(-r*(T-(j+1)*dt)));
    
    if flag==1 %Boundary Conditions
        if U==inf
            BC(end)= ( S0*exp(xmax)-K*exp(-r*(T-j*dt)) );
        end
        BC(1) = 0; % For Call option
    elseif flag==-1
        BC(end) = 0;
        if D==-inf
            BC(1)= - ( S0*exp(xmax)-K*exp(-r*(T-j*dt))  ); % For Put option
        end
    end
    rhs=MatB*V-I;
    rhs(end) = BC(end);
    rhs(1) = BC(1);
    V=MatA\rhs;
end



%% Plot
% figure
% plot(S0*exp(x),V); title('Price'); xlabel('S - spot price');
Price=interp1( S0*exp(pMer),V, S0,'spline');


function I=integral(nu,x,V,ynodes,S0,discK)
% trapezoidal formula
dy=ynodes(2)-ynodes(1);
w=ones(size(ynodes)); w(1)=1/2; w(end)=1/2;
I=zeros(size(x));
for i=2:length(I)-1
    I(i)=sum( w.*valueV(x(i)+ynodes,x,V,S0,discK).*nu(ynodes) )*dy;
end


function f=valueV(y,x,V,S0,discK)
f=zeros(size(y));
index=find( (y>=x(1)) .* (y<=x(end)) );
f(index)=interp1(x,V,y(index));

% What is this? Help me
if (flag==1) & (U==inf)
    index=find( (y>x(end)) );
    f(index)=S0*exp(y(index))-discK;
elseif (flag==-1) & (D==-inf)
    index=find( (y<x(1)) );
    f(index)=discK-S0*exp(y(index));
end

  