function Price =  Kou_implicitEuler(S0,T,r,K,U,D, par, flag, M, N )
% Price a European Call (or Put) Option in the B&S model using Explicit
% Euler
%% Input
if nargin <10
    N=1e3;
end
if nargin <9
    M = 1e3;
end
if nargin <8
    flag = 1;
end

%% For Kou Model
sigma=par(1);
p=par(2);
lambda=par(3);
lambdap=par(4);
lambdam=par(5);
nu=@(y) lambda*(p*lambdap*exp(-lambdap*y).*(y>0)+...
               (1-p)*lambdam*exp(-lambdam*abs(y)).*(y<0));

%% Grids
dt=T/M; 
% xmin=(r-sigma^2/2)*T-6*sigma*sqrt(T);
% xmax=(r-sigma^2/2)*T+6*sigma*sqrt(T);
Smin = max(0.2*S0,D);
Smax = min(3*S0,U);
xmin=log(Smin/S0); xmax=log(Smax/S0);
dx=(xmax-xmin)/N;
x=linspace(xmin, xmax, N+1)';

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


%% Construct the matrix
Au=-(r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2);
A=-1/dt-sigma^2/(dx^2)-(r+lambdaNum);
Ad=(r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2);
% 1nd way
Mhat=spdiags( [Au*ones(N-1,1) A*ones(N-1,1) Ad*ones(N-1,1)],[-1 0 1],N-1,N-1);
% 2nd way
Mat=sparse(N+1,N+1); Mat(2:end-1,2:end-1)=Mhat; clear Mhat
Mat(1,1)=1; Mat(2,1)=Au; Mat(end,end)=1; Mat(end-1,end)=Ad;

%% Backward-in-time procedure
V=max( flag*(S0*exp(x)-K),0); 
BC=zeros(N+1,1);
for j=M-1:-1:0
    I=integral(nu,x,V,ynodes,S0,K*exp(-r*(T-(j+1)*dt)));
    rhs = -I-1/dt*V;
    if flag==1 %Boundary Conditions
        rhs(end)= (S0*exp(xmax)-K*exp(-r*(T-j*dt)) );
        rhs(1) = 0; % For Call option
    elseif flag==-1
        rhs(end) = 0;
        rhs(1)= - (S0*exp(xmax)-K*exp(-r*(T-j*dt)) ); % For Put option
    end
    
    V=Mat\rhs;
end



%% Plot
% figure
% plot(S0*exp(x),V); title('Price'); xlabel('S - spot price');
Price=interp1( S0*exp(x),V, S0,'spline');


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
index=find( (y>x(end)) );
f(index)=S0*exp(y(index))-discK;
  