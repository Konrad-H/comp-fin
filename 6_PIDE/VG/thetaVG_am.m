function Price =  thetaVG_am(S0,T,r,K,  pVG,theta, flag, M, N,D,U,	epsilon )
% Price a European Call (or Put) Option using FD SCHEME - THETA
% Model: VG
% Framework: General Levy - LogPrice
%% Input

%% For Kou Model
sigmaVG=pVG(1); thetaVG=pVG(2);  kVG=pVG(3); 
Anu=thetaVG/(sigmaVG^2); Bnu=sqrt(thetaVG^2+2*sigmaVG^2/kVG)/(sigmaVG^2);
nu=@(y) 1./(kVG*abs(y)).*exp(Anu*y-Bnu*abs(y));
%% TRUNCATION

if epsilon==0
    sigma=0;
else
    %%%%%%%%%%%%%%%%%%%
    %  define sigma from truncation, and define the truncated nu
    ynodes=linspace(-epsilon,epsilon,2*N);
    sigma=sqrt(trapz(ynodes,ynodes.^2.*nu(ynodes)));
    nu=@(y) nu(y).*(abs(y)>epsilon);
    %%%%%%%%%%%%%%%%%%%
end
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
ymin=-0.1-epsilon;
while nu(ymin)>tol
    ymin=ymin-0.5;
end
ymax=0.1+epsilon;
while nu(ymax)>tol
    ymax=ymax+0.5;
end
ynodes=linspace(ymin,ymax,2*N);
% figure
% plot(ynodes,nu(ynodes));

%% Construct the matrix


Ad=(1-theta)*(-(r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));
A=-1/dt+(1-theta)*(-sigma^2/(dx^2)-(r));
Au=(1-theta)*((r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));
% 1nd way
Mhat=spdiags( [Ad*ones(N-1,1) A*ones(N-1,1) Au*ones(N-1,1)],[-1 0 1],N-1,N-1);
% 2nd way
MatA=sparse(N+1,N+1); MatA(2:end-1,2:end-1)=Mhat; clear Mhat
MatA(1,1)=1; MatA(2,1)=Ad; MatA(end,end)=1; MatA(end-1,end)=Au;


Bd=-theta*(-(r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));
B=-1/dt-theta*(-sigma^2/(dx^2)-(r));
Bu=-theta*((r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));
% 1nd way
Mhat=spdiags( [Bd*ones(N-1,1) B*ones(N-1,1) Bu*ones(N-1,1)],[-1 0 1],N-1,N-1);
% 2nd way
MatB=sparse(N+1,N+1); MatB(2:end-1,2:end-1)=Mhat; clear Mhat
MatB(2,1)=Bd; MatB(end-1,end)=Bu;

%% Backward-in-time procedure
V=max( flag*(S0*exp(x)-K),0).*( S0*exp(x)>D).*( S0*exp(x)<U); 
BC=zeros(N+1,1);
for j=M-1:-1:0
    I=integral(nu,x,V,ynodes,S0,K*exp(-r*(T-(j+1)*dt)));
    
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
    b=MatB*V-I;
    b(end) = BC(end);
    b(1) = BC(1);
    %% Exploit SOR to solve the linear system
    tol=1e-4; maxiter=500; omega=1.5;
    % guess solution: V
    for iter=1:maxiter 
        Vold=V; %given Vold -> V
        for i=1:length(V)
            if i==1
                rhs=b(i)-MatA(i,i+1)*Vold(i+1);
            elseif i==length(V)
                rhs=b(i)-MatA(i,i-1)*V(i-1);
            else
                rhs=b(i)-MatA(i,i+1)*Vold(i+1)-MatA(i,i-1)*V(i-1);
            end
            V(i)=max( (1-omega)*Vold(i)+omega*rhs/MatA(i,i),...
                flag*(S0*exp(x(i))-K) );
        end
        error=norm(V-Vold,'Inf');
        if error<tol
            break
        end
    end
    [j, iter, error]
end



%% Plot
% figure
% plot(S0*exp(x),V); title('Price'); xlabel('S - spot price');
Price=interp1( S0*exp(x),V, S0,'spline');


function I=integral(nu,x,V,ynodes,S0,discK)
% trapezoidal formula
dy=ynodes(2)-ynodes(1);
w=ones(size(ynodes)); w(1)=1/2; w(end)=1/2;
% construct the first order derivative
dx=x(2)-x(1);
dV=(V(3:end)-V(1:end-2))/(2*dx);
% compute the integral
I=zeros(size(x));
for i=2:length(I)-1
    
    I(i)=sum( w.*(valueV(x(i)+ynodes,x,V,S0,discK)... %We call this A, does not converge alone
         -V(i)...   %We call this B, converges only for FA
         -(exp(ynodes)-1)*dV(i-1) )... % We call this C, converges for FA/FV
             .*nu(ynodes) )*dy;
end


function f=valueV(y,x,V,S0,discK)
f=zeros(size(y));
index=find( (y>=x(1)) .* (y<=x(end)) );
f(index)=interp1(x,V,y(index));

% % This is never used
% if (flag==1) & (U==inf)
%     index=find( (y>x(end)) );
%     f(index)=S0*exp(y(index))-discK;
% elseif (flag==-1) & (D==-inf)
%     index=find( (y<x(1)) );
%     f(index)=discK-S0*exp(y(index));
% end