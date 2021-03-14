function Price = ThetaAm_BS(S0,T,r,K,  sigma,theta, flag, M, N,D,U )
%THETA Price an European Option with finite difference using Theta Method
%   Detailed explanation goes here

%% Grids
dt=T/M; 

if D==-inf    
    xmin=(r-sigma^2/2)*T-6*sigma*sqrt(T);
else
    xmin = log(D/S0); 
end
if U==inf
    xmax=(r-sigma^2/2)*T+6*sigma*sqrt(T);
else
    xmax = log(U/S0);
end

dx=(xmax-xmin)/N;
x=linspace(xmin, xmax, N+1)';

%% Construct the matrix
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
V=max( flag*(S0-K)*exp(x),0).*( S0*exp(x)>D).*( S0*exp(x)<U); 
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
    b=MatB*V+BC;
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
% figure
% plot(S0*exp(x),V); title('Price'); xlabel('S - spot price');
% hold on
% plot(S0*exp(x),max( K-S0*exp(x),0))
Price=interp1( S0*exp(x),V, S0,'spline');
    
    
end

