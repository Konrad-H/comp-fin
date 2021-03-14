function S=asset_Bates2(S0,T,r,x,Nsim,M)
% Implementation of Euler Scheme for Heston model
% "INEfficient Simulation of the Heston Stochastic Volatility Model"

% PARAMETERS
epsilon=x(1); % vol-of-vol
k=x(2); % mean reversion speed
rho=x(3); % correlation
theta=x(4); % mean
V0=x(5);
lambdaK=x(6);
mu = x(7);
delta = x(8);


psi=@(u) lambdaK.*(exp(-delta^2*(u.^2)/2 + 1i*mu*u)-1);
drift=r -psi(-1i); % Merton Drift

dt=T/M; 

% MC simulation
t=linspace(0,T,M+1); % time grid
S=zeros(Nsim,M+1); S(:,1)=S0;
V=zeros(Nsim,M+1); V(:,1)=V0;
mu_gbm=[0;0]; VC_gbm=[1 rho;rho 1];

%% Jump Part
NT=icdf('Poisson',rand(Nsim,1),lambdaK*T);
maxJump   = (ones(Nsim,1)*(1:max(NT)) )>NT;
JumpTimes = T*(rand(Nsim,max(NT))) + maxJump*T;
JumpTimes = sort(JumpTimes,2);

for i=1:M
    Z=mvnrnd(mu_gbm,VC_gbm,Nsim); 
    
    idx_jumps = find( (JumpTimes>(j-1)*dt).*(JumpTimes<=j*dt));
    u = randn(length(idx_jumps),1);
    Jumps = zeros(Nsim, max(NT));
    Jumps(idx_jumps)  = (mu+delta*u);
    
    
    S(:,i+1)= exp(log(S(:,i))+drift*dt-.5*V(:,i)*dt+sqrt(V(:,i)*dt).*Z(:,1)...
        +sum(Jumps,2));
    
    V(:,i+1)=max(0, V(:,i)+k*(theta-V(:,i))*dt+...
        epsilon*sqrt(V(:,i)*dt).*Z(:,2) );
end
% figure; hold on
% subplot(1,2,1); plot(t,S); title('Price');
% subplot(1,2,2); plot(t,V); title('Variance');


    