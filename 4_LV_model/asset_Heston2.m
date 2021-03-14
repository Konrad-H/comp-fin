function S=asset_Heston2(S0,T,r,x,Nsim,M)
% Implementation of Euler Scheme for Heston model
% "INEfficient Simulation of the Heston Stochastic Volatility Model"

% PARAMETERS
epsilon=x(1); % vol-of-vol
k=x(2); % mean reversion speed
rho=x(3); % correlation
theta=x(4); % mean
V0=x(5);

dt=T/M; 

% MC simulation
t=linspace(0,T,M+1); % time grid
S=zeros(Nsim,M+1); S(:,1)=S0;
V=zeros(Nsim,M+1); V(:,1)=V0;
mu=[0;0]; VC=[1 rho;rho 1];
for i=1:M
    Z=mvnrnd(mu,VC,Nsim); 
    S(:,i+1)= exp(log(S(:,i))+r*dt-.5*V(:,i)*dt+sqrt(V(:,i)*dt).*Z(:,1));
    V(:,i+1)=max(0, V(:,i)+k*(theta-V(:,i))*dt+...
        epsilon*sqrt(V(:,i)*dt).*Z(:,2) );
end
% figure; hold on
% subplot(1,2,1); plot(t,S); title('Price');
% subplot(1,2,2); plot(t,V); title('Variance');


    