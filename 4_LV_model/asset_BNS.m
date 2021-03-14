function S=asset_Heston(X0,T,r,x,Nsim,M)

% PARAMETERS
drift=x(1); % vol-of-vol
pho=x(2); % mean reversion speed
lambda=x(3); % correlation


rn_drift = r-drift;


X = zeros(Nsim, M+1);
V = zeros(Nsim, M+1);
Z = zerps(Nsim, M+1);
for i=1:M

    % Levy Increments
    dZ

    % Main Process
    X(:,i+1) = X(:,i)+drift*dt+sqrt(V(:,i))*rand(Nsim,1)+pho*dZ;
    V(:,i+1) = V(:,i) - lambda*V(:,i)*dt + dZ;  
    

S=exp(X);