function [S, S_AV] = asset_ENIG(S0, T, r, param,Nsim,M)

% Sample VG
dt = T/M;
sigmaNIG=param(1); thetaNIG=param(2); kNIG=param(3);sigmaGBM=param(4);
 
% risk neutral
Psi=@(v) 1/kNIG-sqrt(1+v.^2*sigmaNIG^2*kNIG-2i*thetaNIG*v*kNIG)/kNIG...
    +sigmaGBM^2/2*1i*v; % without drift
drift=r-Psi(-1i);

mu = dt;
lambda = (dt)^2/kNIG;

X=zeros(Nsim,M+1);
X_AV=zeros(Nsim,M+1);
for i=1:M
    % Sample the subordinator
    DeltaS=icdf('Inverse Gaussian',rand(Nsim,1),mu,lambda);
    
    % Sample the process
    Z = randn(Nsim,1);
    bm = sigmaGBM*sqrt(dt)*randn(Nsim,1);
    X(:,i+1)=X(:,i)+bm+...
        drift*dt+thetaNIG*DeltaS+sigmaNIG*sqrt(DeltaS).*Z;
    X_AV(:,i+1)=X_AV(:,i)-bm+...
        drift*dt+thetaNIG*DeltaS+sigmaNIG*sqrt(DeltaS).*-Z;
end

S = S0*exp(X);
S_AV = S0*exp(X_AV);
end