function [S S_AV] = asset_NIG(S0, T, r, param,Nsim,M)

% Sample VG
dt = T/M;
sigmaNIG=param(1); thetaNIG=param(2); kNIG=param(3);
 
% risk neutral
Psi=@(v) 1/kNIG-sqrt(1+v.^2*sigmaNIG^2*kNIG-2i*thetaNIG*v*kNIG)/kNIG; % without drift
drift=r-Psi(-1i);

mu = dt;
lambda = (dt)^2/kNIG;

X=zeros(Nsim,M+1);
X_AV=zeros(Nsim,M+1);
for i=1:M
    % Sample the subordinator
    DeltaS=icdf('Inverse Gaussian',rand(Nsim,1),mu,lambda);
    % Sample the process
%     X(:,i+1)=X(:,i)+thetaVG*DeltaS+sigmaVG*sqrt(DeltaS).*randn(Nsim,1);
    Z = randn(Nsim,1);
    X(:,i+1)=X(:,i)+drift*dt+thetaNIG*DeltaS+sigmaNIG*sqrt(DeltaS).*Z;
    X_AV(:,i+1)=X_AV(:,i)+drift*dt+thetaNIG*DeltaS+sigmaNIG*sqrt(DeltaS).*-Z;
end

S = S0*exp(X);
S_AV = S0*exp(X_AV);
end