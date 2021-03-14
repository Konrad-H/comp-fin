function [S S_AV] = asset_VG(S0, T, r, param,Nsim,M)

% Sample VG

dt = T/M;
sigmaVG=param(1); thetaVG=param(2); kVG=param(3);

% risk neutral
Psi= @(u) -1/kVG*log(1+u.^2*sigmaVG^2*kVG/2-1i*thetaVG*kVG*u);
drift=r-Psi(-1i);

a=dt/kVG;

X=zeros(Nsim,M+1);
X_AV = zeros(Nsim,M+1);
for i=1:M
    % Sample the subordinator
    DeltaS=kVG*icdf('Gamma',rand(Nsim,1),a,1);
    % Sample the process
%     X(:,i+1)=X(:,i)+thetaVG*DeltaS+sigmaVG*sqrt(DeltaS).*randn(Nsim,1);
    Z = randn(Nsim,1);
    X(:,i+1)=X(:,i)+drift*dt+thetaVG*DeltaS+sigmaVG*sqrt(DeltaS).*Z;
    X_AV(:,i+1)=X_AV(:,i)+drift*dt+thetaVG*DeltaS+sigmaVG*sqrt(DeltaS).*-Z;
end

S = S0*exp(X);
S_AV = S0*exp(X_AV);
end
    