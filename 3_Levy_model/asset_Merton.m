function [S,SAV]=matrix_asset_merton_AV(S0, T, r, par,Nsim,M)
% Simulation of the Kou process
sigma=par(1);
lambdaK=par(2);
mu=par(3);
delta=par(4);


% Under Q
dt=T/M;
psi=@(u) -sigma^2/2*u.^2 +lambdaK.*(exp(-delta^2*(u.^2)/2 + 1i*mu*u)-1);
drift=r -psi(-1i); % risk neutral drift

% Simulating -> conditional simulation
NT=icdf('Poisson',rand(Nsim,1),lambdaK*T);
X=zeros(Nsim,M+1); XAV=zeros(Nsim,M+1);
Z=randn(Nsim,M);
maxJump   = (ones(Nsim,1)*(1:max(NT)) )>NT;
JumpTimes = T*(rand(Nsim,max(NT))) + maxJump*T;
JumpTimes = sort(JumpTimes,2);


for j=1:M
        X(:,j+1)=X(:,j)+drift*dt+sigma*sqrt(dt)*Z(:,j);
        XAV(:,j+1)=XAV(:,j)+drift*dt-sigma*sqrt(dt)*Z(:,j);
        
        idx_jumps = find( (JumpTimes>(j-1)*dt).*(JumpTimes<=j*dt));
        
        u = randn(length(idx_jumps),1);
        Jumps = zeros(Nsim, max(NT));
        JumpsAV = zeros(Nsim, max(NT));
        Jumps(idx_jumps)  = (mu+delta*u);
        JumpsAV(idx_jumps)= (mu-delta*u);
        
        X(:,j+1)=X(:,j+1)+sum(Jumps,2);
        XAV(:,j+1)=XAV(:,j+1)+sum(JumpsAV,2);        

end

S=S0*exp(X); SAV=S0*exp(XAV);
      
                
        
        
        
