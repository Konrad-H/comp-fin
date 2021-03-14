function [S,SAV]=matrix_asset_kou_AV(S0, T, r, par,Nsim,M)
% Simulation of the Kou process
sigma=par(1);
p=par(2);
lambdaK=par(3);
lambdap=par(4);
lambdam=par(5);


% Under Q
dt=T/M;
psi=@(u) -sigma^2/2*u.^2+1i*u*lambdaK.*(p./(lambdap-1i*u)-(1-p)./(lambdam+1i*u));
drift=r-psi(-1i); % risk neutral drift

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
        
        which_jumps  = (JumpTimes>(j-1)*dt).*(JumpTimes<=j*dt);
        idx_jumps = find(which_jumps );
        u_jump = rand(length(idx_jumps),1);
%         prob = (u_jump<p);
        idx_pos = idx_jumps(u_jump<p);
        idx_neg = idx_jumps( u_jump>p );        
%         which_pos = which_jumps*(prob);
%         which_neg = which_jumps*(1-prob);
%         idx_pos = find(which_pos);
%         idx_neg = find(which_neg);
        u_pos = rand(length(idx_pos),1);
        u_neg =rand(length(idx_neg),1);
        Jumpp = icdf('Exponential',u_pos,1/lambdap);
        Jumpm = -icdf('Exponential',u_neg,1/lambdam);
        
        JumppAV = icdf('Exponential',1-u_pos,1/lambdap);
        JumpmAV =- icdf('Exponential',1-u_neg,1/lambdam);
        
        Jumps = zeros(Nsim, max(NT));
        JumpsAV = zeros(Nsim, max(NT));
        
        Jumps(idx_pos)  = Jumpp;
        Jumps(idx_neg)= Jumpm;
        JumpsAV(idx_pos)  = JumppAV;
        JumpsAV(idx_neg)= JumpmAV;
        
        X(:,j+1)=X(:,j+1)+sum(Jumps,2);
        XAV(:,j+1)=XAV(:,j+1)+sum(JumpsAV,2);        

end

S=S0*exp(X); SAV=S0*exp(XAV);
      
                
        
        
        
