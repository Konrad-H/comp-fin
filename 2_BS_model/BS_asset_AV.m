function [S, S_AV] = BS_asset_AV(S0,T,r,sigma, Msim, steps )
    

    dt = T/steps;
    
    X = zeros(Msim,steps+1);
    X_AV = zeros(Msim,steps+1);

    for i=2:(steps+1)
        Z = randn(Msim,1);
        X   (:,i)=   X (:,i-1) +(r-sigma^2/2)*dt+sigma*sqrt(dt)*Z;
        X_AV(:,i)= X_AV(:,i-1) +(r-sigma^2/2)*dt+sigma*sqrt(dt)*-Z;
    end

    S = S0*exp(X);
    S_AV = S0*exp(X_AV);
    
end