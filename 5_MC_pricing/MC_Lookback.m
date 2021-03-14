function [Price,CI,DiscPayoff] = MC_Lookback( St, T,r, K,flag)

m = min(St,[],2);
M = max(St,[],2);

if K==0     % Floating Strike
    switch flag 
        case 1 
            Strike=m; 
        case -1 
            Strike=M;
    end
    Payoff = @(St) max( flag*(St(:,end)-Strike), 0);
    
else        % Fixed Strike
    switch flag 
        case 1 
            Payoff = @(St) max( flag*(M-K), 0);
        case -1 
            Payoff = @(St) max( flag*(m-K), 0);
    end
    
end


DiscPayoff=exp(-r*T)*Payoff(St);
[Price,~,CI]=normfit(DiscPayoff);

end