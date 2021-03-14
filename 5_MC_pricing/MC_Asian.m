function [Price,CI,DiscPayoff] = MC_Asian( St, T,r,K, flag)

if K==0     % Floating Strike
    Strike = mean(St,2);
    Payoff = @(St) max( flag*(St(:,end)-Strike), 0);
    
else        % Fixed Strike
    Payoff = @(St) max( flag*(mean(St,2)-K), 0);
end

DiscPayoff=exp(-r*T)*Payoff(St);
[Price,~,CI]=normfit(DiscPayoff);

end