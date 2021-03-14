function [Price,CI,DiscPayoff] = MC_European( St, T,r, K, flag, D, U)

knock = (min(St,[],2)>D).*(max(St,[],2)<U);
DiscPayoff=exp(-r*T)*max( flag*(St(:,end)-K), 0).*knock;
[Price,~,CI]=normfit(DiscPayoff);

end

% in = (max(St,[],2)>U)