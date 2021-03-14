function [Price,CI,DiscPayoff] = MC_BarrierIn( St, T,r, K, flag, D, U)
if nargin <7
    U = inf;
end
if nargin <6
    D = -inf;
end

if D>-inf
    knockin = (min(St,[],2)<D);
elseif U<inf
    knockin = (max(St,[],2)>U);
else
    knockin = 1;
end
DiscPayoff=exp(-r*T)*max( flag*(St(:,end)-K), 0).*knockin;
[Price,~,CI]=normfit(DiscPayoff);

end

% in = (max(St,[],2)>U)