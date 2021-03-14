function price = crrBs_am(S0,K,r,T,sigma, flag, M,D, U,div )
%% Binomial method (Cox Ross Rubinstein 1979)
% Price an European Call Option


%% 1. Parameters

if nargin<8
    D=-inf; % number of time steps
end
if nargin<9
    U=inf; % number of time steps
end
if nargin <10
    div=0;
end

f_Payoff = @(S) max( flag*(S-K),0 ).*(S>D).*(S<U);
dt=T/M;


f_IV = @(S) max(flag*(S-K),0);

%% 2. Tree Generation
u=exp(sigma*sqrt(dt)); d=1/u;
ST=S0*u.^[M:-2:-M];
%% 3. Payoff
Payoff = f_Payoff(ST);

%% 4. Backward-in-time formula to price the derivative
CallTree=Payoff;
r_div = r-div;
q=(exp( r_div*dt)-d)/(u-d);

for j=M:-1:1
    CV=exp(-r*dt)* (q* CallTree((1:j))+(1-q)*CallTree((1:j)+1));
    Sj=S0*u.^( (j-1):-2:-(j-1) );
    IV= f_IV(Sj);
    CallTree=max(IV,CV).*(Sj>D).*(Sj<U);
end
price=CallTree;
   
end