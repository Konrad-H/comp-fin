function price = CrrBs_eu(S0,K,r,T,sigma, flag,  M,D, U )
%% Binomial method (Cox Ross Rubinstein 1979)
% Price an European Call Option

%% 1. Parameters
if nargin<6
    flag = 1; % 1 for call option, -1 for put option
end
if nargin<7
    D=-inf; % number of time steps
end
if nargin<8
    U=inf; % number of time steps
end
if nargin<9
    M=5000; % number of time steps
end

f_Payoff = @(S) max( flag*(S-K),0 );
dt=T/M;

%% 2. Tree Generation
u=exp(sigma*sqrt(dt)); d=1/u;
ST=S0*u.^[M:-2:-M];
%% 3. Payoff

Payoff=f_Payoff(ST).*(ST>D).*(ST<U);

%% 4. Backward-in-time formula to price the derivative
CallTree=Payoff';
q=(exp( r*dt)-d)/(u-d);

for j=M:-1:1
    Sj = S0*u.^[(j-1):-2:-(j-1)]';    
    CallTree=(Sj>D).*(Sj<U).*exp(-r*dt).*...
        (q* CallTree(1:j)+(1-q)*CallTree(2:(j+1)));


end
price=CallTree;
    

end