clear all

%% Parameters
S0 = 1;
% S0 = 25.67;
T = 0.5; r = 0.05; K= 1;
addpath("../1_CarrMadam")
addpath("../3_Levy_model")
addpath("../5_MC_pricing")
addpath("../6_PIDE/Merton")
addpath('../7_Conv')
load("../9_Calibration/Params/pMer.mat")

%% Simulation Parameters
% FD parameters
Mpide=100; Npide=100; theta=.5; 

% MC parameters
Nsim = 1e5;

% CONV
N_conv = 2^13;

% M_mc = round(T*254); % Daily Monitoring
M = round(T*52); % Weekly Monitoring
% M = round(T*12); % Monthly Monitoring



%% Carr&Madam prices
P_CM_call = FFT_Merton(K,S0,T,r,pMer,1);
P_CM_put = FFT_Merton(K,S0,T,r,pMer,-1);
CP_check = P_CM_call-P_CM_put+K*exp(-T*r)-S0
CM_prices = table(P_CM_call, P_CM_put)

%% Monte Carlo Simulation
[Smer, Smer_AV] = asset_Merton(S0,T, r, pMer,Nsim,M);

% Matrix Kou Martingale Check
[Check,~,CI]=normfit( Smer(:,end)*exp(-r*T) -S0) 

%% European
% flag=1; D=-inf; U=inf;    % CALL
% flag=-1; D=-inf; U=inf;   % PUT
flag=1; D=-inf; U=inf;       % DO PUT
% flag=1; D=-inf; U=1.2;    % UP CALL


% Methods for pricing
[P_mc, ~]  = MC_European(Smer, T, r, K, flag, D, U);
P_fd = thetaMerton_eu (S0,T,r,K,pMer, theta, flag, Mpide, Npide,D,U );
% P_conv = CONV(S0,T,r,K,  pKou,4, flag, M, N_conv,D,U);

if flag==1 P_CM=P_CM_call; else P_CM=P_CM_put; end
    
method = ["CarrMadam", "Monte Carlo","fd"]';
prices = [P_CM, P_mc, P_fd]' ;
error = 100*(prices-P_mc)/P_mc;
european = table(method,prices,error)

%% American
flag=1; D=-inf; U=inf;    % CALL
% flag=-1; D=-inf; U=inf;   % PUT


[P_mc, ~]  = MC_American(Smer, T, r, K, flag);
P_fd = thetaMerton_am (S0,T,r,K,pMer, theta, flag, Mpide, Npide,D,U );

if flag==1 P_CM=P_CM_call; else P_CM=P_CM_put; end
    
method = ["CarrMadam", "Monte Carlo","fd"]';
prices = [P_CM, P_mc, P_fd]' ;
error = 100*(prices-P_mc)/P_mc;
american = table(method,prices,error)

%% Asian
flag = 1;
Strike = 0;
% Strike = K;
MC_Asian(Smer, T, r, Strike, flag)

%% Lookback
Strike = 0;
% Strike = K;
MC_Lookback(Smer, T, r, Strike, flag)

