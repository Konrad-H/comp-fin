clear all;
%% Paths
addpath("../1_CarrMadam")
addpath("../3_Levy_model")
addpath("../5_MC_pricing")
addpath("../6_PIDE/NIG")
addpath('../7_Conv')
load("../9_Calibration/Params/pNIG.mat")

%% Parameters
S0 = 1;
% S0 = 25.67;
T = 0.5; r = 0.05; K= 1;

%% Simulation Parameters
% FD parameters
Mpide=100; Npide=100; theta=.5; 
epsilon=0;
% MC parameters
Nsim = 1e4;

% CRR parameters
M_crr = 500;

% CONV
N_conv = 2^13;

% M = round(T*254*8); % Hourly Monitoring
% M = round(T*254); % Daily Monitoring
M = round(T*52); % Weekly Monitoring
% M = round(T*12); % Monthly Monitoring

%% Theoretical Price
P_CM_call = FFT_NIG( K,S0, T, r, pNIG,1);
P_CM_put = FFT_NIG( K,S0, T, r, pNIG,-1);

table(P_CM_call, P_CM_put)
%% Monte Carlo Simulation
% BS Martingale Check
q = 0; % No dividends
r_q = r-q;
[Snig, Snig_AV] =asset_NIG(S0,T, r_q, pNIG,Nsim,M);

[Check,~,CI]=normfit( Snig(:,end)*exp(-r*T) -S0)

%% European
% flag = 1;
flag = -1;
D=-inf; U=inf;    % European
 D=-inf; U=inf;   % PUT
% D=-inf; U=inf;       % DO PUT
% D=-inf; U=1.2;    % UP CALL

% Methods for pricing
[P_mc, ~]  = MC_European(Snig, T, r, K, flag, D, U);
P_fd = thetaENIG_eu (S0,T,r,K,[pNIG 0], theta, flag, Mpide, Npide,D,U, epsilon);

P_conv = CONV(S0,T,r,K,  pNIG,3, flag, M, N_conv,D,U);

if flag==1 P_CM=P_CM_call; else P_CM=P_CM_put; end
    
method = ["CarrMadam", "Monte Carlo","conv","fd"]';
prices = [P_CM, P_mc, P_conv, P_fd]' ;
P_mean = mean([P_mc, P_conv,P_fd]);
error = 100*(prices-P_mean)/P_CM;
european = table(method,prices,error)

%% American
% flag=1; D=-inf; U=inf;    % CALL
flag=-1; D=-inf; U=inf;   % PUT


[P_mc, ~]  = MC_American(Snig, T, r, K, flag);
P_fd = thetaENIG_am (S0,T,r,K,[pNIG 0], theta, flag, Mpide, Npide,D,U,epsilon );

if flag==1 P_CM=P_CM_call; else P_CM=P_CM_put; end
    
method = ["CarrMadam", "Monte Carlo","fd"]';
prices = [P_CM, P_mc, P_fd]' ;
P_mean = mean([P_mc, P_fd]);
error = 100*(prices-P_mean)/P_mean;
european = table(method,prices,error)

%% Asian
flag = 1;
Strike = 0;
% Strike = K;
MC_Asian(Snig, T, r, Strike, flag)

%% Lookback
Strike = 0;
% Strike = K;
MC_Lookback(Snig, T, r, Strike, flag)



%%%%
%%%%
%%%%
%% Confrontiamo COV
thetaIG=pNIG(1); sigmaIG=pNIG(2); kIG=pNIG(3); 
temp=1/(kIG*sigmaIG^2); %alpha^2-beta^2 
beta=thetaIG*kIG*temp; delta=1/(kIG*sqrt(temp)); alpha=sqrt(temp+beta^2);
pNIG2=[alpha, beta, delta]

P_conv = CONV(S0,T,r,K,  pNIG2,2, flag, 50, N_conv,D,U)
P_conv = CONV(S0,T,r,K,  pNIG,3, flag, 50, N_conv,D,U)


