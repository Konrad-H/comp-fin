clear all

%% Parameters
S0 = 1;
% S0 = 25.67;
T = 0.5; r = 0.05; K= 0.3;
addpath("../1_CarrMadam")
addpath("../5_MC_pricing")
addpath("../4_LV_model")
addpath("../6_PIDE/EVG")
load("../9_Calibration/Params/pBat.mat")
%% Heston to Bates
% if we want to use the Heston model with Bates we only need to change
% parameters
% load("../Calibration/Params/pHes.mat")
% pBat = [pHes, 0, 0 1]

%% Simulation Parameters
% FD parameters
Mpide=200; Npide=200; theta=.5; 

% MC parameters
Nsim_mc = 1e5;
% M_mc = round(T*254); % Daily Monitoring
% M_mc = round(T*52); % Weekly Monitoring
M_mc = round(T*12); % Monthly Monitoring

%% Pricing
P_CM_call = FFT_Bates(K,S0,T,r,pBat,1);
P_CM_put = FFT_Bates(K,S0,T,r,pBat,-1);
CP_check = P_CM_call-P_CM_put+K*exp(-T*r)-S0
CM_prices = table(P_CM_call, P_CM_put)

%% Monte Carlo Simulation - EULER SCHEME
Sbat = asset_Bates(S0,T, r, pBat,Nsim_mc,M_mc);
% Bates Martingale Check
[Check,~,IC]=normfit( Shes2(:,end)*exp(-r*T) -S0)

%% European
flag=1; D=-inf; U=inf;    % CALL
% flag=-1; D=-inf; U=inf;       % PUT
% flag=1; D=0.01; U=inf;    % UP CALL
% flag=-1; D=-inf; U=1.1;   % DO PUT

% Methods for pricing
[P_mc, ~]  = MC_European(Sbat, T, r, K, flag, D, U);
% P_conv = CONV(S0,T,r,K,  ~,~, flag, M, N_conv,D,U);
P_conv = 0;
if flag==1 P_CM=P_CM_call; else P_CM=P_CM_put; end
    
method = ["CarrMadam", "Monte Carlo", "conv"]';
prices = [P_CM, P_mc, P_conv]' ;
error = 100*(prices-P_mc)/P_mc;
european = table(method,prices,error)

%% American
flag=1; D=-inf; U=inf;    % CALL
% flag=-1; D=-inf; U=inf;       % PUT
% flag=1; D=0.01; U=inf;    % UP CALL
% flag=-1; D=-inf; U=1.1;   % DO PUT

% Methods for pricing
[P_mc, ~]  = MC_American(Sbat, T, r, K, flag);
if flag==1 P_CM=P_CM_call; else P_CM=P_CM_put; end
    
method = ["CarrMadam", "Monte Carlo"]';
prices = [P_CM, P_mc]' ;
error = 100*(prices-P_mc)/P_mc;
european = table(method,prices,error)





