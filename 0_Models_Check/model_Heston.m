clear all

%% Parameters
S0 = 1;
% S0 = 25.67;
T = 0.5; r = 0.05; K= 0.3;
load("../9_Calibration/Params/pHes.mat")
addpath("../5_MC_pricing")
addpath("../4_LV_model")
addpath("../1_CarrMadam")

%% Simulation Parameters
%FD parameters
Mpide=200; Npide=200; theta=.5; 
%MC parameters
M = round(T*12); % Daily Monitoring
Nsim = 1e5;

%% Carr&Madam prices
P_CM_call = FFT_Heston(K,S0,T,r,pHes,1);
P_CM_put = FFT_Heston(K,S0,T,r,pHes,-1);
CP_check = P_CM_call-P_CM_put+K*exp(-T*r)-S0

%% Monte Carlo Simulation - QE SCHEME
Shes = asset_Heston(S0,T, r, pHes,Nsim,M); %doesnt perform well
% Shes = [S_nAV; S_AV];
% Heston Martingale Check
[Check,~,IC]=normfit( Shes(:,end)*exp(-r*T) -S0) 

%% Monte Carlo Simulation - EULER SCHEME
Shes2 = asset_Heston2(S0,T, r, pHes,Nsim,M);
% Heston Martingale Check
[Check,~,IC]=normfit( Shes2(:,end)*exp(-r*T) -S0)


%% European
flag=1; D=-inf; U=inf;    % CALL
% flag=-1; D=-inf; U=inf;       % PUT
% flag=1; D=0.01; U=inf;    % UP CALL
% flag=-1; D=-inf; U=1.1;   % DO PUT

% Methods for pricing
[P_mc, ~]  = MC_European(Shes, T, r, K, flag, D, U);
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
[P_mc, ~]  = MC_American(Shes, T, r, K, flag);
if flag==1 P_CM=P_CM_call; else P_CM=P_CM_put; end
    
method = ["CarrMadam", "Monte Carlo"]';
prices = [P_CM, P_mc]' ;
error = 100*(prices-P_mc)/P_mc;
european = table(method,prices,error)