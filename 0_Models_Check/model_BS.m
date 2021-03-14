clear all
% 
%% Parameters
S0 = 1;
% S0 = 25.67;
T = 0.5; r = 0.05; K= 1;
addpath("../1_CarrMadam")
addpath("../2_BS_model")
addpath("../5_MC_pricing")
addpath("../6_PIDE/BS")
addpath('../7_Conv')
addpath("../8_CRR")
load("../9_Calibration/Params/pBs.mat")

%% Theoretical Price
[P_CM_call, P_CM_put] = blsprice(S0, K, r, T, sigmaBS);
table(P_CM_call, P_CM_put)
%% Simulation Parameters
% FD parameters
Mpide=100; Npide=100; theta=.5; 

% MC parameters
Nsim_mc = 1e4;

% CRR parameters
M_crr = 5000;

% CONV
N_conv = 2^12;
M = round(T*254*8); % Hourly Monitoring
% M = round(T*254); % Daily Monitoring
% M = round(T*52); % Weekly Monitoring
% M = round(T*12); % Monthly Monitoring

%% Monte Carlo Simulation
% BS Martingale Check
q = 0; % No dividends
r_q = r-q;
[Sbs, Sbs_AV] =BS_asset_AV(S0,T, r_q, sigmaBS,Nsim_mc,M);

[Check,~,CI]=normfit( Sbs(:,end)*exp(-r*T) -S0)

%% European
% flag=1; D=-inf; U=inf;    % CALL
% flag=-1; D=-inf; U=inf;   % PUT
flag=1; D=.98; U=inf;       % DO PUT
% flag=1; D=-inf; U=1.1;    % UP CALL


[P_mc, ~]  = MC_European(Sbs, T, r, K, flag, D, U);
P_fd = thetaBS_eu (S0,T,r,K,sigmaBS, theta, flag, Mpide, Npide,D,U );
P_conv = CONV(S0,T,r,K,  sigmaBS,1, flag, M, N_conv,D,U);
P_crr = crrBs_eu(S0,K,r,T,sigmaBS, flag,  M_crr,D, U );


if flag==1 P_CM=P_CM_call; else P_CM=P_CM_put; end
    
method = ["CarrMadam", "Monte Carlo","conv","fd","crr"]';
prices = [P_CM, P_mc, P_conv, P_fd, P_crr]' ;
error = 100*(prices-P_mc)/P_mc;
european = table(method,prices,error)

%% American
% flag=1; D=-inf; U=inf;    % CALL
flag=-1; D=-inf; U=inf;   % PUT


[P_mc, ~]  = MC_American(Sbs, T, r, K, flag);
P_fd = thetaBS_am (S0,T,r,K,sigmaBS, theta, flag, Mpide, Npide,D,U );
P_crr = crrBs_am(S0,K,r,T,sigmaBS, flag,  M_crr,D, U );

if flag==1 P_CM=P_CM_call; else P_CM=P_CM_put; end
    
method = ["CarrMadam", "Monte Carlo","fd","crr"]';
prices = [P_CM, P_mc, P_fd, P_crr]' ;
error = 100*(prices-P_mc)/P_mc;
american = table(method,prices,error)

%% Asian
flag = 1;
Strike = 0;
% Strike = K;
MC_Asian(Sbs, T, r, Strike, flag)

%% Lookback
Strike = 0;
% Strike = K;
MC_Lookback(Sbs, T, r, Strike, flag)

%% BarrierIn
flag=1;
% D=-inf; U= 1.1;
D=.9; U= inf;
P_in = MC_BarrierIn(Sbs, T, r, K, flag, D, U);
P_out = MC_European(Sbs, T, r, K, flag, D, U);

P_out+P_in
P_CM_call
%% With Dividend
flag = 1; T=.5;
r = .02; q = .01; r_q = r-q;
K=1;

[Sbsq_nAV, Sbs_AV] =BS_asset_AV(S0,T, r_q, sigmaBS,Nsim_mc,M_mc); Sbsq = [Sbsq_nAV; Sbs_AV];

P_bs_CRR_cam = crrBs_am(S0,K,r,T,sigmaBS, flag,  M_crr,D,U,q)
MC_European(Sbsq, T, r, K, flag)
MC_American(Sbsq, T, r, K, flag)
%%