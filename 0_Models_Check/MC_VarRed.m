%% Variance Reduction Asian
clc; clear; close all;
%% Parameters
S0 = 1;
% S0 = 25.67;
T = .5; r = 0.05; K= 1;

load("../9_Calibration/Calibration_Parameters.mat")
addpath("../2_BS_model")
addpath("../3_Levy_model")
addpath("../4_LV_model")
addpath("../5_MC_pricing")
%% MC parameters
Nsim = 1e5;
M = round(T*24); % Monitoring

Sgen = @(N)  BS_asset_AV(S0,T, r, sigmaBS,N,M);
len = @(CI) CI(2)-CI(1);

%% EUROPEAN
flag=1;
Nsm=1e6;
D=-inf; U=inf;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Normal MC
Sbs = Sgen(Nsm);
[~, CI1]  = MC_European(Sbs, T, r, K, flag, D, U);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Control Variable
Svr = Sgen(Nsm/10);
S = Sgen(Nsm);

ECV =  flag*S0;
gvr = flag*(Svr(:,end))*exp(-r*T);
g   = flag*(S(:,end))*exp(-r*T);

[~, ~, fvr] = MC_European(Svr, T, r, K, flag, D, U);
[~, ~, f] = MC_European(S, T, r, K, flag, D, U);
VC = cov(fvr, gvr);
alpha = -VC(1,2)/VC(2,2);

[~, CI_cv]  = MC_CV(fvr,gvr,f,g,ECV);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Antithetic Variable
[Sbs, Sbsav] =  Sgen(Nsm);
[~, ~,f1]  = MC_European(Sbs, T, r, K, flag, D, U);
[~, ~,f2]  = MC_European(Sbsav, T, r, K, flag, D, U);
[~,~, CI_av] = normfit ( (f1+f2)/2);

[[CI1; len(CI1)], [CI_cv; len(CI_cv)], [CI_av; len(CI_av)]]

%% ASIAN
flag=1;
Nsm=1e6;


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Normal MC
Sbs = Sgen(Nsm);
[~, CI1]  = MC_Asian(Sbs, T, r, 0, flag);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Control Variable
Svr = Sgen(Nsm/10);
S = Sgen(Nsm);

% FOR FLOATING STRIKE ASIAN
ECV =  flag*(S0-mean(S0*exp(r*T/M*(0:M)))*exp(-r*T));
gvr = flag*(Svr(:,end)-mean(Svr(:,:),2))*exp(-r*T);
g   = flag*(S(:,end)-mean(S(:,:),2))*exp(-r*T);

[~, ~, fvr] = MC_Asian(Svr, T, r, 0, flag);
[~, ~, f] = MC_Asian(S, T, r, 0, flag);

[~, CI_cv]  = MC_CV(fvr,gvr,f,g,ECV);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Antithetic Variable
[Sbs, Sbsav] =  Sgen(Nsm);
[~, ~,f1]  = MC_Asian(Sbs, T, r, 0, flag);
[~, ~,f2]  = MC_Asian(Sbsav, T, r, 0, flag);
[~,~, CI_av] = normfit ( (f1+f2)/2);
[[CI1; len(CI1)], [CI_cv; len(CI_cv)], [CI_av; len(CI_av)]]
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% American
flag=-1;
Nsm=1e4;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Normal MC
Sbs = Sgen(Nsm);
[~, CI1]  = MC_American(Sbs, T, r, K, flag);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Control Variable
Svr = Sgen(Nsm/10);
S = Sgen(Nsm);

% ECV =  flag*(S0-mean(S0*exp(r*T/M*(0:M)))*exp(-r*T));
% gvr = flag*(Svr(:,end)-mean(Svr(:,:),2))*exp(-r*T);
% g   = flag*(S(:,end)-mean(S(:,:),2))*exp(-r*T);

% This is a much better variable
ECV =  flag*S0;
gvr = flag*(Svr(:,end))*exp(-r*T);
g   = flag*(S(:,end))*exp(-r*T);

[~, ~, fvr] = MC_American(Svr, T, r, K, flag);
[~, ~, f] = MC_American(S, T, r, K, flag);

[~, CI_cv]  = MC_CV(fvr,gvr,f,g,ECV);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Antithetic Variable
[Sbs, Sbsav] =  Sgen(Nsm);
[~, ~,f1]  = MC_American(Sbs, T, r, K, flag);
[~, ~,f2]  = MC_American(Sbsav, T, r, K, flag);
[~,~, CI_av] = normfit ( (f1+f2)/2);

Sbs = Sgen(Nsm*10);
[~, CI10]  = MC_American(Sbs, T, r, K, flag);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% MANY SIMULATIONS
[[CI1; len(CI1)], [CI_cv; len(CI_cv)], [CI_av; len(CI_av)], [CI10; len(CI10)]]
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% LOOKBACK
flag=1;
Nsm=1e5;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Normal MC
Sbs = Sgen(Nsm);
[~, CI1]  = MC_Lookback(Sbs, T, r, 0, flag);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Control Variable
Svr = Sgen(Nsm/10);
S = Sgen(Nsm);

% FOR FLOATING STRIKE ASIAN

ECV =  flag*(S0-mean(S0*exp(r*T/M*(0:M)))*exp(-r*T));
gvr = flag*(Svr(:,end)-mean(Svr(:,:),2))*exp(-r*T);
g   = flag*(S(:,end)-mean(S(:,:),2))*exp(-r*T);

[~, ~, fvr] = MC_Lookback(Svr, T, r, 0, flag);
[~, ~, f] = MC_Lookback(S, T, r, 0, flag);

[~, CI_cv]  = MC_CV(fvr,gvr,f,g,ECV);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Antithetic Variable
[Sbs, Sbsav] =  Sgen(Nsm);
[~, ~,f1]  = MC_Lookback(Sbs, T, r, 0, flag);
[~, ~,f2]  = MC_Lookback(Sbsav, T, r, 0, flag);
[~,~, CI_av] = normfit ( (f1+f2)/2);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% MANY SIMULATIONS
Sbs = Sgen(Nsm*10);
[~, CI10]  = MC_American(Sbs, T, r, K, flag);

[[CI1; len(CI1)], [CI_cv; len(CI_cv)], [CI_av; len(CI_av)], [CI10; len(CI10)]]
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
