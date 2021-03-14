clc; clear; close all;
%% Parameters
S0 = 1;
T = .5; r = 0.05; K= 1;

load("../9_Calibration/Calibration_Parameters.mat")
addpath("../2_BS_model")
addpath("../3_Levy_model")
addpath("../4_LV_model")
addpath("../5_MC_pricing")
%% MC parameters
Nsim = 1e5;
M = round(T*52); % Monitoring
%% Simulate MonteCarlo Assets

% BS Martingale Check
[Sbs, Sbs_AV] =BS_asset_AV(S0,T, r, sigmaBS,Nsim,M);
[Check,~,IC]=normfit( Sbs(:,end)*exp(-r*T) -S0)

% Matrix Kou Martingale Check
[Skou, Skou_AV] = asset_Kou(S0,T, r, pKou,Nsim,M);
[Check,~,IC]=normfit( Skou(:,end)*exp(-r*T) -S0)

% Matrix Merton Martingale Check
[Smer, Smer_AV]=asset_Merton(S0,T, r, pMer,Nsim,M);
[Check,~,IC]=normfit( Smer(:,end)*exp(-r*T) -S0)

% VG Martingale Check
[SVG, SVG_AV]=asset_VG(S0,T, r, pVG,Nsim,M);
[Check,~,IC]=normfit( SVG(:,end)*exp(-r*T) -S0)

% NIG Martinagale Check
[SNIG, SNIG_AV]=asset_NIG(S0,T, r, pNIG,Nsim,M);
[Check,~,IC]=normfit( SNIG(:,end)*exp(-r*T) -S0) 

% ENIG Martinagale Check
[SENIG, SENIG_AV]=asset_ENIG(S0,T, r, pENIG,Nsim,M);
[Check,~,IC]=normfit( SENIG(:,end)*exp(-r*T) -S0) 

% Heston QE Martingale Check
Shes = asset_Heston(S0,T, r, pHes,Nsim,M); %doesnt perform well
[Check,~,IC]=normfit( Shes(:,end)*exp(-r*T) -S0) 

% % Heston Euler Martingale Check
% Shes2 = asset_Heston2(S0,T, r, pHes,Nsim,M);
% [Check,~,IC]=normfit( Shes2(:,end)*exp(-r*T) -S0)

% Bates Martingale Check
Sbat = asset_Bates(S0,T, r, pBat,Nsim,M);
[Check,~,IC]=normfit( Sbat(:,end)*exp(-r*T) -S0)

% % Bates Martingale Check
% Sbat2 = asset_Bates2(S0,T, r, pBat,Nsim,M);
% [Check,~,IC]=normfit( Sbat2(:,end)*exp(-r*T) -S0)

%% Plotting
figure;
hold on
plot(mean(Sbs,1))
plot(mean(Skou,1))
plot(mean(Smer,1))
plot(mean(SVG,1))
plot(mean(SNIG,1))
plot(mean(SENIG,1))
plot(mean(Shes,1))
% plot(mean(Shes2,1))
plot(mean(Sbat,1))
models = ["Black & Scholes","Kou", "Merton","VG","NIG","ENIG","Heston","Bates"]' ;
legend(models)

%% Price European
flag=-1;  D=-inf;U=inf;
[P_sbs, ~]  = MC_European(Sbs, T, r, K, flag,D,U);
[P_kou, ~]  = MC_European(Skou, T, r, K, flag,D,U);
[P_mer, ~]  = MC_European(Smer, T, r, K, flag,D,U);
[P_vg, ~]  = MC_European(SVG, T, r, K, flag,D,U);
[P_nig, ~]  = MC_European(SNIG, T, r, K, flag,D,U);
[P_enig, ~]  = MC_European(SENIG, T, r, K, flag,D,U);
[P_hes, ~]  = MC_European(Shes, T, r, K, flag,D,U);
[P_bat, ~]  = MC_European(Sbat, T, r, K, flag,D,U);

prices = [P_sbs,P_kou ,P_mer , P_vg, P_nig,P_enig P_hes, P_bat]';
P_mean = mean(prices);
error = 100*(prices-P_mean)/P_mean;
table(models,prices,error )
%% Price American
flag=-1; 
[P_sbs, ~]  = MC_American(Sbs, T, r, K, flag);
[P_kou, ~]  = MC_American(Skou, T, r, K, flag);
[P_mer, ~]  = MC_American(Smer, T, r, K, flag);
[P_vg, ~]  = MC_American(SVG, T, r, K, flag);
[P_nig, ~]  = MC_American(SNIG, T, r, K, flag);
[P_enig, ~]  = MC_American(SENIG, T, r, K, flag);
[P_hes, ~]  = MC_American(Shes, T, r, K, flag);
[P_bat, ~]  = MC_American(Sbat, T, r, K, flag);

prices = [P_sbs,P_kou ,P_mer , P_vg, P_nig,P_enig P_hes, P_bat]';
error = 100*(prices-P_nig)/P_nig;
table(models,prices,error )

%% Price Asian
flag = 1;
Strike = 0;
% Strike = K;
MC_Asian(Sbs, T, r, Strike, flag)

[P_sbs, ~]  = MC_Asian(Sbs, T, r, Strike, flag);
[P_kou, ~]  = MC_Asian(Skou, T, r, Strike, flag);
[P_mer, ~]  = MC_Asian(Smer, T, r, Strike, flag);
[P_vg, ~]  = MC_Asian(SVG, T, r, Strike, flag);
[P_nig, ~]  = MC_Asian(SNIG, T, r, Strike, flag);
[P_enig, ~]  = MC_Asian(SENIG, T, r, K, flag);
[P_hes, ~]  = MC_Asian(Shes, T, r, Strike, flag);
[P_bat, ~]  = MC_Asian(Sbat, T, r, Strike, flag);

prices = [P_sbs,P_kou ,P_mer , P_vg, P_nig,P_enig P_hes, P_bat]';
error = 100*(prices-P_nig)/P_nig;
table(models,prices,error )

%% Price Lookback
Strike = 0; flag=-1;
% Strike = K;
[P_sbs, ~]  = MC_Lookback(Sbs, T, r, Strike, flag);
[P_kou, ~]  = MC_Lookback(Skou, T, r, Strike, flag);
[P_mer, ~]  = MC_Lookback(Smer, T, r, Strike, flag);
[P_vg, ~]  = MC_Lookback(SVG, T, r, Strike, flag);
[P_nig, ~]  = MC_Lookback(SNIG, T, r, Strike, flag);
[P_enig, ~]  = MC_Lookback(SENIG, T, r, K, flag);
[P_hes, ~]  = MC_Lookback(Shes, T, r, Strike, flag);
[P_bat, ~]  = MC_Lookback(Sbat, T, r, Strike, flag);

prices = [P_sbs,P_kou ,P_mer , P_vg, P_nig,P_enig P_hes, P_bat]';
error = 100*(prices-P_nig)/P_nig;
table(models,prices,error )

%% Variance Reduction European/DownOut
Sbs = BS_asset_AV(S0,T, r, sigmaBS,Nsim,M);
[P_sbs1, CI1]  = MC_European(Sbs, T, r, K, flag,D,U)
Svr = Sbs(1:Nsim/3,:);
S = Sbs((1+Nsim/3):end , :);

gvr = flag*Svr(:,end)*exp(-r*T);
g   = flag*S(:,end)*exp(-r*T);
ECV = flag*S0;

[~, ~, fvr] = MC_European(Svr, T, r, K, flag, D,U);
[~, ~, f] = MC_European(S, T, r, K, flag, D,U);

[P_sbs2, CI2]  = MC_CV(fvr,gvr,f,g,ECV)

FFT_BS(K,S0,T,r,sigmaBS, flag)

