clear all
%% Load Path
load("../9_Calibration/Calibration_Parameters.mat")
addpath("../1_CarrMadam")

%% Parameters
S0 = 1;
% S0 = 25.67;
T = .5; r = 0.05; K= 1;
models = ["Black & Scholes","Kou", "Merton","VG","EVG","NIG","ENIG","Heston","Bates"]' ;


%% Pricing
flag =-1;
P_bs = FFT_BS (K,S0,T,r,sigmaBS,flag);
P_fft_Kou = FFT_Kou(K,S0,T,r,pKou,flag);
P_fft_Mer = FFT_Merton(K,S0,T,r,pMer,flag);
P_fft_VG = FFT_VG(K,S0,T,r,pVG,flag);
P_fft_EVG = FFT_EVG(K,S0,T,r,pEVG,flag);
P_fft_NIG = FFT_NIG(K,S0,T,r,pNIG,flag);
P_fft_ENIG = FFT_ENIG(K,S0,T,r,pENIG,flag);
P_fft_Hes = FFT_Heston(K,S0,T,r,pHes,flag);
P_fft_Bat = FFT_Bates(K,S0,T,r,pBat,flag);

prices = [P_bs,...
    P_fft_Kou,P_fft_Mer,P_fft_VG,P_fft_EVG ,P_fft_NIG,P_fft_ENIG,...
    P_fft_Hes , P_fft_Bat]';
P_mean = mean(prices);
error = 100*(prices-P_mean)/P_mean;
table(models, prices, error)
