clear; close all;

addpath('7_Conv_Method')
param.rf=0.05; %risk-free rate
param.q= 0; %dividend
param.distr=1; %BS distribution
x = [.1862];
param.x = x; %volatility

S_0=1; K=1; M=256*8; D=-inf; U=1.2; N=2^12;
param.T= 0.5; %maturity
param.dt=param.T/M;
flag=-1;
[S,v] = CONV( S_0, K, M, N, flag,D,U, param);
price=interp1(S,v,S_0,'spline')