function [Price]=CM_Call_VG(Strike,S0,T,r,x, flag)
% Price of a Plain Vanilla Call exploiting the Carr-Madan algorithm
% Model: VG

% parameter
sigmaVG=x(1); thetaVG=x(2); kVG=x(3);

% discretization parameter
Npow=15; N=2^Npow; A=1200;

% v-> compute integral as a summation
eta=A/N;
v=[0:eta:A*(N-1)/N]; v(1)=1e-22;

% lambda-> compute summation via FFT
lambda=2*pi/(N*eta);
k=-lambda*N/2+lambda*(0:N-1);

% Fourier transform of z_k
CharFunc=@(v) exp(T*CharExp(v,sigmaVG,thetaVG,kVG));
Z_k=exp(1i*r*v*T).*(CharFunc(v-1i)-1)./(1i*v.*(1i*v+1));
% Option Price
w=ones(1,N); w(1)=0.5; w(end)=0.5;
x=w.*eta.*Z_k.*exp(1i*pi*(0:N-1));
z_k=real(fft(x)/pi);
C=S0*(z_k+max( flag*(1-exp(k-r*T)),0));
K=S0*exp(k);

% Output
index=find( K>0.1*S0 & K<3*S0 );
C=C(index); K=K(index);
% plot(K,C)
% title( 'Option Price' );
% xlabel('Strike');
Price=interp1(K,C,Strike,'spline');
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function V=CharExp(v,sigma,theta,k)
% risk-neutral characteristic exponent
V=@(v)-log(1+v.^2*sigma^2*k/2-1i*theta*k*v)/k; % without drift
drift_rn=-V(-1i); % Drift Risk_neutral
V=drift_rn*1i*v+V(v);
end















