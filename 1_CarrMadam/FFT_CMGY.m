function [Price]=FFT_CM_Call(Strike,S0,T,r,x,flag)
% Price of a Plain Vanilla Call exploiting the Carr-Madan algorithm
% Model: CMGY

% discretization parameter
Npow=20; N=2^Npow; A=1200;

% v-> compute integral as a summation
eta=A/N; 
v=[0:eta:A*(N-1)/N]; v(1)=1e-22;

% lambda-> compute summation via FFT
lambda=2*pi/(N*eta); 
k=-lambda*N/2+lambda*(0:N-1);

% Fourier transform of z_k
CharFunc=@(v) exp(T*CharExp(v,x));

% temp=1/(kIG*sigmaIG^2); %alpha^2-beta^2 
% beta=thetaIG*kIG*temp; delta=1/(kIG*sqrt(temp)); alpha=sqrt(temp+beta^2);
% CharFunc=@(v) exp(T*CharExp2(v,alpha,beta,delta));

Z_k=exp(1i*r*v*T).*(CharFunc(v-1i)-1)./(1i*v.*(1i*v+1));
% Option Price
w=ones(1,N); w(1)=0.5; w(end)=0.5;
x=w.*eta.*Z_k.*exp(1i*pi*(0:N-1));
z_k=real(fft(x)/pi);
C=S0*(z_k+max(  flag*(1-exp(k-r*T))  ,0));
K=S0*exp(k);
% Output
index=find( K>0.1*S0 & K<3*S0 );
C=C(index); K=K(index);
% plot(K,C)
% title( 'Option Price' );
% xlabel('Strike');
Price=interp1(K,C,Strike,'spline');

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function V=CharExp(v,x)

sigmaIG=x(1); thetaIG=x(2); kIG=x(3); a = x(4);

% risk-neutral characteristic exponent
V=@(v) (1-a)/(a*kIG)*...
    (1-(1+ kIG*(v.^2*sigmaIG^2/2-1i*thetaIG*v)/(1-a)).^(a)); % without drift
drift_rn=-V(-1i); % Drift Risk_neutral
V=drift_rn*1i*v+V(v);

% function V=CharExp2(v,x)
%
% thetaIG=x(1); sigmaIG=x(2); kIG=x(3); 
% temp=1/(kIG*sigmaIG^2); %alpha^2-beta^2 
% beta=thetaIG*kIG*temp; delta=1/(kIG*sqrt(temp)); alpha=sqrt(temp+beta^2);
% 
% % risk-neutral characteristic exponent
% V=@(v) delta*sqrt(alpha^2-beta^2)-delta*sqrt(alpha^2-(beta+1i*v).^2); % without drift
% drift_rn=-V(-1i); % Drift Risk_neutral
% V=drift_rn*1i*v+V(v);





