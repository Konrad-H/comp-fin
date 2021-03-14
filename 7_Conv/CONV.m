function price = CONV(S0,T,r,K,  par,dist, flag, M, N,D,U,q )
if nargin<12
    q=0;
end
param.rf=r; %risk-free rate
param.q= q; %dividend
param.distr=dist; %what_dist
param.x = par; %parameters
param.T= T; %maturity
param.dt=param.T/M;
[S,v] = CONV0( S0, K, M, N, flag,D,U, param);
price=interp1(S,v,S0,'spline');

%%%%%%%%%%%%%%%%%
% DISTRIBUTIONS %
%%%%%%%%%%%%%%%%%
% 1 - BS
% 2 - Merton
% 3 - Kou

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S,v] = CONV0( S_0, K, M, N, flag, D,U, param)

b=2.5;
[x,~,~,H] = kernel(N,-b,b,param,0); % logprice grid x, conjugate of the characteristic function
S = S_0*exp(x);
% Payoff e trasformata di Fourier della densità
v = max(flag*(S-K),0).*(S>D).*(S<U); % Barrier
H=ifftshift(H);
for j = 1:M
    %v(S<=Barriera) = 0; % se in t=0 la barriera NON ha effetto
%     v=real((fft(ifft(v).*H)))*exp(-param.rf*param.dt);
    v=real(fftshift(fft(ifft(ifftshift(v)).*H)))*exp(-param.rf*param.dt);
    v(S<=D) = 0; % se in t=0 la barriera ha effetto
    v(S>=U) = 0;
end
index=find( (S>0.1*S_0).*(S<3*S_0));
S=S(index); v=v(index);
% figure
% plot(S,v)


