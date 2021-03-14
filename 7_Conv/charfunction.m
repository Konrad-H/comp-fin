function F = charfunction(u,parameters,flag)
% flag=0 --> funzione caratteristica per problema backward 
% flag=1 --> funzione caratteristica per problema forward

if nargin==2
    flag=0;
end

drift_r = (parameters.rf-parameters.q)*parameters.dt;
drift_rn = -log(charfunction0(-1i,parameters));

meancorrection = drift_r  + drift_rn;

F = exp(1i*meancorrection*u).*charfunction0(u,parameters);

if flag==0
    F=conj(F);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = charfunction0(u,parameters)
dt=parameters.dt;
x = parameters.x;       
switch parameters.distr

    case 1 % Normal
%         m = parameters.m;
        s = x(1);
	
	% Rearrange parameters (time rescaling)
        m = 0; %         m = m*dt;
        s = s*sqrt(dt);

        F = exp(1i*u*m-0.5*(s*u).^2);

    case 2 % Normal inverse Gaussian (NIG) 
        alpha = x(1);
        beta = x(2);
        delta = x(3);

        % Rearrange parameters (time rescaling)
        Psi = -delta*(sqrt(alpha^2-(beta+1i*u).^2)-sqrt(alpha^2-beta^2));

        F = exp(dt*Psi);
        
    case 3 % Normal inverse Gaussian (NIG) 
        thetaIG = x(1);
        sigmaIG = x(2);
        kIG = x(3);

        % Rearrange parameters (time rescaling)
        thetaIG = thetaIG*dt;
        sigmaIG = sigmaIG*sqrt(dt);
        kIG = kIG/dt;

        % risk-neutral characteristic exponent
        F= exp( 1/kIG-sqrt(1+u.^2*sigmaIG^2*kIG-2i*thetaIG*u*kIG)/kIG ); % without drift
        
    case 4 % Kou model
        sigma=x(1); 
        p=x(2); 
        lambdaK=x(3); 
        lambdap=x(4); 
        lambdam=x(5);
        
        % time scale ??

        Psi = -sigma^2/2*u.^2+1i*u*lambdaK.* ...
            (p./(lambdap-1i*u)-(1-p)./(lambdam+1i*u));
        
        F=  exp( dt*Psi  ); % without drift
end
