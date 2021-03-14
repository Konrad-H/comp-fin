clear all; close all;
%% PLAIN VANILLA CALL calibration
% [Strike/Spot,Time To Maturity,Price]
data=[0.5  0.3  13.0391
      0.5  0.6  13.2413
      0.5  0.9  13.4410
      0.5  1    13.5070
      0.5  1.5  13.8325
      0.5  2    14.1501
      0.6  0.6  10.7861
      0.7  0.6  8.36675
      0.8  0.6  6.01538
      0.9  0.6  3.79936
      1    0.6  1.87872
      1.1  0.6  0.61598
      1.3  0.6  0.05320
      1.4  0.6  0.01986
      1.5  0.6  0.00832];
spot=25.67;   rf=0.05;
strike=spot*data(:,1); 
maturity=data(:,2);
PriceMkt=data(:,3); %market price
options=optimoptions('lsqnonlin','FunctionTolerance',1e-8,'OptimalityTolerance',1e-8,'StepTolerance',1e-8);

%% Calibrating the BS

sigmaBS=lsqnonlin(@(x) funBS(x,spot,strike,rf,maturity,PriceMkt),0.3,[0],[0.5],options)
[Error,PriceBS]=funBS(sigmaBS,spot,strike,rf,maturity,PriceMkt);
ErrorBS=sum(Error.^2)


%% Calibrating the KOU

pKou=lsqnonlin(@(x) funKou(x,spot,strike,rf,maturity,PriceMkt),[0.3 0.5 3 10 10],[0 0 0 1 1], [0.5 1 100 1000 1000],options)
[Error,PriceKou]=funKou(pKou,spot,strike,rf,maturity,PriceMkt);
ErrorKou=sum(Error.^2)


%% Calibrating the MERTON

pMer=lsqnonlin(@(x) funMerton(x,spot,strike,rf,maturity,PriceMkt),[0.3 3 1 .3 ],[0 0 -1 0 ], [0.5 100 100 1],options)
[Error,PriceMer]=funMerton(pMer,spot,strike,rf,maturity,PriceMkt);
ErrorMer=sum(Error.^2)

%% Calibrating the VG
pVG=lsqnonlin(@(x) funVG(x,spot,strike,rf,maturity,PriceMkt),[0.1 0 0.3],[0 -5 0], [1 5 5],options)
[Error,PriceVG]=funVG(pVG,spot,strike,rf,maturity,PriceMkt);
ErrorVG=sum(Error.^2)


%% Calibrating the EVG
pEVG=lsqnonlin(@(x) funEVG(x,spot,strike,rf,maturity,PriceMkt),[0.3 0 0.3 0.01],[0 -5 0 0], [1 5 5 1],options)
[Error,PriceEVG]=funEVG(pEVG,spot,strike,rf,maturity,PriceMkt);
ErrorEVG=sum(Error.^2)


%% Calibrating the NIG

pNIG=lsqnonlin(@(x) funNIG(x,spot,strike,rf,maturity,PriceMkt),[0.3 0 0.3],[0 -5 0], [1 5 5],options)
[Error,PriceNIG]=funNIG(pNIG,spot,strike,rf,maturity,PriceMkt);
ErrorNIG=sum(Error.^2)

%% Calibrating the ENIG

pENIG=lsqnonlin(@(x) funENIG(x,spot,strike,rf,maturity,PriceMkt),[0.3 0 0.3 0.1],[0 -5 0 0], [1 5 5 1],options)
[Error,PriceENIG]=funENIG(pENIG,spot,strike,rf,maturity,PriceMkt);
ErrorENIG=sum(Error.^2)

%% Calibrating the Heston
pHes=lsqnonlin(@(x) funHeston(x,spot,strike,rf,maturity,PriceMkt),[0.2 0.01 -0.1 0.1 0.1],[0 0 -1 0 0], [ 10 100 1 1 100],options)
[Error,PriceHes]=funHeston(pHes,spot,strike,rf,maturity,PriceMkt);
ErrorHeston=sum(Error.^2)

%% Calibrating the Bates
pBat=lsqnonlin(@(x) funBates(x,spot,strike,rf,maturity,PriceMkt),[0.2 0.01 -0.1 0.1 0.1 3 1 .3],[0 0 -1 0 0 0 -1 0], [ 10 100 1 1 100 100 100 1],options)
[Error,PriceBat]=funBates(pBat,spot,strike,rf,maturity,PriceMkt);
ErrorBates=sum(Error.^2)
%% Calibrating the CMGY
pCMGY=lsqnonlin(@(x) funCMGY(x,spot,strike,rf,maturity,PriceMkt),[0.3 0 0.3 .9],[0 -5 0 1e-8], [1 5 5 2],options)
[Error,PriceCMGY]=funCMGY(pCMGY,spot,strike,rf,maturity,PriceMkt);
ErrorCMGY=sum(Error.^2)

%%

save('./Params/pBS','sigmaBS')
save('./Params/pKou','pKou')
save('./Params/pMer','pMer')
save('./Params/pVG','pVG')
save('./Params/pEVG','pEVG')
save('./Params/pNIG','pNIG')
save('./Params/pENIG','pENIG')
save('./Params/pHes','pHes')
save('./Params/pBat','pBat')
save('./Params/pBat','pCMGY')
%% save
save('Calibration_Parameters','sigmaBS','pKou','pMer','pVG','pNIG','pEVG','pENIG','pHes','pBat','pCMGY')

%% load
load('Calibration_Parameters')

[~,PriceBS]=funBS(sigmaBS,spot,strike,rf,maturity,PriceMkt);
[~,PriceKou]=funKou(pKou,spot,strike,rf,maturity,PriceMkt);
[~,PriceMer]=funMerton(pMer,spot,strike,rf,maturity,PriceMkt);
[~,PriceVG]=funVG(pVG,spot,strike,rf,maturity,PriceMkt);
[~,PriceEVG]=funEVG(pEVG,spot,strike,rf,maturity,PriceMkt);
[~,PriceNIG]=funNIG(pNIG,spot,strike,rf,maturity,PriceMkt);
[~,PriceENIG]=funENIG(pENIG,spot,strike,rf,maturity,PriceMkt);
[~,PriceHes]=funHeston(pHes,spot,strike,rf,maturity,PriceMkt);
[~,PriceBat]=funBates(pBat,spot,strike,rf,maturity,PriceMkt);
[~,PriceCMGY]=funBates(pBat,spot,strike,rf,maturity,PriceMkt);

%% PLot
figure
hold on


plot(PriceBS,'d')

plot(PriceKou,'s')

plot(PriceVG,'<')

plot(PriceMer, '*')
plot(PriceNIG, 'o')
plot(PriceHes, '.')
plot(PriceBat, '.')
plot(PriceMkt,'+k')

legend('Calibrated Price (BS)','Calibrated Price (Kou)',...
    'Calibrated Price (VG)','Calibrated Price (Merton)',...
    'Calibrated Price (NIG)', 'Calibrated Price (Heston)',...
    'Calibrated Price (Bates)',...
    'Market Price')
hold off

%% Implied Volatility
options = optimoptions('fsolve','OptimalityTolerance',1e-8,'FunctionTolerance',1e-9,'StepTolerance',1e-8,'MaxFunctionEvaluations',400);
% wrt Strike
T=0.6;
id=find(maturity==T);
sigma_mkt=zeros(size(id)); sigma_BS=zeros(size(id)); sigma_Kou=zeros(size(id));  
sigma_Mer=zeros(size(id));sigma_VG=zeros(size(id));sigma_EVG=zeros(size(id));
sigma_NIG=zeros(size(id));sigma_ENIG=zeros(size(id));sigma_Hes=zeros(size(id));
sigma_Bat=zeros(size(id));
for i=1:length(id)
    sigma_mkt(i)=fsolve(@(s) blsprice(spot, strike(id(i)), rf, T, s)-PriceMkt(id(i))+100*(s>=0.5)*(s<=0), 0.5,options);
    sigma_BS(i) =fsolve(@(s) blsprice(spot, strike(id(i)), rf, T, s)-PriceBS(id(i))+100*(s>=0.5)*(s<=0), sigmaBS,options);
    sigma_Kou(i)=fsolve(@(s) blsprice(spot, strike(id(i)), rf, T, s)-PriceKou(id(i))+100*(s>=0.5)*(s<=0), 0.5,options);
    sigma_Mer(i) =fsolve(@(s) blsprice(spot, strike(id(i)), rf, T, s)-PriceMer(id(i))+100*(s>=0.5)*(s<=0), 0.5,options);    
    sigma_VG(i) =fsolve(@(s) blsprice(spot, strike(id(i)), rf, T, s)-PriceVG(id(i))+100*(s>=0.5)*(s<=0), 0.5,options);
    sigma_EVG(i) =fsolve(@(s) blsprice(spot, strike(id(i)), rf, T, s)-PriceEVG(id(i))+100*(s>=0.5)*(s<=0), 0.5,options);
    sigma_NIG(i) =fsolve(@(s) blsprice(spot, strike(id(i)), rf, T, s)-PriceNIG(id(i))+100*(s>=0.5)*(s<=0), 0.5,options);
    sigma_ENIG(i) =fsolve(@(s) blsprice(spot, strike(id(i)), rf, T, s)-PriceENIG(id(i))+100*(s>=0.5)*(s<=0), 0.5,options);
    sigma_Hes(i) =fsolve(@(s) blsprice(spot, strike(id(i)), rf, T, s)-PriceHes(id(i))+100*(s>=0.5)*(s<=0), 0.5,options);
    sigma_Bat(i) =fsolve(@(s) blsprice(spot, strike(id(i)), rf, T, s)-PriceBat(id(i))+100*(s>=0.5)*(s<=0), 0.5,options);
end

figure
hold on
plot(strike(id),sigma_mkt,'ok-')
plot(strike(id),sigma_BS,'+-')
plot(strike(id),sigma_Kou,'+-')
plot(strike(id),sigma_Mer,'+-')
plot(strike(id),sigma_VG, '+-')
plot(strike(id),sigma_EVG,'+-')
plot(strike(id),sigma_NIG,'+-')
plot(strike(id),sigma_ENIG,'+-')
plot(strike(id),sigma_Hes,'+-')
plot(strike(id),sigma_Bat,'+-')
title("Implied Volatility vs Strike")
legend('Market IV','BS IV','KOU IV','Merton IV','VG IV', 'EVG IV', ...
    'NIG IV', 'ENIG IV', 'Heston IV','Bates IV');
xlabel('Strike'); ylabel('Volatility');
hold off
% wrt Maturity
K=0.5*spot;
id=find(strike==K);
sigma_mkt=zeros(size(id)); sigma_BS=zeros(size(id)); sigma_Kou=zeros(size(id));  
sigma_Mer=zeros(size(id));sigma_VG=zeros(size(id));sigma_EVG=zeros(size(id));
sigma_NIG=zeros(size(id));sigma_ENIG=zeros(size(id));sigma_Hes=zeros(size(id));
sigma_Bat=zeros(size(id));
for i=1:length(id)
    sigma_mkt(i)=fsolve(@(s) blsprice(spot, K, rf, maturity(id(i)), s)-PriceMkt(id(i))+100*(s>=0.5)*(s<=0), 0.5,options);
    sigma_BS(i) =fsolve(@(s) blsprice(spot, K, rf, maturity(id(i)), s)-PriceBS(id(i))+100*(s>=0.5)*(s<=0), sigmaBS,options);
    sigma_Kou(i)=fsolve(@(s) blsprice(spot, K, rf, maturity(id(i)), s)-PriceKou(id(i))+100*(s>=0.5)*(s<=0), 0.5,options);
 sigma_Mer(i)=fsolve(@(s) blsprice(spot, K, rf, maturity(id(i)), s)-PriceMer(id(i))+100*(s>=0.5)*(s<=0), 0.5,options);    
sigma_VG(i) =fsolve(@(s) blsprice(spot, K, rf, maturity(id(i)), s)-PriceVG(id(i))+100*(s>=0.5)*(s<=0), 0.5,options);
sigma_EVG(i) =fsolve(@(s) blsprice(spot, K, rf, maturity(id(i)), s)-PriceEVG(id(i))+100*(s>=0.5)*(s<=0), 0.5,options);
sigma_NIG(i) =fsolve(@(s) blsprice(spot, K, rf, maturity(id(i)), s)-PriceNIG(id(i))+100*(s>=0.5)*(s<=0), 0.5,options);
sigma_ENIG(i) =fsolve(@(s) blsprice(spot, K, rf, maturity(id(i)), s)-PriceENIG(id(i))+100*(s>=0.5)*(s<=0), 0.5,options);
sigma_Hes(i) =fsolve(@(s) blsprice(spot, K, rf, maturity(id(i)), s)-PriceHes(id(i))+100*(s>=0.5)*(s<=0), 0.5,options);
sigma_Bat(i) =fsolve(@(s) blsprice(spot, K, rf, maturity(id(i)), s)-PriceBat(id(i))+100*(s>=0.5)*(s<=0), 0.5,options);

end
figure
hold on
plot(maturity(id),sigma_mkt,'ok-'); 

plot(maturity(id),sigma_BS,'+-')
plot(maturity(id),sigma_Kou,'+-')
plot(maturity(id),sigma_Mer,'+-')
plot(maturity(id),sigma_VG,'+-')
plot(maturity(id),sigma_EVG,'+-')
plot(maturity(id),sigma_NIG,'+-')
plot(maturity(id),sigma_ENIG,'+-')
plot(maturity(id),sigma_Hes,'+-')
plot(maturity(id),sigma_Bat,'+-')
title(" Implied Volatility vs TTM")
legend('Market IV','BS IV','KOU IV','Merton IV','VG IV', 'EVG IV', ...
    'NIG IV', 'ENIG IV', 'Heston IV','Bates IV'); 
xlabel('Time to Maturity'); ylabel('Volatility');
hold off

%% Plot
figure;
hold on
plot(maturity(id),sigma_mkt,'ok-'); 
plot(maturity(id),sigma_BS,'+-');
plot(maturity(id),sigma_Kou,'+-');
plot(maturity(id),sigma_Hes,'+-');
hold off