function [Price,CI, DiscPayoff] = MC_American( St, T, Rate,Strike, flag)

steps = size(St,2);
Nsim = size(St,1);
dt = T/steps;

%% Longstaff-Schwartz
ExerciseTime=steps*ones(Nsim,1);
Cashflow=max( flag*( St(:,end)-Strike),0 );
for step=steps-1:-1:1
    if flag==1
        InMoney=find( St(:,step+1)>Strike );    
    elseif flag==-1
        InMoney=find( St(:,step+1)<Strike );
    end
    Stemp=St(InMoney,step+1);
    %-------- Regression [1, S, S.^2]
    RegressionMatrix=[ones(size(Stemp)), Stemp, Stemp.^2];
    YData=Cashflow(InMoney).*exp( -Rate*dt*(ExerciseTime(InMoney)-step));
    alpha=RegressionMatrix\YData;
    %----------------------------------------------------------------------
    CV=RegressionMatrix*alpha;
    IV=flag*(Stemp-Strike);
    %-------- Early Exercise
    EarlyExercise_inMoney_index=find( IV> CV );
    EarlyExercise_index=InMoney(EarlyExercise_inMoney_index);
    Cashflow(EarlyExercise_index)=IV(EarlyExercise_inMoney_index);
    ExerciseTime(EarlyExercise_index)=step;
end
[Price,~,CI]=normfit( Cashflow.*exp(-Rate*dt*ExerciseTime) );

DiscPayoff = Cashflow.*exp(-Rate*dt*ExerciseTime);


end