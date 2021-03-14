function [price, CI] = MC_CV(f_vr, g_vr, f, g, ECV)
    % Compute variance reduction MC with control var g over payoff f
    alpha = CV(f_vr, g_vr);
    [price, ~,CI]= normfit(f+alpha*(g-ECV)); 



function alpha = CV(fvr, gvr)
    % f is the payoff
    % g is control variable
    VC = cov(fvr, gvr);
    alpha = -VC(1,2)/VC(2,2);