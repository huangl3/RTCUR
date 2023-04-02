addpath(genpath(cd))
for r1 = [3 5 10]
%% Parameters
rng('default')
n1 = 300;
n = [n1 n1 n1]; r = [r1 r1 r1]; 
k_max = 10;
alpha_s = 0:0.02:0.6;
if r1 == 10
    CI_s = 0.2:0.2:1;
else
    CI_s = 0.2:0.2:2;
end
[ff_error, ff_time, rf_error, rf_time, fc_error, fc_time, rc_error, rc_time] = deal(zeros(length(alpha_s),length(CI_s),k_max));
%% k_max trails for each pair of {alpha,CI}
for i = 1:length(alpha_s)
    alpha = alpha_s(i);
    [L,D] = rpca_tensor(n,r,alpha);
    for j = 1:length(CI_s)
        para.CI = CI_s(j);
        round_tim = tic;
        
        func = @RTCUR_ff_morefiber;
        [tim,err] = test_k_round(func, D, L, r, para, k_max);
        ff_error(i,j,:) = err;
        ff_time(i,j,:) = tim;
        
        func = @RTCUR_rf_morefiber;
        [tim,err] = test_k_round(func, D, L, r, para, k_max);
        rf_error(i,j,:) = err;
        rf_time(i,j,:) = tim;
        
        func = @RTCUR_fc;
        [tim,err] = test_k_round(func, D, L, r, para, k_max);
        fc_error(i,j,:) = err;
        fc_time(i,j,:) = tim;
        
        func = @RTCUR_rc;
        [tim,err] = test_k_round(func, D, L, r, para, k_max);
        rc_error(i,j,:) = err;
        rc_time(i,j,:) = tim;

        fprintf('Current setting: alpha = %.2f, CI = %.1f, relative error: %.4e\n', alpha,CI_s(j),err(end))
        fprintf('Time elapsed for recent round: %.3f seconds.\n', toc(round_tim))
        fprintf('Progress for r = %i: %.3f%%\n', r1, 100*((i-1)*length(CI_s)+j)/length(alpha_s)/length(CI_s))
    end
end
save(sprintf('phase_n300_r%i.mat',r1),'fc_error','rc_error','ff_error','rf_error','CI_s','alpha_s')
end
function [times,errors] = test_k_round(func,D,L,r,para,k_max)
    times = zeros(1,k_max);
    errors = zeros(1,k_max);
    k = 1;
    while k <= k_max
        try
            [L_core, X_sub_mat, time, ~] = func(D, r, para);
            L_out = ttm(L_core,X_sub_mat);
            times(k) = sum(time(time>0));
            errorL = norm(L-L_out)/norm(L);
            errors(k) = errorL;
            
        catch
            fprintf('SVD not converge in this trail, resampling...\n')
        end
        k = k+1;
    end
end