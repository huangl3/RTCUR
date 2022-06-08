addpath(genpath(cd))
%% Generate RPCA problem with noise density alpha.
n = [500 500 500]; 
r = [3 3 3];
alpha = 0.3;
CI = 4;
[L,D] = rpca_tensor(n,r,alpha);
para.CI = CI; 
%% Resampled Fiber
ff = tic;
[L_core, X_sub_mat,~, ~] = RTCUR_ff(D, r, para);
L_out = ttm(L_core,X_sub_mat);
fprintf('Relative error for Fixed Fiber: %d\n', norm(L-L_out)/norm(L))
toc(ff)

%% Resampled Fiber_cpp
rf = tic;
% [L_core, X_sub_mat,~, ~] = RTCUR_rf(D, r, para);
% L_out = ttm(L_core,X_sub_mat);
[L_core, X_sub_mat,~, ~] = RTCUR_ff_cpp(D, r, para);
L_out = ttm(L_core,X_sub_mat);
fprintf('Relative error for Fixed Fiber CPP: %d\n', norm(L-L_out)/norm(L))
toc(rf)

% %% Fixed Chidori
% fc = tic;
% [L_core, X_sub_mat,~, ~] = RTCUR_fc(D, r, para);
% L_out = ttm(L_core,X_sub_mat);
% fprintf('Relative error for Fixed Chidori: %d\n', norm(L-L_out)/norm(L))
% toc(fc)
% 
% %% Resampled Chidori
% rc = tic;
% [L_core, X_sub_mat,runtime, errors] = RTCUR_rc(D, r, para);
% L_out = ttm(L_core,X_sub_mat);
% fprintf('Relative error for Resampled Chidori: %d\n', norm(L-L_out)/norm(L))
% toc(rc)
