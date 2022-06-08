% fc - fixed chidori
function [L_core,X_sub_mat,timer,err ] = RTCUR_rc(D, r, para)
d = size(D);
Sub_ten_ind = struct('type','()','subs',{[]});
n = ndims(D);  % number of modes
I = cell(n,1);
%% Default/Inputed parameters
max_iter  = 200;
epsilon   = 1e-6;
zeta      = 1.1*max(D(:));
gamma     = 0.7;    
CI        = 2;
if isfield(para,'max_iter') 
    max_iter = para.max_iter; 
end
if isfield(para,'epsilon') 
    epsilon = para.epsilon; 
end
if isfield(para,'zeta') 
    zeta = para.zeta; 
end
if isfield(para,'gamma') 
    gamma = para.gamma; 
end
if isfield(para,'CI') 
    CI = para.CI; 
end
D_sub_mat = cell(n,1); % tensor
S_sub_mat = cell(n,1); % tensor
L_sub_mat = cell(n,1); % tensor
X_sub_mat = cell(n,1); % double / matrix
err    = -1*ones(max_iter,1);
timer  = -1*ones(max_iter,1);

%% Initialize L core and chidori（to be zeros), 
len = zeros(1,n); % size of L_core
for i = 1:n
    len(i) = min(ceil(CI*r(i)*log(d(i))),d(i));
end
for i = 1:n
    shapeL = len; shapeL(i) = d(i);
    L_sub_mat{i} = tensor(zeros(shapeL));
    L_core = tensor(zeros(len));
end

%% Iteration
for it = 1:max_iter
tic
% Sample I
for i = 1:n
    I{i} = randperm(d(i),len(i));
end
% Extract core and chidori from D corresponding I
norm_of_D = 0;
Sub_ten_ind.subs = I;
D_core = subsref(D,Sub_ten_ind);
norm_of_D = norm_of_D + norm(D_core)^2;
for i = 1:n
    J = I;
    J{i} = 1:d(i);
    Sub_ten_ind.subs = J;
    D_sub_mat{i} = subsref(D,Sub_ten_ind);
    norm_of_D = norm_of_D + norm(D_sub_mat{i})^2;
end

% Compute chidori and core for L^(k+1) (this usually faster than computing the whole L)
if it > 1
    for i = 1:n
        L_sub_mat{i} = L_core;
        for j = 1:n
            if j == i
                L_sub_mat{i} = ttm(L_sub_mat{i},X_sub_mat{j},j);
            else
                L_sub_mat{i} = ttm(L_sub_mat{i},X_sub_mat{j}(I{j},:),j);
            end
        end
    end
    for i = 1:n
        L_core = ttm(L_core,X_sub_mat{i}(I{i},:),i);
    end
end

% Phase 1: update sparse component S
error = 0;
for i = 1:n
    S_sub_mat{i} = tensor(Thres(double(D_sub_mat{i}-L_sub_mat{i}),zeta));
    error = error + norm(D_sub_mat{i}-L_sub_mat{i}-S_sub_mat{i})^2;
end
S_core = tensor(Thres(double(D_core-L_core),zeta));
error = error+ norm(D_core-S_core-L_core)^2;

% phase 2: update low rank component L
L_core = D_core - S_core;

for i = 1:n
    C = double(tenmat(D_sub_mat{i}-S_sub_mat{i},i));
    U = C(I{i},:);
    X_sub_mat{i} = C*rinv(U,r(i));
end

zeta = zeta*gamma;
timer(it) = toc;
err(it) = sqrt(error/norm_of_D);

if err(it) <= epsilon
    if isfield(para,'info') 
        fprintf('Total %d iteration, final error: %e, total time: %f  \n', i, err(i), sum(timer(timer>0)));
    end
    return
end
end
end

function TX = Thres(X,zeta)
TX = X.*(abs(X)>zeta);
end