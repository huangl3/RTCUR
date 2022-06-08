% ff - fixed fiber
function [L_core,X_sub_mat,timer,err ] = RTCUR_ff_sym(D, R, para)
n_dim = size(D);
Sub_ten_ind = struct('type','()','subs',{[]});
n = ndims(D);
I = cell(n,1);
J = cell(n,1);
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
D_sub_mat = cell(n,1);
S_sub_mat = cell(n,1);
L_sub_mat = cell(n,1);
X_sub_mat = cell(n,1);
err    = -1*ones(max_iter,1);
timer  = -1*ones(max_iter,1);

%% Sample I,J, Initialize L,D fibers and core correspondingly, 
r_i = R(1);
d_i = n_dim(1);
l1 = min(ceil(CI*r_i*log(d_i)),d_i); % size of L_core
len1 = [l1 l1 l1];
l2 = min(ceil(10*CI*r_i*log(prod(n_dim)/d_i)),prod(n_dim)/d_i); % width of C_i(number of picked fibers)
len2 = [l2 l2 l2];
Ii = randperm(n_dim(1),len1(1)); I = {Ii,Ii,Ii};
Ji = randperm(prod(n_dim)/n_dim(1),len2(1)); J = {Ji,Ji,Ji};
for i_mod = 1:n
    L_sub_mat{i_mod} = zeros(d_i,len2(i_mod));
    L_core = zeros(len1);
end

% Extract fibers and core from D corresponding I and J.
norm_of_D = 0;
for i_mod = 1:n
   P = vec2ten_ind(J{i_mod},n_dim,i_mod);
   Sub_ten_ind.subs = P;  
   D_sub_mat{i_mod} = reshape(subsref(D,Sub_ten_ind),n_dim(i_mod),[]);
   norm_of_D = norm_of_D + norm(D_sub_mat{i_mod},'fro')^2;
end
Sub_ten_ind.subs = I;
D_core = subsref(D,Sub_ten_ind);
norm_of_D = norm_of_D + norm(D_core)^2;

%% Iteration
for i = 1:max_iter
tic
% Compute fibers and core for L^(k+1) (this usually faster than computing the whole L)
if i > 1
    for i_mod = 1:n
%         L_sub_mat{i_mod} = ttm_reduce(L_core,X_sub_mat,J,n_dim,i_mod);
        L_sub_mat{i_mod} = ttm_cpp(double(L_core),X_sub_mat{1},X_sub_mat{2},X_sub_mat{3},...
            J{i_mod},i_mod);
    end
    for i_mod = 1:n
        L_core = ttm(L_core,X_sub_mat{i_mod}(I{i_mod},:),i_mod);
    end
end

% Phase 1: update sparse component S
error = 0;
for i_mod = 1:n
    S_sub_mat{i_mod} = Thres(D_sub_mat{i_mod}-L_sub_mat{i_mod},zeta);
    error = error + norm(D_sub_mat{i_mod}-L_sub_mat{i_mod}-S_sub_mat{i_mod},'fro')^2;
end
S_core = tensor(Thres(double(D_core-L_core),zeta));
error = error+ norm(D_core-S_core-L_core)^2;

% phase 2: update low rank component L
L_core = D_core - S_core;

for i_mod = 1:n
    C = D_sub_mat{i_mod}-S_sub_mat{i_mod};
    U = C(I{i_mod},:);
    X_sub_mat{i_mod} = C*rinv(U,R(i_mod));
end

zeta = zeta*gamma;
timer(i) = toc;
err(i) = sqrt(error/norm_of_D);

if err(i) <= epsilon
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