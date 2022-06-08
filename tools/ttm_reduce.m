% compute the production of T = ttensor(L_core,X_sub_mat) but only take |J|
% fibers along mode 1,...,i-1,i+1,...,n. (i = it_mod, n = n_dim)
function C = ttm_reduce(core,X_sub_mat,J,n_dim,it_mod)
N = length(J);
s = length(J{it_mod});
p = vecJ2ten_ind(J{it_mod},n_dim,it_mod);
C = zeros(size(core,it_mod),s);
for i = 1:s
    ttm_tube = core;
    for j = setdiff(1:N,it_mod)
        X_tube = X_sub_mat{j}(p(i,j),:);
        ttm_tube = ttm(ttm_tube,X_tube,j);
    end
    C(:,i) = double(ttm_tube);
end
C = X_sub_mat{it_mod}*C;
end