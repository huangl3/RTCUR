function P1 = vecJ2ten_ind(V,n_dims,mod)
%%%% This function is used to transform the scalar indices to tensor indices
% similar to vec2ten_ind, this function convert J(indices between 
% 1 and d1*...*d_{i-1}*d_{i+1}*...d_n)
% into N-1 dimension vectors where each dimension in [1, d_k]
num_mod = length(n_dims);
num_indices = length(V);
if nargin<3
    mod = num_mod+1;
    n_dims = [n_dims 1];
end
prod_all_mod = prod(n_dims)/n_dims(mod);
P1 = zeros(num_indices,num_mod);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = num_mod:-1:mod+1
    prod_all_mod = prod_all_mod/n_dims(i);
    P1(:,i) = ceil(V/prod_all_mod);
    V = rem(V-1,prod_all_mod)+1;
end
for i = mod-1:-1:1
    prod_all_mod = prod_all_mod/n_dims(i);
    P1(:,i) = ceil(V/prod_all_mod);
    V = rem(V-1,prod_all_mod)+1;
end






