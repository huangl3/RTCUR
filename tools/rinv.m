% best rank r approximation of inv(X)
function rinvX = rinv(X, r)
[U,S,V] = svd(X,'econ');
rinvX = V(:,1:r)*pinv(S(1:r,1:r))*U(:,1:r)';
end