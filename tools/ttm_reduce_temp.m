function C = ttm_reduce_temp(core,X1,X2,X3,J,dim)
if dim == 1
    d = size(X2,1);
else
    d = size(X1,1);
end
s = length(J);
p = index2D(J,d);
C = zeros(size(core,dim),s);
if dim == 1
    for i = 1:s
    C(:,i) = ttm(core,{X2(p(i,1),:),X3(p(i,2),:)},[2,3]);
    end
    C = X1*C;
end
if dim == 2
    for i = 1:s
    C(:,i) = ttm(core,{X1(p(i,1),:),X3(p(i,2),:)},[1,3]);
    end
    C = X2*C;
end
if dim == 3
    for i = 1:s
    C(:,i) = ttm(core,{X1(p(i,1),:),X2(p(i,2),:)},[1,2]);
    end
    C = X3*C;
end
end

function p = index2D(J, d1)
sz = size(J,2); p = zeros(sz,2);
p(:,1) = mod(J-1,d1)+1;
p(:,2) = floor((J-1)./d1)+1;
end