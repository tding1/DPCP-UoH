function [groups, t] = main(X,K,r,affine,alpha,outlier,rho)

tic;

if (nargin < 7)
    rho = 1;
end
if (nargin < 6)
    outlier = false;
end
if (nargin < 5)
    alpha = 20;
end
if (nargin < 4)
    affine = false;
end
if (nargin < 3)
    r = 0;
end

[groups] = SSC(X,K,r,affine,alpha,outlier,rho);

t = toc;

end