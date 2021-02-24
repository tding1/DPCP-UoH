function [y] = c_min(U, X)
% Compute N * c_{X,min}
% X: inlier points in subspace S
% U: orthonormal basis of S
[~,y] = util.psgm(U' * X);
end