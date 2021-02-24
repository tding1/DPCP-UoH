function [y] = c_max(U, X)
% Compute N * c_{X,max}
% X: inlier points in subspace S
% U: orthonormal basis of S
[~,y] = util.psgm_max(U' * X);
end

