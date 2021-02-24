function U = CoP(X,r,nSamp,p,G)
% Coherence Pursuit algorithm for Robust Subspace Recovery
%   U = CoP(X,r,nSamp,p,G)
% Input:
%   X: ambient dimension x number of data points
%   r: rank of subspace to estimate
%   nSamp: number of inliers used to estimate subspace
%   p: optional choice of 1-norm or 2-norm outlier rejection
%   G: optional Gram matrix
% Output:
%   U: estimated orthonormal basis for r-dimensional subspace
%--------------------------------------------------------------------------
% Copyright @ John Lipor, 2018
%--------------------------------------------------------------------------

[D,N] = size(X);

% normalize data
X = X./repmat(sqrt(sum(X.^2,1)),D,1);
X(isnan(X)) = 0;

if nargin < 4
    p = 2;
end

% compute Gram matrix
if nargin < 5
    G = abs(X'*X);
end
G = G - diag(diag(G));

if p == 1
    g = sum(abs(G));
elseif p == 2
    g = sum(G.^2);
else
    G = 1./(1 + exp(-10*(G-0.4)));
    g = sum(abs(G));
    g = sum(abs(G));
end
g = g/max(g); 

% take top nSamp inds
[~,inds] = sort(g,'descend');
Y = X(:,inds(1:nSamp)); 

% find subspace basis
[U,~,~] = svd(Y,'econ'); 
U = U(:,1:r); 
