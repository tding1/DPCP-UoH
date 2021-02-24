function [err,estLabels,A] = EKSS(X,K,d,trueLabels,B,q,T,heur,sc,base)
% The Ensemble K-Subspaces algorithm for subspace clustering
%   [trueErr,l2Err,estLabels] = EKSS(X,K,d,trueLabels,B,q,heur,sc,base)
% Input:
%   X: ambient dimension by number of data points
%   K: number of subspaces
%   d: dimension of subspaces (assumed all equal)
%   trueLabels: vector of true labels
%   B: number of base clusterings
%   q: threshold parameter (default is no thresholding)
%   T: number of KSS iterations (default 30, set to 0 for EKSS-0)
%   heur: set to 1 to add a value besides 1 for co-clustered points (default 0)
%   sc: select whether to run spectral clustering or not (default 1)
%   base: base clustering algorithm, choose from 'kss' and 'copKSS' (default 'kss')
% Output:
%   err: clustering error
%   estLabels: estimated labels
%   A: affinity matrix
%--------------------------------------------------------------------------
% Copyright @ John Lipor, 2018
%--------------------------------------------------------------------------

[D,N] = size(X);

if nargin < 6 || length(q) == 0
    q = N;
end

if nargin < 7 || length(T) == 0
    T = 30;
end

if nargin < 8 || length(heur) == 0
    heur = 0;
end

if nargin < 9 || length(sc) == 0
    sc = 1;
end

if nargin < 10 || length(base) == 0
    base = 'kss';
end

% form affinity matrix
A = zeros(N);
for bb = 1:B
    if strcmp(base,'kss')
        normals = normc(rand(D, K));
        for kk = 1:K
            Uinit{kk} = null(normals(:, kk)');
        end
        [trueErr,l2Err,estLabels,state] = KSS(X,K,d,trueLabels,Uinit,T);
    else
        [trueErr,l2Err,estLabels,state] = copKSS(X,K,d,trueLabels,[],T,2);
    end
    for kk = 1:K
        kInds = find(estLabels==kk);
        if heur == 0
            A(kInds,kInds) = A(kInds,kInds) + 1;
        else
            A(kInds,kInds) = A(kInds,kInds) + 1 - l2Err;
        end
    end
end
A = A/B;

% threshold
A = A - diag(diag(A));
if q < N
    Aq = thresh(A,q);
else
    Aq = A;
end

% estimate labels using spectral clustering
if sc == 1
    estLabels = SpectralClustering(Aq,K);
else 
    estLabels = nan;
end

% calculate error
if length(trueLabels) > 1 && sc == 1
    [err,estLabels] = missRate(trueLabels,estLabels);
else
    err = nan;
end
