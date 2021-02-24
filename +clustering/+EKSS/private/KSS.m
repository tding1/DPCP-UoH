function [trueErr,l2Err,estLabels,state] = KSS(X,K,d,trueLabels,Uinit,maxIter)
% The K-Subspaces algorithm for subspace clustering
%   [trueErr,l2Err,estLabels] = KSS(X,K,d,trueLabels,Uinit)
% Input:
%   X: ambient dimension x number of data points
%   K: number of subspaces
%   d: dimension of subspaces (assumed all equal)
%   trueLabels: vector of true labels
%   Uinit: optional initialization of subspaces
% Output:
%   trueErr: true clustering error
%   l2Err: final cost
%   estLabels: estimated labels
%   state: struct array with elements
%       iter: number of iterations run
%       Ufull: length-K cell array of subspace bases
%       fullResids: N x K array of residuals from points to subspaces
%--------------------------------------------------------------------------
% Copyright @ John Lipor, 2018
%--------------------------------------------------------------------------

[D,N] = size(X);

if nargin < 5 || isempty(Uinit)
    % randomly initialize true subspaces
    Uinit = cell(K,1);
    fullResids = zeros(N,K);
    for kk = 1:K
        Uinit{kk} = orth(randn(D,d));
        fullResids(:,kk) = norms(X - Uinit{kk}*(Uinit{kk}'*X));
    end
    [minResids,estLabels] = min(fullResids,[],2);    % label by nearest subspace
else
    % use input initialization
    fullResids = zeros(N,K);
    for kk = 1:K
        fullResids(:,kk) = norms(X - Uinit{kk}*(Uinit{kk}'*X));
    end
    [minResids,estLabels] = min(fullResids,[],2);    % label by nearest subspace
end

if nargin < 6
    maxIter = 1000;
end

% alternate between subspace estimation and assignment
prevLabels = rand(size(estLabels));
iter = 1;
Ufull = cell(K,1);
while sum(estLabels ~= prevLabels) > 0 && iter < maxIter
    for kk = 1:K
        % first update residuals after labels obtained
        Z = X(:,estLabels==kk);  % form subspace
        if isempty(Z)
            U = orth(randn(D,d));
            Ufull{kk} = U;
        else
            [U,V] = eig(Z*Z');
            [~,vind] = sort(diag(V),'descend');
            U = U(:,vind(1:d));
            Ufull{kk} = U;
        end
        fullResids(:,kk) = norms(X - U*(U'*X));   % find distances to subspace
    end
    prevLabels = estLabels;
    [minResids,estLabels] = min(fullResids,[],2);  % label by nearest subspace
    iter = iter + 1;
end

l2Err = norm(min(minResids,[],2))/norm(X,'fro');
% compute misclassification rate if true labels given
if length(trueLabels) > 1
    trueErr = missRate(trueLabels,estLabels);
else
    trueErr = nan;
end

state = struct;
state.iter = iter;
state.Ufull = Ufull;
state.fullResids = fullResids;
