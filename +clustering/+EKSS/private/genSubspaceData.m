function [X,trueLabels,U] = genSubspaceData(D,d,K,Nk,varn,theta)
% Generate data from a union of subspaces.
%   [X,trueLabels] = genSubspaceData(D,d,K,Nk,varn,theta)
% Input:
%   D: ambient dimension
%   d: dimension of subspaces (assumed all equal)
%   K: number of subspaces
%   Nk: points per subspace
%   varn: variance of additive Gaussian noise (default 0)
%   theta: principle angles between all subspaces if K = 3 or pairs of subspaces if K > 2 (default 0, which results in random subspace generation)
% Output:
%   X: data of shape D x N
%   trueLabels: labels of data
%   U: cell array of size K containing true subspace bases
%--------------------------------------------------------------------------
% Copyright @ John Lipor, 2018
%--------------------------------------------------------------------------

if nargin < 5
    varn = 0;
end

if nargin < 6
    theta = 0;
end

if theta > 0 && mod(K,2) > 0 && K ~= 3
    disp('Error: To set principal angles, must have either K = 3 or K even.')
    return
end

if theta > 0
    if K == 3   % CURRENTLY BROKEN
        % set principal angles between all subspaces
        fullU = orth(randn(D));
        for dd = 1:d
            u = fullU(:,(dd-1)*3+1:dd*3);
            U1(:,dd) = u(:,1);
            U2(:,dd) = cos(theta)*u(:,1) + sin(theta)*u(:,2);
            U3(:,dd) = cos(theta)*u(:,1) + sin(theta)*u(:,3);
        end
        U{1} = U1;
        U{2} = U2;
        U{3} = U3;
    else
        % set principal angles to theta for each pair of subspaces
        for kk = 1:K/2
            fullU = orth(randn(D));
            for dd = 1:d
                u = fullU(:,(dd-1)*2+1:dd*2);
                U1(:,dd) = u(:,1);
                U2(:,dd) = cos(theta)*u(:,1) + sin(theta)*u(:,2);
            end
            U{(kk-1)*2+1} = U1;
            U{(kk-1)*2+2} = U2;
        end
    end
else
    for kk = 1:K
        U{kk} = orth(randn(D,d));
    end
end
% generate data
X = [];
trueLabels = [];
for kk = 1:K
    x = U{kk}*randn(d,Nk);
    X = [X x];
    trueLabels = [trueLabels; kk*ones(Nk,1)];
end
[D,N] = size(X);
X = X + sqrt(varn)*randn(size(X));
X = X./repmat(sqrt(sum(X.^2,1)),D,1);
