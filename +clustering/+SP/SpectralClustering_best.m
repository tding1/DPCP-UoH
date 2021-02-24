%--------------------------------------------------------------------------
% This function takes an adjacency matrix of a graph and computes the 
% clustering of the nodes using the spectral clustering algorithm of 
% Ng, Jordan and Weiss.
% CMat: NxN adjacency matrix
% n: number of groups for clustering
% groups: N-dimensional vector containing the memberships of the N points 
% to the n groups obtained by spectral clustering
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------

function [BestGROUPS, GROUPS] = SpectralClustering_best(CKSym,n,X)
normType = 1;GROUPS = [];
warning off;
N = size(CKSym,1);
MAXiter = 100; % Maximum number of iterations for KMeans 
REPlic = 10; % Number of replications for KMeans

% Normalized spectral clustering according to Ng & Jordan & Weiss
% using Normalized Symmetric Laplacian L = I - D^{-1/2} W D^{-1/2}

DN = diag( 1./sqrt(sum(CKSym)+eps) );
LapN = speye(N) - DN * CKSym * DN;
[uN,sN,vN] = svd(LapN);
kerN = vN(:,N-n+1:N);
for i = 1:N
    kerNS(i,:) = kerN(i,:) ./ norm(kerN(i,:)+eps);
end
J_best = inf;
for i = 1:REPlic
    [groups C] = kmeans(kerNS,n,'maxiter',MAXiter,'replicates',1,'EmptyAction','singleton');
    J_new = 0;normals = zeros(size(X,1), n);
   for j = 1 : n
       groupi = find(groups==j);
       if normType == 1
            [~,normals(:,j),~] = solvers.DPCP.rsgm(X(:,groupi));
           J_new = J_new + norm(abs(X(:,groupi)'*normals(:,j)),1);
       elseif normType == 2
           [U_temp,~,~] = svd(X(:,groupi));
           normals(:,j) = U_temp(:, end);
           J_new = J_new + norm(abs(X(:,groupi)'*normals(:,j)))^2;
       end
   end
   GROUPS = [GROUPS,groups];
   if J_new < J_best
        J_best = J_new;
        BestGROUPS = groups;
        NORMALS = normals;
    end
end


% for i = 1 : 10
%     [groups C] = kmeans(kerNS,[],'maxiter',MAXiter,'EmptyAction','singleton','start',C);
% end