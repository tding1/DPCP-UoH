function [normals, groups, elapsed_time] = KSS(method, Xtilde, init_normals, n, t_budget, threshold, maxiter)
% KSS framework

t_start = tic;

normals = init_normals;

% initial clustering
[~,groups] = min(abs(normals'*Xtilde));

% initial cost
J_old = 0;
for i = 1 : n
    mask = groups==i;
    J_old = J_old + norm(Xtilde(:,mask)'*normals(:,i),1);
end

% loop over KH iterations
iter = 0;
DJ = inf;
eps = 1e-10;
while (abs(DJ) > eps) && (toc(t_start) < t_budget) && iter < maxiter
    iter = iter + 1;
    
    % update groups and normals
    [normals, groups] = LOCAL_NORMALS(method, Xtilde,groups,normals, threshold, t_budget/n);
    
    % new cost
    J_new = 0;
    for i = 1 : n
        mask = groups==i;
        J_new = J_new + norm(Xtilde(:,mask)'*normals(:,i),1);
    end
    DJ = (J_old - J_new) / (J_old + 10^(-6));
    J_old = J_new;
end

elapsed_time = toc(t_start);

end


function [NORMALS, GROUPS] = LOCAL_NORMALS(method, X,GROUPS,NORMALS,threshold, t_budget)

% DESCRIPTION: This function takes as input a dataset X together with
% a clustering (GROUPS) and a normal vector for each cluster
% (NORMALS). The function fits afresh a normal vector to each cluster
% and recomputes the clustering based on these normals. This function is
% used internally inside the function Iterative_DPCP_HC.

[D, n] = size(NORMALS);

if strcmp(method, 'rsgm')
    for i = 1 : n
        mask = GROUPS == i;
        if ~any(mask)
            NORMALS(:,i) = normc(randn(D, 1));
        else
            [~, normal] = solvers.DPCP.rsgm(X(:,mask));        
            NORMALS(:,i) = normal;
        end
    end
elseif strcmp(method, 'pca')
    for i = 1 : n
        mask = GROUPS == i;
        if ~any(mask)
            NORMALS(:,i) = normc(randn(D, 1));
        else
            [U,~] = svd(X(:,mask), 'econ');
            NORMALS(:,i) = U(:, end);
        end
    end
elseif strcmp(method, 'cop')
    for i = 1 : n
        mask = GROUPS == i;
        if ~any(mask)
            NORMALS(:,i) = normc(randn(D, 1));
        else
            nSamp = max(D-1,size(X(:,mask),2)-2);
            U = solvers.CoP.CoP(X(:,mask),D-1,nSamp,2);
            NORMALS(:,i) = null(U');
        end
    end
end

% update the clustering
[~,GROUPS] = min(abs(NORMALS'*X));

end