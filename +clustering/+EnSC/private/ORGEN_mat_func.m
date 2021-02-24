function C = ORGEN_mat_func(X, EN_solver, varargin)
%ORGEN_MAT_FUNC Perform elstic net by ORGEN for self representation.
%   This code implements the subspace clustering algorithm described in
% 
%   Chong You, Chun-guang Li, Daniel Robinson, Rene Vidal,
%   "Oracle Based Active Set Algorithm for Scalable Elastic Net Subspace
%   Clustering", CVPR 2016.
% 
% 	It perform elstic net for each column of data X = [x_1, \dots, x_N]
%   using all other columns as a dicitonary. 
%   i.e., for each j = 1, \dots, N, solves the following:
%   min_{c_j} lambda ||c_j||_1 + (1-lambda)/2 ||c_j||_2^2 + nu/2 ||x_j - X
%   c_j||_2^2, c_jj = 0. (*)
%   The output C is given by [c_1, \dots, c_N]. The algorithm for the above
%   optimization problem is ORGEN. See also ORGEN.m

% Input Arguments
% X                 -- data matrix D by N where each column is a data point.
% EN_solver, Nsample, lambda: same as ORGEN.m
% maxiter           -- max number of iterations in ORGEN.
% nu0, nu_method    -- if nu_method == 'fixed', then nu = nu0 for all
%                      columns. If nu_method == 'nonzero' then 
%                      nu(j) = nu0 * k(j) where k(j) is the largest number
%                      such that the solution to (*) above is zero if nu =
%                      k(j).

% Copyright Chong You @ Johns Hopkins University, 2016
% chong.you1987@gmail.com

MEMORY_TOTAL = 0.2 * 10^9; % memory available for double precision.
% Check input parameters
if mod(length(varargin), 2) ~= 0
   error(['' mfilename '' 'needs propertyName/propertyValue pairs']);
end
% Set default 
vararg = {'Nsample', 100, ...
          'maxiter', 5, ...
          'nu0', 20, ...
          'nu_method', 'fixed', ...
          'lambda', 0.9, ...
          'outflag', true};
% Overwrite by input
vararg = vararginParser(vararg, varargin);
% Generate variables
for pair = reshape(vararg, 2, []) % pair is {propName;propValue}
   eval([pair{1} '= pair{2};']);
end

% 
N = size(X, 2);
D = size(X, 1);
blockSize = round(MEMORY_TOTAL / N);
if Nsample > N
    error('Nsample is bigger than N\n')
end

% Preallocation
rows = cell(1, N);
vals = cell(1, N);

% Compute Dic; 
% Compute nu;
if strcmpi(nu_method, 'fixed')
    Dic = X' / (eye(D)/nu0+X*X');
    nu = nu0 * ones(1, N);
else
    Dic = X' / (eye(D)/ (nu0*lambda*100) + X*X');
    nu = zeros(1, N);
%     % The following while loop does the following:
%     I = X' * X;
%     diag(I) = 0;
%     nu = nu0 * lambda / max(abs(I), [], 1);
%     % essentially, X' * X is too big if N is large, so the problem is
%     % decomposed into blocks. Comparing to decomposing into single
%     % columns-wise processing, the computation can be accelarated by
%     % using blocks.
    counter = 0;
    while(1)
        mask = [counter+1 : min(counter + blockSize, N)];
        I = abs(X' * X(:, mask));
        I(counter+1:N+1:end) = 0; % set diagonal = 0
        nu(mask) = nu0 * lambda ./ max(I, [], 1);
        counter = counter + blockSize;
        if counter >= N
            break;
        end
    end
end

%
V = zeros(D, N); % oracle
counter = 0;
while(1)
    mask = [counter+1 : min(counter + blockSize, N)];
    % initialize active set
    I = abs(Dic * X(:, mask));
    I(counter+1:N+1:end) = 0; % set diagonal = 0
    [~, Supp] = sort(I, 1, 'descend');
    Supp = Supp(1:Nsample, :);
    % compute rows, vals, V
    for ii = 1:length(mask)
        index = mask(ii); % index in 1:N
        % S1: Xs -> oracle v 
        x = X(:, index);
        Xs = X(:, Supp(:, ii));
        cs = EN_solver(Xs, x, lambda, nu(index));
        
        V(:, index) = nu(index) * (x - Xs*cs); % oracle 
        supp_mask = (abs(cs) > eps);
        rows{index} = Supp(supp_mask, ii);
        vals{index} = cs(supp_mask);
    end
    
    counter = counter + blockSize;
    if counter >= N
        break;
    end
end

indexing = 1:N;
for iter = 2:maxiter
    counter = 0;
    while(1)
        mask = indexing(counter+1 : min(counter + blockSize, length(indexing)));
    
        % S2: guided neighborhood.
        coherence = abs(X' * V(:, mask));
        coherence(sub2ind(size(coherence), mask, 1:length(mask))) = 0;
    
        for ii = 1:length(mask)
            index = mask(ii);
            Supp = find( coherence(:, ii) > lambda); % a mask for C that out of mask is zero.
            if length(Supp) == length(rows{index})
                indexing(indexing==index) = 0;
                continue;
            end

            if length(Supp) > Nsample % control the size of the problem
                coherence(rows{index}, ii) = 0;
                addedsupp = find(coherence(:, ii) > lambda); %
                [~, ind] = sort( coherence(addedsupp, ii), 'ascend' ); % 
                if length(Supp) - Nsample == length(ind) % all addedsupp should be removed from newsupp
                    fprintf('Warning: not converged termination\n');
                    indexing(indexing==index) = 0;
                    continue;
                end
                ind = ind( 1: length(Supp) - Nsample ); % 
                Supp = setdiff(Supp, addedsupp(ind)); %
            end

            % S1: Xs -> oracle v 
            x = X(:, index);
            Xs = X(:, Supp);
            cs = EN_solver(Xs, x, lambda, nu(index));
            V(:, index) = nu(index) * (x - Xs*cs); % oracle 
            
            supp_mask = (abs(cs) > eps);
            rows{index} = Supp(supp_mask);
            vals{index} = cs(supp_mask);
        end
        
        counter = counter + blockSize;
        if counter >= length(indexing)
            break;
        end
    end
   
    indexing(indexing == 0) = [];
    if outflag
        fprintf('Iteration: %d, No. remaining: %d in %d\n', iter-1, length(indexing), N)
    end
    if isempty(indexing) % all have converged
        break;
    end
end

counter = 0;
for ii = 1:N
    counter = counter + length(rows{ii});
end
C = sparse([],[],[], N, N, counter);
for ii = 1:N
    C(rows{ii}, ii) = vals{ii};
end

% if outflag
%     Obj_val = lambda * sum(sum(abs(C))) + (1 - lambda)/2 * sum(sum(C .^2)) + nu/2 * sum(sum((X - X*C) .^2));
%     fprintf('obj: %2.4f, sparsity: %f\n', Obj_val, counter/N^2);
% end


