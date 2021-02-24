% -------------------------------------------------------------------------
% This function implements Algorithm 1 in [arXiv:1706.03860]

% X: The given data matrix

% mu: The parameter mu in the ADMM algorithm

% gamma: The pa rameter gamma in the algorithm

% itr: The maximum number of iterations allowed for the ADMM algorithm

% nrm: It determines if we use ell_1 norm or ell_2 norm.

% c_neighborhood: Cardinality of a neighborhood set for each data point

% error: clustering error

% labels: obtained labels

% Xlabels: clustering ground-truth
%--------------------------------------------------------------------------
% Copyright @ Mostafa Rahmani, 2017
%--------------------------------------------------------------------------

function [labels, t] = DSC_function(X,mu,gamma,itr,nrm,c_neighborhood,Xlabels) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the optimization problem and constructing the neighborhood matrix
tic;

[~,N2] = size(X) ;
T = Solver_DSC(X, X , mu , gamma , itr , nrm) ;
T = abs(T - diag(diag(T))) ;
Z=zeros(N2) ;

for k = 1:N2
    [~,b] = sort(T(:,k) , 'descend') ;
    b = b(1:c_neighborhood ) ;
    inner_product = abs(X(: , k)'*X(: , b)) ;
    Z(k, b ) = exp(-2*acos(inner_product));
end

Z = Z + Z';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (normalized) spectral clustering step

D = ( 1./sqrt(sum(Z)+eps) ); 
H = Z.*repmat(D' , 1 , N2) ; H = H.*repmat(D , N2 , 1) ;
L = speye(N2) - H; 
[~,~,V] = svd(L);
n_classes =  max(Xlabels) ;
V = V(:,N2-n_classes+1:N2);
V = normr(V);   
warning off;
maxiter = 1000;     % number iterations
replicates = 20;    % number of replications of kmeans algorithm
labels = kmeans(V,n_classes,'maxiter',maxiter,'replicates',replicates,'EmptyAction','singleton');

t = toc;

