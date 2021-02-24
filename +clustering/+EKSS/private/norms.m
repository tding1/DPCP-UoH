function nx = norms(X,p)
% Return p-norm of each column of the matrix X
%   nx = norms(X,p)
% Input:
%   X: (D x N) matrix 
%   p: optional norm to compute, deafult is 2
% Output:
%   nx: vector of size N whose elements are the p-norms of the columns of X
%--------------------------------------------------------------------------
% Copyright @ John Lipor, 2018
%--------------------------------------------------------------------------

if nargin < 2
    p = 2;
end

nx = (sum(abs(X).^p,1)).^(1/p);
