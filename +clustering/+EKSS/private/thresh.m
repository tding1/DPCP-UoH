function Aq = thresh(A,q)
% Threshold affinity matrix to keep top q entries of each row/column
%   Aq = thresh(A,q)
% Input:
%   A: affinity matrix
%   q: number of entries to keep in each row/column
% Output:
%   Aq: thresholded affinity matrix
%--------------------------------------------------------------------------
% Copyright @ John Lipor, 2018
%--------------------------------------------------------------------------

N = size(A,1);
A = A - diag(diag(A));

Aq = zeros(N);
for ii = 1:N
    ai = A(:,ii);
    ai(ii) = 0;
    [val,ind] = sort(ai,'descend');
    Aq(ind(1:q),ii) = val(1:q);
    Aq(ii,ind(1:q)) = val(1:q);
end
Aq = Aq + Aq';
