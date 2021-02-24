function [err,permLabels] = missRate(trueLabels,estLabels)
% Clustering error/misclassification rate
%   [err,estLabels] = missRate(trueLabels,estLabels)
% Input:
%   trueLabels: vector of true labels
%   estLabels: vector of estimated labels to evaluate
% Output:
%   err: true clustering error
%   permLabels: estimated labels permuted to best match true labels
%--------------------------------------------------------------------------
% Copyright @ John Lipor, 2018
%--------------------------------------------------------------------------

permLabels = bestMap(trueLabels,estLabels);
err  = sum(trueLabels ~= permLabels) / length(trueLabels);
