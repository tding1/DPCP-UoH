function permLabels = bestMap(trueLabels,estLabels)
% Best permutation of estimated labels to match true labels
%   permLabels = bestMap(trueLabels,estLabels)
% Input:
%   trueLabels: vector of true labels
%   estLabels: vector of estimated labels to evaluate
% Output:
%   permLabels: estimated labels permuted to best match true labels

trueLabelVals = unique(trueLabels);
kTrue = length(trueLabelVals);
estLabelVals = unique(estLabels);
kEst = length(estLabelVals);

cost_matrix = zeros(kEst,kTrue);
for ii = 1:kEst
    inds = find(estLabels == estLabelVals(ii));
    for jj = 1:kTrue
        cost_matrix(ii,jj) = length(find(trueLabels(inds) == trueLabelVals(jj)));
    end
end

[rInd,cInd] = linear_sum_assignment(-cost_matrix);

permLabels = inf*ones(size(estLabels));
for ii = 1:length(rInd)
    permLabels(estLabels==estLabelVals(rInd(ii))) = trueLabelVals(cInd(ii));
end

outLabelVals = unique(permLabels);
if length(outLabelVals) < max(permLabels)
    lVal = 1;
    for ii = 1:length(outLabelVals)
        permLabels(permLabels==outLabelVals(ii)) = lVal;
        lVal = lVal + 1;
    end
end
