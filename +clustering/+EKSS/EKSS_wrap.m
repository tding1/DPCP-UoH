function [err,groups] = EKSS_wrap(X,n,trueLabels,B,T)

[D, N] = size(X);
[err, groups, ~] = EKSS(X,n,D-1,trueLabels,B, N,T, 0, 1, 'kss');

end