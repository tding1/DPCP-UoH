function [y] = obj_func(X,b)
% Compute the objective function of DPCP problem ||X^T b||_1
y = norm(X'*b, 1);
end

