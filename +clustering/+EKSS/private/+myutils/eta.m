function [y] = eta(P, X)
% Simulate eta_X
% P: orthonormal projection to S
D = size(X, 1);
k = 1e5;
y = -1;
s = normc(P * randn(D, k));
v = X*sign(X'*s);
pv = P*v;
for i = 1:k
    t = norm(pv(:,i) - s(:,i)*s(:,i)'*v(:,i));
    if t > y
        y = t;
    end
end
end

