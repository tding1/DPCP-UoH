function [t, b, f] = rsgm(Xtilde,mu_min,maxiter)
% Riemannian subgradient method.

t_start = tic;

if nargin < 2
    mu_min = 1e-10; maxiter = 200;
end

if nargin <3
    maxiter = 200;
end

mu_0 = 1e-1; alpha = 1e-3; beta = 0.5;

% initialization
[U,~,~] = svd(Xtilde, 'econ'); 
b = U(:,end);

mu = mu_0;
grad = Xtilde*sign(Xtilde'*b);
grad = grad - b*(grad'*b);
grad_norm = norm(grad)^2;
% line search
tmp = b - mu*grad;
obj_old = utils.obj_func(Xtilde, b);
while (utils.obj_func(Xtilde, tmp/norm(tmp)) > obj_old - alpha*mu*grad_norm) && mu>mu_min
    mu = mu*beta;
    tmp = b - mu*grad;
end

b = tmp / norm(tmp);

DJ = inf;
eps = 1e-3;
J_old = utils.obj_func(Xtilde, b);
i = 1;
while (abs(DJ) > eps) && i <= maxiter
    i = i+1;
    grad = Xtilde*sign(Xtilde'*b);
    grad = grad - b*(grad'*b);
   
    tmp = b - mu*grad;
    b = tmp / norm(tmp);
    
    J_new = utils.obj_func(Xtilde, b);
    DJ = (J_old - J_new) / (J_old + 10^(-6));
    J_old = J_new;
    
    mu = mu * 0.9;
end

f = utils.obj_func(Xtilde, b);

t = toc(t_start);
end

