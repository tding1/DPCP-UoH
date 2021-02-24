function [groups, t] = main(X, lambda, nu0, Nsample, K)

tic;

EN_solver =  @(X, y, lambda, nu) rfss( X, y, lambda / nu, (1-lambda) / nu );

R = ORGEN_mat_func(X, EN_solver, 'nu0', nu0, 'nu_method', 'nonzero', 'lambda', lambda, ...
                                                          'Nsample', Nsample, 'maxiter', 5, 'outflag', true); 
N = size(X,2);                                                
R(1:N+1:end) = 0;
A = abs(R) + abs(R)';
groups = SpectralClustering(A, K, 'Eig_Solver', 'eigs');       

t = toc;

end