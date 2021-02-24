function [groups, t] = main(X,K)

tic;

R = OMP_mat_func(X, 5, 1e-8);
N = size(X,2);     

R(1:N+1:end) = 0;
A = abs(R) + abs(R)';
groups = SpectralClustering(A, K, 'Eig_Solver', 'eigs');

t = toc;

end