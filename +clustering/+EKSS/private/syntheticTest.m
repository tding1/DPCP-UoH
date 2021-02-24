% Example usage of Ensemble K-Subspaces algorithm from "Subspace Clustering
% using Ensembles of K-Subspaces," by Lipor, Hong, Tan, and Balzano.
% NOTE: To run, download the SSC author code and place SpectralClustering.m
% in this folder. To use CoP-KSS as a base clustering algorithm, download
% the code at the website of John Lipor or Laura Balzano.

clear

D = 4;
d = D-1;
K = 4;
Nk = 200;
varn = 0.01;

% generate random data
[X,trueLabels,Utrue] = genSubspaceData(D,d,K,Nk,varn);
[D,N] = size(X);

% algorithm parameters
B = 1000;       % number of base clusterings
q = N;          % thresholding parameter
T = 3;          % number of KSS iterations to run
heur = 0;       % weight clustering based on l2-error of KSS
sc = 1;         % run spectral clustering
base = 'kss';   % base clustering algorithm

% run B individual KSS clusterings and take lowest error
% kssErr = inf;
% for bb = 1:B
%     err = KSS(X,K,d,trueLabels);
%     if err < kssErr
%         kssErr = err;
%     end
% end
% kssErr

% compare to EKSS
[ekssErr, gg] = EKSS(X,K,d,trueLabels,B,q,T,heur,sc,base);
ekssErr






