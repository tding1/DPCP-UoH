clc
clear

randn('seed',2021);rand('seed',2021)

D = 4;
K = 2;

acc_scc = 0;                 t_scc = 0;
acc_mkf = 0;                 t_mkf = 0;

acc_dpcpkss_10 = 0;          t_dpcpkss_10 = 0;
acc_pcakss_10 = 0;           t_pcakss_10 = 0;
acc_copkss_10 = 0;           t_copkss_10 = 0;

acc_pcakss_10_core = 0;      t_pcakss_10_core = 0;
acc_copkss_10_core = 0;      t_copkss_10_core = 0;
acc_dpcpkss_10_core = 0;     t_dpcpkss_10_core = 0;

acc_epcakss = 0;             t_epcakss = 0;
acc_ecopkss = 0;             t_ecopkss = 0;
acc_edpcpkss = 0;            t_edpcpkss = 0;

acc_ensc = 0;           t_ensc = 0;
acc_ssc_admm = 0;       t_ssc_admm = 0;
acc_ssc_omp = 0;        t_ssc_omp = 0;

parms = rand_parms_spec();
parms.D = D;
parms.K = K;
parms.alpha = 1;
parms.r = 0.3;

info = info_compute(parms);
N = size(info.Xtilde,2);

%% SCC
[C, ~, t] = clustering.SCC.lscc(info.Xtilde', parms.D-1, parms.K);
t_scc = t_scc + t;
acc_scc = acc_scc + utils.compute_accuracy(info.C, C(1:length(info.C)));

%% MKF
[~,C,~, ~, t] = clustering.MKF.kmcl1sd(info.Xtilde', (parms.D-1)*ones(1, parms.K), 0.001, 10000, 2, 0);
t_mkf = t_mkf + t;
acc_mkf = acc_mkf + utils.compute_accuracy(info.C, C(1:length(info.C)));

%% EnSC
lambda_ensc = 0.95;
alpha_ensc = 3;
Nsample_ensc = 200;
[groups,t] = clustering.EnSC.main(info.Xtilde, lambda_ensc, alpha_ensc, Nsample_ensc, parms.K);
acc_ensc = acc_ensc + utils.compute_accuracy(info.C, groups(1:length(info.C)));
t_ensc = t_ensc + t;

%% SSC-ADMM
[groups,t] = clustering.SSC_ADMM.main(info.Xtilde, parms.K);
acc_ssc_admm = acc_ssc_admm + utils.compute_accuracy(info.C, groups(1:length(info.C)));
t_ssc_admm = t_ssc_admm + t;

%% SSC-OMP
[groups,t] = clustering.SSC_OMP.main(info.Xtilde, parms.K);
acc_ssc_omp = acc_ssc_omp + utils.compute_accuracy(info.C, groups(1:length(info.C)));
t_ssc_omp = t_ssc_omp + t;

%% KSS
kk = 10;
for i = 1:kk
    normals = normc(randn(D, K));

    [normals_pcakss, C0, t_bud] = clustering.KSS.KSS('pca', info.Xtilde, normals, parms.K, inf, -1, 100);
    t_pcakss_10 = t_pcakss_10 + t_bud;
    meta_normals_pcakss{i} = normals_pcakss;
    meta_groups_pcakss{i} = C0;

    [normals_dpcpkss, C0, t] = clustering.KSS.KSS('rsgm',info.Xtilde, normals, parms.K, inf, -1,100);
    t_dpcpkss_10 = t_dpcpkss_10 + t;
    meta_normals_dpcpkss{i} = normals_dpcpkss;
    meta_groups_dpcpkss{i} = C0;

    [normals_copkss, C0, t] = clustering.KSS.KSS('cop', info.Xtilde, normals, parms.K, inf, -1, 100);
    t_copkss_10 = t_copkss_10 + t;
    meta_normals_copkss{i} = normals_copkss;
    meta_groups_copkss{i} = C0;

end
acc_dpcpkss_10 = acc_dpcpkss_10 + compute_acc(meta_groups_dpcpkss, meta_normals_dpcpkss, info);
acc_pcakss_10 = acc_pcakss_10 + compute_acc(meta_groups_pcakss, meta_normals_pcakss, info);
acc_copkss_10 = acc_copkss_10 + compute_acc(meta_groups_copkss, meta_normals_copkss, info);

%% CoRe-KSS
t = tic;
[~, best_loss] = compute_acc(meta_groups_pcakss, meta_normals_pcakss, info);
best_meta_groups = meta_groups_pcakss;
best_meta_normals = meta_normals_pcakss;
for s = 1:10
    inner_best_loss = Inf;
    for i = 1:kk
        [groups, normals] = CoRe(info.Xtilde, i, meta_normals_pcakss, meta_groups_pcakss, ceil(K/2), 1e-3);
        [normals, groups, ~] = clustering.KSS.KSS('pca',info.Xtilde, normals, K, inf, -1, 100);
        l = loss(info.Xtilde, groups, normals);
        meta_groups_pcakss{i} = groups;
        meta_normals_pcakss{i} = normals;
        if l < inner_best_loss
            inner_best_loss = l;
        end
    end
    if inner_best_loss < 0.999 * best_loss
        best_loss = inner_best_loss;
        best_meta_groups = meta_groups_pcakss;
        best_meta_normals = meta_normals_pcakss;
    else
        meta_groups_pcakss = best_meta_groups;
        meta_normals_pcakss = best_meta_normals;
        break;
    end
end
t_pcakss_10_core = t_pcakss_10_core + toc(t);
acc_pcakss_10_core = acc_pcakss_10_core + compute_acc(meta_groups_pcakss, meta_normals_pcakss, info);


t = tic;
[~, best_loss] = compute_acc(meta_groups_copkss, meta_normals_copkss, info);
best_meta_groups = meta_groups_copkss;
best_meta_normals = meta_normals_copkss;
for s = 1:10
    inner_best_loss = Inf;
    for i = 1:kk
        [groups, normals] = CoRe(info.Xtilde, i, meta_normals_copkss, meta_groups_copkss, ceil(K/2), 1e-3);
        [normals, groups, ~] = clustering.KSS.KSS('copkss',info.Xtilde, normals, K, inf, -1, 100);
        l = loss(info.Xtilde, groups, normals);
        meta_groups_copkss{i} = groups;
        meta_normals_copkss{i} = normals;
        if l < inner_best_loss
            inner_best_loss = l;
        end
    end
    if inner_best_loss < 0.999 * best_loss
        best_loss = inner_best_loss;
        best_meta_groups = meta_groups_copkss;
        best_meta_normals = meta_normals_copkss;
    else
        meta_groups_copkss = best_meta_groups;
        meta_normals_copkss = best_meta_normals;
        break;
    end
end
t_copkss_10_core = t_copkss_10_core + toc(t);
acc_copkss_10_core = acc_copkss_10_core + compute_acc(meta_groups_copkss, meta_normals_copkss, info);

t = tic;
[~, best_loss] = compute_acc(meta_groups_dpcpkss, meta_normals_dpcpkss, info);
best_meta_groups = meta_groups_dpcpkss;
best_meta_normals = meta_normals_dpcpkss;
for s = 1:10
    inner_best_loss = Inf;
    for i = 1:kk
        [groups, normals] = CoRe(info.Xtilde, i, meta_normals_dpcpkss, meta_groups_dpcpkss, ceil(parms.K/2), 1e-3);
        [normals, groups, ~] = clustering.KSS.KSS('rsgm',info.Xtilde, normals, parms.K, inf, -1, 100);
        l = loss(info.Xtilde, groups, normals);
        meta_groups_dpcpkss{i} = groups;
        meta_normals_dpcpkss{i} = normals;
        if l < inner_best_loss
            inner_best_loss = l;
        end
    end
    if inner_best_loss < 0.999 * best_loss
        best_loss = inner_best_loss;
        best_meta_groups = meta_groups_dpcpkss;
        best_meta_normals = meta_normals_dpcpkss;
    else
        meta_groups_dpcpkss = best_meta_groups;
        meta_normals_dpcpkss = best_meta_normals;
        break;
    end
end
t_dpcpkss_10_core = t_dpcpkss_10_core + toc(t);
acc_dpcpkss_10_core = acc_dpcpkss_10_core + compute_acc(meta_groups_dpcpkss, meta_normals_dpcpkss, info);

t_pcakss_10_core = t_pcakss_10_core + t_pcakss_10;
t_copkss_10_core = t_copkss_10_core + t_copkss_10;
t_dpcpkss_10_core = t_dpcpkss_10_core + t_dpcpkss_10;

%% EKSS
A_pcakss = zeros(size(info.Xtilde,2));
A_copkss = zeros(size(info.Xtilde,2));
A_dpcpkss = zeros(size(info.Xtilde,2));
for i = 1:1000
    if rem(i, 100) == 0
        fprintf('EKSS iter: %d\n',i)
    end

    normals = normc(randn(D, K));

    t = tic;
    [~, C_kss, ~] = clustering.KSS.KSS('pca',info.Xtilde, normals, K, inf, -1, 10);
    for ii = 1:K
        kInds = find(C_kss==ii);
        A_pcakss(kInds,kInds) = A_pcakss(kInds,kInds) + 1;
    end
    t_epcakss = t_epcakss + toc(t);

    t = tic;
    [~, C_copkss, ~] = clustering.KSS.KSS('cop',info.Xtilde, normals, K, inf, -1, 10);
    for ii = 1:K
        kInds = find(C_copkss==ii);
        A_copkss(kInds,kInds) = A_copkss(kInds,kInds) + 1;
    end
    t_ecopkss = t_ecopkss + toc(t);

    t = tic;
    [~, C_kss, ~] = clustering.KSS.KSS('rsgm',info.Xtilde, normals, K, inf, -1, 10);
    for ii = 1:K
        kInds = find(C_kss==ii);
        A_dpcpkss(kInds,kInds) = A_dpcpkss(kInds,kInds) + 1;
    end
    t_edpcpkss = t_edpcpkss + toc(t);
end

t = tic;
A_pcakss = A_pcakss/1000;
A_pcakss = A_pcakss - diag(diag(A_pcakss));
estLabels = clustering.SP.SpectralClustering(A_pcakss,K);
t_epcakss = t_epcakss + toc(t);
acc_epcakss = acc_epcakss + utils.compute_accuracy(info.C, estLabels(1:length(info.C)));     

t = tic;
A_copkss = A_copkss/1000;
A_copkss = A_copkss - diag(diag(A_copkss));
estLabels = clustering.SP.SpectralClustering(A_copkss,K);
t_ecopkss = t_ecopkss + toc(t);
acc_ecopkss = acc_ecopkss + utils.compute_accuracy(info.C, estLabels(1:length(info.C)));     

t = tic;
A_dpcpkss = A_dpcpkss/1000;
A_dpcpkss = A_dpcpkss - diag(diag(A_dpcpkss));
estLabels = clustering.SP.SpectralClustering(A_dpcpkss,K);
t_edpcpkss = t_edpcpkss + toc(t);
acc_edpcpkss = acc_edpcpkss + utils.compute_accuracy(info.C, estLabels(1:length(info.C)));  

%% Output
fprintf('\n')
fprintf('=======================================================\n')
fprintf('              Summary (D = %d, K = %d)                   \n', D, K)
fprintf('=======================================================\n')
fprintf('Accuracy for MKF:           %.4f, elapsed time: %.2fs\n', acc_mkf, t_mkf)
fprintf('Accuracy for SCC:           %.4f, elapsed time: %.2fs\n', acc_scc, t_scc)
fprintf('Accuracy for EnSC:          %.4f, elapsed time: %.2fs\n', acc_ensc, t_ensc)
fprintf('Accuracy for SSC-ADMM:      %.4f, elapsed time: %.2fs\n', acc_ssc_admm, t_ssc_admm)
fprintf('Accuracy for SSC-OMP:       %.4f, elapsed time: %.2fs\n', acc_ssc_omp, t_ssc_omp)
fprintf('Accuracy for DPCP-KSS:      %.4f, elapsed time: %.2fs\n', acc_dpcpkss_10, t_dpcpkss_10)
fprintf('Accuracy for CoP-KSS:       %.4f, elapsed time: %.2fs\n', acc_copkss_10, t_copkss_10)
fprintf('Accuracy for PCA-KSS:       %.4f, elapsed time: %.2fs\n', acc_pcakss_10, t_pcakss_10)
fprintf('Accuracy for DPCP-EKSS:     %.4f, elapsed time: %.2fs\n', acc_edpcpkss, t_edpcpkss)
fprintf('Accuracy for CoP-EKSS:      %.4f, elapsed time: %.2fs\n', acc_ecopkss, t_ecopkss)
fprintf('Accuracy for PCA-EKSS:      %.4f, elapsed time: %.2fs\n', acc_epcakss, t_epcakss)
fprintf('Accuracy for DPCP-CoRe-KSS: %.4f, elapsed time: %.2fs\n', acc_dpcpkss_10_core, t_dpcpkss_10_core)
fprintf('Accuracy for CoP-CoRe-KSS:  %.4f, elapsed time: %.2fs\n', acc_copkss_10_core, t_copkss_10_core)
fprintf('Accuracy for PCA-CoRe-KSS:  %.4f, elapsed time: %.2fs\n', acc_pcakss_10_core, t_pcakss_10_core)

