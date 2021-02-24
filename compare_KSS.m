clc
clear
close all
randn('seed',2020);rand('seed',2020)

D = 9;
K = 3;
parms = rand_parms_spec();
parms.D = D;
parms.K = K;
parms.alpha = 1;
parms.r = 0.3;

t_dpcp = 0;
t_pca = 0;
t_cop = 0;
T = 100;
maxiter = 100;
for num_datasets = 1:T
    num_datasets
    
    info = info_compute(parms);

    init = normc(randn(D, K));
    normals_dpcp = init;
    normals_pca = init;
    normals_cop = init;
    
    [~,init_groups] = min(abs(init'*info.Xtilde));
    groups_dpcp = init_groups;
    groups_pca = init_groups;
    groups_cop = init_groups;
    
    acc_dpcp = [];
    acc_pca = [];
    acc_cop = [];
    acc_dpcp = [acc_dpcp utils.compute_accuracy(info.C, init_groups(1:length(info.C)))];
    acc_pca = [acc_pca utils.compute_accuracy(info.C, init_groups(1:length(info.C)))];
    acc_cop = [acc_cop utils.compute_accuracy(info.C, init_groups(1:length(info.C)))];

    iter = 0;
    while iter < maxiter
        iter = iter + 1;

        % update groups and normals
        for i = 1 : K
            tt = tic;
            mask = groups_dpcp == i;
            if ~any(mask)
                normals_dpcp(:,i) = normc(randn(D, 1));
            else
                [~, normal] = solvers.DPCP.rsgm(info.Xtilde(:,mask));        
                normals_dpcp(:,i) = normal;
            end
            t_dpcp = t_dpcp + toc(tt);
            
            tt = tic;
            mask = groups_pca == i;
            if ~any(mask)
                normals_pca(:,i) = normc(randn(D, 1));
            else
                [U,~] = svd(info.Xtilde(:,mask), 'econ');
                normals_pca(:,i) = U(:, end);
            end
            t_pca = t_pca + toc(tt);
            
            tt = tic;
            mask = groups_cop == i;
            if ~any(mask)
                normals_cop(:,i) = normc(randn(D, 1));
            else
                nSamp = max(D-1,size(info.Xtilde(:,mask),2)-2);
                U = solvers.CoP.CoP(info.Xtilde(:,mask),D-1,nSamp,2);        
                normals_cop(:,i) = null(U');
            end
            t_cop = t_cop + toc(tt);
        end

        % update the clustering
        [~,groups_dpcp] = min(abs(normals_dpcp'*info.Xtilde));
        [~,groups_pca] = min(abs(normals_pca'*info.Xtilde));
        [~,groups_cop] = min(abs(normals_cop'*info.Xtilde));

        acc_dpcp = [acc_dpcp utils.compute_accuracy(info.C, groups_dpcp(1:length(info.C)))];
        acc_pca = [acc_pca utils.compute_accuracy(info.C, groups_pca(1:length(info.C)))];
        acc_cop = [acc_cop utils.compute_accuracy(info.C, groups_cop(1:length(info.C)))];
    end
    
    if num_datasets == 1
        acc_list_dpcp = acc_dpcp;
        acc_list_pca = acc_pca;
        acc_list_cop = acc_cop;
    else
        acc_list_dpcp = acc_list_dpcp + acc_dpcp;
        acc_list_pca = acc_list_pca + acc_pca;
        acc_list_cop = acc_list_cop + acc_cop;
    end
end

acc_list_dpcp = acc_list_dpcp / T;
acc_list_pca = acc_list_pca / T;
acc_list_cop = acc_list_cop / T;
t_dpcp = t_dpcp / T
t_pca = t_pca / T
t_cop = t_cop / T

plot(acc_list_dpcp, '-r', 'linewidth', 4)
hold on
plot(acc_list_cop, '--g', 'linewidth', 4)
plot(acc_list_pca, '-.b', 'linewidth', 4)

legend('DPCP-KSS', 'CoP-KSS', 'PCA-KSS', 'location','best')
xlabel('iteration')
ylabel('accuracy')
xlim([0 maxiter])
ylim([0 1])

set(gca,'YDir','normal')
set(gca, ...
    'LineWidth' , 2                     , ...
    'FontSize'  , 30              , ...
    'FontName'  , 'Times New Roman'     , ...
    'XTick', [0:50:maxiter]           );

