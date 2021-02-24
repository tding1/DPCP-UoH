clc
clear
close all
randn('seed',2020);rand('seed',2020)

D = 9;
K = 3;
parms = rand_parms_spec();
parms.D = D;
parms.K = K;
parms.alpha = 0.8;
parms.r = 0.3;

info = info_compute(parms);

maxiter = 2e2;
beta_list = [0.3  0.6  0.9];
mu_o = 1e-2;

Xtilde = info.Xtilde;
c = 1;
d = D-1;

% initialization
[U,~,~] = svd(Xtilde, 'econ');
Bo = U(:,end);
for i_beta = 1:length(beta_list)
    beta = beta_list(i_beta);
    B = Bo;
    i = 1;
    
    dist(1,i_beta) = min(norm(B-info.b1), norm(B+info.b1));
    while i<= maxiter
        i = i+1;
        mu = mu_o*beta^(i);
        grad = Xtilde*sign(Xtilde'*B);
        grad = grad - B*(grad'*B);
   
        tmp = B - mu*grad;
        B = tmp / norm(tmp);
        dist(i,i_beta) =  min(norm(B-info.b1), norm(B+info.b1));
    end
end

%%
fontsize = 30;
plotStyle = {'b-','r:','k--','g-','b:','r--','k-','k:','k--'};
figure
for i_beta = 1:length(beta_list)
    semilogy(0:length(find(dist(:,i_beta)>0))-1,dist(find(dist(:,i_beta)>0),i_beta),plotStyle{i_beta},'linewidth',3);
    legendInfo{i_beta} = ['\beta = ' num2str(beta_list(i_beta))];
    
    hold on
end
ylim([1e-12,1])
xlim([0 size(dist,1)])
legend(legendInfo,'Location','Best','fontsize', 30)
xlabel('iteration','FontSize',fontsize,'FontName','Times New Roman','Interpreter','LaTex');
y=ylabel('dist($\widehat{\bf b}_t, \{\pm\bf b_1\}$)','FontSize',fontsize,'FontName','Times New Roman','Interpreter','LaTex');
get(y,'position');
set(y, 'position', get(y,'position')-[10,0,0]);
set(gca,'YDir','normal')
set(gca, ...
    'LineWidth' , 2                     , ...
    'FontSize'  , fontsize-2              , ...
    'FontName'  , 'Times New Roman'         );
set(gcf, 'Color', 'white');
