function [info] = info_compute(parms)

%% basic info
info.N       =  50*(parms.D-1)*parms.K;    % total number of inliers
if parms.alpha == 1
    info.N1  =  ceil(info.N / parms.K);
else
    info.N1  =  ceil(info.N*(1-parms.alpha)/(1-parms.alpha^parms.K));
end

info.b1      =  normc(randn(parms.D,1));     
info.M       =  ceil(parms.r * info.N / (1 - parms.r));
info.C       =  ones(1,info.N1);

for i = 2 : parms.K
    eval(['info.N' num2str(i) '= ceil(parms.alpha*info.N' num2str(i-1) ');']);
    eval(['info.b' num2str(i) '= normc(randn(parms.D,1));']);
    eval(['info.C = [info.C ' num2str(i) '*ones(1,info.N' num2str(i) ')];']);
end

for i = 1:parms.K
    eval(['info.U' num2str(i) ' = null(info.b' num2str(i) ''');']);    % orthonormal basis
    eval(['info.P' num2str(i) ' = info.U' num2str(i) '*info.U' num2str(i) ''';']);   % projection
    eval(['info.X' num2str(i) ' = normc(info.P' num2str(i) '*randn(parms.D, info.N' num2str(i) ')' ...
        '+ (eye(parms.D)-info.P' num2str(i) ')*parms.sigma*randn(parms.D, info.N' num2str(i) '));']);  % inlier points
end

%% dataset creation
assign_data  =  'info.Xtilde = [';
for i = 1:parms.K
    assign_data = strcat(assign_data, ['info.X' num2str(i) ', ']);
end
assign_data = strcat(assign_data, 'info.O];');
info.O = normc(randn(parms.D, info.M));    % outliers
eval(assign_data);

return