clc
clear
close all

files = {'timing_outliers_0_noise_0','timing_outliers_1_noise_0', 'timing_outliers_1_noise_1'};

methods = {'dpcp_rsgm', 'dpcp_d', 'dpcp_irls', ...
           'ransac_0_1', 'ransac_0_01', 'ransac_0_001', ...
           'ransac10_0_01', 'ransac100_0_01', 'fsasc'};

settings1 = {'_D4_n2_alp10', '_D4_n3_alp10', '_D4_n4_alp10'};
settings2 = {'_D9_n2_alp10', '_D9_n3_alp10', '_D9_n4_alp10'};
settings3 = {'_D30_n2_alp10', '_D30_n3_alp10', '_D30_n4_alp10'};

num_settings = length(settings1);

data1 = zeros(num_settings, 9);
data2 = zeros(num_settings, 9);
data3 = zeros(num_settings, 9);

lbls = cell(1,9);
lbls{1}='-rsgm'; lbls{2}='-d'; lbls{3}='-irls'; 
lbls{4}='ransac (10^{-1})'; lbls{5}='ransac (10^{-2})'; lbls{6}='ransac (10^{-3})'; 
lbls{7}='ransac10 (10^{-2})'; lbls{8}='ransac100 (10^{-2})'; lbls{9}='FSASC'; 

for k = 1:3
    file = files{k};
    
    load(['saved_data/', strcat(file, '.mat')]);

    for j = 1:9
        for i = 1:num_settings
            eval(['data1(i, j) = getfield(', file, ',"', methods{j}, settings1{i}, '");'])
            eval(['data2(i, j) = getfield(', file, ',"', methods{j}, settings2{i}, '");']) 
            eval(['data3(i, j) = getfield(', file, ',"', methods{j}, settings3{i}, '");']) 
        end
    end

    figure
    h = bar(data1);
    colors = parula(9);
    for k = 1:length(h)
        h(k).FaceColor = colors(10-k, :);
    end
    legend(h,lbls, 'fontsize', 18, 'Location', 'northwest');
    grid on
    set(gca, ...
    'xticklabel', {'n = 2', 'n = 3', 'n = 4'}, ...
    'FontSize'  , 20             , ...
    'FontName'  , 'Times New Roman');
    ylabel('Running Time (s)');
    xlabel('D = 4, \alpha = 1.0');

    
    figure
    h = bar(data2);
    colors = parula(9);
    for k = 1:length(h)
        h(k).FaceColor = colors(10-k, :);
    end
    legend(h,lbls, 'fontsize', 18, 'Location', 'northwest');
    grid on
    set(gca, ...
    'xticklabel', {'n = 2', 'n = 3', 'n = 4'}, ...
    'FontSize'  , 20             , ...
    'FontName'  , 'Times New Roman');
    ylabel('Running Time (s)');
    xlabel('D = 9, \alpha = 1.0');  
    
    figure
    h = bar(data3(:, 1:end-1));
    colors = parula(9);
    for k = 1:length(h)
        h(k).FaceColor = colors(9-k+1, :);
    end
    legend(h,lbls{:,1:end-1}, 'fontsize', 18, 'Location', 'northwest');
    grid on
    set(gca, ...
    'xticklabel', {'n = 2', 'n = 3', 'n = 4'}, ...
    'FontSize'  , 20             , ...
    'FontName'  , 'Times New Roman');
    ylabel('Running Time (s)');
    xlabel('D = 30, \alpha = 1.0');  
    
end



       
       
       