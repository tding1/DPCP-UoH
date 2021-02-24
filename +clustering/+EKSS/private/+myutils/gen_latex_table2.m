function gen_latex_table2(filename, out)

if nargin < 2
    out =  'output/tex.txt';
end

fid = fopen(out, 'w');

% write table heading
fprintf(fid, '%s\n', '\begin{table}[]');
fprintf(fid, '%s\n', '\tiny\centering');
fprintf(fid,'%s\n', '\begin{tabular}{c|c|c|cccccccc}');
fprintf(fid,'%s\n', '\hline');
fprintf(fid,'%s\n', '\multicolumn{3}{c|}{Settings}  & \multirow{2}{*}{-rsgm} & \multirow{2}{*}{-d} & \multirow{2}{*}{-irls} & \multirow{2}{*}{\begin{tabular}[c]{@{}c@{}}ransac\\ ($ 10^{-1} $)\end{tabular}} & \multirow{2}{*}{\begin{tabular}[c]{@{}c@{}}ransac\\ ($ 10^{-2} $)\end{tabular}} & \multirow{2}{*}{\begin{tabular}[c]{@{}c@{}}ransac\\ ($ 10^{-3} $)\end{tabular}} & \multirow{2}{*}{\begin{tabular}[c]{@{}c@{}}ransac10\\ ($ 10^{-2} $)\end{tabular}} & \multirow{2}{*}{\begin{tabular}[c]{@{}c@{}}ransac100\\ ($ 10^{-2} $)\end{tabular}} \\ \cline{1-3}');
fprintf(fid,'%s\n', '$ D $  & $ n $   & $ \alpha $     &     &      &      &    &       &      &       &    \\ \hline');

% load data
load([strcat(filename, '.mat')]);

% read data
methods = {'dpcp_rsgm', 'dpcp_d', 'dpcp_irls', ...
           'ransac_0_1', 'ransac_0_01', 'ransac_0_001', ...
           'ransac10_0_01', 'ransac100_0_01'};

D = {'_D4', '_D9', '_D30'};
n = {'_n2', '_n3', '_n4', '_n5', '_n6'};
alp = {'_alp6', '_alp8', '_alp10'};

data = -1 * ones(45, 8);
i = 1;
for di = 1:3
    for ni = 1:5
        for ai = 1:3
            Dna = [D{di}, n{ni}, alp{ai}];
            for mi = 1:8
                eval(['data(i, mi) = getfield(', filename, ',"', methods{mi}, Dna, '");']);
            end
            i = i + 1;
        end
    end
end

mx = max(data, [], 2);
mx_ind = data == mx;
tmp = data;
tmp(mx_ind) = -1;
sec_mx = max(tmp, [], 2);
sec_mx_ind = tmp == sec_mx;
sec_mx_ind(mx_ind) = 0;

mask = zeros(45,8);
mask(mx_ind) = 1;
mask(sec_mx_ind) = 2;

% write latex code
data_str = {};
for i = 1:45
    str = '';
    for j = 1:8
        if data(i,j) < 0
            acc = '-';
        else
            acc = round(100*data(i,j), 2);
        end
        if mask(i,j) == 1
            str = [str, '{\color{red}', string(acc), '}'];
        elseif mask(i,j) == 2
            str = [str, '{\color{blue}', string(acc), '}'];
        else
            str = [str, string(acc)];
        end
        if j ~= 8
            str = [str, '  &  '];
        end
    end
    data_str{i} = str;
end

num_best = sum(mask==1, 1);
num_sec_best = sum(mask==2, 1);
num_sum = num_best + num_sec_best;
num = [num_best; num_sec_best; num_sum];
mx = max(num, [], 2);
mx_ind = num == mx;
tmp = num;
tmp(mx_ind) = -1;
sec_mx = max(tmp, [], 2);
sec_mx_ind = tmp == sec_mx;

mask = zeros(3,8);
mask(mx_ind) = 1;
mask(sec_mx_ind) = 2;

for i = [46, 47, 48]
    str = '';
    for j = 1:8
        if mask(i-45,j) == 1
            str = [str, '{\color{red}', string(num(i-45, j)), '}'];
        elseif mask(i-45,j) == 2
            str = [str, '{\color{blue}', string(num(i-45, j)), '}'];
        else
            str = [str, string(num(i-45, j))];
        end
        if j ~= 8
            str = [str, '  &  '];
        end
    end
    data_str{i} = str;
end


fprintf(fid, '\\multirow{15}{*}{4}  & \\multirow{3}{*}{2} & 0.6 & %s  \\\\ \n', sprintf('%s',data_str{1}));
fprintf(fid, '&    & 0.8 &  %s  \\\\ \n', sprintf('%s',data_str{2}));
fprintf(fid, '&    & 1.0 &   %s  \\\\ \\cline{2-11} \n', sprintf('%s',data_str{3}));
fprintf(fid, '& \\multirow{3}{*}{3} & 0.6 & %s  \\\\ \n', sprintf('%s',data_str{4}));
fprintf(fid, '&    & 0.8 & %s  \\\\ \n', sprintf('%s',data_str{5}));
fprintf(fid, '&    & 1.0 &   %s  \\\\ \\cline{2-11} \n', sprintf('%s',data_str{6}));
fprintf(fid, '& \\multirow{3}{*}{4} & 0.6 &   %s  \\\\ \n', sprintf('%s',data_str{7}));
fprintf(fid, '&    & 0.8 & %s  \\\\ \n', sprintf('%s',data_str{8}));
fprintf(fid, '&    & 1.0 & %s  \\\\ \\cline{2-11}\n', sprintf('%s',data_str{9}));
fprintf(fid, '& \\multirow{3}{*}{5} & 0.6 &   %s  \\\\ \n', sprintf('%s',data_str{10}));
fprintf(fid, '&    & 0.8 & %s  \\\\ \n', sprintf('%s',data_str{11}));
fprintf(fid, '&    & 1.0 & %s  \\\\ \\cline{2-11}\n', sprintf('%s',data_str{12}));
fprintf(fid, '& \\multirow{3}{*}{6} & 0.6 &   %s  \\\\ \n', sprintf('%s',data_str{13}));
fprintf(fid, '&    & 0.8 & %s  \\\\ \n', sprintf('%s',data_str{14}));
fprintf(fid, '&    & 1.0 & %s  \\\\ \\hline\n', sprintf('%s',data_str{15}));
fprintf(fid, '\\multirow{15}{*}{9}  & \\multirow{3}{*}{2} & 0.6 &   %s  \\\\ \n', sprintf('%s',data_str{16}));
fprintf(fid, '&    & 0.8 &  %s  \\\\ \n', sprintf('%s',data_str{17}));
fprintf(fid, '&    & 1.0 &  %s  \\\\ \\cline{2-11} \n', sprintf('%s',data_str{18}));
fprintf(fid, '& \\multirow{3}{*}{3} & 0.6 &    %s  \\\\ \n', sprintf('%s',data_str{19}));
fprintf(fid, '&    & 0.8 &  %s  \\\\ \n', sprintf('%s',data_str{20}));
fprintf(fid, '&    & 1.0 &   %s  \\\\ \\cline{2-11} \n', sprintf('%s',data_str{21}));
fprintf(fid, '& \\multirow{3}{*}{4} & 0.6 & %s  \\\\ \n', sprintf('%s',data_str{22}));
fprintf(fid, '&    & 0.8 &  %s  \\\\ \n', sprintf('%s',data_str{23}));
fprintf(fid, '&    & 1.0 &   %s  \\\\ \\cline{2-11} \n', sprintf('%s',data_str{24}));
fprintf(fid, '& \\multirow{3}{*}{5} & 0.6 &   %s  \\\\ \n', sprintf('%s',data_str{25}));
fprintf(fid, '&    & 0.8 & %s  \\\\ \n', sprintf('%s',data_str{26}));
fprintf(fid, '&    & 1.0 & %s  \\\\ \\cline{2-11}\n', sprintf('%s',data_str{27}));
fprintf(fid, '& \\multirow{3}{*}{6} & 0.6 &   %s  \\\\ \n', sprintf('%s',data_str{28}));
fprintf(fid, '&    & 0.8 & %s  \\\\ \n', sprintf('%s',data_str{29}));
fprintf(fid, '&    & 1.0 & %s  \\\\ \\hline\n', sprintf('%s',data_str{30}));
fprintf(fid, '\\multirow{15}{*}{30} & \\multirow{3}{*}{2} & 0.6 & %s  \\\\ \n', sprintf('%s',data_str{31}));
fprintf(fid, '&    & 0.8 &  %s  \\\\ \n', sprintf('%s',data_str{32}));
fprintf(fid, '&    & 1.0 &   %s  \\\\ \\cline{2-11} \n', sprintf('%s',data_str{33}));
fprintf(fid, '& \\multirow{3}{*}{3} & 0.6 & %s  \\\\ \n', sprintf('%s',data_str{34}));
fprintf(fid, '&    & 0.8 &  %s  \\\\ \n', sprintf('%s',data_str{35}));
fprintf(fid, '&    & 1.0 &   %s  \\\\ \\cline{2-11} \n', sprintf('%s',data_str{36}));
fprintf(fid, '& \\multirow{3}{*}{4} & 0.6 &  %s  \\\\ \n', sprintf('%s',data_str{37}));
fprintf(fid, '&    & 0.8 &  %s  \\\\ \n', sprintf('%s',data_str{38}));
fprintf(fid, '&    & 1.0 &   %s  \\\\ \\cline{2-11} \n', sprintf('%s',data_str{39}));
fprintf(fid, '& \\multirow{3}{*}{5} & 0.6 &   %s  \\\\ \n', sprintf('%s',data_str{40}));
fprintf(fid, '&    & 0.8 & %s  \\\\ \n', sprintf('%s',data_str{41}));
fprintf(fid, '&    & 1.0 & %s  \\\\ \\cline{2-11}\n', sprintf('%s',data_str{42}));
fprintf(fid, '& \\multirow{3}{*}{6} & 0.6 &   %s  \\\\ \n', sprintf('%s',data_str{43}));
fprintf(fid, '&    & 0.8 & %s  \\\\ \n', sprintf('%s',data_str{44}));
fprintf(fid, '&    & 1.0 & %s  \\\\ \\hline\n', sprintf('%s',data_str{45}));
fprintf(fid, '\\multicolumn{3}{c|}{best}   &  %s  \\\\ \\hline \n', sprintf('%s',data_str{46}));
fprintf(fid, '\\multicolumn{3}{c|}{second best}   &  %s  \\\\ \\hline \n', sprintf('%s',data_str{47}));
fprintf(fid, '\\multicolumn{3}{c|}{sum}   &  %s  \\\\ \\hline \n', sprintf('%s',data_str{48}));

fprintf(fid,'%s\n', '\end{tabular}');
fprintf(fid,'%s\n', '\end{table}');

fclose(fid);

end
