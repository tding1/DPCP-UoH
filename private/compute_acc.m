function [acc, l] = compute_acc(meta_groups, meta_normals, info)
l = Inf;
for j = 1:length(meta_groups)
    tmp_l = loss(info.Xtilde, meta_groups{j}, meta_normals{j});
    if tmp_l < l
        l = tmp_l;
        groups = meta_groups{j};
    end
end

acc = utils.compute_accuracy(info.C, groups(1:length(info.C)));

end

function val = loss(Xtilde, groups, normals)
    val = 0;
    n = max(groups);
    for i = 1 : n
        mask = groups==i;
        val = val + norm(Xtilde(:,mask)'*normals(:,i), 1);
    end
end
