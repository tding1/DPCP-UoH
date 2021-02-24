function val = loss(Xtilde, groups, normals)
    val = 0;
    n = max(groups);
    for i = 1 : n
        mask = groups==i;
        val = val + norm(Xtilde(:,mask)'*normals(:,i), 1);
    end
end