function [groups, normals] = CoRe(Xtilde, r0, meta_normals, meta_groups, T, delta)

    num_replica = length(meta_normals);
    sibling_index = 1:num_replica;
    sibling_index(r0) = [];
    
    r0_normals = meta_normals{r0};
    r0_groups = meta_groups{r0};
    K = size(r0_normals, 2);
    for s = 1:T
        most_decrease = -Inf;
        obj_old = loss(Xtilde, r0_groups, r0_normals);
        for i = sibling_index
            ri_normals = meta_normals{i};
            for j = 1:K
                tmp_normals = [r0_normals, ri_normals(:, j)];
                dist = abs(tmp_normals' * Xtilde);
                [~, tmp_groups] = min(dist, [], 1);
                obj_new = loss(Xtilde, tmp_groups, tmp_normals);
                decrease = obj_old - obj_new;
                if decrease > most_decrease
                    most_decrease = decrease;
                    adding_index = [i, j];
                    obj_mid = obj_new;
                end
            end
        end
        
        adding_normals = meta_normals{adding_index(1)};
        adding_normal = adding_normals(:, adding_index(2));
        tmp_normals = [r0_normals, adding_normal];
        
        least_increase = Inf;
        for j = 1:K
            tmp_normals_hat = tmp_normals;
            tmp_normals_hat(:, j) = [];
            dist = abs(tmp_normals_hat' * Xtilde);
            [~, tmp_groups] = min(dist, [], 1);
            obj_new = loss(Xtilde, tmp_groups, tmp_normals_hat);
            increase = obj_new - obj_mid;
            if increase < least_increase
                least_increase = increase;
                dropping_index = [r0, j];
            end
        end
        
        normals = tmp_normals;
        normals(:, dropping_index(2)) = [];
        dist = abs(normals' * Xtilde);
        [~, groups] = min(dist, [], 1);
        obj = loss(Xtilde, groups, normals);

        if obj > (1-delta) * obj_old
            normals = r0_normals;
            groups = r0_groups;
        end
    end

end


function val = loss(Xtilde, groups, normals)
    val = 0;
    n = max(groups);
    for i = 1 : n
        mask = groups==i;
        val = val + norm(Xtilde(:,mask)'*normals(:,i), 1);
    end
end

