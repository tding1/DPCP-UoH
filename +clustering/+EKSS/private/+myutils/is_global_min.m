function [out] = is_global_min(X, b)
    out = 1;
    val = norm(X'*b, 1);
    disp('---')
    disp(val)
    D = length(b);
    t = inf;
    for x = linspace(0,1,200)
        for y = linspace(0,1,200)
            for z = linspace(0,1,200)
                bb = [x, y ,z]';
                bb = bb/norm(bb);
                v = norm(X'*bb, 1);
                if v < t
                    t = v; 
                end
                if v < val
                    out = 0;
                    break;
                end
            end
        end
    end
    disp(t)
%     for n = 1:1000
%         bb = randn(D,1);
%         bb = bb/norm(bb);
%         v = norm(X'*bb, 1);
%         if v < val
%             out = 0;
%             break;
%         end
%     end
end