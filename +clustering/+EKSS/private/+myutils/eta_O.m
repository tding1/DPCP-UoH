function out = eta_O(X1, X2)
% Compute eta_O quantity only when n = 3.
    
    s = 0;
    D = size(X1,1);
    for i = 1:1e4
        b = randn(D,1);
        b = b / norm(b);
        t = norm((eye(D)-b*b') * (X1*sign(X1'*b)+X2*sign(X2'*b)));
        if t > s
            s = t;
        end
    end
    out = s;

end

