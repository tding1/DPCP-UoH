function [vec] = gen_vec_with_angle3(u1, u2, inner_prod1, inner_prod2)
% generate a unit vector such that it has the specified inner product with u1,
% and the specified inner product with u2
% note: u1 and u2 can be any dimension

u1 = u1 / norm(u1);
u2 = u2 / norm(u2);

A = [u1'; u2'];
assert(size(A,1) == 2)
b = [inner_prod1; inner_prod2];
t = A\b;

N = null(A);
n = N(:,1);

coeff1 = 1;
coeff2 = 2*n'*t;
coeff3 = t'*t-1;

rts = roots([coeff1, coeff2, coeff3]);

for i = 1:2
    a = rts(i);
    vec = t + a*n;
    vec = real(vec);
    diff1 = abs(dot(vec, u1) - inner_prod1);
    diff2 = abs(dot(vec, u2) - inner_prod2);
    if diff1 < 1e-5 && diff2 < 1e-5
        break;
    end
end

assert(diff1 < 1e-5 && diff2 < 1e-5)

end


