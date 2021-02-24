function [vec] = gen_vec_with_angle2(u1, u2, inner_prod1, inner_prod2)
% generate a 3d unit vector such that it has the specified inner product with u1 (3d),
% and the specified inner product with u2 (3d)

u1 = u1 / norm(u1);
u2 = u2 / norm(u2);

a = u1(1); b = u1(2); c = u1(3);
d = u2(1); e = u2(2); f = u2(3);

det = b*f-c*e;
s1 = (b*inner_prod2-e*inner_prod1) / det;
s2 = (a*e-b*d) / det;
r1 = (f*inner_prod1-c*inner_prod2) / det;
r2 = (c*d-a*f) / det;

coeff1 = 1+r2^2+s2^2;
coeff2 = 2*(r1*r2+s1*s2);
coeff3 = r1^2+s1^2-1;

rts = roots([coeff1, coeff2, coeff3]);

for i = 1:2
    x = rts(i);
    y = r1+r2*x;
    z = s1+s2*x;
    vec = [x,y,z]';
    vec = real(vec);
    diff1 = abs(dot(vec, u1) - inner_prod1);
    diff2 = abs(dot(vec, u2) - inner_prod2);
    if diff1 < 1e-5 && diff2 < 1e-5
        break;
    end
end

assert(diff1 < 1e-5 && diff2 < 1e-5)

end

