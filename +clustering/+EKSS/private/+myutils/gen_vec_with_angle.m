function [vec] = gen_vec_with_angle(u, inner_prod)
% generate an unit vector such that it has the specified inner product with u

% randomly initialize an another vector v to form a plane with u
v = randn(size(u));

u = u / norm(u);
v = v / norm(v);

% suppose vec = u + y * v, solve for y
d = dot(u, v);
t = 1-inner_prod^2;
coeff1 = d^2 - inner_prod^2;
coeff2 = 2 * t * d;
coeff3 = t;

% solve a quadratic equation for y
rts = roots([coeff1, coeff2, coeff3]);

% two roots correspond to different directions of vec
% we want the one with exactly the same inner product (including the sign)
for i = 1:2
    y = rts(i);
    vec = u + y * v;
    vec = vec / norm(vec);
    vec = real(vec);
    diff = abs(dot(vec, u) - inner_prod);
    if diff < 1e-5
        break;
    end
end

assert(diff < 1e-5)

end