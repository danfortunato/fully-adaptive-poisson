function [u, f] = poissonPair()

K = 1e4; % steepness in gaussian test function

% centers of the two gaussians
center_one = [0.3, 0.25];
center_two = [0.8, 0.8];

g1 = @(x,y) exp(-100000*((x-0.2).^2 + (y-0.3).^2));
g = cheb.gallery2('tiltedpeg');
g1 = @(x,y) g(2*x-1, 2*y-1);
ufun = chebfun2(g1, [0 1 0 1]);
ffun = -lap(ufun);

    function u = u_func(x, y)
        u = ufun(x, y);
%         r2a = (x-center_one(1)).^2 + (y-center_one(2)).^2;
%         r2b = (x-center_two(1)).^2 + (y-center_two(2)).^2;
%         ea = 3*exp(-K*r2a);
%         eb = 3*exp(-K*r2b);
%         u = ea + eb + 2*x.^3.*exp(y);
    end

    function f = f_func(x, y)
        fprintf("evaluating f\n")
        f = ffun(x, y);
%         r2a = (x-center_one(1)).^2 + (y-center_one(2)).^2;
%         r2b = (x-center_two(1)).^2 + (y-center_two(2)).^2;
%         ea = 3*exp(-K*r2a);
%         eb = 3*exp(-K*r2b);
%         fa = 4*K*ea.*(K*r2a - 1);
%         fb = 4*K*eb.*(K*r2b - 1);
%         fc = 2*x.*exp(y).*(6 + x.^2);
%         f = -(fa + fb + fc);
    end

u = @u_func;
f = @f_func;

end
