n = 16;
x1 = chebpts(n, [0 2]);
x2 = chebpts(n, [2 3]);
x3 = chebpts(n, [3 3.5]);
x4 = chebpts(n, [3.5 4.5]);
x5 = chebpts(n, [4.5 6.5]);
%x6 = chebpts(n, [6.5 7.5]);
plot(x1, 0*x1, 'ko', x2, 0*x2, 'ko', x3, 0*x3, 'ko', x4, 0*x4, 'ko', ...
     x5, 0*x5, 'ko')%, x6, 0*x6, 'ko')

dist = zeros(5*n, 1);
pts = [x1; x2; x3; x4; x5];
pts_wrapped = [x5-6.5; x1; x2; x3; x4; x5; x1+6.5];
for k = 1:5*n
    dist(k) = mean([abs(pts_wrapped(n+k+n) - pts_wrapped(n+k)), ...
                    abs(pts_wrapped(n+k-n) - pts_wrapped(n+k))]);
end

plot(pts(1:5*n), dist, 'o')
shg

distf = chebfun(mat2cell(dist, repmat(n, 5, 1), 1), [0 2 3 3.5 4.5 6.5]);

%%
blend_width = 12*n;
[step, bump] = util.makeMask(1000, blend_width);
%phi = chebfun(@(x) bump(x/6.5), [0 6.5]);
phi = chebfun(@(x) exp(-20*(x-a).^2), [0 6.5]);
a = 2;
f = distf;
g = @(x) phi(x).*f(a) + (1-phi(x)).*f(x);
plot(pts(1:5*n), g(pts(1:5*n)))
hold on
plot(pts(1:5*n), distf(pts(1:5*n)))
hold off
shg