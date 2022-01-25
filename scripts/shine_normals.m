n = 16;
Gamma = Boundary.circle(n, 'quadrature', 'panel');
%Gamma = refine(Gamma, 2);
%Gamma = Boundary.star(n, 'quadrature', 'panel');
%Gamma = refine(Gamma);

n_re = n;
%n_re = n;
Gamma_re = resample(Gamma, n_re);

% Define Gamma'
beta = 4;
bleed = 10;
width = 0.5;
[x, y] = smoothStrip2(Gamma_re, n_re, beta, bleed, width);
xleg = chebvals2legvals(reshape(x, n_re, Gamma_re.np)); xleg = xleg(:);
yleg = chebvals2legvals(reshape(y, n_re, Gamma_re.np)); yleg = yleg(:);
z1 = mat2cell(xleg + 1i*yleg, repmat(n_re, Gamma_re.np, 1), 1);
Gamma1 = Boundary(z1);

% Build the strip grid
xsem = chebpts(n_re);
t = chebpts(n_re, [0 1]);
xx = cell(Gamma_re.np,1);
yy = cell(Gamma_re.np,1);
gamx_cheb  = cell(Gamma_re.np, 1);
gamy_cheb  = cell(Gamma_re.np, 1);
gam1x_cheb = cell(Gamma_re.np, 1);
gam1y_cheb = cell(Gamma_re.np, 1);
gam_nx_cheb = cell(Gamma_re.np, 1);
gam1_nx_cheb = cell(Gamma_re.np, 1);
for k = 1:Gamma_re.np
    gamx_cheb{k}  = legvals2chebvals(real(Gamma_re.x{k}(:,1)));
    gamy_cheb{k}  = legvals2chebvals(real(Gamma_re.x{k}(:,2)));
    gam1x_cheb{k} = legvals2chebvals(real(Gamma1.x{k}(:,1)));
    gam1y_cheb{k} = legvals2chebvals(real(Gamma1.x{k}(:,2)));
    gam_nx_cheb{k} = legvals2chebvals(Gamma_re.normal{k});
    gam1_nx_cheb{k} = legvals2chebvals(Gamma1.normal{k});
    % Upsample
    xx{k} = t.*gam1x_cheb{k}.' + (1-t).*gamx_cheb{k}.';
    yy{k} = t.*gam1y_cheb{k}.' + (1-t).*gamy_cheb{k}.';
end
dom = cell2struct([xx yy], {'x','y'}, 2);

%%
e = 6;

surf(dom(e).x, dom(e).y, 0*dom(e).x, 'FaceAlpha', 0)
view(2)
axis equal
hold on

xfun = chebfun2(dom(e).x);
yfun = chebfun2(dom(e).y);
ss = 2*rand(1) - 1;
rr = 2*rand(1) - 1;
ss = 1;
rr = 0;
%q = [0.91 0.14];
q = [-0.51 0.7];
%q = [0.58 0.35];
%q = [-0.15 0.45];
%q = [0.05 -0.1];
% q = [0 0.9];
%q = [xfun(ss,rr) yfun(ss,rr)];
%q = [0.696773816326466 -0.464515877550977];
%q = [q; q];
plot(q(1), q(2), 'b.', 'markersize', 30)

kd = KDTree([xx{e}(:) yy{e}(:)]);
nearestIdx = kd.nn(q);
q_nearest = [xx{e}(nearestIdx) yy{e}(nearestIdx)];
%plot(q_nearest(1), q_nearest(2), 'bo', 'markersize', 10)

[i, j] = ind2sub([n_re n_re], nearestIdx);
%nx_cheb = legvals2chebvals(Gamma.normal{e});
nx_j = gam_nx_cheb{e}(j,:);
perpx = -nx_j(2);
perpy =  nx_j(1);
%perp = rand(2,1); perp = perp ./ norm(perp);
%perpx = perp(1);
%perpy = perp(2);
Lx_q = @(t) q(1) + t*perpx;
Ly_q = @(t) q(2) + t*perpy;
tt = linspace(-1,1,100);
plot(Lx_q(tt), Ly_q(tt), 'k--', 'linewidth', 1)

tic
tt = linspace(-1,2,100);
tk = zeros(n_re, 1);
for k = 1:n_re
    sx = gam_nx_cheb{e}(k, 1);
    sy = gam_nx_cheb{e}(k, 2);
    %sx = gamx_cheb{e}(k)-gam1x_cheb{e}(k);
    %sy = gamy_cheb{e}(k)-gam1y_cheb{e}(k);
    Lx_k = @(t) gam1x_cheb{e}(k) + sx*t;
    Ly_k = @(t) gam1y_cheb{e}(k) + sy*t;
    plot(Lx_k(tt), Ly_k(tt), 'k--', 'linewidth', 1)
    A = [perpx -sx;
         perpy -sy];
    t_intersect = A \ [gam1x_cheb{e}(k) - q(1);
                       gam1y_cheb{e}(k) - q(2)];
    xk = Lx_k(t_intersect(2));
    yk = Ly_k(t_intersect(2));
    plot(xk, yk, 'r.', 'linewidth', 2, 'markersize', 30)
    tk(k) = norm(q - [xk yk], 2);
    tk(k) = tk(k) * sign(dot(q - [perpx perpy], q - [xk yk]));
end
toc

tic
sx = reshape(gam_nx_cheb{e}(:,1), 1, 1, []);
sy = reshape(gam_nx_cheb{e}(:,2), 1, 1, []);
e3 = ones(1, 1, n_re); % Tube vector
invA = [-sy sx; -perpy*e3 perpx*e3] ./ (perpy*sx - perpx*sy);
rhs = reshape(([gam1x_cheb{e} gam1y_cheb{e}] - q).', 2, 1, []);
t_intersect = pagemtimes(invA, rhs);
xk = gam1x_cheb{e} + squeeze(sx .* t_intersect(2,:,:));
yk = gam1y_cheb{e} + squeeze(sy .* t_intersect(2,:,:));
tk2 = vecnorm(q - [xk yk], 2, 2);
tk2 = tk2 .* sign((q - [xk yk]) * (q - [perpx perpy]).');
toc

[xcheb, ~, vcheb] = chebpts(n_re);
s_q = util.bary(0, xcheb, tk, vcheb);
gamma_q  = util.bary(s_q, [gamx_cheb{e}  gamy_cheb{e}], xcheb, vcheb);
gamma1_q = util.bary(s_q, [gam1x_cheb{e} gam1y_cheb{e}], xcheb, vcheb);
r_q = 2*norm(q - gamma_q) ./ norm(gamma1_q - gamma_q) - 1;

plot(gamma_q(1),  gamma_q(2),  'm.', 'markersize', 30)
plot(gamma1_q(1), gamma1_q(2), 'm.', 'markersize', 30)
%plot(xfun(s_q,r_q), yfun(s_q,r_q), 'bo')

% max(abs(s_q - ss), [], 'all')
% max(abs(r_q - rr), [], 'all')

max(abs([xfun(s_q,r_q) yfun(s_q,r_q)] - q), [], 'all')

%nxfun = chebfun(nx_cheb);
shg

axis([minmax(xx{e}) minmax(yy{e})])
box on
set(gca, 'xtick', [], 'ytick', [])