%% An adaptive partition-of-unity strip solver III
%
%  Let Omega be a domain, Gamma be its boundary, B be a bounding box
%  containing Omega, S be a strip, and Gamma' be the strip boundary.
%
%  Solve the Poisson problem
%
%        lap(u) = f  in Omega
%             u = g  on Gamma
%
%  in four stages:
%
%    u = u_bulk + u_strip + u_glue + u_bc

%% The problem

n = 16;      % Number of nodes per panel
m = n;       % Number of radial Chebyshev nodes in strip
nbulk = 200; % Number of uniform nodes for bulk FFT solve

%Gamma = Boundary.star(n, 'quadrature', 'panel');
Gamma = Boundary.circle(n, 'quadrature', 'panel');
Gamma = refine(Gamma, 2);

gamxy = cell2mat(Gamma.x);
gamx  = gamxy(:,1);
gamy  = gamxy(:,2);

scale = 1.2;
dom_global = boundingbox(Gamma, scale);

sol = @(x,y) -cos(x).*exp(sin(x)).*sin(y);
f   = @(x,y) (2*cos(x)+3*cos(x).*sin(x)-cos(x).^3).*exp(sin(x)).*sin(y);
g   = @(x,y) sol(x,y);
% sol = randnfun2(1, dom_global, 'trig');
% f = lap(sol);
% g = @(x,y) sol(x,y);

%% Define the fake interior boundary, Gamma'
[x, y] = smoothStrip(Gamma, n, 0.8);
xleg = chebvals2legvals(reshape(x, n, Gamma.np)); xleg = xleg(:);
yleg = chebvals2legvals(reshape(y, n, Gamma.np)); yleg = yleg(:);
z1 = mat2cell(xleg + 1i*yleg, repmat(n, Gamma.np, 1), 1);
Gamma1 = Boundary(z1);

% Build the strip grid
cx = mat2cell(x, repmat(n,Gamma.np,1), 1);
cy = mat2cell(y, repmat(n,Gamma.np,1), 1);
t = chebpts(n, [0 1]);
xx = cell(Gamma.np,1);
yy = cell(Gamma.np,1);
gamx_cheb = cell(Gamma.np, 1);
gamy_cheb = cell(Gamma.np, 1);
gam1x_cheb = cell(Gamma.np, 1);
gam1y_cheb = cell(Gamma.np, 1);
for k = 1:Gamma.np
    gamx_cheb{k} = legvals2chebvals(real(Gamma.x{k}(:,1)));
    gamy_cheb{k} = legvals2chebvals(real(Gamma.x{k}(:,2)));
    gam1x_cheb{k} = legvals2chebvals(real(Gamma1.x{k}(:,1)));
    gam1y_cheb{k} = legvals2chebvals(real(Gamma1.x{k}(:,2)));
    xx{k} = (1-t).*gam1x_cheb{k}.' + t.*gamx_cheb{k}.';
    yy{k} = (1-t).*gam1y_cheb{k}.' + t.*gamy_cheb{k}.';
    %xx{k} = (1-t).*cx{k}.' + t.*gamx_cheb{k}.';
    %yy{k} = (1-t).*cy{k}.' + t.*gamy_cheb{k}.';
end
dom = cell2struct([xx yy], {'x','y'}, 2);

%% Get grids
tic

[xx, yy] = meshgrid(trigpts(nbulk, dom_global(1:2)), ...
                    trigpts(nbulk, dom_global(3:4)));
inGamma  = isinterior(Gamma,  xx, yy);
inGamma1 = isinterior(Gamma1, xx, yy);
inStrip  = inGamma & ~inGamma1;
inStripFunc = util.inStripFactory(Gamma, Gamma1);
%inStrip2 = reshape(inStripFunc([xx(:) yy(:)]), size(xx));

xx_g = xx(inGamma); xx_g1 = xx(inGamma1); xx_s = xx(inStrip);
yy_g = yy(inGamma); yy_g1 = yy(inGamma1); yy_s = yy(inStrip);

% Convert global grid coordinates to local element coordinates.
% This involves (1) determining which element each point is located in and
% (2) inverting the element's coordinate maps using Newton's method.
tt_s = cell(length(dom), 1);
rr_s = cell(length(dom), 1);
rr_s_global = zeros(size(xx_s));
found = false(size(xx_s));
globalIdx = cell(length(dom), 1);
for k = 1:length(dom)
    k
    xfun = chebfun2(dom(k).x);
    yfun = chebfun2(dom(k).y);
    [tt, rr, converged] = cart2curv_element(xfun, yfun, xx_s, yy_s);
    %idx = ~found & converged & abs(tt) <= 1+100*eps & abs(rr) <= 1+100*eps;
    idx = ~found & converged & abs(tt) <= 1.01 & abs(rr) <= 1.01;
    tt_s{k} = tt(idx);
    rr_s{k} = rr(idx);
    globalIdx{k} = idx;
    rr_s_global(idx) = rr_s{k};
    found = found | idx;
end

fprintf('\n');
fprintf('Grid computation took...%.6gs\n', toc);

% Did we miss any points?
idx = globalIdx{1};
for k = 2:length(dom)
    idx = idx | globalIdx{k};
end
if ( ~all(idx) )
    warning('Some grid points were not assigned to an element.');
    figure
    plot(Gamma), hold on, plot(Gamma1), hold on
    scatter(xx_s(idx), yy_s(idx), 'bo'), hold on
    scatter(xx_s(~idx), yy_s(~idx), 'ro'), hold off
end

% Plot the points for each element
plot(Gamma), hold on
plot(Gamma1)
for k = 1:length(dom)
    idx = globalIdx{k};
    scatter(xx_s(idx), yy_s(idx)), hold on
end
hold off
shg

%%

m = 16;
blend_width = 2*m;
[step01, bump] = util.makeMask(1000, blend_width);
step = chebfun(@(x) step01((1-x)/2));
ff_s = zeros(size(xx_s));
for k = 1:length(dom)
    idx = globalIdx{k};
    ff_s(idx) = step(rr_s{k});
    %scatter3(xx_s(idx), yy_s(idx), ff_s(idx), 100, ff_s(idx), '.')
    %hold on
end
%shg

ff = zeros(size(xx));
ff(inStrip) = ff_s;
ff(inGamma1) = 1;
figure(2)
surf(xx, yy, ff)
axis equal
shading interp

%%
tt = zeros(size(xx_s));
rr = zeros(size(xx_s));
xfun = chebfun2(dom(1).x);
yfun = chebfun2(dom(1).y);
for k = 1:numel(xx_s)
    k
    [tt(k), rr(k)] = cart2curv_element(xfun, yfun, xx_s(k), yy_s(k));
end

%% Treefun evaluation
sigma = ones(Gamma.N*Gamma.np, 1);
inGammaFunc   = @(xy) kernels.laplace.dlp(Gamma,  'target', xy, 'density', sigma, 'side', 'i');
outGammaFunc  = @(xy) kernels.laplace.dlp(Gamma,  'target', xy, 'density', sigma, 'side', 'e');
inGamma1Func  = @(xy) kernels.laplace.dlp(Gamma1, 'target', xy, 'density', sigma, 'side', 'i');
outGamma1Func = @(xy) kernels.laplace.dlp(Gamma1, 'target', xy, 'density', sigma, 'side', 'e');
%inStripFunc = @(xy) inGammaFunc(xy) & ~inGamma1Func(xy);
%inStripFunc = @(xy) inGammaFunc(xy) + inGamma1Func(xy);
%inStripFunc = @(xy) outGammaFunc(xy) + inGamma1Func(xy);
%inStripFunc = @(xy) abs(inGammaFunc(xy)+1) < 0.99 & abs(outGamma1Func(xy)) < 0.99;
%inStripFunc = @(xy) abs(outGammaFunc(xy)) > 0.5 & abs(inGamma1Func(xy)+1) > 0.5;
inStripFunc = util.inStripFactory(Gamma, Gamma1);
[xx, yy] = meshgrid(trigpts(nbulk, dom_global(1:2)), ...
                    trigpts(nbulk, dom_global(3:4)));
xy_s2 = inStripFunc([xx(:) yy(:)]);
idx = xy_s2;
%idx = abs(xy_s2 + 1) < 0.5;
scatter(xx(idx),yy(idx))
hold on
plot(Gamma)
plot(Gamma1)
shg

%% (1) Solve the bulk problem
%
%    lap(u_bulk) = phi*f     in B
%         u_bulk = periodic  on B
tic

ff = zeros(nbulk);
ff(inGamma) = f(xx_g, yy_g);

blend_width = 2*m;
[step, bump] = util.makeMask(1000, blend_width);
phi = zeros(nbulk);
phi(inStrip) = step((-rr_s_global+1)/2);
phi(inGamma1) = 1;
bumpf = @(x) (-1<=x&x<=1).*bump((x+1)/2);

% TODO: Make phi*f mean zero.
rhs = phi .* ff;

% unit integral bump function located at the corner of the rectangular domain
% this should not intersect the real domnain!
% easy to guarantee by increasing the size of the rectangular domain a bit
hx = diff(dom_global(1:2))/nbulk;
hy = diff(dom_global(3:4))/nbulk;
h = min(hx,hy);
loc = 0.3;
grr = sqrt((xx-(dom_global(2)-loc)).^2 + (yy-(dom_global(4)-loc)).^2);
bumpy = bumpf( grr./loc );
bumpy = bumpy ./ (sum(sum(bumpy))*h.^2);
rhs = rhs - bumpy * sum(sum(rhs)) * h.^2;
u_bulk = util.periodicPoisson( rhs, dom_global );

fprintf('FFT solve took..........%.6gs\n', toc);

%% (2) Solve the strip problem
%
%    lap(u_strip) = f  in S
%               u = 0  on Gamma & Gamma'
tic

S = StripSolver(dom, f);
build(S);
%bc = zeros(size(S.patches{1}.xy, 1), 1);
bc = u_bulk(S.patches{1}.xy(:,1), S.patches{1}.xy(:,2));
% bc = g(S.patches{1}.xy(:,1), S.patches{1}.xy(:,2));
innerIdx = reshape((1:n-2).' + (n-2)*(0:2:2*length(dom)-1), [], 1);
outerIdx = true(size(S.patches{1}.xy, 1), 1);
outerIdx(innerIdx) = false;
innerx = S.patches{1}.xy(innerIdx, 1);
innery = S.patches{1}.xy(innerIdx, 2);
%bc(innerIdx) = u_bulk(innerx, innery);
%bc(outerIdx) = bc(innerIdx);
%bc(outerIdx) = 0*bc(outerIdx);

ubc = zeros(size(gamx));
for k = 1:length(dom)
    ubc((k-1)*n+(1:n)) = util.modchebvals2legvals(bc(2*(k-1)*(n-2)+(1:n-2)+n-2));
end

u_strip = S \ bc;

fprintf('Strip solve took........%.6gs\n', toc);

%% (3) Solve the Neumann glue problem
%
%     lap(u_glue) = 0             in B \ Gamma'
%        [u_glue] = 0             on Gamma'
%    [du_glue/dn] = -du_strip/dn  on Gamma'

tic

% Compute  Neumann jump on Gamma'
normal1 = cell2mat(Gamma1.normal);
nx1 = normal1(:,1);
ny1 = normal1(:,2);
gam1x = cell2mat(Gamma1.x); gam1x = gam1x(:,1);
gam1y = cell2mat(Gamma1.x); gam1y = gam1y(:,2);
u_bulk_dn = feval(diff(u_bulk,1,2), gam1x, gam1y) .* nx1 + ...
            feval(diff(u_bulk,1,1), gam1x, gam1y) .* ny1;
% strip_dx = reshape(Dx * uu(:), m, n); strip_dx = strip_dx(m,:).';
% strip_dy = reshape(Dy * uu(:), m, n); strip_dy = strip_dy(m,:).';
% u_strip_dn = -strip_dx.*nx1 - strip_dy.*ny1;

u_strip_dn_all = S.patches{1}.D2N * [bc ; 1];
u_strip_dn = zeros(size(u_bulk_dn));
for k = 1:length(dom)
    u_strip_dn((k-1)*n+(1:n)) = util.modchebvals2legvals(u_strip_dn_all(2*(k-1)*(n-2)+(1:n-2)));
end

neu_jump = u_strip_dn - u_bulk_dn;
u_glue_i = @(x,y) kernels.laplace.slp(Gamma1, 'density', neu_jump, 'target', [x y], 'closeeval', true, 'side', 'i');
u_glue_e = @(x,y) kernels.laplace.slp(Gamma1, 'density', neu_jump, 'target', [x y], 'closeeval', true, 'side', 'e');

fprintf('Glue correction took....%.6gs\n', toc);

%% (4) Solve the boundary correction problem
%
%    lap(u_bc) = 0                                  in B \ Gamma
%         u_bc = g - u_bulk|_Gamma - u_glue|_Gamma  on Gamma

bc = g(gamx,gamy) - ubc - u_glue_e(gamx,gamy);
K = kernels.laplace.dlp(Gamma);
I = eye(n * Gamma.np);
sigma = (K - I/2) \ bc;
u_bc = @(x,y) kernels.laplace.dlp(Gamma, 'density', sigma, 'target', [x y], 'closeeval', true, 'side', 'i');

%% Add up the solutions
uu = nan(nbulk);
ub = u_bulk(xx,yy);
uu(inGamma1) = ub(inGamma1) + u_glue_i(xx_g1,yy_g1) + u_bc(xx_g1,yy_g1);

stripvals = zeros(size(xx_s));
for k = 1:length(dom)
    stripvals(globalIdx{k}) = util.bary2d(u_strip{k}, tt_s{k}, rr_s{k});
end
uu(inStrip) = stripvals + u_glue_e(xx_s, yy_s) + u_bc(xx_s, yy_s);

uu_sol = nan(nbulk);
uu_sol(inGamma) = sol(xx_g,yy_g);
uu_sol1 = nan(nbulk);
uu_sol1(inGamma1) = sol(xx_g1,yy_g1);

uu_solstrip = nan(nbulk);
uu_solstrip(inStrip) = sol(xx_s,yy_s);

ug = nan(size(xx));
ug(inGamma1) = uu(inGamma1);

%surf(xx, yy, ug - uu_sol1)
%surf(xx, yy, uu_sol1)
surf(uu-uu_sol)
surf(log10(abs(uu-uu_sol)))

shading interp

%%
plotopts = {'-k.', 'LineWidth', 2, 'MarkerSize', 14};

% subplot(121)
% hold on
% pcolor(xx, yy, uu)
% plot(Gamma,  plotopts{:})
% plot(Gamma1, plotopts{:}, 'Color', 'r')
% hold off
% shading interp
% colorbar, axis equal, axis(dom_global)
% shg
% 
% subplot(122)
% hold on
% pcolor(xx, yy, abs(uu-uu_sol))
% plot(Gamma,  plotopts{:})
% plot(Gamma1, plotopts{:}, 'Color', 'r')
% hold off
% shading interp
% colorbar, axis equal, axis(dom_global)
% shg

%norm(uu(inGamma)-uu_sol(inGamma), inf)
%return

vv_g1 = nan(nbulk);
ub = u_bulk(xx,yy);
vv_g1(inGamma1) = ub(inGamma1) + u_glue_i(xx_g1,yy_g1) + u_bc(xx_g1,yy_g1);
sol_g1 = nan(nbulk);
sol_g1(inGamma1) = sol(xx_g1,yy_g1);
err_g1 = norm(vv_g1(inGamma1)-sol_g1(inGamma1), inf);

t = linspace(dom_strip(1), dom_strip(2), n);
r = chebpts(m, dom_strip(3:4));
[tt,rr] = meshgrid(t,r);
xxs = xfun(tt,rr);
yys = yfun(tt,rr);
vv_s = u_strip(tt(:),rr(:)) + u_glue_e(xxs(:),yys(:)) + u_bc(xxs(:),yys(:));
err_s = norm(vv_s - sol(xxs(:),yys(:)), inf);

pcolor(xx, yy, log10(abs(vv_g1-sol_g1)))
%pcolor(xx, yy, vv_g1)
hold on
pcolor(xxs, yys, log10(abs(reshape(vv_s,m,n)-reshape(sol(xxs(:),yys(:)),m,n))))
%pcolor(xxs, yys, reshape(vv_s,m,n))
hold off
shading interp
hold on, plot(Gamma,'b-','Linewidth',3), plot(Gamma1, 'r-','Linewidth',3), shg
colorbar, axis(dom_global/scale), axis equal
shg

fprintf('\n');
fprintf('Error inside Gamma'' = %g\n', err_g1);
fprintf('Error inside strip  = %g\n', err_s);
fprintf('\n');
