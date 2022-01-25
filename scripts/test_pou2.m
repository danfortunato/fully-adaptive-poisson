%% An adaptive partition-of-unity strip solver
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

n = 200;     % Number of quadrature nodes on Gamma & Gamma'
m = 16;      % Number of radial Chebyshev nodes in strip
nbulk = 200;   % Number of uniform nodes for bulk FFT solve

%Gamma = Boundary.circle(n);
%Gamma = Boundary.squished_circle(n, 2.5, 0.2);
%Gamma = Boundary.squircle(n, 4, 1, 1);
Gamma = Boundary.star(n, 'wobble', 0.2, 'radius', 1.7);
gamxy = cell2mat(Gamma.x);
gamx  = gamxy(:,1);
gamy  = gamxy(:,2);

scale = 1.2;
dom_global = boundingbox(Gamma, scale);
%dom_global = [-pi pi -pi pi];

%sol = chebfun2(@(x,y) -cos(x).*exp(sin(x)).*sin(y).*(1-(x.^2+y.^2)), [-1.2 1.2 -1.2 1.2]);
%f = lap(sol);
sol = @(x,y) -cos(x).*exp(sin(x)).*sin(y);
f   = @(x,y) (2*cos(x)+3*cos(x).*sin(x)-cos(x).^3).*exp(sin(x)).*sin(y);
g   = @(x,y) sol(x,y);

% sol = randnfun2(1, dom_global);
% f = lap(sol);
% g = @(x,y) sol(x,y);

%% Define the fake interior boundary, Gamma'

fac = m;
h_ratio = 0.5;
hx = diff(dom_global(1:2))/nbulk;
hy = diff(dom_global(3:4))/nbulk;
h = min(hx,hy);
bh = h * min(Gamma.speed{1}) * h_ratio;
width = -fac * bh;
%width = -0.2;
blend_width = 2*m;

% Define the strip region
dom_strip = [0 2*pi 0 1];
curv = chebfun(Gamma.curvature, [0 2*pi], 'trig');
dff = chebfun(Gamma.dff, [0 2*pi], 'trig');
%[Gamma1, xyfun] = perturb(Gamma, @(t) width-0.01.*abs(Gamma.dff(t)), dom_strip);
[Gamma1, xyfun] = perturb(Gamma, @(t) width-0.1.*(max(abs(curv))-abs(curv(t))), dom_strip);
%[Gamma1, xyfun] = perturb(Gamma, @(t) width-0.01.*(max(abs(dff))-abs(Gamma.dff(t))), dom_strip);
xfun = real(xyfun);
yfun = imag(xyfun);
gam1xy = cell2mat(Gamma1.x);
gam1x  = gam1xy(:,1);
gam1y  = gam1xy(:,2);

%% Get grids
tic

[xx, yy] = meshgrid(trigpts(nbulk, dom_global(1:2)), ...
                    trigpts(nbulk, dom_global(3:4)));
%[xx, yy] = meshgrid(linspace(dom_global(1), dom_global(2), nbulk), ...
%                    linspace(dom_global(3), dom_global(4), nbulk));

inGamma  = isinterior(Gamma,  xx, yy);
inGamma1 = isinterior(Gamma1, xx, yy);
inStrip  = inGamma & ~inGamma1;

xx_g = xx(inGamma); xx_g1 = xx(inGamma1); xx_s = xx(inStrip);
yy_g = yy(inGamma); yy_g1 = yy(inGamma1); yy_s = yy(inStrip);
[tt_s, rr_s] = util.cart2curv(Gamma, xfun, yfun, xx_s, yy_s);

fprintf('\n');
fprintf('Grid computation took...%.6gs\n', toc);

%% (1) Solve the bulk problem
%
%    lap(u_bulk) = phi*f     in B
%         u_bulk = periodic  on B
tic

ff = zeros(nbulk);
ff(inGamma) = f(xx_g, yy_g);

[step, bump] = util.makeMask(1000, blend_width);
phi = zeros(nbulk);
phi(inStrip) = step(rr_s);
phi(inGamma1) = 1;
bumpf = @(x) (-1<=x&x<=1).*bump((x+1)/2);

% TODO: Make phi*f mean zero.
rhs = phi .* ff;

% unit integral bump function located at the corner of the rectangular domain
% this should not intersect the real domnain!
% easy to guarantee by increasing the size of the rectangular domain a bit
loc = 0.5;
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

ubc = u_bulk(gam1x, gam1y);
dbc = ubc;
[uu, Dx, Dy] = util.curvedPoisson(f, xfun, yfun, [m n], 'periodic', dbc, ubc);
us = uu;
u_strip = chebfun2(uu, xfun.domain, 'trigx');

fprintf('Strip solve took........%.6gs\n', toc);

%% (3) Solve the glue problem
%
%     lap(u_glue) = 0             in B \ Gamma'
%        [u_glue] = 0             on Gamma'
%    [du_glue/dn] = -du_strip/dn  on Gamma'
tic

% Compute Dirichlet and Neumann jumps on Gamma'
%dir_jump = 0*u_bulk(xy(:,1),xy(:,2));
normal1 = cell2mat(Gamma1.normal);
nx1 = normal1(:,1);
ny1 = normal1(:,2);
u_bulk_dn =  feval(diff(u_bulk,1,2), gam1x, gam1y) .* nx1 + ...
             feval(diff(u_bulk,1,1), gam1x, gam1y) .* ny1;
strip_dx = reshape(Dx * uu(:), m, n); strip_dx = strip_dx(m,:).';
strip_dy = reshape(Dy * uu(:), m, n); strip_dy = strip_dy(m,:).';
u_strip_dn = -strip_dx.*nx1 - strip_dy.*ny1;

neu_jump = u_strip_dn - u_bulk_dn;
u_glue_i = @(x,y) kernels.laplace.slp(Gamma1, 'density', neu_jump, 'target', [x y], 'closeeval', true, 'side', 'i');
u_glue_e = @(x,y) kernels.laplace.slp(Gamma1, 'density', neu_jump, 'target', [x y], 'closeeval', true, 'side', 'e');

fprintf('Glue correction took....%.6gs\n', toc);

%% (4) Solve the boundary correction problem
%
%    lap(u_bc) = 0                                  in B \ Gamma
%         u_bc = g - u_bulk|_Gamma - u_glue|_Gamma  on Gamma

bc = g(gamx,gamy) - (ubc + u_glue_e(gamx,gamy));
K = kernels.laplace.dlp(Gamma);
I = eye(Gamma.N);
sigma = (K - I/2) \ bc;
u_bc = @(x,y) kernels.laplace.dlp(Gamma, 'density', sigma, 'target', [x y], 'closeeval', true, 'side', 'i');

%% Add up the solutions
uu = nan(nbulk);
ub = u_bulk(xx,yy);
uu(inGamma1) = ub(inGamma1) + u_glue_i(xx_g1,yy_g1) + u_bc(xx_g1,yy_g1);
uu(inStrip)  = u_strip(tt_s,rr_s)  + u_glue_e(xx_s, yy_s)  + u_bc(xx_s, yy_s);

uu_sol = nan(nbulk);
uu_sol(inGamma) = sol(xx_g,yy_g);

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
