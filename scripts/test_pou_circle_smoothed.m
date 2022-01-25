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
psem = n;

Gamma = Boundary.circle(n, 'quadrature', 'panel');
Gamma = refine(Gamma, 2);
%Gamma = Boundary.star(n, 'quadrature', 'panel');
%Gamma = refine(Gamma);
%Gamma = Boundary.multiscale(n, 'quadrature', 'panel');

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
tic
beta = 4;
bleed = 5;
width = 0.5;
%[x, y] = smoothStrip(Gamma, n, 0.4);
[x, y] = smoothStrip2(Gamma, n, beta, bleed, width);
xleg = chebvals2legvals(reshape(x, n, Gamma.np)); xleg = xleg(:);
yleg = chebvals2legvals(reshape(y, n, Gamma.np)); yleg = yleg(:);
z1 = mat2cell(xleg + 1i*yleg, repmat(n, Gamma.np, 1), 1);
Gamma1 = Boundary(z1);

% zz = cell2mat(Gamma.z);
% normals = cell2mat(Gamma.normal);
% th = 0;
% R = [cos(th) -sin(th);
%      sin(th)  cos(th)];
% normals = (R*normals.').';
% cn = normals(:,1) + 1i*normals(:,2);
% zz1 = zz - 0.1*cn;
% z1 = mat2cell(zz1, repmat(n, Gamma.np, 1), 1);
% Gamma1 = Boundary(z1);

gam1xy = cell2mat(Gamma1.x);
gam1x  = gam1xy(:,1);
gam1y  = gam1xy(:,2);
% x = gam1xy(:,1);
% y = gam1xy(:,2);

% Build the strip grid
cx = mat2cell(x, repmat(n,Gamma.np,1), 1);
cy = mat2cell(y, repmat(n,Gamma.np,1), 1);
t = chebpts(psem, [0 1]);
xsem = chebpts(psem);
xx = cell(Gamma.np,1);
yy = cell(Gamma.np,1);
gamx_cheb  = cell(Gamma.np, 1);
gamy_cheb  = cell(Gamma.np, 1);
gam1x_cheb = cell(Gamma.np, 1);
gam1y_cheb = cell(Gamma.np, 1);
for k = 1:Gamma.np
    gamx_cheb{k}  = legvals2chebvals(real(Gamma.x{k}(:,1)));
    gamy_cheb{k}  = legvals2chebvals(real(Gamma.x{k}(:,2)));
    gam1x_cheb{k} = legvals2chebvals(real(Gamma1.x{k}(:,1)));
    gam1y_cheb{k} = legvals2chebvals(real(Gamma1.x{k}(:,2)));
    % Upsample
    gamx_cheb{k}  = bary(xsem, gamx_cheb{k});
    gamy_cheb{k}  = bary(xsem, gamy_cheb{k});
    gam1x_cheb{k} = bary(xsem, gam1x_cheb{k});
    gam1y_cheb{k} = bary(xsem, gam1y_cheb{k});
    %xx{k} = (1-t).*gam1x_cheb{k}.' + t.*gamx_cheb{k}.';
    %yy{k} = (1-t).*gam1y_cheb{k}.' + t.*gamy_cheb{k}.';
    xx{k} = t.*gam1x_cheb{k}.' + (1-t).*gamx_cheb{k}.';
    yy{k} = t.*gam1y_cheb{k}.' + (1-t).*gamy_cheb{k}.';
end
dom = cell2struct([xx yy], {'x','y'}, 2);
toc
%% Get grids - kd-tree
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
fprintf('Grid split took...%.6gs\n', toc);

tic
p = [xx_s yy_s];
kd = KDTree(p);
tt_s = cell(length(dom), 1);
rr_s = cell(length(dom), 1);
rr_s_global = zeros(size(xx_s));
found = false(size(xx_s));
globalIdx = cell(length(dom), 1);
for k = 1:length(dom)
    xfun = chebfun2(dom(k).x);
    yfun = chebfun2(dom(k).y);
    % Filter out the points that are not in the bounding box of this element
    pad = 0.01;
    box = [min(dom(k).x(:))-pad max(dom(k).x(:))+pad;
           min(dom(k).y(:))-pad max(dom(k).y(:))+pad];
    testIdx = kd.range(box);
    [tt, rr, converged] = cart2curv_element(xfun, yfun, xx_s(testIdx), yy_s(testIdx), [0 0]);
    goodIdx = ~found(testIdx) & converged & abs(tt) <= 1.01 & abs(rr) <= 1.01;
    tt_s{k} = tt(goodIdx);
    rr_s{k} = rr(goodIdx);
    gidx = testIdx(goodIdx);
    globalIdx{k} = gidx;
    %globalIdx{k} = false(size(xx_s));
    %globalIdx{k}(gidx) = true;
    rr_s_global(gidx) = rr_s{k};
    found(gidx) = true;
    %found = found | globalIdx{k};
end
fprintf('Assigning elements took...%.6gs\n', toc);

% Did we miss any points?
%idx = globalIdx{1};
idx = false(size(xx_s(:)));
for k = 1:length(dom)
    %idx = idx | globalIdx{k};
    idx(globalIdx{k}) = true;
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

%% Get grids - old way
% tic
% 
% [xx, yy] = meshgrid(trigpts(nbulk, dom_global(1:2)), ...
%                     trigpts(nbulk, dom_global(3:4)));
% inGamma  = isinterior(Gamma,  xx, yy);
% inGamma1 = isinterior(Gamma1, xx, yy);
% inStrip  = inGamma & ~inGamma1;
% inStripFunc = util.inStripFactory(Gamma, Gamma1);
% %inStrip2 = reshape(inStripFunc([xx(:) yy(:)]), size(xx));
% 
% xx_g = xx(inGamma); xx_g1 = xx(inGamma1); xx_s = xx(inStrip);
% yy_g = yy(inGamma); yy_g1 = yy(inGamma1); yy_s = yy(inStrip);
% 
% % Convert global grid coordinates to local element coordinates.
% % This involves (1) determining which element each point is located in and
% % (2) inverting the element's coordinate maps using Newton's method.
% tt_s = cell(length(dom), 1);
% rr_s = cell(length(dom), 1);
% rr_s_global = zeros(size(xx_s));
% found = false(size(xx_s));
% globalIdx = cell(length(dom), 1);
% for k = 1:length(dom)
%     k
%     xfun = chebfun2(dom(k).x);
%     yfun = chebfun2(dom(k).y);
%     [tt, rr, converged] = cart2curv_element(xfun, yfun, xx_s, yy_s);
%     %idx = ~found & converged & abs(tt) <= 1+100*eps & abs(rr) <= 1+100*eps;
%     idx = ~found & converged & abs(tt) <= 1.01 & abs(rr) <= 1.01;
%     tt_s{k} = tt(idx);
%     rr_s{k} = rr(idx);
%     globalIdx{k} = idx;
%     rr_s_global(idx) = rr_s{k};
%     found = found | idx;
% end
% 
% fprintf('\n');
% fprintf('Grid computation took...%.6gs\n', toc);
% 
% % Did we miss any points?
% idx = globalIdx{1};
% for k = 2:length(dom)
%     idx = idx | globalIdx{k};
% end
% if ( ~all(idx) )
%     warning('Some grid points were not assigned to an element.');
%     figure
%     plot(Gamma), hold on, plot(Gamma1), hold on
%     scatter(xx_s(idx), yy_s(idx), 'bo'), hold on
%     scatter(xx_s(~idx), yy_s(~idx), 'ro'), hold off
% end
% 
% % Plot the points for each element
% plot(Gamma), hold on
% plot(Gamma1)
% for k = 1:length(dom)
%     idx = globalIdx{k};
%     scatter(xx_s(idx), yy_s(idx)), hold on
% end
% hold off
% shg

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
%phi(inStrip) = step((-rr_s_global+1)/2);
phi(inStrip) = step((rr_s_global+1)/2);
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
[u_bulk, X_bulk] = util.periodicPoisson( rhs, dom_global );
%u_bulk = chebfun2(sol, dom_global);
%X_bulk = u_bulk(xx,yy);

fprintf('FFT solve took..........%.6gs\n', toc);

%% (2) Solve the strip problem
%
%    lap(u_strip) = f  in S
%               u = 0  on Gamma & Gamma'
tic

S = StripSolver(dom, f);
build(S);

%innerIdx = reshape((1:n-2).' + (n-2)*(0:2:2*length(dom)-1), [], 1);
%outerIdx = true(size(S.patches{1}.xy, 1), 1);
%outerIdx(innerIdx) = false;
outerIdx = false(size(S.patches{1}.xy, 1), 1);
outerIdx((1:n-2).' + (n-2)*(0:2:2*length(dom)-1)) = true;
innerIdx = ~outerIdx;
innerx = S.patches{1}.xy(innerIdx, 1);
innery = S.patches{1}.xy(innerIdx, 2);
outerx = S.patches{1}.xy(outerIdx, 1);
outery = S.patches{1}.xy(outerIdx, 2);

%bc = zeros(size(S.patches{1}.xy, 1), 1);
opts = struct();
%opts.modeord = 1;
%opts.nthreads = 8;
X_bulk_cfs = trigtech.vals2coeffs( trigtech.vals2coeffs( X_bulk ).' ).';
xy = S.patches{1}.xy;
xy(:,1) = util.scale_interval(xy(:,1), dom_global(1:2), [-pi pi]);
xy(:,2) = util.scale_interval(xy(:,2), dom_global(3:4), [-pi pi]);
plan = finufft_plan(2, [nbulk nbulk], 1, 1, 1e-14, opts);
plan.setpts(xy(:,1), xy(:,2));
u_bulk_on_strip = real(finufft2d2(xy(:,1), xy(:,2), 1, 1e-14, X_bulk_cfs.'));
u_bulk_on_strip = real(plan.execute(X_bulk_cfs.'));
bc = u_bulk_on_strip;
%u_bulk_on_strip = finufft2d2(xy(:,1), xy(:,2), 1, 1e-14, fftshift(fftshift(fft2(X_bulk), 1),2));
%bc = u_bulk_on_strip;
bc = u_bulk(S.patches{1}.xy(:,1), S.patches{1}.xy(:,2));
%bc = g(S.patches{1}.xy(:,1), S.patches{1}.xy(:,2));
%bc(innerIdx) = u_bulk(innerx, innery);
%bc(outerIdx) = bc(innerIdx);
%bc(outerIdx) = 0*bc(outerIdx);

% ubc = zeros(size(gamx));
% for k = 1:length(dom)
%     ubc((k-1)*n+(1:n)) = util.modchebvals2legvals(bc(2*(k-1)*(n-2)+(1:n-2)+n-2));
% end
ubc = reshape(bc(outerIdx), n-2, length(dom));
ubc = util.modchebvals2legvals(ubc);
ubc = ubc(:);

u_strip = S \ bc;

fprintf('Strip solve took........%.6gs\n', toc);

%% (3) Solve the Neumann glue problem
%
%     lap(u_glue) = 0             in B \ Gamma'
%        [u_glue] = 0             on Gamma'
%    [du_glue/dn] = -du_strip/dn  on Gamma'

tic

X_bulk_cfs = trigtech.vals2coeffs( trigtech.vals2coeffs( X_bulk ).' ).';
X1x = real(trigtech.coeffs2vals( trigtech.coeffs2vals( (2*pi)/diff(dom_global(3:4))*trigspec.diffmat(nbulk, 1, 1) * X_bulk_cfs ).' ).');
u1 = chebfun2(X1x, dom_global, 'trig');
X2x = trigtech.vals2coeffs( trigtech.vals2coeffs( feval(diff(u_bulk,1,1), xx, yy) ).' ).';

xy = [gam1x gam1y];
xy(:,1) = util.scale_interval(xy(:,1), dom_global(1:2), [-pi pi]);
xy(:,2) = util.scale_interval(xy(:,2), dom_global(3:4), [-pi pi]);
D = trigspec.diffmat(nbulk, 1);
X_bulk_dx = (2*pi)/diff(dom_global(1:2)) * X_bulk_cfs * D;
X_bulk_dy = (2*pi)/diff(dom_global(3:4)) * D * X_bulk_cfs;
X_bulk_dx_gam1 = real(plan.execute(X_bulk_dx.'));
X_bulk_dy_gam1 = real(plan.execute(X_bulk_dy.'));
%feval(diff(u_bulk,1,2), gam1x, gam1y)

% Compute Neumann jump on Gamma'
normal1 = cell2mat(Gamma1.normal);
nx1 = normal1(:,1);
ny1 = normal1(:,2);
gam1x = cell2mat(Gamma1.x); gam1x = gam1x(:,1);
gam1y = cell2mat(Gamma1.x); gam1y = gam1y(:,2);

% Slow way: call chebfun/feval
u_bulk_dn = feval(diff(u_bulk,1,2), gam1x, gam1y) .* nx1 + ...
            feval(diff(u_bulk,1,1), gam1x, gam1y) .* ny1;
% Fast way: use NUFFT
% xy = [gam1x gam1y];
% xy(:,1) = util.scale_interval(xy(:,1), dom_global(1:2), [-pi pi]);
% xy(:,2) = util.scale_interval(xy(:,2), dom_global(3:4), [-pi pi]);
% plan.setpts(xy(:,1), xy(:,2));
% D = trigspec.diffmat(nbulk, 1);
% X_bulk_dx = (2*pi)/diff(dom_global(1:2)) * X_bulk_cfs * D;
% X_bulk_dy = (2*pi)/diff(dom_global(3:4)) * D * X_bulk_cfs;
% X_bulk_dx_gam1 = real(plan.execute(X_bulk_dx.'));
% X_bulk_dy_gam1 = real(plan.execute(X_bulk_dy.'));
% u_bulk_dn2 = X_bulk_dx_gam1 .* nx1 + X_bulk_dy_gam1 .* ny1;
      
% strip_dx = reshape(Dx * uu(:), m, n); strip_dx = strip_dx(m,:).';
% strip_dy = reshape(Dy * uu(:), m, n); strip_dy = strip_dy(m,:).';
% u_strip_dn = -strip_dx.*nx1 - strip_dy.*ny1;
u_strip_dn_all = S.patches{1}.D2N * [bc ; 1];
u_strip_dn = reshape(u_strip_dn_all(innerIdx), n-2, length(dom));
u_strip_dn = util.modchebvals2legvals(u_strip_dn);
u_strip_dn = reshape(u_strip_dn, [], 1);

% Downsample

% u_strip_dn = zeros(size(u_bulk_dn));
% for k = 1:length(dom)
%     u_strip_dn((k-1)*n+(1:n)) = util.modchebvals2legvals(u_strip_dn_all(2*(k-1)*(n-2)+(1:n-2)));
% end

%neu_jump = (u_strip_dn - u_bulk_dn);
neu_jump = -(u_strip_dn + u_bulk_dn);
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

uu = real(uu);

%surf(xx, yy, ug - uu_sol1)
%surf(xx, yy, uu_sol1)
%surf(uu)
surf(uu-uu_sol)
%surf(log10(abs(uu-uu_sol)))

shading interp
shg

%%
% for k = 1:length(dom)
%     domx = dom(k).x;
%     domy = dom(k).y;
%     ugluek = reshape(u_glue_e(domx(:),domy(:)), size(domx));
%     ubck = reshape(u_bc(domx(:),domy(:)), size(domx));
%     ustripk = u_strip{k} + ugluek + ubck;
%     err = abs(sol(domx,domy) - ustripk);
%     %err = log10(abs(sol(domx,domy) - ustripk));
%     %err(isinf(err)) = log10(eps);
%     surf(domx,domy,err)
%     hold on
% end
% shading interp
% shg