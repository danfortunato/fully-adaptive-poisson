n = 16;      % Number of quadrature nodes per panel
nbulk = 100; % Number of uniform nodes for bulk FFT solve

%Gamma = Boundary.circle(n, 'quadrature', 'panel');
%Gamma = refine(Gamma, 1);
%Gamma = Boundary(Gamma.z);
Gamma = Gamma1;

scale = 1.2;
dom_global = boundingbox(Gamma, scale);
sol = randnfun2(1, dom_global, 'trig');
f = lap(sol);
g = @(x,y) sol(x,y);

[xx, yy] = meshgrid(trigpts(nbulk, dom_global(1:2)), ...
                    trigpts(nbulk, dom_global(3:4)));
inGamma  = isinterior(Gamma, xx, yy);
xx_g = xx(inGamma);
yy_g = yy(inGamma);
gamxy = cell2mat(Gamma.x);
gamx  = gamxy(:,1);
gamy  = gamxy(:,2);

% Solve the bulk problem
ff = f(xx,yy);
u_bulk = util.periodicPoisson(ff, dom_global);

% Correct the boundary conditions
bc = g(gamx,gamy) - u_bulk(gamx,gamy);
K = kernels.laplace.dlp(Gamma);
I = eye(n * Gamma.np);
sigma = (K - I/2) \ bc;
u_bc = @(x,y) kernels.laplace.dlp(Gamma, 'density', sigma, 'target', [x y], 'closeeval', false, 'side', 'i');

% Add up the solutions
uu = nan(nbulk);
ub = u_bulk(xx,yy);
tic
uu(inGamma) = ub(inGamma) + u_bc(xx_g,yy_g);
toc

err = log10(abs(uu-sol(xx,yy)));
err(isinf(err)) = log10(eps);
surf(xx, yy, err)
view(2)
shading interp
colorbar
shg

%% Now try with blending
n = 16;      % Number of quadrature nodes per panel
nbulk = 200; % Number of uniform nodes for bulk FFT solve

Gamma = Boundary.circle(n, 'quadrature', 'panel');
Gamma = refine(Gamma, 1);

scale = 1.2;
dom_global = boundingbox(Gamma, scale);
% sol = randnfun2(2, dom_global);
% f = lap(sol);
% g = @(x,y) sol(x,y);
sol = @(x,y) -cos(x).*exp(sin(x)).*sin(y);
f   = @(x,y) (2*cos(x)+3*cos(x).*sin(x)-cos(x).^3).*exp(sin(x)).*sin(y);
g   = @(x,y) sol(x,y);

% Define the fake interior boundary, Gamma'
[x, y] = smoothStrip(Gamma, n, 0.5);
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
for k = 1:Gamma.np
    Gamma.x{k}(:,1) = legvals2chebvals(real(Gamma.x{k}(:,1)));
    Gamma.x{k}(:,2) = legvals2chebvals(real(Gamma.x{k}(:,2)));
    xx{k} = (1-t).*cx{k}.' + t.*Gamma.x{k}(:,1).';
    yy{k} = (1-t).*cy{k}.' + t.*Gamma.x{k}(:,2).';
end
dom = cell2struct([xx yy], {'x','y'}, 2);

% Get grids
tic

[xx, yy] = meshgrid(trigpts(nbulk, dom_global(1:2)), ...
                    trigpts(nbulk, dom_global(3:4)));
inGamma  = isinterior(Gamma,  xx, yy);
inGamma1 = isinterior(Gamma1, xx, yy);
inStrip  = inGamma & ~inGamma1;

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
    idx = ~found & converged & abs(tt) <= 1+100*eps & abs(rr) <= 1+100*eps;
    tt_s{k} = tt(idx);
    rr_s{k} = rr(idx);
    globalIdx{k} = idx;
    rr_s_global(idx) = rr_s{k};
    found = found | idx;
end

% Did we miss any points?
idx = globalIdx{1};
for k = 2:length(dom)
    idx = idx | globalIdx{k};
end
if ( ~all(idx) )
    warning('Some grid points were not assigned to an element.');
end

fprintf('Grid computation took...%.6gs\n', toc);

tic
% Solve the bulk problem
ff = zeros(nbulk);
ff(inGamma) = f(xx_g, yy_g);

blend_width = 2*n;
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

% Correct the boundary conditions
gam1xy = cell2mat(Gamma1.x);
gam1x  = gam1xy(:,1);
gam1y  = gam1xy(:,2);
bc = g(gam1x,gam1y) - u_bulk(gam1x,gam1y);
K = kernels.laplace.dlp(Gamma1);
I = eye(n * Gamma1.np);
sigma = (K - I/2) \ bc;
u_bc = @(x,y) kernels.laplace.dlp(Gamma1, 'density', sigma, 'target', [x y], 'closeeval', true, 'side', 'i');

% Add up the solutions
uu = nan(nbulk);
ub = u_bulk(xx,yy);
tic
uu(inGamma1) = ub(inGamma1) + u_bc(xx_g1,yy_g1);
toc

err = log10(abs(uu-sol(xx,yy)));
err(isinf(err)) = log10(eps);
surf(xx, yy, err)
view(2)
shading interp
colorbar
shg