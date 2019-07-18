%% Periodic problem
m = 10;
n = 100;

% Define a geometry and exact solution
dom = [0 2*pi 0 1];
a1 = 0.5; b1 = 1;
a2 = 2;   b2 = 1.5;
bdy_in  = @(t) a1*b1./sqrt((b1*cos(t)).^2 + (a1*sin(t)).^2);
bdy_out = @(t) a2*b2./sqrt((b2*cos(t)).^2 + (a2*sin(t)).^2);
tfun = chebfun2(@(t,r) t, dom);
rfun = chebfun2(@(t,r) (1-r).*bdy_in(t) + r.*bdy_out(t), dom);
[x,y] = pol2cart(tfun,rfun);
sol = chebfun2(@(x,y) (1-((x/a1).^2+(y/b1).^2)).*(1-((x/a2).^2+(y/b2).^2)), [-1 1 -1 1] * max(max(rfun)));
f = lap(sol);

u = curvedPoisson(f, x, y, m, n, 'periodic');
norm(u-sol(x,y))

subplot(121)
surf(x,y,u)
colorbar, view(0,90), shading interp, axis equal tight
subplot(122)
surf(x,y,sol(x,y)-u)
colorbar, view(0,90), shading interp, axis equal tight
shg

%% Patch problem
m = 25;
n = 11;

% Define a geometry and exact solution
dom = [0 2*pi 0 1];
bdy_up   = chebfun(@(x) sin(x),     dom(1:2)) + 2;
bdy_down = chebfun(@(x) 0.2*sin(x), dom(1:2));
x = chebfun2(@(t,r) t, dom);
y = chebfun2(@(t,r) (1-r).*bdy_down(t) + r.*bdy_up(t), dom);
sol = chebfun2(@(x,y) (y-bdy_down(x)).*(y-bdy_up(x)).*x.*(x-2*pi), [dom(1:2) min(bdy_down) max(bdy_up)]);
f = lap(sol);

u = curvedPoisson(f, x, y, m, n);
norm(u - sol(x,y))

subplot(121)
surf(x,y,u)
colorbar, view(0,90), shading interp, axis equal tight
subplot(122)
surf(x,y,u-sol(x,y))
colorbar, view(0,90), shading interp, axis equal tight
shg
