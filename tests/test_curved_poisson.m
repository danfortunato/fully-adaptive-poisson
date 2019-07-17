%% Periodic problem
m = 21;
n = 21;

% Define a geometry and exact solution
a = 0.25;
R = 1;
bdy = @(t) 1/2+0*t;
bdy2 = @(t) a*sin(t) + sqrt(R^2 - a^2*cos(t).^2);
tfun = chebfun2(@(t,r) t, [0 2*pi 0 1]);
rfun = chebfun2(@(t,r) r.*bdy(t) + (1-r).*bdy2(t), [0 2*pi 0 1]);
[x,y] = pol2cart(tfun,rfun);
sol = chebfun2(@(x,y) (R^2-(x.^2+(y-a).^2)).*(1/4-(x.^2+y.^2)));
f = lap(sol);

[uu, xx, yy] = curvedPoisson(f, x, y, m, n, 'periodic');
norm(uu - sol(xx,yy), 'fro')

subplot(121)
surf(xx,yy,uu)
colorbar, view(0,90), shading interp, axis equal tight
subplot(122)
surf(xx,yy,sol(xx,yy))
colorbar, view(0,90), shading interp, axis equal tight
shg

%% Patch problem
m = 25;
n = 11;

% Define a geometry and exact solution
v = chebfun(@(x) sin(x), [0 2*pi]) + 2;
x = chebfun2(@(t,r) t, [0 2*pi 0 1]);
y = chebfun2(@(t,r) (1-r).*v(t), [0 2*pi 0 1]);
sol = chebfun2(@(x,y) y.*(y-v(x)).*x.*(x-2*pi), [0 2*pi 0 max(v)]);
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
