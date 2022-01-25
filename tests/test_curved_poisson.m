function pass = test_curved_poisson( show )

if ( nargin == 0 )
    show = false;
end

tol = 1e-9;
pass = [];

%% Periodic problem
m = 20;
n = 120;
%m = 12;
%n = 400;

% Define a geometry and exact solution
dom = [0 2*pi 0 1];
a1 = 0.7; b1 = 0.7;
a2 = 1; b2 = 1;
%a1 = 0.5; b1 = 1;
%a2 = 2;   b2 = 1.5;
bdy_in  = @(t) a1*b1./sqrt((b1*cos(t)).^2 + (a1*sin(t)).^2);
bdy_out = @(t) a2*b2./sqrt((b2*cos(t)).^2 + (a2*sin(t)).^2);
tfun = chebfun2(@(t,r) t, dom);
rfun = chebfun2(@(t,r) (1-r).*bdy_in(t) + r.*bdy_out(t), dom);
[x,y] = pol2cart(tfun,rfun);
rng(0)
%bc = randnfun2(1, [-1 1 -1 1] * max(max(rfun)));
%sol = chebfun2(@(x,y) -(1-((x/a1).^2+(y/b1).^2)).*(1-((x/a2).^2+(y/b2).^2))+bc(x,y), [-1 1 -1 1] * max(max(rfun)));
k = 2*pi/3;
k = 6;
sol = chebfun2(@(x,y) exp(sin(k*x)).*sin(k*y), [-1 1 -1 1]);
%sol = randnfun2(1);
f = lap(sol);

tt = trigpts(n, [0 2*pi]);
[x_in, y_in]   = pol2cart(tfun(tt,0),rfun(tt,0));
[x_out, y_out] = pol2cart(tfun(tt,1),rfun(tt,1));
dbc = sol(x_out,y_out);
ubc = sol(x_in,y_in);
u = util.curvedPoisson(f, x, y, [m n], 'periodic', dbc, ubc);

[uu, ~, ~, xx, yy] = util.curvedPoisson(f, x, y, [m n], 'periodic', dbc, ubc);
uu = [uu uu(:,1)];
xx = [xx xx(:,1)];
yy = [yy yy(:,1)];
%mesh(xx,yy,uu,'edgecolor','k'), view(0,90)
%axis equal off

[xxx,yyy] = meshgrid(linspace(0,2*pi,n),linspace(0,1,m));
xxx = [xxx xxx(:,1)];
yyy = [yyy yyy(:,1)];
%norm(u(xxx,yyy)-sol(x(xxx,yyy),y(xxx,yyy)), inf)

pass(end+1) = ( norm(u-sol(x,y), inf) < tol );

if ( show )
    figure(numel(pass))
    subplot(121)
    surf(x,y,u)
    colorbar, view(0,90), shading interp, axis equal tight
    subplot(122)
    %surf(xx,yy,log10(abs(uu-sol(xx,yy))))
    surf(x(xxx,yyy),y(xxx,yyy),log10(abs(u(xxx,yyy)-sol(x(xxx,yyy),y(xxx,yyy)))))
    colorbar, view(0,90), shading interp, axis equal tight
    set(gca, 'FontSize', 18)
    shg
end

%% Patch problem
m = 30;
n = 100;

% Define a geometry and exact solution
dom = [0 2*pi 0 1];
bdy_up   = chebfun(@(x) sin(x),     dom(1:2)) + 2;
bdy_down = chebfun(@(x) 0.2*sin(x), dom(1:2));
%bdy_up = chebfun(@(x) 0*x+1, dom(1:2));
%bdy_down = chebfun(@(x) 0*x, dom(1:2));
x = chebfun2(@(t,r) t, dom);
y = chebfun2(@(t,r) (1-r).*bdy_down(t) + r.*bdy_up(t), dom);
%sol = randnfun2(1, [dom(1:2) min(bdy_down) max(bdy_up)]);
sol = chebfun2(@(x,y) sin(2*pi*x).*sin(2*pi*y), [dom(1:2) min(bdy_down) max(bdy_up)]);
f = lap(sol);

tt = chebpts(n, dom(1:2));
rr = chebpts(m, dom(3:4));
x_l = x(dom(1),rr); y_l = y(dom(1),rr);
x_r = x(dom(2),rr); y_r = y(dom(2),rr);
x_u = x(tt,dom(3)); y_u = y(tt,dom(3));
x_d = x(tt,dom(4)); y_d = y(tt,dom(4));

lbc = sol(x_l,y_l);
rbc = sol(x_r,y_r);
ubc = sol(x_u,y_u);
dbc = sol(x_d,y_d);
u = util.curvedPoisson(f, x, y, [m n], lbc, rbc, dbc, ubc);

[uu, ~, ~, xx, yy] = util.curvedPoisson(f, x, y, [m n], lbc, rbc, dbc, ubc);

pass(end+1) = ( norm(u-sol(x,y), inf) < tol );

if ( show )
    figure(numel(pass))
    subplot(121)
    surf(x,y,u)
    colorbar, view(0,90), shading interp, axis equal tight
    subplot(122)
    surf(xx,yy,log10(abs(uu-sol(xx,yy))))
    colorbar, view(0,90), shading interp, axis equal tight
    shg
end

end
