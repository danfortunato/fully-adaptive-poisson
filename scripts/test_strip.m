% Global coordinate system is given by:
%
%    x(t) = X(t) + r*f(t)*n(t)

N = 100;

S = Boundary.star(N, 'wobble', 0.15);
%S = Boundary.ellipse(N, 2, 1);
speed = [S.speed{1}; S.speed{1}(1,:)];
curvature = [S.curvature{1}; S.curvature{1}(1,:)];
x = [S.x{1}; S.x{1}(1,:)];
nx = [S.normal{1}; S.normal{1}(1,:)];

%f = @(k) sign(speed(k)) .* min(abs(1./speed(k)), 0.2); % f(t) = speed at t
f = @(k) 0.3*speed(k); % f(t) = speed at t
ff = f(1:N+1);

m = 10;
h_ratio = 0.5;
dt = S.w{1}(1)/S.speed{1}(1);
bh = dt*min(speed)*h_ratio;
r = m*bh;

%strip = x + 0.1.*min(speed,1./max(abs(curvature),1)).*nx;
rs = linspace(0.1,0.5,10);
rs = chebpts(10,[0,0.5])';

clf
plot(S,'-o')
hold on
for r = rs
    strip = x + r.*ff.*nx;
    plot(strip(:,1),strip(:,2),'-k')
end
plot([x(:,1), strip(:,1)].', [x(:,2), strip(:,2)].', '-k')
hold off
shg

%%

m = 61; % Number of Chebyshev nodes/modes
n = 61; % Number of Fourier nodes/modes
mn = m*n;
periodic = true;
coeffs   = false;
r = chebpts(m, [0 1]);
if ( periodic )
    %t = linspace(0, 2*pi, n).';
    t = trigpts(n, [0 2*pi]);
else
    t = chebpts(n, [0 2*pi]);
end
[tt,rr] = meshgrid(t,r);

if ( ~exist('s','var') )
    s = rng;
end
rng(s)

%v = randnfun2(2, [0 2*pi -1 2], 'trig');
%v = randnfun2(2, [0 2*pi 0 1], 'trig');
%sol = chebfun2(@(x,y) y.*(y-1).*sin(x), [0 2*pi 0 1]);
%sol = chebfun2(@(x,y) y.*(y-1).*x.*(x-2*pi).*v(x,y), [0 2*pi 0 1]);
%sol = chebfun2(@(x,y) y.*(y-1).*x.*(x-2*pi), [0 2*pi -1 2]);
%sol = chebfun2(@(x,y) y.*(y-1).*cos(y).*v(x,y), [0 2*pi -1 2]);
%sol = chebfun2(@(x,y) (y-sin(x)).*(y-(sin(x)+1)), [0 2*pi -1 2]);
%surf(xx,yy,sol(xx,yy))

% Coordinate change: (r,t) -> (x,y)
%v = 0.1*randnfun(2, [0 2*pi], 'trig') + 2;
v = chebfun(@(x) sin(x), [0 2*pi]) + 2;
%v = chebfun(@(t) sin(t), [0 2*pi]);
x = chebfun2(@(r,t) t, [0 1 0 2*pi]);
%y = chebfun2(@(r,t) r, [0 1 0 2*pi]);
%y = chebfun2(@(r,t) r + v(t), [0 1 0 2*pi]);
y = chebfun2(@(r,t) (1-r).*v(t), [0 1 0 2*pi]);
%rcoord = chebfun2(@(x,y) y - v(x), [0 2*pi -1 2]);
%tcoord = chebfun2(@(x,y) x, [0 2*pi -1 2]);
rcoord = chebfun2(@(x,y) 1-y./v(x), [0 2*pi 0 3]);
tcoord = chebfun2(@(x,y) x, [0 2*pi 0 3]);
bdy = @(t) 1 + 0.2*cos(5*t);
bdy = @(t) 1/2+0*t;
a = 0.25; R = 1;
bdy2 = @(t) a*sin(t) + sqrt(R^2 - a^2*cos(t).^2);
tfun = chebfun2(@(r,t) t, [0 1 0 2*pi]);
%rfun = chebfun2(@(r,t) bdy(t) + 0.5*r + 0.2*r.*cos(5*t), [0 1 0 2*pi]);
%rfun = chebfun2(@(r,t) bdy(t) + 0.5*(2+sin(t)).*r, [0 1 0 2*pi]);
rfun = chebfun2(@(r,t) r.*bdy(t) + (1-r).*bdy2(t), [0 1 0 2*pi]);
[x,y] = pol2cart(tfun,rfun);
xx = x(rr,tt);
yy = y(rr,tt);

%sol = chebfun2(@(x,y) (y-v(x)).*(y-(v(x)+1)).*x.*(x-2*pi), [0 2*pi -1 2]);
%sol = chebfun2(@(x,y) y.*(y-v(x)).*x.*(x-2*pi), [0 2*pi 0 max(v)]);

sol = chebfun2(@(x,y) ((R)^2-(x.^2+(y-a).^2)).*(1/4-(x.^2+y.^2)), [-1 1 -1 1]);
%g = chebfun2(@(r,t) r.*(1-r).*sin(t), [0 1 0 2*pi]);
% sol = chebfun2(@(x,y) g(rcoord(x,y),tcoord(x,y)), [0 2*pi -1 2]);
%sol = chebfun2(@(x,y) g(rcoord(x,y),tcoord(x,y)), [0 2*pi 0 3]);
f = lap(sol);

ff = f(xx,yy);
if ( periodic && coeffs )
    ff = trigtech.vals2coeffs(ff.').';
end

% Compute the Jacobian
dxdr = feval(diff(x,1,2), rr, tt);
dxdt = feval(diff(x,1,1), rr, tt);
dydr = feval(diff(y,1,2), rr, tt);
dydt = feval(diff(y,1,1), rr, tt);
J = dxdr.*dydt - dxdt.*dydr;
% Compute d[r,t]/d[x,y] as the inverse of the Jacobian
drdx =  dydt ./ J;
drdy = -dxdt ./ J;
dtdx = -dydr ./ J;
dtdy =  dxdr ./ J;

% A = [diag(dxdr(:)) diag(dxdt(:)) ;
%      diag(dydr(:)) diag(dydt(:)) ];
% Ji = inv(A);
% drdx = Ji(1:mn,1:mn); drdx = reshape(diag(drdx), m, n);
% drdy = Ji(1:mn,mn+1:end); drdy = reshape(diag(drdy), m, n);
% dtdx = Ji(mn+1:end,1:mn); dtdx = reshape(diag(dtdx), m, n);
% dtdy = Ji(mn+1:end,mn+1:end); dtdy = reshape(diag(dtdy), m, n);

% Chebyshev differentiation on [0 1]
dom = [0 1];
scl = 2/(dom(end)-dom(1));
Dr = -scl * cheb(m-1);

% Fourier differentiation
dom = [0 2*pi];
scl = 2/(dom(end)-dom(1));
if ( periodic && coeffs )
    Dt = trigspec.diffmat(n, 1, 1);
elseif ( periodic )
    Dt = scl * trigtech.diffmat(n, 1);
else
    Dt = -scl * cheb(n-1);
end

Im = speye(m); In = speye(n);
Dx = spdiags(drdx(:),0,mn,mn) * kron(In,Dr) + spdiags(dtdx(:),0,mn,mn) * kron(Dt,Im);
Dy = spdiags(drdy(:),0,mn,mn) * kron(In,Dr) + spdiags(dtdy(:),0,mn,mn) * kron(Dt,Im);
L = Dx^2 + Dy^2;
%gr = rfun(r,0).';
%L = kron(In, Dr^2 + Dr./gr) + kron(Dt^2./gr.^2, Im);

% Impose boundary conditions
if ( periodic )
    b = find(rr == 0 | rr == 1);
else
    b = find(rr == 0 | rr == 1 | tt == 0 | tt == t(end));
end
%L(b,:) = 0;
%ubc = exp(1i*(-n/2:n/2-1));
%dbc = exp(1i*(-n/2:n/2-1));
%ubc = feval(trigpoly(0:n-1), t); ubc = ubc(:).';
%dbc = ubc;
%L(b,:) = kron([ubc; dbc], Im);
%L(b,:) = [ubc ; dbc];
if ( periodic && coeffs )
    L(b,b) = diag( reshape([exp(-1i*(-n/2:n/2-1).*t.'); exp(-1i*(-n/2:n/2-1).*t.')], 2*n, 1) );
elseif ( periodic )
    %L(b,b) = speye(2*n);
    L(b,:) = [];
    L(:,b) = [];
else
    %L(b,b) = speye(2*m+2*n-4);
    L(b,:) = [];
    L(:,b) = [];
end
%ff(b) = 0;
ff(b) = [];

% Solve
u = L \ ff(:);

% Plot
%uu = reshape(u,m,n);
uu = zeros(m,n);
if ( periodic )
    uu(2:m-1,:) = reshape(u,m-2,n);
else
    uu(2:m-1,2:n-1) = reshape(u,m-2,n-2);
end
if ( periodic && coeffs )
    uu = real( trigtech.coeffs2vals(uu.').' );
end
vv = u;

[uu, ~] = util.curvedPoisson(f, x, y, [m n], 'periodic');

u = chebfun2(uu, [0 2*pi 0 1]);

norm(uu - sol(xx,yy), 'fro')

subplot(121)
%xx = [xx xx(:,1)+2*pi];
xx = [xx xx(:,1)];
yy = [yy yy(:,1)];
uu = [uu uu(:,1)];
%dd = f(xx,yy);
%gg = dd ./ feval(lap(u),tt,rr);
%surf(reshape(L*vv,m-2,n)-dd(2:m-1,:))
%surf(rr,tt,feval(lap(u),xx,yy))
%surf(xx,yy,feval(lap(u),tt,rr))
% urr = feval(diff(u,2,1),tt,rr);
% ur  = feval(diff(u,1,1),tt,rr);
% utt = feval(diff(u,2,2),tt,rr);
% surf(rr,tt, urr + ur./rr + utt./rr.^2);
mesh(xx,yy,uu)
%surf(uu)
%gg = f(xx,yy);
%gg(2:m-1,2:n-1) = reshape(L*vv,m-2,n-2);
%surf(xx,yy,gg-f(xx,yy))
colorbar
view(0,90)
shading interp
axis equal tight

subplot(122)
%surf(xx,yy,f(xx,yy))
%surf(xx,yy,log10(abs(uu-sol(xx,yy))))
%surf(tt,rr,uu-sol(xx,yy))
surf(xx,yy,sol(xx,yy))
%[xx,yy] = meshgrid(linspace(0,2*pi,100),linspace(0,1,100));
%surf(xx,yy,log10(abs(uu-sol(xx,yy))))
colorbar
view(0,90)
shading interp
axis equal tight
shg

%% Chebyshev Laplacian

m = 61; % Number of Chebyshev nodes/modes
n = 61; % Number of Fourier nodes/modes
mn = m*n;
r = chebpts(m, [0 1]);
t = chebpts(n, [0 2*pi]);
[tt,rr] = meshgrid(t,r);

if ( ~exist('s','var') )
    s = rng;
end
rng(s)

% Coordinate change: (r,t) -> (x,y)
v = chebfun(@(x) 1*sin(x), [0 2*pi]) + 2;
x = chebfun2(@(r,t) t, [0 1 0 2*pi]);
y = chebfun2(@(r,t) (1-r).*v(t), [0 1 0 2*pi]);
xx = x(rr,tt);
yy = y(rr,tt);

u = randnfun2(2, [0 2*pi 0 3]);
f = lap(u);
uu = u(xx,yy);
ff = f(xx,yy);

%drdx = y.*cos(x/2)./(2*sin(x/2).^2);

% Compute the Jacobian
dxdr = feval(diff(x,1,2), rr, tt);
dxdt = feval(diff(x,1,1), rr, tt);
dydr = feval(diff(y,1,2), rr, tt);
dydt = feval(diff(y,1,1), rr, tt);
J = dxdr.*dydt - dxdt.*dydr;
% Compute d[r,t]/d[x,y] as the inverse of the Jacobian
drdx =  dydt ./ J;
drdy = -dxdt ./ J;
dtdx = -dydr ./ J;
dtdy =  dxdr ./ J;

% Chebyshev differentiation on [0 1]
dom = [0 1];
scl = 2/(dom(end)-dom(1));
Dr = -scl * cheb(m-1);

% Fourier differentiation
dom = [0 2*pi];
scl = 2/(dom(end)-dom(1));
Dt = -scl * cheb(n-1);

Im = speye(m); In = speye(n);
Dx = spdiags(drdx(:),0,mn,mn) * kron(In,Dr) + spdiags(dtdx(:),0,mn,mn) * kron(Dt,Im);
Dy = spdiags(drdy(:),0,mn,mn) * kron(In,Dr) + spdiags(dtdy(:),0,mn,mn) * kron(Dt,Im);
L = Dx^2 + Dy^2;

% Impose boundary conditions
% b = find(rr == 0 | rr == 1 | tt == 0 | tt == t(end));
% L(b,:) = [];
% L(:,b) = [];
% ff(b) = [];

vv = L * uu(:); % Take Laplacian
vv = reshape(vv, m, n);

norm(vv - ff, 'fro')

% subplot(121)
% surf(xx,yy,vv)
% colorbar
% view(0,90)
% shading interp
% axis equal tight

%subplot(122)
surf(xx,yy,log10(abs(vv-ff)))
colorbar
view(0,90)
shading interp
axis equal tight

%% Periodic Laplacian

m = 31; % Number of Chebyshev nodes/modes
n = 31; % Number of Fourier nodes/modes
mn = m*n;
r = chebpts(m, [0 1]);
%t = trigpts(n, [0 2*pi]);
t = chebpts(n, [0 2*pi]);
[tt,rr] = meshgrid(t,r);

if ( ~exist('s','var') )
    s = rng;
end
rng(s)

% Coordinate change: (r,t) -> (x,y)
%v = chebfun(@(x) sin(x), [0 2*pi]) + 2;
v = chebfun(@(x) 0.02*sin(2*pi*x), [0 1]) + 2;
x = chebfun2(@(r,t) t/2/pi, [0 1 0 2*pi]);
y = chebfun2(@(r,t) (1-r).*v(t/2/pi), [0 1 0 2*pi]);
%rcoord = chebfun2(@(x,y) 1-y./v(x), [0 1 0 3]);
%tcoord = chebfun2(@(x,y) 2*pi*x, [0 1 0 3]);
xx = x(rr,tt);
yy = y(rr,tt);

%g = chebfun2(@(r,t) r.*(1-r).*sin(t), [0 1 0 2*pi]);
%u = chebfun2(@(x,y) g(rcoord(x,y),tcoord(x,y)), [0 2*pi 0 3]);
u = chebfun2(@(x,y) sin(2*pi*x)+cos(2*pi*y), [0 1 0 3], 'trig');
%f = lap(u);
f = diff(u,1,2);
uu = u(xx,yy);
ff = f(xx,yy);

% Compute the Jacobian
dxdr = feval(diff(x,1,2), rr, tt);
dxdt = feval(diff(x,1,1), rr, tt);
dydr = feval(diff(y,1,2), rr, tt);
dydt = feval(diff(y,1,1), rr, tt);
J = dxdr.*dydt - dxdt.*dydr;
% Compute d[r,t]/d[x,y] as the inverse of the Jacobian
drdx =  dydt ./ J;
drdy = -dxdt ./ J;
dtdx = -dydr ./ J;
dtdy =  dxdr ./ J;

% Chebyshev differentiation on [0 1]
dom = [0 1];
scl = 2/(dom(end)-dom(1));
Dr = -scl * cheb(m-1);

% Fourier differentiation
dom = [0 2*pi];
scl = 2/(dom(end)-dom(1));
Dt = scl * trigtech.diffmat(n, 1);
Dt = -scl * cheb(n-1);

Im = speye(m); In = speye(n);
Dx = spdiags(drdx(:),0,mn,mn) * kron(In,Dr) + spdiags(dtdx(:),0,mn,mn) * kron(Dt,Im);
Dy = spdiags(drdy(:),0,mn,mn) * kron(In,Dr) + spdiags(dtdy(:),0,mn,mn) * kron(Dt,Im);
L = Dx^2 + Dy^2;
L = Dx;

% Impose boundary conditions
% b = find(rr == 0 | rr == 1);
% L(b,:) = [];
% L(:,b) = [];
% ff(b) = [];

vv = L * uu(:); % Take Laplacian
vv = reshape(vv, m, n);

norm(vv - ff, 'fro')

subplot(121)
surf(xx,yy,vv)
colorbar
view(0,90)
shading interp
axis square tight

subplot(122)
surf(xx,yy,ff)
colorbar
view(0,90)
shading interp
axis square tight
shg