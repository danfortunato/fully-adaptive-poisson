n = 60;
pref = chebfunpref;
pref.tech = @chebtech1;
sol = chebfun2(@(x,y) sin(10*pi*x).*sin(10*pi*y));
f = lap(sol);

% Rectangular double differentiation matrix: (n-2) x n
D2 = diffmat([n-2 n], 2, 'chebkind1');

% Rectangular double conversion matrix: (n-2) x n
[x0, ~, v0] = chebpts(n, [-1 1], 1);
[y0, ~, w0] = chebpts(n, [-1 1], 2);
x2 = chebpts(n-2, [-1 1], 1);
S02 = barymat(x2, x0, v0);
S00 = barymat(y0, x0, v0);

% Rectangular Laplacian matrix: (n-2)^2 x n^2
L = kron(D2,S02) + kron(S02,D2);

% Set the boundary conditions: (4n-4) x n^2
bb = zeros(n); bb(:,[1,n]) = 1; bb([1,n],:) = 1; bb = find(bb);
B = zeros(4*n-4,n^2);
% ii = sub2ind(size(B), 1:4*n-4, bb.');
% B(ii) = 1;

[xx, yy] = meshgrid(chebpts(n, [-1 1], 2));
for k = 1:4*n-4
    idx = bb(k);
    B(k,:) = kron(barymat(xx(idx), x0, v0), barymat(yy(idx), x0, v0));
end

% Solve: n^2 x n^2
[xx, yy] = meshgrid(x2);
ff = f(xx,yy); ff = ff(:);
bc = zeros(4*n-4,1);
uu = [L; B] \ [ff; bc];

uu = reshape(uu, n, n);
vv = S00*uu*S00.';
u = chebfun2(vv);
plot(u-sol)
norm(u-sol)

%%
n = 50;
sol = chebfun2(@(x,y) sin(5*pi*x).*sin(5*pi*y));
f = lap(sol);
dom = specdd.rectangle([-1 1 -1 1]);
pdo = {1,0,0};
bc = 0;
rhs = @(x,y) f(x,y);
S = specdd(dom, pdo, rhs, n);
u = S \ bc;

xx = u.x{1}; yy = u.y{1}; uu = u.u{1};
surf(xx,yy,uu-sol(xx,yy))
shg
