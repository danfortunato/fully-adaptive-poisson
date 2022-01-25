N = 32;
[D,x] = cheb(N);
D2 = D^2;

f = exp(4*x(1:N+1));

%D2 = D2(2:N,2:N); % boundary conditions
b = find(x == -1 | x == 1);
D2(b,:) = 0;
D2(b,b) = speye(2);
f(b) = 0;

u = D2\f; % Poisson eq. solved here
%u = [0;u;0];

clf, subplot('position',[.1 .4 .8 .5])
plot(x,u,'.','markersize',16)
%xx = -1:.01:1;
xx = x;
%uu = polyval(polyfit(x,u,N),xx); % interpolate grid data
uu = u;
line(xx,uu,'linewidth',.8)
grid on
exact = ( exp(4*xx) - sinh(4)*xx - cosh(4) )/16;
title(['max err = ' num2str(norm(uu-exact,inf))],'fontsize',12)

%%

N = 30;
D = cheb(N);
x = chebpts(N+1);
y = chebpts(N+1);
%[xx,yy] = meshgrid(x(2:N),y(2:N));
[xx,yy] = meshgrid(x,y);

%exact = chebfun2(@(x,y) (1-y.^2).*(1-x.^2).*sin(x+y).*cos(y));
exact = randnfun2;
exact = chebfun2(@(x,y) (1-y.^2).*(1-x.^2).*exact(x,y));
f = lap(exact);
ff = f(xx,yy);

D2 = D^2;
%D2 = D2(2:N,2:N);
%I = speye(N-1);
I = speye(N+1);
L = kron(I,D2) + kron(D2,I);

b = find(xx == -1 | xx == 1 | yy == -1 | yy == 1);
%L(b,:) = [];
%L(:,b) = [];
L(b,b) = speye(size(b,1));
%ff(b) = [];
%ff(b) = exact(xx(b),yy(b));

tic, u = L\ff(:); toc

% Reshape long 1D results onto 2D grid:
uu = reshape(u,N+1,N+1);
%uu = zeros(N+1,N+1);
%uu(2:N,2:N) = reshape(u,N-1,N-1);

u = chebfun2(uu);
plot(u-exact)
[xx,yy] = meshgrid(linspace(-1,1,100));
surf(log10(abs(u(xx,yy) - exact(xx,yy))))
%zlim([-16 0])
title(['max err = ' num2str(norm(uu-exact,inf))],'fontsize',12)
shg
