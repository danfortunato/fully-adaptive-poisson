function [uu, xx, yy] = curvedPoisson(f, x, y, varargin)
%CURVEDPOISSON   Spectral Poisson solver in curvilinear coordinates.
%
%   CURVEDPOISSON(F, X, Y) solves Poisson's equation with zero Dirichlet
%   boundary conditions using the coordinate transformations given by
%   X = X(t,r) and Y(t,r), specified as a function handle or chebfun2.
%   The default discretization is tensor-product Chebyshev.
%
%   CURVEDPOISSON(F, X, Y, 'periodic') solves using a tensor-product
%   Fourier-Chebyshev discretization.
%
%   CURVEDPOISSON(F, X, Y, N) solves using an N x N discretization.
%
%   CURVEDPOISSON(F, X, Y, M, N) solves using an M x N discretization.

%% Parse arguments
periodic = false;
n = 21;
m = n;
dom = [0 2*pi 0 1];

if ( isa(x, 'function_handle') )
    x = chebfun2(x, dom);
end
if ( isa(y, 'function_handle') )
    y = chebfun2(y, dom);
end

if ( nargin == 4 )
    if ( isscalar(varargin{1}) )
        n = varargin{1};
        m = n;
    else
        periodic = strcmpi(varargin{1}, 'periodic');
    end
elseif ( nargin == 5 )
    n = varargin{1};
    m = n;
    if ( isscalar(varargin{2}) )
        m = varargin{2};
    else
        periodic = strcmpi(varargin{2}, 'periodic');
    end
elseif ( nargin == 6 )
    m = varargin{1};
    n = varargin{2};
    periodic = strcmpi(varargin{3}, 'periodic');
end

%% Set up grids
if ( periodic )
    n = n-1;
end

if ( periodic )
    t = trigpts(n, dom(1:2));
else
    t = chebpts(n, dom(1:2));
end
r = chebpts(m, dom(3:4));
[tt,rr] = meshgrid(t,r);

xx = x(tt,rr);
yy = y(tt,rr);
ff = f(xx,yy);

%% Compute the Jacobian
dxdt = feval(diff(x,1,2), tt, rr);
dxdr = feval(diff(x,1,1), tt, rr);
dydt = feval(diff(y,1,2), tt, rr);
dydr = feval(diff(y,1,1), tt, rr);
J = dxdr.*dydt - dxdt.*dydr;

% Compute d[t,r]/d[x,y] as the inverse of the Jacobian
dtdx = -dydr ./ J;
dtdy =  dxdr ./ J;
drdx =  dydt ./ J;
drdy = -dxdt ./ J;

%% Construct operators

sclt = 2/diff(dom(1:2));
sclr = 2/diff(dom(3:4));

if ( periodic )
    Dt = sclt * diffmat(n, 'periodic');
else
    Dt = sclt * diffmat(n);
end
Dr = sclr * diffmat(m);

Im = eye(m);
In = eye(n);
Dt2d = kron(Dt,Im);
Dr2d = kron(In,Dr);
Dx = diag(dtdx(:)) * Dt2d + diag(drdx(:)) * Dr2d;
Dy = diag(dtdy(:)) * Dt2d + diag(drdy(:)) * Dr2d;
L = Dx^2 + Dy^2;

%% Impose boundary conditions

% Find boundary indices
if ( periodic )
    b = find(rr == r(1) | rr == r(m));
else
    b = find(rr == r(1) | rr == r(m) | tt == t(1) | tt == t(n));
end

% Remove boundary DOFs
L(b,:) = [];
L(:,b) = [];
ff(b) = [];

%% Solve
u = L \ ff(:);

% Add in the boundary conditions
uu = zeros(m,n);
if ( periodic )
    uu(2:m-1,:) = reshape(u,m-2,n);
    if ( nargout == 1 )
        uu = chebfun2(uu, dom, {'trig', []});
    else
        % Wrap the data
        uu = [uu uu(:,1)];
        xx = [xx xx(:,1)];
        yy = [yy yy(:,1)];
    end
else
    uu(2:m-1,2:n-1) = reshape(u,m-2,n-2);
    if ( nargout == 1 )
        uu = chebfun2(uu, dom);
    end
end

end
