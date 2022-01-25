function [uu, Dx, Dy, xx, yy] = curvedPoisson(f, x, y, n, varargin)
%CURVEDPOISSON   Spectral Poisson solver in curvilinear coordinates.
%
%   CURVEDPOISSON(F, X, Y, N) solves Poisson's equation with righthand side
%   F and zero Dirichlet boundary conditions using the coordinate
%   transformations given by X = X(t,r) and Y(t,r), specified as a function
%   handle or chebfun2. The discretization is N x N tensor-product
%   Chebyshev. The coordinate transformations X and Y are assumed to be
%   defined on the (t,r)-domain [0 2*pi]x[0 1] if specified as function
%   handles.
%
%   CURVEDPOISSON(F, X, Y, [M N]) solves using an M x N discretization.
%
%   CURVEDPOISSON(F, X, Y, [M N], 'periodic') solves using a tensor-product
%   Fourier-Chebyshev discretization.
%
%   CURVEDPOISSON(F, X, Y, [M N], 'periodic', DBC, UBC) applies Dirichlet
%   boundary conditions on the bottom and top given by the values in DBC
%   and UBC.
%
%   CURVEDPOISSON(F, X, Y, [M N], LBC, RBC, DBC, UBC) applies Dirichlet
%   boundary conditions on the left, right, bottom, and top given by the
%   values in LBC, RBC, DBC, and UBC.

%% Parse arguments
dom = [0 2*pi 0 1];
periodic = false;

if ( isa(x, 'chebfun2') && isa(y, 'chebfun2') )
    if ( domainCheck(x, y) )
        dom = x.domain;
    else
        error('X and Y must have the same domain.');
    end
end

if ( isscalar(n) )
    % Call is CURVEDPOISSON(F, X, Y, N, ...)
    m = n;
else
    % Call is CURVEDPOISSON(F, X, Y, [M N], ...)
    m = n(1);
    n = n(2);
end

lbc = zeros(m,1);
rbc = zeros(m,1);
dbc = zeros(n,1);
ubc = zeros(n,1);

if ( nargin == 5 )
    % Call is CURVEDPOISSON(F, X, Y, [M N], 'periodic')
    periodic = strcmpi(varargin{1}, 'periodic');
elseif ( nargin == 6 )
    % Call is CURVEDPOISSON(F, X, Y, [M N], 'periodic', DBC)
    periodic = strcmpi(varargin{1}, 'periodic');
    dbc = varargin{2};
    ubc = dbc;
elseif ( nargin == 7 )
    % Call is CURVEDPOISSON(F, X, Y, [M N], 'periodic', DBC, UBC)
    periodic = strcmpi(varargin{1}, 'periodic');
    dbc = varargin{2};
    ubc = varargin{3};
elseif ( nargin == 8 )
    % Call is CURVEDPOISSON(F, X, Y, [M N], LBC, RBC, DBC, UBC)
    lbc = varargin{1};
    rbc = varargin{2};
    dbc = varargin{3};
    ubc = varargin{4};
end

% Check compatibility at corners
if ( ~periodic && ~(lbc(1) == ubc(1) && lbc(m) == dbc(1) && ...
                    rbc(1) == ubc(n) && rbc(m) == dbc(n)) )
   error('Boundary conditions are incompatible.');
end

%% Set up grids
if ( periodic )
    t = trigpts(n, dom(1:2));
else
    t = chebpts(n, dom(1:2));
end
r = chebpts(m, dom(3:4));
[tt,rr] = meshgrid(t,r);

%% Construct operators
if ( periodic )
    [Dx, Dy, xx, yy] = util.curvedDiffmat(x, y, [m n], 'periodic');
else
    [Dx, Dy, xx, yy] = util.curvedDiffmat(x, y, [m n]);
end
L = Dx^2 + Dy^2;
ff = f(xx,yy);

%% Impose boundary conditions

% Find boundary indices
if ( periodic )
    b = find(rr == r(1) | rr == r(m));
else
    b = find(rr == r(1) | rr == r(m) | tt == t(1) | tt == t(n));
end

% Remove boundary DOFs
II = kron(eye(n),eye(m));
BC = zeros(m,n);
BC(:,1) = lbc;
BC(:,n) = rbc;
BC(1,:) = ubc;
BC(m,:) = dbc;
L(b,:) = II(b,:);
ff(b) = BC(b);

%% Solve
u = L \ ff(:);

uu = reshape(u,m,n);
if ( nargout == 1 )
    if ( periodic )
        uu = chebfun2(uu, dom, 'trigx');
    else
        uu = chebfun2(uu, dom);
    end
end

end
