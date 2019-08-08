function [Dx, Dy, xx, yy] = curvedDiffmat(x, y, n, periodic)
%CURVEDDIFFMAT   Construct differentiation matrices in curved coordinates.
%
%   [DX, DY] = CURVEDDIFFMAT(X, Y, N) returns the N^2 x N^2 Cartesian
%   differentiation matrices DX and DY associated with the Chebyshev
%   spectral collocation method at tensor-product second-kind Chebyshev
%   points XX and YY, using the coordinate transformations given by
%   X = X(t,r) and Y(t,r), specified as a function handle or chebfun2. The
%   coordinate transformations X and Y are assumed to be defined on the
%   (t,r)-domain [0 2*pi]x[0 1] if specified as function handles.
%
%   [DX, DY, XX, YY] = CURVEDDIFFMAT(X, Y, [M N]) constructs M*N x M*N
%   matrices.
%
%   [DX, DY, XX, YY] = CURVEDDIFFMAT(X, Y, N, 'periodic') constructs
%   differentiation matrices associated with Fourier-Chebyshev collocation.

%% Parse arguments
dom = [0 2*pi 0 1];

if ( isa(x, 'function_handle') )
    x = chebfun2(x, dom);
end
if ( isa(y, 'function_handle') )
    y = chebfun2(y, dom);
end
if ( ~domainCheck(x, y) )
    error('X and Y must have the same domain.');
end
dom = x.domain;

if ( isscalar(n) )
    % Call is CURVEDDIFFMAT(X, Y, N, ...)
    m = n;
else
    % Call is CURVEDDIFFMAT(X, Y, [M N], ...)
    m = n(1);
    n = n(2);
end

if ( nargin == 3 )
    periodic = false;
else
    periodic = strcmpi(periodic, 'periodic');
end

%% Set up grids
if ( periodic )
    t = trigpts(n, dom(1:2));
else
    t = chebpts(n, dom(1:2));
end
r = chebpts(m, dom(3:4));
[tt,rr] = meshgrid(t,r);
xx = x(tt,rr);
yy = y(tt,rr);

%% Compute the Jacobian

%%% The chebfun2 way:

% [dxdt, dxdr] = gradient(x);
% [dydt, dydr] = gradient(y);
% J = jacobian(x,y);
% dtdx = feval(-dydr ./ J, tt, rr);
% dtdy = feval( dxdr ./ J, tt, rr);
% drdx = feval( dydt ./ J, tt, rr);
% drdy = feval(-dxdt ./ J, tt, rr);

%%% The fast way:

% Compute the Jacobian
dxdt = feval(diff(x,1,2), tt, rr);
dxdr = feval(diff(x,1,1), tt, rr);
dydt = feval(diff(y,1,2), tt, rr);
dydr = feval(diff(y,1,1), tt, rr);
J = -(dxdr.*dydt - dxdt.*dydr);
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
Dx = dtdx(:).*Dt2d + drdx(:).*Dr2d;
Dy = dtdy(:).*Dt2d + drdy(:).*Dr2d;

end
