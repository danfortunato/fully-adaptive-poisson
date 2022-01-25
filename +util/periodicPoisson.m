function [u, X_cfs] = periodicPoisson(f, varargin)
%PERIODICPOISSON   Periodic spectral Poisson solver.
%
%   PERIODICPOISSON(F) solves Poisson's equation with periodic boundary
%   conditions on [-1,1]x[-1,1], with the righthand side F specified as
%   either a matrix of function values at equispaced points, a function
%   handle, or a chebfun2.
%
%   PERIODICPOISSON(F, [A B C D]) specifies a domain [A,B]x[C,D] on which
%   to solve.
%
%   PERIODICPOISSON(F, N) solves using an N x N discretization.
%
%   PERIODICPOISSON(F, M, N) solves using an M x N discretization.

%% Parse arguments
n = 21;
m = n;
dom = [-1 1 -1 1];

if ( isa(f, 'chebfun2') )
    dom = f.domain;
end

if ( nargin == 2 )
    if ( isscalar(varargin{1}) )
        n = varargin{1};
        m = n;
    else
        dom = varargin{1};
    end
elseif ( nargin == 3 )
    m = varargin{1};
    n = m;
    if ( isscalar(varargin{2}) )
        n = varargin{2};
    else
        dom = varargin{2};
    end
elseif ( nargin == 6 )
    m = varargin{1};
    n = varargin{2};
    dom = varargin{3};
end

if ( isa(f, 'function_handle') || isa(f, 'chebfun2') )
    [xx,yy] = meshgrid(trigpts(n, dom(1:2)), trigpts(m, dom(3:4)));
    F = f(xx,yy);
else
    F = f;
    [m,n] = size(F);
end

%% Solve using FFTs
% k = [ 0:ceil(m/2-1) -floor(m/2):-1 ] .* (2*pi)/diff(dom(3:4)); k = k.';
% l = [ 0:ceil(n/2-1) -floor(n/2):-1 ] .* (2*pi)/diff(dom(1:2));
% ilap = -1./(k.^2+l.^2);
% ilap(1,1) = 0;
% X = ifft2( fft2( F ) .* ilap );

k = [ -floor(m/2):-1 0:ceil(m/2-1) ] .* (2*pi)/diff(dom(3:4)); k = k.';
l = [ -floor(n/2):-1 0:ceil(n/2-1) ] .* (2*pi)/diff(dom(1:2));
ilap = -1./(k.^2+l.^2);
ilap(floor(m/2)+1, floor(n/2)+1) = 0;
F_cfs = trigtech.vals2coeffs( trigtech.vals2coeffs( F ).' ).';
X_cfs = F_cfs .* ilap;
X = trigtech.coeffs2vals( trigtech.coeffs2vals( X_cfs ).' ).';

u = chebfun2(X, dom, 'trig');

end
