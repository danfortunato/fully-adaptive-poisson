function L = initialize(dom, rhs)
%INITIALIZE   Initialize an array of LEAF objects.
%   L = STRIPSOLVER.LEAF.INITIALIZE(DOM) returns a cell array L of LEAF
%   objects which contain the solution and D2N operators for Poisson's
%   equation on the domain DOM with zero righthand side.
%
%   L = STRIPSOLVER.LEAF.INITIALIZE(DOM, RHS) is as above, but with the
%   righthand side RHS, which may be a scalar or a function handle.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(isstruct(dom), 'Invalid domain.');

if ( nargin < 2 )
    % Default to homogeneous problem:
    rhs = 0;
end

numPatches = length(dom);
n = size(dom(1).x, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% DEFINE REFERENCE GRID %%%%%%%%%%%%%%%%%%%%%%%%%%

[X, Y] = chebpts2(n);             % Chebyshev points and grid.
ii = abs(X) < 1 & abs(Y) < 1;     % Interior indices.
ii([1,end],[1,end]) = true;       % Treat corners as interior.
ee = ~ii;                         % Boundary indices.
leftIdx  = 1:n-2;
rightIdx = n-1:2*(n-2);
downIdx  = 2*(n-2)+1:3*(n-2);
upIdx    = 3*(n-2)+1:4*(n-2);
numBdyPts = sum(ee(:)); 
numIntPts = sum(ii(:));
ibc = 3*(n-2)+1;
ss = [1:n-2, ibc:4*(n-2), n-1:2:ibc-1, n:2:ibc-1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%% DEFINE OPERATORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Interpolation operator for corner values:
Xii = X(1,2:(n-1)).';
B = [Xii-1, -1-Xii].'; B(:,1:2:end) = -B(:,1:2:end);
if ( mod(n-1, 2) )
    B(2,:) = -B(2,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANT RHS? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define scalar RHSs:
if ( isnumeric(rhs) && isscalar(rhs) )
    rhs = repmat(rhs, numIntPts, 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% SOLVE LOCAL PROBLEMS %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Each patch needs a different solution operator.

% Initialize
L = cell(numPatches, 1);

% Loop over each patch:
for k = 1:numPatches

    % Define the left and right edges for this patch:
    domk = dom(k);
    x = domk.x;
    y = domk.y;
    edges = [ x(n,1) y(n,1) x(1,1) y(1,1) n ;  % "Left" side
              x(n,n) y(n,n) x(1,n) y(1,n) n ;  % "Right" side
              x(n,1) y(n,1) x(n,n) y(n,n) n ;  % "Down" side
              x(1,1) y(1,1) x(1,n) y(1,n) n ]; % "Up" side

    % Evaluate non-constant RHSs if required:
    if ( isa(rhs, 'function_handle') || isa(rhs, 'chebfun2') )
        rhs_eval = feval(rhs, x(ii), y(ii));
    elseif ( iscell(rhs) )
        rhs_eval = rhs{k};
    else
        rhs_eval = rhs;
    end

    [Dx, Dy] = curvedDiffmat(x, y);
    A = Dx^2 + Dy^2;

    % Construct solution operator:
    S = A(ii,ii) \ [-A(ii,ee), rhs_eval];
    Ainv = @(u) A(ii,ii) \ u;

    % Replace solution operator for corners with interp conditions:
    S([1:2,end-1:end],:) = 0;
    S(1:2,1:n-2) = B;
    S([end-1,end],end-n+2:end-1) = B;
    
    % Append boundary points to solution operator:
    tmpS = zeros(n^2, size(S, 2));
    tmpS(ii,:) = S;
    tmpS(ee,:) = eye(numBdyPts, numBdyPts + 1);
    S = tmpS;
    S = S(:,[ss end]);

    % Construct the D2N map:
    dx = Dx(ee,:) * S; dx = dx(ss,:);
    dy = Dy(ee,:) * S; dy = dy(ss,:);
    [nl, nr, nd, nu] = normals(x, y);
    D2N = zeros(numBdyPts, numBdyPts + 1);
    D2N(leftIdx,:)  = nl(:,1).*dx(leftIdx,:)  + nl(:,2).*dy(leftIdx,:);
    D2N(rightIdx,:) = nr(:,1).*dx(rightIdx,:) + nr(:,2).*dy(rightIdx,:);
    D2N(downIdx,:)  = nd(:,1).*dx(downIdx,:)  + nd(:,2).*dy(downIdx,:);
    D2N(upIdx,:)    = nu(:,1).*dx(upIdx,:)    + nu(:,2).*dy(upIdx,:);

    % Assemble the patch:
    xee = x(ee); xee = xee(ss);
    yee = y(ee); yee = yee(ss);
    xy = [xee yee];
    L{k} = StripSolver.Leaf(domk, S, D2N, edges, xy, Ainv);

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE OPERATORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Dx, Dy] = curvedDiffmat(x, y)
%CURVEDDIFFMAT   Construct differentiation matrices in curved coordinates.
%
%   [DX, DY] = CURVEDDIFFMAT(X, Y) returns the M*N x M*N Cartesian
%   differentiation matrices DX and DY associated with spectral collocation
%   on an M x N grid of tensor-product second-kind Chebyshev points (X, Y).

[m, n] = size(x);

persistent Dt Dr
if ( size(Dt,1) ~= n )
    Dt = diffmat(n);
end
if ( size(Dr,1) ~= m )
    Dr = diffmat(m);
end

% Compute the Jacobian:
dxdt = x * Dt.'; dxdr = Dr * x;
dydt = y * Dt.'; dydr = Dr * y;
J = -(dxdr.*dydt - dxdt.*dydr);

% Compute d[t,r]/d[x,y] as the inverse of the Jacobian:
dtdx = -dydr ./ J;
dtdy =  dxdr ./ J;
drdx =  dydt ./ J;
drdy = -dxdt ./ J;

% Construct operators:
Im = eye(m);
In = eye(n);
Dt2d = kron(Dt, Im);
Dr2d = kron(In, Dr);
Dx = dtdx(:).*Dt2d + drdx(:).*Dr2d;
Dy = dtdy(:).*Dt2d + drdy(:).*Dr2d;

end

function [nl, nr, nd, nu] = normals(x, y)
%NORMALS   Outward pointing normal vectors to the edges of a mapping.

persistent D

n = size(x, 1);
if ( size(D,1) ~= n )
    D = diffmat(n);
end

dx = D*flipud(x(:,1)); % flipud?
dy = D*flipud(y(:,1));
nl = -[dy, -dx] ./ sqrt(dx.^2 + dy.^2);
nl([1,n],:) = [];

dx = D*flipud(x(:,n));
dy = D*flipud(y(:,n));
nr = [dy, -dx] ./ sqrt(dx.^2 + dy.^2);
nr([1,n],:) = [];

dx = D*x(1,:).';
dy = D*y(1,:).';
nd = -[dy, -dx] ./ sqrt(dx.^2 + dy.^2);
nd([1,n],:) = [];

dx = D*x(n,:).';
dy = D*y(n,:).';
nu = [dy, -dx] ./ sqrt(dx.^2 + dy.^2);
nu([1,n],:) = [];

% v = [ x([-1 1 1 -1],[-1 -1 1 1]) ;
%       y([-1 1 1 -1],[-1 -1 1 1]) ].';
% v = [ x(n,1) x(n,n) x(1,n) x(1,1) ;
%       y(n,1) y(n,n) y(1,n) y(1,1) ].';
% v = [ v ; v(1,:) ];
% dx = diff(v(:,1));
% dy = diff(v(:,2));
% vn = [dy.' ; -dx.' ];
% vn = vn(:,[4 2 1 3]); % Reorder to left, right, down, up
% vn = vn ./ sqrt(sum(vn.^2));  % Normalize

end

function D2N = transformD2N2(x, y, V, W, ee)

% dsdx = T.dinvT1{1}( T.T1(XX,YY), T.T2(XX,YY) );
% dsdy = T.dinvT1{2}( T.T1(XX,YY), T.T2(XX,YY) );
% dtdx = T.dinvT2{1}( T.T1(XX,YY), T.T2(XX,YY) );
% dtdy = T.dinvT2{2}( T.T1(XX,YY), T.T2(XX,YY) );

[m, n] = size(x);

% Compute the Jacobian:
Dt = diffmat(n);
Dr = diffmat(m);
dxdt = x * Dt.'; dxdr = Dr * x;
dydt = y * Dt.'; dydr = Dr * y;
J = -(dxdr.*dydt - dxdt.*dydr);

% Compute d[t,r]/d[x,y] as the inverse of the Jacobian:
dtdx = -dydr ./ J; dtdx = dtdx(ee);
dtdy =  dxdr ./ J; dtdy = dtdy(ee);
drdx =  dydt ./ J; drdx = drdx(ee);
drdy = -dxdt ./ J; drdy = drdy(ee);

dtdx = diag(dtdx);
dtdy = diag(dtdy);
drdx = diag(drdx);
drdy = diag(drdy);

dx = dtdx*V + drdx*W; % d/dx
dy = dtdy*V + drdy*W; % d/dy

% Now make normal derivatives:
ibc = (3*(n-2)+1);
leftIdx  = 1:n-2;
rightIdx = n-1:2*(n-2);       % ibc:4*(n-2)
downIdx  = 2*(n-2)+1:3*(n-2); % n-1:2:ibc-1
upIdx    = 3*(n-2)+1:4*(n-2); % n:2:ibc-1
V = dx;
[nl, nr, nd, nu] = normals(x, y);
V(leftIdx,:)  = diag(nl(:,1))*dx(leftIdx,:)  + diag(nl(:,2))*dy(leftIdx,:);
V(rightIdx,:) = diag(nr(:,1))*dx(rightIdx,:) + diag(nr(:,2))*dy(rightIdx,:);
V(downIdx,:)  = diag(nd(:,1))*dx(downIdx,:)  + diag(nd(:,2))*dy(downIdx,:);
V(upIdx,:)    = diag(nu(:,1))*dx(upIdx,:)    + diag(nu(:,2))*dy(upIdx,:);

D2N = V;

end
