function [u, du] = solve(L, rhs, opts)

defaults = [];
defaults.debug = true;
if ( nargin < 3 )
    opts = defaults;
else
    opts = setDefaults(opts, defaults);
end

dom = L.dom;
numPatches = length(dom);
[ny, nx] = size(dom(1).x);

[X, Y] = chebpts2(nx, ny);        % Chebyshev points and grid.
ii = ( abs(X) < 1 & abs(Y) < 1 ); % Interior indices.
ee = ( abs(X) == 1 );             % Boundary indices.
innerIdx = ( Y ==  1 );           % Inner indices.
outerIdx = ( Y == -1 );           % Outer indices.
numBdyPts = sum(ee(:)); 
numIntPts = sum(ii(:));

% Skeleton mappings
nxskel = nx-2;
nyskel = ny-2;
numSkelPts = 2*nyskel;
S2L = (skel2leaf(ny, nyskel));
L2S = (leaf2skel(nyskel, ny));
xskel = chebpts(nxskel, 1); xleaf = chebpts(nx, 2);
yskel = chebpts(nyskel, 1); yleaf = chebpts(ny, 2);
Bx = barymat(xskel, xleaf);
By = barymat(yskel, yleaf);

% Skeleton indices for each side
leftSkel  = 1:nyskel;
rightSkel = (nyskel+1):(2*nyskel);

% Define scalar RHSs:
if ( isnumeric(rhs) && isscalar(rhs) )
    rhs = repmat(rhs, numIntPts, 1);
end

%% Initialize
Timer.tic();

S   = zeros(nx*ny, numSkelPts, numPatches);
D2N = zeros(numSkelPts, numSkelPts, numPatches);
u_part  = zeros(nx*ny, numPatches);
du_part = zeros(numSkelPts, numPatches);
inner_d = zeros(nx, nx*ny, numPatches);
outer_d = zeros(nx, nx*ny, numPatches);

% Loop over each patch:
for k = 1:numPatches
    % Define the left and right edges for this patch:
    domk = dom(k);
    x = domk.x;
    y = domk.y;

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
    Sk = A(ii,ii) \ [-A(ii,ee) rhs_eval];

    % Append boundary points to solution operator:
    tmpS = zeros(nx*ny, numBdyPts+1);
    tmpS(ii,:) = Sk;
    tmpS(ee,:) = eye(numBdyPts, numBdyPts+1);
    Sk = [tmpS(:,1:end-1) * S2L, tmpS(:,end)];

    [nl, nr, nd, nu] = normals(x, y);
    nl = -By * nl;
    nr = -By * nr;

    % Construct the D2N map:
    dx = L2S * Dx(ee,:) * Sk;
    dy = L2S * Dy(ee,:) * Sk;
    D2Nk = zeros(numSkelPts, numSkelPts+1);
    D2Nk(leftSkel,:)  = nl(:,1).*dx(leftSkel,:)  + nl(:,2).*dy(leftSkel,:);
    D2Nk(rightSkel,:) = nr(:,1).*dx(rightSkel,:) + nr(:,2).*dy(rightSkel,:);

    S(:,:,k)   = Sk(:,1:end-1);   u_part(:,k)  = Sk(:,end);
    D2N(:,:,k) = D2Nk(:,1:end-1); du_part(:,k) = D2Nk(:,end);
    
    % Normal derivative for the top and bottom:
    inner_d(:,:,k) = nu(:,1).*Dx(innerIdx,:) + nu(:,2).*Dy(innerIdx,:);
    outer_d(:,:,k) = nd(:,1).*Dx(outerIdx,:) + nd(:,2).*Dy(outerIdx,:);
end

Timer.toc('Initialize');

%% Merge
Timer.tic();

i1 = leftSkel;  s1 = rightSkel;
i2 = rightSkel; s2 = leftSkel;
Z12 = zeros(numel(i1), numel(i2));

% Save this data
S_m = zeros(nyskel, numSkelPts, numPatches-1);
u_part_m = zeros(nyskel, numPatches);

% Overwrite this data
D2N_m = D2N(:,:,1);
du_part_m = du_part(:,1);

for k = 1:numPatches-1
    % Merge patches k-1 and k
    A = -( D2N_m(s1,s1) + D2N(s2,s2,k+1) );
    z = [ D2N_m(s1,i1) D2N(s2,i2,k+1) ];
    z_part = du_part_m(s1) + du_part(s2,k+1);
    temp = A \ [z z_part];
    Sk = temp(:,1:end-1);
    u_partk = temp(:,end);

    % Compute new D2N maps:
    M = [ D2N_m(i1,s1) ; D2N(i2,s2,k+1) ];
    D2N_m = [ D2N_m(i1,i1) Z12 ; Z12.' D2N(i2,i2,k+1) ] + M * Sk;
    du_part_m = [ du_part_m(i1) ; du_part(i2,k+1) ] + M * u_partk;

    S_m(:,:,k)    = Sk;
    u_part_m(:,k) = u_partk;
end

% Final merge is periodic
A = -( D2N_m(s1,s1) + D2N_m(s2,s2) );
z_part = du_part_m(s1) + du_part_m(s2);
u_part_m(:,end) = A \ z_part;

Timer.toc('Merge');

%% Solve
Timer.tic();

% Skeleton solve
u_skel = zeros(nyskel, numPatches+1);
u_skel(:,1)   = u_part_m(:,end);
u_skel(:,end) = u_part_m(:,end);
bc = u_skel(:,[1 end]);
bc = bc(:);
for k = numPatches:-1:2
    u_skel(:,k) = u_part_m(:,k-1) + S_m(:,:,k-1) * bc;
    bc(rightSkel) = u_skel(:,k);
end

% Leaf solve
u = zeros(ny, nx, numPatches);
du = zeros(nx, numPatches, 2);
for k = 1:numPatches
    bc = u_skel(:,k:k+1);
    sol = u_part(:,k) + S(:,:,k) * bc(:);
    u(:,:,k) = reshape(sol, ny, nx);
    du(:,k,1) = inner_d(:,:,k)*sol(:);
    du(:,k,2) = outer_d(:,:,k)*sol(:);
end

Timer.toc('Solve');

end

%% curvedDiffmat
function [Dx, Dy] = curvedDiffmat(x, y)
%CURVEDDIFFMAT   Construct differentiation matrices in curved coordinates.
%
%   [DX, DY] = CURVEDDIFFMAT(X, Y) returns the M*N x M*N Cartesian
%   differentiation matrices DX and DY associated with spectral collocation
%   on an M x N grid of tensor-product second-kind Chebyshev points (X, Y).

[ny, nx] = size(x);

persistent Dt Dr Dt2d Dr2d
if ( size(Dt,1) ~= nx || size(Dr,1) ~= ny )
    Dt = diffmat(nx);
    Dr = diffmat(ny);
    Im = eye(ny);
    In = eye(nx);
    Dt2d = kron(Dt, Im);
    Dr2d = kron(In, Dr);
end

% Im = eye(ny);
% In = eye(nx);
% Dr2d = kron(Dt, Im);
% Dt2d = kron(In, Dr);

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
Dx = dtdx(:).*Dt2d + drdx(:).*Dr2d;
Dy = dtdy(:).*Dt2d + drdy(:).*Dr2d;

end

%% normals
function [nl, nr, nd, nu] = normals(x, y)
% Outward-pointing normal vectors on a patch.

[ny, nx] = size(x);

persistent Dt Dr
if ( size(Dt,1) ~= nx || size(Dr,1) ~= ny )
    Dt = diffmat(nx);
    Dr = diffmat(ny);
end

dx = Dr*flipud(x(:,1)); % flipud?
dy = Dr*flipud(y(:,1));
nl = [dy, -dx] ./ sqrt(dx.^2 + dy.^2);

dx = Dr*flipud(x(:,nx));
dy = Dr*flipud(y(:,nx));
nr = -[dy, -dx] ./ sqrt(dx.^2 + dy.^2);

dx = Dt*x(1,:).';
dy = Dt*y(1,:).';
nd = [dy, -dx] ./ sqrt(dx.^2 + dy.^2);

dx = Dt*x(ny,:).';
dy = Dt*y(ny,:).';
nu = -[dy, -dx] ./ sqrt(dx.^2 + dy.^2);

end

%% skel2leaf
function P = skel2leaf(nleaf, nskel)
%SKEL2LEAF   Boundary interpolation matrix.
%   SKEL2LEAF(NLEAF, NSKEL) returns the 2*NLEAF x 2*NSKEL matrix that maps
%   2 pieces of length-NSKEL first-kind boundary values to 2*NLEAF second-
%   kind boundary values, including the corners. At each corner, the
%   average of the two interpolated values is used.

[xskel, ~, wskel] = chebpts(nskel, 1);
[xleaf, ~, wleaf] = chebpts(nleaf, 2);
B = barymat(xleaf, xskel, wskel);
B([1 end],:) = 0.5*B([1 end],:);
P = blkdiag(B, B);

end

%% leaf2skel
function P = leaf2skel(nskel, nleaf)
%LEAF2SKEL   Boundary interpolation matrix.
%   LEAF2SKEL(NSKEL, NLEAF) returns the 2*NSKEL x 2*NLEAF matrix that maps
%   2*NLEAF second-kind boundary values to 2 pieces of length-NSKEL first-
%   kind boundary values.

[xskel, ~, wskel] = chebpts(nskel, 1);
[xleaf, ~, wleaf] = chebpts(nleaf, 2);
B = barymat(xskel, xleaf, wleaf);
P = blkdiag(B, B);

end
