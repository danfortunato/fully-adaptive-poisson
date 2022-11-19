function L = parinitialize(dom, rhs)
%INITIALIZE   Initialize an array of LEAF objects.
%   L = STRIPSOLVER.LEAF.PATINITIALIZE(DOM) returns a cell array L of LEAF
%   objects which contain the solution and D2N operators for Poisson's
%   equation on the domain DOM with zero righthand side.
%
%   L = STRIPSOLVER.LEAF.PARINITIALIZE(DOM, RHS) is as above, but with the
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

x = zeros(n, n, numPatches);
y = zeros(n, n, numPatches);
rhs_eval = zeros(numIntPts, numPatches);

% Loop over each patch:
for k = 1:numPatches
    x(:,:,k) = dom(k).x;
    y(:,:,k) = dom(k).y;
    % Evaluate non-constant RHSs if required:
    if ( isa(rhs, 'function_handle') || isa(rhs, 'chebfun2') )
        rhs_eval(:,k) = feval(rhs, x(ii), y(ii));
    elseif ( iscell(rhs) )
        rhs_eval(:,k) = rhs{k};
    else
        rhs_eval(:,k) = rhs;
    end
end

nthreads = 1;
[S, D2N, Aii] = constructOperators(x, y, rhs_eval, B, nthreads);

% Assemble the patches:
for k = 1:numPatches
    % Define the edges for this patch:
    edges = [ x(n,1,k) y(n,1,k) x(1,1,k) y(1,1,k) n ;  % "Left" side
              x(n,n,k) y(n,n,k) x(1,n,k) y(1,n,k) n ;  % "Right" side
              x(n,1,k) y(n,1,k) x(n,n,k) y(n,n,k) n ;  % "Down" side
              x(1,1,k) y(1,1,k) x(1,n,k) y(1,n,k) n ]; % "Up" side
    xee = x(ee); xee = xee(ss);
    yee = y(ee); yee = yee(ss);
    xy = [xee yee];
    Ainv = @(u) Aii(:,:,k) \ u;
    L{k} = StripSolver.Leaf(dom(k), S(:,:,k), D2N(:,:,k), edges, xy, Ainv);
end

end
