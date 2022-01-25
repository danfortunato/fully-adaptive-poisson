function [u, d] = solve(P, bc)
%SOLVE   Solve a parent patch.
%   [U, D] = SOLVE(P, BC) returns a cell array U of function values
%   representing the PDE solution on the parent P with Dirichlet boundary
%   data given by BC, along with a vector D containing the domains of the
%   subpatches.

if ( ~isnumeric(bc) )
    % Evaluate the RHS if given a function handle:
    bc = feval(bc, P.xy(:,1), P.xy(:,2));
elseif ( isscalar(bc) )
    % Convert a scalar to a constant vector:
    bc = repmat(bc, size(P.xy, 1), 1);
end

% Evaluate the solution operator for the parent:
u = P.S * [bc ; 1]; % The 1 accounts for the particular part.

% Construct boundary conditions for children and recurse.

% Construct boundary indices:
i1 = 1:numel(P.idx1{1});
i2 = 1:numel(P.idx2{1});
if ( ~isempty(i1) )
    i2 = i2 + i1(end);
end
% CAT() is 10x faster than CELL2MAT().
idx1 = cat(1, P.idx1{:}); % idx1 = cell2mat(p.idx1.');
idx2 = cat(1, P.idx2{:}); % idx2 = cell2mat(p.idx2.');

% Assemble boundary conditions for child patches:
ubc1 = ones(size(P.child1.S, 2)-1, 1);
ubc1(idx1) = [bc(i1) ; u];
ubc2 = ones(size(P.child2.S, 2)-1, 1);
ubc2(idx2) = [bc(i2) ; u];

% Solve for the child patches:
[u1, d1] = solve(P.child1, ubc1);
[u2, d2] = solve(P.child2, ubc2);

% Concatenate for output:
u = [u1 ; u2]; 
d = [d1 ; d2];

end
