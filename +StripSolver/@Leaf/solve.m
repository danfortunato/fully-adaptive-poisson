function [u, d] = solve(P, bc)
%SOLVE   Solve a leaf patch.
%   [U, D] = SOLVE(P, BC) returns instead a cell array containing the
%   solution coefficients U and a vector D containing the domain of the
%   patch.

% Extract the domain from the patch:
d = P.domain;
n = size(d.x, 1);

if ( ~isnumeric(bc) )
    % Evaluate the RHS if given a function handle:
    bc = feval(bc, P.xy(:,1), P.xy(:,2));
elseif ( isscalar(bc) )
    % Convert a scalar to a constant vector:
    bc = repmat(bc, size(P.xy, 1), 1);
end

% Evaluate the solution operator for the patch:
u = P.S * [bc ; 1]; % The 1 accounts for the particular part.
u = reshape(u, n, n);

% Return cell output for consistency with PARENT/SOLVE():
u = {u};

end
