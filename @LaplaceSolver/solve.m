function [u, sigma] = solve(L, bc, tol)
%SOLVE   Solve a Laplace problem.

if ( nargin < 3 )
    tol = L.tol;
end

if ( size(bc, 1) ~= L.n )
    error('Incorrect dimensions.');
end

% Solve for the unknown density on the boundary
if ( strcmpi(L.method, 'direct') )
    sigma = L.op \ bc;
else
    maxiter = min(L.n, 500);
    %[sigma, ~] = gmres(L.op, bc, [], tol, maxiter);
    sigma = gmres(L.op, bc, [], tol, maxiter);
end

% Instantiate the layer potential representation with the computed density
u = L.rep(sigma);

end