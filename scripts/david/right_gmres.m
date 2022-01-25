function [x, flag, relres, iter, resvec] = right_gmres(A, b, restart, tol, maxit, M, x)
%RIGHT_GMRES   GMRES with right preconditioning.

if ( nargin == 6 )
    x = [];
end

if ( isa(A, 'function_handle') )
    afun = A;
else
    afun = @(x) A*x;
end

if ( isa(M, 'function_handle') )
    mfun = M;
else
    mfun = @(x) M\x;
end

% Define new LinearOperator A*P^{-1}
APinv = @(x) afun(mfun(x));

% Solve system A*P^{-1} * y = b
[y, flag, relres, iter, resvec] = gmres(APinv, b, restart, tol, maxit, x);

% Solve system P*x = y
x = mfun(y);

end
