function t = split(f, a, b, N, tol)
%SPLIT   Adaptively split interval [a,b] according to function(s).

if ( nargin < 5 )
    tol = 1e-8;
end

if ( isa(f, 'function_handle') || isa(f, 'chebfun') )
    f = {f};
end

nf = numel(f);
t = cell(nf,1);
for k = 1:nf
    % Check resolution of function f{k}
    resolved = isResolved(f{k}, a, b, N, tol);
    if ( ~resolved )
        t1 = split(f(k), a, (a+b)/2, N, tol);
        t2 = split(f(k), (a+b)/2, b, N, tol);
        t = unique([a t1 t2 b]);
    else
        t = [a b];
    end
    
end

end

function resolved = isResolved(f, a, b, N, tol)
%ISRESOLVED   Check if a function F is resolved to tolerance TOL on the
% interval [A, B] using N nodes.

pts = gauss(N,a,b);
vals = abs(f(pts));
coeffs = legvals2legcoeffs(vals);
est = max(abs(coeffs(end-1:end))) / max(abs(vals));
resolved = ( est < tol );

end
