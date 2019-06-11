function t = split(f, a, b, N, tol)
%SPLIT   Adaptively split interval [a,b] according to function(s).

if ( nargin < 5 )
    tol = 1e-10;
end

if ( isa(f, 'function_handle') )
    f = {f};
end

nf = numel(f);
t = cell(nf,1);
for k = 1:nf
    % Check resolution of function f{k}
    vals = f{k}(chebpts(N, [a,b]).').';
    coeffs = chebtech2.vals2coeffs(vals);
    resolved = max(abs(coeffs(end-1:end)))*(b-a) < N*tol;

    if ( ~resolved )
        t1 = split(f(k), a, (a+b)/2, N, tol);
        t2 = split(f(k), (a+b)/2, b, N, tol);
        t = unique([a t1 t2 b]);
    else
        t = [a b];
    end
    
end

end
