function t = split(f, a, b, N, tol)
%SPLIT   Adaptively split interval [a,b] according to function(s).

if ( nargin < 5 )
    %tol = 1e-7;
    tol = 1e-12;
    %tol = 10.^(1-N);
end

if ( isa(f, 'function_handle') || isa(f, 'chebfun') )
    f = {f};
end

%resolved = isResolved(f, a, b, N, tol);
resolved = isResolvedCoeffs(f, a, b, N, tol);
bent = isBent(f, a, b, N, tol);
if ( ~resolved || bent )
% if ( ~resolved )
    t1 = util.split(f, a, (a+b)/2, N, tol);
    t2 = util.split(f, (a+b)/2, b, N, tol);
    t = unique([a t1 t2 b]);
else
    t = [a b];
end

% nf = numel(f);
% t = cell(nf,1);
% for k = 1:nf
%     % Check resolution of function f{k}
%     resolved = isResolved(f{k}, a, b, N, tol);
%     if ( ~resolved )
%         t1 = util.split(f(k), a, (a+b)/2, N, tol);
%         t2 = util.split(f(k), (a+b)/2, b, N, tol);
%         t = unique([a t1 t2 b]);
%     else
%         t = [a b];
%     end
% end

end

function resolved = isResolved(f, a, b, N, tol)
%ISRESOLVED   Check if a function F is resolved to tolerance TOL on the
% interval [A, B] using N nodes.

M = 2*N;
%M = 10;
pts = util.gauss(N, a, b);
xx = linspace(a, b, M).';

persistent R savedN
if ( isempty(R) || savedN~=N)
    R = util.interpmat(util.gauss(N), linspace(-1,1,M));
    savedN = N;
end

nf = numel(f);
vals = zeros(N,nf);
tvals = zeros(M,nf);
for k = 1:nf
    vals(:,k) = f{k}(pts);
    tvals(:,k) = f{k}(xx);
end

% Interpolation error
err = tvals - R * vals;
relerr = err ./ sum(abs(tvals));
resolved = min(norm(relerr,inf), norm(err,inf)) < 10*M*max(1e-14, tol)/(b-a);

if ( (b-a) > 0.1 || abs(f{1}(b)-f{1}(a)) > pi/8 )
    resolved = false;
end

% if ( ~resolved )
%     disp([norm(relerr,inf) norm(err,inf) a, b, b-a])
% end

end

% function resolved = isResolved(f, a, b, N, tol)
% %ISRESOLVED   Check if a function F is resolved to tolerance TOL on the
% % interval [A, B] using N nodes.
% 
% m = 2*N;
% pts = util.gauss(N,a,b);
% vals = f(pts);
% fun = chebfun( legvals2chebvals(vals), [a b] );
% xx = linspace(a, b, m).';
% est = max(abs(f(xx) - fun(xx))) / max(abs(vals))
% %coeffs = legvals2legcoeffs(vals);
% %est = max(abs(coeffs(end-1:end))) / max(abs(vals));
% resolved = ( est < tol );
% 
% end

function resolved = isResolvedCoeffs(f, a, b, N, tol)
%ISRESOLVEDCOEFFS   Check if a function F is resolved to tolerance TOL on
% the interval [A, B] using N nodes.

% pts = util.gauss(N, a, b);
% vals = abs(f(pts));
% coeffs = legvals2legcoeffs(vals);
% est = max(abs(coeffs(end-1:end))) / max(abs(vals));
% resolved = ( est < tol );

pts = util.gauss(2*N, a, b);
nf = numel(f);
vals = zeros(2*N, nf);
for k = 1:nf
    vals(:,k) = abs(f{k}(pts));
end
coeffs = legvals2legcoeffs(vals);
errs0 = sum(coeffs(1:N,:).^2,1);
errs = sum(coeffs(N+1:2*N,:).^2,1);
rmsemax = max(sqrt(errs./errs0/N));
resolved = ( rmsemax < tol );

end

function bent = isBent(f, a, b, N, tol)

[pts, w, D] = util.gauss(N, a, b);
zp  = D   * f{1}(pts);
zpp = D^2 * f{1}(pts);
% Bending energy (normalized by 2*pi)
bending = sum(w .* imag(zpp./zp).^2 ./ abs(zp));% * (b-a);% / (2*pi);
arclen = sum(abs(zp) .* w);
bending = bending * arclen;
bent = ( bending > 0.5 );
%bent = ( bending > 5/log10(1/tol) );

end
