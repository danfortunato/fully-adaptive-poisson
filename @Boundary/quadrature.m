function [x, w, breaks] = quadrature(N, rule, np, resfun)
%QUADRATURE   Quadrature rule on [0, 2pi].

if ( nargin < 3 )
    np = 1;
    if ( nargin < 2 )
        rule = 'ptr';
    end
end

switch rule
    case 'ptr'
        x{1} = 2*pi/N*(0:N-1).';
        w{1} = 2*pi/N*ones(N,1);
        breaks = [0 2*pi];
    case 'panel'
        if ( strcmp(np, 'auto') )
            % Adaptive panels
            if ( isa(resfun{1}, 'chebfun') )
                dom = domain(resfun{1});
                a = dom(1);
                b = dom(end);
            else
                a = 0;
                b = 2*pi;
            end
            breaks = split(resfun, a, b, N);
        else
            % Uniform panels
            breaks = 2*pi*(0:np)/np;
        end
        [x01, w01] = util.gauss(N, 0, 1);
        np = numel(breaks)-1;
        x = cell(np,1);
        w = cell(np,1);
        for k = 1:np
            x{k} = (breaks(k+1)-breaks(k))*x01 + breaks(k);
            w{k} = (breaks(k+1)-breaks(k))*w01;
        end
    otherwise
        error('Unknown quadrature rule.');
end

end

%%
function t = split(f, a, b, N, tol)
%SPLIT   Adaptively split interval [a,b] according to function(s).

if ( nargin < 5 )
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
%if ( ~resolved )
    t1 = split(f, a, (a+b)/2, N, tol);
    t2 = split(f, (a+b)/2, b, N, tol);
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
%         t1 = split(f(k), a, (a+b)/2, N, tol);
%         t2 = split(f(k), (a+b)/2, b, N, tol);
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

min(norm(relerr,inf), norm(err,inf))

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

persistent x01 V2C nstored
if ( isempty(x01) || N ~= nstored )
    nstored = N;
    x01 = legpts(2*N, [0 1]);
    V2C = legvals2legcoeffs(eye(2*N));
end
x = a + x01*(b-a);

nf = length(f);
vals = zeros(2*N, nf);
for k = 1:nf
    vals(:,k) = abs(f{k}(x));
end
coeffs = V2C * vals;

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
