function h = mycircconv(f, g)
%CIRCCONV   Circular convolution of CHEBFUN objects.
%   H = CONV(F, G) produces the convolution of CHEBFUN objects F and G:
%                     -
%                    /
%           H(x) =   |    F(t) G(x-t, F(x-t)) dt,  x in [-1, 1]
%                    /
%                   -
%   The integral is taken over all t in [-1, 1]. The breakpoints of H are
%   all pairwise sums of the breakpoints of F and G.

% Return empty for an empty input:
if ( isempty(f) || isempty(g) )
    h = chebfun();
    return
end

% Return a warning if F and G have too many pieces (the computation is
% probably going to be very slow):
if ( numel(f.funs) > 50 ) 
   warning('CHEBFUN:CHEBFUN:circconv:piecewise',...
       ['Convolving CHEBFUNs with many pieces can be very slow.\n', ...
        'Try calling MERGE() on the inputs before calling CIRCCONV().']);
end

% Extract the domain:
[a, b] = domain(f);
c = a;
d = b;
if ( isa(g, 'chebfun') )
    [c, d] = domain(g);
    if ( abs(a-c) > eps ||  abs(b-d) > eps )
        error('CHEBFUN:circconv:domain', 'Domains of f and g must match.');
    end
end

% No support for unbounded domains:
if ( any(isinf([a b c d])) )
    error('CHEBFUN:CHEBFUN:circconv:bounded', ...
        'CIRCCONV only supports CHEBFUN objects on bounded domains.');
end

if ( isa(g, 'function_handle') )
    g = periodize(g, [a b]);
end

h = oldConv(f, g);

% Simplify the result:
%h = simplify(h);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g = periodize(g, dom)

t = @(x) mod(x+dom(2), dom(2)-dom(1)) + dom(1);
g = @(x, varargin) feval(g, t(x), varargin{:});

end

function h = oldConv(f, g)
% The old convolution algorithm based on quadrature.

% Set preferences:
p = chebfunpref();
p.splitting = false;
p.blowup = false;
p.tech = @trigtech;
p.techPrefs.extrapolate = true;
p.techPrefs.refinementFunction = 'nested';
p.techPrefs.sampleTest = false;
p.techPrefs.chebfuneps = 1e-8;

% Construct CHEBFUN:
h = chebfun(@(x) convIntegral(x, f, g), p);

%xx = chebpts(200);
%h = convIntegral(xx, f, g);

end

function out = convIntegral(x, f, g)
%CONVINTEGRAL   Evaluate convolution integral.
%   Y = CONVINTEGRAL(X, F, G) evaluates the convolution of the CHEBFUNs F
%   and G at the points X.

dom = f.domain;
out = 0*x;
ff = f(x);
for k = 1:length(x)
    %integrand = @(t) feval(f, t) .* feval(g, x(k)-t, ff(k));
    integrand = @(t) feval(f, t) .* feval(g, x(k)-t);
    for j = 1:length(dom)-1
        out(k) = out(k) + integral(integrand, dom(j), dom(j+1), ...
            'AbsTol', 1e-15, 'RelTol', 1e-15);
    end
end

end
