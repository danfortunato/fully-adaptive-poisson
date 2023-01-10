function C = spikey(N, varargin)
%SPIKEY   Set up a spikey domain.

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'N', @isfloat);
parse(p, N, varargin{:});

r = randnfun(0.1, [0 2*pi], 'trig');
r = r - min(r);
r = r ./ max(r);
%r = chebfun(@(x) 100*cos(x), [0 2*pi], 'trig');
%r = r/5 + 1;
%r = r+1;

a = 0.3;
k = 5;
th = 0;
rad = 1;
star = chebfun(@(t) rad + a*cos(k*(t-th)), [0 2*pi], 'trig');
r = 0.4*r + star;

dr = diff(r);
drr = diff(r, 2);

z   = @(t) r(t).*cos(t) + r(t).*sin(t)*1i;
dz  = @(t) (dr(t).*cos(t) - r(t).*sin(t)) + ...
           (dr(t).*sin(t) + r(t).*cos(t))*1i;
dzz = @(t) (drr(t).*cos(t) - 2*dr(t).*sin(t) - r(t).*cos(t)) + ...
           (drr(t).*sin(t) + 2*dr(t).*cos(t) - r(t).*sin(t))*1i;

levelset = @(x,y) sqrt(x.^2+y.^2) - r(atan2(y,x));

C = Boundary(N, z, dz, dzz, levelset, varargin{:});

end
