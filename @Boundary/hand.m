function C = hand(N, varargin)
%HAND   Set up a smooth hand.

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'N', @isfloat);
addParameter(p, 'radius', 1, @isfloat);
parse(p, N, varargin{:});
rad = p.Results.radius;

r = @(t) rad + (1+cos(20*t)).* ((1 + sin(t))/2).^6;
syms t
dr  = matlabFunction(diff(r(t)));
drr = matlabFunction(diff(dr(t)));

z   = @(t) r(t).*cos(t) + r(t).*sin(t)*1i;
dz  = @(t) (dr(t).*cos(t) - r(t).*sin(t)) + ...
           (dr(t).*sin(t) + r(t).*cos(t))*1i;
dzz = @(t) (drr(t).*cos(t) - 2*dr(t).*sin(t) - r(t).*cos(t)) + ...
           (drr(t).*sin(t) + 2*dr(t).*cos(t) - r(t).*sin(t))*1i;

levelset = @(x,y) sqrt(x.^2+y.^2) - r(atan2(y,x));

C = Boundary(N, z, dz, dzz, levelset, varargin{:});

end
