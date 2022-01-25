function C = star(N, varargin)
%STAR   Set up a smooth star.

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'N', @isfloat);
addParameter(p, 'arms', 5, @isfloat);
addParameter(p, 'wobble', 0.3, @isfloat);
addParameter(p, 'rotate', 0, @isfloat);
addParameter(p, 'radius', 1, @isfloat);
parse(p, N, varargin{:});
k  = p.Results.arms;
a  = p.Results.wobble;
th = p.Results.rotate;
rad = p.Results.radius;
a = a*rad;

%t = np.linspace(0.0, 2.0*np.pi, N, endpoint=False)
%c = (x+1j*y) + (r + r*a*np.cos(f*(t-rot)))*np.exp(1j*t)

r   = @(t)  rad + a*cos(k*(t-th));
dr  = @(t)   -k*a*sin(k*(t-th));
drr = @(t) -k^2*a*cos(k*(t-th));

z   = @(t) r(t).*cos(t) + r(t).*sin(t)*1i;
dz  = @(t) (dr(t).*cos(t) - r(t).*sin(t)) + ...
           (dr(t).*sin(t) + r(t).*cos(t))*1i;
dzz = @(t) (drr(t).*cos(t) - 2*dr(t).*sin(t) - r(t).*cos(t)) + ...
           (drr(t).*sin(t) + 2*dr(t).*cos(t) - r(t).*sin(t))*1i;

levelset = @(x,y) sqrt(x.^2+y.^2) - r(atan2(y,x));

C = Boundary(N, z, dz, dzz, levelset, varargin{:});

end
