function C = raindrop(N, varargin)

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'N', @isfloat);
addParameter(p, 'm', 3, @isfloat);
addParameter(p, 's', 1/7, @isfloat);
parse(p, N, varargin{:});
m = p.Results.m;
s = p.Results.s;

z = chebfun(@(t) 1i*cos(t) - sin(t).*sin(0.5*t).^m, [0 2*pi], 'trig');
g = chebfun(@(x) exp(-0.5*(x-pi).^2/(pi*s).^2)/((pi*s)*sqrt(2*pi)), [0 2*pi], 'trig');
z = circconv(z, g);

%data = load('data/raindrop.mat', 'h');
%z = chebfun(@(t) 1i*data.h(t/pi-1), [0 2*pi], 'trig');
dz = diff(z);
dzz = diff(dz);

levelset = [];
C = Boundary(N, z, dz, dzz, levelset, varargin{:});

end
