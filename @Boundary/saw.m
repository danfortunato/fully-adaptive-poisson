function C = saw(N, varargin)
%SAW   Set up a saw domain.

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'N', @isfloat);
parse(p, N, varargin{:});

R = 1; % radius (note that it is not normalized)
k = 5; % number of arms
c = 0.5; % depth and length of arms
s = 1; % orientation
b = 1; % Set b = 1 for simplicity
center = 0+0i;

z = @(t) center + R * ((b+c*sin(s*k*t)).*cos(s*t+c*sin(k*s*t)) + 1i*((b+c*sin(s*k*t)).*sin(s*t + c*sin(s*k*t))));
syms t
dz  = matlabFunction(diff(z(t)));
dzz = matlabFunction(diff(dz(t)));

r = randnfun(0.05, [0 2*pi], 'trig');
r = r - min(r);
r = r ./ max(r);

% z = chebfun(z, [0 2*pi], 'trig');
% zn = normal(z, 'unit');
% zn = zn(:,1) + zn(:,2)*1i;
% z = z + 0.1*r.*zn;

%z = chebfun(@(t) 0.4*r(t) + z(t), [0 2*pi], 'trig');

levelset = [];
C = Boundary(N, z, dz, dzz, levelset, varargin{:});

end
