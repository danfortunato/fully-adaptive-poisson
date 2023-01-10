function C = larrycup(N, varargin)
%LARRYCUP

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'N', @isfloat);
parse(p, N, varargin{:});

a = 0.1;
b = 0.1;

nhalf = ceil(N/2);
s = ((1:nhalf)-0.5)/nhalf * pi;  % note half-offset, needed for easy reflection abt z
r = 1 - a*erf((s-pi/2)/a);  % radius: starts at 1+a, ends at 1-a
c = a; %*(1-b/pi);  % is theta rounding scale
sabs = @(x) exp(-(x/c).^2)*c/sqrt(pi)+x.*erf(x/c); % c-smoothed absval
th = b-a + 2*(1-(b-a)/pi)*sabs(s-pi/2);
rho = r.*sin(th); z = r.*cos(th);  % theta down from z axis as in 3D cyl coords
%z = z*1.2;  % vert stretch! makes ellipse cavity
Z = [rho -rho(end:-1:1)] + 1i*[z z(end:-1:1)]; % complex coords of full curve

z = chebfun(Z(:), [0 2*pi], 'trig');
dz = diff(z);
dzz = diff(dz);

levelset = [];
C = Boundary(N, z, dz, dzz, levelset, varargin{:});

end
