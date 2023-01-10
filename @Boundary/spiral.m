function C = spiral(N, varargin)
%SPIRAL   Set up a spiral domain.

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'N', @isfloat);
parse(p, N, varargin{:});

aa = 0.2;
bb = 0;
%r = @(t) 1 - aa*erf((t-pi)/aa);
r = @(t) 1 - aa*erf((pi/2-abs(mod(t-pi/2,2*pi)-pi))/aa);
%r = @(t) r(pi-t);

c = 1;
sabs = @(x) exp(-(x/c).^2)*c/sqrt(pi)+x.*erf(x/c); % c-smoothed absval
%sabs = @(x) sqrt(x.^2 + c^2);
%th = @(t) bb - aa + 2*(1-(bb-aa)/pi)*sabs(pi/2-abs(mod(t,2*pi)-pi));
th0 = @(t) bb - aa + 2*(1-(bb-aa)/pi)*sabs((t-pi)/2);
thL = @(t) pi - th0(t+pi);
thR = @(t) pi - th0(t-pi);
%th = @(t) th(t).*(mod(t,2*pi)<=pi) + (-th(t) + 2*th(pi)).*(mod(t,2*pi)>pi);
%r = @(t) r(pi-t);
%th = @(t) th(pi-t);

dist = 0.1;
th = blend(blend(thL, th0, pi/2, dist), thR, 3*pi/2, dist);

tt = linspace(0, 2*pi, 100).';

b = 0.5;
k = 20;
m = 4;
c = 0;
a = 1./exp(k*b*2 + c);

z = @(t) r(t).*a*exp(k*(b+m*1i)*th(t)+c);

%syms t
%dz  = matlabFunction(diff(z(t)));
%dzz = matlabFunction(diff(dz(t)));
z = chebfun(z, [0 2*pi]);
dz = diff(z);
dzz = diff(dz);

levelset = [];
C = Boundary(N, z, dz, dzz, levelset, varargin{:});

end
