function C = c_shape(N, varargin)

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'N', @isfloat);
parse(p, N, varargin{:});

a = 0.2;
b = 0;
c = 0.1 - 0.003;
sabs = @(x) sqrt(x.^2 + c^2);
rho = @(t) 1 - a*erf((pi/2-abs(mod(t,2*pi)-pi))/a);
th = @(t) b - a + 2*(1-(b-a)/pi)*sabs(pi/2-abs(mod(t,2*pi)-pi));
th = @(t) th(t).*(mod(t,2*pi)<=pi) + (-th(t) + 2*th(pi)).*(mod(t,2*pi)>pi);
rho = @(t) rho(pi-t);
th = @(t) th(pi-t);

m = 101;
z = chebfun(@(t) rho(t).*cos(th(t)) + rho(t).*sin(th(t))*1i, [0 2*pi], 'trig');
cfs = trigcoeffs(z, m);
%cfs = [cfs(1:50); cfs(51); cfs(50:-1:1)];
%cfs = [cfs(end:-1:52); cfs(51); cfs(52:end)];
%scl = (1./(1:(m-1)/2).^2).';
%cfs = cfs .* [fliplr(scl); 1; scl];
cfs = [cfs(1)*cfs(1:(m-1)/2); cfs; cfs((m-1)/2+2:end)*cfs(end)];
z = chebfun(cfs, [0 2*pi], 'trig', 'coeffs');
dz = diff(z);
dzz = diff(dz);

%levelset = @(x,y) sqrt(x.^2+y.^2) - r(atan2(y,x));
levelset = [];
C = Boundary(N, z, dz, dzz, levelset, varargin{:});

end
