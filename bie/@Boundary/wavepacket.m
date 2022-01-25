function C = wavepacket(N, varargin)

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'N', @isfloat);
parse(p, N, varargin{:});

%w = 1e6;
w = 50;

syms x
packet_f = @(x) 1/sqrt(w) * cos(4*pi*sqrt(w)*x) .* exp(-w*x.^2);
r = @(x) 1 + packet_f((x-pi/2)/pi);
dr = matlabFunction(diff(r(x)));
drr = matlabFunction(diff(dr(x)));

%packet = @(x) 1/sqrt(w) * cos(2*pi*sqrt(w)*x) .* exp(-w*x.^2);
%packet = @(x) packet((x-pi/2)/pi);
%r = chebfun(@(x) 1 + packet(x), [0 2*pi], 'trig');
%dr = diff(r);
%drr = diff(r, 2);

z   = @(t) r(t).*cos(t) + r(t).*sin(t)*1i;
dz  = @(t) (dr(t).*cos(t) - r(t).*sin(t)) + ...
           (dr(t).*sin(t) + r(t).*cos(t))*1i;
dzz = @(t) (drr(t).*cos(t) - 2*dr(t).*sin(t) - r(t).*cos(t)) + ...
           (drr(t).*sin(t) + 2*dr(t).*cos(t) - r(t).*sin(t))*1i;
levelset = @(x,y) sqrt(x.^2+y.^2) - r(atan2(y,x));
C = Boundary(N, z, dz, dzz, levelset, varargin{:});

end
