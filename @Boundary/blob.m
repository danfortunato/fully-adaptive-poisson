function C = blob(N, varargin)

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'N', @isfloat);
parse(p, N, varargin{:});

r = randnfun(1.5, [0 2*pi], 'trig');
%r = r - min(r);
r = 1 + 0.3 * r ./ max(min(r),max(r));
z = chebfun(@(t) r(t).*cos(t) + r(t).*sin(t)*1i, [0 2*pi], 'trig');

% vars = load('data/blob.mat');
% z = chebfun(@(t) vars.zfun_conv(t/pi-1), [0 2*pi], 'trig');
% z = z - mean(z);
% z = z ./ max(abs(z));

dz = diff(z);
dzz = diff(dz);
levelset = [];
C = Boundary(N, z, dz, dzz, levelset, varargin{:});

end
