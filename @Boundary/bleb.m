function C = bleb(N, varargin)
%BLEB

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'N', @isfloat);
parse(p, N, varargin{:});

vars = load('data/bleb_z.mat');
z = chebfun(@(t) vars.zfun_conv(t/pi-1), [0 2*pi], 'trig');
z = z - mean(z);
z = z ./ max(abs(z));

dz = diff(z);
dzz = diff(dz);

levelset = [];
C = Boundary(N, z, dz, dzz, levelset, varargin{:});

end
