function C = squircle(N, v, a, b, varargin)
%SQUIRCLE   Set up a rounded rectangle (squircle).

if ( nargin == 1 )
    v = 4;
    a = 1;
    b = 1;
end

if ( mod(v,2) == 1 )
    error('Squircle parameter must be even.');
end

r = @(t) 1 ./ ( (cos(t)/a).^v + (sin(t)/b).^v ).^(1/v);
z = chebfun(@(t) r(t).*exp(1i*t), [0 2*pi], 'trig');
dz  = diff(z);
dzz = diff(dz);

levelset = @(x,y) ((x/a).^v+(y/b).^v).^(1/v) - 1;

C = Boundary(N, z, dz, dzz, levelset, varargin{:});

end
