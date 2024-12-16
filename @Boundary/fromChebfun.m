function C = fromChebfun(n, z, varargin)

dz  = diff(z);
dzz = diff(dz);
levelset = [];
C = Boundary(n, z, dz, dzz, levelset, varargin{:});

end
