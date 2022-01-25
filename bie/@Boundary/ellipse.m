function C = ellipse(N, a, b, varargin)
%ELLIPSE   Set up an ellipse.

if ( nargin == 1 )
    a = 1;
    b = 1;
end

z   = @(t)  a.*cos(t) + b.*sin(t)*1i;
dz  = @(t) -a.*sin(t) + b.*cos(t)*1i;
dzz = @(t) -a.*cos(t) - b.*sin(t)*1i;

levelset = @(x,y) sqrt((x/a).^2+(y/b).^2) - 1;

C = Boundary(N, z, dz, dzz, levelset, varargin{:});

end
