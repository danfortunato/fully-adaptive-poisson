function C = raindrop2(N, varargin)

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'N', @isfloat);
addParameter(p, 'eta',   0.1, @isfloat);
addParameter(p, 'alpha', 0.5, @isfloat);
parse(p, N, varargin{:});
eta   = p.Results.eta;
alpha = p.Results.alpha;

z = @(t) alpha*sin(2*pi-t)/2 - 1i*sqrt(sin((2*pi-t)/2).^2 + eta.^2);
%z = @(t) z(t-pi/exp(pi)); % Shift so that the panels "see" the corner

t = sym('t');
dz  = matlabFunction(diff(z(t), 1));
dzz = matlabFunction(diff(z(t), 2));

levelset = [];
C = Boundary(N, z, dz, dzz, levelset, varargin{:});

end
