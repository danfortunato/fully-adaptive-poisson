function C = gobbler(N, varargin)
%GOBBLER   Set up a gobbler.

parser = inputParser;
parser.KeepUnmatched = true;
addRequired(parser, 'N', @isfloat);
addParameter(parser, 'delta', 1e-3, @isfloat);
addParameter(parser, 'align', 1,    @isfloat);
addParameter(parser, 'amp',   1,    @isfloat);
addParameter(parser, 'freq',  3,    @isfloat);
addParameter(parser, 'gap',   0.1,  @isfloat);
parse(parser, N, varargin{:});

p = struct();
p.delta = parser.Results.delta;
p.align = parser.Results.align;
p.a     = parser.Results.amp;
p.f     = parser.Results.freq;
p.g     = parser.Results.gap;
z = chebfun(@(t) puxgobbler(t,p), [0 2*pi], 'trig');
dz = diff(z);
dzz = diff(dz);

levelset = [];
C = Boundary(N, z, dz, dzz, levelset, varargin{:});

end
