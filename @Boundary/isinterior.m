function ii = isinterior(S, x, y, varargin) %#ok<*UNRCH>
%ISINTERIOR   Find interior points.

if ( ~all(size(x) == size(y)) )
    error('Number of X and Y points must match.');
end

if ( isempty(x) || isempty(y) )
    ii = zeros(size(x));
    return
end

parser = inputParser;
parser.KeepUnmatched = true;
addRequired(parser, 'S');
addRequired(parser, 'x');
addRequired(parser, 'y');
addParameter(parser, 'fmm', false, @islogical);
parse(parser, S, x, y, varargin{:});
useFMM = parser.Results.fmm;

if ( ~isempty(S.levelset) )
    ii = S.levelset(x,y) < 10*eps;
elseif ( useFMM )
    sigma = ones(S.N*S.np, 1);
    fmm = kernels.laplace.dlp(S, 'target', [x(:) y(:)], 'density', sigma);
    ii = abs(fmm) > 0.5;
    ii = reshape(ii, size(x));
else
    ii = inpoly2([x(:) y(:)], S.polynodes, S.polyedges);
    ii = reshape(ii, size(x));
end

end
