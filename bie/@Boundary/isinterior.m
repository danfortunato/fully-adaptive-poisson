function ii = isinterior(S, x, y) %#ok<*UNRCH>
%ISINTERIOR   Find interior points.

if ( ~all(size(x) == size(y)) )
    error('Number of X and Y points must match.');
end

if ( isempty(x) || isempty(y) )
    ii = zeros(size(x));
    return
end

useFMM = false;

if ( ~isempty(S.levelset) )
    ii = S.levelset(x,y) < 10*eps;
elseif ( useFMM )
    %inS  = kernels.laplace.dlp(S, 'target', [x(:) y(:)], 'density', sigma, 'closeeval', true, 'side', 'i');
    %outS = kernels.laplace.dlp(S, 'target', [x(:) y(:)], 'density', sigma, 'closeeval', true, 'side', 'e');
    %inS = abs(inS+1) < 1e-8;
    %outS = abs(outS) < 1e-8;
    %ii = inS & ~outS;
    %S = refine(S);

    sigma = ones(S.N*S.np, 1);
    fmm = kernels.laplace.dlp(S, 'target', [x(:) y(:)], 'density', sigma);
    ii = abs(fmm) > 0.5;
    ii = reshape(ii, size(x));

    %sigma = ones(S.N*S.np, 1);
    %fmm = kernels.laplace.dlp(S, 'target', [x(:) y(:)], 'density', sigma, 'closeeval', true, 'side', 'i');
    %ii = abs(fmm+1) < 1e-12;
    %ii = reshape(ii, size(x));
else
    %warning('Level set not found. Defaulting to polygon.');
    %ii = isinterior(S.polygon, x(:), y(:));
    %ii = reshape(ii, size(x));

    ii = inpoly2([x(:) y(:)], S.polynodes, S.polyedges);
    ii = reshape(ii, size(x));
end

end
