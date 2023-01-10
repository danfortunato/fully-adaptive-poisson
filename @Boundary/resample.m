function dom = resample(dom, n)

if ( nargin < 1 )
    n = 2*dom.N;
end

if ( n == dom.N )
    return
end

if ( ~strcmpi(dom.rule, 'panel') )
    error('Resampling is not implemented for this domain.');
end

[x, w] = legpts(n);
[xk, ~, vk] = legpts(dom.N);
w = w(:);
for k = 1:dom.np
    dom.x{k}         = bary(x, dom.x{k}, xk, vk);
    dom.dx{k}        = bary(x, dom.dx{k}, xk, vk);
    dom.dxx{k}       = bary(x, dom.dxx{k}, xk, vk);
    dom.z{k}         = bary(x, dom.z{k}, xk, vk);
    dom.dz{k}        = bary(x, dom.dz{k}, xk, vk);
    dom.dzz{k}       = bary(x, dom.dzz{k}, xk, vk);
    dom.s{k}         = bary(x, dom.s{k}, xk, vk);
    dom.speed{k}     = bary(x, dom.speed{k}, xk, vk);
    dom.normal{k}    = bary(x, dom.normal{k}, xk, vk);
    dom.curvature{k} = bary(x, dom.curvature{k}, xk, vk);
    wk = (dom.breaks(k+1)-dom.breaks(k))/2 * w;
    dom.w{k}         = wk .* dom.speed{k};
    dom.cw{k}        = wk .* dom.dz{k};
end
dom.N = n;

%nodes = cell2mat(dom.x);
%warnstate = warning('off', 'MATLAB:polyshape:repairedBySimplify');
%dom.polygon = polyshape(nodes(:,1), nodes(:,2));
%warning(warnstate);

nodes = cell2mat(dom.x);
npts = size(nodes, 1);
dom.polynodes = nodes;
dom.polyedges = [1:npts ; 2:npts 1].';

end
