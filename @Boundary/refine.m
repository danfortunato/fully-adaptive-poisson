function newdom = refine(dom, m, ids)

if ( nargin < 2 )
    m = 1;
end

if ( nargin < 3 )
    ids = 1:dom.np;
end

newdom = dom;

if ( m == 0 )
    return
elseif ( m > 1 )
    for k = 1:m
        newdom = refine(newdom, 1, ids);
    end
    return
end

% % Refine the breaks by two
% newdom.np = 2*dom.np;
% newdom.breaks = zeros(1, newdom.np+1);
% for k = 1:dom.np
%     newdom.breaks(2*k-1) = dom.breaks(k);
%     newdom.breaks(2*k)   = (dom.breaks(k) + dom.breaks(k+1)) / 2;
% end
% newdom.breaks(end) = dom.breaks(end);

newdom.np = dom.np + length(ids);
newdom.breaks = zeros(1, newdom.np+1);
needToRefine = false(dom.np);
needToRefine(ids) = true;
j = 1;
for k = 1:dom.np
    if ( needToRefine(k) )
        newdom.breaks(j) = dom.breaks(k);
        newdom.breaks(j+1) = (dom.breaks(k) + dom.breaks(k+1)) / 2;
        j = j+2;
    else
        newdom.breaks(j) = dom.breaks(k);
        j = j+1;
    end
end
newdom.breaks(end) = dom.breaks(end);

[s01, w01] = util.gauss(dom.N, 0, 1);
for k = 1:newdom.np
    w = (newdom.breaks(k+1)-newdom.breaks(k))*w01;
    newdom.s{k}   = (newdom.breaks(k+1)-newdom.breaks(k))*s01 + newdom.breaks(k);
    newdom.z{k}   = newdom.f(newdom.s{k});
    newdom.dz{k}  = newdom.df(newdom.s{k});
    newdom.dzz{k} = newdom.dff(newdom.s{k});
    newdom.x{k}   = [ real(newdom.z{k})   imag(newdom.z{k})   ];
    newdom.dx{k}  = [ real(newdom.dz{k})  imag(newdom.dz{k})  ];
    newdom.dxx{k} = [ real(newdom.dzz{k}) imag(newdom.dzz{k}) ];
    newdom.speed{k} = sqrt(sum(newdom.dx{k}.^2,2));
    newdom.accel{k} = real(newdom.dz{k} ./ newdom.speed{k} .* conj(newdom.dzz{k}));
    newdom.normal{k} = [ newdom.dx{k}(:,2), -newdom.dx{k}(:,1) ] ./ newdom.speed{k};
    newdom.curvature{k} = -dot(newdom.dxx{k}, newdom.normal{k}, 2) ./ newdom.speed{k}.^2;
    newdom.w{k}  = w .* newdom.speed{k};
    newdom.cw{k} = w .* newdom.dz{k};
end
newdom.zbreaks = newdom.f(newdom.breaks(:)).';

%nodes = cell2mat(newdom.x);
%warnstate = warning('off', 'MATLAB:polyshape:repairedBySimplify');
%newdom.polygon = polyshape(nodes(:,1), nodes(:,2));
%warning(warnstate);

nodes = cell2mat(dom.x);
npts = size(nodes, 1);
dom.polynodes = nodes;
dom.polyedges = [1:npts ; 2:npts 1].';

end
