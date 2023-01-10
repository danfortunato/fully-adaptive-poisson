function domnew = reparametrize(dom)

n = dom.N;
np = dom.np;
[s01, w01, v01] = legpts(n, [0 1]);
w01 = w01(:);

speed = chebfun(legvals2chebvals([dom.speed{:}]), [0 1]);
total = sum([dom.w{:}]);
scls = diff(dom.breaks);
f = cumsum(speed .* scls) ./ total;

% Find nodes on the curve that are Legendre in arc length:
t = zeros(n, np);
for k = 1:n
    t(k,:) = roots(f - s01(k));
end

domnew = dom;
for k = 1:np
    T = barymat(t(:,k), s01, v01);
    domnew.z{k}   = T * dom.z{k};
    domnew.dz{k}  = T * (dom.dz{k} ./ dom.speed{k});
    domnew.dzz{k} = T * ((dom.dzz{k} - dom.dz{k}./dom.speed{k} .* dom.accel{k}) ./ dom.speed{k}.^2);

    domnew.x{k}   = [ real(domnew.z{k})   imag(domnew.z{k})   ];
    domnew.dx{k}  = [ real(domnew.dz{k})  imag(domnew.dz{k})  ];
    domnew.dxx{k} = [ real(domnew.dzz{k}) imag(domnew.dzz{k}) ];

    % Arc length parametrizations have unit speed:
    domnew.speed{k} = ones(n, 1);
    domnew.accel{k} = zeros(n, 1);
    domnew.normal{k} = [domnew.dx{k}(:,2), -domnew.dx{k}(:,1)];
    domnew.curvature{k} = -dot(domnew.dxx{k}, domnew.normal{k}, 2);

    domnew.breaks(k+1) = domnew.breaks(k) + sum(dom.w{k});
    scl = domnew.breaks(k+1) - domnew.breaks(k);
    domnew.s{k}  = scl*s01 + domnew.breaks(k);
    domnew.w{k}  = scl*w01;
    domnew.cw{k} = scl*w01 .* domnew.dz{k};
end

z   = cellfun(@legvals2chebvals, domnew.z,   'UniformOutput', false);
dz  = cellfun(@legvals2chebvals, domnew.dz,  'UniformOutput', false);
dzz = cellfun(@legvals2chebvals, domnew.dzz, 'UniformOutput', false);
domnew.f   = chebfun(z,   domnew.breaks);
domnew.df  = chebfun(dz,  domnew.breaks);
domnew.dff = chebfun(dzz, domnew.breaks);

% Unnecessary:
% domnew.zbreaks = domnew.f(domnew.breaks(:)).';
% But, domnew.f/df/dff are now wrong!

nodes = cell2mat(domnew.x);
npts = size(nodes, 1);
domnew.polynodes = nodes;
domnew.polyedges = [1:npts ; 2:npts 1].';

end
