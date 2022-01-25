boxes = struct();
boxes(1).domain = [-1 1 -1 1];
boxes(1).level = 0;
boxes(1).height = 1;
boxes(1).id = 1;
boxes(1).parent = 0;
boxes(1).children = [2 3 4 5];
boxes(1).coeffs = [];
boxes(1).col = 1;
boxes(1).row = 1;

n = 16;
dom = [-1 0 -1 0;
        0 1 -1 0;
       -1 0  0 1;
        0 1  0 1];
for k = 1:4
    [xx, yy] = chebpts2(n, n, dom(k,:));
    vals = sin(2*xx).*sin(2*yy)-yy.^2.*exp(xx);
    coeffs = chebtech2.vals2coeffs( chebtech2.vals2coeffs(vals).' ).';
    boxes(k+1).domain = dom(k,:);
    boxes(k+1).level = 1;
    boxes(k+1).height = 0;
    boxes(k+1).id = k+1;
    boxes(k+1).parent = 1;
    boxes(k+1).coeffs = coeffs;
    boxes(k+1).col = mod(k-1,2)+1;
    boxes(k+1).row = double(k>2)+1;
end

f = treefun2(boxes);
plot(f)

%%
% Convert second-kind Chebyshev points to Legendre points on each box
leaf = leaves(f);
xx = cell(length(leaf), 1);
yy = cell(length(leaf), 1);
ww = cell(length(leaf), 1);
ff = leafvals(f);
[x0, w0] = legpts(f.n, [0 1]);
[xx0, yy0] = meshgrid(x0);
ww0 = w0(:) * w0(:).';
for k = 1:length(leaf)
    dom = leaf(k).domain;
    sclx = diff(dom(1:2));
    scly = diff(dom(3:4));
    xx{k} = sclx*xx0 + dom(1);
    yy{k} = scly*yy0 + dom(3);
    ww{k} = sclx*scly*ww0;
    ff{k} = chebvals2legvals( chebvals2legvals(ff{k}).' ).';
end


xx_ = [xx{:}]; xx_ = xx_(:);
yy_ = [yy{:}]; yy_ = yy_(:);
ww_ = [ww{:}]; ww_ = ww_(:);
ff_ = [ff{:}]; ff_ = ff_(:);

iprec = 4;
nsource = length(f) * f.n^2;
source = [xx_ yy_].';
ifcharge = 1;
charge = -ww_.*ff_/(2*pi);
ifdipole = 0;
dipstr = zeros(1, nsource);
dipvec = zeros(2, nsource);
ifpot  = 1;
ifgrad = 0;
ifhess = 0;
target = [];
ntarget = length(target);
ifpottarg  = 0;
ifgradtarg = 0;
ifhesstarg = 0;

% Call a point FMM
% Note: This will not work if box points overlap (e.g., if we use
% second-kind Chebyshev points)
out = rfmm2dpart(iprec,                    ...
                 nsource, source,          ...
                 ifcharge, charge,         ...
                 ifdipole, dipstr, dipvec, ...
                 ifpot, ifgrad, ifhess,    ...
                 ntarget, target,          ...
                 ifpottarg, ifgradtarg, ifhesstarg);

uu = reshape(out.pot, f.n, f.n, length(leaf));
uu = squeeze(mat2cell(uu, f.n, f.n, ones(length(leaf), 1)));
ff_cfs = cell(length(leaf), 1);
for k = 1:length(leaf)
    ff_cfs{k} = legvals2legcoeffs( legvals2legcoeffs(ff{k}).' ).';
end

u = f;
for k = 1:length(leaf)
    uu{k} = legvals2chebcoeffs( legvals2chebcoeffs(uu{k}).' ).';
    id = leaf(k).id;
    u.boxes(id).coeffs = uu{k};
end
