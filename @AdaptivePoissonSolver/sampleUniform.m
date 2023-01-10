function [uu, xx, yy, inGamma, inGamma1, inStrip] = sampleUniform(S, u, nbulk, dom)

if ( nargin < 3 )
    nbulk = 400;
end

if ( nargin < 4 )
    dom = S.domain;
end

uu = nan(nbulk);

% Sample on a uniform grid of size nbulk x nbulk
[xx, yy] = meshgrid(linspace(dom(1), dom(2), nbulk).', ...
                    linspace(dom(3), dom(4), nbulk).');

% Partition the grid
[inGamma, inGamma1, inStrip] = S.partition(xx, yy);

% Sample the solution inside Gamma'
if ( any(inGamma1, 'all') )
    xx_g1 = xx(inGamma1);
    yy_g1 = yy(inGamma1);
    uu(inGamma1) = u.bulk(xx_g1, yy_g1) + u.glue.int(xx_g1, yy_g1) + u.bc(xx_g1, yy_g1);
end

% Sample the solution in the strip
if ( any(inStrip, 'all') )
    xx_s = xx(inStrip);
    yy_s = yy(inStrip);
    stripvals = zeros(size(xx_s));
    [tt_s, rr_s, globalIdx, notFound] = assignStripPoints(S, xx_s, yy_s);
    for k = 1:length(S.strip_dom)
        stripvals(globalIdx{k}) = util.bary2d(u.strip(:,:,k), tt_s{k}, rr_s{k});
    end
    uu(inStrip) = stripvals + u.glue.ext(xx_s, yy_s) + u.bc(xx_s, yy_s);
    uu_s = uu(inStrip);
    uu_s(notFound) = nan;
    uu(inStrip) = uu_s;
end

uu = real(uu);

end

function [tt_s, rr_s, globalIdx, notFound] = assignStripPoints(S, xx_s, yy_s)

nthreads = 8;
dom = S.strip_dom;
n = size(dom(1).x, 2);
[xcheb, ~, vcheb] = chebpts(n);
getElemCoordinates = util.getElemCoordinates_(n);

p = [xx_s yy_s];
kd = KDTree(p);
tt_s = cell(length(dom), 1);
rr_s = cell(length(dom), 1);
rr_s_global = zeros(size(xx_s));
found = false(size(xx_s));
globalIdx = cell(length(dom), 1);
for k = 1:length(dom)
    % Filter out the points that are not in the bounding box of this element
    box = [min(dom(k).x(:))-1e-12 max(dom(k).x(:))+1e-12;
           min(dom(k).y(:))-1e-12 max(dom(k).y(:))+1e-12];
    testIdx = kd.range(box);
    nx = S.gam_nx_cheb{k}(:,1);
    ny = S.gam_nx_cheb{k}(:,2);
    [tt, rr, found_k] = getElemCoordinates(xx_s(testIdx), yy_s(testIdx), ...
      S.gamx_cheb{k}, S.gamy_cheb{k}, S.gam1x_cheb{k}, S.gam1y_cheb{k}, ...
      nx, ny, xcheb, vcheb, nthreads);
    %xfun = chebfun2(dom(k).x);
    %yfun = chebfun2(dom(k).y);
    %[tt, rr, found_k] = util.cart2curv_element(xfun, yfun, xx_s(testIdx), yy_s(testIdx), [0 0]);
    goodIdx = ~found(testIdx) & found_k & abs(tt) <= 1.01 & abs(rr) <= 1.01;
    tt_s{k} = tt(goodIdx);
    rr_s{k} = rr(goodIdx);
    gidx = testIdx(goodIdx);
    globalIdx{k} = gidx;
    rr_s_global(gidx) = rr_s{k};
    found(gidx) = true;
end

% Did we miss any points?
idx = false(size(xx_s(:)));
for k = 1:length(dom)
    idx(globalIdx{k}) = true;
end

notFound = [];
if ( ~all(idx) )
    warning('Some grid points were not assigned to an element.');
    notFound = find(~idx);
%     figure
%     plot(S.Gamma), hold on, plot(S.Gamma1), hold on
%     scatter(xx_s(idx), yy_s(idx), 'bo'), hold on
%     scatter(xx_s(~idx), yy_s(~idx), 'ro'), hold off
%     keyboard
end

% % Plot the points for each element
% plot(Gamma), hold on
% plot(Gamma1)
% for k = 1:length(dom)
%     idx = globalIdx{k};
%     scatter(xx_s(idx), yy_s(idx)), hold on
% end
% hold off
% shg

end
