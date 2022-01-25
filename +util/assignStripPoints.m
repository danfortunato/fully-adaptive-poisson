function [tt_s, rr_s, globalIdx] = assignStripPoints(dom, xx_s, yy_s, ...
    gamx_cheb, gamy_cheb, gam1x_cheb, gam1y_cheb, gam_nx_cheb)

nthreads = 8;
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
    xfun = chebfun2(dom(k).x);
    yfun = chebfun2(dom(k).y);
    % Filter out the points that are not in the bounding box of this element
    box = [min(dom(k).x(:))-1e-12 max(dom(k).x(:))+1e-12;
           min(dom(k).y(:))-1e-12 max(dom(k).y(:))+1e-12];
    testIdx = kd.range(box);
    %[tt, rr, converged] = cart2curv_element(xfun, yfun, xx_s(testIdx), yy_s(testIdx), [0 0]);
    nx = gam_nx_cheb{k}(:,1);
    ny = gam_nx_cheb{k}(:,2);
    [tt, rr, found_k] = getElemCoordinates(xx_s(testIdx), yy_s(testIdx), ...
                    gamx_cheb{k}, gamy_cheb{k}, gam1x_cheb{k}, gam1y_cheb{k}, ...
                    nx, ny, xcheb, vcheb, nthreads);
    %goodIdx = ~found(testIdx) & converged & abs(tt) <= 1.01 & abs(rr) <= 1.01;
    goodIdx = ~found(testIdx) & found_k & abs(tt) <= 1.01 & abs(rr) <= 1.01;
    tt_s{k} = tt(goodIdx);
    rr_s{k} = rr(goodIdx);
    gidx = testIdx(goodIdx);
    globalIdx{k} = gidx;
    %globalIdx{k} = false(size(xx_s));
    %globalIdx{k}(gidx) = true;
    rr_s_global(gidx) = rr_s{k};
    found(gidx) = true;
    %found = found | globalIdx{k};
end

% Did we miss any points?
%idx = globalIdx{1};
idx = false(size(xx_s(:)));
for k = 1:length(dom)
    %idx = idx | globalIdx{k};
    idx(globalIdx{k}) = true;
end
if ( ~all(idx) )
    warning('Some grid points were not assigned to an element.');
    % figure
    % plot(Gamma), hold on, plot(Gamma1), hold on
    % scatter(xx_s(idx), yy_s(idx), 'bo'), hold on
    % scatter(xx_s(~idx), yy_s(~idx), 'ro'), hold off
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
