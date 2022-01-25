function [tt_s, rr_s, globalIdx] = assignStripPoints(S, xx_s, yy_s)

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
