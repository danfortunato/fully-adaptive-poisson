doExport = false; %#ok<*UNRCH>

n = 16;
n_sem = 12;
gray = 0.6*[1 1 1];
ms = 25;
lw = 3;
strip_lw = 1;

B = Boundary.wavepacket(n, quadrature='panel', w=150, k=10, amp=1);
B = levelrestrict(B);

%% Over-resolved strip
clf

width = 0.005;
z = cell(B.np, 1);
for k = 1:B.np
    nx = B.normal{k};
    nx = nx(:,1)+nx(:,2)*1i;
    nx = nx./abs(nx);
    z{k} = B.z{k} - width*nx;
end
B1 = Boundary(z);

figure(1)
clf

% Plot strip grid
strip_dom = buildStripGrid(B, B1, n_sem);
for k = 1:length(strip_dom)
    mesh(strip_dom(k).x, strip_dom(k).y, 0*strip_dom(k).x, 'LineWidth', strip_lw, 'EdgeColor', gray, 'FaceColor', 'none')
    view(2)
    hold on
end

plot(B, 'k-', 'LineWidth', lw, endpoints=false), hold on
plot(B.zbreaks, 'k.', 'MarkerSize', ms)
xx = cat(1, B1.x{:});
plot(xx(:,1), xx(:,2), 'k--', 'LineWidth', lw), hold on
%plot(B1, 'k--', endpoints=false), hold on
%plot(B1.zbreaks, 'r.', 'MarkerSize', 15)
hold off

axis off

%% Under-resolved strip

width = 0.1;
z = cell(B.np, 1);
for k = 1:B.np
    nx = B.normal{k};
    nx = nx(:,1)+nx(:,2)*1i;
    nx = nx./abs(nx);
    z{k} = B.z{k} - width*nx;
end
B2 = Boundary(z);

figure(2)
clf

% Plot strip grid
strip_dom = buildStripGrid(B, B2, n_sem);
for k = 1:length(strip_dom)
    mesh(strip_dom(k).x, strip_dom(k).y, 0*strip_dom(k).x, 'LineWidth', strip_lw, 'EdgeColor', gray, 'FaceColor', 'none')
    view(2)
    hold on
end

plot(B, 'k-', 'LineWidth', lw, endpoints=false), hold on
plot(B.zbreaks, 'k.', 'MarkerSize', ms)
xx = cat(1, B2.x{:});
plot(xx(:,1), xx(:,2), 'k--', 'LineWidth', lw), hold on
%plot(B2, 'k--', endpoints=false), hold on
%plot(B2.zbreaks, 'k.', 'MarkerSize', 15)
hold off

axis off

%% Adaptive strip

opts = [];
opts.n_sem = n_sem;
f = @(x,y) 0*x;
S = AdaptivePoissonSolver(B, f, opts);
B3 = S.Gamma1;

figure(3)
clf

% Plot strip grid
strip_dom = S.strip_dom;
for k = 1:length(S.strip_dom)
    mesh(strip_dom(k).x, strip_dom(k).y, 0*strip_dom(k).x, 'LineWidth', strip_lw, 'EdgeColor', gray, 'FaceColor', 'none')
    view(2)
    hold on
end

plot(B, 'k-', 'LineWidth', lw, endpoints=false), hold on
plot(B.zbreaks, 'k.', 'MarkerSize', ms)
xx = cat(1, B3.x{:});
plot(xx(:,1), xx(:,2), 'k--', 'LineWidth', lw), hold on
%plot(B3, 'k--', endpoints=false), hold on
%plot(B3.zbreaks, 'k.', 'MarkerSize', 15)
hold off
axis off

%% Export
alignfigs

if ( doExport )
    figure(1)
    exportgraphics(gcf, 'strip_overresolved.eps', 'BackgroundColor', 'none', 'ContentType', 'vector')
    figure(2)
    exportgraphics(gcf, 'strip_underresolved.eps', 'BackgroundColor', 'none', 'ContentType', 'vector')
    figure(3)
    exportgraphics(gcf, 'strip_adaptive.eps', 'BackgroundColor', 'none', 'ContentType', 'vector')
end

%% Now do insets

%ax = [-0.1 0.1 0.915 1.083];
ax = [-0.025 0.025 1.06 1.085];
%ax = [-0.05 0.05 1.06 1.085];
pos = [500 500 1200 600];

figure(1), axis(ax)
%set(gcf, 'Position', pos)
figure(2), axis(ax)
%set(gcf, 'Position', pos)
figure(3), axis(ax)
%set(gcf, 'Position', pos)
alignfigs

if ( doExport )
    figure(1)
    exportgraphics(gcf, 'strip_overresolved_inset.eps', 'BackgroundColor', 'none', 'ContentType', 'vector')
    figure(2)
    exportgraphics(gcf, 'strip_underresolved_inset.eps', 'BackgroundColor', 'none', 'ContentType', 'vector')
    figure(3)
    exportgraphics(gcf, 'strip_adaptive_inset.eps', 'BackgroundColor', 'none', 'ContentType', 'vector')
end

%% Helper functions

function [strip_dom, gamx_cheb, gamy_cheb, gam_nx_cheb, gam1x_cheb, ...
    gam1y_cheb, gam1_nx_cheb] = buildStripGrid(Gamma, Gamma1, n_sem)

if ( ~isscalar(n_sem) )
    nt = n_sem(1);
    nr = n_sem(2);
else
    nt = n_sem;
    nr = n_sem;
end

n  = Gamma.N;
np = Gamma.np;
t = chebpts(nt, [0 1]);
LV2CV = legvals2chebvals(eye(n));
r = chebpts(nr);
x = chebpts(n);
RE = barymat(r, x);
xx = cell(np, 1);
yy = cell(np, 1);
gamx_cheb    = cell(np, 1);
gamy_cheb    = cell(np, 1);
gam1x_cheb   = cell(np, 1);
gam1y_cheb   = cell(np, 1);
gam_nx_cheb  = cell(np, 1);
gam1_nx_cheb = cell(np, 1);
for k = 1:np
    gamx_cheb{k}    = RE * LV2CV * real(Gamma.x{k}(:,1));
    gamy_cheb{k}    = RE * LV2CV * real(Gamma.x{k}(:,2));
    gam1x_cheb{k}   = RE * LV2CV * real(Gamma1.x{k}(:,1));
    gam1y_cheb{k}   = RE * LV2CV * real(Gamma1.x{k}(:,2));
    gam_nx_cheb{k}  = RE * LV2CV * Gamma.normal{k};
    gam1_nx_cheb{k} = RE * LV2CV * Gamma1.normal{k};
    xx{k} = t.*gam1x_cheb{k}.' + (1-t).*gamx_cheb{k}.';
    yy{k} = t.*gam1y_cheb{k}.' + (1-t).*gamy_cheb{k}.';
end
strip_dom = cell2struct([xx yy], {'x','y'}, 2);

end
