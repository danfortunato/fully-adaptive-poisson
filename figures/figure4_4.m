doExport = false; %#ok<*UNRCH>

n = 16;
Gamma = Boundary.star(n, quadrature='panel');
Gamma = levelrestrict(Gamma);
%Gamma = reparametrize(Gamma);

% Initialize the solver
opts = [];
opts.debug = false;
opts.n_re  = 16;
opts.n_sem = [16 24];
opts.n_box = 16;
opts.beta = 1;
opts.stripWidth = 0.46;
opts.tol = 1e-12;
opts.boxToStripRatio = 0.5;
f = @(x,y) 1+0*x;
S = AdaptivePoissonSolver(Gamma, f, opts);

% Plot the strip grid
gray = 0.6*[1 1 1];
ms = 25;
lw = 2;
strip_lw = 1;
n_sem = 16;

strip_dom = buildStripGrid(S.Gamma, S.Gamma1, n_sem);

clf
hold on
k = 19;
mesh(strip_dom(k).x, strip_dom(k).y, 0*strip_dom(k).x, 'LineWidth', strip_lw, 'EdgeColor', gray, 'FaceColor', 'none')
plot(strip_dom(k).x(1,:),   strip_dom(k).y(1,:),   'k.', 'MarkerSize', ms)
plot(strip_dom(k).x(end,:), strip_dom(k).y(end,:), 'k.', 'MarkerSize', ms)
dx = strip_dom(k).x(1,1) - strip_dom(k).x(end,1);
dy = strip_dom(k).y(1,1) - strip_dom(k).y(end,1);
nn = [-dy dx] ./ sqrt(dx.^2 + dy.^2);
cx = mean(strip_dom(k).x(:,1));
cy = mean(strip_dom(k).y(:,1));
view(2)

plot(S.Gamma, 'k-', 'LineWidth', lw, endpoints=false), hold on
xx = cat(1, S.Gamma1.x{:});
plot(xx(:,1), xx(:,2), 'k--', 'LineWidth', lw)
hold off

axis([0.33 0.513 1.104 1.26])
axis off
shg

% Export
if ( doExport )
    %export_fig -transparent -r300 newfigs/sem_grid2.png
    exportgraphics(gcf, 'newfigs/sem_grid.eps', 'BackgroundColor', 'none', 'ContentType', 'vector')
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
