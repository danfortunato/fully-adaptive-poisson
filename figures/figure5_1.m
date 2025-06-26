doExport = false; %#ok<*UNRCH>

%% Set parameters

n = 16;
tol = 1e-10;
nplot = 400;
eta = 1e-3;  % Corner rounding length scale
alpha = 0.5; % Opening angle of raindrop tip
%inset = [-0.01 0.01 -0.02 -0.005];
inset = [-0.0025/2 0.0025/2 -0.0023 -0.001/2];
cmap = 'inferno';
Timer.reset();

%% Construct exponentially clustered Gaussians for reference solution

ngauss = 30;
ref = 1;
etamin = 1e-8;

ymax = 0;
ymin = -0.8;

% Exponential clustering into the corner
cy = logspace(log10(etamin), 0, ceil(ref*ngauss));
cy = ymax - cy*(ymax-ymin);
cx = zeros(size(cy));

s = -1e2./logspace(log10(etamin), 0, ceil(ref*ngauss));
s = -0.003*s.^2;

cx = reshape(cx, 1, 1, 1, []);
cy = reshape(cy, 1, 1, 1, []);
s  = reshape(s,  1, 1, 1, []);

bg = @(x,y) 2*cos(3*pi*y).*x;
sol = @(x,y) sum(exp(s.*((x-cx).^2+(y-cy).^2)), 4) + bg(x,y);
syms x y
f = matlabFunction(laplacian(sol(x,y)), 'vars', [x y]);
f = @(x,y) f(x,y) + 0*x;
bc = @(x,y) sol(x,y);

%% Create raindrop geometry

Gamma = Boundary.raindrop2(n, quadrature='panel', eta=eta, alpha=alpha, tol=tol);
Gamma = levelrestrict(Gamma, fac=2.1);

%% Compute solution

opts = [];
opts.debug = true;
opts.n_re  = 2*n;
opts.n_sem = [n 2*n];
opts.n_box = n;
opts.beta = 1;
opts.stripWidth = 0.5;
opts.boxToStripRatio = 0.5;
opts.tol = tol;
solve_opts = [];
solve_opts.tol = tol;

Timer.tic();
S = AdaptivePoissonSolver(Gamma, f, opts);
Timer.toc('Initialization');

Timer.tic();
u = S.solve(bc, solve_opts);
Timer.toc('Solve');

Timer.tic();
[uu, xx, yy, inGamma, inGamma1, inStrip] = S.sampleUniform(u, nplot);
Timer.toc('Evaluate plot');

Timer.tic();
[uuu, xxx, yyy, iinGamma, iinGamma1, iinStrip] = S.sampleUniform(u, nplot, inset);
Timer.toc('Evaluate inset');

Timer.tic();
[u_bulk, u_strip, err_bulk, err_strip] = S.evaluatePotentials(u, sol);
Timer.toc('Evaluate back on self');

uu_sol = nan(nplot);
uu_sol(inGamma) = sol(xx(inGamma),yy(inGamma));

max(max(abs(uu - uu_sol)))

%% Figure 5.1a

figure(1), clf
ff = nan(nplot);
ff(inGamma) = sol(xx(inGamma),yy(inGamma));

surf(xx, yy, 0*ff, ff)
view(2)

hold on
plot(S.Gamma,  'k-', markersize=18)
plot(S.Gamma1, 'k--', endpoints=false)
hold off

colormap(cmap)
axis equal off
shading interp

cb = colorbar;
clim([-1 1])
cb.Ticks = [-1 0 1];
cb.TickLabels = cellfun(@(s) ['$' s '$'], cb.TickLabels, 'UniformOutput', false);
cb.TickLabelInterpreter = 'latex';
cb.TickDirection = 'out';
cb.Location = 'WestOutside';
cb.FontSize = 18;
cb.LineWidth = 1;
cb.Color = [0 0 0];
cb.Label.Color = [0 0 0];
% pos = cb.Position;
% %pos(1) = 0.7;
% pos(3) = 0.04;
% cb.Position = pos;

pos = cb.Position;
pos(1) = 0.28;
pos(2) = 0.25;
pos(3) = 0.02;
pos(4) = 0.4;
cb.Position = pos;

shg

if ( doExport )
    export_fig -transparent -r400 raindrop_solution.png
end

%% Figure 5.1a inset

figure(2), clf

fff = nan(nplot);
fff(iinGamma) = sol(xxx(iinGamma),yyy(iinGamma));

pcolor(xxx, yyy, fff)
shading interp
view(2)

style = {'markersize', 50, 'linewidth', 5};
%style = {};

hold on
plot(S.Gamma,  'k-', style{:})
plot(S.Gamma1, 'k--', style{:}, endpoints=false)
hold off

colormap(cmap)
clim([-1 1])
axis equal
axis(inset)
axis off
shg

if ( doExport )
    export_fig -transparent -r400 raindrop_solution_inset.png
end

%% Figure 5.1b

figure(3), clf

ids = leaves(S.tf);
keep = [];
for id = ids(:).'
    if ( ~all(S.tf.coeffs{id} == 0, 'all') )
        keep = [keep id];
    end
end

plotOnly(S.tf, keep, boxes=true, values=false)
axis equal off

if ( doExport )
    exportgraphics(gcf, 'raindrop_tree.eps', 'BackgroundColor', 'none', 'ContentType', 'vector')
end

%% Figure 5.1b inset

figure(4), clf
plotOnly(S.tf, keep, boxes=true, values=false, linewidth=2)
axis equal
axis(inset)
axis off

if ( doExport )
    exportgraphics(gcf, 'raindrop_tree_inset.eps', 'BackgroundColor', 'none', 'ContentType', 'vector')
end

%% Figure 5.1c

stylefig = @(cb) chain(...
    @() set(gca, 'fontsize', fs), ...
    @() set(gca, 'linewidth', lw), ...
    @() set(gcf, 'Color', [1 1 1]), ...
    @() colormap('inferno'), ...
    @() axis('equal'), ...
    @() axis('tight'), ...
    @() axis(dom), ...
    @() axis('off'), ...
    @() shading('interp'), ...
    @() view(2), ...
    @() set(cb, 'location', 'EastOutside'), ...
    @() set(cb, 'fontsize', fs), ...
    @() set(cb, 'linewidth', 2), ...
    @() set(cb, 'TickLabelInterpreter', 'latex'), ...
    @() set(cb, 'TickLabels', ticks2tex(cb.TickLabels)), ...
    @() set(cb, 'TickDirection', 'out'), ...
    @() set(cb, 'Color', [0 0 0]), ...
    @() set(cb.Label, 'Color', [0 0 0]));

figure(5), clf
err = abs(uu-uu_sol);
err(abs(err) < 1e-16) = 1e-16;
err = log10(err ./ max(max(abs(uu_sol))));
surf(xx, yy, 0*uu, err)
view(2)
shading interp
axis equal
%axis(S.domain)
axis off
colormap inferno

hold on
plot(Gamma, 'k-', endpoints=false, linewidth=1.5)
hold off
cb = colorbar;
clim([-16 -10])
cb.Ticks = -16:2:-10;
cb.TickLabels = cellfun(@(s) ['$10^{' s '}$'], cb.TickLabels, 'UniformOutput', false);
cb.TickLabelInterpreter = 'latex';
cb.TickDirection = 'out';
cb.Location = 'EastOutside';
cb.FontSize = 18;
cb.LineWidth = 1;
cb.Color = [0 0 0];
cb.Label.Color = [0 0 0];
pos = cb.Position;
pos(1) = 0.72;
pos(2) = 0.25;
pos(3) = 0.02;
pos(4) = 0.4;
cb.Position = pos;

shg

if ( doExport )
    export_fig -transparent -r400 raindrop_error.png
end

%% Figure 5.1c inset

uuu_sol = nan(nplot);
uuu_sol(iinGamma) = sol(xxx(iinGamma),yyy(iinGamma));

figure(6), clf
err = abs(uuu-uuu_sol);
err(abs(err) < 1e-16) = 1e-16;
err = log10(err ./ max(max(abs(uu_sol))));
surf(xxx, yyy, 0*uuu, err)
view(2)
hold on
plot(Gamma, 'k-', endpoints=false, linewidth=3)
hold off
shading interp
axis equal
axis(inset)
axis off
colormap inferno
clim([-16 -10])

shg

if ( doExport )
    export_fig -transparent -r400 raindrop_error_inset.png
end

alignfigs