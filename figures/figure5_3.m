doExport = false; %#ok<*UNRCH>

%% Set parameters

tol = 1e-9;
n = 10; % Scale n with 1/eps
nplot = 400;
geom_tol = 1e-12;
align = 1;
amp   = 1;
delta = 1e-4;
freq  = 3;
gap   = 0.3;
inset = [-1e-4 3e-4 -2e-4 2e-4];

cmap = inferno(512);
cmap = cmap(60:end,:);

%% Make the geometry

% Exact solution
L = log(1/delta);
be = L*freq;
t1 = @(r) pi*(1 + log(r)/L);
mouthang = @(r) pi + 0.5*gap*amp/freq - amp/freq*sin(be*t1(r)).^2;
r = @(x,y) sqrt(x.^2 + y.^2);
real_scl = @(x,y) x.*cos(mouthang(r(x,y))) - y.*sin(mouthang(r(x,y)));
imag_scl = @(x,y) x.*sin(mouthang(r(x,y))) + y.*cos(mouthang(r(x,y)));
theta = @(x,y) atan2(imag_scl(x,y), real_scl(x,y)) - mouthang(r(x,y));
sol = @(x,y) 1/16*(x.^2+y.^2) .* (theta(x,y) + pi);
f = @(x,y) 1/16*4*(theta(x,y) + pi);
bc = @(x,y) sol(x,y);

Gamma = Boundary.gobbler(n, quadrature='panel', align=align, amp=amp, ...
    delta=delta, freq=freq, gap=gap, tol=geom_tol);

% Further refinement
Gamma = kdrefine(Gamma);
Gamma = levelrestrict(Gamma, fac=2, maxsize=0.5);

%% Do a run for a given tolerance

Gamma = resample(Gamma, n);

Timer.reset();

% Initialize the solver
opts = [];
opts.debug = true;
opts.n_sem = [n 2*n];
opts.n_re  = 2*n;
opts.n_box = n;
opts.beta = 1;
opts.average = 5;
opts.stripWidth = 0.5;
opts.tol = tol;
opts.boxToStripRatio = 0.5;

% Initialize the solver
Timer.tic();
S = AdaptivePoissonSolver(Gamma, f, opts);
Timer.toc('Initialize');

% Solve
Timer.tic();
u = S.solve(bc);
Timer.toc('Solve');

% Evaluate
Timer.tic();
dom = boundingbox(S.Gamma);
[uu, xx, yy, inGamma, inGamma1, inStrip] = S.sampleUniform(u, nplot, dom);
Timer.toc('Evaluate on uniform grid');

Timer.tic();
[uuu, xxx, yyy, iinGamma, iinGamma1, iinStrip] = S.sampleUniform(u, nplot, inset);
Timer.toc('Evaluate inset');

% Calculate error
uu_sol = nan(nplot);
uu_sol(inGamma) = sol(xx(inGamma),yy(inGamma));
max(abs(uu(:) - uu_sol(:)))

%% Figure 5.3a

figure(1), clf
ff = nan(nplot);
ff(inGamma) = f(xx(inGamma),yy(inGamma));

surf(xx, yy, 0*ff, ff)
view(2)

hold on
plot(S.Gamma,  'k-',  linewidth=2.5, markersize=20)
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
cb.FontSize = 24;
cb.LineWidth = 1.5;
cb.Color = [0 0 0];
cb.Label.Color = [0 0 0];
% pos = cb.Position;
% %pos(1) = 0.7;
% pos(3) = 0.04;
% cb.Position = pos;

pos = cb.Position;
pos(1) = 0.16;
pos(2) = 0.25;
pos(3) = 0.028;
pos(4) = 0.5;
cb.Position = pos;

shg

if ( doExport )
    export_fig -transparent -r400 gobbler_solution.png
end

%% Figure 5.3a inset

figure(2), clf

fff = nan(nplot);
fff(iinGamma) = f(xxx(iinGamma),yyy(iinGamma));

pcolor(xxx, yyy, fff)
shading interp
view(2)

style = {'markersize', 24, 'linewidth', 3};

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
    export_fig -transparent -r400 gobbler_solution_inset.png
end

%% Figure 5.3b

figure(3), clf

ids = leaves(S.tf);
keep = [];
for id = ids(:).'
    if ( ~all(S.tf.coeffs{id} == 0, 'all') )
        keep = [keep id]; %#ok<*AGROW>
    end
end

plotOnly(S.tf, keep, boxes=true, values=false);

axis equal off
shg

if ( doExport )
    exportgraphics(gcf, 'gobbler_tree.eps', 'BackgroundColor', 'none', 'ContentType', 'vector')
end

%% Figure 5.3b inset

figure(4), clf
plotOnly(S.tf, keep, boxes=true, values=false, linewidth=2)
axis equal
axis(inset)
axis off

if ( doExport )
    exportgraphics(gcf, 'gobbler_tree_inset.eps', 'BackgroundColor', 'none', 'ContentType', 'vector')
end

%% Figure 5.3c

figure(5), clf

surf(xx, yy, 0*uu, log10(abs(uu-uu_sol)))
view(2)

hold on
plot(S.Gamma, 'k-', endpoints=false, linewidth=2)
hold off

colormap(cmap)

axis equal off
shading interp

cb = colorbar;
clim([-16 -10])
cb.Ticks = -16:2:-10;
cb.TickLabels = cellfun(@(s) ['$10^{' s '}$'], cb.TickLabels, 'UniformOutput', false);
cb.TickLabelInterpreter = 'latex';
cb.TickDirection = 'out';
cb.Location = 'EastOutside';
cb.FontSize = 24;
cb.LineWidth = 1.5;
cb.Color = [0 0 0];
cb.Label.Color = [0 0 0];
pos = cb.Position;
pos(1) = 0.85;
pos(2) = 0.25;
pos(3) = 0.028;
pos(4) = 0.5;
cb.Position = pos;

shg

if ( doExport )
    export_fig -transparent -r400 gobbler_error.png
end

%% Figure 5.3c inset

figure(6), clf

uuu_sol = nan(nplot);
uuu_sol(iinGamma) = sol(xxx(iinGamma),yyy(iinGamma));

err = abs(uuu-uuu_sol);
err(abs(err) < 1e-16) = 1e-16;
surf(xxx, yyy, 0*uuu, log10(err))
view(2)
hold on
plot(Gamma, 'k-', endpoints=false, linewidth=3)
hold off
shading interp
axis equal
axis(inset)
axis off
colormap(cmap)
clim([-16 -10])

shg

if ( doExport )
    export_fig -transparent -r400 gobbler_error_inset.png
end

alignfigs
