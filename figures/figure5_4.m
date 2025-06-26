doExport = false; %#ok<*UNRCH>

%% Set parameters
rng(2)

% Random Gaussians in unit box
m = 200;
cx = 2*rand(m,1) - 1;
cy = 2*rand(m,1) - 1;

% Random strengths in [1 1e4]
s = -10.^(5*rand(m,1));
cx = reshape(cx, 1, 1, []);
cy = reshape(cy, 1, 1, []);
s  = reshape(s,  1, 1, []);
sol = @(x,y) sum(exp(s.*((x-cx).^2+(y-cy).^2)), 3);
sol = @(x,y) (sol(x,y)-4)/10;

syms x y
f = matlabFunction(laplacian(sol(x,y)), 'vars', [x y]);
f = @(x,y) f(x,y) + 0*x;
bc = @(x,y) sol(x,y);

%% Construct geometry

n = 16;
tol = 1e-10;
Gamma = Boundary.bleb(n, quadrature='panel', tol=tol);
Gamma = kdrefine(Gamma);
Gamma = levelrestrict(Gamma, fac=2);

%% Initialize the solver

opts = [];
opts.debug = true;
opts.n_sem = [n 2*n];
opts.n_re  = 2*n;
opts.n_box = 16;
opts.beta = 2;
opts.average = 5;
opts.stripWidth = 0.5;
opts.tol = tol;
opts.boxToStripRatio = 0.5;

Timer.reset();

Timer.tic();
S = AdaptivePoissonSolver(Gamma, f, opts);
Timer.toc('Initialize');

%% Solve
Timer.tic();
u = S.solve(bc);
Timer.toc('Solve');

%% Evaluate
nbulk = 800;
Timer.tic();
dom = S.domain;
[uu, xx, yy, inGamma, inGamma1, inStrip] = S.sampleUniform(u, nbulk, dom);
Timer.toc('Plot');

uu_sol = nan(nbulk);
uu_sol(inGamma) = sol(xx(inGamma),yy(inGamma));

normU = max(abs(uu_sol(:)));
max(abs(uu(:) - uu_sol(:))) ./ normU

ff = nan(nbulk);
ff(inGamma) = f(xx(inGamma),yy(inGamma));

%% Plot the solution
figure(1)
clf
surf(xx, yy, 0*uu, uu_sol)
view(2)
shading interp

hold on
plot(Gamma, 'k-', linewidth=2, markersize=18)
hold off

box = boundingbox(Gamma);
box(2) = 1.1;
axis(box)
axis off

colormap inferno
clim([0 1])
cb = colorbar;
cb.Ticks = [0 1];
cb.TickLabels = cellfun(@(s) ['$' s '$'], cb.TickLabels, 'UniformOutput', false);
cb.TickLabelInterpreter = 'latex';
cb.TickDirection = 'out';
cb.Location = 'WestOutside';
cb.FontSize = 24;
cb.LineWidth = 1.5;
cb.Color = [0 0 0];
cb.Label.Color = [0 0 0];
pos = cb.Position;
pos(1) = 0.08;
pos(2) = 0.27;
pos(3) = 0.03;
pos(4) = 0.5;
cb.Position = pos;

if ( doExport )
    export_fig -r400 -transparent bleb_solution.png
end

%% Plot the tree

figure(2), clf

ids = leaves(S.tf);
keep = [];
for id = ids(:).'
    if ( ~all(S.tf.coeffs{id} == 0, 'all') )
        keep = [keep id]; %#ok<*AGROW>
    end
end

plotOnly(S.tf, keep, boxes=true, values=false, linewidth=0.5)
axis equal tight off

if ( doExport )
    exportgraphics(gcf, 'bleb_tree.eps', 'BackgroundColor', 'none', 'ContentType', 'vector')
end

%% Plot the error

figure(3), clf

surf(xx, yy, 0*uu, log10(abs(uu-uu_sol) ./ normU))
view(2)
hold on
plot(S.Gamma, 'k-', endpoints=false)
hold off
colormap inferno

cb = colorbar;
clim([-16 -10])
cb.Ticks = -16:2:-10;
cb.TickLabels = cellfun(@(s) ['$10^{' s '}$'], cb.TickLabels, 'UniformOutput', false);
cb.TickLabelInterpreter = 'latex';
cb.TickDirection = 'out';
cb.Location = 'EastOutside';
cb.FontSize = 22;
cb.LineWidth = 1.5;
cb.Color = [0 0 0];
cb.Label.Color = [0 0 0];
pos = cb.Position;
pos(1) = 0.87;
pos(2) = 0.27;
pos(3) = 0.03;
pos(4) = 0.5;
cb.Position = pos;

axis equal off
shading interp

set(gcf, 'position', [1280 1020 620 420])

if ( doExport )
    export_fig -r400 -transparent bleb_error.png
end
